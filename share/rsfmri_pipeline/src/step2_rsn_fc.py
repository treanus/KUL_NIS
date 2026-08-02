#!/usr/bin/env python3
# Run with:  python3 step2_rsn_fc.py [options]
# NOT with bash — this is a Python script.

# Suppress matplotlib GUI before any library import to prevent display freeze.
import matplotlib
matplotlib.use("Agg")

# Limit BLAS/OpenMP threads before numpy is imported to prevent CPU saturation.
import os
os.environ.setdefault("OMP_NUM_THREADS",     "4")
os.environ.setdefault("MKL_NUM_THREADS",     "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS","4")

"""
Step 2 — RSN-level Functional Connectivity Matrix + Whole-Brain FC Maps

For each subject:
  1. Load the subject-specific Yeo17 atlas (produced by step0) and extract
     one binary mask per network (up to 17).
  2. For every run:
       a. Extract mean timeseries per RSN and compute the full
          N×N inter-network Pearson-r matrix → Fisher-Z.
       b. Treat each RSN mask as a seed and compute a whole-brain
          seed-to-voxel FC map (Pearson-r → Fisher-Z), identical in
          method to Step 1 SBA.
  3. Average matrices and maps across runs → subject-level outputs.
  4. Save both raw and thresholded (r≥0.25, min 20-voxel clusters) FC maps.

Outputs (per subject):
  analysis/rsn_fc/sub-<ID>/sub-<ID>_rsn_fc_matrix.csv
  analysis/rsn_fc/sub-<ID>/sub-<ID>_rsn_fc_matrix.npy
  analysis/rsn_fc/sub-<ID>/sub-<ID>_rsn-<name>_desc-sba_statmap.nii.gz      (× N)
  analysis/rsn_fc/sub-<ID>/sub-<ID>_rsn-<name>_desc-sbaThresh_statmap.nii.gz (× N)

Usage:
  python step2_rsn_fc.py [--subjects HV01 PT01 ...] [--hv-only] [--force]
"""

import argparse
import logging
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import ALL_SUBJECTS, HV_SUBJECTS, RSN_FC_DIR, get_profile, get_seed
from utils import (
    average_maps, extract_mean_ts, get_auditory_network_mask,
    get_language_network_mask, get_runs, get_seed_mask,
    get_yeo17_subject_masks, is_rest_run, load_bold, load_brain_mask,
    stat_maps, threshold_fc_map, whole_brain_fc, write_provenance_json,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------

def _rtag(run_type: str | None) -> str:
    return f"_runtype-{run_type}" if run_type else ""


def fc_matrix_path(subject: str, ext: str, run_type: str | None = None) -> Path:
    out_dir = RSN_FC_DIR / f"sub-{subject}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"sub-{subject}{_rtag(run_type)}_rsn_fc_matrix.{ext}"


def coupling_matrix_path(subject: str, profile: str, run_type: str | None = None) -> Path:
    out_dir = RSN_FC_DIR / f"sub-{subject}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"sub-{subject}{_rtag(run_type)}_profile-{profile}_seed_rsn_coupling.csv"


def rsn_sba_map_path(subject: str, rsn_name: str, thresholded: bool = False,
                     run_type: str | None = None) -> Path:
    out_dir = RSN_FC_DIR / f"sub-{subject}"
    out_dir.mkdir(parents=True, exist_ok=True)
    desc = "sbaThresh" if thresholded else "sba"
    return out_dir / f"sub-{subject}_rsn-{rsn_name}{_rtag(run_type)}_desc-{desc}_statmap.nii.gz"


def rsn_zstat_path(subject: str, rsn_name: str, run_type: str | None = None) -> Path:
    out_dir = RSN_FC_DIR / f"sub-{subject}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"sub-{subject}_rsn-{rsn_name}{_rtag(run_type)}_desc-zstat_statmap.nii.gz"


def rsn_neglog10p_path(subject: str, rsn_name: str, run_type: str | None = None) -> Path:
    out_dir = RSN_FC_DIR / f"sub-{subject}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"sub-{subject}_rsn-{rsn_name}{_rtag(run_type)}_desc-neglog10p_statmap.nii.gz"


def rsn_mask_path(subject: str, rsn_name: str) -> Path:
    out_dir = RSN_FC_DIR / f"sub-{subject}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"sub-{subject}_rsn-{rsn_name}_mask.nii.gz"


# ---------------------------------------------------------------------------
# FC matrix for one run
# ---------------------------------------------------------------------------

def compute_run_fc_matrix(bold_data: np.ndarray,
                           bold_ref: nib.Nifti1Image,
                           rsn_masks: list[tuple[str, nib.Nifti1Image]]
                           ) -> np.ndarray:
    """
    Compute the N×N inter-RSN Fisher-Z FC matrix for a single run.
    bold_data: pre-loaded float32 array (x, y, z, t).
    bold_ref:  the source NIfTI image, used as spatial reference for resampling.
    rsn_masks: [(name, binary_mask_img), ...] already in BOLD space.
    Returns array of shape (N, N).
    """
    n = len(rsn_masks)

    # Extract mean timeseries for every RSN (skip empty masks)
    timeseries = []
    valid_idx  = []
    for i, (name, mask_img) in enumerate(rsn_masks):
        try:
            ts = extract_mean_ts(bold_data, mask_img, bold_ref)
            timeseries.append(ts)
            valid_idx.append(i)
        except ValueError:
            log.warning("    RSN '%s' mask is empty after resampling — skipped", name)

    if len(timeseries) < 2:
        raise RuntimeError("Fewer than 2 valid RSN masks — cannot compute matrix.")

    # Stack → (N_valid, T) and compute full correlation matrix
    X   = np.stack(timeseries, axis=0)               # (N_valid, T)
    X  -= X.mean(axis=1, keepdims=True)              # demean each RSN ts
    std = np.linalg.norm(X, axis=1, keepdims=True) + 1e-12
    Xn  = X / std                                     # unit L2-norm rows (N_valid, T)
    r_mat = Xn @ Xn.T                                  # (N_valid, N_valid) Pearson r
    np.fill_diagonal(r_mat, 0.0)                      # zero diagonal before clip
    r_mat = np.clip(r_mat, -0.9999, 0.9999)
    z_mat = np.arctanh(r_mat)                         # Fisher Z

    # Place results into full N×N matrix (NaN for any skipped RSN)
    full = np.full((n, n), np.nan)
    for ii, gi in enumerate(valid_idx):
        for jj, gj in enumerate(valid_idx):
            full[gi, gj] = z_mat[ii, jj]

    return full


# ---------------------------------------------------------------------------
# Coupling helpers
# ---------------------------------------------------------------------------

def _compute_run_coupling(
    bold_data:  np.ndarray,
    bold_img:   "nib.Nifti1Image",
    seed_items: list[tuple[str, dict]],
    rsn_ts:     dict[str, np.ndarray],
    rsn_names:  list[str],
    subject:    str,
    run_id:     str,
) -> "np.ndarray | None":
    """
    Compute (n_seeds × n_rsns) Fisher-Z coupling matrix for one run.
    rsn_ts: {rsn_name: timeseries} for RSNs that succeeded this run.
    Returns array or None if all seed timeseries failed.
    """
    n_seeds = len(seed_items)
    n_rsns  = len(rsn_names)
    mat     = np.full((n_seeds, n_rsns), np.nan, dtype=np.float32)

    for i, (seed_name, seed_cfg) in enumerate(seed_items):
        try:
            seed_img = get_seed_mask(seed_name, seed_cfg, subject=subject)
            s_ts     = extract_mean_ts(bold_data, seed_img, bold_img)
        except Exception as exc:
            log.warning("    [%s|%s|%s] coupling seed TS: %s", subject, seed_name, run_id, exc)
            continue
        s_d = s_ts - s_ts.mean()
        s_n = np.linalg.norm(s_d) + 1e-12
        for j, rsn_name in enumerate(rsn_names):
            if rsn_name not in rsn_ts:
                continue
            r_ts = rsn_ts[rsn_name]
            r_d  = r_ts - r_ts.mean()
            r_n  = np.linalg.norm(r_d) + 1e-12
            r_val = np.clip(np.dot(s_d, r_d) / (s_n * r_n), -0.9999, 0.9999)
            mat[i, j] = float(np.arctanh(r_val))

    return None if np.all(np.isnan(mat)) else mat


# ---------------------------------------------------------------------------
# Per-subject processing
# ---------------------------------------------------------------------------

def _compute_group(subject: str, group_runs: list[str],
                   rsn_masks: list[tuple[str, nib.Nifti1Image]],
                   rsn_names: list[str],
                   seed_items: list[tuple[str, dict]] | None,
                   ) -> tuple[list[np.ndarray], dict[str, list], dict[str, list],
                              list[np.ndarray], list[str]]:
    """
    Run the per-run FC computation over exactly this group's runs.
    Returns (run_matrices, run_rsn_maps, run_rsn_dofs, run_coupling_mats, runs_used).
    """
    run_matrices: list[np.ndarray] = []
    run_rsn_maps: dict[str, list] = {name: [] for name, _ in rsn_masks}
    run_rsn_dofs: dict[str, list] = {name: [] for name, _ in rsn_masks}
    run_coupling_mats: list[np.ndarray] = []
    runs_used: list[str] = []

    for run_id in group_runs:
        log.info("  [%s | %s] Computing RSN FC", subject, run_id)
        try:
            bold_img  = load_bold(subject, run_id)
            mask_img  = load_brain_mask(subject, run_id)
            bold_data = bold_img.get_fdata(dtype=np.float32)
            run_dof   = bold_data.shape[-1] - 3

            fc_mat = compute_run_fc_matrix(bold_data, bold_img, rsn_masks)
            run_matrices.append(fc_mat)

            run_rsn_ts: dict[str, np.ndarray] = {}
            for rsn_name, rsn_mask_img in rsn_masks:
                try:
                    rsn_ts  = extract_mean_ts(bold_data, rsn_mask_img, bold_img)
                    run_rsn_ts[rsn_name] = rsn_ts
                    fc_map  = whole_brain_fc(bold_data, rsn_ts, mask_img, bold_img)
                    run_rsn_maps[rsn_name].append(fc_map)
                    run_rsn_dofs[rsn_name].append(run_dof)
                except Exception as exc:
                    log.warning("    [%s | %s | %s] SBA map failed: %s",
                                subject, rsn_name, run_id, exc)

            if seed_items and run_rsn_ts:
                coup = _compute_run_coupling(
                    bold_data, bold_img, seed_items, run_rsn_ts, rsn_names, subject, run_id
                )
                if coup is not None:
                    run_coupling_mats.append(coup)

            del bold_data, bold_img, mask_img
            runs_used.append(run_id)

        except Exception as exc:
            log.warning("  [%s | %s] Failed (%s) — skipping run", subject, run_id, exc)

    return run_matrices, run_rsn_maps, run_rsn_dofs, run_coupling_mats, runs_used


def _save_group_rsn_maps(subject: str, run_type: str | None,
                         run_rsn_maps: dict[str, list], run_rsn_dofs: dict[str, list],
                         runs_used: list[str],
                         fdr_q: float, min_cluster: int, bilateral: bool,
                         min_r_relative: float) -> dict[str, tuple[nib.Nifti1Image, int]]:
    """Average + save (raw/thresh/zstat/neglog10p/provenance) per-RSN maps
    for one run-type group (or the combined result when run_type=None).
    Returns {rsn_name: (avg_map, pooled_dof)} for use in group recombination."""
    out: dict[str, tuple[nib.Nifti1Image, int]] = {}
    for rsn_name, maps in run_rsn_maps.items():
        if not maps:
            log.warning("  [%s | %s | %s] No SBA maps produced — skipping",
                        subject, rsn_name, run_type or "combined")
            continue
        avg_map, pooled_dof = average_maps(maps, run_rsn_dofs[rsn_name])
        out[rsn_name] = (avg_map, pooled_dof)

        map_path = rsn_sba_map_path(subject, rsn_name, run_type=run_type)
        nib.save(avg_map, str(map_path))

        thresh_map  = threshold_fc_map(avg_map, dof=pooled_dof, fdr_q=fdr_q,
                                       min_cluster=min_cluster, bilateral=bilateral,
                                       min_r_relative=min_r_relative)
        thresh_path = rsn_sba_map_path(subject, rsn_name, thresholded=True, run_type=run_type)
        nib.save(thresh_map, str(thresh_path))

        zstat_img, neglog10p_img = stat_maps(avg_map, dof=pooled_dof, bilateral=bilateral)
        zstat_path = rsn_zstat_path(subject, rsn_name, run_type=run_type)
        neglog10p_path = rsn_neglog10p_path(subject, rsn_name, run_type=run_type)
        nib.save(zstat_img, str(zstat_path))
        nib.save(neglog10p_img, str(neglog10p_path))

        provenance = dict(subject=subject, rsn=rsn_name, run_type=run_type or "combined",
                          runs_used=runs_used, dof=pooled_dof, fdr_q=fdr_q,
                          min_cluster=min_cluster, bilateral=bilateral,
                          min_r_relative=min_r_relative)
        for p in (map_path, thresh_path, zstat_path, neglog10p_path):
            write_provenance_json(p, **provenance)

    log.info("  [%s | %s] RSN maps saved (%d networks, dof pooled over %d run(s))",
             subject, run_type or "combined", len(out), len(runs_used))
    return out


def run_rsn_fc_subject(subject: str, force: bool = False,
                       fdr_q: float = 0.05, min_cluster: int = 50,
                       bilateral: bool = False, min_r_relative: float = 2.0,
                       profile_name: str | None = None,
                       seed_items: list[tuple[str, dict]] | None = None) -> Path | None:
    csv_path = fc_matrix_path(subject, "csv")
    if csv_path.exists() and not force:
        log.info("  [%s] Already exists — skipping", subject)
        return csv_path

    # Load subject-specific Yeo17 network masks (already in MNI 2mm, no resampling needed)
    rsn_masks = get_yeo17_subject_masks(subject)

    # Append language network from SynthSeg (not present in Yeo17)
    lang = get_language_network_mask(subject)
    if lang is not None:
        rsn_masks.append(lang)

    # Append auditory network from Smith 2009 RSN10 (intersected with subject GM)
    aud = get_auditory_network_mask(subject)
    if aud is not None:
        rsn_masks.append(aud)

    # Persist each RSN's mask alongside its SBA map so downstream consumers
    # (masked ICA) can discover "what is RSN <name>" purely from files this
    # step wrote, without importing this step's atlas-loading logic. Keeps
    # step3 atlas-agnostic and unable to drift out of sync with whichever
    # atlas step2 actually uses (this loop is cheap -- no BOLD involved).
    for rsn_name, rsn_mask_img in rsn_masks:
        mpath = rsn_mask_path(subject, rsn_name)
        if not mpath.exists() or force:
            nib.save(rsn_mask_img, str(mpath))

    rsn_names = [name for name, _ in rsn_masks]
    log.info("  [%s] RSN masks (%d): %s", subject, len(rsn_masks), rsn_names)
    seed_names = [n for n, _ in seed_items] if seed_items else []

    runs      = get_runs(subject)
    rest_runs = [r for r in runs if is_rest_run(r)]
    task_runs = [r for r in runs if not is_rest_run(r)]

    # rest/task runs are never blended into one number (see module docstring):
    # each run-type is pooled separately, both kept as distinct outputs, and
    # only combined afterward at the group level.
    group_data: dict[str, dict] = {}  # run_type -> {matrix, coupling, rsn_maps}

    for run_type, group_runs in (("rest", rest_runs), ("task", task_runs)):
        if not group_runs:
            continue

        run_matrices, run_rsn_maps, run_rsn_dofs, run_coupling_mats, runs_used = _compute_group(
            subject, group_runs, rsn_masks, rsn_names, seed_items
        )
        if not run_matrices:
            log.warning("  [%s | %s] No runs succeeded in this group — skipping", subject, run_type)
            continue

        group_mat = np.nanmean(np.stack(run_matrices, axis=0), axis=0)
        df = pd.DataFrame(group_mat, index=rsn_names, columns=rsn_names)
        df.to_csv(str(fc_matrix_path(subject, "csv", run_type=run_type)), float_format="%.4f")
        np.save(str(fc_matrix_path(subject, "npy", run_type=run_type)), group_mat)

        group_coup = None
        if run_coupling_mats and profile_name and seed_items:
            group_coup = np.nanmean(np.stack(run_coupling_mats, axis=0), axis=0)
            coup_df = pd.DataFrame(group_coup, index=seed_names, columns=rsn_names)
            coup_df.to_csv(str(coupling_matrix_path(subject, profile_name, run_type=run_type)),
                          float_format="%.4f")

        rsn_map_results = _save_group_rsn_maps(
            subject, run_type, run_rsn_maps, run_rsn_dofs, runs_used,
            fdr_q, min_cluster, bilateral, min_r_relative,
        )

        log.info("  [%s | %s] FC matrix saved (%d run(s) pooled)", subject, run_type, len(runs_used))
        group_data[run_type] = dict(matrix=group_mat, coupling=group_coup,
                                    rsn_maps=rsn_map_results, runs_used=runs_used)

    if not group_data:
        log.error("  [%s] No runs succeeded — no output written", subject)
        return None

    # Combine group-level results into the final "combined" (untagged) output.
    # Unweighted mean between the two groups, matching average_maps()'s own
    # point-estimate philosophy (dof is pooled for significance only, never
    # used to let one group outweigh the other in the estimate itself).
    if len(group_data) == 1:
        (run_type, only), = group_data.items()
        combined_mat  = only["matrix"]
        combined_coup = only["coupling"]
        combined_rsn_maps = only["rsn_maps"]
        combined_runs = only["runs_used"]
    else:
        rest_d, task_d = group_data["rest"], group_data["task"]
        combined_mat  = np.nanmean(np.stack([rest_d["matrix"], task_d["matrix"]]), axis=0)
        combined_coup = None
        if rest_d["coupling"] is not None and task_d["coupling"] is not None:
            combined_coup = np.nanmean(np.stack([rest_d["coupling"], task_d["coupling"]]), axis=0)
        combined_runs = rest_d["runs_used"] + task_d["runs_used"]

        combined_rsn_maps = {}
        rsn_names_present = set(rest_d["rsn_maps"]) | set(task_d["rsn_maps"])
        for rsn_name in rsn_names_present:
            in_rest = rsn_name in rest_d["rsn_maps"]
            in_task = rsn_name in task_d["rsn_maps"]
            if in_rest and in_task:
                m1, d1 = rest_d["rsn_maps"][rsn_name]
                m2, d2 = task_d["rsn_maps"][rsn_name]
                combined_rsn_maps[rsn_name] = average_maps([m1, m2], [d1, d2])
            elif in_rest:
                combined_rsn_maps[rsn_name] = rest_d["rsn_maps"][rsn_name]
            else:
                combined_rsn_maps[rsn_name] = task_d["rsn_maps"][rsn_name]

    df = pd.DataFrame(combined_mat, index=rsn_names, columns=rsn_names)
    df.to_csv(str(fc_matrix_path(subject, "csv")), float_format="%.4f")
    np.save(str(fc_matrix_path(subject, "npy")), combined_mat)
    log.info("  [%s] Combined FC matrix saved → %s", subject, fc_matrix_path(subject, "csv").name)

    if combined_coup is not None and profile_name and seed_items:
        coup_df = pd.DataFrame(combined_coup, index=seed_names, columns=rsn_names)
        coup_path = coupling_matrix_path(subject, profile_name)
        coup_df.to_csv(str(coup_path), float_format="%.4f")
        log.info("  [%s] Combined coupling matrix → %s", subject, coup_path.name)
    elif profile_name and seed_items:
        log.warning("  [%s] No coupling matrix produced", subject)

    # Save combined per-RSN maps: reconstruct (maps, dofs) lists per RSN from
    # whichever single map(s) contribute, then reuse the same save helper.
    combined_run_rsn_maps = {name: [m] for name, (m, _) in combined_rsn_maps.items()}
    combined_run_rsn_dofs = {name: [d] for name, (_, d) in combined_rsn_maps.items()}
    _save_group_rsn_maps(subject, None, combined_run_rsn_maps, combined_run_rsn_dofs,
                        combined_runs, fdr_q, min_cluster, bilateral, min_r_relative)

    return csv_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profile",      default=None, metavar="NAME",
                   help="Condition profile (from profiles.yaml); enables seed–RSN "
                        "coupling matrix computation alongside RSN FC maps")
    p.add_argument("--subjects",     nargs="+", default=None, metavar="ID")
    p.add_argument("--hv-only",      action="store_true")
    p.add_argument("--force",        action="store_true",
                   help="Recompute even if output already exists")
    p.add_argument("--fdr-q",        type=float, default=0.05, metavar="Q",
                   help="Benjamini-Hochberg FDR level, DOF-aware z-test on "
                        "the Fisher-Z map (default: 0.05)")
    p.add_argument("--min-cluster",  type=int,   default=50,   metavar="N",
                   help="Minimum cluster size in voxels (default: 50)")
    p.add_argument("--bilateral",    action="store_true",
                   help="Keep negative FC clusters (|r|≥thresh) in addition to positive")
    p.add_argument("--min-r-relative", type=float, default=2.0, metavar="X",
                   help="Effect-size floor as a multiple of each map's own "
                        "robust noise scale (default: 3.0x); 0 disables it")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.subjects:
        subjects = args.subjects
    elif args.hv_only:
        subjects = HV_SUBJECTS
    else:
        subjects = ALL_SUBJECTS

    # Load seed items if profile given (for coupling matrix)
    profile_name: str | None = None
    seed_items: list[tuple[str, dict]] | None = None
    if args.profile:
        profile_name = args.profile
        seed_names   = get_profile(profile_name)["seeds"]
        seed_items   = [(n, get_seed(n)) for n in seed_names]
        log.info("Profile   : %s | Seeds: %s", profile_name, seed_names)

    log.info("Subjects: %s", subjects)
    log.info("FDR q     : %.3f  |  min r relative: %.1fx  |  min cluster: %d  |  bilateral: %s",
             args.fdr_q, args.min_r_relative, args.min_cluster, args.bilateral)

    n_ok = n_fail = 0
    for subject in subjects:
        result = run_rsn_fc_subject(subject, force=args.force,
                                    fdr_q=args.fdr_q,
                                    min_cluster=args.min_cluster,
                                    bilateral=args.bilateral,
                                    min_r_relative=args.min_r_relative,
                                    profile_name=profile_name,
                                    seed_items=seed_items)
        if result:
            n_ok += 1
        else:
            n_fail += 1

    log.info("=" * 60)
    log.info("RSN FC complete: %d succeeded  |  %d failed", n_ok, n_fail)


if __name__ == "__main__":
    main()
