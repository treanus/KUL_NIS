#!/usr/bin/env python3
# Run with:  python3 step1_sba.py --profile <name> [options]
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
Step 1 — Seed-Based Analysis (SBA)

For each subject × seed defined in the chosen profile:
  1. Load the denoised BOLD and brain mask for every run.
  2. Resample the atlas seed mask to BOLD space.
  3. Extract the mean seed timeseries.
  4. Compute a whole-brain Pearson-r map → Fisher-Z transform, per run.
  5. Pool runs WITHIN each run-type (rest, task) separately -- never blended
     -- into a rest-only and a task-only subject-level map, each with its
     own dof-weighted average (see average_maps()). Mixing rest and tbfMRI
     task runs into one number was silently blending genuine resting-state
     connectivity with task-block co-activation, with the task group's
     larger raw timepoint count able to dominate simply by outnumbering
     rest, not because it's a better estimate of anything.
  6. If both run-types are present, combine the two GROUP-level maps (again
     via average_maps(), now weighted by each group's own pooled dof rather
     than by individual run count) into a "combined" map -- mirrors the
     rest/task recombination already used in step3's masked ICA.

Outputs (per subject, per seed; "" or "_runtype-rest" / "_runtype-task"):
  analysis/sba/sub-<ID>/sub-<ID>_seed-<name>[_runtype-X]_desc-sba_statmap.nii.gz
  ...same with desc-sbaThresh / desc-zstat / desc-neglog10p
  ...plus a <name>.json provenance sidecar next to every statmap above

Usage:
  python step1_sba.py --profile MDD [--subjects HV01 PT01 ...] [--save-runs]
"""

import argparse
import logging
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    ALL_SUBJECTS, HV_SUBJECTS, SBA_DIR,
    get_profile, get_seed, seed_subdir,
)
from utils import (
    average_maps, extract_mean_ts, get_runs, get_seed_mask, is_rest_run,
    load_bold, load_brain_mask, pooled_dof_for_runs, stat_maps,
    threshold_fc_map, whole_brain_fc, write_provenance_json,
)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Output path
# ---------------------------------------------------------------------------

def _tag(subject: str, seed_name: str, run_type: str | None) -> str:
    tag = f"sub-{subject}_seed-{seed_name}"
    if run_type:
        tag += f"_runtype-{run_type}"
    return tag


def sba_map_path(subject: str, seed_name: str, run_id: str | None = None,
                 thresholded: bool = False, run_type: str | None = None) -> Path:
    """
    Subject-level map (run_id=None) or per-run map (run_id='run-01' etc.).
    run_type=None is the combined (rest+task) map; "rest"/"task" are the
    per-run-type pooled maps. thresholded=True returns the cluster-filtered
    variant path.
    """
    out_dir = SBA_DIR / f"sub-{subject}" / seed_subdir(seed_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag  = _tag(subject, seed_name, run_type)
    desc = "sbaThresh" if thresholded else "sba"
    if run_id:
        return out_dir / f"{tag}_{run_id}_desc-{desc}_statmap.nii.gz"
    return out_dir / f"{tag}_desc-{desc}_statmap.nii.gz"


def sba_zstat_path(subject: str, seed_name: str, run_type: str | None = None) -> Path:
    out_dir = SBA_DIR / f"sub-{subject}" / seed_subdir(seed_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"{_tag(subject, seed_name, run_type)}_desc-zstat_statmap.nii.gz"


def sba_neglog10p_path(subject: str, seed_name: str, run_type: str | None = None) -> Path:
    out_dir = SBA_DIR / f"sub-{subject}" / seed_subdir(seed_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"{_tag(subject, seed_name, run_type)}_desc-neglog10p_statmap.nii.gz"


# ---------------------------------------------------------------------------
# Per-subject, per-seed processing
# ---------------------------------------------------------------------------

def _save_group_outputs(subject: str, seed_name: str, run_type: str | None,
                        group_map: nib.Nifti1Image, group_dof: float,
                        runs_used: list[str],
                        fdr_q: float, min_cluster: int, bilateral: bool,
                        min_r_relative: float) -> None:
    """Save raw + thresholded + zstat + -log10(p) + provenance for one map
    (either a per-run-type group map, or the final combined map)."""
    raw_path = sba_map_path(subject, seed_name, run_type=run_type)
    nib.save(group_map, str(raw_path))

    thresh_map = threshold_fc_map(group_map, dof=group_dof, fdr_q=fdr_q,
                                  min_cluster=min_cluster, bilateral=bilateral,
                                  min_r_relative=min_r_relative)
    thresh_path = sba_map_path(subject, seed_name, thresholded=True, run_type=run_type)
    nib.save(thresh_map, str(thresh_path))

    zstat_img, neglog10p_img = stat_maps(group_map, dof=group_dof, bilateral=bilateral)
    zstat_path = sba_zstat_path(subject, seed_name, run_type=run_type)
    neglog10p_path = sba_neglog10p_path(subject, seed_name, run_type=run_type)
    nib.save(zstat_img, str(zstat_path))
    nib.save(neglog10p_img, str(neglog10p_path))

    provenance = dict(
        subject=subject, seed=seed_name, run_type=run_type or "combined",
        runs_used=runs_used, dof=group_dof, fdr_q=fdr_q,
        min_cluster=min_cluster, bilateral=bilateral,
        min_r_relative=min_r_relative,
    )
    for p in (raw_path, thresh_path, zstat_path, neglog10p_path):
        write_provenance_json(p, **provenance)

    log.info("  [%s | %s | %s] Saved raw + thresholded + zstat + -log10(p) (dof=%d, %d run(s))",
             subject, seed_name, run_type or "combined", group_dof, len(runs_used))


def run_sba_subject_seed(subject: str,
                         seed_name: str,
                         seed_cfg: dict,
                         save_runs: bool = False,
                         fdr_q: float = 0.05,
                         min_cluster: int = 50,
                         bilateral: bool = False,
                         min_r_relative: float = 2.0,
                         force: bool = False) -> Path | None:
    """
    Compute subject-level SBA map for one seed: rest-only and task-only
    pooled maps kept as separate outputs, plus a combined map when both
    run-types are present. See module docstring.
    Returns path to the saved combined map, or None if processing failed.
    """
    out_path = sba_map_path(subject, seed_name)
    if out_path.exists() and not force:
        log.info("  [%s | %s] Already exists — skipping", subject, seed_name)
        return out_path

    try:
        seed_img = get_seed_mask(seed_name, seed_cfg, subject=subject)
    except Exception as exc:
        log.error("  [%s | %s] Failed to build seed mask: %s", subject, seed_name, exc)
        return None

    runs      = get_runs(subject)
    rest_runs = [r for r in runs if is_rest_run(r)]
    task_runs = [r for r in runs if not is_rest_run(r)]

    group_results: list[tuple[str, nib.Nifti1Image, float, list[str]]] = []

    for run_type, group_runs in (("rest", rest_runs), ("task", task_runs)):
        if not group_runs:
            continue

        run_maps, run_dofs, runs_used = [], [], []
        for run_id in group_runs:
            log.info("  [%s | %s | %s] Computing SBA", subject, seed_name, run_id)
            try:
                bold_img  = load_bold(subject, run_id)
                mask_img  = load_brain_mask(subject, run_id)
                bold_data = bold_img.get_fdata(dtype=np.float32)
                seed_ts   = extract_mean_ts(bold_data, seed_img, bold_img)
                fc_map    = whole_brain_fc(bold_data, seed_ts, mask_img, bold_img)
                run_dof   = bold_data.shape[-1] - 3
                del bold_data, bold_img, mask_img

                if save_runs:
                    run_path = sba_map_path(subject, seed_name, run_id, run_type=run_type)
                    nib.save(fc_map, str(run_path))
                    log.info("    Saved per-run map → %s", run_path.name)

                run_maps.append(fc_map)
                run_dofs.append(run_dof)
                runs_used.append(run_id)

            except Exception as exc:
                log.warning(
                    "  [%s | %s | %s] Run failed (%s) — skipping this run",
                    subject, seed_name, run_id, exc,
                )

        if not run_maps:
            log.warning("  [%s | %s | %s] No runs succeeded in this group — skipping",
                        subject, seed_name, run_type)
            continue

        group_map, group_dof = average_maps(run_maps, run_dofs)
        _save_group_outputs(subject, seed_name, run_type, group_map, group_dof, runs_used,
                            fdr_q, min_cluster, bilateral, min_r_relative)
        group_results.append((run_type, group_map, group_dof, runs_used))

    if not group_results:
        log.error("  [%s | %s] No runs succeeded — no output written", subject, seed_name)
        return None

    if len(group_results) == 1:
        _, combined_map, combined_dof, combined_runs = group_results[0]
    else:
        (_, m1, d1, r1), (_, m2, d2, r2) = group_results
        combined_map, combined_dof = average_maps([m1, m2], [d1, d2])
        combined_runs = r1 + r2

    _save_group_outputs(subject, seed_name, None, combined_map, combined_dof, combined_runs,
                        fdr_q, min_cluster, bilateral, min_r_relative)

    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profile",   required=True,
                   help="Condition profile name (defined in profiles.yaml)")
    p.add_argument("--subjects",  nargs="+", default=None,
                   metavar="ID",
                   help="Subject IDs to process (default: all subjects in config)")
    p.add_argument("--hv-only",   action="store_true",
                   help="Process healthy volunteers only")
    p.add_argument("--save-runs",    action="store_true",
                   help="Also save per-run SBA maps (in addition to averaged map)")
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

    profile = get_profile(args.profile)
    seed_names = profile["seeds"]

    if args.subjects:
        subjects = args.subjects
    elif args.hv_only:
        subjects = HV_SUBJECTS
    else:
        subjects = ALL_SUBJECTS

    log.info("Profile   : %s", args.profile)
    log.info("Seeds     : %s", seed_names)
    log.info("Subjects  : %s", subjects)
    log.info("Save runs : %s", args.save_runs)
    log.info("FDR q     : %.3f  |  min r relative: %.1fx  |  min cluster: %d  |  bilateral: %s",
             args.fdr_q, args.min_r_relative, args.min_cluster, args.bilateral)

    n_ok = n_fail = 0

    for subject in subjects:
        for seed_name in seed_names:
            seed_cfg = get_seed(seed_name)

            result = run_sba_subject_seed(
                subject, seed_name, seed_cfg,
                save_runs=args.save_runs,
                fdr_q=args.fdr_q,
                min_cluster=args.min_cluster,
                bilateral=args.bilateral,
                min_r_relative=args.min_r_relative,
                force=args.force,
            )
            if result:
                n_ok += 1
            else:
                n_fail += 1

    log.info("=" * 60)
    log.info("SBA complete: %d succeeded  |  %d failed", n_ok, n_fail)


if __name__ == "__main__":
    main()
