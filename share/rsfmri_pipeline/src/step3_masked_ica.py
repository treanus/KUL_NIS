#!/usr/bin/env python3
# Limit BLAS/OpenMP threads before numpy is imported to prevent CPU saturation.
import os
os.environ.setdefault("OMP_NUM_THREADS",     "4")
os.environ.setdefault("MKL_NUM_THREADS",     "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS","4")

"""
Step 3 — RSN-masked ICA, decomposed separately per run-type and recombined

Rest and tbfMRI-task runs are never concatenated together for a single
MELODIC decomposition: they routinely have different native TRs (a
resting-state protocol vs a motor/language localizer protocol), and even
when TR happens to match, task-block activation dominates a joint temporal
ICA and produces spurious "which run is this" components rather than
genuine networks (confirmed on real data: a component whose entire
timecourse was a two-level step function aligned with the task-run/rest-run
boundary).

For each subject × RSN that step 2 actually produced (discovered from the
sub-<ID>_rsn-<name>_desc-sba_statmap.nii.gz / _mask.nii.gz files step 2
writes -- not from any hardcoded atlas here, so this step automatically
follows whatever atlas step2 uses):
  1. Load the RSN mask and whole-brain FC map produced by step 2.
  2. Build the union mask: RSN mask ∪ FDR+effect-size-thresholded SBA map.
  3. Split runs into a "rest" group and a "task" group; within each group
     (only), temporally concatenate and run one MELODIC decomposition, using
     that group's own real per-run TR (read from the fmriprep BOLD sidecar,
     not a single hardcoded config value).
  4. In each group's decomposition, pick the single component with the best
     spatial correlation to the network's SBA map (within the union mask).
  5. If both groups produced a component, Stouffer-combine the two
     (already-in-z-units) maps, weighted by sqrt(dof) per group. If only one
     group is present (the common case for subjects with no tbfMRI runs),
     use it directly.
  6. Save the combined map, plus an FDR+cluster-extent thresholded version.

Runs are concatenated within subject and within run-type — never across
subjects, and never across rest/task.

Outputs (per subject, per RSN):
  analysis/masked_ica/sub-<ID>/rsn-<name>/
    sub-<ID>_rsn-<name>_union_mask.nii.gz
    sub-<ID>_rsn-<name>_rest_melodic/      FSL MELODIC output (if rest runs exist)
    sub-<ID>_rsn-<name>_task_melodic/      FSL MELODIC output (if task runs exist)
    sub-<ID>_rsn-<name>_ica_zstat.nii.gz          combined, unthresholded
    sub-<ID>_rsn-<name>_ica_zstatThresh.nii.gz    FDR + cluster-extent thresholded

Usage:
  python step3_masked_ica.py [--subjects HV01 PT01 ...] [--hv-only] [--dims 5] [--force]
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    ALL_SUBJECTS, HV_SUBJECTS, ICA_DIR,
    RSN_FC_DIR,
)
from utils import (
    build_union_mask, combine_zmaps, concatenate_bold_runs, get_run_tr,
    get_runs, is_rest_run, load_bold, load_melodic_components,
    pooled_dof_for_runs, resample_to, run_melodic, t1w_path,
    threshold_zstat_map, write_provenance_json,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def rsn_out_dir(subject: str, rsn_name: str) -> Path:
    d = ICA_DIR / f"sub-{subject}" / f"rsn-{rsn_name}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def rsn_sba_map_path(subject: str, rsn_name: str) -> Path:
    return (
        RSN_FC_DIR / f"sub-{subject}"
        / f"sub-{subject}_rsn-{rsn_name}_desc-sba_statmap.nii.gz"
    )


def rsn_mask_path(subject: str, rsn_name: str) -> Path:
    return RSN_FC_DIR / f"sub-{subject}" / f"sub-{subject}_rsn-{rsn_name}_mask.nii.gz"


def rsn_melodic_dir(subject: str, rsn_name: str, label: str) -> Path:
    return rsn_out_dir(subject, rsn_name) / f"sub-{subject}_rsn-{rsn_name}_{label}_melodic"


def discover_rsn_masks(subject: str) -> list[tuple[str, nib.Nifti1Image]]:
    """
    Find every RSN step 2 actually produced for this subject, by globbing its
    SBA statmap outputs -- never a hardcoded atlas/name list here, so step 3
    stays in lockstep with whatever atlas step 2 uses.
    """
    out_dir = RSN_FC_DIR / f"sub-{subject}"
    sba_paths = sorted(out_dir.glob(f"sub-{subject}_rsn-*_desc-sba_statmap.nii.gz"))

    masks = []
    for sba_path in sba_paths:
        rsn_name = sba_path.name.split("_rsn-")[1].split("_desc-")[0]
        mpath = rsn_mask_path(subject, rsn_name)
        if not mpath.exists():
            log.warning(
                "  [%s | rsn:%s] No mask file (%s) -- re-run step2 to regenerate it; skipping",
                subject, rsn_name, mpath.name,
            )
            continue
        masks.append((rsn_name, nib.load(str(mpath))))
    return masks


# ---------------------------------------------------------------------------
# Component selection
# ---------------------------------------------------------------------------

def select_best_component(comp_img: nib.Nifti1Image,
                          sba_img: nib.Nifti1Image,
                          mask_img: nib.Nifti1Image) -> tuple[np.ndarray, int, float]:
    """
    Correlate every IC's spatial map against this network's SBA statmap,
    restricted to voxels inside mask_img (already on comp_img's grid), and
    return the single best-|r|-matching component: (z_map, 0-based index, r).

    MELODIC component sign is arbitrary, so the returned z_map is flipped
    to positive-going if its best-matching r was negative -- otherwise a
    genuinely-matching but anti-correlated component would cancel out real
    signal when later combined with the other run-type's component.
    """
    sba_r = resample_to(sba_img, comp_img, interpolation="nearest").get_fdata(dtype=np.float32)
    mask  = mask_img.get_fdata(dtype=np.float32) > 0
    sba_vals = sba_r[mask]

    comp_data = comp_img.get_fdata(dtype=np.float32)
    best_idx, best_r, best_absr = 0, 0.0, -1.0
    for k in range(comp_data.shape[-1]):
        comp_vals = comp_data[..., k][mask]
        if comp_vals.std() == 0 or sba_vals.std() == 0:
            r = 0.0
        else:
            r = float(np.corrcoef(comp_vals, sba_vals)[0, 1])
        if abs(r) > best_absr:
            best_idx, best_r, best_absr = k, r, abs(r)

    z_map = comp_data[..., best_idx].copy()
    if best_r < 0:
        z_map = -z_map
    return z_map, best_idx, best_r


# ---------------------------------------------------------------------------
# Per-subject, per-RSN processing
# ---------------------------------------------------------------------------

def run_ica_subject_rsn(subject: str,
                         rsn_name: str,
                         rsn_mask_img: nib.Nifti1Image,
                         n_dims: int,
                         force: bool = False,
                         fdr_q: float = 0.05,
                         min_cluster: int = 15) -> bool:
    out_dir   = rsn_out_dir(subject, rsn_name)
    done_flag = out_dir / ".done"
    if done_flag.exists() and not force:
        log.info("  [%s | rsn:%s] Already complete — skipping", subject, rsn_name)
        return True

    # Clean up outputs from the old single-decomposition, dump-every-IC scheme
    for stale in out_dir.glob(f"sub-{subject}_rsn-{rsn_name}_IC-*_zstat.nii.gz"):
        stale.unlink()
    old_mel_dir = rsn_out_dir(subject, rsn_name) / f"sub-{subject}_rsn-{rsn_name}_melodic"
    if old_mel_dir.exists():
        shutil.rmtree(str(old_mel_dir))

    # Load RSN SBA map from step 2
    sba_path = rsn_sba_map_path(subject, rsn_name)
    if not sba_path.exists():
        log.error(
            "  [%s | rsn:%s] RSN SBA map not found: %s  (run step2 first)",
            subject, rsn_name, sba_path,
        )
        return False

    sba_img = nib.load(str(sba_path))
    t1_img  = nib.load(str(t1w_path(subject)))

    # Build union mask using first run as spatial reference
    runs     = get_runs(subject)
    ref_bold = load_bold(subject, runs[0])
    dof      = pooled_dof_for_runs(subject, runs)

    try:
        union_img = build_union_mask(rsn_mask_img, sba_img, ref_bold, dof=dof)
    except ValueError as exc:
        log.error("  [%s | rsn:%s] Union mask failed: %s", subject, rsn_name, exc)
        return False

    mask_save = out_dir / f"sub-{subject}_rsn-{rsn_name}_union_mask.nii.gz"
    nib.save(union_img, str(mask_save))
    log.info("  [%s | rsn:%s] Union mask: %d voxels → %s",
             subject, rsn_name, int(union_img.get_fdata().sum()), mask_save.name)

    # Decompose rest and task runs separately -- never concatenated together
    # (different native TRs, and task-block activation dominates a joint ICA).
    rest_runs = [r for r in runs if is_rest_run(r)]
    task_runs = [r for r in runs if not is_rest_run(r)]

    common_grid_img: nib.Nifti1Image | None = None
    group_results: list[tuple[np.ndarray, float, str, int, float]] = []

    for label, group_runs in (("rest", rest_runs), ("task", task_runs)):
        if not group_runs:
            continue

        mel_dir = rsn_melodic_dir(subject, rsn_name, label)

        if mel_dir.exists() and not force:
            log.info("  [%s | rsn:%s | %s] MELODIC output exists — reusing",
                     subject, rsn_name, label)
        else:
            if mel_dir.exists():
                shutil.rmtree(str(mel_dir))

            log.info("  [%s | rsn:%s | %s] Concatenating %d run(s)",
                     subject, rsn_name, label, len(group_runs))
            bold_concat = concatenate_bold_runs(subject, group_runs)
            if bold_concat is None:
                log.warning("  [%s | rsn:%s | %s] No BOLD runs could be loaded — skipping group",
                            subject, rsn_name, label)
                continue

            group_mask_img = resample_to(union_img, bold_concat, interpolation="nearest")
            group_bg_img   = resample_to(t1_img, bold_concat, interpolation="continuous")
            group_tr = get_run_tr(subject, group_runs[0])
            result = run_melodic(bold_concat, group_mask_img, mel_dir, n_dims, group_tr,
                                 bg_img=group_bg_img)
            del bold_concat
            if result is None:
                continue

        comp_img = load_melodic_components(mel_dir)
        if comp_img is None:
            continue

        comp_mask_img = resample_to(union_img, comp_img, interpolation="nearest")
        z_map, comp_idx, r = select_best_component(comp_img, sba_img, comp_mask_img)
        log.info("  [%s | rsn:%s | %s] Best IC %d  r=%.3f",
                 subject, rsn_name, label, comp_idx + 1, r)

        z_img_native = nib.Nifti1Image(z_map, comp_img.affine)
        if common_grid_img is None:
            common_grid_img = z_img_native
            z_common = z_map
        else:
            z_common = resample_to(z_img_native, common_grid_img,
                                   interpolation="nearest").get_fdata(dtype=np.float32)

        group_dof = pooled_dof_for_runs(subject, group_runs)
        group_results.append((z_common, group_dof, label, comp_idx, r))

    if not group_results:
        log.error("  [%s | rsn:%s] No group produced a usable ICA decomposition",
                  subject, rsn_name)
        return False

    if len(group_results) == 1:
        combined = group_results[0][0]
    else:
        (z1, dof1, *_), (z2, dof2, *_) = group_results
        combined = combine_zmaps(z1, dof1, z2, dof2)

    combined_img = nib.Nifti1Image(combined.astype(np.float32), common_grid_img.affine)
    raw_path = out_dir / f"sub-{subject}_rsn-{rsn_name}_ica_zstat.nii.gz"
    nib.save(combined_img, str(raw_path))

    # min_cluster default (15) is deliberately much smaller than steps 1-2's
    # whole-brain 50: the union mask here is only ~0.5-1.5% of whole-brain
    # voxel count, so real focal signal is naturally far smaller in absolute
    # voxel terms -- confirmed on real data where a genuine z=5.9 peak's
    # largest contiguous cluster was only 22 voxels, entirely killed by a
    # flat 50-voxel floor calibrated for full-brain-resolution maps.
    thresh_img = threshold_zstat_map(combined_img, fdr_q=fdr_q, min_cluster=min_cluster,
                                     bilateral=True)
    thresh_path = out_dir / f"sub-{subject}_rsn-{rsn_name}_ica_zstatThresh.nii.gz"
    nib.save(thresh_img, str(thresh_path))

    provenance = dict(
        subject=subject, rsn=rsn_name, run_type="+".join(g[2] for g in group_results),
        runs_used=runs, groups=[dict(run_type=g[2], best_ic=g[3] + 1, r=g[4], dof=g[1])
                                for g in group_results],
        fdr_q=fdr_q, min_cluster=min_cluster,
    )
    for p in (raw_path, thresh_path):
        write_provenance_json(p, **provenance)

    log.info("  [%s | rsn:%s] Combined ICA map saved (%s)",
             subject, rsn_name, " + ".join(g[2] for g in group_results))
    done_flag.touch()
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--subjects", nargs="+", default=None, metavar="ID")
    p.add_argument("--hv-only",  action="store_true")
    p.add_argument("--dims",     type=int, default=5, metavar="N",
                   help="ICA dimensionality (default: 5)")
    p.add_argument("--force",    action="store_true",
                   help="Re-run MELODIC even if output already exists")
    p.add_argument("--fdr-q",       type=float, default=0.05, metavar="Q",
                   help="Benjamini-Hochberg FDR threshold for the combined ICA map (default: 0.05)")
    p.add_argument("--min-cluster", type=int,   default=15,   metavar="N",
                   help="Minimum cluster size (voxels) to survive thresholding (default: 15 -- "
                        "the union mask is a small fraction of whole-brain, so real clusters "
                        "here are naturally much smaller than a whole-brain-calibrated floor)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.subjects:
        subjects = args.subjects
    elif args.hv_only:
        subjects = HV_SUBJECTS
    else:
        subjects = ALL_SUBJECTS

    log.info("Subjects : %s", subjects)
    log.info("ICA dims : %d", args.dims)

    n_ok = n_fail = 0

    for subject in subjects:
        get_runs(subject)  # fail fast with a clear error if no runs exist
        rsn_masks = discover_rsn_masks(subject)
        if not rsn_masks:
            log.error("  [%s] No step2 RSN outputs found -- run step2 first", subject)
            continue

        for rsn_name, rsn_mask_img in rsn_masks:
            log.info("[%s | rsn:%s]", subject, rsn_name)
            ok = run_ica_subject_rsn(
                subject, rsn_name, rsn_mask_img, args.dims, force=args.force,
                fdr_q=args.fdr_q, min_cluster=args.min_cluster,
            )
            n_ok, n_fail = (n_ok + 1, n_fail) if ok else (n_ok, n_fail + 1)

    log.info("=" * 60)
    log.info("Masked ICA complete: %d succeeded  |  %d failed", n_ok, n_fail)


if __name__ == "__main__":
    main()
