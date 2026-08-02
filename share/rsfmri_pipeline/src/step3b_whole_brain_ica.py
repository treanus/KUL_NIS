#!/usr/bin/env python3
# Limit BLAS/OpenMP threads before numpy is imported to prevent CPU saturation.
import os
os.environ.setdefault("OMP_NUM_THREADS",     "4")
os.environ.setdefault("MKL_NUM_THREADS",     "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS","4")

"""
Step 3b — Whole-brain, unconstrained ICA

Unlike step3 (masked ICA, one small union-mask search space per RSN, single
best-matching component picked out), this runs one MELODIC decomposition
over the WHOLE brain mask at a higher dimensionality, purely exploratory --
there is no target network to match components against, so every component
is kept and shown (like a standard group-ICA report), not reduced to a
single "best" map.

Same rest/task-separation and per-run z-scoring as step3, for the same
reasons (different native TRs; task-block activation dominating a joint
decomposition otherwise) -- see step3_masked_ica.py's module docstring for
the full rationale and the real-data numbers behind it. Rest and task are
decomposed independently and NEVER recombined here: there is no natural
component-to-component correspondence between two independent whole-brain
decompositions the way there is a specific network to match against in
step3, so recombination would just be pairing up arbitrary components.

Outputs (per subject):
  analysis/whole_brain_ica/sub-<ID>/
    sub-<ID>_wbica_rest_melodic/         FSL MELODIC output (if rest runs exist)
    sub-<ID>_wbica_task_melodic/         FSL MELODIC output (if task runs exist)
    sub-<ID>_wbica_rest_IC-<k>_zstat.nii.gz   (× dims)
    sub-<ID>_wbica_task_IC-<k>_zstat.nii.gz   (× dims)
    sub-<ID>_wbica_<rest|task>.json      provenance (dof, dims, runs used, TR)

Usage:
  python step3b_whole_brain_ica.py [--subjects HV01 PT01 ...] [--hv-only] [--dims 25] [--force]
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from config import ALL_SUBJECTS, HV_SUBJECTS, WBICA_DIR
from utils import (
    concatenate_bold_runs, get_run_tr, get_runs, is_rest_run, load_bold,
    load_brain_mask, load_melodic_components, pooled_dof_for_runs,
    resample_to, run_melodic, t1w_path, write_provenance_json,
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

def subject_out_dir(subject: str) -> Path:
    d = WBICA_DIR / f"sub-{subject}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def melodic_dir(subject: str, run_type: str) -> Path:
    return subject_out_dir(subject) / f"sub-{subject}_wbica_{run_type}_melodic"


def ic_path(subject: str, run_type: str, k: int) -> Path:
    return subject_out_dir(subject) / f"sub-{subject}_wbica_{run_type}_IC-{k:02d}_zstat.nii.gz"


def done_flag(subject: str, run_type: str) -> Path:
    return subject_out_dir(subject) / f".done_{run_type}"


# ---------------------------------------------------------------------------
# Per-subject, per-run-type processing
# ---------------------------------------------------------------------------

def run_wbica_subject_group(subject: str, run_type: str, group_runs: list[str],
                            n_dims: int, force: bool = False) -> bool:
    if not group_runs:
        return False

    flag = done_flag(subject, run_type)
    if flag.exists() and not force:
        log.info("  [%s | %s] Already complete — skipping", subject, run_type)
        return True

    out_dir = subject_out_dir(subject)
    for stale in out_dir.glob(f"sub-{subject}_wbica_{run_type}_IC-*_zstat.nii.gz"):
        stale.unlink()

    mel_dir = melodic_dir(subject, run_type)
    if mel_dir.exists():
        shutil.rmtree(str(mel_dir))

    log.info("  [%s | %s] Concatenating %d run(s)", subject, run_type, len(group_runs))
    bold_concat = concatenate_bold_runs(subject, group_runs)
    if bold_concat is None:
        log.error("  [%s | %s] No BOLD runs could be loaded", subject, run_type)
        return False

    mask_img = resample_to(load_brain_mask(subject, group_runs[0]), bold_concat,
                           interpolation="nearest")
    bg_img   = resample_to(nib.load(str(t1w_path(subject))), bold_concat,
                           interpolation="continuous")
    tr = get_run_tr(subject, group_runs[0])

    result = run_melodic(bold_concat, mask_img, mel_dir, n_dims, tr, bg_img=bg_img)
    del bold_concat
    if result is None:
        return False

    comp_img = load_melodic_components(mel_dir)
    if comp_img is None:
        return False

    for k in range(comp_img.shape[-1]):
        comp_vol = comp_img.get_fdata(dtype=np.float32)[..., k]
        comp_nii = nib.Nifti1Image(comp_vol, comp_img.affine)
        path = ic_path(subject, run_type, k + 1)
        nib.save(comp_nii, str(path))
        write_provenance_json(
            path, subject=subject, run_type=run_type, ic=k + 1,
            n_dims=comp_img.shape[-1], tr=tr, runs_used=group_runs,
            dof=pooled_dof_for_runs(subject, group_runs),
        )

    log.info("  [%s | %s] Saved %d IC maps", subject, run_type, comp_img.shape[-1])
    flag.touch()
    return True


def run_wbica_subject(subject: str, n_dims: int, force: bool = False) -> bool:
    runs      = get_runs(subject)
    rest_runs = [r for r in runs if is_rest_run(r)]
    task_runs = [r for r in runs if not is_rest_run(r)]

    ok_any = False
    for run_type, group_runs in (("rest", rest_runs), ("task", task_runs)):
        if not group_runs:
            continue
        ok = run_wbica_subject_group(subject, run_type, group_runs, n_dims, force=force)
        ok_any = ok_any or ok

    return ok_any


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--subjects", nargs="+", default=None, metavar="ID")
    p.add_argument("--hv-only",  action="store_true")
    p.add_argument("--dims",     type=int, default=25, metavar="N",
                   help="ICA dimensionality (default: 25 -- whole-brain exploratory "
                        "decomposition, higher than masked-ICA's per-network default "
                        "since it isn't targeting one known network)")
    p.add_argument("--force",    action="store_true",
                   help="Re-run MELODIC even if output already exists")
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
        log.info("[%s]", subject)
        ok = run_wbica_subject(subject, args.dims, force=args.force)
        n_ok, n_fail = (n_ok + 1, n_fail) if ok else (n_ok, n_fail + 1)

    log.info("=" * 60)
    log.info("Whole-brain ICA complete: %d succeeded  |  %d failed", n_ok, n_fail)


if __name__ == "__main__":
    main()
