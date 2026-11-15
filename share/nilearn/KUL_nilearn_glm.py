#!/usr/bin/env python3
"""
KUL_nilearn_glm.py

MATLAB/SPM-free first-level task-fMRI GLM using Nilearn, designed as a drop-in
replacement for the SPM12 step in KUL_fmriproc_*.sh.

Assumes a block design: <rest-dur> s REST followed by <task-dur> s TASK epochs,
repeated for the length of the run. Models a single "task" regressor; the
implicit baseline is REST, so the task>rest contrast == spmT_0001.

Inputs are fmriprep MNI-space, already brain-masked and SUSAN-smoothed by the
bash driver. Supports one or more runs (fixed-effects combination across runs).

Outputs (written to --output-dir, SPM-compatible names so the bash warp-back
step works unchanged):
    spmT_0001.nii              unthresholded t map (task>rest)
    spmZ_0001.nii              unthresholded z map
    spmT_0001_p001unc_k50.nii  voxelwise p<0.001 uncorrected, cluster k>=--cluster-k
    spmT_0001_FWE<tag>_k50.nii  voxelwise FWE (Bonferroni) p<--pfwe, cluster k>=--cluster-k
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import nibabel as nib

from nilearn.glm.first_level import FirstLevelModel
from nilearn.glm import threshold_stats_img


HIGH_PASS_HZ = 1.0 / 128.0   # matches SPM's default 128 s high-pass filter


def build_events(n_scans, t_r, rest_dur, task_dur, first):
    """Block design: alternating rest/task; only 'task' is modelled."""
    run_dur = n_scans * t_r
    cycle = rest_dur + task_dur
    start = 0.0 if first == "task" else float(rest_dur)
    onsets = np.arange(start, run_dur, cycle)
    onsets = onsets[onsets < run_dur]
    durations = np.minimum(task_dur, run_dur - onsets)
    return pd.DataFrame(
        {"onset": onsets, "duration": durations, "trial_type": "ON"}
    )


def load_confounds(path):
    """Load the (header-less) filtered confounds .txt produced by KUL_tsv_filter."""
    df = pd.read_csv(path, sep="\t", header=None)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df.columns = [f"c{i}" for i in range(df.shape[1])]
    return df


def n_scans_of(img_path):
    img = nib.load(img_path)
    return img.shape[3] if img.ndim == 4 else 1


def pfwe_to_tag(pfwe):
    """Mirror the bash tag logic: strip '0.' then trailing zeros. 0.01 -> '01'."""
    s = ("%g" % pfwe)
    if s.startswith("0."):
        s = s[2:]
    s = s.rstrip("0")
    return s if s else "0"


def main():
    ap = argparse.ArgumentParser(description="Nilearn first-level task GLM")
    ap.add_argument("--bold", nargs="+", required=True,
                    help="one or more (masked, smoothed) MNI BOLD runs")
    ap.add_argument("--confounds", nargs="*", default=None,
                    help="filtered confounds .txt per run (order matches --bold)")
    ap.add_argument("--mask", required=True, help="brain mask (fmriprep, MNI)")
    ap.add_argument("--tr", type=float, required=True)
    ap.add_argument("--task-name", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--rest-dur", type=float, default=30.0)
    ap.add_argument("--task-dur", type=float, default=30.0)
    ap.add_argument("--first", choices=["rest", "task"], default="rest")
    ap.add_argument("--hrf-derivative", type=int, default=0,
                    help="0=canonical HRF, 1=HRF + temporal derivative")
    ap.add_argument("--pfwe", type=float, default=0.01)
    ap.add_argument("--cluster-k", type=int, default=50)
    ap.add_argument("--n-jobs", type=int, default=1,
                    help="Nilearn joblib workers for the voxelwise GLM fit "
                         "(parallelizes a single run's GLM; -1 = all cores)")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    bolds = args.bold
    confounds_paths = args.confounds if args.confounds else None
    if confounds_paths is not None and len(confounds_paths) != len(bolds):
        sys.exit("ERROR: number of --confounds (%d) != number of --bold (%d)"
                 % (len(confounds_paths), len(bolds)))

    hrf_model = "spm + derivative" if args.hrf_derivative else "spm"

    # per-run events (and confounds)
    events_list = [build_events(n_scans_of(b), args.tr,
                                args.rest_dur, args.task_dur, args.first)
                   for b in bolds]
    confounds_list = ([load_confounds(c) for c in confounds_paths]
                      if confounds_paths is not None else None)

    print("  [nilearn] runs=%d  TR=%.4f  hrf='%s'  confounds=%s  n_jobs=%d"
          % (len(bolds), args.tr, hrf_model, confounds_paths is not None, args.n_jobs))

    # Parameters matched to the SPM12 job files:
    #   timing.RT=TR, hpf=128 s, cvi='AR(1)', bases.hrf.derivs -> hrf_model,
    #   fmri_t=16 / fmri_t0=8  ->  slice_time_ref = 8/16 = 0.5,
    #   global='None' + grand-mean scaling (t-stats are invariant to the scale factor).
    # Intentional difference: SPM used implicit masking (mthresh=0.8, mask='');
    # here the fmriprep brain mask is passed explicitly (the whole point of this port).
    model = FirstLevelModel(
        t_r=args.tr,
        slice_time_ref=0.5,          # SPM fmri_t0=8 of fmri_t=16
        hrf_model=hrf_model,
        drift_model="cosine",
        high_pass=HIGH_PASS_HZ,
        noise_model="ar1",           # SPM cvi='AR(1)'
        mask_img=args.mask,          # <-- the fmriprep brain mask enters the GLM here
        smoothing_fwhm=None,         # already SUSAN-smoothed upstream
        signal_scaling=0,            # grand-mean scaling (SPM-like); t-map unaffected
        n_jobs=args.n_jobs,          # parallelize the voxelwise fit within this GLM
        minimize_memory=False,
    )
    model.fit(bolds if len(bolds) > 1 else bolds[0],
              events=events_list if len(bolds) > 1 else events_list[0],
              confounds=confounds_list if confounds_list is None or len(bolds) > 1
              else confounds_list[0])

    tmap = model.compute_contrast("ON", output_type="stat")     # t (== spmT_0001)
    zmap = model.compute_contrast("ON", output_type="z_score")  # z (for thresholding)

    out = args.output_dir
    nib.save(tmap, os.path.join(out, "spmT_0001.nii"))
    nib.save(zmap, os.path.join(out, "spmZ_0001.nii"))
    print("  [nilearn] wrote spmT_0001.nii / spmZ_0001.nii")

    fwe_tag = "FWE%s_k%d" % (pfwe_to_tag(args.pfwe), args.cluster_k)
    unc_tag = "p001unc_k%d" % args.cluster_k

    tdata = tmap.get_fdata()

    def save_t_in_clusters(thr_img, tag):
        """Write SPM-style output: t-values inside surviving clusters, 0 elsewhere."""
        cluster_mask = (np.asanyarray(thr_img.get_fdata()) != 0).astype(tdata.dtype)
        out_img = nib.Nifti1Image(tdata * cluster_mask, tmap.affine, tmap.header)
        nib.save(out_img, os.path.join(out, "spmT_0001_%s.nii" % tag))
        return int(np.count_nonzero(cluster_mask))

    # Primary: voxelwise p<0.001 uncorrected, cluster k>=cluster_k, positive only.
    # (voxel set == SPM's tinv(1-0.001, df) height threshold; we store t-values.)
    try:
        thr_unc, z_unc = threshold_stats_img(
            zmap, alpha=0.001, height_control="fpr",
            cluster_threshold=args.cluster_k, two_sided=False)
        nvox = save_t_in_clusters(thr_unc, unc_tag)
        print("  [nilearn] p<0.001 unc (z>=%.3f) k>=%d -> %s (%d voxels)"
              % (z_unc, args.cluster_k, unc_tag, nvox))
    except Exception as e:
        print("  [nilearn] WARNING uncorrected thresholding failed: %s" % e)

    # Secondary: FWE via Bonferroni. NB SPM used RFT (spm_uc); see notes. Bonferroni
    # is a valid but more conservative FWE control than SPM's smoothness-based RFT.
    try:
        thr_fwe, z_fwe = threshold_stats_img(
            zmap, alpha=args.pfwe, height_control="bonferroni",
            cluster_threshold=args.cluster_k, two_sided=False)
        nvox = save_t_in_clusters(thr_fwe, fwe_tag)
        print("  [nilearn] FWE(Bonferroni) p<%g (z>=%.3f) k>=%d -> %s (%d voxels)"
              % (args.pfwe, z_fwe, args.cluster_k, fwe_tag, nvox))
    except Exception as e:
        print("  [nilearn] WARNING FWE thresholding failed: %s" % e)


if __name__ == "__main__":
    main()
