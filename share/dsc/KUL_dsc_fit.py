#!/usr/bin/env python3
"""DSC-MRI quantification: dR2*, BSW-style leakage correction and SVD deconvolution.

Reads a preprocessed 4D DSC series plus a brain mask and writes the parametric
maps (rCBV corrected/uncorrected, K1, K2, rCBF, MTT, TTP, TT0) that
KUL_dsc_perfusion.sh registers into the participant's T1w space.

Derived from tools/KUL_DSC_analysis/Good_DSCLC_fit5.py (A. Radwan). Relative to
that prototype:
  * the removed scipy.integrate.cumtrapz/simps are replaced by
    cumulative_trapezoid/simpson (the originals are gone in scipy >= 1.14),
  * the per-voxel leakage-fit and deconvolution loops are vectorised, which is
    what makes this usable at clinical turnaround (minutes, not hours),
  * deconvolution defaults to a dt-scaled AIF with truncated-SVD
    regularisation, so MTT comes out in seconds. --legacy reproduces the
    prototype's z-scored AIF and 1/(S+1e-3) damping exactly,
  * the pre-contrast baseline window is detected from the data instead of being
    fixed at frames 5-10, and a bolus that arrives before frame 10 no longer
    silently inverts the sign of every map.

@ Ahmed Radwan - KU Leuven, Translational MRI - ahmed.radwan@kuleuven.be
"""

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")  # no DISPLAY on the compute nodes

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy.linalg import svd
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.signal import savgol_filter

try:  # scipy >= 1.6
    from scipy.integrate import cumulative_trapezoid, simpson
except ImportError:  # pragma: no cover - very old scipy
    from scipy.integrate import cumtrapz as cumulative_trapezoid
    from scipy.integrate import simps as simpson

try:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "scikit-learn is required for AIF detection. Install it in the conda env "
        "used for this step (the KUL_NIS 'pyfMRI' env gets it via nilearn)."
    ) from exc


# --- I/O -----------------------------------------------------------------


def load_nifti(path):
    img = nib.load(path)
    return img.get_fdata(dtype=np.float32), img.affine


def save_nifti(data, affine, path):
    """Save 'data' as float32 NIfTI; NaNs (out-of-mask, with --apply_mask) survive."""
    nib.save(nib.Nifti1Image(np.asarray(data, dtype=np.float32), affine), path)


# --- signal -> dR2* ------------------------------------------------------


def compute_delta_r2star(signal, TE, baseline_start=5, baseline_end=10, eps=1e-3):
    """dR2*(t) = -(1/TE) * ln(S(t) / S0), with S0 the mean over the baseline frames.

    Both S0 and S are clipped at 'eps' so the log argument stays positive; any
    residual non-finite value becomes 0.
    """
    baseline = np.mean(signal[..., baseline_start:baseline_end], axis=-1, keepdims=True)
    baseline = np.clip(baseline, eps, None)
    signal = np.clip(signal, eps, None)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        delta = -np.log(signal / baseline) / TE
        delta[~np.isfinite(delta)] = 0.0

    return delta


def compute_r2star_no_negatives(signal, TE, baseline_start=5, baseline_end=10, eps=1e-3):
    """R2*(t) = (1/TE) * ln(S0 / S(t)), floored at 0 where S(t) > S0."""
    baseline = np.mean(signal[..., baseline_start:baseline_end], axis=-1, keepdims=True)
    baseline = np.clip(baseline, eps, None)
    signal = np.clip(signal, eps, None)

    with np.errstate(divide="ignore", invalid="ignore"):
        r2star = np.log(baseline / signal) / TE
        r2star[r2star < 0] = 0
        r2star[~np.isfinite(r2star)] = 0.0

    return r2star


def detect_baseline_window(dsc, mask, min_frames=3, fallback=(5, 10), frac=0.05, guard=2, edge=3):
    """Find the pre-contrast frames from the mask-average signal time course.

    The bolus peak is the deepest point of the mean signal, searched away from
    the first and last 'edge' frames. The baseline is then the run of frames
    immediately preceding it that are still within 'frac' of the plateau level,
    ending 'guard' frames before the signal starts to fall.

    The prototype hardcoded frames 5-10. That is fine when the bolus arrives
    late, but if it arrives earlier the "baseline" is measured part-way down the
    bolus, S0 comes out too low, dR2* goes negative and every downstream map --
    MTT in particular -- silently changes sign.

    Two details keep this honest on real series. The peak search ignores the
    edge frames, because non-steady-state at the start and filter roll-off at
    the end can both dip below the true bolus. And the baseline is taken as the
    *contiguous* plateau run before the bolus, so a depressed leading frame is
    excluded rather than dragging S0 down. Falls back to the fixed window when
    no onset is detectable.
    """
    n_frames = dsc.shape[-1]
    curve = dsc[mask > 0].reshape(-1, n_frames).mean(axis=0)

    lo = min(edge, max(0, n_frames // 4))
    hi = max(lo + 1, n_frames - lo)
    peak_idx = lo + int(np.argmin(curve[lo:hi]))

    if peak_idx > min_frames:
        plateau = curve[:peak_idx].max()
        bottom = curve[peak_idx]
        if plateau > bottom:
            threshold = plateau - frac * (plateau - bottom)
            pre = np.where(curve[:peak_idx] >= threshold)[0]
            if pre.size:
                # walk back over the contiguous run ending at pre[-1]
                run_start = int(pre[-1])
                pre_set = set(int(i) for i in pre)
                while run_start - 1 in pre_set:
                    run_start -= 1
                start = max(run_start, 1)  # frame 0 is often not at steady state
                end = int(pre[-1]) + 1 - guard
                if end - start >= min_frames:
                    return start, end

    start, end = fallback
    return start, min(end, n_frames)


# --- AIF -----------------------------------------------------------------


def detect_aif(dsc, mask, num_candidates=200, mode="single", baseline=(5, 10)):
    """Pick an arterial input function from the voxels with the largest signal drop.

    The candidate pool is the 'num_candidates' voxels with the biggest
    baseline-minus-minimum drop. PCA + KMeans then locates the most typical
    curve in that pool; mode='single' returns it, mode='mean' returns the pool
    average (less noisy, but broadens the bolus).

    The prototype additionally required an absolute drop > 100 intensity units
    before ranking. That threshold is scanner- and scaling-dependent: where more
    than num_candidates voxels cleared it the ranking was identical anyway, and
    where none did it left an empty pool and crashed PCA. It is gone.
    """
    ts = dsc[mask > 0].reshape(-1, dsc.shape[-1])
    if ts.shape[0] == 0:
        raise SystemExit("Brain mask is empty - cannot detect an AIF.")

    pre_contrast = np.mean(ts[:, baseline[0] : baseline[1]], axis=1)
    drop = pre_contrast - np.min(ts, axis=1)

    n_take = int(min(num_candidates, ts.shape[0]))
    top = np.argsort(drop)[-n_take:]
    candidates = ts[top]

    if mode == "mean":
        return candidates.mean(axis=0), n_take

    reduced = PCA(n_components=min(2, candidates.shape[0])).fit_transform(candidates)
    kmeans = KMeans(n_clusters=1, n_init=10, random_state=42).fit(reduced)
    return candidates[np.argmin(kmeans.transform(reduced))], n_take


def smooth_signal(sig):
    """Savitzky-Golay smoothing of the 1D AIF (window 7, order 2)."""
    window = 7 if len(sig) >= 7 else (len(sig) // 2) * 2 + 1
    if window < 3:
        return sig.astype(np.float64)
    return savgol_filter(sig.astype(np.float64), window, 2)


def enforce_flat_baseline_tail(sig, N=5):
    """Flatten the first/last N points, so deconvolution sees a stable baseline and tail."""
    sig = np.array(sig, dtype=np.float64, copy=True)
    if len(sig) > 2 * N + 3:
        sig[:N] = np.mean(sig[N : N + 3])
        sig[-N:] = np.mean(sig[-N - 3 : -N])
    return sig


def save_aif_plots(raw, smoothed, deltaR2, output_dir):
    for data, title, ylabel, fname in (
        (raw, "Raw AIF", "Signal intensity", "aif_raw.png"),
        (smoothed, "Smoothed AIF", "Signal intensity", "aif_smoothed.png"),
        (deltaR2, "dR2* AIF", "dR2* (1/s)", "aif_deltaR2star.png"),
    ):
        plt.figure()
        plt.plot(data)
        plt.title(title)
        plt.xlabel("Timepoints")
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, fname))
        plt.close()


# --- leakage correction --------------------------------------------------


def bsw_leakage_correction(tissue, aif_deltaR2, time):
    """Vectorised leakage correction over an [N_vox, T] array of dR2* curves.

    Each voxel curve is regressed on the cumulative AIF, and the fitted
    component is subtracted before integrating:

        dR2*_tissue(t) ~ K2 * INT(AIF) + K1
        rCBV_corr      = INT[ dR2*_tissue(t) - K2 * INT(AIF) ]

    Same model and same K1/K2 convention as the prototype's per-voxel
    scipy.stats.linregress loop, computed in closed form: because the regressor
    INT(AIF) is shared by every voxel, the slope reduces to
    cov(x, y) / var(x) and the whole volume becomes one matrix-vector product.

    NOTE: this regresses on the integral of the *arterial* input function.
    Textbook Boxerman-Schmainda-Weisskoff instead regresses tissue on a
    whole-brain *non-enhancing reference tissue* curve and its integral. The
    prototype's formulation is kept here deliberately - see
    docs/KUL_dsc_perfusion/KUL_dsc_perfusion.md, "Known deviations".
    """
    aif_int = cumulative_trapezoid(aif_deltaR2, time, initial=0)

    x = aif_int - aif_int.mean()
    denom = float(x @ x)
    if denom <= 0:
        raise SystemExit(
            "The cumulative AIF is constant - leakage correction is undefined. "
            "Check the AIF diagnostic plots in the output directory."
        )

    y_mean = tissue.mean(axis=1)
    K2 = (tissue - y_mean[:, None]) @ x / denom
    K1 = y_mean - K2 * aif_int.mean()

    corrected = tissue - K2[:, None] * aif_int[None, :]

    rCBV_unc = simpson(tissue, x=time, axis=-1)
    rCBV_corr = simpson(corrected, x=time, axis=-1)

    return rCBV_unc, rCBV_corr, K2, K1


# --- deconvolution -------------------------------------------------------


def build_deconvolution_operator(aif_deltaR2, dt, legacy=False, svd_threshold=0.2):
    """Pseudo-inverse of the lower-triangular AIF convolution matrix.

    Default: A[i, j] = dt * AIF(t_i - t_j), regularised by truncated SVD
    (singular values below svd_threshold * S_max are discarded, the standard
    sSVD approach). max(A^-1 c) is then CBF in 1/s.

    legacy=True reproduces the prototype: the AIF is z-scored, dt is dropped
    from the operator, and S_inv = 1/(S + 1e-3). rCBF is then in arbitrary
    units and MTT is not in seconds.
    """
    if legacy:
        aif = (aif_deltaR2 - np.mean(aif_deltaR2)) / (np.std(aif_deltaR2) + 1e-8)
        scale = 1.0
    else:
        aif = np.asarray(aif_deltaR2, dtype=np.float64)
        scale = dt

    nT = len(aif)
    A = np.zeros((nT, nT), dtype=np.float64)
    for i in range(nT):
        A[i:, i] = aif[: nT - i]
    A *= scale

    U, S, Vt = svd(A, full_matrices=False)
    if legacy:
        S_inv = 1.0 / (S + 1e-3)
    else:
        S_inv = np.zeros_like(S)
        keep = S > (svd_threshold * S.max())
        S_inv[keep] = 1.0 / S[keep]

    return Vt.T @ np.diag(S_inv) @ U.T


# --- main ----------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="DSC-MRI leakage correction and perfusion map estimation."
    )
    parser.add_argument("dsc_nifti", help="preprocessed 4D DSC image")
    parser.add_argument("mask_nifti", help="3D brain mask")
    parser.add_argument("output_dir", help="directory to save outputs")
    parser.add_argument("--te", type=float, required=True, help="echo time (s)")
    parser.add_argument("--tr", type=float, required=True, help="TR / sampling interval (s)")

    parser.add_argument(
        "--baseline",
        default="auto",
        help="pre-contrast frames as START:END (end exclusive), or 'auto' (default) to "
        "detect the bolus onset from the mask-average signal. Use --baseline 5:10 to "
        "reproduce the prototype's fixed window.",
    )

    parser.add_argument(
        "--apply_mask",
        action="store_true",
        help="set voxels outside the brain mask to NaN in every output",
    )
    parser.add_argument(
        "--spatial_sigma",
        type=float,
        default=0.0,
        help="spatial smoothing sigma in voxels (3D Gaussian); 0=off",
    )
    parser.add_argument(
        "--temporal_sigma",
        type=float,
        default=0.0,
        help="temporal smoothing sigma in timepoints (1D Gaussian); 0=off",
    )

    parser.add_argument(
        "--aif_mode",
        choices=("single", "mean"),
        default="single",
        help="AIF = most typical candidate voxel (default) or the candidate-pool mean",
    )
    parser.add_argument(
        "--n_aif_candidates", type=int, default=200, help="AIF candidate pool size (default 200)"
    )
    parser.add_argument(
        "--svd_threshold",
        type=float,
        default=0.2,
        help="truncated-SVD cutoff as a fraction of the largest singular value (default 0.2)",
    )
    parser.add_argument(
        "--legacy",
        action="store_true",
        help="reproduce the prototype's z-scored AIF and 1/(S+1e-3) damping "
        "(rCBF in arbitrary units, MTT not in seconds)",
    )
    parser.add_argument(
        "--save_timeseries",
        action="store_true",
        help="also write the 4D R2*(t) and dR2*(t) series (large files)",
    )
    parser.add_argument(
        "--chunk", type=int, default=50000, help="voxels processed per chunk (default 50000)"
    )

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    dsc, aff = load_nifti(args.dsc_nifti)
    mask, _ = load_nifti(args.mask_nifti)
    if dsc.ndim != 4:
        raise SystemExit(f"{args.dsc_nifti} is {dsc.ndim}D - a 4D DSC series is required.")
    if mask.shape != dsc.shape[:3]:
        raise SystemExit(
            f"Mask shape {mask.shape} does not match DSC spatial shape {dsc.shape[:3]}."
        )

    n_frames = dsc.shape[-1]
    time = np.arange(n_frames) * args.tr

    if args.baseline == "auto":
        bl_start, bl_end = detect_baseline_window(dsc, mask)
        print(f"Baseline frames detected: [{bl_start}:{bl_end}]")
    else:
        try:
            bl_start, bl_end = (int(v) for v in args.baseline.split(":"))
        except ValueError:
            raise SystemExit(f"--baseline must be 'auto' or START:END, got '{args.baseline}'.")
        if bl_end > n_frames or bl_start < 0 or bl_end - bl_start < 2:
            raise SystemExit(
                f"Baseline window [{bl_start}:{bl_end}] is invalid for a {n_frames}-frame series."
            )

    if args.spatial_sigma > 0:
        print(f"Applying spatial smoothing (sigma={args.spatial_sigma} voxels)...")
        s = args.spatial_sigma
        dsc = gaussian_filter(dsc, sigma=(s, s, s, 0.0))
    if args.temporal_sigma > 0:
        print(f"Applying temporal smoothing (sigma={args.temporal_sigma} frames)...")
        dsc = gaussian_filter1d(dsc, sigma=args.temporal_sigma, axis=-1)

    if args.save_timeseries:
        print("Computing R2*(t) and dR2*(t) series...")
        save_nifti(
            compute_r2star_no_negatives(dsc, args.te, bl_start, bl_end),
            aff,
            os.path.join(args.output_dir, "R2star_timeseries.nii.gz"),
        )
        save_nifti(
            compute_delta_r2star(dsc, args.te, bl_start, bl_end),
            aff,
            os.path.join(args.output_dir, "deltaR2star_timeseries.nii.gz"),
        )

    print("Detecting AIF...")
    aif_raw, n_candidates = detect_aif(
        dsc, mask, args.n_aif_candidates, args.aif_mode, (bl_start, bl_end)
    )
    aif_smooth = enforce_flat_baseline_tail(smooth_signal(aif_raw))
    aif_deltaR2 = compute_delta_r2star(aif_smooth, args.te, bl_start, bl_end)
    save_aif_plots(aif_raw, aif_smooth, aif_deltaR2, args.output_dir)

    aif_auc = float(simpson(aif_deltaR2, x=time))
    if aif_auc <= 0:
        raise SystemExit(
            f"The AIF dR2* curve has a non-positive area ({aif_auc:.4g}), which means no "
            "contrast passage was measured relative to the chosen baseline. Usually the "
            f"baseline window [{bl_start}:{bl_end}] overlaps the bolus, or the series is "
            "truncated. Inspect aif_raw.png / aif_deltaR2star.png in the output directory "
            "and set --baseline explicitly."
        )

    A_inv = build_deconvolution_operator(
        aif_deltaR2, args.tr, legacy=args.legacy, svd_threshold=args.svd_threshold
    )

    shape_3d = dsc.shape[:3]
    rCBV_unc = np.zeros(shape_3d, dtype=np.float64)
    rCBV_corr = np.zeros(shape_3d, dtype=np.float64)
    K2_map = np.zeros(shape_3d, dtype=np.float64)
    K1_map = np.zeros(shape_3d, dtype=np.float64)
    rCBF = np.zeros(shape_3d, dtype=np.float64)
    MTT = np.zeros(shape_3d, dtype=np.float64)
    TTP = np.zeros(shape_3d, dtype=np.float64)
    TT0 = np.zeros(shape_3d, dtype=np.float64)

    coords = np.where(mask > 0)
    n_vox = coords[0].size
    print(f"Fitting {n_vox} voxels (AIF from {n_candidates} candidates)...")

    for start in range(0, n_vox, args.chunk):
        stop = min(start + args.chunk, n_vox)
        sel = tuple(c[start:stop] for c in coords)

        curves = compute_delta_r2star(dsc[sel], args.te, bl_start, bl_end).astype(np.float64)

        unc, corr, k2, k1 = bsw_leakage_correction(curves, aif_deltaR2, time)
        rCBV_unc[sel] = unc
        rCBV_corr[sel] = corr
        K2_map[sel] = k2
        K1_map[sel] = k1

        # baseline-shift before deconvolution, as in the prototype
        residues = (curves - curves[:, bl_start:bl_end].mean(axis=1, keepdims=True)) @ A_inv.T
        cbf = residues.max(axis=1)
        rCBF[sel] = cbf

        peak = curves.max(axis=1)
        TTP[sel] = time[curves.argmax(axis=1)]

        # TT0: first frame at or above 10% of the peak (0 where the curve never rises)
        above = curves >= (0.1 * peak)[:, None]
        has_peak = peak > 0
        first_above = np.where(has_peak, above.argmax(axis=1), 0)
        TT0[sel] = np.where(has_peak, time[first_above], 0.0)

        with np.errstate(divide="ignore", invalid="ignore"):
            valid = cbf > 1e-6
            if args.legacy:
                mtt = np.where(valid, corr / np.where(valid, cbf, 1.0), 0.0)
            else:
                # CBV normalised by the AIF area is dimensionless, so CBV/CBF is in seconds
                mtt = np.where(valid, (corr / aif_auc) / np.where(valid, cbf, 1.0), 0.0)
            mtt[~np.isfinite(mtt)] = 0.0
        MTT[sel] = mtt

    outputs = {
        "rCBV_uncorrected.nii.gz": rCBV_unc,
        "rCBV_corrected.nii.gz": rCBV_corr,
        "K2.nii.gz": K2_map,
        "K1.nii.gz": K1_map,
        "rCBF.nii.gz": rCBF,
        "MTT.nii.gz": MTT,
        "TTP.nii.gz": TTP,
        "TT0.nii.gz": TT0,
    }

    if args.apply_mask:
        print("Masking outputs outside the brain mask (-> NaN)...")
        outside = mask <= 0
        for data in outputs.values():
            data[outside] = np.nan

    print("Saving parametric maps...")
    for name, data in outputs.items():
        save_nifti(data, aff, os.path.join(args.output_dir, name))

    provenance = {
        "TE_s": args.te,
        "TR_s": args.tr,
        "baseline_frames": [bl_start, bl_end],
        "baseline_mode": args.baseline,
        "n_frames": int(n_frames),
        "n_voxels_fitted": int(n_vox),
        "aif_mode": args.aif_mode,
        "n_aif_candidates": int(n_candidates),
        "aif_auc": aif_auc,
        "deconvolution": "legacy_zscored_tikhonov" if args.legacy else "dt_scaled_truncated_svd",
        "svd_threshold": None if args.legacy else args.svd_threshold,
        "spatial_sigma": args.spatial_sigma,
        "temporal_sigma": args.temporal_sigma,
        "units": {
            "rCBV_corrected": "arbitrary (integral of dR2*, normalise against NAWM)",
            "rCBV_uncorrected": "arbitrary (integral of dR2*, normalise against NAWM)",
            "rCBF": "arbitrary" if args.legacy else "1/s",
            "MTT": "arbitrary" if args.legacy else "s",
            "TTP": "s",
            "TT0": "s",
            "K1": "arbitrary",
            "K2": "arbitrary (leakage coefficient)",
        },
    }
    with open(os.path.join(args.output_dir, "KUL_dsc_fit.json"), "w") as f:
        json.dump(provenance, f, indent=4)

    print(f"All maps saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
