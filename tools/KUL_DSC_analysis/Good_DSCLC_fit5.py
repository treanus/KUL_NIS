import numpy as np
import nibabel as nib
import argparse
import os
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.signal import savgol_filter
from scipy.integrate import cumtrapz, simps
from scipy.stats import linregress
import matplotlib.pyplot as plt
from scipy.linalg import svd

# For optional smoothing:
from scipy.ndimage import gaussian_filter, gaussian_filter1d


def load_nifti(path):
    img = nib.load(path)
    return img.get_fdata(), img.affine

def save_nifti(data, affine, path):
    """
    Save 'data' to path as float32 NIfTI. 
    If 'data' contains NaNs, they will be preserved in float32 format.
    """
    nib.save(nib.Nifti1Image(data.astype(np.float32), affine), path)

def extract_time_series(dsc, mask):
    """
    Flatten all masked voxels into [N_vox, T].
    """
    return dsc[mask > 0].reshape(-1, dsc.shape[-1])

def detect_aif(dsc, mask, num_candidates=200):
    """
    Heuristic approach:
    1) Find voxels with signal drop > 100 between baseline(0..4) and min of time-series.
    2) Take the top 'num_candidates' of those voxels.
    3) Run PCA->KMeans(1 cluster) on that set; pick voxel closest to cluster center.
    """
    ts = extract_time_series(dsc, mask)
    baseline = np.mean(ts[:, :5], axis=1)
    drop = baseline - np.min(ts, axis=1)
    # Filter for significant drop
    valid = np.where(drop > 100)[0]
    top = valid[np.argsort(drop[valid])[-num_candidates:]]

    pca = PCA(n_components=2)
    reduced = pca.fit_transform(ts[top])
    kmeans = KMeans(n_clusters=1, n_init=10, random_state=42)
    kmeans.fit(reduced)

    # Pick the voxel closest to the cluster centroid
    best_idx = top[np.argmin(kmeans.transform(reduced))]
    return ts[best_idx]

def smooth_signal(sig):
    """
    Smooth a 1D signal (AIF) using Savitzky-Golay filter.
    """
    return savgol_filter(sig, 7, 2)

def enforce_flat_baseline_tail(sig, N=5):
    """
    Force the first N and last N points to be flat (mean of next 3 or previous 3).
    Helps keep the baseline and tail stable for deconvolution.
    """
    sig[:N] = np.mean(sig[N:N+3])
    sig[-N:] = np.mean(sig[-N-3:-N])
    return sig

def compute_delta_r2star(signal, TE, baseline=None, eps=1e-3):
    """
    ΔR2*(t) = -(1/TE) * ln(signal / baseline).
    If baseline=None, use frames 5..9. We clip both baseline & signal to 'eps'
    to avoid zero or negative log arguments. Then we suppress warnings
    and fix any infinite or NaN results.
    """
    if baseline is None:
        baseline = np.mean(signal[..., 5:10], axis=-1, keepdims=True)
    
    # Ensure no zero or near-zero values
    baseline = np.clip(baseline, eps, None)
    signal   = np.clip(signal,   eps, None)

    with np.errstate(divide='ignore', invalid='ignore', over='ignore', under='ignore'):
        ratio = signal / baseline
        delta = -np.log(ratio) / TE
        # Replace any inf, -inf, or NaN with 0 (or np.nan if you prefer)
        delta[~np.isfinite(delta)] = 0.0

    return delta


def compute_r2star_no_negatives(signal, TE, baseline_start=5, baseline_end=10):
    """
    Compute R2*(t) = (1/TE) * ln( baseline / signal(t) ),
    using frames [baseline_start..baseline_end-1] as the baseline.
    
    If the formula would produce negative values (signal>baseline),
    we set R2*(t)=0 instead of clipping arbitrarily. 
    """
    # Baseline from the second 5 volumes
    baseline = np.mean(signal[..., baseline_start:baseline_end], axis=-1, keepdims=True)
    # Prevent zero or negative baseline
    baseline = np.clip(baseline, 1e-3, None)
    # Prevent zero or negative signals
    signal   = np.clip(signal,   1e-3, None)

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = baseline / signal
        r2star = (1.0 / TE) * np.log(ratio)
        # Where ratio < 1 => log(ratio) negative => set to 0
        r2star[r2star < 0] = 0
        # Replace inf or NaN with 0
        r2star[~np.isfinite(r2star)] = 0.0

    return r2star

def bsw_leakage_correction(deltaR2, aif_deltaR2, mask, time):
    """
    BSW leakage correction (voxelwise).
      - For each masked voxel, do linear regression of tissue_curve on aif_int(t).
      - Subtract the leakage component, then integrate (simps).
      - slope (K2) is the leakage coefficient, intercept (K1).
    """
    dt = time[1] - time[0]
    aif_int = cumtrapz(aif_deltaR2, time, initial=0)

    shape = deltaR2.shape[:3]
    rCBV_unc = np.zeros(shape)
    rCBV_corr = np.zeros(shape)
    K2_map = np.zeros(shape)
    K1_map = np.zeros(shape)

    coords = np.where(mask > 0)
    for x, y, z in zip(*coords):
        tissue_curve = deltaR2[x, y, z, :]
        tissue_auc = simps(tissue_curve, time)

        # Fit: tissue_curve(t) ~ slope * aif_int(t) + intercept
        slope, intercept, _, _, _ = linregress(aif_int, tissue_curve)
        leakage_component = slope * aif_int
        corrected_curve = tissue_curve - leakage_component
        corrected_auc = simps(corrected_curve, time)

        rCBV_unc[x, y, z] = tissue_auc
        rCBV_corr[x, y, z] = corrected_auc
        K2_map[x, y, z] = slope
        K1_map[x, y, z] = intercept

    return rCBV_unc, rCBV_corr, K2_map, K1_map

def save_aif_plots(raw, smoothed, deltaR2, output_dir):
    """
    Diagnostic plots for the AIF:
      - raw signal
      - smoothed signal
      - ΔR2*
    """
    plt.figure()
    plt.plot(raw)
    plt.title("Raw AIF")
    plt.xlabel("Timepoints")
    plt.ylabel("Signal Intensity")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "aif_raw.png"))
    plt.close()

    plt.figure()
    plt.plot(smoothed)
    plt.title("Smoothed AIF")
    plt.xlabel("Timepoints")
    plt.ylabel("Signal Intensity")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "aif_smoothed.png"))
    plt.close()

    plt.figure()
    plt.plot(deltaR2)
    plt.title("ΔR2* AIF")
    plt.xlabel("Timepoints")
    plt.ylabel("ΔR2* (1/s)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "aif_deltaR2star.png"))
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dsc_nifti", help="4D DSC image")
    parser.add_argument("mask_nifti", help="3D brain mask")
    parser.add_argument("output_dir", help="Directory to save outputs")
    parser.add_argument("--te", type=float, default=0.03, help="Echo time (s)")
    parser.add_argument("--tr", type=float, default=1.5, help="TR / sampling interval (s)")
    
    # 1) Optional masking
    parser.add_argument(
        "--apply_mask",
        action="store_true",
        help="If set, all output voxels outside the mask will be set to NaN."
    )
    
    # 2) Optional smoothing
    parser.add_argument(
        "--spatial_sigma",
        type=float,
        default=0.0,
        help="Spatial smoothing sigma in voxels (3D Gaussian). 0=off."
    )
    parser.add_argument(
        "--temporal_sigma",
        type=float,
        default=0.0,
        help="Temporal smoothing sigma in timepoints (1D Gaussian). 0=off."
    )
    
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load data
    dsc, aff = load_nifti(args.dsc_nifti)
    mask, _ = load_nifti(args.mask_nifti)
    time = np.arange(dsc.shape[-1]) * args.tr

    # --- (Optional) apply spatial and/or temporal smoothing to dsc
    #     1) Spatial smoothing with sigma=(sx, sy, sz, 0)
    if args.spatial_sigma > 0:
        print(f"Applying spatial smoothing with sigma={args.spatial_sigma}...")
        dsc = gaussian_filter(
            dsc, sigma=(args.spatial_sigma, args.spatial_sigma, args.spatial_sigma, 0.0)
        )

    #     2) Temporal smoothing along last axis
    if args.temporal_sigma > 0:
        print(f"Applying temporal smoothing with sigma={args.temporal_sigma}...")
        dsc = gaussian_filter1d(dsc, sigma=args.temporal_sigma, axis=-1)

    # --- Compute and save 4D R2*(t) and 4D ΔR2*(t).
    #     Use our 'compute_r2star_positive' to avoid negative R2* values.
    print("Computing R2* (4D) with positivity enforced...")
    r2star_4d = compute_r2star_no_negatives(dsc, TE=args.te, 
                                        baseline_start=5, baseline_end=10)

    print("Computing ΔR2* (4D) from baseline points...")
    delta_r2star_4d = compute_delta_r2star(dsc, TE=args.te)

    # (Optional) Mask out-of-brain => NaN
    if args.apply_mask:
        print("Masking R2* and ΔR2* timeseries outside the brain mask (-> NaN).")
        r2star_4d[mask <= 0, ...] = np.nan
        delta_r2star_4d[mask <= 0, ...] = np.nan

    # Save them
    save_nifti(r2star_4d, aff, os.path.join(args.output_dir, "R2star_timeseries.nii.gz"))
    save_nifti(delta_r2star_4d, aff, os.path.join(args.output_dir, "deltaR2star_timeseries.nii.gz"))

    # --- Detect AIF
    print("📌 Detecting AIF...")
    aif_raw = detect_aif(dsc, mask)  # raw intensities
    aif_smooth = enforce_flat_baseline_tail(smooth_signal(aif_raw))
    aif_deltaR2 = compute_delta_r2star(aif_smooth, args.te)  # 1D

    print("💾 Saving AIF plots...")
    save_aif_plots(aif_raw, aif_smooth, aif_deltaR2, args.output_dir)

    # --- Compute ΔR2* for all voxels (already done above, but we keep it for clarity)
    print("📌 Computing ΔR2* for all voxels...")
    deltaR2 = compute_delta_r2star(dsc, args.te)

    # --- BSW leakage correction
    print("📌 Performing BSW leakage correction...")
    rCBV_unc, rCBV_corr, K2, K1 = bsw_leakage_correction(deltaR2, aif_deltaR2, mask, time)

    # --- Deconvolution for rCBF, MTT, TTP, TT0
    print("📌 Preparing AIF convolution matrix (standardized)...")
    aif_std = (aif_deltaR2 - np.mean(aif_deltaR2)) / (np.std(aif_deltaR2) + 1e-8)
    nT = len(aif_std)

    A = np.zeros((nT, nT), dtype=np.float32)
    for i in range(nT):
        A[i:, i] = aif_std[:(nT - i)]

    U, S, Vt = svd(A, full_matrices=False)
    S_inv = np.diag(1 / (S + 1e-3))
    A_inv = Vt.T @ S_inv @ U.T

    shape_3d = deltaR2.shape[:3]
    rCBF = np.zeros(shape_3d)
    MTT = np.zeros(shape_3d)
    TTP = np.zeros(shape_3d)
    TT0 = np.zeros(shape_3d)

    print("🚀 Computing rCBF, MTT, TTP, TT0...")
    coords = np.where(mask > 0)
    for x, y, z in zip(*coords):
        curve = deltaR2[x, y, z, :]
        if np.all(curve == 0):
            continue

        # Baseline shift
        curve_bs = curve - np.mean(curve[:5])

        residue = A_inv @ curve_bs
        rCBF[x, y, z] = np.max(residue)

        # TTP from raw ΔR2*(t) curve
        TTP[x, y, z] = time[np.argmax(curve)]

        # TT0: first time above 10% of max
        peak_val = np.max(curve)
        if peak_val > 0:
            threshold = 0.1 * peak_val
            idxs = np.where(curve >= threshold)[0]
            if len(idxs) > 0:
                TT0[x, y, z] = time[idxs[0]]

    # MTT = rCBV_corr / rCBF (avoid division by 0)
    valid_mask = (rCBF > 1e-6)
    MTT[valid_mask] = rCBV_corr[valid_mask] / rCBF[valid_mask]

    # --- (Optional) apply mask => out-of-mask = NaN
    if args.apply_mask:
        print("Masking final parametric maps outside the brain mask (-> NaN).")
        rCBV_unc[mask <= 0] = np.nan
        rCBV_corr[mask <= 0] = np.nan
        K2[mask <= 0] = np.nan
        K1[mask <= 0] = np.nan
        rCBF[mask <= 0] = np.nan
        MTT[mask <= 0] = np.nan
        TTP[mask <= 0] = np.nan
        TT0[mask <= 0] = np.nan

    print("💾 Saving leakage correction outputs and final maps...")
    save_nifti(rCBV_unc, aff, os.path.join(args.output_dir, "rCBV_uncorrected.nii.gz"))
    save_nifti(rCBV_corr, aff, os.path.join(args.output_dir, "rCBV_corrected.nii.gz"))
    save_nifti(K2, aff, os.path.join(args.output_dir, "K2.nii.gz"))
    save_nifti(K1, aff, os.path.join(args.output_dir, "K1.nii.gz"))
    save_nifti(rCBF, aff, os.path.join(args.output_dir, "rCBF.nii.gz"))
    save_nifti(MTT, aff, os.path.join(args.output_dir, "MTT.nii.gz"))
    save_nifti(TTP, aff, os.path.join(args.output_dir, "TTP.nii.gz"))
    save_nifti(TT0, aff, os.path.join(args.output_dir, "TT0.nii.gz"))

    print("✅ All maps saved to:", args.output_dir)


if __name__ == "__main__":
    main()
