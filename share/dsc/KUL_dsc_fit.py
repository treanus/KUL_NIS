#!/usr/bin/env python3
"""DSC-MRI quantification: dR2*, BSW leakage correction, oSVD deconvolution.

Reads a motion-corrected 4D DSC series ON ITS NATIVE GRID plus a brain mask, and
writes rCBV (leakage-corrected and uncorrected), rCBF, MTT, Tmax, TTP, TT0, the
BSW coefficients K1/K2 and a bolus-CNR quality map. KUL_dsc_perfusion.sh then
moves the 3D maps to T1w space with one interpolation.

Method (version 2), with the evidence each choice rests on:

  timing    Baseline, bolus arrival and end of first pass come from the
            brain-average signal. The baseline is every steady-state frame up to
            2 frames before arrival; arrival is where the signal leaves the early
            pre-bolus level by more than max(3 sigma, 10% of the bolus depth).
            This tolerates the slow pre-bolus drift that cut the v1 window to 7
            frames. More baseline frames = less CBV noise (Boxerman 1997, JMRI, PMID 9170038;
            ASFNR, Welker 2015 AJNR, PMID 25907520).

  leakage   Boxerman-Schmainda-Weisskoff (Boxerman 2006 AJNR, PMID 16611779), the
            correction with the most simulation and clinical evidence (2020
            consensus, Boxerman 2020 Neuro-Oncol, PMID 32516388):
                dR2*(t) = K1*ref(t) - K2*INT_0^t ref
                dR2*_corr(t) = dR2*(t) + K2*INT_0^t ref
            ref = mean dR2* of non-enhancing brain: outside --exclude_mask, not
            T1-leaky, not CSF/noise, not vessels. K2 > 0 = T1-dominant leakage
            (rCBV underestimated before correction), K2 < 0 = T2*-dominant.
            v1 regressed on INT(AIF) of one noisy voxel without the K1*ref term,
            so K2 absorbed recirculation and depended on the AIF voxel.

  rCBV      trapezoidal integral of dR2* from 2 frames before arrival to the end
            of the series ('full', the ASFNR default) or to the end of the first
            pass plus margin ('firstpass'). Relative units: read as a ratio to
            NAWM. Not AIF-normalised, so it does not depend on the AIF at all.

  AIF       automatic cluster selection (Mouridsen 2006 MRM, PMID 16453314): the 10%
            of voxels with the largest first-pass area, minus the 25% roughest,
            k-means on area-normalised curves, keep the cluster with the
            earliest first moment, twice. AIF = mean of that cluster (several
            voxels, no Savitzky-Golay over the peak, no forced flat tail).

  rCBF      block-circulant SVD with per-voxel oscillation-index regularisation
            (oSVD; Wu 2003 MRM, PMID 12815691), which is insensitive to tracer delay.
            Deconvolves the leakage-corrected curve. Tmax = time of the residue
            maximum. MTT = CBV/CBF with CBV = INT(corrected)/INT(AIF) over the
            same window (central volume theorem), in seconds.

  timing    TTP and TT0 come from the temporally smoothed (sigma 1 frame)
  maps      corrected curve, searched only inside the first-pass window, with
            sub-frame interpolation; TT0 walks back from the peak, so pre-bolus
            noise cannot trigger it (v1 put TT0 in pre-bolus frames in 96% of
            voxels). Voxels without a detectable bolus (CNR < 3) get 0.

  smoothing none on the primary maps. --smooth_fwhm_mm writes additional
            *_smooth maps (normalised Gaussian convolution inside the mask),
            labelled as such (ASFNR: "a small kernel may be considered").

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
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.ndimage import gaussian_filter, gaussian_filter1d, median_filter

try:
    from sklearn.cluster import KMeans
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "scikit-learn is required for AIF selection. Install it in the conda env "
        "used for this step (the KUL_NIS 'pyfMRI' env gets it via nilearn)."
    ) from exc

METHOD_VERSION = "2.0"


def log(msg):
    print(msg, flush=True)


# --- I/O -----------------------------------------------------------------


def load_nifti(path):
    img = nib.load(path)
    return np.nan_to_num(img.get_fdata(dtype=np.float32)), img


def save_map(values, mask, ref_img, path, fill=0.0):
    out = np.full(mask.shape, fill, dtype=np.float32)
    out[mask] = values
    hdr = ref_img.header.copy()
    hdr.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(out, ref_img.affine, hdr), path)


def mad_sigma(x):
    x = np.asarray(x, dtype=np.float64)
    return 1.4826 * np.median(np.abs(x - np.median(x)))


# --- timing from the brain-average curve -----------------------------------


def detect_timing(curve, edge=3, guard=2, max_skip=3):
    """Baseline window, bolus arrival and bolus peak from the mean signal curve.

    Returns (ss, bl_end, arrival, peak, sigma, info). Frames [ss, bl_end) are the
    baseline. Raises ValueError when no bolus can be found.
    """
    curve = np.asarray(curve, dtype=np.float64)
    n = curve.size
    if n < 2 * edge + 8:
        raise ValueError(f"only {n} frames - too short for a DSC series")

    peak = edge + int(np.argmin(curve[edge : n - edge]))
    if peak < 8:
        raise ValueError(f"the bolus peaks at frame {peak}, leaving no usable baseline")

    pre = curve[:peak]
    sigma = mad_sigma(np.diff(pre)) / np.sqrt(2.0)
    early = curve[1 : max(4, peak // 2)]
    level = float(np.median(early))
    depth = level - curve[peak]
    if depth <= 0:
        raise ValueError("the mean signal never drops below its pre-bolus level")

    thr = level - max(3.0 * sigma, 0.10 * depth)
    i = peak
    while i > 0 and curve[i] < thr:
        i -= 1
    arrival = i + 1

    # frame 0 is never at steady state; skip a few more only if clearly brighter
    ss = 1
    while ss < 1 + max_skip and curve[ss] - level > max(4.0 * sigma, 0.25 * depth):
        ss += 1

    bl_end = arrival - guard
    if bl_end - ss < 3:
        raise ValueError(
            f"only {max(bl_end - ss, 0)} baseline frames before the bolus (arrival frame {arrival})"
        )

    # report (not correct) the baseline trend: a drift over the baseline that is
    # a sizeable fraction of the bolus depth makes S0 ambiguous
    tb = np.arange(ss, bl_end)
    slope = np.polyfit(tb, curve[ss:bl_end], 1)[0]
    info = {
        "bolus_depth_fraction": float(depth / level),
        "baseline_drift_fraction_of_depth": float(slope * (bl_end - ss) / depth),
        "baseline_noise_fraction": float(sigma / level),
    }
    return ss, bl_end, arrival, peak, sigma, info


def first_pass_end(ref, arrival, peak_ref, n):
    """End of the first pass on the reference dR2* curve: the recirculation trough.

    First local minimum after the peak once the curve has fallen below half its
    peak, capped at peak + 2*(peak - arrival) + 2 frames.
    """
    r = gaussian_filter1d(np.asarray(ref, dtype=np.float64), 1.0)
    cap = min(n - 1, peak_ref + 2 * max(peak_ref - arrival, 2) + 2)
    half = 0.5 * r[peak_ref]
    for k in range(peak_ref + 1, cap):
        if r[k] < half and r[k] <= r[k + 1]:
            return k
    return cap


# --- dR2* ----------------------------------------------------------------


def delta_r2star(S, S0, TE, eps=1e-3):
    with np.errstate(divide="ignore", invalid="ignore"):
        d = -np.log(np.clip(S, eps, None) / np.clip(S0, eps, None)[:, None]) / TE
    d[~np.isfinite(d)] = 0.0
    return d


# --- reference tissue and BSW ----------------------------------------------


def reference_tissue(d, S0, excl, fp, late):
    """Mask of non-enhancing brain voxels for the BSW reference curve."""
    ok = ~excl & (S0 > 0.5 * np.median(S0[~excl]))
    auc = trapezoid(d[:, fp], axis=1)
    lo, hi = np.percentile(auc[ok], [20, 95])  # drop CSF/noise and vessels
    ok &= (auc > lo) & (auc < hi)
    tail = d[:, late].mean(axis=1)
    med, sig = np.median(tail[ok]), mad_sigma(tail[ok])
    ok &= (tail > med - 2.0 * sig) & (tail < med + 3.0 * sig)  # T1-leaky / T2*-leaky out
    return ok


def shift_curve(c, t, delay):
    """c(t - delay) by linear interpolation; zero before the start."""
    return np.interp(t - delay, t, c, left=0.0, right=c[-1])


def estimate_delay(d, ref, t, fp, cnr_ok, shifts):
    """Per-voxel bolus delay relative to the reference: the shift of ref that best
    correlates with the voxel's first pass. Voxels without a detectable bolus get 0."""
    x = d[:, fp] - d[:, fp].mean(axis=1, keepdims=True)
    xn = np.sqrt((x ** 2).sum(axis=1)) + 1e-12
    best = np.full(d.shape[0], -np.inf)
    delay = np.zeros(d.shape[0])
    for s in shifts:
        r = shift_curve(ref, t, s)[fp]
        r = r - r.mean()
        cc = (x @ r) / (xn * (np.sqrt(r @ r) + 1e-12))
        better = cc > best
        best[better] = cc[better]
        delay[better] = s
    return np.where(cnr_ok, delay, 0.0)


def bsw_fit(d, ref, t, fit_idx, delay=None):
    """Per-voxel OLS of dR2* on [ref, -INT ref]. Returns K1, K2 and the correction term
    K2*INT ref(t) per voxel.

    With 'delay' (s, per voxel), the reference is shifted to each voxel's bolus
    arrival before fitting (arrival-time correction, Leigh 2012 PLoS One,
    PMID 23285132): an unshifted reference makes a late-arriving bolus look like
    T2*-leakage and an early one like T1-leakage.
    """
    if delay is None:
        delay = np.zeros(d.shape[0])
    K1 = np.zeros(d.shape[0])
    K2 = np.zeros(d.shape[0])
    term = np.zeros_like(d)
    for s in np.unique(delay):
        v = delay == s
        r = shift_curve(ref, t, s)
        ir = cumulative_trapezoid(r, t, initial=0.0)
        X = np.stack([r[fit_idx], -ir[fit_idx]], axis=1)
        coef = np.linalg.pinv(X) @ d[v][:, fit_idx].T
        K1[v], K2[v] = coef[0], coef[1]
        term[v] = coef[1][:, None] * ir[None, :]
    return K1, K2, term


# --- AIF -------------------------------------------------------------------


def select_aif(d, S0, excl, fp, t, n_clusters=5, pool_frac=0.10, rough_frac=0.25, min_vox=3,
               max_vox=20):
    """Cluster-based automatic AIF (Mouridsen 2006). Returns (aif, voxel_indices, info)."""
    cand = np.where(~excl & (S0 > 0.25 * np.median(S0[~excl])))[0]
    dfp = d[cand][:, fp]
    auc = trapezoid(dfp, axis=1)
    pool = cand[auc >= np.percentile(auc, 100 * (1 - pool_frac))]

    c = d[pool][:, fp]
    peak = c.max(axis=1)
    rough = np.sum(np.diff(c, 2, axis=1) ** 2, axis=1) / np.maximum(peak, 1e-6) ** 2
    pool = pool[rough <= np.percentile(rough, 100 * (1 - rough_frac))]

    tfp = t[fp]
    sel = pool
    for _ in range(2):
        if sel.size < max(min_vox, n_clusters) * 2:
            break
        c = np.clip(d[sel][:, fp], 0, None)
        area = trapezoid(c, tfp, axis=1)
        keep = area > 0
        sel, c, area = sel[keep], c[keep], area[keep]
        shape = c / area[:, None]
        km = KMeans(n_clusters=n_clusters, n_init=10, random_state=0).fit(shape)
        m1 = (km.cluster_centers_ * tfp).sum(axis=1) / km.cluster_centers_.sum(axis=1)
        best = None
        for k in np.argsort(m1):  # earliest first moment with enough voxels
            if np.sum(km.labels_ == k) >= min_vox:
                best = k
                break
        if best is None:
            break
        sel = sel[km.labels_ == best]

    # within the earliest cluster, keep the highest-peaking voxels: the cluster
    # also holds partial-volumed and smaller vessels, which lower and broaden
    # the AIF (peak criterion as in Bjornerud & Emblem 2010 JCBFM, PMID 20087370)
    n_cluster = int(sel.size)
    if sel.size > max_vox:
        pk = d[sel][:, fp].max(axis=1)
        sel = sel[np.argsort(pk)[-max_vox:]]
    aif = d[sel].mean(axis=0)
    info = {"n_pool": int(pool.size), "n_cluster": n_cluster, "n_voxels": int(sel.size)}
    return aif, sel, info


# --- oSVD deconvolution ------------------------------------------------------


def osvd_deconvolve(curves, aif, dt, oi_target=0.095, chunk=8000):
    """Block-circulant SVD with per-voxel oscillation-index thresholding (Wu 2003).

    Returns CBF (max of the residue, in 1/s with a dt-scaled AIF), Tmax (s) and
    the chosen SVD threshold per voxel.
    """
    T = aif.size
    L = 2 * T
    a = np.concatenate([aif, np.zeros(T)])
    idx = (np.arange(L)[:, None] - np.arange(L)[None, :]) % L
    D = dt * a[idx]
    U, S, Vt = np.linalg.svd(D)
    thresholds = np.geomspace(0.005, 0.6, 24)  # ascending
    inverses = []
    for p in thresholds:
        s_inv = np.where(S > p * S.max(), 1.0 / S, 0.0)
        inverses.append((Vt.T * s_inv) @ U.T)
    inverses = np.stack(inverses)  # (P, L, L)

    n = curves.shape[0]
    cbf = np.zeros(n)
    tmax = np.zeros(n)
    used = np.full(n, thresholds[-1])
    for s in range(0, n, chunk):
        c = np.concatenate([curves[s : s + chunk], np.zeros((min(chunk, n - s), T))], axis=1)
        r = np.einsum("pij,nj->pni", inverses, c)  # (P, N, L)
        rmax = r.max(axis=2)
        second = np.abs(r[:, :, 2:] - 2 * r[:, :, 1:-1] + r[:, :, :-2]).sum(axis=2)
        with np.errstate(divide="ignore", invalid="ignore"):
            oi = second / (L * rmax)
        ok = (oi <= oi_target) & (rmax > 0)
        # smallest threshold (least regularisation) that is not oscillating
        first = np.where(ok.any(axis=0), ok.argmax(axis=0), len(thresholds) - 1)
        cols = np.arange(first.size)
        rs = r[first, cols]  # (N, L)
        cbf[s : s + chunk] = rs.max(axis=1)
        k = rs.argmax(axis=1)
        tmax[s : s + chunk] = np.where(k < T, k, k - L) * dt  # wrap = tissue before AIF
        used[s : s + chunk] = thresholds[first]
    return np.clip(cbf, 0, None), tmax, used


# --- per-voxel timing --------------------------------------------------------


def smoothed_peak_cnr(curves, bl, search):
    """Peak of the temporally smoothed curve inside 'search' over the SD of the
    smoothed baseline - the same smoothing on both sides, so CNR >= 3 means a
    bolus is detectable."""
    cs = gaussian_filter1d(curves, 1.0, axis=1)
    sig = cs[:, bl].std(axis=1, ddof=1)
    pk = cs[:, search].max(axis=1)
    return cs, sig, pk / np.maximum(sig, 1e-9)


def bolus_timing(corr, bl, search, t, dt):
    """TTP, TT0 (s from series start) and bolus CNR from the corrected curves."""
    cs, sig_s, cnr = smoothed_peak_cnr(corr, bl, search)
    w = cs[:, search]
    k = w.argmax(axis=1)
    pk = w[np.arange(w.shape[0]), k]

    # parabolic sub-frame peak
    km = np.clip(k, 1, w.shape[1] - 2)
    y0, y1, y2 = (w[np.arange(w.shape[0]), km + o] for o in (-1, 0, 1))
    den = y0 - 2 * y1 + y2
    off = np.where(np.abs(den) > 1e-9, 0.5 * (y0 - y2) / den, 0.0)
    off = np.where((k > 0) & (k < w.shape[1] - 1), np.clip(off, -0.5, 0.5), 0.0)
    ttp = t[search][0] + (k + off) * dt

    # TT0: walk back from the peak to the last sample below the arrival threshold
    thr = np.maximum(3.0 * sig_s, 0.15 * pk)
    below = w < thr[:, None]
    pos = np.arange(w.shape[1])[None, :]
    cand = np.where(below & (pos < k[:, None]), pos, -1)
    j = cand.max(axis=1)  # last below-threshold frame before the peak
    j_ok = j >= 0
    jj = np.clip(j, 0, w.shape[1] - 2)
    ya = w[np.arange(w.shape[0]), jj]
    yb = w[np.arange(w.shape[0]), jj + 1]
    frac = np.where(yb > ya, (thr - ya) / np.maximum(yb - ya, 1e-9), 0.0)
    tt0 = t[search][0] + (jj + np.clip(frac, 0, 1)) * dt
    tt0 = np.where(j_ok, tt0, t[search][0])

    detect = cnr >= 3.0
    return np.where(detect, ttp, 0.0), np.where(detect, tt0, 0.0), cnr


# --- smoothing ---------------------------------------------------------------


def smooth_in_mask(vol, mask, fwhm_mm, zooms):
    sig = [fwhm_mm / (2.3548 * z) for z in zooms[:3]]
    num = gaussian_filter(np.where(mask, vol, 0.0), sig)
    den = gaussian_filter(mask.astype(np.float64), sig)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(mask & (den > 1e-3), num / den, 0.0)
    return out


# --- QC ------------------------------------------------------------------------


def qc_figure(path, t, mean_sig, ss, bl_end, arrival, fp_end, win, ref, aif, aif_curves):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    ax[0].plot(t, mean_sig / mean_sig[ss:bl_end].mean(), "k")
    ax[0].axvspan(t[ss], t[bl_end - 1], color="tab:blue", alpha=0.2, label="baseline")
    ax[0].axvline(t[arrival], color="tab:red", ls="--", label="arrival")
    ax[0].axvline(t[fp_end], color="tab:green", ls="--", label="end first pass")
    ax[0].axvspan(t[win[0]], t[win[-1]], color="tab:orange", alpha=0.08, label="integration")
    ax[0].set_title("brain-average S/S0")
    ax[0].legend(fontsize=7)
    for c in aif_curves[:50]:
        ax[1].plot(t, c, color="0.75", lw=0.5)
    ax[1].plot(t, aif, "r", lw=2, label=f"AIF (mean of {len(aif_curves)} voxels)")
    ax[1].plot(t, ref, "b", lw=2, label="BSW reference tissue")
    ax[1].set_title("dR2* (1/s)")
    ax[1].legend(fontsize=7)
    ax[2].plot(t, ref / max(ref.max(), 1e-6), "b", label="reference (peak-normalised)")
    ax[2].plot(t, aif / max(aif.max(), 1e-6), "r", label="AIF (peak-normalised)")
    ax[2].set_title("shape: the AIF should be earlier and narrower")
    ax[2].legend(fontsize=7)
    for a in ax:
        a.set_xlabel("time (s)")
        a.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=90)
    plt.close()


# --- main ------------------------------------------------------------------------


def parse_frames(spec, name, n):
    try:
        a, b = (int(v) for v in spec.split(":"))
    except ValueError:
        raise SystemExit(f"--{name} must be 'auto' or START:END, got '{spec}'.")
    if a < 0 or b > n or b - a < 3:
        raise SystemExit(f"--{name} [{a}:{b}] is invalid for a {n}-frame series.")
    return a, b


def main():
    p = argparse.ArgumentParser(description="DSC-MRI perfusion fit (BSW leakage correction, oSVD).")
    p.add_argument("dsc_nifti", help="motion-corrected 4D DSC series, native grid")
    p.add_argument("mask_nifti", help="3D brain mask on the same grid")
    p.add_argument("output_dir")
    p.add_argument("--te", type=float, required=True, help="echo time (s)")
    p.add_argument("--tr", type=float, required=True, help="frame interval (s)")
    p.add_argument("--exclude_mask", help="3D mask (e.g. dilated lesion) kept out of the "
                   "BSW reference and the AIF search; same grid")
    p.add_argument("--baseline", default="auto", help="'auto' or START:END (end exclusive)")
    p.add_argument("--window", default="full",
                   help="rCBV integration: 'full' (arrival-2 to end; default), 'firstpass' "
                   "(arrival-2 to end of first pass + 3 frames), or START:END")
    p.add_argument("--oi_target", type=float, default=0.095, help="oSVD oscillation index (Wu 2003)")
    p.add_argument("--max_delay_s", type=float, default=6.0,
                   help="largest bolus delay searched for the arrival-time-corrected BSW")
    p.add_argument("--no_delay_correction", action="store_true",
                   help="plain BSW with an unshifted reference (Boxerman 2006 as published)")
    p.add_argument("--no_drift_adjust", action="store_true",
                   help="S0 = plain baseline mean, no levelling to the pre-arrival signal")
    p.add_argument("--aif_clusters", type=int, default=5)
    p.add_argument("--smooth_fwhm_mm", type=float, default=0.0,
                   help="also write *_smooth maps with this 3D Gaussian FWHM; 0 = off")
    p.add_argument("--save_timeseries", action="store_true",
                   help="also write the 4D dR2* and leakage-corrected dR2* series")
    args = p.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    out = lambda name: os.path.join(args.output_dir, name)

    dsc, img = load_nifti(args.dsc_nifti)
    maskvol, mimg = load_nifti(args.mask_nifti)
    if dsc.ndim != 4:
        raise SystemExit(f"{args.dsc_nifti} is {dsc.ndim}D - a 4D DSC series is required.")
    if maskvol.shape != dsc.shape[:3] or not np.allclose(mimg.affine, img.affine, atol=1e-3):
        raise SystemExit("The mask is not on the DSC grid (shape or affine differ).")
    mask = maskvol > 0
    excl_vol = np.zeros(mask.shape, bool)
    if args.exclude_mask:
        ev, eimg = load_nifti(args.exclude_mask)
        if ev.shape != mask.shape:
            raise SystemExit("--exclude_mask is not on the DSC grid.")
        excl_vol = ev > 0

    T = dsc.shape[-1]
    dt = args.tr
    t = np.arange(T) * dt
    S = dsc[mask].astype(np.float64)
    excl = excl_vol[mask]
    log(f"{mask.sum()} brain voxels, {excl.sum()} excluded from reference/AIF, {T} frames")

    # timing
    mean_sig = S.mean(axis=0)
    try:
        ss, bl_end, arrival, peak, _, tinfo = detect_timing(mean_sig)
    except ValueError as e:
        if args.baseline == "auto":
            raise SystemExit(f"Automatic timing failed: {e}. Give --baseline START:END.")
        ss = bl_end = arrival = peak = None
        tinfo = {}
    if args.baseline != "auto":
        ss, bl_end = parse_frames(args.baseline, "baseline", T)
        if arrival is None or arrival < bl_end:
            arrival = bl_end + 2
            peak = arrival + int(np.argmin(mean_sig[arrival:]))
    log(f"Baseline frames [{ss}:{bl_end}], bolus arrival frame {arrival}, mean-signal peak frame {peak}")
    if abs(tinfo.get("baseline_drift_fraction_of_depth", 0.0)) > 0.25:
        log("WARNING: the pre-bolus signal drifts by more than 25% of the bolus depth; "
            "S0 (and so every map) is less certain. See qc_timing.png.")

    # S0: averaged over the whole baseline for noise, but levelled to the signal
    # just before arrival using the (high-SNR) brain-average curve, so a slow
    # pre-bolus drift does not leave a constant dR2* offset that the integral
    # then multiplies by the window length. g = 1 when there is no drift.
    g = 1.0
    if not args.no_drift_adjust:
        n_last = min(5, bl_end - ss)
        g = float(mean_sig[bl_end - n_last : bl_end].mean() / mean_sig[ss:bl_end].mean())
    S0 = g * S[:, ss:bl_end].mean(axis=1)
    log(f"S0 drift adjustment factor {g:.4f}")
    d = delta_r2star(S, S0, args.te)
    d[:, :ss] = 0.0  # non-steady-state frames carry no information
    bl = np.arange(ss, bl_end)

    start = max(arrival - 2, ss)
    fp_guess = np.arange(start, min(T, peak + 2 * max(peak - arrival, 3) + 3))
    late = np.arange(max(T - max(T // 6, 5), fp_guess[-1] + 1), T)
    if late.size < 3:
        late = np.arange(T - 3, T)

    # reference tissue and end of first pass
    ref_ok = reference_tissue(d, S0, excl, fp_guess, late)
    if ref_ok.sum() < 200:
        raise SystemExit(f"Only {ref_ok.sum()} reference-tissue voxels - check the mask.")
    ref = d[ref_ok].mean(axis=0)
    peak_ref = start + int(np.argmax(ref[start:]))
    fp_end = first_pass_end(ref, arrival, peak_ref, T)
    fp = np.arange(start, fp_end + 1)
    log(f"Reference tissue: {ref_ok.sum()} voxels, peak frame {peak_ref}, end of first pass frame {fp_end}")

    if args.window == "full":
        win = np.arange(start, T)
    elif args.window == "firstpass":
        win = np.arange(start, min(T, fp_end + 4))
    else:
        a, b = parse_frames(args.window, "window", T)
        win = np.arange(a, b)
    log(f"rCBV integration frames [{win[0]}:{win[-1] + 1}] ({args.window})")

    # bolus delay relative to the reference, for the arrival-time-corrected BSW.
    # Searched in half-frame steps; spatially median-filtered because delays are
    # smooth in space and the per-voxel estimate is noisy.
    search = np.arange(max(start - 2, ss), min(T, fp_end + 4))
    _, _, cnr_unc = smoothed_peak_cnr(d, bl, search)
    delay = np.zeros(d.shape[0])
    if not args.no_delay_correction:
        shifts = np.arange(-2.0, args.max_delay_s + 1e-6, dt / 2)
        delay = estimate_delay(d, ref, t, fp, cnr_unc >= 3, shifts)
        dv = np.zeros(mask.shape)
        dv[mask] = delay
        dv = median_filter(dv, size=(3, 3, 1))
        delay = np.round(dv[mask] / (dt / 2)) * (dt / 2)
        log(f"Bolus delay vs reference: median {np.median(delay):.2f} s, "
            f"{np.mean(np.abs(delay) > 0):.0%} of voxels shifted")

    # BSW leakage correction, fitted over every steady-state frame
    fit_idx = np.arange(ss, T)
    K1, K2, leak = bsw_fit(d, ref, t, fit_idx, delay)
    corr = d + leak
    rcbv_unc = trapezoid(d[:, win], t[win], axis=1)
    rcbv_corr = trapezoid(corr[:, win], t[win], axis=1)

    # AIF and oSVD
    aif, aif_idx, ainfo = select_aif(d, S0, excl, fp, t, n_clusters=args.aif_clusters)
    aif_ttp = int(np.argmax(aif[start:])) + start
    aif_auc = float(trapezoid(aif[win], t[win]))
    checks = {
        "aif_peak_over_ref_peak": float(aif.max() / max(ref.max(), 1e-6)),
        "aif_ttp_frame": aif_ttp,
        "ref_ttp_frame": peak_ref,
    }
    log(f"AIF: mean of {ainfo['n_voxels']} voxels (pool {ainfo['n_pool']}), peak {aif.max():.1f}/s "
        f"= {checks['aif_peak_over_ref_peak']:.1f}x reference, TTP frame {aif_ttp} vs tissue {peak_ref}")
    if aif_auc <= 0:
        raise SystemExit("The AIF has a non-positive area; inspect qc_timing.png.")
    if checks["aif_peak_over_ref_peak"] < 2 or aif_ttp > peak_ref:
        log("WARNING: the selected AIF is not clearly arterial (low peak or late). rCBF/MTT/Tmax "
            "are unreliable for this series; rCBV is unaffected (it does not use the AIF).")

    dec = np.zeros_like(corr)
    dec[:, start:] = corr[:, start:]  # baseline frames are zero by construction
    aif_dec = np.zeros_like(aif)
    aif_dec[start:] = aif[start:]
    log("oSVD deconvolution...")
    cbf, tmax, svd_thr = osvd_deconvolve(dec, aif_dec, dt, args.oi_target)
    with np.errstate(divide="ignore", invalid="ignore"):
        mtt = np.where(cbf > 1e-6, (rcbv_corr / aif_auc) / cbf, 0.0)
    mtt = np.where(np.isfinite(mtt) & (mtt > 0) & (mtt < 60), mtt, 0.0)

    # timing maps, searched inside the first pass (+margin)
    ttp, tt0, cnr = bolus_timing(corr, bl, search, t, dt)

    maps = {
        "rCBV_corrected": rcbv_corr,
        "rCBV_uncorrected": rcbv_unc,
        "rCBF": cbf,
        "MTT": mtt,
        "Tmax": tmax,
        "TTP": ttp,
        "TT0": tt0,
        "K1": K1,
        "K2": K2,
        "CNR": cnr,
    }
    log("Saving maps...")
    for name, v in maps.items():
        save_map(v, mask, img, out(f"{name}.nii.gz"))
    save_map(ref_ok.astype(np.float32), mask, img, out("reference_tissue_mask.nii.gz"))
    am = np.zeros(mask.sum(), np.float32)
    am[aif_idx] = 1
    save_map(am, mask, img, out("aif_mask.nii.gz"))

    if args.smooth_fwhm_mm > 0:
        zooms = img.header.get_zooms()
        for name in ("rCBV_corrected", "rCBV_uncorrected", "rCBF"):
            vol = np.zeros(mask.shape)
            vol[mask] = maps[name]
            sm = smooth_in_mask(vol, mask, args.smooth_fwhm_mm, zooms)
            save_map(sm[mask], mask, img, out(f"{name}_smooth.nii.gz"))

    if args.save_timeseries:
        for name, v in (("deltaR2star", d), ("deltaR2star_corrected", corr)):
            vol = np.zeros(mask.shape + (T,), np.float32)
            vol[mask] = v
            nib.save(nib.Nifti1Image(vol, img.affine), out(f"{name}_timeseries.nii.gz"))

    qc_figure(out("qc_timing.png"), t, mean_sig, ss, bl_end, arrival, fp_end, win, ref, aif, d[aif_idx])
    np.savetxt(out("curves.tsv"), np.column_stack([t, mean_sig, ref, aif]), delimiter="\t",
               header="time_s\tmean_signal\tref_dR2star\taif_dR2star", comments="")

    k2_brain = K2[ref_ok]
    prov = {
        "method_version": METHOD_VERSION,
        "TE_s": args.te,
        "TR_s": args.tr,
        "n_frames": T,
        "n_voxels": int(mask.sum()),
        "baseline_frames": [int(ss), int(bl_end)],
        "baseline_mode": args.baseline,
        "bolus_arrival_frame": int(arrival),
        "first_pass_end_frame": int(fp_end),
        "integration_frames": [int(win[0]), int(win[-1]) + 1],
        "integration_mode": args.window,
        "timing_qc": tinfo,
        "leakage_correction": "BSW reference tissue (Boxerman 2006)"
        + ("" if args.no_delay_correction else " + arrival-time correction (Leigh 2012)"),
        "s0_drift_factor": g,
        "delay_median_s": float(np.median(delay)),
        "reference_tissue_voxels": int(ref_ok.sum()),
        "K2_reference_median": float(np.median(k2_brain)),
        "aif": {"method": "cluster (Mouridsen 2006)", **ainfo, **checks, "auc": aif_auc},
        "deconvolution": "oSVD block-circulant (Wu 2003)",
        "oi_target": args.oi_target,
        "svd_threshold_median": float(np.median(svd_thr)),
        "bolus_detected_fraction": float(np.mean(cnr >= 3)),
        "smooth_fwhm_mm": args.smooth_fwhm_mm,
        "units": {
            "rCBV_corrected": "arbitrary (integral of dR2*); read as a ratio to NAWM",
            "rCBV_uncorrected": "arbitrary (integral of dR2*); read as a ratio to NAWM",
            "rCBF": "1/s relative to the AIF (partial-volume scaled); read as a ratio to NAWM",
            "MTT": "s", "Tmax": "s", "TTP": "s from series start", "TT0": "s from series start",
            "K1": "dimensionless (1 = reference tissue)",
            "K2": "1/s; >0 T1-dominant leakage, <0 T2*-dominant",
            "CNR": "peak corrected dR2* / baseline dR2* SD; <3 = no detectable bolus, TTP/TT0 = 0",
        },
    }
    with open(out("KUL_dsc_fit.json"), "w") as f:
        json.dump(prov, f, indent=4)
    log(f"All maps saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
