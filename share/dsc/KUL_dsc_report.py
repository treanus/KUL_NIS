#!/usr/bin/env python3
"""ROI-based DSC perfusion measures and the DSC report.

Called by KUL_dsc_perfusion.sh after the maps have been delivered on the T1w
grid. It is the single place where perfusion numbers are produced:

  * NAWM-normalised maps (n<map>.nii.gz) for rCBV/rCBF (and their *_smooth)
  * per-ROI statistics for every map            -> *_perfusion_roi_stats.tsv
  * one row per ROI with the measures you read  -> *_perfusion_report.tsv
  * the legacy lesion-vs-NAWM tables            -> *_perfusion_stats.tsv,
                                                   *_perfusion_summary.tsv
  * a PDF report (tables, montages, fit QC)     -> REPORT/*_06_DSC_report.pdf
  * the nrCBV QC montage                        -> REPORT/*_06_DSC_rCBV.png

What is measured, and why
  * Only voxels inside the DSC slab (the fit mask brought to T1w) count: the
    maps are 0 outside it, and a lesion that extends beyond the slab would
    otherwise have its median pulled towards 0.
  * MTT, Tmax, TTP and TT0 are only meaningful where a bolus was detected
    (CNR >= 3); MTT/TTP/TT0 are also 0 where they could not be estimated. Those
    voxels are excluded from these maps' statistics, not counted as 0.
  * The headline numbers are ratios to the contralesional NAWM median. Besides
    the median, the 90th/95th percentiles are given as a reproducible
    "hotspot": histogram measures were more reproducible than hand-placed
    maximum ROIs (Wetzel 2002, PMID 12202717).
  * The % of voxels above 1.0 and above 1.75 is a threshold-based volume
    fraction in the spirit of fractional tumour burden (Hu 2012, PMID
    22561797); 1.75 is the conventional high-grade glioma rCBV cut-off (Law
    2003, PMID 14625221).
  * Statistics always use the UNSMOOTHED maps; *_smooth maps are display only.
  * A DSC voxel is ~15x larger than a T1w voxel, so an ROI covering fewer than
    10 native DSC voxels is flagged as too small to trust.

ROIs are passed as name=path. Masks on another grid are resampled (nearest
neighbour) to the maps. --roi_dir adds every NIfTI in a folder: a binary mask
becomes one ROI; an integer label map becomes one ROI per label, named from an
ITK-SNAP label description file next to it (<stem>.txt or <stem>_labels.txt)
when there is one.
"""

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import nibabel as nib
from nibabel.processing import resample_from_to

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

VERSION = "1.0"

# every map the fit can deliver, in report order
MAPS = ["rCBV_corrected", "rCBV_uncorrected", "rCBF", "MTT", "Tmax", "TTP",
        "TT0", "K1", "K2", "CNR",
        "rCBV_corrected_smooth", "rCBV_uncorrected_smooth", "rCBF_smooth"]
# maps that are read as a ratio to NAWM
NORMALISE = {"rCBV_corrected", "rCBV_uncorrected", "rCBF",
             "rCBV_corrected_smooth", "rCBV_uncorrected_smooth", "rCBF_smooth"}
# maps that need a detected bolus, and where 0 means "not estimated"
NEED_BOLUS = {"MTT", "Tmax", "TTP", "TT0"}
ZERO_INVALID = {"MTT", "TTP", "TT0"}
CNR_MIN = 3.0
MIN_NATIVE_VOXELS = 10
MIN_COVERAGE = 90.0
MIN_BOLUS = 50.0

# descriptions and contour colours for the ROIs the pipeline passes by name
KNOWN_ROIS = {
    "enhancing":    ("contrast-enhancing tumour, holes closed (includes enclosed necrosis) - HD-GLIO", "magenta", "-"),
    "nonenhancing": ("non-enhancing FLAIR abnormality - HD-GLIO", "cyan", "-"),
    "cavity":       ("resection cavity - ResSeg", "yellow", ":"),
    "lesion":       ("whole lesion mask as given (may include the cavity)", "magenta", "-"),
    "NAWM":         ("contralesional normal-appearing WM (reference)", "lime", "-"),
    "contralesional_GM": ("contralesional cortex (QC reference)", None, None),
}
REFERENCE_ROIS = {"NAWM", "contralesional_GM"}


def log(msg):
    print(msg, flush=True)


# ----------------------------------------------------------------- loading ---

def load_on(path, ref_img, order):
    """Load an image onto ref_img's grid (resampled when the grid differs)."""
    img = nib.load(path)
    if img.ndim > 3:
        img = img.slicer[..., 0]
    if img.shape[:3] != ref_img.shape[:3] or not np.allclose(img.affine, ref_img.affine, atol=1e-3):
        img = resample_from_to(img, (ref_img.shape[:3], ref_img.affine), order=order, cval=0)
    return np.asanyarray(img.dataobj).astype(np.float32)


def read_itksnap_labels(stem):
    # <stem>_itksnap.txt is what KUL_lesion_compartments.py writes
    for cand in (stem + ".txt", stem + "_itksnap.txt", stem + "_labels.txt", stem + ".label.txt"):
        if os.path.isfile(cand):
            names = {}
            pat = re.compile(r'^\s*(\d+)\s+\d+\s+\d+\s+\d+\s+[\d.]+\s+\d+\s+\d+\s+"(.*)"')
            for line in open(cand):
                m = pat.match(line)
                if m:
                    names[int(m.group(1))] = m.group(2)
            return names
    return {}


def clean_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_") or "roi"


def gather_rois(args, ref_img):
    """[(name, mask bool array, description, colour, linestyle)]"""
    rois = []

    def add(name, mask, desc, colour="orange", ls="-"):
        base, k = name, 2
        while any(r[0] == name for r in rois):
            name, k = f"{base}_{k}", k + 1
        if mask.any():
            rois.append((name, mask, desc, colour, ls))
        else:
            log(f"  ROI {name}: empty on the map grid, skipped")

    for spec in args.roi:
        name, _, path = spec.partition("=")
        if not path or not os.path.isfile(path):
            log(f"  ROI {spec}: file not found, skipped")
            continue
        desc, colour, ls = KNOWN_ROIS.get(name, (f"user ROI ({os.path.basename(path)})", "orange", "-"))
        add(name, load_on(path, ref_img, 0) > 0.5, desc, colour, ls)

    if args.roi_dir and os.path.isdir(args.roi_dir):
        for path in sorted(glob.glob(os.path.join(args.roi_dir, "*.nii*"))):
            stem = re.sub(r"\.nii(\.gz)?$", "", path)
            data = load_on(path, ref_img, 0)
            vals = np.unique(np.round(data[data > 0.5]))
            base = clean_name(os.path.basename(stem))
            if len(vals) <= 1:
                add(base, data > 0.5, f"user ROI ({os.path.basename(path)})")
            else:
                names = read_itksnap_labels(stem)
                for v in vals.astype(int):
                    label = clean_name(names.get(v, f"label{v}"))
                    add(f"{base}_{label}", np.round(data) == v,
                        f"user label map {os.path.basename(path)}, label {v}")

    if args.nawm and os.path.isfile(args.nawm):
        add("NAWM", load_on(args.nawm, ref_img, 0) > 0.5, *KNOWN_ROIS["NAWM"])
    if args.ref_gm and os.path.isfile(args.ref_gm):
        add("contralesional_GM", load_on(args.ref_gm, ref_img, 0) > 0.5, *KNOWN_ROIS["contralesional_GM"])
    return rois


# -------------------------------------------------------------- statistics ---

def describe(v):
    if v.size == 0:
        return dict(n=0, mean=np.nan, std=np.nan, p10=np.nan, p25=np.nan,
                    median=np.nan, p75=np.nan, p90=np.nan, p95=np.nan)
    p = np.percentile(v, [10, 25, 50, 75, 90, 95])
    return dict(n=int(v.size), mean=float(v.mean()), std=float(v.std()),
                p10=p[0], p25=p[1], median=p[2], p75=p[3], p90=p[4], p95=p[5])


def valid_mask(map_name, vol, coverage, cnr):
    ok = coverage & np.isfinite(vol)
    if map_name in NEED_BOLUS and cnr is not None:
        ok &= cnr >= CNR_MIN
    if map_name in ZERO_INVALID:
        ok &= vol != 0
    return ok


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{x:.{nd}g}"


# ----------------------------------------------------------------- figures ---

def to_canonical(arr, img):
    """Reorient a grid-aligned array the way nib.as_closest_canonical would."""
    ornt = nib.orientations.io_orientation(img.affine)
    return nib.orientations.apply_orientation(arr, ornt)


def pick_slices(lesion_union, coverage, n):
    zc = np.where(coverage.any(axis=(0, 1)))[0]
    zl = np.where(lesion_union.any(axis=(0, 1)))[0] if lesion_union is not None else np.array([])
    if zl.size:
        zl = zl[(zl >= zc.min()) & (zl <= zc.max())] if zc.size else zl
    if zl.size:
        lo, hi = zl.min(), zl.max()
        if hi - lo + 1 < n:
            mid = (lo + hi) // 2
            lo, hi = mid - n, mid + n
    elif zc.size:
        lo, hi = np.percentile(zc, 20), np.percentile(zc, 85)
    else:
        lo, hi = 0, coverage.shape[2] - 1
    lo, hi = int(max(lo, 0)), int(min(hi, coverage.shape[2] - 1))
    return sorted(set(np.linspace(lo, hi, n).round().astype(int)))


def bbox2d(mask2d_stack, shape, margin=6):
    xy = mask2d_stack.any(axis=2)
    xs, ys = np.where(xy)
    if xs.size == 0:
        return slice(0, shape[0]), slice(0, shape[1])
    return (slice(max(xs.min() - margin, 0), min(xs.max() + margin + 1, shape[0])),
            slice(max(ys.min() - margin, 0), min(ys.max() + margin + 1, shape[1])))


def panel(ax, under, over, cov, contours, z, crop, aspect, cmap, vmin, vmax, alpha=0.6):
    """One axial slice, radiological convention (patient right on image left)."""
    def sl(a):
        return a[crop[0], crop[1], z].T[:, ::-1]
    u = sl(under)
    lo, hi = np.percentile(under[under > 0], [1, 99.5]) if (under > 0).any() else (0, 1)
    ax.imshow(u, cmap="gray", vmin=lo, vmax=hi, origin="lower", aspect=aspect, interpolation="bilinear")
    im = None
    if over is not None:
        o = np.ma.masked_where(~sl(cov), sl(over))
        im = ax.imshow(o, cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha, origin="lower",
                       aspect=aspect, interpolation="nearest")
    for mask, colour, ls in contours:
        m = sl(mask)
        if m.any():
            ax.contour(m.astype(float), levels=[0.5], colors=[colour], linestyles=[ls],
                       linewidths=0.9, origin="lower")
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(0.02, 0.97, f"z={z}", color="w", fontsize=6, va="top", transform=ax.transAxes)
    return im


# -------------------------------------------------------------------- main ---

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perf_dir", required=True, help="RESULTS/sub-X/Perfusion")
    ap.add_argument("--participant", required=True, help="participant label, without sub-")
    ap.add_argument("--anat", required=True, help="underlay on the maps' grid")
    ap.add_argument("--coverage", help="DSC fit mask on the maps' grid")
    ap.add_argument("--nawm", help="contralesional NAWM mask (the normalisation reference)")
    ap.add_argument("--ref_gm", help="contralesional cortex mask (QC: GM/WM ratio)")
    ap.add_argument("--roi", action="append", default=[], help="name=path, repeatable")
    ap.add_argument("--roi_dir", help="folder of extra ROI masks / label maps")
    ap.add_argument("--native_ref", help="a map on the native DSC grid (for the voxel size)")
    ap.add_argument("--fit_json", help="KUL_dsc_fit.json")
    ap.add_argument("--qc_png", help="qc_timing.png from the fit")
    ap.add_argument("--report_dir", help="where the PDF and the montage go (default: no figures)")
    ap.add_argument("--thresholds", type=float, nargs="+", default=[1.0, 1.75],
                    help="nrCBV thresholds for the %%-above columns (default 1.0 1.75)")
    ap.add_argument("--n_slices", type=int, default=6)
    args = ap.parse_args()

    p = args.participant
    pre = os.path.join(args.perf_dir, f"sub-{p}_")

    # ---- maps
    maps, map_img = {}, None
    for m in MAPS:
        f = pre + m + ".nii.gz"
        if os.path.isfile(f):
            img = nib.load(f)
            if map_img is None:
                map_img = img
            maps[m] = load_on(f, map_img, 1)
    if map_img is None:
        sys.exit(f"ERROR: no perfusion maps found as {pre}<map>.nii.gz")
    zooms = map_img.header.get_zooms()[:3]
    vox_ml = float(np.prod(zooms)) / 1000.0
    native_mm3 = None
    if args.native_ref and os.path.isfile(args.native_ref):
        native_mm3 = float(np.prod(nib.load(args.native_ref).header.get_zooms()[:3]))

    if args.coverage and os.path.isfile(args.coverage):
        coverage = load_on(args.coverage, map_img, 0) > 0.5
    else:  # fall back to "any map non-zero"
        log("  no coverage mask given: using voxels where any map is non-zero")
        coverage = np.zeros(map_img.shape[:3], bool)
        for v in maps.values():
            coverage |= (v != 0) & np.isfinite(v)
    cnr = maps.get("CNR")

    rois = gather_rois(args, map_img)
    names = [r[0] for r in rois]
    log(f"  ROIs: {', '.join(names) if names else 'none'}")

    # ---- NAWM reference and normalised maps
    ref_median = {}
    nawm = next((r[1] for r in rois if r[0] == "NAWM"), None)
    if nawm is not None:
        for m in NORMALISE & set(maps):
            v = maps[m][valid_mask(m, maps[m], coverage, cnr) & nawm]
            med = float(np.median(v)) if v.size else np.nan
            if np.isfinite(med) and med > 1e-6:
                ref_median[m] = med
                out = nib.Nifti1Image((maps[m] / med).astype(np.float32), map_img.affine, map_img.header)
                out.set_data_dtype(np.float32)
                nib.save(out, pre + "n" + m + ".nii.gz")
        log(f"  NAWM reference median rCBV_corrected = {fmt(ref_median.get('rCBV_corrected'))}")
    else:
        log("  no NAWM reference: absolute (arbitrary-unit) values only, no ratios")
        # ratio maps from an earlier run no longer match any reference, and the
        # PACS export would pick them up
        for m in NORMALISE:
            stale = pre + "n" + m + ".nii.gz"
            if os.path.isfile(stale):
                os.remove(stale)
                log(f"  removed stale {os.path.basename(stale)}")

    # ---- per ROI, per map statistics
    long_rows, stats = [], {}
    for name, mask, *_ in rois:
        for m, vol in maps.items():
            ok = valid_mask(m, vol, coverage, cnr) & mask
            d = describe(vol[ok])
            stats[(name, m)] = d
            long_rows.append((m, name, "value", d))
            if m in ref_median:
                dn = describe(vol[ok] / ref_median[m])
                stats[(name, "n" + m)] = dn
                long_rows.append((m, name, "ratio_to_NAWM", dn))

    cols = ["n", "mean", "std", "p10", "p25", "median", "p75", "p90", "p95"]
    with open(pre + "perfusion_roi_stats.tsv", "w") as fh:
        fh.write("map\troi\tscale\tvalid_voxels\t" + "\t".join(cols[1:]) + "\n")
        for m, name, scale, d in long_rows:
            fh.write(f"{m}\t{name}\t{scale}\t" + "\t".join(fmt(d[c], 6) for c in cols) + "\n")

    # ---- one row per ROI: the numbers you read
    thr = args.thresholds
    report_rows = []
    for name, mask, desc, *_ in rois:
        n_roi = int(mask.sum())
        n_cov = int((mask & coverage).sum())
        cov_pct = 100.0 * n_cov / n_roi if n_roi else np.nan
        bolus_pct = (100.0 * float(((cnr >= CNR_MIN) & mask & coverage).sum()) / n_cov
                     if (cnr is not None and n_cov) else np.nan)
        native = n_cov * vox_ml * 1000.0 / native_mm3 if native_mm3 else np.nan
        row = dict(roi=name, description=desc, volume_ml=n_roi * vox_ml,
                   dsc_coverage_pct=cov_pct, native_dsc_voxels=native, bolus_detected_pct=bolus_pct)
        nr = stats.get((name, "nrCBV_corrected"))
        if nr:
            row.update(nrCBV_corr_median=nr["median"], nrCBV_corr_p90=nr["p90"], nrCBV_corr_p95=nr["p95"])
            v = maps["rCBV_corrected"][valid_mask("rCBV_corrected", maps["rCBV_corrected"], coverage, cnr) & mask]
            v = v / ref_median["rCBV_corrected"]
            for t in thr:
                row[f"pct_nrCBV_corr_gt_{t:g}"] = 100.0 * float((v > t).mean()) if v.size else np.nan
        else:
            row["rCBV_corr_median_au"] = stats.get((name, "rCBV_corrected"), {}).get("median", np.nan)
        if stats.get((name, "nrCBV_uncorrected")):
            row["nrCBV_uncorr_median"] = stats[(name, "nrCBV_uncorrected")]["median"]
        if stats.get((name, "nrCBF")):
            row["nrCBF_median"] = stats[(name, "nrCBF")]["median"]
            row["nrCBF_p90"] = stats[(name, "nrCBF")]["p90"]
        for m, key in (("MTT", "MTT_median_s"), ("Tmax", "Tmax_median_s"), ("K2", "K2_median")):
            if (name, m) in stats:
                row[key] = stats[(name, m)]["median"]
        flags = []
        if name in REFERENCE_ROIS:
            flags.append("reference")
        if np.isfinite(native) and native < MIN_NATIVE_VOXELS:
            flags.append(f"small(<{MIN_NATIVE_VOXELS} DSC voxels)")
        if np.isfinite(cov_pct) and cov_pct < MIN_COVERAGE:
            flags.append(f"partly outside DSC slab ({cov_pct:.0f}% covered)")
        if np.isfinite(bolus_pct) and bolus_pct < MIN_BOLUS and name not in ("cavity",):
            flags.append(f"bolus detected in {bolus_pct:.0f}% only")
        row["flags"] = "; ".join(flags)
        report_rows.append(row)

    keys = []
    for r in report_rows:
        keys += [k for k in r if k not in keys]
    keys = [k for k in keys if k != "flags"] + ["flags"]
    with open(pre + "perfusion_report.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in report_rows:
            fh.write("\t".join(r.get(k, "") if isinstance(r.get(k, ""), str) else fmt(r.get(k), 4)
                               for k in keys) + "\n")

    # ---- the legacy lesion-vs-NAWM tables (same columns as before)
    if "lesion" in names:
        with open(pre + "perfusion_stats.tsv", "w") as fh:
            fh.write("map\troi\tmean\tmedian\tstd\tvoxels\n")
            for m in maps:
                for name, tag in (("lesion", "lesion"), ("NAWM", "nawm")):
                    if (name, m) in stats:
                        d = stats[(name, m)]
                        fh.write(f"{m}\t{tag}\t{fmt(d['mean'], 6)}\t{fmt(d['median'], 6)}\t"
                                 f"{fmt(d['std'], 6)}\t{d['n']}\n")
        with open(pre + "perfusion_summary.tsv", "w") as fh:
            fh.write("map\tlesion_median\tnawm_median\tnormalised_ratio\n")
            for m in maps:
                les = stats[("lesion", m)]["median"]
                ref = stats.get(("NAWM", m), {}).get("median", np.nan)
                ratio = les / ref if (np.isfinite(ref) and ref > 1e-6) else np.nan
                fh.write(f"{m}\t{fmt(les, 6)}\t{fmt(ref, 6)}\t{fmt(ratio, 4)}\n")

    # ---- QC numbers
    fit = {}
    if args.fit_json and os.path.isfile(args.fit_json):
        fit = json.load(open(args.fit_json))
    qc, warnings = [], []
    if fit:
        bl = fit.get("baseline_frames", [None, None])
        nbl = (bl[1] - bl[0] + 1) if None not in bl else None
        qc += [("frames / baseline frames", f"{fit.get('n_frames')} / {nbl} (frames {bl[0]}-{bl[1]})"),
               ("bolus arrival / first-pass end (frame)", f"{fit.get('bolus_arrival_frame')} / {fit.get('first_pass_end_frame')}"),
               ("AIF voxels, peak AIF / peak tissue", f"{fit.get('aif', {}).get('n_voxels')}, {fmt(fit.get('aif', {}).get('aif_peak_over_ref_peak'), 3)}"),
               ("voxels with a detectable bolus", f"{100 * fit.get('bolus_detected_fraction', np.nan):.0f}%"),
               ("baseline drift factor (S0 levelling)", fmt(fit.get("s0_drift_factor"), 4)),
               ("K2 of the reference tissue (should be ~0)", fmt(fit.get("K2_reference_median"), 3))]
        if nbl is not None and nbl < 30:
            warnings.append(f"only {nbl} baseline frames (2020 consensus: 30-50; Boxerman 2020, PMID 32516388)")
        if fit.get("n_frames", 999) < 120:
            warnings.append(f"{fit.get('n_frames')} frames: shorter than the 2020 consensus (>=120), the tail is truncated")
        if fit.get("aif", {}).get("aif_peak_over_ref_peak", 99) < 3:
            warnings.append("AIF peak < 3x tissue peak: partial-volumed AIF, absolute rCBF/MTT unreliable")
    gm = stats.get(("contralesional_GM", "rCBV_corrected"), {}).get("median", np.nan)
    wm = stats.get(("NAWM", "rCBV_corrected"), {}).get("median", np.nan)
    if np.isfinite(gm) and np.isfinite(wm) and wm > 0:
        gmwm = gm / wm
        qc.append(("contralesional cortex / NAWM rCBV", f"{gmwm:.2f}  (cortex has roughly twice the blood volume of WM; "
                   "partial volume with 5 mm slices lowers this)"))
        if gmwm < 1.3:
            warnings.append(f"cortex/NAWM rCBV only {gmwm:.2f}: poor GM/WM contrast - check registration, AIF and the NAWM mask")
    if nawm is not None:
        qc.append(("NAWM reference voxels (in DSC slab)", str(int((nawm & coverage).sum()))))

    # ---- log the headline table
    log("")
    w = max([22] + [len(r["roi"]) + 2 for r in report_rows])
    if ref_median:
        log(f"  {'ROI':<{w}}{'vol mL':>8}{'nrCBV med':>11}{'p90':>7}{'%>1.75':>8}{'nrCBF med':>11}  flags")
        for r in report_rows:
            log(f"  {r['roi']:<{w}}{fmt(r['volume_ml'], 3):>8}{fmt(r.get('nrCBV_corr_median'), 3):>11}"
                f"{fmt(r.get('nrCBV_corr_p90'), 3):>7}{fmt(r.get('pct_nrCBV_corr_gt_1.75'), 3):>8}"
                f"{fmt(r.get('nrCBF_median'), 3):>11}  {r['flags']}")
    else:
        log(f"  {'ROI':<{w}}{'vol mL':>8}{'rCBV a.u. med':>15}  flags")
        for r in report_rows:
            log(f"  {r['roi']:<{w}}{fmt(r['volume_ml'], 3):>8}{fmt(r.get('rCBV_corr_median_au'), 4):>15}  {r['flags']}")
    for w in warnings:
        log(f"  WARNING: {w}")

    if not args.report_dir:
        return
    os.makedirs(args.report_dir, exist_ok=True)

    # ---- figures (canonical RAS, shown radiologically)
    anat = to_canonical(load_on(args.anat, map_img, 1), map_img)
    cov_c = to_canonical(coverage, map_img)
    ornt = nib.orientations.io_orientation(map_img.affine)
    cz = [zooms[int(i)] for i in ornt[:, 0].argsort()]  # voxel size along canonical x, y, z
    aspect = cz[1] / cz[0]

    target = [r for r in rois if r[0] not in REFERENCE_ROIS]
    compartments = [r for r in target if r[0] != "lesion"]
    drawn = compartments if compartments else target
    contours = [(to_canonical(r[1], map_img), r[3], r[4]) for r in drawn]
    contours += [(to_canonical(r[1], map_img), r[3], r[4]) for r in rois if r[0] == "NAWM"]
    union = np.zeros(cov_c.shape, bool)
    for r in target:
        union |= to_canonical(r[1], map_img)
    slices = pick_slices(union if union.any() else None, cov_c, args.n_slices)
    crop = bbox2d(cov_c, cov_c.shape)
    legend = [Line2D([0], [0], color=r[3], ls=r[4], lw=1.5, label=r[0]) for r in drawn]
    if nawm is not None:
        legend.append(Line2D([0], [0], color="lime", lw=1.5, label="NAWM (reference)"))

    def disp(name):
        for m in (name + "_smooth", name):
            if m in maps:
                return m
        return None

    rows = []
    if ref_median:
        m = disp("rCBV_corrected")
        rows.append((f"nrCBV corrected ({'3D smoothed, display' if m.endswith('smooth') else 'unsmoothed'})",
                     to_canonical(maps[m] / ref_median[m], map_img), "turbo", 0, 4))
        m = disp("rCBF")
        if m in ref_median:
            rows.append(("nrCBF", to_canonical(maps[m] / ref_median[m], map_img), "turbo", 0, 4))
    elif "rCBV_corrected" in maps:
        v = maps["rCBV_corrected"][coverage]
        rows.append(("rCBV corrected (a.u., no NAWM reference)", to_canonical(maps["rCBV_corrected"], map_img),
                     "turbo", 0, float(np.percentile(v, 99)) if v.size else 1))
    if "K2" in maps:
        v = np.abs(maps["K2"][coverage])
        k = float(np.percentile(v, 98)) if v.size else 1
        rows.append(("K2 leakage (>0 T1-dominant, <0 T2*-dominant)", to_canonical(maps["K2"], map_img),
                     "RdBu_r", -k, k))

    def montage(fig, row_defs, n_rows_per_map=1):
        ncol = len(slices) if n_rows_per_map == 1 else int(np.ceil(len(slices) / n_rows_per_map))
        nrow = len(row_defs) * n_rows_per_map
        gs = fig.add_gridspec(nrow, ncol + 1, width_ratios=[1] * ncol + [0.05],
                              left=0.02, right=0.95, top=0.90, bottom=0.06, wspace=0.03, hspace=0.12)
        for i, (title, vol, cmap, vmin, vmax) in enumerate(row_defs):
            im = None
            for j, z in enumerate(slices):
                rr, cc = i * n_rows_per_map + j // ncol, j % ncol
                ax = fig.add_subplot(gs[rr, cc])
                im = panel(ax, anat, vol, cov_c, contours, z, crop, aspect, cmap, vmin, vmax) or im
                if cc == 0 and j // ncol == 0:
                    ax.set_title(title, fontsize=8, loc="left")
                if rr == 0 and cc == 0:
                    ax.text(0.02, 0.03, "R", color="w", fontsize=7, transform=ax.transAxes)
                    ax.text(0.95, 0.03, "L", color="w", fontsize=7, transform=ax.transAxes)
            if im is not None:
                cax = fig.add_subplot(gs[i * n_rows_per_map:(i + 1) * n_rows_per_map, -1])
                fig.colorbar(im, cax=cax).ax.tick_params(labelsize=6)
        if legend:
            fig.legend(handles=legend, loc="lower center", ncol=len(legend), fontsize=7, frameon=False)

    # the REPORT montage (replaces the mrview render)
    if rows:
        fig = plt.figure(figsize=(12, 5.2), facecolor="white")
        fig.suptitle(f"sub-{p}  DSC {rows[0][0]}  - ratio to contralesional NAWM, 1.75 = high-grade cut-off"
                     if ref_median else f"sub-{p}  DSC {rows[0][0]}", fontsize=10)
        slices_saved = slices
        slices = pick_slices(union if union.any() else None, cov_c, 2 * args.n_slices)
        montage(fig, rows[:1], n_rows_per_map=2)
        slices = slices_saved
        fig.savefig(os.path.join(args.report_dir, f"sub-{p}_06_DSC_rCBV.png"), dpi=150)
        plt.close(fig)

    pdf_path = os.path.join(args.report_dir, f"sub-{p}_06_DSC_report.pdf")
    with PdfPages(pdf_path) as pdf:
        # page 1: the numbers
        fig = plt.figure(figsize=(11.69, 8.27), facecolor="white")
        fig.text(0.04, 0.95, f"DSC perfusion report - sub-{p}", fontsize=15, weight="bold")
        method = (f"KUL_dsc_perfusion method v{fit.get('method_version', '?')}: BSW reference-tissue leakage correction "
                  f"+ arrival-time correction, oSVD, cluster AIF. Values are ratios to the contralesional NAWM median "
                  f"(unsmoothed maps; voxels outside the DSC slab excluded). Report v{VERSION}.")
        fig.text(0.04, 0.915, method, fontsize=7.5, wrap=True)
        if ref_median:
            head = ["ROI", "vol\nmL", "DSC\ncover %", "bolus\n%", "nrCBV\nmedian", "nrCBV\np90", "nrCBV\np95"] + \
                   [f"% nrCBV\n> {t:g}" for t in thr] + ["nrCBF\nmedian", "nrCBF\np90", "MTT\ns", "K2\nmedian", "flags"]
            keys_t = ["roi", "volume_ml", "dsc_coverage_pct", "bolus_detected_pct", "nrCBV_corr_median",
                      "nrCBV_corr_p90", "nrCBV_corr_p95"] + [f"pct_nrCBV_corr_gt_{t:g}" for t in thr] + \
                     ["nrCBF_median", "nrCBF_p90", "MTT_median_s", "K2_median", "flags"]
        else:
            head = ["ROI", "vol\nmL", "DSC\ncover %", "bolus\n%", "rCBV a.u.\nmedian", "MTT\ns", "K2\nmedian", "flags"]
            keys_t = ["roi", "volume_ml", "dsc_coverage_pct", "bolus_detected_pct", "rCBV_corr_median_au",
                      "MTT_median_s", "K2_median", "flags"]
        cell = []
        for r in report_rows:
            line = []
            for k in keys_t:
                v = r.get(k, "")
                if k in ("roi", "flags"):
                    line.append(v)
                elif k in ("volume_ml",):
                    line.append(fmt(v, 3))
                elif k.startswith("pct") or k.endswith("_pct"):
                    line.append("n/a" if not np.isfinite(v) else f"{v:.0f}")
                elif k == "K2_median":
                    line.append(fmt(v, 2))
                else:
                    line.append("n/a" if v == "" or not np.isfinite(v) else f"{v:.2f}")
            cell.append(line)
        ax = fig.add_axes([0.03, 0.47, 0.94, 0.40]); ax.axis("off")
        if cell:
            ncols = len(head)
            widths = [0.12] + [(1 - 0.12 - 0.20) / (ncols - 2)] * (ncols - 2) + [0.20]
            tb = ax.table(cellText=cell, colLabels=head, colWidths=widths, loc="upper center", cellLoc="center")
            tb.auto_set_font_size(False); tb.set_fontsize(7.5); tb.scale(1, 1.9)
            for (i, j), c in tb.get_celld().items():
                if i == 0:
                    c.set_facecolor("#e8e8e8"); c.set_height(c.get_height() * 1.4)
                elif j in (0, ncols - 1):
                    c.get_text().set_horizontalalignment("left"); c.PAD = 0.03
                if i > 0 and cell[i - 1][0] in REFERENCE_ROIS:
                    c.set_facecolor("#f4f4f4")
        else:
            ax.text(0.5, 0.8, "no ROIs available (no lesion mask, no NAWM)", ha="center")

        y = 0.44
        fig.text(0.04, y, "ROIs", fontsize=9, weight="bold"); y -= 0.022
        for r in report_rows:
            fig.text(0.05, y, f"{r['roi']}: {r['description']}", fontsize=7); y -= 0.018
        y -= 0.01
        fig.text(0.04, y, "Acquisition / fit QC", fontsize=9, weight="bold"); y -= 0.022
        for k, v in qc:
            fig.text(0.05, y, f"{k}: {v}", fontsize=7); y -= 0.018
        for w in warnings:
            fig.text(0.05, y, f"WARNING: {w}", fontsize=7, color="#b00000"); y -= 0.018
        notes = ("Reading: nrCBV median describes the whole ROI; p90/p95 are reproducible 'hotspot' measures (Wetzel 2002, "
                 "PMID 12202717); % > 1.75 is the volume fraction above the conventional high-grade cut-off (Law 2003, "
                 "PMID 14625221), cf. fractional tumour burden (Hu 2012, PMID 22561797). The enhancing ROI includes "
                 "enclosed necrosis, which lowers its median - read the percentiles. The NAWM row's % > cut-off is the "
                 "background rate from voxel noise; a lesion value is meaningful only well above it. Cut-offs depend on protocol and "
                 "software; compare within a patient over time with the same protocol. Research output, not a "
                 "diagnostic device.")
        fig.text(0.04, 0.035, notes, fontsize=6.5, wrap=True, color="#333333")
        pdf.savefig(fig); plt.close(fig)

        # page 2: montages
        if rows:
            fig = plt.figure(figsize=(11.69, 8.27), facecolor="white")
            fig.suptitle(f"sub-{p}  perfusion maps through the lesion (radiological: R on the left)", fontsize=11)
            montage(fig, rows)
            pdf.savefig(fig); plt.close(fig)

        # page 3: the fit QC figure
        if args.qc_png and os.path.isfile(args.qc_png):
            fig = plt.figure(figsize=(11.69, 8.27), facecolor="white")
            ax = fig.add_axes([0.02, 0.06, 0.96, 0.86]); ax.axis("off")
            ax.imshow(plt.imread(args.qc_png))
            fig.suptitle(f"sub-{p}  fit QC: baseline, bolus arrival, AIF and reference-tissue curves", fontsize=11)
            pdf.savefig(fig); plt.close(fig)

    log(f"  report: {pdf_path}")


if __name__ == "__main__":
    main()
