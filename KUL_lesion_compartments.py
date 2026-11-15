#!/usr/bin/env python3
#
# KUL_lesion_compartments.py, 2026-10-20
# Split a binary lesion mask into intensity-based compartments.
#
# The IMAGE intensities inside the (dilated) lesion mask are clustered, not the
# mask itself: every available channel on the T1 grid (T1, FLAIR, T2, cT1-T1
# subtraction, SWI, ADC) is z-scored against normal-appearing white matter
# (FreeSurfer aseg WM away from the lesion), a Gaussian mixture is fitted with
# k chosen from the BIC curve, and the classes are named from their intensity
# profile with a lesion-type preset. Spatial regularisation is deliberately
# very light (sub-voxel posterior smoothing): only voxels whose own evidence is
# ambiguous change class; confident voxels, even isolated ones, never do.
#
# The result is a SUGGESTION: a label image + ITK-SNAP label description file
# (rename/edit there), a QA png and a JSON with every number used.
#
import argparse
import glob
import json
import os
import sys

import numpy as np
import nibabel as nib
from nibabel.processing import resample_from_to
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
from sklearn.mixture import GaussianMixture

VERSION = "0.2"

# channel name -> path relative to the study dir ({p} = participant)
CHANNELS = [
    ("T1",     "RESULTS/sub-{p}/Anat/T1w.nii.gz"),
    ("FLAIR",  "RESULTS/sub-{p}/Anat/FLAIR_reg2_T1w.nii.gz"),
    ("T2",     "RESULTS/sub-{p}/Anat/T2w_reg2_T1w.nii.gz"),
    ("cT1sub", "RESULTS/sub-{p}/Anat/cT1w_T1w_subtracted.nii.gz"),
    ("SWI",    "RESULTS/sub-{p}/Anat/SWI_reg2_T1w.nii.gz"),
    ("ADC",    "dwiprep/sub-{p}/sub-{p}/qa/adc_reg2T1w.nii.gz"),
]
# original acquisitions, only to report the native voxel size of each channel
NATIVE = {
    "T1": "sub-{p}_T1w.nii.gz", "FLAIR": "sub-{p}_FLAIR.nii.gz",
    "T2": "sub-{p}_T2w.nii.gz", "SWI": "sub-{p}_SWI.nii.gz",
    "cT1sub": "sub-{p}_ce-gadolinium_T1w.nii.gz",
}
# first existing file wins
MASKS = [
    "RESULTS/sub-{p}/Lesion/sub-{p}_lesion_and_cavity.nii.gz",
    "RESULTS/sub-{p}/Lesion/sub-{p}_hdglio_lesion_total.nii.gz",
    "RESULTS/sub-{p}/Lesion/lesion.nii.gz",
    "RESULTS/sub-{p}/Anat/sub-{p}_lesion.nii.gz",
]
# reference segmentations used for validation when present
REFS = [
    ("hdglio_solid",        "RESULTS/sub-{p}/Lesion/sub-{p}_hdglio_lesion_solid_tissue.nii.gz"),
    ("hdglio_perilesional", "RESULTS/sub-{p}/Lesion/sub-{p}_hdglio_lesion_perilesional_tissue.nii.gz"),
    ("resseg_cavity",       "RESULTS/sub-{p}/Lesion/sub-{p}_resseg_cavity_only.nii.gz"),
]
ASEG = "BIDS/derivatives/freesurfer/sub-{p}/mri/aseg.mgz"
# types 1/2: adopted as the compartments instead of clustering (earlier wins on overlap)
HDGLIO_CLASSES = [
    ("enhancing_hdglio",          "RESULTS/sub-{p}/Lesion/sub-{p}_hdglio_lesion_solid_tissue.nii.gz"),
    ("cavity_resseg",             "RESULTS/sub-{p}/Lesion/sub-{p}_resseg_cavity_only.nii.gz"),
    ("nonenhancing_FLAIR_hdglio", "RESULTS/sub-{p}/Lesion/sub-{p}_hdglio_lesion_perilesional_tissue.nii.gz"),
]
# a class is accepted when its split-half Dice, coherence and volume clear these
QC = {"repro": 0.7, "coherence": 0.5, "min_ml": 0.05, "min_frac": 0.01}
NON_LESION = ("margin_partial_volume", "normal_appearing")

ASEG_WM = (2, 41)
ASEG_CORTEX = (3, 42)
# ventricles, CSF, choroid plexus: the dilation ring is not grown into these
ASEG_CSF = (4, 5, 14, 15, 24, 31, 43, 44, 63, 72)
# channels where 0 is a legitimate value (all others: 0 = outside FOV/mask)
ZERO_OK = ("cT1sub",)

# name priority per preset (also fixes label numbering, so colours follow names)
PRESETS = {
    "glioma": ["enhancing", "cyst_cavity_necrosis", "nonenhancing_tumour_like",
               "FLAIR_abnormal", "oedema_like", "margin_partial_volume", "normal_appearing", "unassigned"],
    "cavernoma": ["blood_products", "hemosiderin", "core_blood_products", "susceptibility_blooming",
                  "oedema_gliosis", "CSF_like", "margin_partial_volume", "normal_appearing", "unassigned"],
    "fcd": ["FLAIR_hyperintense", "T1_abnormal", "CSF_like", "margin_partial_volume", "normal_appearing", "unassigned"],
    "generic": ["class", "CSF_like", "margin_partial_volume", "normal_appearing"],
}

# categorical slots 1-8 (dark-surface steps; validated on #1a1a19 and on black)
COLOURS = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181",
           "#008300", "#9085e9", "#e66767"]
SURFACE, INK, INK2, GRID = "#1a1a19", "#ffffff", "#c3c2b7", "#3a3a38"

Z_ABN = 1.5   # |median z| above which a channel counts as abnormal in a class


def log(msg):
    print(msg, flush=True)


def load_on_grid(path, ref, order):
    """Load image data on the grid of ref (resample only when needed)."""
    img = nib.load(path)
    if img.ndim > 3:
        img = img.slicer[..., 0]
    if img.shape[:3] == ref.shape[:3] and np.allclose(img.affine, ref.affine, atol=1e-3):
        return np.asanyarray(img.dataobj).astype(np.float32), False
    out = resample_from_to(img, (ref.shape[:3], ref.affine), order=order, cval=0)
    return np.asanyarray(out.dataobj).astype(np.float32), True


def robust_z(x, ref_vals):
    med = float(np.median(ref_vals))
    mad = float(np.median(np.abs(ref_vals - med))) * 1.4826
    if mad <= 0:
        mad = float(np.std(ref_vals)) or 1.0
    return (x - med) / mad, med, mad


def descriptor(med):
    parts = []
    for ch, z in med.items():
        if z > 4:
            parts.append(ch + "++")
        elif z > Z_ABN:
            parts.append(ch + "+")
        elif z < -4:
            parts.append(ch + "--")
        elif z < -Z_ABN:
            parts.append(ch + "-")
    return " ".join(parts) if parts else "all channels ~NAWM"


def name_classes(meds, frac_ring, preset):
    """Name each class from its median-z profile. Returns list of names."""
    k = len(meds)
    names = [None] * k
    g = lambda c, ch, default=0.0: meds[c].get(ch, default)
    has = lambda ch: ch in meds[0]

    def normal(c):
        # one mildly abnormal channel (typically the coarse, resampled ADC
        # smearing CSF/oedema into its neighbours) does not make tissue abnormal
        z = sorted((abs(v) for v in meds[c].values()), reverse=True)
        return z[0] < 3 and (len(z) < 2 or z[1] < Z_ABN)

    # an abnormal class that lives mostly in the dilation ring is the lesion
    # border: partial volume / transition, not a tissue type of its own
    for c in range(k):
        if frac_ring[c] >= 0.5 and not normal(c):
            names[c] = "margin_partial_volume"
    # fluid inside the mask (sulci, enlarged spaces); glioma has its own cyst/cavity rule
    if preset != "glioma":
        for c in range(k):
            # (blood products are T1/FLAIR-dark with high ADC too, but CSF is never SWI-dark)
            if names[c] is None and g(c, "T1") < -2 and g(c, "FLAIR", -9) < -Z_ABN and g(c, "SWI") > -3 and \
                    (g(c, "ADC") > 3 or g(c, "T2") > 3 or not (has("ADC") or has("T2"))):
                names[c] = "CSF_like"

    if preset == "glioma":
        if has("cT1sub"):
            c = max(range(k), key=lambda i: g(i, "cT1sub"))
            if g(c, "cT1sub") > 3:
                names[c] = "enhancing"
        for c in range(k):
            if names[c] is None and g(c, "T1") < -2 and g(c, "FLAIR") < Z_ABN and \
                    (g(c, "T2") > 2 or g(c, "ADC") > 2 or not (has("T2") or has("ADC"))):
                names[c] = "cyst_cavity_necrosis"
        for c in range(k):
            if names[c] is None and normal(c):
                names[c] = "normal_appearing"
        key = "FLAIR" if has("FLAIR") else "T2"
        rest = [c for c in range(k) if names[c] is None and g(c, key) > Z_ABN and frac_ring[c] < 0.5]
        if len(rest) >= 2 and has("ADC"):
            rest.sort(key=lambda c: g(c, "ADC"))
            names[rest[0]] = "nonenhancing_tumour_like"
            names[rest[-1]] = "oedema_like"
            for c in rest[1:-1]:
                names[c] = "FLAIR_abnormal"
        else:
            for c in rest:
                names[c] = "FLAIR_abnormal"
    elif preset == "cavernoma":
        dark = "SWI" if has("SWI") else "T2"
        n_dark = sum(g(c, dark) < -3 for c in range(k))
        for c in range(k):
            if g(c, dark) < -3 and n_dark == 1:
                # one dark class = the whole cavernoma (hemosiderin + core together)
                names[c] = "blood_products"
            elif g(c, dark) < -3:
                if dark == "T2" or g(c, "T2") < 0:
                    names[c] = "hemosiderin"
                elif frac_ring[c] >= 0.5:
                    # SWI-dark but T2-bright, outside the drawn lesion: blooming
                    names[c] = "susceptibility_blooming"
                else:
                    names[c] = "core_blood_products"
        for c in range(k):
            if names[c] is None and \
                    (g(c, "T1") > Z_ABN or g(c, "T2") > 2 or g(c, dark) < -Z_ABN):
                names[c] = "core_blood_products"
        for c in range(k):
            if names[c] is None and g(c, "FLAIR") > Z_ABN:
                names[c] = "oedema_gliosis"
        for c in range(k):
            if names[c] is None and normal(c):
                names[c] = "normal_appearing"
    elif preset == "fcd":
        for c in range(k):
            if g(c, "FLAIR") > Z_ABN:
                names[c] = "FLAIR_hyperintense"
            elif abs(g(c, "T1")) > Z_ABN:
                names[c] = "T1_abnormal"
            elif normal(c):
                names[c] = "normal_appearing"
    else:
        for c in range(k):
            if names[c] is None and normal(c):
                names[c] = "normal_appearing"
        names = [n or "class" for n in names]
    return [n or "unassigned" for n in names]


def split_half_dice(X, fit, k, full_lab, gmm, seed, nrep=2):
    """Refit on two disjoint halves of the fitting voxels, label every voxel
    with each half-model, match classes to the full model (Hungarian) and
    return the per-class Dice between the two halves, shape (k, nrep)."""
    rng = np.random.default_rng(seed + 1000 * k)
    out = np.zeros((k, nrep))
    for r in range(nrep):
        perm = rng.permutation(len(fit))
        h = len(fit) // 2
        labs = []
        for half in (perm[:h], perm[h:2 * h]):
            lab = gmm(k, fit[half]).predict(X)
            conf = np.zeros((k, k))
            np.add.at(conf, (full_lab, lab), 1)
            ri, ci = linear_sum_assignment(-conf)
            m = np.empty(k, int)
            m[ci] = ri
            labs.append(m[lab])
        for c in range(k):
            a, b = labs[0] == c, labs[1] == c
            n = a.sum() + b.sum()
            out[c, r] = 2 * (a & b).sum() / n if n else 0.0
    return out


def bhattacharyya(gm):
    """Pairwise Bhattacharyya coefficient (0 = separate, 1 = identical) of the mixture components."""
    k = gm.n_components
    bc = np.eye(k)
    for i in range(k):
        for j in range(i + 1, k):
            ci, cj = gm.covariances_[i], gm.covariances_[j]
            s = (ci + cj) / 2
            d = gm.means_[i] - gm.means_[j]
            db = d @ np.linalg.solve(s, d) / 8 + 0.5 * (np.linalg.slogdet(s)[1] - 0.5 * (
                np.linalg.slogdet(ci)[1] + np.linalg.slogdet(cj)[1]))
            bc[i, j] = bc[j, i] = np.exp(-db)
    return bc


def coherence(lab, zone, lv):
    """Mean fraction of a class voxel's 26 in-zone neighbours that share its class
    (salt-and-pepper ~0.3, thin shell ~0.6, compact region > 0.8)."""
    m = (lab == lv).astype(np.float32)
    same = ndi.uniform_filter(m, 3) * 27 - 1
    nb = ndi.uniform_filter(zone.astype(np.float32), 3) * 27 - 1
    sel = lab == lv
    return float(np.mean(same[sel] / np.maximum(nb[sel], 1))) if sel.any() else 0.0


def unique_names(names):
    out, seen = [], {}
    for n in names:
        if names.count(n) > 1:
            seen[n] = seen.get(n, 0) + 1
            out.append(f"{n}_{seen[n]}")
        else:
            out.append(n)
    return out


def to_ras(a, affine):
    ornt = nib.orientations.io_orientation(affine)
    return nib.orientations.apply_orientation(a, ornt)


def qa_figure(path, title, chans, ras, zooms_ras, mask_ras, lab_ras, comps_ras,
              classes, bic, k, colours):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, to_rgba

    plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                         "savefig.facecolor": SURFACE, "text.color": INK,
                         "axes.labelcolor": INK2, "axes.edgecolor": GRID,
                         "xtick.color": INK2, "ytick.color": INK2, "font.size": 9})
    comp_ids = comps_ras["ids"][:4]
    ncol = len(chans) + 2
    nrow = len(comp_ids) + 1
    fig = plt.figure(figsize=(2.3 * ncol, 2.5 * nrow + 0.6))
    gs = fig.add_gridspec(nrow, ncol, height_ratios=[1] * len(comp_ids) + [1.25])
    cmap = ListedColormap([(0, 0, 0, 0)] + [to_rgba(c, 0.55) for c in colours])
    margin = np.ceil(15 / np.array(zooms_ras)).astype(int)

    def crop(a, bb, ax2):
        sl = [slice(max(bb[d][0] - margin[d], 0), min(bb[d][1] + margin[d], a.shape[d]))
              for d in range(3)]
        sl[ax2] = bb["slice"][ax2]
        return a[tuple(sl)]

    for r, cid in enumerate(comp_ids):
        comp = comps_ras["lab"] == cid
        idx = np.argwhere(comp)
        bb = {d: (idx[:, d].min(), idx[:, d].max() + 1) for d in range(3)}
        bb["slice"] = np.round(idx.mean(0)).astype(int)
        vol = comp.sum() * np.prod(zooms_ras) / 1000
        for j, (ch, vol_ch) in enumerate(chans + [("classes, axial", None), ("classes, coronal", None)]):
            ax = fig.add_subplot(gs[r, j])
            ax2 = 1 if ch == "classes, coronal" else 2
            base = ras["FLAIR" if "FLAIR" in ras else "T1"] if vol_ch is None else vol_ch
            img = crop(base, bb, ax2)
            lo, hi = np.percentile(img[img != 0], [1, 99]) if np.any(img != 0) else (0, 1)
            # radiological display: patient right on screen left
            show = lambda a: a.T[::-1, ::-1] if ax2 == 2 else a.T[::-1, ::-1]
            asp = (zooms_ras[1] / zooms_ras[0]) if ax2 == 2 else (zooms_ras[2] / zooms_ras[0])
            ax.imshow(show(img), cmap="gray", vmin=lo, vmax=hi, aspect=asp, interpolation="nearest")
            if vol_ch is None:
                ax.imshow(show(crop(lab_ras, bb, ax2)), cmap=cmap, vmin=0, vmax=len(colours),
                          aspect=asp, interpolation="nearest")
                ax.contour(show(crop(mask_ras.astype(float), bb, ax2)), [0.5],
                           colors=INK, linewidths=0.7)
            ax.set_xticks([]), ax.set_yticks([])
            ax.text(0.03, 0.97, "R", transform=ax.transAxes, va="top", color=INK2, fontsize=7)
            if r == 0:
                ax.set_title(ch, color=INK, fontsize=9)
            if j == 0:
                ax.set_ylabel(f"component {cid}\n{vol:.2f} mL", color=INK2)

    # class profiles: median z per channel, IQR whiskers
    ax = fig.add_subplot(gs[-1, : ncol - 2])
    ax_prof, direct = ax, []
    chn = [c for c, _ in chans]
    x = np.arange(len(chn))
    for i, cl in enumerate(classes):
        off = (i - (len(classes) - 1) / 2) * 0.08
        med = [cl["median_z"][c] for c in chn]
        q1 = [cl["iqr_z"][c][0] for c in chn]
        q3 = [cl["iqr_z"][c][1] for c in chn]
        ax.errorbar(x + off, med, yerr=[np.subtract(med, q1), np.subtract(q3, med)],
                    color=colours[i], lw=2, marker="o", ms=8, mec=SURFACE, mew=2,
                    elinewidth=1, capsize=0,
                    label=f"{cl['label']}  {cl['name']}  ({cl['vol_mask_ml']:.2f} + {cl['vol_ring_ml']:.2f} mL)")
        direct.append((x[-1] + off, med[-1], cl["name"]))
    ax.axhline(0, color=INK2, lw=0.8, ls="--")
    ax.set_yscale("symlog", linthresh=2)
    ax.set_xticks(x, chn)
    ax.set_xlim(-0.4, len(chn) - 0.4 + 0.9)
    ax.set_ylabel("median z vs NAWM (IQR)")
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.12), ncol=2, frameon=False,
              fontsize=8, title="label  name  (in mask + ring)", title_fontsize=8, labelcolor=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # class checks (cluster mode) or the source of the labels (HD-GLIO mode)
    ax = fig.add_subplot(gs[-1, ncol - 2:])
    ax.axis("off")
    if bic is None:
        ax.text(0, 0.95, "labels adopted from HD-GLIO / resseg\nnot clustered", color=INK,
                fontsize=8, va="top", transform=ax.transAxes)
    else:
        ax.set_title(f"class checks, k={k}", color=INK, fontsize=9, loc="left")
        lines = [f"{'':3s}{'split-half':>11s}{'coher-':>8s}{'vol':>8s}", f"{'':3s}{'Dice':>11s}{'ence':>8s}{'mL':>8s}"]
        for cl in classes:
            lines.append(f"{cl['label']:<3d}{cl['reproducibility_dice']:>11.2f}{cl['coherence']:>8.2f}"
                         f"{cl['vol_ml']:>8.2f}" + ("  " + ", ".join(cl["checks_failed"]) if cl["checks_failed"] else "  ok"))
        lines.append("")
        lines.append(f"pass: Dice >= {QC['repro']}, coherence >= {QC['coherence']}")
        lines.append(f"      vol >= max({QC['min_ml']} mL, {QC['min_frac']:.0%} of zone)")
        ax.text(0, 0.88, "\n".join(lines), color=INK, fontsize=8, va="top", family="monospace",
                transform=ax.transAxes)

    fig.suptitle(title, color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    # direct labels at the line ends, pushed apart so they never overlap
    fig.canvas.draw()
    tr = ax_prof.transData
    pts = sorted(((tr.transform((xx, yy))[1], xx, yy, nm) for xx, yy, nm in direct))
    gap = 11 * fig.dpi / 72
    placed = []
    for ydisp, xx, yy, nm in pts:
        placed.append(max(ydisp, placed[-1] + gap) if placed else ydisp)
    # keep the stack inside the axes
    top = ax_prof.get_window_extent().y1 - gap / 2
    shift = max(0.0, placed[-1] - top) if placed else 0.0
    for (ydisp, xx, yy, nm), y_new in zip(pts, placed):
        ax_prof.annotate(nm, (xx, yy), xytext=(9, (y_new - shift - ydisp) * 72 / fig.dpi),
                         textcoords="offset points", color=INK2, fontsize=7, va="center")
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Split a lesion mask into intensity-based compartments (suggestion + QA).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (run from the study folder):
  KUL_lesion_compartments.py -p JohnDoe --preset glioma           # type 1/2: adopts HD-GLIO
  KUL_lesion_compartments.py -p JohnDoe --preset glioma --cluster
  KUL_lesion_compartments.py -p JohnDoe --preset cavernoma --exclude ADC
  KUL_lesion_compartments.py -p JohnDoe --preset generic --k 3 --mask my_mask.nii.gz
  KUL_lesion_compartments.py --mask m.nii.gz --channel T1=t1.nii.gz --channel FLAIR=fl.nii.gz \\
                             --aseg aseg.mgz -o outdir --prefix test

Types 1/2: when HD-GLIO compartments exist in RESULTS/sub-<p>/Lesion/ they are
adopted as-is (enhancing > resseg cavity > non-enhancing FLAIR) and nothing is
clustered, unless --cluster is given.

Choosing k: every k in --ks (default 3,4) is fitted and written. Each class is
checked: split-half reproducibility (refit on two disjoint halves of the voxels,
Dice between them), spatial coherence (fraction of neighbours sharing the class)
and volume. The largest k whose classes all pass is the main output; the other
k is written with a _k<k> suffix, with a cross-table in the JSON. A larger k is
only preferred when it adds a lesion class (not just another margin/normal one).

Outputs (default BIDS/derivatives/KUL_compute/sub-<p>/KUL_lesion_compartments/):
  <prefix>_compartments.nii.gz        labels after light spatial regularisation
  <prefix>_compartments_raw.nii.gz    labels straight from the mixture model
  <prefix>_compartments_posteriors.nii.gz
  <prefix>_compartments_zone.nii.gz   1 = original mask, 2 = dilation ring
  <prefix>_compartments_itksnap.txt   ITK-SNAP label description (rename there)
  <prefix>_compartments.json          everything used and measured
  <prefix>_compartments_QA.png
  <prefix>_compartments_validation.csv  (only when reference masks exist)
  <prefix>_compartments_k<k>*          the alternative k, same set of files
""")
    ap.add_argument("-p", "--participant", help="participant (without sub-), study-dir mode")
    ap.add_argument("--preset", default="generic", choices=sorted(PRESETS),
                    help="lesion type used to NAME the classes (default generic)")
    ap.add_argument("--mask", help="lesion mask (any value > 0 is lesion)")
    ap.add_argument("--channel", action="append", default=[], metavar="NAME=PATH",
                    help="add/override a channel (T1, FLAIR, T2, cT1sub, SWI, ADC or any name)")
    ap.add_argument("--exclude", action="append", default=[], metavar="NAME", help="drop a channel")
    ap.add_argument("--aseg", help="FreeSurfer aseg (NAWM reference + CSF exclusion)")
    ap.add_argument("--ref", action="append", default=[], metavar="NAME=PATH",
                    help="reference mask for validation (repeatable)")
    ap.add_argument("-o", "--outdir")
    ap.add_argument("--prefix")
    ap.add_argument("--dilate", type=float, default=4.0, help="ring around the mask, mm (default 4)")
    ap.add_argument("--nawm-margin", type=float, default=10.0,
                    help="NAWM excludes this many mm around the lesion (default 10)")
    ap.add_argument("--ks", type=lambda s: [int(v) for v in s.split(",")], default=[3, 4],
                    help="candidate numbers of classes, all fitted and written (default 3,4); the "
                         "largest k whose classes all pass the checks becomes the main output")
    ap.add_argument("--k", type=int, help="force the number of classes")
    ap.add_argument("--cluster", action="store_true",
                    help="cluster even when HD-GLIO compartments exist (types 1/2)")
    ap.add_argument("--smooth", type=float, default=0.4,
                    help="posterior smoothing sigma in voxels, 0 = off (default 0.4: only "
                         "ambiguous voxels can change class)")
    ap.add_argument("--max-fit", type=int, default=100000, help="max voxels used to fit (default 1e5)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(1)

    p = args.participant
    if not p and not (args.mask and args.channel and args.outdir):
        sys.exit("ERROR: give -p <participant> or --mask, --channel(s) and -o")

    # ---------- inputs
    chan_paths, native = {}, {}
    if p:
        for name, tmpl in CHANNELS:
            f = tmpl.format(p=p)
            if os.path.isfile(f):
                chan_paths[name] = f
        for name, tmpl in NATIVE.items():
            hits = sorted(glob.glob(f"BIDS/sub-{p}/anat/{tmpl.format(p=p)}") +
                          glob.glob(f"BIDS/sub-{p}/ses-*/anat/{tmpl.format(p=p)}") +
                          glob.glob(f"BIDS/derivatives/KUL_compute/sub-{p}/cT1w/{tmpl.format(p=p)}"))
            if hits:
                native[name] = [round(float(z), 2) for z in nib.load(hits[0]).header.get_zooms()[:3]]
    for c in args.channel:
        name, f = c.split("=", 1)
        chan_paths[name] = f
    for name in args.exclude:
        chan_paths.pop(name, None)
    if "ADC" in chan_paths:
        native["ADC"] = [round(float(z), 2) for z in nib.load(chan_paths["ADC"]).header.get_zooms()[:3]]

    mask_path = args.mask or next((m.format(p=p) for m in MASKS if os.path.isfile(m.format(p=p))), None)
    if not mask_path or not os.path.isfile(mask_path):
        sys.exit("ERROR: no lesion mask found (tried: " + ", ".join(m.format(p=p) for m in MASKS) + ")")
    aseg_path = args.aseg or (ASEG.format(p=p) if p else None)
    if not aseg_path or not os.path.isfile(aseg_path):
        sys.exit(f"ERROR: aseg needed for the NAWM reference, not found: {aseg_path}")
    refs = {}
    if p:
        for name, tmpl in REFS:
            if os.path.isfile(tmpl.format(p=p)):
                refs[name] = tmpl.format(p=p)
    for r in args.ref:
        name, f = r.split("=", 1)
        refs[name] = f

    outdir = args.outdir or f"BIDS/derivatives/KUL_compute/sub-{p}/KUL_lesion_compartments"
    prefix = args.prefix or (f"sub-{p}" if p else "lesion")
    os.makedirs(outdir, exist_ok=True)
    out = lambda s: os.path.join(outdir, f"{prefix}_compartments{s}")

    ref_img = nib.load(mask_path)
    zooms = np.array(ref_img.header.get_zooms()[:3], dtype=float)
    vox_ml = float(np.prod(zooms)) / 1000
    mask = np.asanyarray(ref_img.dataobj) > 0
    log(f"mask    {mask_path}: {mask.sum() * vox_ml:.2f} mL")

    # ---------- channels on the mask grid
    data, warn = {}, []
    for name, f in chan_paths.items():
        a, resampled = load_on_grid(f, ref_img, order=1)
        data[name] = a
        log(f"channel {name:7s}{f}" + ("  (resampled to mask grid)" if resampled else ""))
        if resampled and name != "ADC":
            warn.append(f"{name} was not on the mask grid and was resampled (linear)")
    if not data:
        sys.exit("ERROR: no channels found")

    aseg, _ = load_on_grid(aseg_path, ref_img, order=0)
    aseg = np.rint(aseg).astype(int)
    csf = np.isin(aseg, ASEG_CSF)
    brain = (aseg > 0) & ~csf

    # ---------- zone = mask + ring (ring restricted to non-CSF brain)
    dist = ndi.distance_transform_edt(~mask, sampling=zooms)
    ring = (dist > 0) & (dist <= args.dilate) & brain
    zone = mask | ring

    valid = {n: np.isfinite(a) & ((a != 0) if n not in ZERO_OK else True) for n, a in data.items()}
    dropped = {}
    for n in list(data):
        miss = 1 - valid[n][zone].mean()
        if miss > 0.10:
            dropped[n] = f"missing in {miss:.0%} of the lesion zone"
            log(f"WARNING: dropping {n}: {dropped[n]}")
            del data[n], valid[n]
    chn = list(data)
    zv = zone.copy()
    for n in chn:
        zv &= valid[n]
    n_missing = int((zone & ~zv).sum())

    # ---------- NAWM reference and z-scores
    wm = ndi.binary_erosion(np.isin(aseg, ASEG_WM), iterations=1)
    nawm = wm & (dist > args.nawm_margin)
    for n in chn:
        nawm &= valid[n]
    if nawm.sum() < 5000:
        sys.exit(f"ERROR: only {nawm.sum()} NAWM voxels; check aseg alignment")
    cortex = np.isin(aseg, ASEG_CORTEX) & valid.get("T1", True)
    aseg_check = None
    if "T1" in chn and cortex.sum() > 0:
        aseg_check = float(np.median(data["T1"][nawm]) / np.median(data["T1"][cortex]))
        if aseg_check < 1.1:
            warn.append(f"aseg may be misaligned: T1 WM/cortex median ratio {aseg_check:.2f} (expect >1.2)")

    # Z = z vs NAWM (reported); X = signed-log compressed z (fitted). NAWM is so
    # homogeneous that lesion z runs to 10-50: compression keeps those
    # differences instead of clipping them away, and stops one channel's tail
    # from dominating the mixture.
    Z = np.empty((int(zv.sum()), len(chn)), dtype=np.float32)
    nawm_stats = {}
    for j, n in enumerate(chn):
        z, med, mad = robust_z(data[n][zv], data[n][nawm])
        Z[:, j] = np.clip(z, -100, 100)
        nawm_stats[n] = {"median": med, "mad_sd": mad}
    X = np.sign(Z) * np.log1p(np.abs(Z))
    log(f"zone    {zv.sum() * vox_ml:.2f} mL usable ({mask.sum() * vox_ml:.2f} mask + "
        f"{ring.sum() * vox_ml:.2f} ring; {n_missing} voxels lack a channel); NAWM {nawm.sum()} voxels")

    # ---------- components of the original mask (26-connected), 1 = largest
    comp_lab, ncomp = ndi.label(mask, structure=np.ones((3, 3, 3)))
    sizes = ndi.sum(mask, comp_lab, range(1, ncomp + 1))
    rank = np.zeros(ncomp + 1, int)
    rank[np.argsort(-sizes) + 1] = np.arange(1, ncomp + 1)
    comp_lab = rank[comp_lab]
    sizes = np.sort(sizes)[::-1]
    comp_order = list(range(1, ncomp + 1))
    # ring voxels belong to the nearest component
    _, near = ndi.distance_transform_edt(~mask, sampling=zooms, return_indices=True)
    comp_zone = np.where(mask, comp_lab, comp_lab[tuple(near)]) * zone
    in_mask = mask[zv]
    zone_ml = float(zv.sum() * vox_ml)

    def profile(sel):
        return ({n: float(np.median(Z[sel, j])) for j, n in enumerate(chn)},
                {n: [float(q) for q in np.percentile(Z[sel, j], [25, 75])] for j, n in enumerate(chn)})

    def describe(lab, names, extra=None):
        """Per-class volumes, per-component volumes and z profiles of a label image."""
        classes = []
        lab_zv = lab[zv]
        for i, name in enumerate(names):
            lv = i + 1
            cm = lab == lv
            sel = lab_zv == lv
            med, iqr = profile(sel) if sel.any() else ({n: 0.0 for n in chn}, {n: [0.0, 0.0] for n in chn})
            per_comp = [{"component": int(cid),
                         "vol_mask_ml": round(float((cm & mask & (comp_lab == cid)).sum() * vox_ml), 4),
                         "vol_ring_ml": round(float((cm & ring & (comp_zone == cid)).sum() * vox_ml), 4)}
                        for cid in comp_order]
            cl = {"label": lv, "name": name, "descriptor": descriptor(med),
                  "colour": COLOURS[i % len(COLOURS)],
                  "median_z": {n: round(v, 2) for n, v in med.items()},
                  "iqr_z": {n: [round(q, 2) for q in v] for n, v in iqr.items()},
                  "vol_mask_ml": round(float((cm & mask).sum() * vox_ml), 4),
                  "vol_ring_ml": round(float((cm & ring).sum() * vox_ml), 4),
                  "per_component": per_comp}
            if extra:
                cl.update(extra[i])
            classes.append(cl)
        return classes

    rng = np.random.default_rng(args.seed)
    fit_idx = np.arange(len(X)) if len(X) <= args.max_fit else rng.choice(len(X), args.max_fit, replace=False)
    fit = X[fit_idx]

    def gmm(k, data_):
        return GaussianMixture(k, covariance_type="full", n_init=3, reg_covar=1e-3,
                               random_state=args.seed, max_iter=500).fit(data_)

    def solve(k):
        """Fit k classes, regularise lightly, name, and measure how real each class is."""
        gm = gmm(k, fit)
        post = gm.predict_proba(X)
        raw = post.argmax(1)
        if args.smooth > 0:
            w = ndi.gaussian_filter(zv.astype(np.float32), args.smooth)
            post_s = np.empty_like(post)
            for c in range(k):
                vol = np.zeros(mask.shape, np.float32)
                vol[zv] = post[:, c]
                post_s[:, c] = ndi.gaussian_filter(vol, args.smooth)[zv] / w[zv]
            reg = post_s.argmax(1)
        else:
            reg = raw
        flipped = raw != reg
        flip_info = {"sigma_vox": args.smooth, "voxels_changed": int(flipped.sum()),
                     "fraction_changed": float(flipped.mean()),
                     "max_own_posterior_of_changed_voxel":
                         float(post[flipped].max(1).max()) if flipped.any() else None}

        meds, frac_ring = [], []
        for c in range(k):
            sel = reg == c
            if not sel.any():
                sel = raw == c
            meds.append(profile(sel)[0])
            frac_ring.append(float(1 - in_mask[sel].mean()) if sel.any() else 1.0)
        names = name_classes(meds, frac_ring, args.preset)
        prio = PRESETS[args.preset]
        order = sorted(range(k), key=lambda c: (prio.index(names[c]) if names[c] in prio else len(prio),
                                                frac_ring[c]))
        remap = np.empty(k, int)
        remap[order] = np.arange(1, k + 1)
        lab = np.zeros(mask.shape, np.uint8)
        lab[zv] = remap[reg]
        lab_raw = np.zeros(mask.shape, np.uint8)
        lab_raw[zv] = remap[raw]
        post_vol = np.zeros(mask.shape + (k,), np.float32)
        post_vol[zv] = post[:, order]

        # quality: split-half reproducibility, overlap with other classes,
        # spatial coherence, size
        rep = split_half_dice(X, fit, k, raw, gmm, args.seed, nrep=2)
        bc = bhattacharyya(gm)
        extra = []
        for i, c in enumerate(order):
            lv = i + 1
            coh = coherence(lab, zv, lv)
            vol = float((lab == lv).sum() * vox_ml)
            q = {"reproducibility_dice": round(float(rep[c].mean()), 3),
                 "reproducibility_dice_min": round(float(rep[c].min()), 3),
                 "max_overlap_bhattacharyya": round(float(np.delete(bc[c], c).max()), 3) if k > 1 else 0.0,
                 "coherence": round(coh, 3),
                 "vol_ml": round(vol, 4), "fraction_of_zone": round(vol / zone_ml, 4),
                 "fraction_in_ring": round(frac_ring[c], 3),
                 "mixture_weight": round(float(gm.weights_[c]), 4)}
            fails = []
            if q["reproducibility_dice"] < QC["repro"]:
                fails.append("not reproducible")
            if q["coherence"] < QC["coherence"]:
                fails.append("scattered")
            if vol < max(QC["min_ml"], QC["min_frac"] * zone_ml):
                fails.append("too small")
            q["checks_failed"] = fails
            extra.append(q)
        classes = describe(lab, unique_names([names[c] for c in order]), extra)
        log(f"k={k}    " + "; ".join(
            f"{cl['name']} repro {cl['reproducibility_dice']:.2f} coh {cl['coherence']:.2f} "
            f"{cl['vol_ml']:.2f} mL" + (" FAIL: " + ",".join(cl['checks_failed']) if cl['checks_failed'] else "")
            for cl in classes))
        return {"k": k, "lab": lab, "lab_raw": lab_raw, "post_vol": post_vol, "classes": classes,
                "flip_info": flip_info, "passes": all(not cl["checks_failed"] for cl in classes)}

    def write(sol, tag, title, bic, extra_json):
        o = lambda s: out(tag + s)
        hdr = ref_img.header.copy()
        hdr.set_data_dtype(np.uint8)
        nib.save(nib.Nifti1Image(sol["lab"], ref_img.affine, hdr), o(".nii.gz"))
        if "lab_raw" in sol:
            nib.save(nib.Nifti1Image(sol["lab_raw"], ref_img.affine, hdr), o("_raw.nii.gz"))
            nib.save(nib.Nifti1Image(sol["post_vol"], ref_img.affine), o("_posteriors.nii.gz"))
        with open(o("_itksnap.txt"), "w") as fh:
            fh.write('    0     0    0    0        0  0  0    "Clear Label"\n')
            for cl in sol["classes"]:
                r, g, b = (int(cl["colour"][i:i + 2], 16) for i in (1, 3, 5))
                fh.write(f'{cl["label"]:5d} {r:5d} {g:4d} {b:4d}        1  1  1    '
                         f'"{cl["name"]} [{cl["descriptor"]}]"\n')
        validation = []
        if refs:
            rows = ["class_label,class_name,ref,ref_vol_ml,class_vol_ml,overlap_ml,dice,frac_of_ref,frac_of_class"]
            for rname, rm in ref_masks.items():
                for cl in sol["classes"]:
                    cm = sol["lab"] == cl["label"]
                    ov = (cm & rm).sum()
                    d = 2 * ov / (cm.sum() + rm.sum()) if (cm.sum() + rm.sum()) else 0
                    row = {"class_label": cl["label"], "class_name": cl["name"], "ref": rname,
                           "ref_vol_ml": round(float(rm.sum() * vox_ml), 4),
                           "class_vol_ml": round(float(cm.sum() * vox_ml), 4),
                           "overlap_ml": round(float(ov * vox_ml), 4), "dice": round(float(d), 3),
                           "frac_of_ref": round(float(ov / rm.sum()), 3) if rm.sum() else None,
                           "frac_of_class": round(float(ov / cm.sum()), 3) if cm.sum() else None}
                    validation.append(row)
                    rows.append(",".join(str(row[h]) for h in rows[0].split(",")))
            with open(o("_validation.csv"), "w") as fh:
                fh.write("\n".join(rows) + "\n")
        summary = dict(common, **extra_json, classes=sol["classes"], validation=validation)
        if "flip_info" in sol:
            summary["spatial_regularisation"] = sol["flip_info"]
        with open(o(".json"), "w") as fh:
            json.dump(summary, fh, indent=2)
        qa_figure(o("_QA.png"), title, [(n, ras[n]) for n in chn], ras, zooms_ras, mask_ras,
                  to_ras(sol["lab"], ref_img.affine), comps_ras, sol["classes"], bic, sol["k"],
                  [cl["colour"] for cl in sol["classes"]])
        return o("")

    # shared by every output
    ref_masks = {}
    for rname, rpath in refs.items():
        rm, _ = load_on_grid(rpath, ref_img, order=0)
        ref_masks[rname] = (rm > 0.5) & zone
    zone_img = (mask * 1 + ring * 2).astype(np.uint8)
    hdr = ref_img.header.copy()
    hdr.set_data_dtype(np.uint8)
    nib.save(nib.Nifti1Image(zone_img, ref_img.affine, hdr), out("_zone.nii.gz"))
    common = {
        "tool": "KUL_lesion_compartments.py", "version": VERSION,
        "participant": p, "preset": args.preset, "mask": mask_path, "aseg": aseg_path,
        "channels": {n: {"path": chan_paths[n], "native_voxel_mm": native.get(n),
                         "nawm_median": float(f'{nawm_stats[n]["median"]:.6g}'),
                         "nawm_robust_sd": float(f'{nawm_stats[n]["mad_sd"]:.6g}')} for n in chn},
        "dropped_channels": dropped, "warnings": warn,
        "aseg_check_T1_wm_over_cortex": round(aseg_check, 3) if aseg_check else None,
        "zone": {"dilate_mm": args.dilate, "mask_ml": round(float(mask.sum() * vox_ml), 4),
                 "ring_ml": round(float(ring.sum() * vox_ml), 4),
                 "voxels_missing_a_channel": n_missing, "nawm_voxels": int(nawm.sum()),
                 "nawm_margin_mm": args.nawm_margin},
        "components": [{"component": int(cid), "vol_ml": round(float(sizes[cid - 1] * vox_ml), 4),
                        "centroid_mm": [round(float(v), 1) for v in nib.affines.apply_affine(
                            ref_img.affine, ndi.center_of_mass(comp_lab == cid))]} for cid in comp_order],
    }
    ras = {n: to_ras(data[n], ref_img.affine) for n in chn}
    ornt = nib.orientations.io_orientation(ref_img.affine)
    zooms_ras = [zooms[int(ax)] for ax in ornt[:, 0]]
    comps_ras = {"lab": to_ras(comp_lab, ref_img.affine), "ids": comp_order}
    mask_ras = to_ras(mask, ref_img.affine)

    # ---------- types 1/2: HD-GLIO (+ resseg cavity) already gives the compartments
    hd = {name: tmpl.format(p=p) for name, tmpl in HDGLIO_CLASSES} if p else {}
    hd = {n: f for n, f in hd.items() if os.path.isfile(f)}
    if p and not args.cluster and "nonenhancing_FLAIR_hdglio" in hd:
        lab = np.zeros(mask.shape, np.uint8)
        names = []
        for name, f in hd.items():          # earlier entries win where they overlap
            m, _ = load_on_grid(f, ref_img, order=0)
            m = (m > 0.5) & (lab == 0)
            names.append(name)
            lab[m] = len(names)
        refs.clear()                        # the references ARE these labels
        sol = {"k": len(names), "lab": lab, "classes": describe(lab, names)}
        sol["classes"] = [dict(cl, fraction_in_ring=0.0) for cl in sol["classes"]]
        path = write(sol, "", f"{prefix}  lesion compartments from HD-GLIO/resseg (not clustered)",
                     None, {"status": "HD-GLIO/resseg compartments adopted as-is (use --cluster to cluster)",
                            "source": hd})
        log("")
        for cl in sol["classes"]:
            log(f"  {cl['label']}  {cl['name']:26s} {cl['vol_mask_ml']:8.3f} mL   [{cl['descriptor']}]")
        log(f"\nHD-GLIO compartments adopted; written to {outdir}/")
        return

    # ---------- cluster: fit every candidate k, keep the largest k whose
    # classes are all reproducible, coherent and not tiny
    ks = [args.k] if args.k else sorted(set(args.ks))
    bic = {kk: float(gmm(kk, fit).bic(fit)) for kk in range(1, max(ks) + 2)}
    log("BIC     " + "  ".join(f"k={kk}: {b:.4g}" for kk, b in bic.items()))
    sols = {kk: solve(kk) for kk in ks}
    # among the k whose classes all pass, step up only when the larger k adds a
    # lesion class (margin / normal subdivisions alone are not worth a class)
    passing = [kk for kk in ks if sols[kk]["passes"]]
    n_lesion = {kk: sum(not cl["name"].startswith(NON_LESION) for cl in sols[kk]["classes"]) for kk in ks}
    chosen = min(passing) if passing else min(ks)
    for kk in passing:
        if kk > chosen and n_lesion[kk] > n_lesion[chosen]:
            chosen = kk
    if not passing:
        warn.append(f"no k in {ks} passed all class checks; kept k={chosen}, review the QA")
    if args.k:
        chosen = args.k

    # how the classes of each k map onto those of the chosen k
    crosstab = {}
    for kk in ks:
        if kk == chosen:
            continue
        a, b = sols[chosen], sols[kk]
        crosstab[f"k{chosen}_vs_k{kk}"] = [
            {"k%d_class" % chosen: ca["name"],
             **{cb["name"]: round(float(((a["lab"] == ca["label"]) & (b["lab"] == cb["label"])).sum() * vox_ml), 3)
                for cb in b["classes"]}} for ca in a["classes"]]

    k_info = {"type": "GaussianMixture full covariance on signed-log z", "bic": bic, "candidates": ks,
              "chosen_k": chosen, "fit_voxels": int(len(fit)), "seed": args.seed,
              "k_rule": "forced (--k)" if args.k else
              (f"largest k in {ks} whose classes all pass: split-half Dice >= {QC['repro']}, "
               f"coherence >= {QC['coherence']}, volume >= max({QC['min_ml']} mL, {QC['min_frac']:.0%} of zone); "
               f"a larger k only when it adds a lesion class (not margin/normal)"),
              "passing_k": passing, "lesion_classes_per_k": n_lesion, "crosstab_ml": crosstab}
    for kk in ks:
        tag = "" if kk == chosen else f"_k{kk}"
        status = ("suggestion - review the QA png" if kk == chosen else
                  f"alternative (k={kk}); the chosen solution is k={chosen}")
        write(sols[kk], tag,
              f"{prefix}  lesion compartments ({args.preset}, k={kk}"
              f"{', chosen' if kk == chosen else ', alternative'})  -  suggestion, review before use",
              bic, {"status": status, "model": dict(k_info, k=kk)})

    log("")
    for kk in ks:
        log(f"k={kk}{'  <- chosen' if kk == chosen else ''}")
        for cl in sols[kk]["classes"]:
            log(f"  {cl['label']}  {cl['name']:26s} {cl['vol_mask_ml']:8.3f} mL in mask  "
                f"{cl['vol_ring_ml']:8.3f} mL ring   [{cl['descriptor']}]"
                + (f"   FAIL: {', '.join(cl['checks_failed'])}" if cl["checks_failed"] else ""))
    for key, rows in crosstab.items():
        log(f"\n{key} (mL; rows = chosen classes)")
        for r in rows:
            first = list(r)[0]
            log(f"  {r[first]:26s} " + "  ".join(f"{n}: {v:.2f}" for n, v in list(r.items())[1:] if v > 0))
    for w_ in warn:
        log("WARNING: " + w_)
    log(f"\nwritten to {outdir}/")


if __name__ == "__main__":
    main()
