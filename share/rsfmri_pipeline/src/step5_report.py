#!/usr/bin/env python3
# Run with:  python3 step5_report.py [options]
"""
Step 5 — PDF Report Generator

Reads step4_compare.py outputs and generates a per-patient PDF with:
  • Brain maps in orthogonal views (sagittal / coronal / axial):
      - HV group mean FC map
      - Patient FC map (average over runs)
      - Patient z-score map
      - Patient LOO-significant map  (empirical threshold)
  • ROI z-score summary table

Background image: patient's own T1w in MNI space (fmriprep output),
falling back to MNI152 template if not found.

Usage:
  python step5_report.py --profile MDD --patients PT02
  python step5_report.py --rsn --patients PT02
  python step5_report.py --profile MDD --rsn --patients PT01 PT02
"""

import argparse
import io
import logging
import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.image import resample_to_img

sys.path.insert(0, str(Path(__file__).parent))
from config import ANALYSIS_DIR, FSLDIR, SBA_DIR, RSN_FC_DIR, PT_SUBJECTS, get_profile, seed_subdir
from utils import t1w_path, YEO17_NETWORK_NAMES

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

NORM_DIR     = SBA_DIR    / "normative"
RSN_NORM_DIR = RSN_FC_DIR / "normative"
COMP_DIR     = ANALYSIS_DIR / "comparison"
REPORT_DIR   = ANALYSIS_DIR / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

FSL_TEMPLATE = FSLDIR / "data" / "standard" / "MNI152_T1_2mm_brain.nii.gz"

_ALL_RSN_NAMES = YEO17_NETWORK_NAMES + ["Language", "Auditory"]

_CMAP = "RdBu_r"

_MNI_BG: nib.Nifti1Image | None = None


def _mni_bg() -> nib.Nifti1Image:
    global _MNI_BG
    if _MNI_BG is None:
        if FSL_TEMPLATE.exists():
            _MNI_BG = nib.load(str(FSL_TEMPLATE))
        else:
            from nilearn import datasets
            try:
                _MNI_BG = datasets.load_mni152_template(resolution=2)
            except TypeError:
                _MNI_BG = datasets.load_mni152_template()
    return _MNI_BG


def _patient_bg(subject: str) -> nib.Nifti1Image:
    """Load patient T1w skull-stripped to MNI space, or fall back to MNI template."""
    t1w_p = t1w_path(subject)
    if not t1w_p.exists():
        log.warning("  T1w not found for sub-%s — using MNI template", subject)
        return _mni_bg()
    # Apply the fmriprep brain mask so bg_data > 0 equals the true brain mask.
    # desc-preproc_T1w includes skull; without masking, skull voxels pass the
    # brain_mask check in _plot_ortho and FC interpolation bleed shows there.
    mask_p = t1w_p.parent / t1w_p.name.replace("desc-preproc_T1w", "desc-brain_mask")
    t1w_img = nib.load(str(t1w_p))
    if mask_p.exists():
        log.info("  Background: patient T1w + brain mask (%s)", t1w_p.name)
        mask_data = nib.load(str(mask_p)).get_fdata(dtype=np.float32)
        t1w_data  = t1w_img.get_fdata(dtype=np.float32) * mask_data
        return nib.Nifti1Image(t1w_data, t1w_img.affine, t1w_img.header)
    log.warning("  Brain mask not found for sub-%s — T1w may include skull", subject)
    return t1w_img


# ─────────────────────────────────────────────────────────────────────────────
# Path helpers  (mirrors step4)
# ─────────────────────────────────────────────────────────────────────────────

def _seed_paths(patient: str, seed: str, thresholded: bool = False) -> dict[str, Path]:
    d    = SBA_DIR / f"sub-{patient}" / seed_subdir(seed)   # "" for non-Lausanne seeds
    t    = "Thresh" if thresholded else ""
    mean_tag = "_thresh" if thresholded else ""
    return {
        "mean":   NORM_DIR / f"seed-{seed}_HV_GLM_mean{mean_tag}.nii.gz",
        "pt_fc":  d / f"sub-{patient}_seed-{seed}_desc-sba{t}_statmap.nii.gz",
        "zscore": d / f"sub-{patient}_seed-{seed}_desc-zscore{t}_statmap.nii.gz",
        "pval":   d / f"sub-{patient}_seed-{seed}_desc-pval{t}_statmap.nii.gz",
        "loosig": d / f"sub-{patient}_seed-{seed}_desc-loosig{t}_statmap.nii.gz",
    }


def _rsn_paths(patient: str, rsn: str, thresholded: bool = False) -> dict[str, Path]:
    d    = RSN_FC_DIR / f"sub-{patient}"
    t    = "Thresh" if thresholded else ""
    mean_tag = "_thresh" if thresholded else ""
    return {
        "mean":   RSN_NORM_DIR / f"rsn-{rsn}_HV_GLM_mean{mean_tag}.nii.gz",
        "pt_fc":  d / f"sub-{patient}_rsn-{rsn}_desc-sba{t}_statmap.nii.gz",
        "zscore": d / f"sub-{patient}_rsn-{rsn}_desc-zscore{t}_statmap.nii.gz",
        "pval":   d / f"sub-{patient}_rsn-{rsn}_desc-pval{t}_statmap.nii.gz",
        "loosig": d / f"sub-{patient}_rsn-{rsn}_desc-loosig{t}_statmap.nii.gz",
    }


# ─────────────────────────────────────────────────────────────────────────────
# Brain slice rendering
# ─────────────────────────────────────────────────────────────────────────────

def _peak_mm(
    img:    nib.Nifti1Image,
    bg_img: nib.Nifti1Image,
) -> tuple[float, float, float]:
    """MNI mm coordinates of peak |value| voxel, restricted to brain voxels of bg_img."""
    try:
        img_r  = resample_to_img(img, bg_img, interpolation="linear", force_resample=True)
        data   = img_r.get_fdata(dtype=np.float32)
        mask   = bg_img.get_fdata(dtype=np.float32) > 0
        masked = np.abs(data) * mask
        if masked.max() < 1e-6:
            return 0.0, 0.0, 0.0
        idx = np.unravel_index(np.argmax(masked), data.shape)
        mm  = (bg_img.affine @ np.array([*idx, 1.0]))[:3]
        return float(mm[0]), float(mm[1]), float(mm[2])
    except Exception:
        return 0.0, 0.0, 0.0


def _mm_to_vox(affine: np.ndarray, x: float, y: float, z: float) -> tuple[int, int, int]:
    inv = np.linalg.inv(affine)
    v   = inv @ np.array([x, y, z, 1.0])
    return int(round(v[0])), int(round(v[1])), int(round(v[2]))


def _prepare_overlay_rgba(
    stat_img:  nib.Nifti1Image,
    threshold: float | None,
    vmax:      float | None,
    bg_img:    nib.Nifti1Image,
    cmap:      str = _CMAP,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Resample stat_img onto bg_img's grid and colour-map it into an RGBA
    volume. Shared by _plot_ortho() (single peak-slice view) and
    _plot_ortho_and_axial_montage() (adds a multi-slice lightbox) so both
    stay visually consistent -- same thresholding, same vmax, same colours.
    Returns (bg_norm, stat_rgba, vmax).
    """
    bg_data = bg_img.get_fdata(dtype=np.float32)

    stat_r    = resample_to_img(stat_img, bg_img, interpolation="linear", force_resample=True)
    stat_data = stat_r.get_fdata(dtype=np.float32)

    # Mask stat overlay to brain voxels
    brain_mask = bg_data > 0
    stat_data  = stat_data * brain_mask

    bg_max  = bg_data.max()
    bg_norm = bg_data / bg_max if bg_max > 0 else bg_data

    # Zero sub-threshold voxels before colour scaling so they become fully
    # transparent and do not bias vmax toward weak (near-white) values.
    if threshold is not None:
        stat_data[np.abs(stat_data) < threshold] = 0.0

    if vmax is None or vmax <= 0:
        nonzero = np.abs(stat_data[stat_data != 0])
        vmax    = float(np.percentile(nonzero, 95)) if len(nonzero) > 0 else 1.0
    vmax = max(vmax, 0.01)

    cmap_obj = plt.get_cmap(cmap)
    thr = threshold if (threshold is not None and threshold > 0) else 0.0
    if cmap == _CMAP:                                      # diverging (RdBu_r)
        if thr > 0 and vmax > thr:
            # Map supra-threshold values to the outer coloured bands of RdBu_r,
            # staying well away from the white centre (0.5).
            # In RdBu_r: [0.7, 1.0] = clearly blue, [0.3, 0.0] = clearly red.
            # Below-threshold voxels are already zeroed → alpha = 0.
            scale = vmax - thr
            stat_norm = np.where(
                stat_data > 0,
                0.7 + 0.3 * np.clip((stat_data - thr) / scale, 0.0, 1.0),
                np.where(
                    stat_data < 0,
                    0.3 - 0.3 * np.clip((-stat_data - thr) / scale, 0.0, 1.0),
                    0.5,
                ),
            )
        else:
            stat_norm = np.clip(stat_data / vmax, -1.0, 1.0) * 0.5 + 0.5
    else:                                                  # sequential (positive-only)
        if thr > 0 and vmax > thr:
            # Map [threshold, vmax] to [0.35, 1.0] so even the weakest visible
            # voxel gets a clearly non-white colour in the Reds colormap.
            scale = vmax - thr
            stat_norm = np.where(
                stat_data > 0,
                0.35 + 0.65 * np.clip((stat_data - thr) / scale, 0.0, 1.0),
                0.0,
            )
        else:
            stat_norm = np.clip(stat_data / vmax, 0.0, 1.0)
    stat_rgba = cmap_obj(stat_norm).astype(np.float32)

    alpha = (stat_data != 0).astype(np.float32) * 0.88
    stat_rgba[..., 3] = alpha

    return bg_norm, stat_rgba, vmax


def _add_lr_labels(ax: plt.Axes) -> None:
    """
    Radiological convention label (screen-left = patient's R, screen-right =
    patient's L) for coronal/axial panels -- pairs with the R-L axis flip
    applied to those slices below. Sagittal slices show no L/R axis, so
    they don't get this label.
    """
    ax.text(0.02, 0.02, "R", color="yellow", fontsize=9, fontweight="bold",
            transform=ax.transAxes, ha="left", va="bottom")
    ax.text(0.98, 0.02, "L", color="yellow", fontsize=9, fontweight="bold",
            transform=ax.transAxes, ha="right", va="bottom")


def _plot_ortho(
    stat_img:  nib.Nifti1Image,
    threshold: float | None,
    vmax:      float | None,
    title:     str,
    cut_mm:    tuple[float, float, float],
    bg_img:    nib.Nifti1Image,
    cmap:      str = _CMAP,
) -> plt.Figure:
    """
    Render sagittal / coronal / axial slices of stat_img onto bg_img.
    bg_img is expected to be brain-extracted (nonzero = brain).
    Coronal/axial are shown in radiological convention (screen-left =
    patient's right), labelled explicitly -- see _add_lr_labels().
    Returns a matplotlib Figure for embedding in PDF.
    """
    bg_data = bg_img.get_fdata(dtype=np.float32)
    shape   = bg_data.shape[:3]

    bg_norm, stat_rgba, vmax = _prepare_overlay_rgba(stat_img, threshold, vmax, bg_img, cmap)

    xi, yi, zi = _mm_to_vox(bg_img.affine, *cut_mm)
    xi = int(np.clip(xi, 0, shape[0] - 1))
    yi = int(np.clip(yi, 0, shape[1] - 1))
    zi = int(np.clip(zi, 0, shape[2] - 1))

    # bg_img/stat data are RAS-oriented (axis 0 increases toward patient's
    # right); reversing axis 0 before slicing puts patient's right on
    # screen-left for coronal/axial, i.e. radiological convention.
    slices = [
        (bg_norm[xi, :, :].T,
         stat_rgba[xi, :, :, :].transpose(1, 0, 2),
         "Sagittal", f"x = {cut_mm[0]:+.0f} mm", False),
        (bg_norm[::-1, yi, :].T,
         stat_rgba[::-1, yi, :, :].transpose(1, 0, 2),
         "Coronal",  f"y = {cut_mm[1]:+.0f} mm", True),
        (bg_norm[::-1, :, zi].T,
         stat_rgba[::-1, :, zi, :].transpose(1, 0, 2),
         "Axial",    f"z = {cut_mm[2]:+.0f} mm", True),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5),
                             gridspec_kw={"wspace": 0.04})
    fig.patch.set_facecolor("#111111")

    for ax, (bg_sl, stat_sl, plane, coord_lbl, needs_lr) in zip(axes, slices):
        ax.set_facecolor("black")
        ax.imshow(bg_sl,   origin="lower", cmap="gray", aspect="auto",
                  interpolation="bilinear", vmin=0, vmax=1)
        ax.imshow(stat_sl, origin="lower", aspect="auto",
                  interpolation="nearest")
        ax.set_title(f"{plane}\n{coord_lbl}", color="white", fontsize=8, pad=3)
        ax.axis("off")
        if needs_lr:
            _add_lr_labels(ax)

    cbar_norm = mcolors.Normalize(0, vmax) if cmap != _CMAP else mcolors.Normalize(-vmax, vmax)
    sm   = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=cbar_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(),
                        orientation="vertical", fraction=0.018, pad=0.02, shrink=0.85)
    cbar.ax.tick_params(colors="white", labelsize=7)
    cbar.set_label("", color="white")

    fig.suptitle(title, color="white", fontsize=9, y=1.0, va="bottom")
    return fig


# Standard MNI z-coordinates for the axial lightbox montage, spanning from
# below the cerebellum to the vertex -- fixed regardless of where any given
# map's peak voxel is, so the montage always shows the same anatomical
# levels and is comparable across seeds/RSNs and across subjects.
_LIGHTBOX_Z_MM = (-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60)


def _plot_ortho_and_axial_montage(
    stat_img:     nib.Nifti1Image,
    threshold:    float | None,
    vmax:         float | None,
    title:        str,
    cut_mm:       tuple[float, float, float],
    bg_img:       nib.Nifti1Image,
    cmap:         str = _CMAP,
    z_coords_mm:  tuple[float, ...] = _LIGHTBOX_Z_MM,
    montage_cols: int = 6,
) -> plt.Figure:
    """
    Combined page: the orthogonal view (sagittal/coronal/axial through the
    peak voxel, as in _plot_ortho()) on top, plus an axial lightbox montage
    below -- a fixed grid of evenly-spaced axial slices across standard MNI
    z-levels. The single peak slice shown by the ortho view alone can miss
    how far a cluster actually extends; the montage shows the full spatial
    footprint at a glance, at the cost of using more page space per page.
    """
    bg_data = bg_img.get_fdata(dtype=np.float32)
    shape   = bg_data.shape[:3]

    bg_norm, stat_rgba, vmax = _prepare_overlay_rgba(stat_img, threshold, vmax, bg_img, cmap)

    xi, yi, zi = _mm_to_vox(bg_img.affine, *cut_mm)
    xi = int(np.clip(xi, 0, shape[0] - 1))
    yi = int(np.clip(yi, 0, shape[1] - 1))
    zi = int(np.clip(zi, 0, shape[2] - 1))

    n_slices     = len(z_coords_mm)
    montage_rows = int(np.ceil(n_slices / montage_cols))

    fig = plt.figure(figsize=(15, 4.6 + 2.5 * montage_rows))
    fig.patch.set_facecolor("#111111")
    gs = fig.add_gridspec(2, 1, height_ratios=[4.6, 2.5 * montage_rows], hspace=0.1)

    # -- Top: orthogonal view through the peak voxel --
    # (radiological convention on coronal/axial -- see _plot_ortho())
    gs_top = gs[0].subgridspec(1, 3, wspace=0.04)
    ortho_slices = [
        (bg_norm[xi, :, :].T, stat_rgba[xi, :, :, :].transpose(1, 0, 2),
         "Sagittal", f"x = {cut_mm[0]:+.0f} mm", False),
        (bg_norm[::-1, yi, :].T, stat_rgba[::-1, yi, :, :].transpose(1, 0, 2),
         "Coronal",  f"y = {cut_mm[1]:+.0f} mm", True),
        (bg_norm[::-1, :, zi].T, stat_rgba[::-1, :, zi, :].transpose(1, 0, 2),
         "Axial",    f"z = {cut_mm[2]:+.0f} mm", True),
    ]
    for i, (bg_sl, stat_sl, plane, coord_lbl, needs_lr) in enumerate(ortho_slices):
        ax = fig.add_subplot(gs_top[i])
        ax.set_facecolor("black")
        ax.imshow(bg_sl,   origin="lower", cmap="gray", aspect="auto",
                  interpolation="bilinear", vmin=0, vmax=1)
        ax.imshow(stat_sl, origin="lower", aspect="auto", interpolation="nearest")
        ax.set_title(f"{plane}\n{coord_lbl}", color="white", fontsize=8, pad=3)
        ax.axis("off")
        if needs_lr:
            _add_lr_labels(ax)

    # -- Bottom: axial lightbox, fixed MNI z-levels, radiological convention --
    gs_bot = gs[1].subgridspec(montage_rows, montage_cols, wspace=0.02, hspace=0.3)
    for i, z_mm in enumerate(z_coords_mm):
        ax = fig.add_subplot(gs_bot[i // montage_cols, i % montage_cols])
        ax.set_facecolor("black")
        _, _, zi_m = _mm_to_vox(bg_img.affine, 0.0, 0.0, z_mm)
        if 0 <= zi_m < shape[2]:
            ax.imshow(bg_norm[::-1, :, zi_m].T, origin="lower", cmap="gray", aspect="auto",
                      interpolation="bilinear", vmin=0, vmax=1)
            ax.imshow(stat_rgba[::-1, :, zi_m, :].transpose(1, 0, 2), origin="lower",
                      aspect="auto", interpolation="nearest")
        ax.set_title(f"z={z_mm:+.0f}", color="white", fontsize=6.5, pad=2)
        ax.axis("off")
        _add_lr_labels(ax)

    cbar_norm = mcolors.Normalize(0, vmax) if cmap != _CMAP else mcolors.Normalize(-vmax, vmax)
    sm  = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=cbar_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=fig.axes, orientation="vertical",
                        fraction=0.015, pad=0.02, shrink=0.7)
    cbar.ax.tick_params(colors="white", labelsize=7)

    fig.suptitle(title, color="white", fontsize=10, y=0.995, va="bottom")
    return fig


def _fig_to_buf(fig: plt.Figure) -> io.BytesIO:
    buf = io.BytesIO()
    fig.savefig(buf, dpi=130, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    return buf


# ─────────────────────────────────────────────────────────────────────────────
# PDF page builders
# ─────────────────────────────────────────────────────────────────────────────

def _add_title_page(pdf: PdfPages, patient: str, modes: list[str]) -> None:
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    kw = dict(transform=ax.transAxes, ha="center", va="center")
    ax.text(0.5, 0.78, "Normative Connectivity Report",
            fontsize=22, fontweight="bold", **kw)
    ax.text(0.5, 0.67, f"Patient:  {patient}",
            fontsize=15, **kw)
    ax.text(0.5, 0.59, f"Modes:  {',  '.join(modes)}",
            fontsize=12, color="#444444", **kw)
    ax.text(0.5, 0.52, f"Generated:  {date.today().isoformat()}",
            fontsize=10, color="#666666", **kw)

    legend = (
        "HV mean FC map  — group average functional connectivity (Fisher Z)\n"
        "Patient FC map   — patient average over runs (Fisher Z)\n"
        "Z-score map       — deviation from HV mean in units of HV SD\n"
        "LOO sig map       — voxels exceeding LOO-calibrated empirical threshold\n"
        "Colormap           — blue = decreased FC,  red = increased FC\n"
        "Background         — patient's own T1w in MNI space"
    )
    ax.text(0.5, 0.32, legend,
            fontsize=10, color="#222222", linespacing=1.9,
            bbox=dict(facecolor="#f2f2f2", edgecolor="#cccccc",
                      boxstyle="round,pad=0.6"),
            **kw)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _add_section_divider(pdf: PdfPages, section_title: str) -> None:
    fig = plt.figure(figsize=(11, 2))
    fig.patch.set_facecolor("#2c5f8a")
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(0.5, 0.5, section_title, color="white",
            fontsize=16, fontweight="bold",
            ha="center", va="center", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _add_maps_page(
    pdf:          PdfPages,
    patient:      str,
    label:        str,
    type_str:     str,
    paths:        dict[str, Path],
    bg_img:       nib.Nifti1Image,
    fc_cmap:      str   = _CMAP,
    fc_threshold: float = 0.1,
) -> None:
    """
    One PDF page per seed/RSN — 2×2 grid:
      [0,0] HV group mean FC   |  [0,1] Patient FC (avg over runs)
      [1,0] Patient z-score    |  [1,1] Patient LOO significant
    Cut coordinates derived from peak of patient z-score map.
    """
    def _load(key) -> nib.Nifti1Image | None:
        p = paths.get(key)
        if p and p.exists():
            return nib.load(str(p))
        return None

    mean_img  = _load("mean")
    pt_fc_img = _load("pt_fc")
    z_img     = _load("zscore")
    loo_img   = _load("loosig")

    if z_img is None and mean_img is None:
        if pt_fc_img is not None:
            log.warning("    %s '%s': FC map exists but no z-score/mean — run step4 first",
                        type_str, label)
        else:
            log.warning("    %s '%s': no maps found — run step1/step2 and step4 first",
                        type_str, label)
        return

    cut_mm = _peak_mm(z_img, bg_img) if z_img is not None else (0.0, 0.0, 0.0)

    def _render(img, title, thr, vmax, cmap=_CMAP):
        if img is None:
            return None
        data = img.get_fdata(dtype=np.float32)
        if np.abs(data).max() < 1e-6:
            return None
        try:
            f = _plot_ortho(img, threshold=thr, vmax=vmax,
                            title=title, cut_mm=cut_mm, bg_img=bg_img, cmap=cmap)
            return _fig_to_buf(f)
        except Exception as exc:
            log.warning("    Render failed for '%s' %s: %s", label, title, exc)
            return None

    buf_mean  = _render(mean_img,  "HV Group Mean FC",             thr=fc_threshold, vmax=None, cmap=fc_cmap)
    buf_pt_fc = _render(pt_fc_img, f"{patient} — FC (avg runs)",   thr=fc_threshold, vmax=None, cmap=fc_cmap)
    buf_z     = _render(z_img,     f"{patient} — Z-score",         thr=1.0,  vmax=None)
    buf_loo   = _render(loo_img,   f"{patient} — LOO significant",  thr=0.1,  vmax=None)

    fig = plt.figure(figsize=(16, 10))
    fig.patch.set_facecolor("white")
    fig.suptitle(
        f"{type_str.upper()}:  {label}      Patient: {patient}",
        fontsize=13, fontweight="bold", y=0.995,
    )

    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        hspace=0.22, wspace=0.04,
        left=0.01, right=0.99, top=0.96, bottom=0.01,
    )

    def _place(buf: io.BytesIO | None, gs_slice, label_str: str) -> None:
        ax = fig.add_subplot(gs_slice)
        ax.axis("off")
        if buf is not None:
            img_arr = plt.imread(buf)
            ax.imshow(img_arr)
            ax.set_title(label_str, fontsize=9, pad=4, loc="left")
        else:
            ax.set_facecolor("#f0f0f0")
            ax.text(0.5, 0.5, f"{label_str}\n(not available)",
                    ha="center", va="center", transform=ax.transAxes,
                    color="#aaaaaa", fontsize=9)

    _place(buf_mean,  gs[0, 0], "HV Group Mean FC")
    _place(buf_pt_fc, gs[0, 1], f"{patient}  FC (avg over runs)")
    _place(buf_z,     gs[1, 0], f"{patient}  Z-score (uncorrected)")
    _place(buf_loo,   gs[1, 1], f"{patient}  LOO significant")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _add_table_page(
    pdf:          PdfPages,
    patient:      str,
    profile_name: str | None,
    has_rsn:      bool,
) -> None:
    """Summary table from ROI z-score CSVs."""
    frames: list[pd.DataFrame] = []

    if profile_name:
        csv = COMP_DIR / f"sub-{patient}_profile-{profile_name}_roi_zscores.csv"
        if csv.exists():
            df = pd.read_csv(csv)
            df.insert(0, "mode", "seed")
            frames.append(df)

    if has_rsn:
        csv = COMP_DIR / f"sub-{patient}_rsn_roi_zscores.csv"
        if csv.exists():
            df = pd.read_csv(csv)
            df.insert(0, "mode", "RSN")
            frames.append(df)

    if not frames:
        log.warning("  No ROI CSVs found for %s — skipping table", patient)
        return

    combined = pd.concat(frames, ignore_index=True)

    # Coalesce seed / rsn columns into a single 'network' column
    if "seed" in combined.columns and "rsn" in combined.columns:
        combined["network"] = combined["seed"].fillna(combined["rsn"])
        combined = combined.drop(columns=["seed", "rsn"])
    elif "seed" in combined.columns:
        combined = combined.rename(columns={"seed": "network"})
    elif "rsn" in combined.columns:
        combined = combined.rename(columns={"rsn": "network"})

    col_order = [
        "mode", "network", "patient_mean_FC",
        "ROI_z_score", "ROI_t_CH", "p_value", "significant_p05",
        "t_critical", "df", "method", "n_HV", "interpretation",
        "patient_network_mm3", "HV_mean_network_mm3",
        "network_extent_z", "network_extent_p", "network_extent_sig",
        "fdr_sig_mm3", "loo_sig_mm3",
    ]
    col_order  = [c for c in col_order if c in combined.columns]
    display_df = combined[col_order].copy()

    for c in ["patient_mean_FC", "ROI_z_score", "ROI_t_CH", "p_value", "t_critical"]:
        if c in display_df:
            display_df[c] = display_df[c].round(3)

    display_df["_mode_order"] = display_df["mode"].map({"seed": 0, "RSN": 1}).fillna(2)
    display_df["_abs_z"]      = display_df["ROI_z_score"].abs()
    display_df = (display_df
                  .sort_values(["_mode_order", "_abs_z"], ascending=[True, False])
                  .drop(columns=["_mode_order", "_abs_z"]))

    n_rows = len(display_df)
    fig_h  = max(5.0, 0.45 * n_rows + 3.0)
    fig    = plt.figure(figsize=(18, fig_h))
    fig.patch.set_facecolor("white")
    fig.suptitle(
        f"ROI Summary Table — Patient: {patient}",
        fontsize=13, fontweight="bold", y=0.98,
    )

    ax = fig.add_axes([0.01, 0.02, 0.98, 0.90])
    ax.axis("off")

    col_labels = [c.replace("_", " ").title() for c in display_df.columns]
    cell_text  = [[str(v) for v in row] for _, row in display_df.iterrows()]

    tbl = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc="center",
        loc="upper center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.auto_set_column_width(col=list(range(len(col_labels))))

    header_color = "#2c5f8a"
    for j in range(len(col_labels)):
        cell = tbl[(0, j)]
        cell.set_facecolor(header_color)
        cell.set_text_props(color="white", fontweight="bold")

    sig_col = col_labels.index("Significant P05") if "Significant P05" in col_labels else None
    for i in range(1, n_rows + 1):
        is_sig = (
            sig_col is not None and
            str(display_df.iloc[i - 1].get("significant_p05", "False"))
            in ("True", "true", "1")
        )
        for j in range(len(col_labels)):
            cell = tbl[(i, j)]
            if is_sig:
                cell.set_facecolor("#fff3cd")
            elif i % 2 == 0:
                cell.set_facecolor("#f7f7f7")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Coupling heatmap page
# ─────────────────────────────────────────────────────────────────────────────

def _coupling_zscore_path(sub: str, profile: str) -> Path:
    return COMP_DIR / f"sub-{sub}_profile-{profile}_seed_rsn_coupling_zscore.csv"


def _coupling_sig_path(sub: str, profile: str) -> Path:
    return COMP_DIR / f"sub-{sub}_profile-{profile}_seed_rsn_coupling_sig.csv"


def _add_coupling_heatmap_page(
    pdf:          PdfPages,
    patient:      str,
    profile_name: str,
) -> None:
    """Seed × RSN temporal coupling heatmap with Crawford-Howell z-scores."""
    z_path   = _coupling_zscore_path(patient, profile_name)
    sig_path = _coupling_sig_path(patient, profile_name)

    if not z_path.exists():
        log.warning(
            "  Coupling heatmap: no z-score CSV for %s — run step4 --profile %s first",
            patient, profile_name,
        )
        return

    z_df   = pd.read_csv(str(z_path),   index_col=0)
    sig_df = pd.read_csv(str(sig_path), index_col=0) if sig_path.exists() else None

    n_seeds = len(z_df.index)
    n_rsns  = len(z_df.columns)

    fig_w = max(14, n_rsns * 0.65 + 3)
    fig_h = max(5,  n_seeds * 0.55 + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")

    z_vals = z_df.values.astype(np.float64)
    finite = z_vals[np.isfinite(z_vals)]
    vmax   = max(float(np.abs(finite).max()), 0.1) if len(finite) > 0 else 3.0

    im = ax.imshow(z_vals, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(n_rsns))
    ax.set_xticklabels(z_df.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_seeds))
    ax.set_yticklabels(z_df.index, fontsize=8)

    for i in range(n_seeds):
        for j in range(n_rsns):
            val = z_vals[i, j]
            if not np.isfinite(val):
                continue
            is_sig = sig_df is not None and bool(sig_df.iloc[i, j])
            text   = f"{val:.2f}" + ("*" if is_sig else "")
            # Choose text colour that is readable on the heatmap background
            norm_val       = (val + vmax) / (2.0 * vmax)
            text_color     = "white" if (norm_val < 0.25 or norm_val > 0.75) else "black"
            ax.text(j, i, text, ha="center", va="center", fontsize=6.5,
                    color=text_color, fontweight="bold" if is_sig else "normal")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Crawford-Howell Z", fontsize=9)

    ax.set_title(
        f"Seed × RSN Temporal Coupling — Patient: {patient}  (Profile: {profile_name})\n"
        "* = LOO-significant (|z| > 95th percentile of HV LOO distribution)",
        fontsize=11, fontweight="bold", pad=10,
    )
    ax.set_xlabel("RSN", fontsize=9)
    ax.set_ylabel("Seed", fontsize=9)

    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Main report orchestrator
# ─────────────────────────────────────────────────────────────────────────────

def build_report(
    patient:       str,
    profile_name:  str | None  = None,
    has_rsn:       bool        = False,
    seeds:         list[str] | None = None,
    positive_only: bool        = False,
    fc_threshold:  float       = 0.1,
    thresholded:   bool        = False,
) -> Path:
    modes: list[str] = []
    if profile_name:
        modes.append(f"Seed ({profile_name})")
    if has_rsn:
        modes.append("RSN")

    fc_cmap  = "Reds" if positive_only else _CMAP
    bg_img   = _patient_bg(patient)
    pdf_path = REPORT_DIR / f"sub-{patient}_normative_report.pdf"
    log.info("Building report → %s", pdf_path)

    with PdfPages(str(pdf_path)) as pdf:

        _add_title_page(pdf, patient, modes)

        # ── Seed pages ───────────────────────────────────────────────────────
        if profile_name and seeds:
            _add_section_divider(pdf, f"Seed-Based Analysis — Profile: {profile_name}")
            log.info("  Seed pages (%d seeds)…", len(seeds))
            for seed_name in seeds:
                paths = _seed_paths(patient, seed_name, thresholded)
                has_any = any(p.exists() for k, p in paths.items() if k != "pval")
                if not has_any:
                    log.warning("    Seed '%s': no maps — skipping", seed_name)
                    continue
                log.info("    Seed: %s", seed_name)
                _add_maps_page(pdf, patient, seed_name, "seed", paths, bg_img, fc_cmap=fc_cmap, fc_threshold=fc_threshold)

        # ── RSN pages ────────────────────────────────────────────────────────
        if has_rsn:
            _add_section_divider(pdf, "RSN Normative Analysis")
            log.info("  RSN pages (%d RSNs)…", len(_ALL_RSN_NAMES))
            for rsn_name in _ALL_RSN_NAMES:
                paths = _rsn_paths(patient, rsn_name, thresholded)
                has_any = any(p.exists() for k, p in paths.items() if k != "pval")
                if not has_any:
                    log.warning("    RSN '%s': no maps — skipping", rsn_name)
                    continue
                log.info("    RSN: %s", rsn_name)
                _add_maps_page(pdf, patient, rsn_name, "RSN", paths, bg_img, fc_cmap=fc_cmap, fc_threshold=fc_threshold)

        # ── Coupling heatmap ─────────────────────────────────────────────────
        if profile_name:
            _add_section_divider(pdf, "Seed × RSN Temporal Coupling")
            _add_coupling_heatmap_page(pdf, patient, profile_name)

        # ── Summary table ────────────────────────────────────────────────────
        _add_section_divider(pdf, "ROI Summary Table")
        _add_table_page(pdf, patient, profile_name=profile_name, has_rsn=has_rsn)

    log.info("  Saved → %s", pdf_path)
    return pdf_path


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--profile",        default=None,
                   help="Profile name used in step4 seed mode (e.g. MDD)")
    p.add_argument("--rsn",            action="store_true",
                   help="Include RSN normative pages")
    p.add_argument("--positive-only",  action="store_true",
                   help="Use Reds colormap for FC maps (step4 must have been run with --positive-only)")
    p.add_argument("--fc-threshold",   type=float, default=0.1, metavar="Z",
                   help="Fisher-Z threshold for displaying FC maps (default: 0.1)")
    p.add_argument("--thresholded",     action="store_true",
                   help="Read thresholded+cluster-filtered FC maps (must match --thresholded "
                        "used in step4)")
    p.add_argument("--patients",       nargs="+", default=PT_SUBJECTS, metavar="ID")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.profile and not args.rsn:
        print("Provide --profile, --rsn, or both.")
        sys.exit(1)

    seeds: list[str] | None = None
    if args.profile:
        try:
            seeds = get_profile(args.profile)["seeds"]
        except Exception as exc:
            log.error("Could not load profile '%s': %s", args.profile, exc)
            sys.exit(1)

    for patient in args.patients:
        build_report(
            patient,
            profile_name=args.profile,
            has_rsn=args.rsn,
            seeds=seeds,
            positive_only=args.positive_only,
            fc_threshold=args.fc_threshold,
            thresholded=args.thresholded,
        )

    log.info("Done.")


if __name__ == "__main__":
    main()
