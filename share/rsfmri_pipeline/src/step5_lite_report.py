#!/usr/bin/env python3
# Run with:  python3 step5_lite_report.py [options]
"""
Step 5 (lite) — Patient-own-data PDF Report Generator

Unlike step5_report.py, this report does NOT depend on step4's normative
comparison (HV group mean / z-score / LOO-significance). It renders the
patient's own seed-based and RSN-level connectivity estimates directly from
step1 (SBA) / step2 (RSN-FC) / step3 (masked ICA) outputs. Intended for
opt-in integration paths (e.g. the Presurgical profile) where step4 is out
of scope — its normative comparison is calibrated against a 9-subject pilot
cohort and is not yet a clinically defensible baseline.

Content:
  • One ortho-view page per seed  (patient's own SBA map, no HV comparison)
  • One ortho-view page per RSN   (patient's own RSN-FC map)
  • RSN x RSN FC matrix heatmap
  • Seed <-> RSN coupling heatmap (raw Fisher-Z, requires --profile)
  • Optional masked-ICA thumbnail pages (--include-ica), gated on step3's
    own per-RSN .done flag

Usage:
  python step5_lite_report.py --profile Presurgical --rsn --patients PT01
  python step5_lite_report.py --profile Presurgical --rsn --patients PT01 PT02 --include-ica
"""

import argparse
import base64
import io
import json
import logging
import re
import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.image import resample_to_img

sys.path.insert(0, str(Path(__file__).parent))
from config import ANALYSIS_DIR, DENOISED_DIR, SBA_DIR, RSN_FC_DIR, ICA_DIR, WBICA_DIR, get_profile
from utils import YEO17_NETWORK_NAMES

# Reused verbatim from step5_report.py — generic rendering helpers with no
# dependency on step4's normative-comparison outputs.
from step5_report import (
    _mni_bg, _patient_bg, _peak_mm, _mm_to_vox, _plot_ortho,
    _plot_ortho_and_axial_montage, _fig_to_buf, _add_section_divider,
    _add_lr_labels,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPORT_DIR = ANALYSIS_DIR / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

_ALL_RSN_NAMES = YEO17_NETWORK_NAMES + ["Language", "Auditory"]
_CMAP = "RdBu_r"


# ─────────────────────────────────────────────────────────────────────────────
# Path helpers (patient-own-data only — step1/step2/step3 outputs directly)
# ─────────────────────────────────────────────────────────────────────────────

def _seed_map_path(patient: str, seed: str) -> Path:
    # Thresholded (p<p_thresh & |r|>=min_r & cluster-extent), not the raw
    # signed map -- the report must show what actually survived significance
    # testing, not every voxel including negative FC below any real threshold.
    return SBA_DIR / f"sub-{patient}" / f"sub-{patient}_seed-{seed}_desc-sbaThresh_statmap.nii.gz"


def _rsn_map_path(patient: str, rsn: str) -> Path:
    return RSN_FC_DIR / f"sub-{patient}" / f"sub-{patient}_rsn-{rsn}_desc-sbaThresh_statmap.nii.gz"


def _rsn_matrix_path(patient: str) -> Path:
    return RSN_FC_DIR / f"sub-{patient}" / f"sub-{patient}_rsn_fc_matrix.csv"


def _coupling_path(patient: str, profile_name: str) -> Path:
    return RSN_FC_DIR / f"sub-{patient}" / f"sub-{patient}_profile-{profile_name}_seed_rsn_coupling.csv"


def _ica_rsn_dir(patient: str, rsn: str) -> Path:
    return ICA_DIR / f"sub-{patient}" / f"rsn-{rsn}"


def _pooled_runs_info(patient: str) -> list[tuple[str, str]]:
    """
    [(run_id, task_label), ...] discovered from the denoised func directory
    (e.g. "sub-PT01_task-rest_run-01_..._postproc_nilearn" -> ("run-01", "rest")).
    Used only for the title page's "N runs pooled across task labels" summary.
    """
    func_dir = DENOISED_DIR / f"sub-{patient}" / "func"
    runs: list[tuple[str, str]] = []
    if not func_dir.is_dir():
        return runs
    for d in sorted(func_dir.iterdir()):
        if not d.is_dir():
            continue
        m_run  = re.search(r"(run-\d+)", d.name)
        m_task = re.search(r"_task-([^_]+)", d.name)
        if m_run and m_task:
            runs.append((m_run.group(1), m_task.group(1)))
    return runs


# ─────────────────────────────────────────────────────────────────────────────
# PDF page builders
# ─────────────────────────────────────────────────────────────────────────────

def _add_title_page(
    pdf: PdfPages,
    patient: str,
    profile_name: str | None,
    has_rsn: bool,
    runs_info: list[tuple[str, str]],
) -> None:
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    kw = dict(transform=ax.transAxes, ha="center", va="center")

    modes: list[str] = []
    if profile_name:
        modes.append(f"Seed ({profile_name})")
    if has_rsn:
        modes.append("RSN")

    ax.text(0.5, 0.86, "rsfMRI Network Connectivity Report", fontsize=22, fontweight="bold", **kw)
    ax.text(0.5, 0.78, f"Patient:  {patient}", fontsize=15, **kw)
    ax.text(0.5, 0.71, f"Modes:  {',  '.join(modes)}", fontsize=12, color="#444444", **kw)
    ax.text(0.5, 0.65, f"Generated:  {date.today().isoformat()}", fontsize=10, color="#666666", **kw)

    task_labels = sorted({t for _, t in runs_info}) or ["(none found)"]
    pooling_txt = (
        f"This report pools {len(runs_info)} BOLD run(s) across task label(s):\n"
        f"{', '.join(task_labels)}.\n\n"
        "This is a patient-own-data connectivity estimate — there is NO\n"
        "normative / healthy-volunteer comparison in this report."
    )
    if profile_name == "Presurgical":
        pooling_txt += (
            "\n\nPresurgical profile: pooling task-fMRI runs alongside resting-state\n"
            "is deliberate for motor/language eloquent-cortex mapping, and includes\n"
            "runs where the task paradigm could not be reliably performed — both\n"
            "intrinsic connectivity and task-driven co-activation inform localisation\n"
            "here, not classical resting-state FC alone."
        )

    ax.text(0.5, 0.38, pooling_txt,
            fontsize=10, color="#222222", linespacing=1.8,
            bbox=dict(facecolor="#f2f2f2", edgecolor="#cccccc", boxstyle="round,pad=0.6"),
            **kw)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _read_provenance_caption(nii_path: Path) -> str | None:
    """
    One-line "where did this map come from" summary read from the .json
    sidecar write_provenance_json() writes next to a statmap -- answers the
    question directly on the report page instead of only in a file on disk.
    """
    name = nii_path.name
    if name.endswith(".nii.gz"):
        name = name[: -len(".nii.gz")]
    json_path = nii_path.parent / f"{name}.json"
    if not json_path.exists():
        return None
    try:
        info = json.loads(json_path.read_text())
    except Exception:
        return None
    runs = ", ".join(info.get("runs_used", [])) or "?"
    return f"run-type: {info.get('run_type', '?')}  |  dof: {info.get('dof', '?')}  |  runs pooled: {runs}"


def _add_own_map_page(
    pdf:          PdfPages,
    patient:      str,
    label:        str,
    type_str:     str,
    path:         Path,
    bg_img:       nib.Nifti1Image,
    fc_threshold: float = 0.1,
    cmap:         str   = _CMAP,
) -> None:
    """One page: patient's own FC map only (no HV comparison), ortho views."""
    if not path.exists():
        log.warning("    %s '%s': no map found (%s) — skipping", type_str, label, path.name)
        return

    img  = nib.load(str(path))
    data = img.get_fdata(dtype=np.float32)
    if np.abs(data).max() < 1e-6:
        log.warning("    %s '%s': map is all-zero — skipping", type_str, label)
        return

    caption = _read_provenance_caption(path)
    title = f"{type_str.upper()}: {label}      Patient: {patient}"
    if caption:
        title += f"\n{caption}"

    cut_mm = _peak_mm(img, bg_img)
    try:
        fig = _plot_ortho_and_axial_montage(
            img, threshold=fc_threshold, vmax=None,
            title=title,
            cut_mm=cut_mm, bg_img=bg_img, cmap=cmap,
        )
    except Exception as exc:
        log.warning("    Render failed for '%s' %s: %s", label, type_str, exc)
        return

    pdf.savefig(fig, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def _add_heatmap_page(
    pdf:      PdfPages,
    csv_path: Path,
    title:    str,
    xlabel:   str,
    ylabel:   str,
) -> None:
    """
    Generic Fisher-Z matrix heatmap page (imshow + per-cell annotation +
    colorbar) — rendering structure carried over from step5_report.py's
    _add_coupling_heatmap_page, applied to the raw values step2 already
    writes directly (no step4 z-scores/significance markers involved).
    """
    if not csv_path.exists():
        log.warning("  Heatmap: no CSV found (%s) — skipping", csv_path.name)
        return

    df = pd.read_csv(str(csv_path), index_col=0)
    n_rows, n_cols = df.shape

    fig_w = max(10, n_cols * 0.6 + 3)
    fig_h = max(5,  n_rows * 0.5 + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")

    vals   = df.values.astype(np.float64)
    finite = vals[np.isfinite(vals)]
    vmax   = max(float(np.abs(finite).max()), 0.1) if len(finite) > 0 else 1.0

    im = ax.imshow(vals, cmap=_CMAP, aspect="auto", vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(df.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(df.index, fontsize=8)

    for i in range(n_rows):
        for j in range(n_cols):
            val = vals[i, j]
            if not np.isfinite(val):
                continue
            norm_val   = (val + vmax) / (2.0 * vmax)
            text_color = "white" if (norm_val < 0.25 or norm_val > 0.75) else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=6.5, color=text_color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Fisher Z", fontsize=9)

    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)

    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _ica_combined_map_path(patient: str, rsn: str, thresholded: bool = True) -> Path:
    suffix = "ica_zstatThresh" if thresholded else "ica_zstat"
    return _ica_rsn_dir(patient, rsn) / f"sub-{patient}_rsn-{rsn}_{suffix}.nii.gz"


def _add_ica_thumbnail_page(
    pdf:              PdfPages,
    patient:          str,
    rsn:              str,
    bg_img:           nib.Nifti1Image,
    slice_offsets_mm: tuple[float, ...] = (-12.0, 0.0, 12.0),
) -> None:
    """
    One page per RSN: the single rest/task-combined, FDR+cluster-extent
    thresholded masked-ICA map (see step3_masked_ica.py -- one Stouffer-
    combined component per network, not a dump of every raw IC), rendered
    as a small axial mini-lightbox centred on its own peak |z| voxel.
    Radiological convention with L/R labels, consistent with the rest of
    the report.
    """
    out_dir = _ica_rsn_dir(patient, rsn)
    if not (out_dir / ".done").exists():
        return  # step3 never completed for this RSN — nothing to show

    map_path = _ica_combined_map_path(patient, rsn, thresholded=True)
    if not map_path.exists():
        log.warning("  ICA thumbnails: no combined map for RSN '%s' — skipping", rsn)
        return

    bg_data = bg_img.get_fdata(dtype=np.float32)
    shape   = bg_data.shape[:3]
    bg_max  = bg_data.max()
    bg_norm = bg_data / bg_max if bg_max > 0 else bg_data

    n_slices = len(slice_offsets_mm)
    fig, axes = plt.subplots(1, n_slices, figsize=(2.6 * n_slices, 2.8),
                             squeeze=False, gridspec_kw={"wspace": 0.03})
    fig.patch.set_facecolor("#111111")
    fig.suptitle(f"Masked ICA — RSN: {rsn}      Patient: {patient}",
                 color="white", fontsize=11)

    try:
        comp_img  = nib.load(str(map_path))
        comp_r    = resample_to_img(comp_img, bg_img, interpolation="linear", force_resample=True)
        comp_data = comp_r.get_fdata(dtype=np.float32)
    except Exception as exc:
        log.warning("    ICA thumbnail render failed for RSN '%s': %s", rsn, exc)
        for ax in axes[0]:
            ax.axis("off")
        pdf.savefig(fig, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)
        return

    comp_data = np.where(comp_data != 0, comp_data, np.nan)
    finite = comp_data[np.isfinite(comp_data)]
    vmax = float(np.max(np.abs(finite))) if finite.size > 0 else 1.0
    vmax = max(vmax, 0.01)

    # Centre the mini-lightbox on the thresholded map's own peak |z| voxel.
    _, _, peak_z = _peak_mm(comp_img, bg_img)

    for col, offset in enumerate(slice_offsets_mm):
        ax = axes[0][col]
        ax.set_facecolor("black")
        _, _, zi = _mm_to_vox(bg_img.affine, 0.0, 0.0, peak_z + offset)
        zi = int(np.clip(zi, 0, shape[2] - 1))
        ax.imshow(bg_norm[::-1, :, zi].T, origin="lower", cmap="gray",
                  vmin=0, vmax=1, aspect="auto")
        ax.imshow(comp_data[::-1, :, zi].T, origin="lower", cmap=_CMAP,
                  vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_title(f"z={peak_z + offset:+.0f}", color="white", fontsize=8)
        ax.axis("off")
        _add_lr_labels(ax)

    pdf.savefig(fig, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def _wbica_dir(patient: str) -> Path:
    return WBICA_DIR / f"sub-{patient}"


def _add_wbica_page(pdf: PdfPages, patient: str, run_type: str,
                    bg_img: nib.Nifti1Image, z_thr: float = 2.3) -> None:
    """
    One contact-sheet page per run-type: every whole-brain-ICA component
    (see step3b_whole_brain_ica.py) as a small single-slice thumbnail at its
    own peak |z| voxel -- exploratory/QC, so a compact grid of every
    component beats a full multi-page render per component.
    """
    out_dir = _wbica_dir(patient)
    if not (out_dir / f".done_{run_type}").exists():
        return

    comp_paths = sorted(out_dir.glob(f"sub-{patient}_wbica_{run_type}_IC-*_zstat.nii.gz"))
    if not comp_paths:
        log.warning("  Whole-brain ICA thumbnails: no components for run-type '%s' — skipping",
                    run_type)
        return

    bg_data = bg_img.get_fdata(dtype=np.float32)
    shape   = bg_data.shape[:3]
    bg_max  = bg_data.max()
    bg_norm = bg_data / bg_max if bg_max > 0 else bg_data

    n = len(comp_paths)
    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.4 * ncols, 2.6 * nrows),
                             squeeze=False, gridspec_kw={"wspace": 0.05, "hspace": 0.3})
    fig.patch.set_facecolor("#111111")
    fig.suptitle(f"Whole-brain ICA — {run_type}      Patient: {patient}      "
                f"({n} components, exploratory)",
                color="white", fontsize=11)

    for idx, comp_path in enumerate(comp_paths):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        ax.set_facecolor("black")
        ic_label = comp_path.stem.split("_IC-")[1].split("_")[0]
        try:
            comp_img  = nib.load(str(comp_path))
            comp_r    = resample_to_img(comp_img, bg_img, interpolation="linear",
                                        force_resample=True)
            comp_data = comp_r.get_fdata(dtype=np.float32)
        except Exception as exc:
            log.warning("    Whole-brain ICA thumbnail failed for IC-%s: %s", ic_label, exc)
            ax.axis("off")
            continue

        comp_data = np.where(np.abs(comp_data) >= z_thr, comp_data, np.nan)
        finite = comp_data[np.isfinite(comp_data)]
        vmax = float(np.max(np.abs(finite))) if finite.size > 0 else 1.0
        vmax = max(vmax, 0.01)

        _, _, peak_z = _peak_mm(comp_img, bg_img)
        _, _, zi = _mm_to_vox(bg_img.affine, 0.0, 0.0, peak_z)
        zi = int(np.clip(zi, 0, shape[2] - 1))

        ax.imshow(bg_norm[::-1, :, zi].T, origin="lower", cmap="gray",
                  vmin=0, vmax=1, aspect="auto")
        ax.imshow(comp_data[::-1, :, zi].T, origin="lower", cmap=_CMAP,
                  vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_title(f"IC {ic_label}  z={peak_z:+.0f}", color="white", fontsize=8)
        ax.axis("off")
        _add_lr_labels(ax)

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].axis("off")

    pdf.savefig(fig, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# HTML output (same content as the PDF, reusing the exact same figures)
# ─────────────────────────────────────────────────────────────────────────────

class _MultiSink:
    """
    Duck-types PdfPages' .savefig() so every existing _add_*_page() call
    site works completely unchanged -- forwards each figure to the real
    PdfPages AND snapshots a PNG for the HTML report in the same call,
    before the caller's subsequent plt.close(fig). Both outputs therefore
    come from the exact same rendering pass; there's no separate code path
    that could let the PDF and HTML drift out of sync with each other.
    """
    def __init__(self, pdf_pages: PdfPages):
        self._pdf = pdf_pages
        self.png_b64_pages: list[str] = []

    def savefig(self, fig: plt.Figure, **kwargs) -> None:
        self._pdf.savefig(fig, **kwargs)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=130,
                    bbox_inches=kwargs.get("bbox_inches"),
                    facecolor=fig.get_facecolor())
        buf.seek(0)
        self.png_b64_pages.append(base64.b64encode(buf.read()).decode("ascii"))


def _write_html_report(html_path: Path, png_b64_pages: list[str], patient: str) -> None:
    """
    Self-contained HTML report: the same sequence of rendered pages as the
    PDF, each embedded as a base64 PNG -- no external image files, so the
    single .html file is just as portable/standalone as the PDF (openable
    offline in any browser, nothing to lose track of).
    """
    imgs_html = "\n".join(
        f'  <img src="data:image/png;base64,{b64}" alt="report page {i + 1}">'
        for i, b64 in enumerate(png_b64_pages)
    )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>rsfMRI Network Connectivity Report — sub-{patient}</title>
<style>
  body {{
    background: #1a1a1a; margin: 0; padding: 24px 0;
    display: flex; flex-direction: column; align-items: center; gap: 18px;
  }}
  img {{ max-width: min(1400px, 96vw); box-shadow: 0 2px 10px rgba(0,0,0,0.5); }}
</style>
</head>
<body>
{imgs_html}
</body>
</html>
"""
    html_path.write_text(html, encoding="utf-8")
    log.info("  Saved → %s", html_path)


# ─────────────────────────────────────────────────────────────────────────────
# Main report orchestrator
# ─────────────────────────────────────────────────────────────────────────────

def build_report(
    patient:       str,
    profile_name:  str | None = None,
    has_rsn:       bool       = False,
    seeds:         list[str] | None = None,
    include_ica:   bool       = False,
    include_wbica: bool       = False,
    fc_threshold:  float      = 0.1,
) -> Path:
    bg_img    = _patient_bg(patient)
    runs_info = _pooled_runs_info(patient)
    pdf_path  = REPORT_DIR / f"sub-{patient}_rsfmri_networks_report.pdf"
    html_path = REPORT_DIR / f"sub-{patient}_rsfmri_networks_report.html"
    log.info("Building lite report → %s (+ .html)", pdf_path)

    with PdfPages(str(pdf_path)) as real_pdf:
        sink = _MultiSink(real_pdf)

        _add_title_page(sink, patient, profile_name, has_rsn, runs_info)

        # ── Seed pages ───────────────────────────────────────────────────────
        if profile_name and seeds:
            _add_section_divider(sink, f"Seed-Based Connectivity — Profile: {profile_name}")
            log.info("  Seed pages (%d seeds)…", len(seeds))
            for seed_name in seeds:
                _add_own_map_page(sink, patient, seed_name, "seed",
                                   _seed_map_path(patient, seed_name), bg_img,
                                   fc_threshold=fc_threshold)

        # ── RSN pages ────────────────────────────────────────────────────────
        if has_rsn:
            _add_section_divider(sink, "RSN Connectivity")
            log.info("  RSN pages (%d RSNs)…", len(_ALL_RSN_NAMES))
            for rsn_name in _ALL_RSN_NAMES:
                _add_own_map_page(sink, patient, rsn_name, "RSN",
                                   _rsn_map_path(patient, rsn_name), bg_img,
                                   fc_threshold=fc_threshold)

            _add_section_divider(sink, "RSN x RSN FC Matrix")
            _add_heatmap_page(sink, _rsn_matrix_path(patient),
                               f"RSN x RSN FC Matrix — Patient: {patient}",
                               xlabel="RSN", ylabel="RSN")

        # ── Seed x RSN coupling ──────────────────────────────────────────────
        if profile_name and has_rsn:
            _add_section_divider(sink, "Seed x RSN Coupling")
            _add_heatmap_page(sink, _coupling_path(patient, profile_name),
                               f"Seed x RSN Coupling — Patient: {patient}  (Profile: {profile_name})",
                               xlabel="RSN", ylabel="Seed")

        # ── Masked ICA (opt-in) ──────────────────────────────────────────────
        if include_ica and has_rsn:
            _add_section_divider(sink, "Masked ICA (opt-in)")
            log.info("  ICA thumbnail pages…")
            for rsn_name in _ALL_RSN_NAMES:
                _add_ica_thumbnail_page(sink, patient, rsn_name, bg_img)

        # ── Whole-brain ICA (opt-in) ─────────────────────────────────────────
        if include_wbica:
            _add_section_divider(sink, "Whole-brain ICA (opt-in, exploratory)")
            log.info("  Whole-brain ICA pages…")
            for run_type in ("rest", "task"):
                _add_wbica_page(sink, patient, run_type, bg_img)

    log.info("  Saved → %s", pdf_path)
    _write_html_report(html_path, sink.png_b64_pages, patient)
    return pdf_path


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--profile",       default=None,
                   help="Condition profile for seed pages + seed<->RSN coupling (e.g. Presurgical)")
    p.add_argument("--rsn",           action="store_true",
                   help="Include RSN pages and the RSN x RSN FC matrix")
    p.add_argument("--include-ica",   action="store_true",
                   help="Include masked-ICA thumbnail pages (only for RSNs where step3 completed)")
    p.add_argument("--include-wbica", action="store_true",
                   help="Include whole-brain ICA contact-sheet pages (step3b output)")
    p.add_argument("--fc-threshold",  type=float, default=0.1, metavar="Z",
                   help="Fisher-Z threshold for displaying FC maps (default: 0.1)")
    p.add_argument("--patients",      nargs="+", required=True, metavar="ID")
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
            include_ica=args.include_ica,
            include_wbica=args.include_wbica,
            fc_threshold=args.fc_threshold,
        )

    log.info("Done.")


if __name__ == "__main__":
    main()
