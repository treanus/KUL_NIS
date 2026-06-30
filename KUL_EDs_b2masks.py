#!/usr/bin/env python3
# KUL_EDs_b2masks.py — Swiss-knife distance tool for binary NIfTI masks
#
# Author: Ahmed Radwan, UZ/KU Leuven — ahmed.radwan@kuleuven.be, radwanphd@gmail.com
# v1.1 — 2026-06
#
# Supported metrics : min Euclidean · Hausdorff · 95th-pct Hausdorff · ASSD
# Modes             : one-to-one · one-to-many · all-pairs
# Multi-class       : --classes 1 2 --class-labels core edema
# Parallel          : ProcessPoolExecutor for multi-pair; ThreadPoolExecutor for EDT maps
# Output            : per-pair subdirectory · PNG snapshot · HTML report · CSV


import os
import sys
import re
import csv
import base64
import io
import argparse
import itertools
import datetime
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.spatial.distance import cdist
from scipy.ndimage import (
    distance_transform_edt,
    binary_erosion,
    binary_dilation,
    center_of_mass,
)

# Per-worker cache: populated by _init_worker for one-to-many parallel mode.
# Stores pre-computed mask A data so EDT-A is only paid once per worker process.
_WORKER_CACHE: dict = {}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _subject_tag(path: str) -> str:
    """Extract BIDS subject ID from path, or fall back to filename stem."""
    m = re.search(r'sub-([^/_]+)', path)
    return '_' + m.group(1) if m else '_' + Path(path).name.split('.')[0]


def _pair_tag(mask_a_path: str, mask_b_path: str, class_label: str = '') -> str:
    """Unique subdirectory name for one pair (with optional class suffix)."""
    tag_a = _subject_tag(mask_a_path)
    if class_label:
        tag_a = tag_a + '_' + class_label
    stem_b = Path(mask_b_path).name.split('.')[0]
    return f'{tag_a}_vs_{stem_b}'


def _voxel_size(affine: np.ndarray) -> np.ndarray:
    """Voxel dimensions (mm) — correct for both diagonal and rotated affines."""
    return nib.affines.voxel_sizes(affine).astype(np.float64)


def _dilated_marker(shape: tuple, vv, iterations: int = 5) -> np.ndarray:
    """Single-voxel uint16 map dilated for visibility; center voxel set to 10."""
    arr = np.zeros(shape, np.uint16)
    arr[int(vv[0]), int(vv[1]), int(vv[2])] = 1
    dil = np.uint16(binary_dilation(arr.astype(bool), iterations=iterations))
    dil[int(vv[0]), int(vv[1]), int(vv[2])] = 10
    return dil


def _init_worker(mask_a_path: str, class_val, use_surf: bool):
    """
    Initialiser for ProcessPoolExecutor workers in one-to-many mode.
    Pre-computes all mask-A data (including EDT) once per worker process so
    that each subsequent pair only needs to load and process mask B.
    """
    try:
        img = nib.load(mask_a_path)
        aff = img.affine
        vox = _voxel_size(aff)
        raw = img.get_fdata()
        if class_val is not None:
            im = np.uint16(np.round(raw).astype(np.int64) == int(class_val))
        else:
            im = np.uint16(raw)
        mb = im.astype(bool)
        er = np.uint16(binary_erosion(mb))
        ol = np.uint16(np.abs(np.subtract(im, er)))
        cog = np.array(center_of_mass(im))
        cog_xyz = nib.affines.apply_affine(aff, cog)
        edt = _edt_full(mb, vox)
        _WORKER_CACHE['a'] = dict(
            path=mask_a_path, class_val=class_val,
            aff=aff, vox=vox, im=im, mb=mb, er=er, ol=ol,
            cog=cog, cog_xyz=cog_xyz, edt=edt,
        )
    except Exception as exc:
        _WORKER_CACHE['a'] = {'error': str(exc)}


def _edt_full(mask_bool: np.ndarray, vox_size: np.ndarray) -> np.ndarray:
    """EDT: distance (mm) from every voxel to the nearest True voxel in mask_bool."""
    return distance_transform_edt(~mask_bool, sampling=vox_size).astype(np.float32)


def _make_vtk_snapshot(
    bg_image: str,
    mask_a_path: str,
    mask_b_path: str,
    snapshot_path: str,
    focus_vox: tuple,       # (i, j, k) voxel index — used to set camera focal point
) -> str:
    """
    True off-screen 3-panel volume render (superior | lateral | posterior) using VTK.
    Background rendered as a glass-brain volume (low opacity) so both mask surfaces
    are visible through it.  No display or xvfb required — VTK uses software ray casting.
    Returns base64-encoded PNG on success, '' if vtk or PIL are unavailable.
    """
    try:
        import vtk
        from PIL import Image as _PIL
    except ImportError:
        return ''

    def _import(data: np.ndarray, vox: np.ndarray) -> 'vtk.vtkImageImport':
        imp = vtk.vtkImageImport()
        arr = data.flatten(order='F')
        imp.CopyImportVoidPointer(arr, arr.nbytes)
        imp.SetDataScalarTypeToFloat()
        imp.SetNumberOfScalarComponents(1)
        sx, sy, sz = data.shape
        imp.SetDataExtent(0, sx-1, 0, sy-1, 0, sz-1)
        imp.SetWholeExtent(0, sx-1, 0, sy-1, 0, sz-1)
        imp.SetDataSpacing(*vox)
        imp.Update()
        return imp

    def _glass_vol(imp, data):
        vmax = float(np.percentile(data[data > 0], 99)) if np.any(data > 0) else 1.0
        otf = vtk.vtkPiecewiseFunction()
        otf.AddPoint(0,          0.00)
        otf.AddPoint(vmax * 0.1, 0.00)
        otf.AddPoint(vmax * 0.4, 0.01)
        otf.AddPoint(vmax,       0.015)
        ctf = vtk.vtkColorTransferFunction()
        ctf.AddRGBPoint(0,          0.0, 0.0, 0.0)
        ctf.AddRGBPoint(vmax * 0.3, 0.7, 0.7, 0.8)
        ctf.AddRGBPoint(vmax,       1.0, 1.0, 1.0)
        vp = vtk.vtkVolumeProperty()
        vp.SetScalarOpacity(otf); vp.SetColor(ctf)
        vp.ShadeOff(); vp.SetInterpolationTypeToLinear()
        mapper = vtk.vtkSmartVolumeMapper()
        mapper.SetInputConnection(imp.GetOutputPort())
        vol = vtk.vtkVolume()
        vol.SetMapper(mapper); vol.SetProperty(vp)
        return vol

    def _surface(mask_path, rgb, opacity):
        img  = nib.load(mask_path)
        data = img.get_fdata().astype(np.float32)
        vox  = nib.affines.voxel_sizes(img.affine).astype(float)
        imp  = _import(data, vox)
        mc   = vtk.vtkMarchingCubes()
        mc.SetInputConnection(imp.GetOutputPort()); mc.SetValue(0, 0.5); mc.Update()
        sm   = vtk.vtkSmoothPolyDataFilter()
        sm.SetInputConnection(mc.GetOutputPort())
        sm.SetNumberOfIterations(20); sm.SetRelaxationFactor(0.1); sm.Update()
        pm = vtk.vtkPolyDataMapper()
        pm.SetInputConnection(sm.GetOutputPort()); pm.ScalarVisibilityOff()
        ac = vtk.vtkActor(); ac.SetMapper(pm)
        ac.GetProperty().SetColor(*rgb); ac.GetProperty().SetOpacity(opacity)
        return ac

    def _panel(vol, surf_a, surf_b, cam_offset, view_up, focal_mm, D, W=500, H=500):
        ren = vtk.vtkRenderer(); ren.SetBackground(0.05, 0.05, 0.08)
        ren.AddVolume(vol); ren.AddActor(surf_a); ren.AddActor(surf_b)
        cam = ren.GetActiveCamera()
        cam.SetFocalPoint(*focal_mm)
        cam.SetPosition(*(focal_mm + np.array(cam_offset) * D))
        cam.SetViewUp(*view_up)
        ren.ResetCamera(); ren.ResetCameraClippingRange()
        rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, H)
        rw.AddRenderer(ren); rw.Render()
        wif = vtk.vtkWindowToImageFilter(); wif.SetInput(rw); wif.Update()
        out = wif.GetOutput(); w2, h2, _ = out.GetDimensions()
        arr = np.frombuffer(out.GetPointData().GetScalars(), dtype=np.uint8).reshape(h2, w2, 3)
        return arr[::-1]   # VTK stores bottom-up

    try:
        bg    = nib.load(bg_image)
        bgd   = bg.get_fdata().astype(np.float32)
        vox   = nib.affines.voxel_sizes(bg.affine).astype(float)
        imp   = _import(bgd, vox)
        vol   = _glass_vol(imp, bgd)
        sa    = _surface(mask_a_path, (1.0, 0.00, 0.00), 0.70)
        sb    = _surface(mask_b_path, (0.00, 1.00, 1.00), 0.70)

        # Focal point in VTK image space (voxel × spacing)
        focal = np.array(focus_vox, dtype=float) * vox
        D     = float(max(np.array(bgd.shape) * vox) * 1.6)

        panels = [
            _panel(vol, sa, sb, ( 0,  0,  1), (0, 1, 0), focal, D),   # superior
            _panel(vol, sa, sb, (-1,  0,  0), (0, 0, 1), focal, D),   # left lateral
            _panel(vol, sa, sb, ( 0, -1,  0), (0, 0, 1), focal, D),   # posterior
        ]
        strip = np.concatenate(panels, axis=1)
        _PIL.fromarray(strip).save(snapshot_path)
        with open(snapshot_path, 'rb') as fh:
            return base64.b64encode(fh.read()).decode('ascii')
    except Exception:
        return ''


def _make_mrview_snapshot(
    bg_image: str,
    mask_a_path: str,
    mask_b_path: str,
    dist_a2b_path,          # str path to nii.gz or None
    snapshot_path: str,
    focus_vox: tuple,       # (i, j, k) voxel index for camera focus
    mask_b_edge_path: str = None,
) -> str:
    """
    Generate a 3-panel orthographic QC snapshot via mrview (mode 2) using the
    anatomical bg_image as context.  Always rendered off-screen via xvfb-run.
    Returns base64-encoded PNG on success, '' on failure (missing mrview/xvfb-run).
    Note: mrview volume rendering (mode 3) requires hardware OpenGL and cannot
    run through xvfb, so only orthographic mode is supported here.
    """
    import shutil
    import subprocess

    mrview_bin = shutil.which('mrview')
    if not mrview_bin:
        return ''

    # Always render off-screen via xvfb-run so mrview never pops up on the
    # user's display.  Fall back to using the existing DISPLAY only when
    # xvfb-run is not installed.
    xvfb_bin = shutil.which('xvfb-run')
    if not xvfb_bin and not (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
        return ''   # no display and no xvfb-run: cannot render

    fov  = f'{int(focus_vox[0])},{int(focus_vox[1])},{int(focus_vox[2])}'

    cmd: list = []
    if xvfb_bin:
        cmd += [xvfb_bin, '-a', '--server-args=-screen 0 1920x1080x24']

    cmd += [mrview_bin, bg_image]

    # Distance map first — rendered under masks so colours remain visible
    if dist_a2b_path and os.path.isfile(dist_a2b_path):
        cmd += [
            '-overlay.load',          dist_a2b_path,
            '-overlay.colourmap',     '5',    # viridis — distinct from red/cyan masks
            '-overlay.opacity',       '0.45',
            '-overlay.threshold_min', '0.001',
        ]

    # Mask A — red (slightly transparent)
    cmd += [
        '-overlay.load',          mask_a_path,
        '-overlay.colour',        '1,0,0',
        '-overlay.opacity',       '0.52',
        '-overlay.threshold_min', '0.5',
    ]

    # Mask B — cyan fill (low opacity) + bright edge overlay
    cmd += [
        '-overlay.load',          mask_b_path,
        '-overlay.colour',        '0,1,1',
        '-overlay.opacity',       '0.40',
        '-overlay.threshold_min', '0.5',
    ]
    if mask_b_edge_path and os.path.isfile(mask_b_edge_path):
        cmd += [
            '-overlay.load',          mask_b_edge_path,
            '-overlay.colour',        '1,1,0',
            '-overlay.opacity',       '1.0',
            '-overlay.threshold_min', '0.5',
        ]

    # -capture.grab saves to {folder}/{prefix}0001.png; rename to snapshot_path after.
    cap_dir    = os.path.dirname(snapshot_path)
    cap_prefix = '_mrvcap_'
    cap_file   = os.path.join(cap_dir, cap_prefix + '0000.png')

    cmd += [
        '-voxel', fov,
        '-mode',  '2',
        '-size',  '1200,1200',
        '-noannotations',
        '-capture.folder',  cap_dir,
        '-capture.prefix',  cap_prefix,
        '-capture.grab',
        '-exit',
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        if os.path.isfile(cap_file):
            _mrview_crop_quadrant(cap_file, snapshot_path)
            with open(snapshot_path, 'rb') as fh:
                return base64.b64encode(fh.read()).decode('ascii')
    except Exception:
        pass
    return ''


def _mrview_crop_quadrant(src: str, dst: str):
    """
    mrview orthographic mode always emits a 2×2 grid; the 4th panel (bottom-right)
    is the 3-D overview which renders black for volume data.  Detect and remove it:
    reshape to a 3-panel horizontal strip (coronal | sagittal | axial).
    Falls back to a plain rename if PIL is unavailable.
    """
    try:
        from PIL import Image
        img = Image.open(src)
        arr = np.array(img)
        h, w = arr.shape[:2]
        mh, mw = h // 2, w // 2

        # Is the bottom-right quadrant essentially black?
        br_mean = float(arr[mh:, mw:, :3].mean())
        if br_mean < 8:
            # Build a 3-panel horizontal strip: coronal (top-left) | sagittal (top-right) | axial (bottom-left)
            panel_h = mh
            panel_w = mw
            out = np.zeros((panel_h, panel_w * 3, arr.shape[2]), dtype=arr.dtype)
            out[:, :panel_w]            = arr[:mh, :mw]   # coronal
            out[:, panel_w:2*panel_w]   = arr[:mh, mw:]   # sagittal
            out[:, 2*panel_w:]          = arr[mh:, :mw]   # axial
            Image.fromarray(out).save(dst)
        else:
            os.replace(src, dst)
    except Exception:
        os.replace(src, dst)


def _make_snapshot(im1: np.ndarray, im2: np.ndarray,
                   dist_a2b,          # np.ndarray or None (when --no-maps)
                   a_vox_vv, b_vox_vv, cog1, cog2,
                   mask_a_name: str, mask_b_name: str,
                   class_label: str = '') -> str:
    """
    3-panel axial/coronal/sagittal QC snapshot centred on the closest-voxel pair.
    Returns a base64-encoded PNG string, or '' if matplotlib is unavailable.
    Fallback used when mrview is unavailable or --bg-image is not provided.
    Must be called with the Agg backend (non-interactive); safe inside subprocesses.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        return ''

    def _cl(v, lo, hi): return max(lo, min(int(v), hi))
    si = _cl(a_vox_vv[0], 0, im1.shape[0] - 1)
    sj = _cl(a_vox_vv[1], 0, im1.shape[1] - 1)
    sk = _cl(a_vox_vv[2], 0, im1.shape[2] - 1)

    dm_vmax = 1.0
    if dist_a2b is not None:
        vals = dist_a2b[dist_a2b > 0]
        if vals.size > 0:
            dm_vmax = float(np.percentile(vals, 95))

    def _panel(ax, sl_a, sl_b, sl_dm, pt_a, pt_b, c1, c2, title):
        ax.set_facecolor('black')
        rgba_b  = np.zeros((*sl_b.shape, 4), np.float32)
        rgba_a  = np.zeros((*sl_a.shape, 4), np.float32)
        rgba_ov = np.zeros((*sl_a.shape, 4), np.float32)
        mb = sl_b > 0
        sl_be = mb.astype(np.uint8) - binary_erosion(mb, iterations=3).astype(np.uint8)
        rgba_b [mb]                          = [0.00, 1.00, 1.00, 0.40]
        rgba_a [sl_a > 0]                    = [1.00, 0.00, 0.00, 0.56]
        rgba_ov[(sl_a > 0) & mb]             = [1.00, 1.00, 0.00, 0.80]
        rgba_be = np.zeros((*sl_b.shape, 4), np.float32)
        rgba_be[sl_be > 0]                   = [1.00, 1.00, 0.00, 1.00]

        ax.imshow(np.zeros_like(sl_a), cmap='gray', origin='lower', aspect='auto')
        ax.imshow(rgba_b,  origin='lower', aspect='auto')
        ax.imshow(rgba_be, origin='lower', aspect='auto')
        if sl_dm is not None:
            dm_m = np.ma.masked_where(sl_dm == 0, sl_dm)
            if dm_m.count() > 0:
                ax.imshow(dm_m, cmap='hot', origin='lower', aspect='auto',
                          alpha=0.75, vmin=0, vmax=dm_vmax)
        ax.imshow(rgba_a,  origin='lower', aspect='auto')
        ax.imshow(rgba_ov, origin='lower', aspect='auto')

        ax.plot(pt_a[0], pt_a[1], 'r*', ms=12, mec='white', mew=0.5)
        ax.plot(pt_b[0], pt_b[1], '*', ms=12, mec='white', mew=0.5, color='cyan')
        ax.plot(c1[0],   c1[1],   'r+', ms=9,  mew=1.5)
        ax.plot(c2[0],   c2[1],   '+', ms=9,  mew=1.5, color='cyan')
        ax.set_title(title, color='white', fontsize=8, pad=2)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        for sp in ax.spines.values():
            sp.set_visible(False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), facecolor='black')
    a_label = Path(mask_a_name).stem + (f' [{class_label}]' if class_label else '')
    fig.suptitle(f'{a_label}  ↔  {Path(mask_b_name).stem}', color='white', fontsize=9)

    # Axial   (fix k): im[:,:,k] → display (j=col, i=row) → .T
    _panel(axes[0],
           im1[:, :, sk].T, im2[:, :, sk].T,
           dist_a2b[:, :, sk].T if dist_a2b is not None else None,
           (int(a_vox_vv[1]), int(a_vox_vv[0])),
           (int(b_vox_vv[1]), int(b_vox_vv[0])),
           (int(cog1[1]), int(cog1[0])), (int(cog2[1]), int(cog2[0])),
           f'Axial  z={sk}')

    # Coronal (fix j): im[:,j,:] → display (k=col, i=row) → .T
    _panel(axes[1],
           im1[:, sj, :].T, im2[:, sj, :].T,
           dist_a2b[:, sj, :].T if dist_a2b is not None else None,
           (int(a_vox_vv[2]), int(a_vox_vv[0])),
           (int(b_vox_vv[2]), int(b_vox_vv[0])),
           (int(cog1[2]), int(cog1[0])), (int(cog2[2]), int(cog2[0])),
           f'Coronal  y={sj}')

    # Sagittal (fix i): im[i,:,:] → display (k=col, j=row) → .T
    _panel(axes[2],
           im1[si, :, :].T, im2[si, :, :].T,
           dist_a2b[si, :, :].T if dist_a2b is not None else None,
           (int(a_vox_vv[2]), int(a_vox_vv[1])),
           (int(b_vox_vv[2]), int(b_vox_vv[1])),
           (int(cog1[2]), int(cog1[1])), (int(cog2[2]), int(cog2[1])),
           f'Sagittal  x={si}')

    from matplotlib.lines import Line2D
    legend = [
        Patch(facecolor=(1.0, 0.0, 0.0, 0.8),
              label='Mask A' + (f' [{class_label}]' if class_label else '')),
        Patch(facecolor=(0.0, 1.0, 1.0, 0.8), label='Mask B'),
        Patch(facecolor=(1.00, 1.00, 0.00, 0.9), label='Overlap'),
    ]
    if dist_a2b is not None:
        legend.append(Patch(facecolor='orangered', label='A→B dist (mm)'))
    legend += [
        Line2D([0], [0], marker='*', color='w', markerfacecolor='red',
               markersize=9, label='Closest voxel A', linestyle='None'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='cyan',
               markersize=9, label='Closest voxel B', linestyle='None'),
        Line2D([0], [0], marker='+', color='red',
               markersize=8, markeredgewidth=1.5, label='COG A', linestyle='None'),
        Line2D([0], [0], marker='+', color='cyan',
               markersize=8, markeredgewidth=1.5, label='COG B', linestyle='None'),
    ]
    axes[2].legend(handles=legend, loc='lower right', facecolor='#1a1a1a',
                   labelcolor='white', fontsize=7, framealpha=0.9, edgecolor='#555')

    plt.tight_layout(pad=0.4)
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100, bbox_inches='tight', facecolor='black')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('ascii')


# ---------------------------------------------------------------------------
# Core per-pair computation (module-level → picklable by ProcessPoolExecutor)
# ---------------------------------------------------------------------------

def compute_pair(mask_a_path: str, mask_b_path: str, args,
                 class_val: int = None, class_label: str = '') -> dict:
    """
    Load one pair of binary (or class-extracted) NIfTI masks and compute all metrics.
    Uses EDT-based surface metrics (O(n_voxels)) instead of O(n_A × n_B) cdist.
    Checks _WORKER_CACHE for pre-computed mask-A data (one-to-many parallel mode).
    Saves outputs to a per-pair subdirectory.  Returns a result dict.
    """
    out_n     = args.out
    use_surf  = args.surface
    save_maps = not args.no_maps
    verbose   = getattr(args, 'verbose', False)

    # --- Check per-worker cache for mask A (populated by _init_worker) ---
    ca = _WORKER_CACHE.get('a', {})
    use_cache = (ca.get('path') == mask_a_path and
                 ca.get('class_val') == class_val and
                 'error' not in ca)

    img2 = nib.load(mask_b_path)
    aff2 = img2.affine

    if use_cache:
        aff1      = ca['aff'];  vox_size  = ca['vox']
        im1       = ca['im'];   mask1_bool = ca['mb']
        eroded1   = ca['er'];   outline1   = ca['ol']
        cog1      = ca['cog'];  cog1_xyz   = ca['cog_xyz']
        edt_full_a = ca['edt']
        if not np.allclose(aff1, aff2):
            raise ValueError(
                f'Affines do not match:\n  A: {mask_a_path}\n  B: {mask_b_path}\n'
                'Both masks must be in the same space with identical dimensions.'
            )
    else:
        img1 = nib.load(mask_a_path)
        aff1 = img1.affine
        if not np.allclose(aff1, aff2):
            raise ValueError(
                f'Affines do not match:\n  A: {mask_a_path}\n  B: {mask_b_path}\n'
                'Both masks must be in the same space with identical dimensions.'
            )
        vox_size = _voxel_size(aff1)
        raw1 = img1.get_fdata()
        if class_val is not None:
            im1 = np.uint16(np.round(raw1).astype(np.int64) == int(class_val))
            if not np.any(im1):
                avail = np.unique(np.round(raw1).astype(np.int64)[np.round(raw1) != 0]).tolist()
                raise ValueError(
                    f'No voxels for class_val={class_val} in {mask_a_path}.\n'
                    f'Non-zero values present: {avail}'
                )
        else:
            im1 = np.uint16(raw1)
        mask1_bool = im1.astype(bool)
        if not np.any(mask1_bool):
            raise ValueError(f'Mask A is empty (all zeros): {mask_a_path}')
        eroded1  = np.uint16(binary_erosion(mask1_bool))
        outline1 = np.uint16(np.abs(np.subtract(im1, eroded1)))
        cog1     = np.array(center_of_mass(im1))
        cog1_xyz = nib.affines.apply_affine(aff1, cog1)
        edt_full_a = None  # computed below jointly with edt_full_b

    raw2 = img2.get_fdata()
    im2  = np.uint16(raw2)
    mask2_bool = im2.astype(bool)
    if not np.any(mask2_bool):
        raise ValueError(f'Mask B is empty (all zeros): {mask_b_path}')
    eroded2  = np.uint16(binary_erosion(mask2_bool))
    outline2 = np.uint16(np.abs(np.subtract(im2, eroded2)))
    cog2     = np.array(center_of_mass(im2))
    cog2_xyz = nib.affines.apply_affine(aff2, cog2)
    cogs_d   = float(np.linalg.norm(cog1_xyz - cog2_xyz))

    # --- Overlap detection ---
    in_overlap  = np.uint16(np.multiply(im1, im2))
    has_overlap = bool(np.amax(in_overlap) != 0)
    if verbose:
        print(f'  Overlap: {"found" if has_overlap else "not found"}')

    # --- EDT: always computed (used for metrics AND optional maps) ---
    # If mask A is cached, only compute EDT-B; otherwise compute both in parallel.
    if edt_full_a is None:
        with ThreadPoolExecutor(max_workers=2) as tex:
            fut_a = tex.submit(_edt_full, mask1_bool, vox_size)
            fut_b = tex.submit(_edt_full, mask2_bool, vox_size)
            edt_full_a = fut_a.result()
            edt_full_b = fut_b.result()
    else:
        edt_full_b = _edt_full(mask2_bool, vox_size)

    # --- Surface metrics via EDT lookup: O(n_voxels) — no O(n_A × n_B) cdist ---
    sel1 = outline1.astype(bool) if use_surf else mask1_bool
    sel2 = outline2.astype(bool) if use_surf else mask2_bool
    d_a2b = edt_full_b[sel1]   # distance from each selected A voxel to nearest B
    d_b2a = edt_full_a[sel2]   # distance from each selected B voxel to nearest A

    min_dist         = float(min(d_a2b.min(), d_b2a.min()))
    hausdorff        = float(max(d_a2b.max(), d_b2a.max()))
    hausdorff_95     = float(max(np.percentile(d_a2b, 95), np.percentile(d_b2a, 95)))
    mean_surface_dist = float((d_a2b.mean() + d_b2a.mean()) / 2.0)

    # Closest voxel pair (argmin into EDT-selected index arrays)
    ijk1_sel = np.array(np.where(sel1)).T
    ijk2_sel = np.array(np.where(sel2)).T
    a_vox_vv = ijk1_sel[int(d_a2b.argmin())]
    b_vox_vv = ijk2_sel[int(d_b2a.argmin())]
    a_vox_mm = nib.affines.apply_affine(aff1, a_vox_vv)
    b_vox_mm = nib.affines.apply_affine(aff2, b_vox_vv)

    # --- COG distances ---
    # Scalar: O(1) EDT lookup at the rounded COG voxel
    def _ci(cog, shape):
        return tuple(int(np.clip(round(c), 0, s - 1)) for c, s in zip(cog, shape))
    c1i = _ci(cog1, im1.shape);  c2i = _ci(cog2, im2.shape)
    coga_2b = float(edt_full_b[c1i])
    cogb_2a = float(edt_full_a[c2i])

    # Nearest mask voxel to each COG: single-row O(n) cdist
    ijk2_all = np.array(np.where(mask2_bool)).T
    xyz2_all = nib.affines.apply_affine(aff2, ijk2_all)
    ijk1_all = np.array(np.where(mask1_bool)).T
    xyz1_all = nib.affines.apply_affine(aff1, ijk1_all)
    cog1_ds  = cdist(np.array([cog1_xyz]), xyz2_all)[0]
    cog2_ds  = cdist(np.array([cog2_xyz]), xyz1_all)[0]
    ca2b_vv  = ijk2_all[int(cog1_ds.argmin())];  ca2b_mm = xyz2_all[int(cog1_ds.argmin())]
    cb2a_vv  = ijk1_all[int(cog2_ds.argmin())];  cb2a_mm = xyz1_all[int(cog2_ds.argmin())]

    # --- Per-pair output subdirectory ---
    out_dir  = getattr(args, 'out_dir', os.path.abspath(out_n + '_output'))
    pair_stem = _pair_tag(mask_a_path, mask_b_path, class_label).lstrip('_')
    pair_dir  = os.path.join(out_dir, pair_stem)
    os.makedirs(pair_dir, exist_ok=True)

    def _p(name): return os.path.join(pair_dir, name)

    # --- Optional distance map NIfTIs (EDT already computed above) ---
    dist_a2b = None
    if save_maps:
        dist_a2b = (edt_full_b * mask1_bool).astype(np.float32)
        dist_b2a = (edt_full_a * mask2_bool).astype(np.float32)
        nib.save(nib.Nifti1Image(dist_a2b,   aff1), _p('dist_map_A_to_B.nii.gz'))
        nib.save(nib.Nifti1Image(dist_b2a,   aff2), _p('dist_map_B_to_A.nii.gz'))
        nib.save(nib.Nifti1Image(edt_full_b, aff2), _p('dist_map_full_to_B.nii.gz'))
        nib.save(nib.Nifti1Image(edt_full_a, aff1), _p('dist_map_full_to_A.nii.gz'))

    # --- Snapshot ---
    # Snapshot dispatch:
    #   --bg-image --volume-render → VTK glass-brain 3-D (no display needed)
    #   --bg-image                 → mrview orthographic 3-panel (needs xvfb-run)
    #   (neither)                  → matplotlib
    snap_file    = _p('snapshot.png')
    snap_b64     = ''
    bg_image     = getattr(args, 'bg_image', None)
    vol_render   = getattr(args, 'volume_render', False)

    # Save thick (3-voxel) edge B early so mrview can use it as an outline overlay
    b_edge_path  = _p('mask_B_edge.nii.gz')
    outline2_thick = np.uint16(im2 - binary_erosion(mask2_bool, iterations=3).astype(np.uint16))
    nib.save(nib.Nifti1Image(outline2_thick, aff2), b_edge_path)

    if bg_image:
        focus_vox = tuple(int(round((cog1[i] + cog2[i]) / 2)) for i in range(3))
        dist_path = _p('dist_map_A_to_B.nii.gz') if save_maps else None
        if vol_render:
            snap_b64 = _make_vtk_snapshot(
                bg_image, mask_a_path, mask_b_path, snap_file, focus_vox,
            )
            if not snap_b64:
                print('  [vtk snapshot failed — falling back to matplotlib]', flush=True)
        else:
            snap_b64 = _make_mrview_snapshot(
                bg_image, mask_a_path, mask_b_path,
                dist_path, snap_file, focus_vox, b_edge_path,
            )
            if not snap_b64:
                print('  [mrview snapshot failed — falling back to matplotlib]', flush=True)

    if not snap_b64:
        snap_b64 = _make_snapshot(
            im1, im2, dist_a2b, a_vox_vv, b_vox_vv, cog1, cog2,
            mask_a_path, mask_b_path, class_label,
        )
        if snap_b64:
            with open(snap_file, 'wb') as fh:
                fh.write(base64.b64decode(snap_b64))

    # --- Mask / edge NIfTIs ---
    nib.save(nib.Nifti1Image(im1,      aff1), _p('mask_A.nii.gz'))
    nib.save(nib.Nifti1Image(im2,      aff2), _p('mask_B.nii.gz'))
    nib.save(nib.Nifti1Image(outline1, aff1), _p('mask_A_edge.nii.gz'))
    nib.save(nib.Nifti1Image(outline2, aff2), _p('mask_B_edge.nii.gz'))
    nib.save(nib.Nifti1Image(eroded1,  aff1), _p('mask_A_eroded.nii.gz'))
    nib.save(nib.Nifti1Image(eroded2,  aff2), _p('mask_B_eroded.nii.gz'))

    # --- Point-marker NIfTIs ---
    nib.save(nib.Nifti1Image(_dilated_marker(im1.shape, a_vox_vv), aff1),
             _p('mask_A_vox_mindist.nii.gz'))
    nib.save(nib.Nifti1Image(_dilated_marker(im2.shape, b_vox_vv), aff2),
             _p('mask_B_vox_mindist.nii.gz'))
    nib.save(nib.Nifti1Image(_dilated_marker(im1.shape, np.round(cog1).astype(int)), aff1),
             _p('mask_A_COG.nii.gz'))
    nib.save(nib.Nifti1Image(_dilated_marker(im2.shape, np.round(cog2).astype(int)), aff2),
             _p('mask_B_COG.nii.gz'))

    # --- Overlap workflow ---
    ov_results: dict = dict(has_overlap=False, ov_count=None,
                            ov_perc_mA=None, ov_perc_mB=None,
                            min_AvsOVCOG=None, min_OV2ACOG=None, min_OV2BCOG=None)

    if has_overlap:
        ov_count  = int(np.count_nonzero(in_overlap))
        ov_idx    = np.where(in_overlap)
        ov_ijk    = np.vstack(ov_idx).T
        ov_xyz    = nib.affines.apply_affine(aff2, ov_ijk)
        ov_cog    = np.array(center_of_mass(in_overlap))
        ov_cog_xyz = nib.affines.apply_affine(aff2, ov_cog)

        ov_cog_2_Av = cdist(np.array([ov_cog_xyz]), xyz1_all)[0]
        OVv_2_ACOG  = cdist(ov_xyz, np.array([cog1_xyz]))[:, 0]
        OVv_2_BCOG  = cdist(ov_xyz, np.array([cog2_xyz]))[:, 0]

        n_A = int(np.count_nonzero(im1))
        n_B = int(np.count_nonzero(im2))

        ca_v = ijk1_all[int(ov_cog_2_Av.argmin())]
        cb_v = ov_ijk[int(OVv_2_ACOG.argmin())]
        cc_v = ov_ijk[int(OVv_2_BCOG.argmin())]

        ov_Vox_map = np.zeros(in_overlap.shape, np.uint16)
        ov_Vox_map[ov_ijk[:, 0], ov_ijk[:, 1], ov_ijk[:, 2]] = 1

        ov_cog_map = np.zeros(in_overlap.shape, np.uint16)
        ov_cog_map[int(ov_cog[0]), int(ov_cog[1]), int(ov_cog[2])] = 1
        dil_cogOV  = np.uint16(binary_dilation(ov_cog_map.astype(bool), iterations=5))
        dil_cogOV[int(ov_cog[0]), int(ov_cog[1]), int(ov_cog[2])] = 10

        nib.save(nib.Nifti1Image(ov_Vox_map, aff2), _p('overlap_voxels.nii.gz'))
        nib.save(nib.Nifti1Image(dil_cogOV,  aff2), _p('overlap_COG.nii.gz'))
        nib.save(nib.Nifti1Image(_dilated_marker(im1.shape, ca_v), aff1),
                 _p('maskA_mindist_2_overlap_COG.nii.gz'))
        nib.save(nib.Nifti1Image(_dilated_marker(in_overlap.shape, cb_v), aff2),
                 _p('overlap_mindist_2_maskA_COG.nii.gz'))
        nib.save(nib.Nifti1Image(_dilated_marker(in_overlap.shape, cc_v), aff2),
                 _p('overlap_mindist_2_maskB_COG.nii.gz'))

        ov_results = dict(
            has_overlap=True, ov_count=ov_count,
            ov_perc_mA=100.0 * ov_count / n_A if n_A else 0.0,
            ov_perc_mB=100.0 * ov_count / n_B if n_B else 0.0,
            min_AvsOVCOG=float(ov_cog_2_Av.min()),
            min_OV2ACOG=float(OVv_2_ACOG.min()),
            min_OV2BCOG=float(OVv_2_BCOG.min()),
        )

    # --- Result dict ---
    result = {
        'mask_a':            mask_a_path,
        'mask_b':            mask_b_path,
        'class_val':         class_val,
        'class_label':       class_label,
        'pair_dir':          pair_dir,
        'cog1_vox':          cog1.tolist(),
        'cog1_xyz':          cog1_xyz.tolist(),
        'cog2_vox':          cog2.tolist(),
        'cog2_xyz':          cog2_xyz.tolist(),
        'cogs_d':            cogs_d,
        'min_dist':          min_dist,
        'hausdorff':         hausdorff,
        'hausdorff_95':      hausdorff_95,
        'mean_surface_dist': mean_surface_dist,
        'coga_2b':           coga_2b,
        'cogb_2a':           cogb_2a,
        'a_vox_vv':          a_vox_vv.tolist(),
        'a_vox_mm':          a_vox_mm.tolist(),
        'b_vox_vv':          b_vox_vv.tolist(),
        'b_vox_mm':          b_vox_mm.tolist(),
        'ca2b_vv':           ca2b_vv.tolist(),
        'ca2b_mm':           ca2b_mm.tolist(),
        'cb2a_vv':           cb2a_vv.tolist(),
        'cb2a_mm':           cb2a_mm.tolist(),
        'snapshot_b64':      snap_b64,
        **ov_results,
    }

    _write_text(pair_dir, result)
    return result


# ---------------------------------------------------------------------------
# Text / CSV / HTML output
# ---------------------------------------------------------------------------

def _fmt(vals) -> str:
    return '[' + ', '.join(f'{v:.2f}' for v in vals) + ']'


def _write_text(pair_dir: str, r: dict):
    with open(os.path.join(pair_dir, 'output_measures.txt'), 'w') as fh:
        fh.write('KUL_EDs_b2masks — distance results\n')
        fh.write('=' * 60 + '\n')
        fh.write(f'Mask A : {r["mask_a"]}')
        if r.get('class_label'):
            fh.write(f'  [class {r["class_val"]} — {r["class_label"]}]')
        fh.write(f'\nMask B : {r["mask_b"]}\n\n')

        if r.get('has_overlap'):
            fh.write('OVERLAP DETECTED\n')
            fh.write(f'  Overlapping voxels            : {r["ov_count"]}\n')
            fh.write(f'  Overlap % relative to mask A  : {r["ov_perc_mA"]:.2f}%\n')
            fh.write(f'  Overlap % relative to mask B  : {r["ov_perc_mB"]:.2f}%\n')
            fh.write(f'  Min dist overlap COG → A vox  : {r["min_AvsOVCOG"]:.4f} mm\n')
            fh.write(f'  Min dist A COG → overlap vox  : {r["min_OV2ACOG"]:.4f} mm\n')
            fh.write(f'  Min dist B COG → overlap vox  : {r["min_OV2BCOG"]:.4f} mm\n\n')
        else:
            fh.write('No overlap between masks.\n\n')

        fh.write('DISTANCE METRICS\n')
        fh.write(f'  Min voxel-to-voxel             : {r["min_dist"]:.4f} mm\n')
        fh.write(f'    Mask A voxel (vox)           : {r["a_vox_vv"]}\n')
        fh.write(f'    Mask A voxel (mm)            : {_fmt(r["a_vox_mm"])}\n')
        fh.write(f'    Mask B voxel (vox)           : {r["b_vox_vv"]}\n')
        fh.write(f'    Mask B voxel (mm)            : {_fmt(r["b_vox_mm"])}\n\n')
        fh.write(f'  Hausdorff distance             : {r["hausdorff"]:.4f} mm\n')
        fh.write(f'  95th-pct Hausdorff             : {r["hausdorff_95"]:.4f} mm\n')
        fh.write(f'  Mean surface distance (ASSD)   : {r["mean_surface_dist"]:.4f} mm\n\n')
        fh.write(f'  COG A → COG B distance         : {r["cogs_d"]:.4f} mm\n')
        fh.write(f'    COG A (vox)                  : {_fmt(r["cog1_vox"])}\n')
        fh.write(f'    COG A (mm)                   : {_fmt(r["cog1_xyz"])}\n')
        fh.write(f'    COG B (vox)                  : {_fmt(r["cog2_vox"])}\n')
        fh.write(f'    COG B (mm)                   : {_fmt(r["cog2_xyz"])}\n\n')
        fh.write(f'  Min dist COG A → B surface     : {r["coga_2b"]:.4f} mm\n')
        fh.write(f'    Nearest B voxel (vox)        : {r["ca2b_vv"]}\n')
        fh.write(f'    Nearest B voxel (mm)         : {_fmt(r["ca2b_mm"])}\n\n')
        fh.write(f'  Min dist COG B → A surface     : {r["cogb_2a"]:.4f} mm\n')
        fh.write(f'    Nearest A voxel (vox)        : {r["cb2a_vv"]}\n')
        fh.write(f'    Nearest A voxel (mm)         : {_fmt(r["cb2a_mm"])}\n')


_CSV_FIELDS_BASE = [
    'mask_a', 'mask_b',
    'min_dist', 'hausdorff', 'hausdorff_95', 'mean_surface_dist',
    'cogs_d', 'coga_2b', 'cogb_2a',
    'has_overlap', 'ov_count', 'ov_perc_mA', 'ov_perc_mB',
]
_CSV_FIELDS_CLASS = ['class_val', 'class_label'] + _CSV_FIELDS_BASE


def _write_csv(out_dir: str, out_n: str, results: list):
    has_cls = any(r.get('class_label') for r in results)
    fields  = _CSV_FIELDS_CLASS if has_cls else _CSV_FIELDS_BASE
    csv_path = os.path.join(out_dir, out_n + '_metrics.csv')
    with open(csv_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f'CSV    : {csv_path}')


_HTML = '''\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>KUL_EDs — {out_n}</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{background:#111;color:#ddd;font-family:'Courier New',monospace;font-size:12px;padding:20px}}
h1{{color:#8cf;font-size:17px;margin-bottom:5px}}
.meta{{color:#777;margin-bottom:16px;font-size:11px}}
.meta b{{color:#aaa}}
table{{border-collapse:collapse;width:100%}}
thead th{{background:#1b2232;color:#8cf;padding:8px 11px;cursor:pointer;
          position:sticky;top:0;z-index:1;text-align:left;
          border-bottom:2px solid #2a3550;white-space:nowrap;user-select:none}}
thead th:hover{{background:#252f42}}
thead th.asc::after{{content:" ▲";font-size:9px;color:#aef}}
thead th.desc::after{{content:" ▼";font-size:9px;color:#aef}}
tbody tr:nth-child(even){{background:#161616}}
tbody tr:hover{{background:#1a2336}}
td{{padding:5px 11px;border-bottom:1px solid #222;vertical-align:middle}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
td.tract{{color:#8f8;max-width:200px;word-break:break-word}}
td.cls{{color:#fa8;font-weight:bold}}
img.thumb{{width:190px;cursor:zoom-in;border-radius:3px;display:block;transition:transform .1s}}
img.thumb:hover{{transform:scale(1.04)}}
#modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.93);
        z-index:9999;align-items:center;justify-content:center;cursor:zoom-out}}
#modal.on{{display:flex}}
#modal img{{max-width:96vw;max-height:96vh;border-radius:4px}}
</style>
</head>
<body>
<h1>KUL_EDs — Distance Report</h1>
<p class="meta">Output: <b>{out_n}</b> &nbsp;|&nbsp; Pairs: <b>{n}</b> &nbsp;|&nbsp; Date: {date}</p>
<table id="T">
<thead><tr>{th}</tr></thead>
<tbody>{rows}</tbody>
</table>
<div id="modal" onclick="this.classList.remove('on')">
  <img id="mi" src="" alt="">
</div>
<script>
var D={{}};
function sort(c){{
  var T=document.getElementById('T'),ths=T.querySelectorAll('thead th');
  ths.forEach(function(t){{t.classList.remove('asc','desc')}});
  D[c]=!D[c];ths[c].classList.add(D[c]?'asc':'desc');
  var tb=T.tBodies[0],rows=Array.from(tb.rows);
  rows.sort(function(a,b){{
    var av=a.cells[c].getAttribute('data-v')||a.cells[c].textContent.trim();
    var bv=b.cells[c].getAttribute('data-v')||b.cells[c].textContent.trim();
    var an=parseFloat(av),bn=parseFloat(bv);
    if(!isNaN(an)&&!isNaN(bn)) return D[c]?an-bn:bn-an;
    return D[c]?av.localeCompare(bv):bv.localeCompare(av);
  }});
  rows.forEach(function(r){{tb.appendChild(r)}});
}}
function zoom(img){{
  document.getElementById('mi').src=img.src;
  document.getElementById('modal').classList.add('on');
  event.stopPropagation();
}}
window.onload=function(){{sort({sort_col})}};
</script>
</body>
</html>'''


def _write_html_report(out_dir: str, out_n: str, results: list):
    has_cls = any(r.get('class_label') for r in results)

    # Column definitions: (header, data-key or callable, css-class, sortable?)
    cols = []
    if has_cls:
        cols.append(('Class', lambda r: r.get('class_label', ''), 'cls', True))
    cols += [
        ('Mask A',         lambda r: Path(r['mask_a']).stem, '',      True),
        ('Mask B',         lambda r: Path(r['mask_b']).stem, 'tract', True),
        ('Min (mm)',        'min_dist',          'num', True),
        ('Hausdorff (mm)', 'hausdorff',          'num', True),
        ('H95 (mm)',       'hausdorff_95',       'num', True),
        ('ASSD (mm)',      'mean_surface_dist',  'num', True),
        ('COG↔COG (mm)',  'cogs_d',             'num', True),
        ('Overlap',        None,                 '',    False),
        ('Snapshot',       None,                 '',    False),
    ]
    assd_col = next(i for i, c in enumerate(cols) if c[0] == 'ASSD (mm)')

    th = ''.join(
        f'<th onclick="sort({i})">{lbl}</th>' if sortable else f'<th>{lbl}</th>'
        for i, (lbl, _, __, sortable) in enumerate(cols)
    )

    rows = []
    for r in results:
        cells = []
        for lbl, key, cls, _ in cols:
            if lbl == 'Overlap':
                if r.get('has_overlap'):
                    v = (f'Yes<br><small>{r["ov_count"]} vox · '
                         f'{r["ov_perc_mA"]:.1f}%/A · {r["ov_perc_mB"]:.1f}%/B</small>')
                else:
                    v = 'No'
                cells.append(f'<td>{v}</td>')
            elif lbl == 'Snapshot':
                b64 = r.get('snapshot_b64', '')
                cells.append(
                    f'<td><img class="thumb" src="data:image/png;base64,{b64}" '
                    f'onclick="zoom(this)"></td>' if b64 else '<td>—</td>'
                )
            elif callable(key):
                cells.append(f'<td class="{cls}">{key(r)}</td>')
            else:
                v = r.get(key, '—')
                if isinstance(v, float):
                    cells.append(f'<td class="{cls}" data-v="{v:.6f}">{v:.3f}</td>')
                else:
                    cells.append(f'<td class="{cls}">{v}</td>')
        rows.append('<tr>' + ''.join(cells) + '</tr>')

    html = _HTML.format(
        out_n=out_n, n=len(results),
        date=datetime.date.today().isoformat(),
        th=th, rows='\n'.join(rows),
        sort_col=assd_col,
    )
    path = os.path.join(out_dir, 'summary_report.html')
    with open(path, 'w') as fh:
        fh.write(html)
    print(f'Report : {path}')


def _print_pair_result(r: dict):
    if r.get('class_label'):
        print(f'  Class                         : {r["class_label"]} (val={r["class_val"]})')
    print(f'  COG A (mm)                    : {_fmt(r["cog1_xyz"])}')
    print(f'  COG B (mm)                    : {_fmt(r["cog2_xyz"])}')
    print(f'  Min voxel-to-voxel distance   : {r["min_dist"]:.4f} mm')
    print(f'  Hausdorff distance            : {r["hausdorff"]:.4f} mm')
    print(f'  95th-pct Hausdorff            : {r["hausdorff_95"]:.4f} mm')
    print(f'  Mean surface distance (ASSD)  : {r["mean_surface_dist"]:.4f} mm')
    print(f'  COG-to-COG distance           : {r["cogs_d"]:.4f} mm')
    if r.get('has_overlap'):
        print(f'  Overlap voxels                : {r["ov_count"]} '
              f'({r["ov_perc_mA"]:.1f}% of A, {r["ov_perc_mB"]:.1f}% of B)')


def _print_summary_table(results: list):
    w = 35
    has_cls = any(r.get('class_label') for r in results)
    hdr = (
        (f'{"Class":<12}  ' if has_cls else '') +
        f'{"Mask A":<{w}}  {"Mask B":<{w}}'
        f'  {"Min":>8}  {"Hdorff":>8}  {"H95":>8}  {"ASSD":>8}'
    )
    sep = '=' * len(hdr)
    print(f'\n{sep}\nSUMMARY\n{sep}\n{hdr}\n{"-"*len(hdr)}')
    for r in sorted(results, key=lambda x: x['mean_surface_dist']):
        cls = f'{r.get("class_label",""):<12}  ' if has_cls else ''
        print(
            f'{cls}{Path(r["mask_a"]).name[:w]:<{w}}  {Path(r["mask_b"]).name[:w]:<{w}}'
            f'  {r["min_dist"]:8.3f}  {r["hausdorff"]:8.3f}'
            f'  {r["hausdorff_95"]:8.3f}  {r["mean_surface_dist"]:8.3f}'
        )
    print(sep)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='KUL_EDs_b2masks.py',
        description=(
            'Swiss-knife distance tool for binary NIfTI masks.\n\n'
            'Metrics : min Euclidean · Hausdorff · 95th-pct Hausdorff · ASSD\n'
            'Modes   : one-to-one · one-to-many · all-pairs\n\n'
            'Examples:\n'
            '  One-to-one  : KUL_EDs_b2masks.py -a A.nii.gz -b B.nii.gz -o out\n'
            '  One-to-many : KUL_EDs_b2masks.py -a A.nii.gz -b B1.nii.gz B2.nii.gz -o out -n 4\n'
            '  All-pairs   : KUL_EDs_b2masks.py -a A.nii.gz B.nii.gz C.nii.gz -o out\n'
            '  From list   : KUL_EDs_b2masks.py -l masks.txt -o out -n -1 --csv\n'
            '  Multi-class : KUL_EDs_b2masks.py -a seg.nii.gz -b tract.nii.gz -o out\n'
            '                  --classes 1 2 --class-labels core edema'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('-a', '--mask-a', nargs='+', metavar='MASK',
                     help='Source mask(s). With -b: one-to-one or one-to-many. '
                          'Without -b (≥2 masks): all-pairs mode.')
    src.add_argument('-l', '--list', metavar='FILE',
                     help='Text file (one mask path per line) for all-pairs mode.')
    p.add_argument('-b', '--mask-b', nargs='+', metavar='MASK',
                   help='Target mask(s). Omit for all-pairs mode.')
    p.add_argument('-o', '--out', default='KUL_EDs',
                   help='Output prefix / path (default: KUL_EDs).')
    p.add_argument('--classes', nargs='+', type=int, metavar='VAL',
                   help='Integer class values to extract from mask A separately '
                        '(e.g. --classes 1 2).')
    p.add_argument('--class-labels', nargs='+', metavar='LABEL',
                   help='Human-readable labels for --classes (e.g. core edema). '
                        'Must match the number of --classes values.')
    p.add_argument('--metric', choices=['euclidean', 'hausdorff', 'mean_surface', 'all'],
                   default='all',
                   help='Distance metric to report (default: all; all are always computed).')
    p.add_argument('--surface', action='store_true',
                   help='Restrict Hausdorff/ASSD to 1-voxel erosion edge voxels '
                        '(default: all mask voxels used via EDT lookup).')
    p.add_argument('--no-maps', action='store_true',
                   help='Skip float distance map NIfTIs (faster; snapshot still generated).')
    p.add_argument('--bg-image', metavar='NII',
                   help='Anatomical reference image (T1w, FA, …) used as snapshot background. '
                        'Without --volume-render: mrview orthographic 3-panel (needs xvfb-run). '
                        'With --volume-render: VTK off-screen glass-brain 3-D render.')
    p.add_argument('--volume-render', action='store_true',
                   help='Use VTK off-screen glass-brain 3-D volume rendering for snapshots '
                        '(superior | lateral | posterior views). Requires --bg-image. '
                        'No display or xvfb-run needed.')
    p.add_argument('--csv', action='store_true', help='Write _metrics.csv.')
    p.add_argument('-n', '--workers', type=int, default=1,
                   help='Parallel workers for multi-pair mode (-1 = all CPUs, default: 1).')
    p.add_argument('-v', '--verbose', action='store_true')
    return p


def main():
    parser = _build_parser()
    args   = parser.parse_args()

    # Normalize -o to a plain stem + resolve output dir
    _op = Path(args.out.rstrip('/\\'))
    args.out     = _op.name
    args.out_dir = os.path.abspath(str(_op.parent / (_op.name + '_output')))

    # Validate
    if args.list and args.mask_b:
        parser.error('-b / --mask-b cannot be combined with -l / --list.')
    if args.class_labels and args.classes and len(args.class_labels) != len(args.classes):
        parser.error('--class-labels must have the same number of entries as --classes.')
    if args.volume_render and not args.bg_image:
        parser.error('--volume-render requires --bg-image.')
    if args.bg_image and not os.path.isfile(args.bg_image):
        parser.error(f'--bg-image not found: {args.bg_image}')

    # --- Build base mask-pair list ---
    if args.list:
        masks = [ln.strip() for ln in open(args.list) if ln.strip()]
        if len(masks) < 2:
            parser.error(f'List file needs ≥2 paths; found {len(masks)}.')
        base_pairs = list(itertools.combinations(masks, 2))
        mode = 'all-pairs'
    elif args.mask_a and not args.mask_b:
        if len(args.mask_a) < 2:
            parser.error('All-pairs mode requires ≥2 masks with -a (or use -l).')
        base_pairs = list(itertools.combinations(args.mask_a, 2))
        mode = 'all-pairs'
    else:
        if len(args.mask_a) > 1:
            parser.error('One-to-many needs exactly one -a mask and one or more -b masks.')
        if len(args.mask_b) == 1:
            base_pairs = [(args.mask_a[0], args.mask_b[0])]
            mode = 'one-to-one'
        else:
            base_pairs = [(args.mask_a[0], b) for b in args.mask_b]
            mode = 'one-to-many'

    # --- Expand by classes ---
    if args.classes:
        labels = args.class_labels or [str(v) for v in args.classes]
        full_pairs = [(a, b, cv, cl)
                      for a, b in base_pairs
                      for cv, cl in zip(args.classes, labels)]
    else:
        full_pairs = [(a, b, None, '') for a, b in base_pairs]

    print(f'KUL_EDs_b2masks | mode: {mode} | pairs: {len(full_pairs)} '
          f'| metric: {args.metric} | workers: {args.workers}')

    os.makedirs(args.out_dir, exist_ok=True)
    all_results = []

    if len(full_pairs) == 1 or args.workers == 1:
        for a, b, cv, cl in full_pairs:
            print(f'\nProcessing: {Path(a).name}  ↔  {Path(b).name}'
                  + (f'  [{cl}]' if cl else ''))
            r = compute_pair(a, b, args, cv, cl)
            all_results.append(r)
            _print_pair_result(r)
    else:
        n_workers = args.workers if args.workers > 0 else None

        # In one-to-many mode all pairs share the same mask A and class_val.
        # Pre-compute mask-A data once per worker via initializer (saves EDT-A per pair).
        unique_a = {(a, cv) for a, b, cv, cl in full_pairs}
        if mode == 'one-to-many' and len(unique_a) == 1:
            (init_a_path, init_cv), = unique_a
            pool_kwargs = dict(
                max_workers=n_workers,
                initializer=_init_worker,
                initargs=(init_a_path, init_cv, args.surface),
            )
        else:
            pool_kwargs = dict(max_workers=n_workers)

        with ProcessPoolExecutor(**pool_kwargs) as pool:
            futures = {pool.submit(compute_pair, a, b, args, cv, cl): (a, b, cl)
                       for a, b, cv, cl in full_pairs}
            for fut in as_completed(futures):
                a, b, cl = futures[fut]
                try:
                    r = fut.result()
                    all_results.append(r)
                    print(f'Done: {Path(a).name} ↔ {Path(b).name}'
                          + (f' [{cl}]' if cl else '')
                          + f'  ASSD={r["mean_surface_dist"]:.2f} mm'
                          + f'  H={r["hausdorff"]:.2f} mm')
                except Exception as exc:
                    print(f'ERROR: {Path(a).name} ↔ {Path(b).name}: {exc}',
                          file=sys.stderr)

    if len(all_results) > 1:
        _print_summary_table(all_results)

    _write_html_report(args.out_dir, args.out, all_results)

    if args.csv:
        _write_csv(args.out_dir, args.out, all_results)


if __name__ == '__main__':
    main()
