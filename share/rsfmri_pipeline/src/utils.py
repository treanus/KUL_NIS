"""
Shared utilities: file discovery, atlas loading, BOLD I/O, FC computation,
mask building. All other pipeline steps import from here.
"""

import json
import logging
import os
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import nibabel as nib
import numpy as np
from nilearn import datasets, image
from scipy.ndimage import label as ndimage_label
from scipy.stats import norm

from config import (
    ATLAS_DIR, DENOISED_DIR, FMRIPREP_DIR, FSL_ATLAS_DIR, MELODIC_BIN,
    MNI_SPACE, MIN_MASK_VOXELS,
    SYNTHSEG_DIR, ATLASES,
)

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# File path helpers
# ---------------------------------------------------------------------------

def bold_path(subject: str, run_id: str) -> Path:
    func_dir = DENOISED_DIR / f"sub-{subject}" / "func"
    # Allow optional extra BIDS entities (e.g. acq-singleTE) between task and run
    dirs = sorted(func_dir.glob(f"*_{run_id}_*_postproc_nilearn"))
    if dirs:
        hits = sorted(dirs[0].glob(f"*_{run_id}_*_desc-denoised_bold.nii.gz"))
        if hits:
            return hits[0]
    tag = f"sub-{subject}_task-rest_{run_id}_space-{MNI_SPACE}_res-2"
    return func_dir / f"{tag}_postproc_nilearn" / f"{tag}_desc-denoised_bold.nii.gz"


def mask_path(subject: str, run_id: str) -> Path:
    # Prefer an explicit-resolution (res-2) fmriprep output if present; fall
    # back to whatever resolution tag (or none) fmriprep actually produced --
    # not every fmriprep run is invoked with --output-spaces ...:res-2.
    func_dir = FMRIPREP_DIR / f"sub-{subject}" / "func"
    hits = sorted(func_dir.glob(f"*_{run_id}_*space-{MNI_SPACE}_res-2_desc-brain_mask.nii.gz"))
    if hits:
        return hits[0]
    hits = sorted(func_dir.glob(f"*_{run_id}_*space-{MNI_SPACE}*_desc-brain_mask.nii.gz"))
    if hits:
        return hits[0]
    tag = f"sub-{subject}_{run_id}_space-{MNI_SPACE}_res-2"
    return func_dir / f"{tag}_desc-brain_mask.nii.gz"


def is_rest_run(run_id: str) -> bool:
    """True iff run_id's task label is exactly 'rest' (vs a tbfMRI task)."""
    return bool(re.match(r"task-rest(_run-|$)", run_id))


def get_run_tr(subject: str, run_id: str) -> float:
    """
    Real RepetitionTime (seconds) for one run, read from the fmriprep
    preproc BOLD JSON sidecar -- NOT from the denoised NIfTI header,
    whose pixdim[4] is unreliable (nilearn/SUSAN postprocessing resets it
    to 1.0 rather than preserving the true TR). Different tasks can have
    genuinely different TRs (e.g. a tbfMRI task protocol vs resting-state),
    so this must be looked up per-run rather than assumed constant.
    """
    func_dir = FMRIPREP_DIR / f"sub-{subject}" / "func"
    hits = sorted(func_dir.glob(f"*_{run_id}_*space-{MNI_SPACE}*_desc-preproc_bold.json"))
    if not hits:
        hits = sorted(func_dir.glob(f"*_{run_id}_*_desc-preproc_bold.json"))
    if not hits:
        raise FileNotFoundError(f"No preproc BOLD sidecar found for sub-{subject}/{run_id}")
    with open(hits[0]) as fh:
        sidecar = json.load(fh)
    return float(sidecar["RepetitionTime"])


def t1w_path(subject: str) -> Path:
    # Same res-2-preferred, resolution-agnostic-fallback policy as mask_path().
    anat_dir = FMRIPREP_DIR / f"sub-{subject}" / "anat"
    res2 = anat_dir / f"sub-{subject}_space-{MNI_SPACE}_res-2_desc-preproc_T1w.nii.gz"
    if res2.exists():
        return res2
    hits = sorted(anat_dir.glob(f"sub-{subject}_space-{MNI_SPACE}*_desc-preproc_T1w.nii.gz"))
    return hits[0] if hits else res2


def gm_prob_path(subject: str) -> Path:
    # Same res-2-preferred, resolution-agnostic-fallback policy as mask_path().
    anat_dir = FMRIPREP_DIR / f"sub-{subject}" / "anat"
    res2 = anat_dir / f"sub-{subject}_space-{MNI_SPACE}_res-2_label-GM_probseg.nii.gz"
    if res2.exists():
        return res2
    hits = sorted(anat_dir.glob(f"sub-{subject}_space-{MNI_SPACE}*_label-GM_probseg.nii.gz"))
    return hits[0] if hits else res2


def synthseg_path(subject: str) -> Path:
    return SYNTHSEG_DIR / f"sub-{subject}" / f"sub-{subject}_synthseg_MNI.nii.gz"


def language_mask_path(subject: str) -> Path:
    return SYNTHSEG_DIR / f"sub-{subject}" / f"sub-{subject}_Language_MNI.nii.gz"


def auditory_mask_path(subject: str) -> Path:
    return SYNTHSEG_DIR / f"sub-{subject}" / f"sub-{subject}_Auditory_MNI.nii.gz"


# ---------------------------------------------------------------------------
# Run discovery
# ---------------------------------------------------------------------------

_RUN_RE  = re.compile(r"run-\d+")
_TASK_RE = re.compile(r"task-[A-Za-z0-9]+")


def get_runs(subject: str) -> list[str]:
    """
    Return sorted, task-scoped run identifiers present for a subject in the
    denoised directory, e.g. "task-rest_run-01", "task-rest_run-02",
    "task-HAND" (BIDS omits run- when a task has exactly one run).

    Identifiers are always task-scoped ("task-<label>[_run-N]"), never a
    bare "run-N" -- BIDS run-numbering is conventionally per-task, so two
    different tasks commonly both have a "run-01". A bare "run-N" identifier
    would collide across tasks in that case: get_runs() would silently
    deduplicate them to one entry, and bold_path()/mask_path()'s substring
    glob would then match both tasks' directories, arbitrarily picking one
    and dropping the other with no error.

    bold_path()/mask_path() key off this identifier purely as a literal
    glob substring, so the task- prefix resolves the same way -- no changes
    needed there.
    """
    func_dir = DENOISED_DIR / f"sub-{subject}" / "func"
    runs: set[str] = set()
    for p in func_dir.iterdir():
        if not p.is_dir():
            continue
        m_task = _TASK_RE.search(p.name)
        if not m_task:
            continue  # no task- entity -- not a valid BIDS func acquisition, skip
        m_run = _RUN_RE.search(p.name)
        runs.add(f"{m_task.group(0)}_{m_run.group(0)}" if m_run else m_task.group(0))
    if not runs:
        raise FileNotFoundError(f"No runs found for sub-{subject} in {func_dir}")
    return sorted(runs)


# ---------------------------------------------------------------------------
# BOLD / mask loading
# ---------------------------------------------------------------------------

def load_bold(subject: str, run_id: str) -> nib.Nifti1Image:
    p = bold_path(subject, run_id)
    if not p.exists():
        raise FileNotFoundError(f"Denoised BOLD not found: {p}")
    return nib.load(str(p))


def load_brain_mask(subject: str, run_id: str) -> nib.Nifti1Image:
    p = mask_path(subject, run_id)
    if not p.exists():
        raise FileNotFoundError(f"Brain mask not found: {p}")
    return nib.load(str(p))


# ---------------------------------------------------------------------------
# Shared MELODIC helpers (masked-ICA and whole-brain-ICA both use these)
# ---------------------------------------------------------------------------

def _zscore_run(arr: np.ndarray) -> np.ndarray:
    """
    Per-voxel z-score along the time axis, independently for this one run.

    MELODIC only demeans the concatenated series once, globally -- it has
    no notion of "run boundary". Separately-acquired runs routinely carry a
    per-voxel baseline offset from each other (registration drift, gain,
    whatever) that is huge relative to genuine BOLD variance: confirmed on
    real data (two same-TR task runs from the same session) at std=2230
    for the per-voxel baseline difference between runs vs std=50 for a
    voxel's own real within-run temporal variability -- 44x. Left
    unaddressed, every component of a joint ICA becomes a "which run is
    this" step function instead of any real spatial network (confirmed:
    all 5 components' timecourses were a clean +-1 step at the run
    boundary with near-zero within-run noise). Per-run normalisation
    before concatenation is standard practice for multi-run/session
    temporal-concatenation ICA specifically to prevent this.
    """
    mean = arr.mean(axis=-1, keepdims=True)
    std  = arr.std(axis=-1, keepdims=True)
    std  = np.where(std == 0, 1.0, std)
    return (arr - mean) / std


def concatenate_bold_runs(subject: str, runs: list[str]) -> nib.Nifti1Image | None:
    """
    Load all BOLD runs for a subject, per-voxel z-score each independently
    (see _zscore_run), and concatenate along the time axis.
    Runs are usually already on the same grid, but tasks acquired at
    different native resolutions (e.g. a resting-state run with different
    slice thickness than task runs, when fmriprep wasn't pinned to a
    common res-2 grid -- see average_maps() for the same issue) can differ
    in shape even within the same named MNI space, so each run is
    resampled onto the first successfully-loaded run's grid before
    concatenating.
    Returns a single 4D NIfTI, or None if no runs could be loaded.
    """
    imgs = []

    for run_id in runs:
        try:
            img = load_bold(subject, run_id)
            imgs.append(img)
            log.info("    Loaded %s (%d vols)", run_id, img.shape[-1])
        except Exception as exc:
            log.warning("    Could not load %s for sub-%s: %s", run_id, subject, exc)

    if not imgs:
        return None

    ref = imgs[0]
    arrays = []
    for img in imgs:
        if img.shape[:3] != ref.shape[:3] or not np.allclose(img.affine, ref.affine):
            img = image.resample_to_img(img, ref, interpolation="continuous")
        arrays.append(_zscore_run(img.get_fdata(dtype=np.float32)))

    concat = np.concatenate(arrays, axis=-1) if len(arrays) > 1 else arrays[0]
    log.info("    Concatenated: %d total volumes across %d run(s), per-run z-scored",
             concat.shape[-1], len(arrays))
    return nib.Nifti1Image(concat, ref.affine)


def run_melodic(bold_img: nib.Nifti1Image,
                mask_img: nib.Nifti1Image,
                out_dir: Path,
                n_dims: int,
                tr: float,
                bg_img: nib.Nifti1Image | None = None) -> Path | None:
    """
    Write BOLD (and mask, and optionally a background image) to temp files
    and run FSL MELODIC. tr must be this concatenated group's own real
    RepetitionTime (see get_run_tr()) -- rest and task runs can have
    genuinely different TRs.

    bg_img, if given, is passed as --bgimage so MELODIC's own native HTML
    report has a real anatomical underlay. Without it, MELODIC defaults to
    the mean of its input -- degenerate here, since the input has already
    been per-run z-scored to zero mean by concatenate_bold_runs().

    Returns the output directory, or None on failure.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        bold_path = Path(tmp) / "bold.nii.gz"
        mask_path = Path(tmp) / "mask.nii.gz"
        nib.save(bold_img, str(bold_path))
        nib.save(mask_img, str(mask_path))

        cmd = [
            MELODIC_BIN,
            "-i",           str(bold_path),
            "-o",           str(out_dir),
            "--mask="     + str(mask_path),
            "--dim="      + str(n_dims),
            "--tr="       + str(tr),
            "--nobet",
            "--report",
            "--Oall",
            "-v",
        ]

        if bg_img is not None:
            bg_path = Path(tmp) / "bg.nii.gz"
            nib.save(bg_img, str(bg_path))
            cmd.append("--bgimage=" + str(bg_path))

        log.info("    Running MELODIC (dim=%d) → %s", n_dims, out_dir.name)
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            log.error("    MELODIC failed:\n%s", result.stderr[-2000:])
            return None

    return out_dir


def load_melodic_components(melodic_out: Path) -> nib.Nifti1Image | None:
    """
    Load IC maps from a MELODIC output directory.
    Returns a 4D NIfTI (x, y, z, n_components), or None if not found.
    """
    candidates = sorted(melodic_out.glob("melodic_IC.nii.gz"))
    if not candidates:
        candidates = sorted(melodic_out.glob("*.nii.gz"))

    if not candidates:
        log.error("    No NIfTI output found in %s", melodic_out)
        return None

    img = nib.load(str(candidates[0]))
    if img.ndim == 3:
        data = img.get_fdata(dtype=np.float32)[..., np.newaxis]
        img  = nib.Nifti1Image(data, img.affine)
    return img


# ---------------------------------------------------------------------------
# Atlas loading — internal helpers
# ---------------------------------------------------------------------------

def _parse_fsl_xml(xml_path: Path) -> dict[str, int]:
    """
    Parse an FSL-format atlas XML file.
    Returns {label_name: nifti_value} where nifti_value = xml_index + 1
    (NIfTI background = 0; first label index=0 → value 1).
    """
    tree = ET.parse(str(xml_path))
    return {
        label.text.strip(): int(label.attrib["index"]) + 1
        for label in tree.findall(".//label")
    }


def _label_matches(label: str, queries: list[str]) -> bool:
    label_l = label.lower()
    return any(q.lower() in label_l for q in queries)


def _hemi_matches(label: str, hemi: str) -> bool:
    if hemi == "both":
        return True
    label_l = label.lower()
    if hemi == "L":
        return "left" in label_l or "lh_" in label_l or "lh " in label_l
    if hemi == "R":
        return "right" in label_l or "rh_" in label_l or "rh " in label_l
    return True


def _binary_mask(parc_img: nib.Nifti1Image,
                 selected_ids: list[int]) -> nib.Nifti1Image:
    data = np.round(parc_img.get_fdata()).astype(np.int32)
    mask = np.isin(data, selected_ids).astype(np.uint8)
    return nib.Nifti1Image(mask, parc_img.affine)


def _load_fsl_atlas(atlas_cfg: dict,
                    labels: list[str],
                    hemi: str) -> nib.Nifti1Image:
    """Load an FSL atlas (XML + NIfTI under $FSLDIR/data/atlases/)."""
    nii_path = FSL_ATLAS_DIR / atlas_cfg["path"]
    xml_path = FSL_ATLAS_DIR / atlas_cfg["xml"]

    if not nii_path.exists():
        raise FileNotFoundError(f"FSL atlas NIfTI not found: {nii_path}")
    if not xml_path.exists():
        raise FileNotFoundError(f"FSL atlas XML not found: {xml_path}")

    label_map = _parse_fsl_xml(xml_path)   # {name: nifti_value}
    selected = [
        v for name, v in label_map.items()
        if _label_matches(name, labels) and _hemi_matches(name, hemi)
    ]
    if not selected:
        raise ValueError(
            f"No labels matched {labels} (hemi={hemi}) in {xml_path.name}. "
            f"Available: {sorted(label_map)}"
        )
    log.debug("FSL atlas %s: matched %d regions", atlas_cfg["path"], len(selected))
    return _binary_mask(nib.load(str(nii_path)), selected)


def _load_nilearn_atlas(atlas_cfg: dict,
                        labels: list[str],
                        hemi: str) -> nib.Nifti1Image:
    """Load a nilearn-fetched atlas (Schaefer or Harvard-Oxford)."""
    fn_name = atlas_cfg["fetch_fn"]
    fn_args = atlas_cfg.get("fetch_args", {})
    fetch_fn = getattr(datasets, fn_name)
    atlas = fetch_fn(**fn_args)

    # atlas.maps may be a path string or a NiftiImage depending on nilearn version
    maps = atlas.maps
    parc_img = nib.load(maps) if isinstance(maps, str) else maps

    raw_labels = atlas.labels
    # Decode bytes if needed (older nilearn)
    str_labels = [
        lb.decode() if isinstance(lb, bytes) else lb
        for lb in raw_labels
    ]
    # str_labels[i] corresponds to NIfTI parcel value i
    selected = [
        i for i, lb in enumerate(str_labels)
        if _label_matches(lb, labels) and _hemi_matches(lb, hemi)
    ]
    if not selected:
        raise ValueError(
            f"No nilearn labels matched {labels} (hemi={hemi}) "
            f"in {fn_name}. Sample labels: {str_labels[:6]}"
        )
    log.debug("nilearn %s: matched %d parcels", fn_name, len(selected))
    return _binary_mask(parc_img, selected)


def _load_local_atlas(atlas_cfg: dict,
                      labels: list[str] | None,
                      hemi: str) -> nib.Nifti1Image:
    """
    Load a local atlas from ATLAS_DIR.
    If the atlas has an XML, use it for label-based selection.
    If labels is None, binarise the whole image.
    """
    nii_path = ATLAS_DIR / atlas_cfg["path"]
    if not nii_path.exists():
        raise FileNotFoundError(f"Local atlas not found: {nii_path}")

    img = nib.load(str(nii_path))
    data = img.get_fdata(dtype=np.float32)

    if labels is None:
        # Use entire non-zero volume as mask
        mask = (data != 0).astype(np.uint8)
    else:
        xml_rel = atlas_cfg.get("xml")
        if xml_rel:
            xml_path = ATLAS_DIR / xml_rel
            label_map = _parse_fsl_xml(xml_path)
            selected = [
                v for name, v in label_map.items()
                if _label_matches(name, labels) and _hemi_matches(name, hemi)
            ]
            if not selected:
                raise ValueError(
                    f"No labels matched {labels} in local atlas {nii_path.name}"
                )
            mask = np.isin(np.round(data).astype(np.int32), selected).astype(np.uint8)
        else:
            # No XML: threshold the image and optionally restrict hemisphere
            mask = (data > 0).astype(np.uint8)

    # Hemisphere masking in voxel space (x-axis: low=right, high=left in RAS)
    if hemi == "L":
        mid = mask.shape[0] // 2
        mask[mid:] = 0
    elif hemi == "R":
        mid = mask.shape[0] // 2
        mask[:mid] = 0

    return nib.Nifti1Image(mask, img.affine)


YEO17_NETWORK_NAMES = [
    "Visual_A", "Visual_B",
    "SomatoMotor_A", "SomatoMotor_B",
    "DorsalAttention_A", "DorsalAttention_B",
    "Salience_VentAttn_A", "Salience_VentAttn_B",
    "Limbic_A", "Limbic_B",
    "Control_A", "Control_B", "Control_C",
    "Default_A", "Default_B", "Default_C",
    "TemporalParietal",
]


def yeo17_subject_path(subject: str) -> Path:
    return SYNTHSEG_DIR / f"sub-{subject}" / f"sub-{subject}_Yeo17_subject_specific_MNI.nii.gz"


def get_yeo17_subject_masks(subject: str) -> list[tuple[str, nib.Nifti1Image]]:
    """
    Load the subject-specific Yeo17 atlas (produced by step0) and return
    one binary mask per network as [(name, mask_img), ...].
    Networks with zero voxels are logged and skipped.
    """
    path = yeo17_subject_path(subject)
    if not path.exists():
        raise FileNotFoundError(
            f"Subject-specific Yeo17 atlas not found for sub-{subject}: {path}\n"
            "Run step0_synthseg.sh first."
        )
    yeo_img  = nib.load(str(path))
    yeo_data = np.round(yeo_img.get_fdata()).astype(np.int32)

    masks = []
    for label_idx, name in enumerate(YEO17_NETWORK_NAMES, start=1):
        mask = (yeo_data == label_idx).astype(np.uint8)
        if not mask.any():
            log.warning("sub-%s: Yeo17 network '%s' (label=%d) is empty — skipped",
                        subject, name, label_idx)
            continue
        masks.append((name, nib.Nifti1Image(mask, yeo_img.affine)))

    log.info("sub-%s: loaded %d / %d Yeo17 networks", subject, len(masks), len(YEO17_NETWORK_NAMES))
    return masks


# DKT cortical labels that constitute the language network (bilateral).
# Left:  1018=pars opercularis, 1020=pars triangularis (Broca)
#        1030=superior temporal, 1031=supramarginal (Wernicke)
# Right: homologous labels (2018, 2020, 2030, 2031) included because
#        RSN masks are bilateral even for lateralised functions.
_LANGUAGE_SYNTHSEG_LABELS = [1018, 1020, 1030, 1031, 2018, 2020, 2030, 2031]


# Smith 2009 RSN10: volume index 6 (0-based) = auditory network; threshold z≥3
_SMITH_AUDITORY_IDX   = 6
_SMITH_AUDITORY_Z_THR = 3.0

# Harvard-Oxford label substrings that define the language network
# (Broca: IFG pars opercularis + triangularis; Wernicke: posterior STG + supramarginal)
_LANGUAGE_HO_QUERIES = [
    "pars triangularis",
    "pars opercularis",
    "superior temporal gyrus, posterior",
    "supramarginal gyrus",
]


def _smith_rsn10_img() -> nib.Nifti1Image:
    """Load Smith 2009 RSN10 4-D atlas (nilearn 0.13+ compatible)."""
    smith    = datasets.fetch_atlas_smith_2009()
    rsn_path = getattr(smith, "maps", None) or getattr(smith, "rsn10", None)
    if rsn_path is None:
        raise RuntimeError("fetch_atlas_smith_2009() has no 'maps' or 'rsn10' attribute")
    return nib.load(str(rsn_path))


def get_auditory_network_mask(subject: str) -> tuple[str, nib.Nifti1Image] | None:
    """
    Subject-specific auditory mask: Smith RSN10 component 6 (z≥3) ∩ subject GM.
    Saved to disk on first call; loaded from cache on subsequent calls.
    Returns ("Auditory", mask_img) or None on failure.
    """
    cached = auditory_mask_path(subject)
    if cached.exists():
        return ("Auditory", nib.load(str(cached)))

    try:
        rsn_img = _smith_rsn10_img()
    except Exception as exc:
        log.warning("Could not load Smith 2009 RSN atlas: %s", exc)
        return None

    aud_data   = rsn_img.get_fdata(dtype=np.float32)[..., _SMITH_AUDITORY_IDX]
    aud_binary = (aud_data >= _SMITH_AUDITORY_Z_THR).astype(np.uint8)
    if not aud_binary.any():
        log.warning("Smith auditory component empty after z≥%.1f", _SMITH_AUDITORY_Z_THR)
        return None

    aud_img = nib.Nifti1Image(aud_binary, rsn_img.affine)

    # Intersect with subject union GM for subject-specificity (from step0)
    result_img = aud_img
    gm_path    = SYNTHSEG_DIR / f"sub-{subject}" / f"sub-{subject}_union_GM.nii.gz"
    if gm_path.exists():
        gm_img = nib.load(str(gm_path))
        aud_r  = resample_to(aud_img, gm_img)
        masked = ((aud_r.get_fdata(dtype=np.float32) > 0)
                  & (gm_img.get_fdata(dtype=np.float32) > 0))
        if masked.any():
            result_img = nib.Nifti1Image(masked.astype(np.uint8), gm_img.affine)
        else:
            log.warning("sub-%s: Auditory mask empty after GM intersection — using group mask", subject)

    nib.save(result_img, str(cached))
    log.info("sub-%s: Auditory mask %d vox → %s",
             subject, int((result_img.get_fdata() > 0).sum()), cached.name)
    return ("Auditory", result_img)


def get_language_network_mask(subject: str) -> tuple[str, nib.Nifti1Image] | None:
    """
    Subject-specific language mask from SynthSeg DKT cortical labels (Broca + Wernicke).
    Saved to disk on first call; loaded from cache on subsequent calls.
    Returns ("Language", mask_img) or None if SynthSeg output is missing/empty.
    """
    cached = language_mask_path(subject)
    if cached.exists():
        return ("Language", nib.load(str(cached)))

    mask = _load_synthseg_mask(subject, _LANGUAGE_SYNTHSEG_LABELS)
    if mask is None:
        log.warning("sub-%s: language network mask could not be built — skipped", subject)
        return None

    nib.save(mask, str(cached))
    log.info("sub-%s: Language mask %d vox → %s",
             subject, int(mask.get_fdata().sum()), cached.name)
    return ("Language", mask)


def get_language_group_mask() -> nib.Nifti1Image | None:
    """
    Group-level bilateral language ROI mask from Harvard-Oxford cortical atlas
    (Broca: IFG pars opercularis + triangularis; Wernicke: posterior STG + supramarginal).
    Used as the ROI mask in step4 normative comparison — analogous to Yeo17 group atlas.
    """
    try:
        atlas_cfg = {
            "fetch_fn":   "fetch_atlas_harvard_oxford",
            "fetch_args": {"atlas_name": "cort-maxprob-thr25-2mm"},
        }
        return _load_nilearn_atlas(atlas_cfg, _LANGUAGE_HO_QUERIES, "both")
    except Exception as exc:
        log.warning("Could not build Language group ROI mask: %s", exc)
        return None


def get_auditory_group_mask() -> nib.Nifti1Image | None:
    """
    Group-level auditory ROI mask: Smith RSN10 component 6 thresholded at z≥3,
    without subject-specific GM intersection.
    Used as the ROI mask in step4 normative comparison.
    """
    try:
        rsn_img    = _smith_rsn10_img()
        aud_data   = rsn_img.get_fdata(dtype=np.float32)[..., _SMITH_AUDITORY_IDX]
        aud_binary = (aud_data >= _SMITH_AUDITORY_Z_THR).astype(np.uint8)
        if not aud_binary.any():
            raise RuntimeError("Smith auditory component empty after thresholding")
        return nib.Nifti1Image(aud_binary, rsn_img.affine)
    except Exception as exc:
        log.warning("Could not build Auditory group ROI mask: %s", exc)
        return None


def _load_propagated_subject_atlas(atlas_cfg: dict,
                                   subject: str,
                                   labels: list[str],
                                   hemi: str) -> nib.Nifti1Image:
    """
    Load a subject-specific propagated atlas produced by step0.
    Uses the base FSL atlas XML for label name → integer value lookup,
    then extracts matching labels from the subject's propagated image.
    """
    filename  = atlas_cfg["filename"].format(subject=subject)
    subj_path = SYNTHSEG_DIR / f"sub-{subject}" / filename
    if not subj_path.exists():
        raise FileNotFoundError(
            f"Subject-specific propagated atlas not found: {subj_path}\n"
            "Run step0_synthseg.sh first."
        )

    base_key = atlas_cfg["base_atlas"]
    base_cfg = ATLASES[base_key]
    xml_path = FSL_ATLAS_DIR / base_cfg["xml"]
    if not xml_path.exists():
        raise FileNotFoundError(f"Base atlas XML not found: {xml_path}")

    label_map = _parse_fsl_xml(xml_path)
    selected  = [
        v for name, v in label_map.items()
        if _label_matches(name, labels) and _hemi_matches(name, hemi)
    ]
    if not selected:
        raise ValueError(
            f"No labels matched {labels} (hemi={hemi}) in {xml_path.name}. "
            f"Available: {sorted(label_map)}"
        )

    img  = nib.load(str(subj_path))
    mask = _binary_mask(img, selected)

    # Hemisphere masking (x-axis: low=right, high=left in RAS)
    mask_data = mask.get_fdata(dtype=np.float32)
    if hemi == "L":
        mid = mask_data.shape[0] // 2
        mask_data[mid:] = 0
    elif hemi == "R":
        mid = mask_data.shape[0] // 2
        mask_data[:mid] = 0
    return nib.Nifti1Image(mask_data.astype(np.uint8), img.affine)


def _load_synthseg_mask(subject: str,
                        synthseg_label_ids: list[int]) -> nib.Nifti1Image | None:
    p = synthseg_path(subject)
    if not p.exists():
        log.warning("SynthSeg output missing for sub-%s — skipping", subject)
        return None
    seg = nib.load(str(p))
    data = np.round(seg.get_fdata()).astype(np.int32)
    mask = np.isin(data, synthseg_label_ids).astype(np.uint8)
    if not mask.any():
        log.warning("SynthSeg labels %s empty for sub-%s", synthseg_label_ids, subject)
        return None
    return nib.Nifti1Image(mask, seg.affine)


# ---------------------------------------------------------------------------
# Public: get seed mask
# ---------------------------------------------------------------------------

def resample_to(src: nib.Nifti1Image,
                ref: nib.Nifti1Image,
                interpolation: str = "nearest") -> nib.Nifti1Image:
    return image.resample_to_img(src, ref, interpolation=interpolation, copy=True)


def get_seed_mask(seed_name: str,
                  seed_cfg: dict,
                  subject: str | None = None) -> nib.Nifti1Image:
    """
    Build a binary seed mask for one entry from the seed catalog.

    seed_cfg keys (from profiles.yaml):
      atlas   : atlas key in atlases.yaml
      labels  : list of label substrings to match
      hemi    : "both" | "L" | "R"

    If subject is provided and the atlas config contains synthseg_labels,
    the SynthSeg parcellation is unioned with the atlas mask (subcortical only).

    No GM masking is applied here — seeds are kept as fixed anatomical
    definitions for cross-subject comparability.
    """
    atlas_key = seed_cfg["atlas"]
    labels    = seed_cfg.get("labels")
    hemi      = seed_cfg.get("hemi", "both")

    atlas_cfg = ATLASES[atlas_key]
    source    = atlas_cfg["source"]

    if source == "fsl":
        mask_img = _load_fsl_atlas(atlas_cfg, labels, hemi)
    elif source == "nilearn":
        mask_img = _load_nilearn_atlas(atlas_cfg, labels, hemi)
    elif source == "local":
        mask_img = _load_local_atlas(atlas_cfg, labels, hemi)
    elif source == "subject_synthseg":
        if subject is None:
            raise ValueError(f"Atlas '{atlas_key}' (subject_synthseg) requires a subject ID")
        ss_labels = seed_cfg.get("synthseg_labels")
        if not ss_labels:
            raise ValueError(f"Seed using atlas '{atlas_key}' must specify synthseg_labels")
        mask_img = _load_synthseg_mask(subject, ss_labels)
        if mask_img is None:
            raise ValueError(f"SynthSeg mask empty for sub-{subject} labels={ss_labels}")
    elif source == "subject_propagated":
        if subject is None:
            raise ValueError(f"Atlas '{atlas_key}' (subject_propagated) requires a subject ID")
        mask_img = _load_propagated_subject_atlas(atlas_cfg, subject, labels, hemi)
    else:
        raise ValueError(f"Unknown atlas source '{source}' for '{atlas_key}'")

    # Optionally augment with subject-specific SynthSeg mask (subcortical seeds)
    if subject and "synthseg_labels" in seed_cfg:
        ss = _load_synthseg_mask(subject, seed_cfg["synthseg_labels"])
        if ss is not None:
            ss_r  = resample_to(ss, mask_img)
            union = ((mask_img.get_fdata() > 0) | (ss_r.get_fdata() > 0)).astype(np.uint8)
            mask_img = nib.Nifti1Image(union, mask_img.affine)

    n_vox = int((mask_img.get_fdata() > 0).sum())
    log.info("Seed '%s': %d voxels (atlas=%s, hemi=%s)", seed_name, n_vox, atlas_key, hemi)
    return mask_img


# ---------------------------------------------------------------------------
# Timeseries extraction & FC
# ---------------------------------------------------------------------------

def extract_mean_ts(bold_data: np.ndarray,
                    roi_img: nib.Nifti1Image,
                    bold_ref: nib.Nifti1Image) -> np.ndarray:
    """
    Mean BOLD timeseries within roi_img (resampled to BOLD space).
    bold_data: pre-loaded float32 array (x, y, z, t).
    bold_ref:  the source NIfTI image, used only as spatial reference for resampling.
    Returns array of shape (n_trs,).
    """
    roi_r    = resample_to(roi_img, bold_ref)
    roi_mask = roi_r.get_fdata(dtype=np.float32) > 0

    if not roi_mask.any():
        raise ValueError("ROI mask is empty after resampling to BOLD space.")

    return bold_data[roi_mask].mean(axis=0)   # (n_trs,)


def whole_brain_fc(bold_data: np.ndarray,
                   seed_ts: np.ndarray,
                   brain_mask_img: nib.Nifti1Image,
                   bold_ref: nib.Nifti1Image) -> nib.Nifti1Image:
    """
    Pearson r between seed_ts and every in-mask voxel, Fisher-Z transformed.
    bold_data: pre-loaded float32 array (x, y, z, t).
    bold_ref:  the source NIfTI image, used as spatial reference and for output affine.
    Returns a 3D NIfTI map in the same space as bold_ref.
    """
    mask_data = brain_mask_img.get_fdata(dtype=np.float32) > 0
    shape3d   = bold_data.shape[:3]

    X  = bold_data[mask_data]                           # (n_vox, t)
    X -= X.mean(axis=1, keepdims=True)                  # demean voxels
    y  = seed_ts - seed_ts.mean()                       # demean seed

    num   = X @ y                                       # (n_vox,)
    denom = np.linalg.norm(X, axis=1) * np.linalg.norm(y) + 1e-12
    r     = np.clip(num / denom, -0.9999, 0.9999)
    z     = np.arctanh(r)                               # Fisher Z

    z_vol          = np.zeros(shape3d, dtype=np.float32)
    z_vol[mask_data] = z
    return nib.Nifti1Image(z_vol, bold_ref.affine)


def average_maps(imgs: list[nib.Nifti1Image],
                 dofs: list[int]) -> tuple[nib.Nifti1Image, int]:
    """Voxelwise mean of a list of Fisher-Z FC maps, plus their pooled
    degrees of freedom for downstream significance testing.

    dofs[i] = n_timepoints[i] - 3, one per image in `imgs` (same order).
    Pooled dof is their sum -- treating the runs as independent estimates
    of the same underlying correlation (a fixed-effects approximation;
    the point estimate itself stays an unweighted mean of the per-run
    Fisher-Z values, so this is not a full inverse-variance-weighted
    meta-analytic combination, just a DOF-aware significance test on top
    of the existing averaging).

    Runs are usually already on the same grid (fmriprep run with a pinned
    res-2 output resolution puts every task on one common grid). Without
    res-2, tasks acquired at different native resolutions (e.g. a
    resting-state run with different slice thickness than task runs) can
    differ in shape even within the same named MNI space, so resample
    onto the first image's grid whenever geometry doesn't already match.

    Uses nearest-neighbor interpolation, not continuous/trilinear: a Fisher-Z
    FC map is a per-voxel independent statistical estimate, not a smoothly
    varying field, so there's no basis for blending it across voxels in the
    first place. Confirmed empirically that "continuous" is actively harmful
    here -- resampling a mostly-exact-zero background with trilinear
    interpolation introduces floating-point noise (values as small as ~1e-25,
    not real weak correlation) across most of the volume, not just at the
    true signal boundary, contaminating any downstream statistic computed
    over "nonzero" voxels (e.g. FDR family size, robust noise-scale estimates).
    """
    ref = imgs[0]
    resampled = [
        im if (im.shape == ref.shape and np.allclose(im.affine, ref.affine))
        else image.resample_to_img(im, ref, interpolation="nearest")
        for im in imgs
    ]
    stack = np.stack([im.get_fdata(dtype=np.float32) for im in resampled], axis=0)
    pooled_dof = int(sum(dofs))
    return nib.Nifti1Image(stack.mean(axis=0), ref.affine), pooled_dof


def pooled_dof_for_runs(subject: str, run_ids: list[str]) -> int:
    """
    Sum of (n_timepoints - 3) across a subject's denoised BOLD runs, for
    recomputing the pooled dof used by threshold_fc_map() when only a
    cached (already-averaged) Fisher-Z map is available -- e.g. when
    regenerating a missing thresholded map without re-running the full
    per-run FC computation. nib.load() is lazy, so reading .shape here
    doesn't load actual voxel data.
    """
    dof = 0
    for run_id in run_ids:
        try:
            dof += load_bold(subject, run_id).shape[-1] - 3
        except Exception as exc:
            log.warning("pooled_dof_for_runs: could not read %s/%s (%s)",
                        subject, run_id, exc)
    return dof


def write_provenance_json(nii_path: Path, **fields) -> None:
    """
    Write a small <name>.json sidecar next to a statistical-map NIfTI,
    recording how it was made (dof, fdr_q, min_cluster, min_r_relative,
    which runs/run-type contributed, when) -- so "where did this map come
    from" is answerable by looking at the file next to it, not by reading
    pipeline source or trusting memory of what parameters were last used.
    """
    name = nii_path.name
    if name.endswith(".nii.gz"):
        name = name[: -len(".nii.gz")]
    json_path = nii_path.parent / f"{name}.json"
    with open(json_path, "w") as fh:
        json.dump(fields, fh, indent=2, default=str)


def _zstat_and_p(data: np.ndarray, dof: int, bilateral: bool) -> tuple[np.ndarray, np.ndarray]:
    """Shared by threshold_fc_map() and stat_maps() so both always agree."""
    se = 1.0 / np.sqrt(max(dof, 1))
    z_stat = data / se
    if bilateral:
        p = 2.0 * norm.sf(np.abs(z_stat))
    else:
        p = norm.sf(z_stat)
    return z_stat, p


def stat_maps(img: nib.Nifti1Image, dof: int,
             bilateral: bool = False) -> tuple[nib.Nifti1Image, nib.Nifti1Image]:
    """
    Return (z-statistic map, -log10(p) map) for a Fisher-Z FC map, unthresholded
    -- the full statistical map for direct inspection (e.g. in mrview), the way
    the tbfMRI GLM path saves its unthresholded T-map alongside the thresholded
    variants, rather than only ever exposing the binary pass/fail result.
    Uses the same DOF-aware test as threshold_fc_map() (see there for the
    FDR/min_r_relative rationale); this just skips the pass/fail step.
    """
    data = img.get_fdata(dtype=np.float32)
    z_stat, p = _zstat_and_p(data, dof, bilateral)
    neglog10p = -np.log10(np.clip(p, 1e-300, 1.0)).astype(np.float32)
    return (
        nib.Nifti1Image(z_stat.astype(np.float32), img.affine, img.header),
        nib.Nifti1Image(neglog10p, img.affine, img.header),
    )


def fdr_threshold(p: np.ndarray, mask: np.ndarray, q: float = 0.05) -> float:
    """
    Benjamini-Hochberg FDR-adjusted p-value cutoff, computed only over
    voxels within `mask` (the actual tested search volume -- untested
    background voxels aren't part of the multiple-comparisons family and
    would dilute it if included).

    Standard BH procedure: sort p ascending, find the largest k with
    p_(k) <= (k/m)*q, reject everything at or below that p-value. Unlike a
    fixed uncorrected threshold, the resulting cutoff adapts to how much
    real signal this specific map actually has -- a map with a lot of true
    connectivity gets to keep more of it at the same false-discovery rate;
    a mostly-null map gets pushed tighter automatically instead of passing
    the same fixed bar regardless of what's actually in it.

    Returns the raw-p cutoff to threshold at, or 0.0 if nothing survives at
    the requested FDR level (in which case no voxel should be kept).
    """
    pvals = np.sort(p[mask].ravel())
    m = pvals.size
    if m == 0:
        return 0.0
    ranks = np.arange(1, m + 1)
    below = pvals <= (ranks / m) * q
    if not below.any():
        return 0.0
    return float(pvals[below].max())


def threshold_fc_map(img: nib.Nifti1Image,
                     dof: int,
                     fdr_q: float = 0.05,
                     min_cluster: int = 50,
                     bilateral: bool = False,
                     min_r_relative: float = 2.0) -> nib.Nifti1Image:
    """
    Threshold a Fisher-Z FC map by voxelwise significance + cluster extent,
    mirroring the tbfMRI GLM convention (corrected voxel significance +
    minimum cluster size) rather than a fixed effect-size cutoff.

    Under H0 (true correlation = 0), Fisher-Z has SE = 1/sqrt(dof), so
    z_stat = data / SE is a standard normal deviate; p-values follow from
    the normal survival function. dof must be the pooled degrees of
    freedom behind `img` (see average_maps()/pooled_dof_for_runs()) --
    this is what makes the test DOF-aware instead of an arbitrary fixed
    |r| cutoff applied identically regardless of how many timepoints went
    into the estimate.

    Significance is corrected via Benjamini-Hochberg FDR (see
    fdr_threshold()) rather than a fixed uncorrected p, since a single fixed
    p bar performs very differently depending on how much real signal a
    given map has (the exact problem raised: a fixed r or p threshold that
    looks right for one seed/RSN can be far too permissive or too strict for
    another).

    Significance alone still isn't sufficient on its own: with enough pooled
    timepoints (this pipeline commonly pools 300-700+ TRs across runs), even
    a weak r~0.1-0.15 can be statistically reliable without being a
    clinically meaningful connection. Rather than a single fixed |r| floor
    across every map (same problem as above), `min_r_relative` expresses the
    effect-size floor as a multiple of THIS map's own robust noise scale
    (1.4826 * MAD of its Fisher-Z values, a normal-equivalent robust SD) --
    a map with a naturally broad spread of background correlation needs
    proportionally more effect size to count as real; a cleaner map's
    absolute-r floor ends up lower. BOLD's temporal autocorrelation still
    means the true independent sample size is smaller than raw TR count
    suggests (not corrected for here).

    dof:            pooled degrees of freedom (sum of n_timepoints-3 across
                    the runs averaged into `img`).
    fdr_q:          Benjamini-Hochberg FDR level (default 0.05).
    min_cluster:    minimum cluster size in voxels; smaller clusters are zeroed.
    bilateral:      if True, two-tailed test (keep |z_stat| significant in
                    either direction); if False (default), one-tailed test for
                    positive correlation only (negative FC is never significant
                    under this test regardless of magnitude).
    min_r_relative: effect-size floor as a multiple of this map's own robust
                    noise scale (default 2.0); set to 0 to disable and use
                    significance alone.
    Returns a new NIfTI with sub-threshold and small-cluster voxels set to 0.
    """
    data = img.get_fdata(dtype=np.float32)
    z_stat, p = _zstat_and_p(data, dof, bilateral)

    # NOTE: `data != 0` includes a large population of numerically-degenerate
    # near-zero values (float32 subnormals up through ~1e-3) that are not
    # real weak correlation -- see conversation/investigation notes. A naive
    # `abs(data) > eps` exclusion was tried and overcorrected (0 survivors
    # everywhere), so reverted to this pending root-cause investigation
    # rather than guessing another cutoff.
    brain = data != 0

    p_cut = fdr_threshold(p, brain, q=fdr_q)

    # MAD anchored at 0 (the null-hypothesis center for "no correlation"),
    # NOT at this map's own observed median -- a map whose median sits well
    # away from 0 (e.g. widespread real background correlation, not noise)
    # would otherwise have its "noise scale" measured as spread around that
    # already-elevated median, silently producing a much higher effective
    # floor than intended (confirmed on real data: median r=0.31 pushed a 5x
    # floor up to r=0.90, passing almost nothing).
    z_vals = data[brain]
    if z_vals.size > 0:
        robust_sigma = 1.4826 * float(np.median(np.abs(z_vals)))
    else:
        robust_sigma = 0.0
    min_r_z = min_r_relative * robust_sigma

    effect_ok = np.abs(data) >= min_r_z
    above = (p <= p_cut) & effect_ok & brain

    labeled, n_clusters = ndimage_label(above)
    keep = np.zeros_like(above)
    for i in range(1, n_clusters + 1):
        if (labeled == i).sum() >= min_cluster:
            keep[labeled == i] = True

    out = np.where(keep, data, 0.0).astype(np.float32)
    log.debug(
        "threshold_fc_map: FDR q=%.3f (p<=%.2e), dof=%d, min_r_relative=%.1fx robust-sigma (z>=%.4f, r>=%.3f), "
        "min_cluster=%d, bilateral=%s → %d/%d clusters kept",
        fdr_q, p_cut, dof, min_r_relative, min_r_z, np.tanh(min_r_z), min_cluster, bilateral,
        int(keep.any(axis=(0, 1)).sum()),  # rough cluster count proxy
        n_clusters,
    )
    return nib.Nifti1Image(out, img.affine, img.header)


# ---------------------------------------------------------------------------
# Combining independent z-statistic maps (e.g. masked-ICA rest vs task)
# ---------------------------------------------------------------------------

def combine_zmaps(z1: np.ndarray, dof1: float,
                  z2: np.ndarray | None = None, dof2: float | None = None) -> np.ndarray:
    """
    Stouffer-combine one or two independent z-statistic maps into a single
    z-map, weighting each by sqrt(dof) -- the standard way to pool z-scores
    from independent analyses that test the same hypothesis with different
    sample sizes. Falls back to z1 unchanged if only one map is available
    (e.g. a subject with no tbfMRI runs).

    Unlike threshold_fc_map()'s inputs (Fisher-Z of a Pearson r, needing the
    SE=1/sqrt(dof) rescale in _zstat_and_p), z1/z2 here are assumed to
    already be z-statistics (e.g. MELODIC's own spatial IC maps), so no
    rescale is applied before combining.
    """
    if z2 is None:
        return z1
    w1 = np.sqrt(max(dof1, 1.0))
    w2 = np.sqrt(max(dof2, 1.0))
    return (w1 * z1 + w2 * z2) / np.sqrt(w1 ** 2 + w2 ** 2)


def threshold_zstat_map(img: nib.Nifti1Image, fdr_q: float = 0.05,
                        min_cluster: int = 50, bilateral: bool = True) -> nib.Nifti1Image:
    """
    FDR + cluster-extent threshold for a map that is already in z-statistic
    units (e.g. a MELODIC IC map or a combine_zmaps() output) -- as opposed
    to threshold_fc_map(), whose input is Fisher-Z of a correlation and
    needs a dof-driven rescale to get to z-statistic units first.
    """
    data = img.get_fdata(dtype=np.float32)
    brain = data != 0

    if bilateral:
        p = 2.0 * norm.sf(np.abs(data))
    else:
        p = norm.sf(data)
    p_cut = fdr_threshold(p, brain, q=fdr_q)

    above = (p <= p_cut) & brain
    labeled, n_clusters = ndimage_label(above)
    keep = np.zeros_like(above)
    for i in range(1, n_clusters + 1):
        if (labeled == i).sum() >= min_cluster:
            keep[labeled == i] = True

    out = np.where(keep, data, 0.0).astype(np.float32)
    return nib.Nifti1Image(out, img.affine, img.header)


# ---------------------------------------------------------------------------
# Union mask (seed ∪ thresholded SBA map)
# ---------------------------------------------------------------------------

def build_union_mask(seed_img: nib.Nifti1Image,
                     sba_img: nib.Nifti1Image,
                     bold_ref: nib.Nifti1Image,
                     dof: float,
                     fdr_q: float = 0.05,
                     min_cluster: int = 50,
                     bilateral: bool = False,
                     min_r_relative: float = 2.0) -> nib.Nifti1Image:
    """
    Construct the ICA search-space mask as:
        union(seed_mask, thresholded SBA map)
    where the SBA map is thresholded with the same FDR + relative-effect-size
    framework used elsewhere in the pipeline (see threshold_fc_map).
    Both inputs are resampled to bold_ref space before combining.
    Raises ValueError if the result has fewer than MIN_MASK_VOXELS voxels.
    """
    seed_r = resample_to(seed_img, bold_ref)
    sba_r  = resample_to(sba_img,  bold_ref)

    sba_thresh = threshold_fc_map(sba_r, dof=dof, fdr_q=fdr_q, min_cluster=min_cluster,
                                  bilateral=bilateral, min_r_relative=min_r_relative)

    seed_bin = seed_r.get_fdata(dtype=np.float32) > 0
    sba_bin  = sba_thresh.get_fdata(dtype=np.float32) != 0
    union    = (seed_bin | sba_bin).astype(np.uint8)

    n_vox = int(union.sum())
    if n_vox < MIN_MASK_VOXELS:
        raise ValueError(
            f"Union mask has only {n_vox} voxels "
            f"(min={MIN_MASK_VOXELS}). "
            "Consider lowering min_r_relative or fdr_q."
        )
    log.info(
        "Union mask: %d voxels total  (seed=%d | SBA=%d)",
        n_vox, int(seed_bin.sum()), int(sba_bin.sum()),
    )
    return nib.Nifti1Image(union, bold_ref.affine)
