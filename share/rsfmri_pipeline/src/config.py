"""
Pipeline configuration — reads atlases.yaml and profiles.yaml, sets all paths.
Edit the YAML files under config/ to change atlases, seeds, or profiles.
Edit the constants below to change paths or imaging parameters.
"""

import os
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
def _env_path(var: str, default: Path) -> Path:
    val = os.environ.get(var)
    return Path(val) if val else default

# Fallback for anyone still running this tree standalone (outside KUL_NIS);
# KUL_run_rsfMRI_networks.sh always sets the RSFMRI_* env vars explicitly.
BASE_DIR     = Path(os.environ.get("RSFMRI_BASE_DIR", "/home/radwan/Documents/Work_for_Nada"))
FMRIPREP_DIR = _env_path("RSFMRI_FMRIPREP_DIR", BASE_DIR / "fmriprep")
DENOISED_DIR = _env_path("RSFMRI_DENOISED_DIR", BASE_DIR / "denoised")
ANALYSIS_DIR = _env_path("RSFMRI_ANALYSIS_DIR", BASE_DIR / "analysis")

PIPELINE_DIR = Path(__file__).resolve().parent.parent   # rsfmri_pipeline/
CONFIG_DIR   = PIPELINE_DIR / "config"
ATLAS_DIR    = PIPELINE_DIR / "templates"

def _resolve_fsldir() -> Path:
    fsldir = os.environ.get("FSLDIR")
    if fsldir:
        return Path(fsldir)
    import shutil, subprocess
    flirt = shutil.which("flirt")
    if flirt:
        return Path(flirt).resolve().parent.parent
    raise EnvironmentError(
        "FSLDIR is not set and flirt not found on PATH. "
        "Source your FSL setup script before running the pipeline."
    )

FSLDIR        = _resolve_fsldir()
FSL_ATLAS_DIR = FSLDIR / "data" / "atlases"

# Runtime output directories (created on import)
SYNTHSEG_DIR  = ANALYSIS_DIR / "synthseg"
SBA_DIR       = ANALYSIS_DIR / "sba"
RSN_FC_DIR    = ANALYSIS_DIR / "rsn_fc"
ICA_DIR       = ANALYSIS_DIR / "masked_ica"
WBICA_DIR     = ANALYSIS_DIR / "whole_brain_ica"

for _d in [SYNTHSEG_DIR, SBA_DIR, RSN_FC_DIR, ICA_DIR, WBICA_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Subjects
# ---------------------------------------------------------------------------
# HV10 excluded: run-02 is 3D (single volume), caused fMRIPrep to crash.
HV_SUBJECTS  = [f"HV{i:02d}" for i in range(1, 10)]   # HV01–HV09
PT_SUBJECTS  = ["PT01", "PT02"]
ALL_SUBJECTS = HV_SUBJECTS + PT_SUBJECTS

# ---------------------------------------------------------------------------
# Imaging parameters
# ---------------------------------------------------------------------------
RES       = 2      # mm isotropic
MNI_SPACE = "MNI152NLin2009cAsym"

# ---------------------------------------------------------------------------
# Analysis thresholds
# ---------------------------------------------------------------------------
SMITH_Z_THRESHOLD = 1.5   # threshold applied to Smith 2009 soft RSN volumes
MIN_MASK_VOXELS  = 100    # minimum voxels for a valid union mask

# ---------------------------------------------------------------------------
# External binaries (must be on PATH)
# ---------------------------------------------------------------------------
MELODIC_BIN   = "melodic"
SYNTHSEG_BIN  = "mri_synthseg"
SYNTHSEG_THREADS = 4

# ---------------------------------------------------------------------------
# Load YAML configs
# ---------------------------------------------------------------------------
def _load_yaml(filename: str) -> dict:
    with open(CONFIG_DIR / filename) as fh:
        return yaml.safe_load(fh)

_atlas_data    = _load_yaml("atlases.yaml")
_profile_data  = _load_yaml("profiles.yaml")

ATLASES           = _atlas_data["atlases"]
SEED_CATALOG      = _profile_data["seeds"]
CONDITION_PROFILES = _profile_data["profiles"]

# ---------------------------------------------------------------------------
# Convenience accessors
# ---------------------------------------------------------------------------
def get_atlas(name: str) -> dict:
    if name not in ATLASES:
        raise ValueError(f"Unknown atlas '{name}'. Available: {sorted(ATLASES)}")
    return ATLASES[name]


def get_profile(name: str) -> dict:
    if name not in CONDITION_PROFILES:
        raise ValueError(
            f"Unknown profile '{name}'. "
            f"Available: {sorted(CONDITION_PROFILES)}"
        )
    return CONDITION_PROFILES[name]


def get_seed(name: str) -> dict:
    if name not in SEED_CATALOG:
        raise ValueError(f"Unknown seed '{name}'. Available: {sorted(SEED_CATALOG)}")
    return SEED_CATALOG[name]
