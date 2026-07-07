# KUL_clinical_fmridti

## Purpose

`KUL_clinical_fmridti.sh` is the **clinical batch pipeline** of KUL_NIS. From a single DICOM input it runs the full chain needed for a presurgical / DBS fMRI–dMRI work-up and produces review figures and (optionally) DICOMs for PACS and Brainlab.

It orchestrates the other KUL_NIS tools and the sibling repositories (KUL_VBG, KUL_FWT) so that, for a typical glioma case, one command takes you from dicom to tractography figures.

> **Clinical disclaimer:** any use in a clinical environment is off-label, not FDA-approved, not CE-labeled. See the License file.

## What it does (per processing type)

The `-t` option selects the processing stream:

| `-t` | Stream | Brain extraction / lesion handling |
|---|---|---|
| 1 (default) | Glioma, **intra-axial** — automatic tumor segmentation + VBG | hd-glio-auto segmentation, then KUL_VBG |
| 2 | Glioma, **extra-axial** — automatic tumor segmentation + VBG | hd-glio-auto segmentation, then KUL_VBG (extra-axial) |
| 3 | Tumor with a **manual mask** + VBG | put `lesion.nii.gz` in `RESULTS/sub-{participant}/Lesion` |
| 4 | **No glioma** (cavernoma, epilepsy, …) using cT1w | no VBG lesion filling |
| 5 | **DBS – essential tremor** (DRT tract) | dMRI DBS stream |
| 6 | **DBS – Parkinson's disease** (CSHD pathway) | dMRI DBS stream |
| 7 | **DTI-ALPS** (diffusion along perivascular spaces) processing stream | no tumor/VBG/tractography; see note below |

### DTI-ALPS (`-t 7`)

This stream computes the **DTI-ALPS index** (analysis along the perivascular space). It disables the tumor/VBG/DBS/tractography steps and instead:

1. copies the dedicated config from `study_config/DTI_ALPS_proc/` (sequences, dwiprep, fmriprep, freesurfer, gadolinium BIDS filter),
2. runs dMRI preprocessing and `KUL_dwiprep_MNI` (the ALPS index needs the data in a common space),
3. runs the ALPS computation via the external **`KUL_calc_DTIALPS.sh`**.

> **Requirement:** `KUL_calc_DTIALPS.sh` must be installed and on the `PATH`. It is **not** part of KUL_NIS_unified — provide it separately, otherwise `-t 7` will stop at the ALPS computation step.

Internally the chain (type-dependent) is roughly:

```
dcm2bids  →  anat registration  →  tumor segmentation (hd-glio-auto / resseg)
          →  KUL_VBG (lesion filling so FreeSurfer/FastSurfer can run)
          →  fmriprep (fMRI + anat)  →  SPM12 / FSL melodic activation maps
          →  KUL_dwiprep (+ synb0, MRtrix3)  →  KUL_FWT tractography
          →  mrview rendering of overlays/tracts  →  figures (and PACS/Karawun DICOMs)
```

## Usage

Run from the **study root** (the directory that holds `DICOM/`, `BIDS/`, `study_config/`, `RESULTS/`).

```
KUL_clinical_fmridti.sh -p JohnDoe -d DICOM/JohnDoe.zip
```

`-p` is the **only required** argument. On first use, `-s` will scaffold a default `DICOM/` and `study_config/` for you.

### Two-phase figure → PACS workflow

DICOM export to PACS is deliberately **not** automatic. The intended workflow is:

1. **Run the pipeline** (no `-R`) to generate the review figures (PNG screenshots of tracts and activation overlays).
2. **Review the figures** — check thresholds, orientations, that the right tracts/activations are shown.
3. **Re-run with `-R <underlay>`** to convert the reviewed figures into DICOMs for PACS and Karawun/Brainlab.

This is why `-R` is documented as "run this AFTER reviewing figures."

## Options

```
Required:
  -p   participant name (BIDS name, no underscores)

Optional:
  -t   processing type (1-7, see table above; default 1)
  -d   dicom zip file (or directory)
  -s   scaffold a default DICOM and study_config
  -B   make a backup and cleanup
  -r   redo certain steps (the program will ask)
  -R   generate DICOMs for PACS and Karawun (run AFTER reviewing figures)
         underlay choice: 1=cT1w  2=FLAIR  3=SWI  4=T1w  5=FGATIR  6=DIR  7=MP2RAGE(INV2)
  -O   orientations to render, comma-separated (default: TRA,SAG,COR)
  -e   add an edge outline to SPM/Melodic overlays (dark blue contour at the threshold boundary)
  -T   fixed threshold for ALL SPM/Melodic overlays (default: auto = max/3 per map;
         if omitted and interactive, you are prompted for one threshold per map)
  -a   opacity of SPM/Melodic (fMRI) activation overlays (0=transparent, 1=opaque; default 0.7)
         lower values let underlying anatomy show through on figures AND PACS DICOMs
  -D   dwiprep config file to use from study_config/ (default: run_dwiprep.txt)
         use e.g. -D run_dwiprep_lore_sd.txt to run lore-sd based FOD estimation
  -S   fMRI SUSAN smoothing FWHM in mm (default: adaptive = mean voxel size)
  -P   FWE-corrected p-value for Bizzi fMRI thresholding (default: 0.01)
  -n   number of threads (default 48)
  -v   verbosity (0=silent, 1=normal, 2=verbose; default 1)
  -X   use FastSurfer instead of plain recon-all for the reconstruction step
         in types 4, 5, 6 (faster, requires GPU; default is FreeSurfer 8.2.0 recon-all)
  -f   name of the conda environment to activate for scilpy-based tractography
         post-processing (KUL_FWT step). There is no default — if omitted, an
         empty environment name is passed to conda activate, which will fail;
         always pass -f (e.g. -f scilpy) unless your environment is already active.
```

### Figure / overlay appearance (`-a`, `-e`, `-T`)

These three flags control how the fMRI activation overlays look. Because each screenshot is rendered **once** and that same PNG feeds both the review figures and the PACS DICOMs, these settings apply to **both** outputs:

- **`-a` opacity** — default `0.7`. Lowering it (e.g. `-a 0.3`) makes activation more transparent so the underlying anatomy is visible. Be aware that very low opacity can make weak / threshold-level activation hard to read on PACS.
- **`-e` edge outline** — draws a solid dark-blue contour at the threshold boundary, which keeps cluster extent legible even at low opacity.
- **`-T` threshold** — fixes the activation threshold for every map (default is automatic, `max/3` per map).

> Tracts are rendered as 3D tractography lines (not an image overlay), so anatomy always shows through them regardless of `-a`.

## Outputs

- `RESULTS/sub-{participant}/` — anat, figures and intermediate results
- `*_figures_*/` — PNG screenshots for review (Tracto and SPM/fMRI), per underlay and orientation
- `RESULTS/.../PACS/` — DICOM series for PACS, created only with `-R` (via `KUL_nii2dcm.py`, using a donor DICOM for correct study/series linkage)
- `Karawun/sub-{participant}/` — Brainlab-compatible export

To push the PACS DICOMs to an Orthanc/PACS node, see `tools/send_2_orthanc.sh`.

## Dependencies

This pipeline ties together most of KUL_NIS, so it needs the full software stack — see the **Requirements** section of the main [README](/README.md). The key external tools it invokes are: dcm2bids/dcm2niix, ANTs, FSL, FreeSurfer (8.2.0, or FastSurfer with `-X`), fmriprep, MRtrix3 (**3.0.4-543-g86eb1ea8**, `dev` branch, 2023 build), SPM12 (MATLAB), synb0-disco, hd-glio-auto, resseg, MSBP, **KUL_VBG**, **KUL_FWT** (via the scilpy conda env, see `-f`), Karawun, and `xvfb-run` (**required** for headless `mrview` screenshots — see below). `KUL_nii2dcm.py` additionally needs python3 with SimpleITK, Pillow and numpy.

### Headless `mrview` screenshots (`xvfb-run`)

Every figure/PACS screenshot in this pipeline is captured by running `mrview` off-screen via `xvfb-run`. `xvfb-run` must be installed and on the `PATH` (`sudo apt install -y xvfb libgl1-mesa-dri`); the script checks for it at the point screenshots are generated and prints a warning if it is missing, but will otherwise fail when it tries to render figures. On systems where MATLAB/MCR (used by SPM12) is also installed, the `mrview` calls are additionally run through a filtered environment that strips the MCR's bundled Qt5 from `LD_LIBRARY_PATH` and forces software rendering (`LIBGL_ALWAYS_SOFTWARE=1`), otherwise `mrview` can abort with `Could not find the Qt platform plugin xcb/offscreen`. Override the `mrview` binary with `MRVIEW_BIN=/path/to/mrview` if `PATH` resolution is ambiguous.
