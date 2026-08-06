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
  -E   fMRI GLM engine: spm or nilearn (default: spm). nilearn needs no
         MATLAB/SPM12 license; it runs in the conda env $KUL_PYFMRI_ENV
         (default 'pyfMRI') — see "Conda environments" below.
  -N   opt-in: also run presurgical/eloquent-cortex rsfMRI network mapping
         (KUL_run_rsfMRI_networks.sh). Off by default. Boolean — takes no
         argument. Also runs in $KUL_PYFMRI_ENV.
  -C   condition profile for -N, from share/rsfmri_pipeline/config/profiles.yaml
         (default: Presurgical)
  -U   EXPERIMENTAL opt-in: prefer the rfa-modulated lore_sd FOD in KUL_FWT,
         if found, over the plain lore_sd ODF (only relevant with -D
         run_dwiprep_lore_sd.txt). Off by default.
  -Q   opt-in: run KUL_FWT's per-bundle tractometry (along-tract scalar
         profiles). Off by default — adds real per-bundle runtime. Tracts
         themselves are generated either way.
  -f   conda env to use instead of $KUL_SCILPY_ENV (default 'scilpy') for FWT
         (automated tractography). You shouldn't normally need this.
  -y   conda env to use instead of $KUL_PYFMRI_ENV (default 'pyfMRI') for -N
         and -E nilearn. You shouldn't normally need this either.
```

## Conda environments

Several steps need a conda env with specific Python/tractography packages. As
of the installer's `env-pyfmri`/`env-scilpy`/`env-lore-sd` sections, these are
created under **fixed, predictable names** — `pyfMRI`, `scilpy`, `lore_sd` —
so you do **not** need to tell the pipeline which env to use for normal runs.

> **The pipeline will quit if the expected env is missing.** If FWT (`-f`),
> `-N`, `-E nilearn`, or a `-D` config requesting `lore_sd` can't find their
> conda env (`scilpy`, `pyfMRI`, `pyfMRI`, `lore_sd` respectively) via `conda
> env list`, the script exits immediately with an error naming the missing
> env — it does **not** silently fall back to a bare `python3` or hang. Run
> the corresponding `setup_environment.sh` section (e.g. `--only env-pyfmri`)
> to create it, or override with `-f`/`-y`/`loresd_env:` in the `-D` config
> if you're intentionally using a differently-named env.

## Examples

Paths and participant names below are fictitious — substitute your own study
root layout and DICOM/BIDS naming.

**0. First time for this patient (no study-root folder yet):**

```
mkdir -p /data/studies/glioma_batch3
cd /data/studies/glioma_batch3
KUL_clinical_fmridti.sh -p JaneDoe -s -t 1
```
`-s` does **not** run the pipeline — it creates a new
`clinical_sub-JaneDoe_type1/` folder (named from `-p`/`-t`) containing a
`DICOM/` and a `study_config/` pre-filled with the template configs for that
processing type, then exits. Pass `-t` explicitly if you're scaffolding for
anything other than the default type 1 — it selects which template configs
get copied in (see the type table above). Drop your DICOMs into
`clinical_sub-JaneDoe_type1/DICOM/`, adjust `study_config/` if needed, `cd`
into that folder, then run again without `-s` (example 1 below).

**1. Basic glioma work-up (type 1, default SPM engine), from a zip archive:**

```
cd /data/studies/glioma_batch3/clinical_sub-JaneDoe_type1
KUL_clinical_fmridti.sh -p JaneDoe -d DICOM/JaneDoe.zip -n 32
```
Uses defaults throughout: `-t 1` (intra-axial glioma stream), `-E spm`
(MATLAB/SPM12 GLM), FreeSurfer 8.2.0 recon-all, `scilpy`/`lore_sd` envs only
pulled in if the relevant steps need them.

**2. DBS case (type 3, manual lesion mask) with the nilearn GLM engine,
lore_sd-based FOD estimation, rsfMRI network mapping, and per-bundle
tractometry — a fuller "everything on" run:**

```
cd /data/studies/dbs_cohort/clinical_sub-M0012_type3
mkdir -p RESULTS/sub-M0012/Lesion
cp /data/segmentations/M0012_lesion.nii.gz RESULTS/sub-M0012/Lesion/lesion.nii.gz
KUL_clinical_fmridti.sh -p M0012 -d ./DICOM/M0012 -t 3 \
    -E nilearn -D run_dwiprep_lore_sd.txt -U -N -Q -n 32 -v 2
```
- `-t 3`: tumor with manual mask (the `lesion.nii.gz` you copied in above)
- `-E nilearn`: task-fMRI GLM via nilearn instead of SPM12/MATLAB — runs in
  the `pyfMRI` conda env automatically, no `-y` needed
- `-D run_dwiprep_lore_sd.txt`: use the lore_sd dwiprep config (runs in the
  `lore_sd` conda env automatically — the config's `loresd_env:` field can
  stay blank)
- `-U`: prefer the rfa-modulated lore_sd FOD in tractography, if present
- `-N`: also run rsfMRI network mapping (also uses `pyfMRI` automatically)
- `-Q`: run KUL_FWT's per-bundle tractometry (adds runtime)
- `-v 2`: verbose logging

Neither example passes `-f`, `-y`, or a `loresd_env:` override — the
`scilpy`/`pyfMRI`/`lore_sd` envs are found automatically by their fixed
names. You'd only add those if you deliberately installed one under a
different name.

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

This pipeline ties together most of KUL_NIS, so it needs the full software stack — see the **Requirements** section of the main [README](/README.md). The key external tools it invokes are: dcm2bids/dcm2niix, ANTs, FSL, FreeSurfer (8.2.0, or FastSurfer with `-X`), fmriprep, MRtrix3 (**3.0.8-2097-g99963980**, `dev` branch, built with CMake+Ninja), SPM12 (MATLAB), synb0-disco, hd-glio-auto, resseg, MSBP, **KUL_VBG**, **KUL_FWT** (via the `scilpy` conda env, see "Conda environments" above), Karawun, and `xvfb-run` (**required** for headless `mrview` screenshots — see below). `KUL_nii2dcm.py` additionally needs python3 with SimpleITK, Pillow and numpy.

### Headless `mrview` screenshots (`xvfb-run`)

Every figure/PACS screenshot in this pipeline is captured by running `mrview` off-screen via `xvfb-run`. `xvfb-run` must be installed and on the `PATH` (`sudo apt install -y xvfb libgl1-mesa-dri`); the script checks for it at the point screenshots are generated and prints a warning if it is missing, but will otherwise fail when it tries to render figures. On systems where MATLAB/MCR (used by SPM12) is also installed, the `mrview` calls are additionally run through a filtered environment that strips the MCR's bundled Qt5 from `LD_LIBRARY_PATH` and forces software rendering (`LIBGL_ALWAYS_SOFTWARE=1`), otherwise `mrview` can abort with `Could not find the Qt platform plugin xcb/offscreen`. Override the `mrview` binary with `MRVIEW_BIN=/path/to/mrview` if `PATH` resolution is ambiguous.
