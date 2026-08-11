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

**Donor DICOM.** `-R` needs one donor DICOM — a single slice from a high-resolution
anatomical series is enough — to inherit study/series metadata. Two locations are
searched, `Karawun/sub-{participant}/DICOM/` first and `RESULTS/sub-{participant}/DICOM/`
second; both are created for you, so drop one file into either. `Karawun/` is the better
choice because the `importTractography` command printed at the end reads from there.

**Repeating `-R`.** Running `-R` again with a different underlay (e.g. `-R 4` then `-R 2`)
regenerates the PACS DICOMs for that underlay. Karawun prep is **not** repeated: its output
(tracts, T1w, labels) does not depend on the underlay, and it is gated by
`KUL_LOG/sub-{participant}_karawun_prepare.done`. Delete that marker to force a rebuild.
The marker is only written when Karawun prep actually succeeds.

**`-R` does not push to Brainlab.** When Karawun prep finishes it prints an
`importTractography` command to run manually (after `conda activate KarawunDev`), which
writes `Karawun/sub-{participant}/sub-{participant}_for_elements`.

### DSC perfusion

If a DSC series is present in `BIDS/sub-{participant}/perf/`, `KUL_dsc_perfusion.sh`
runs automatically — the same way fMRI and dMRI data are picked up — and `-W` skips it.
It is scheduled after VBG/multiparc, because the contralesional NAWM reference needs a
FreeSurfer `aseg`. See [KUL_dsc_perfusion](/docs/KUL_dsc_perfusion/KUL_dsc_perfusion.md).

With `-R`, the NAWM-normalised `nrCBV_corrected` and `nrCBF` maps are also exported as
PACS series into `PACS/Clinical_*`. These use **fixed** thresholds rather than the auto
`max/3` applied to the fMRI maps, because they are normalised ratios whose thresholds
carry a fixed clinical meaning: **1.75** on nrCBV (the conventional high-grade glioma
cutoff) and **1.0** on nrCBF ("above contralesional normal white matter"). `-T` overrides
both.

### Backup and cleanup (`-B`)

`-B` is a closing-out step, not something to run mid-analysis: it deletes the
regenerable intermediates, prompts for a password, archives everything left into
`../Finished_<date>_sub-{participant}_type{type}.7z`, and exits.

| removed | kept |
|---|---|
| `fmriprep_work*` | `BIDS/`, `fmriprep/`, `dwiprep/` derivatives |
| `BIDS/tmp_dcm2bids` | `RESULTS/`, `REPORT/` |
| dwiprep `raw/`, `dwi_orig*`, `*dwifsl*tmp*` | `Karawun/`, `KUL_LOG/` |
| the contents of dwiprep `dwi/` … | … **except `geomcorr.mif`** (+ its `.b`) |
| the denoised BOLD variants and SUSAN-smoothed GLM inputs | |

The denoised variants matter here: the task-melodic, resting-state and
pseudo-rest paths each keep their own copy of every 4D series (~6 GB on a
three-run subject), and the smoothed GLM inputs add ~2 GB. All are rebuilt from
the fmriprep output plus the confounds TSVs, so they are pure intermediates — but
before they were listed here they were silently archived into the `.7z`. On a
typical subject the cleanup removes roughly a third of the working directory
(~43 GB of 140 GB on a three-run clinical subject).

**Why `geomcorr.mif` survives.** `dwi/` is otherwise scratch, but `geomcorr.mif`
is what `dwifslpreproc` produced — the output of topup + eddy, the most expensive
step in dwiprep — and its own working directory (`*dwifsl*tmp*`, ~19 GB) is
deleted. Keeping it (1.3 GB) means bias correction, gradient checks and anything
downstream can be redone without re-running eddy. Note `dwi_preproced.mif` one
level up is the *post*-bias-correction volume, so it is not a substitute.
`geomcorr_grad_checked.b` — the gradient table `dwigradcheck` corrected post-eddy
— is kept alongside it, since the volume is not usable without it.

It is destructive and one-way: re-running afterwards has to redo the
denoise/smoothing steps.

### fMRI denoising: task vs resting-state

The three consumers of fMRI data want different things from denoising, and they now say so
explicitly rather than sharing one recipe (`KUL_fmri_denoise.sh`):

| consumer | recipe | why |
|---|---|---|
| task GLM (`-E nilearn` / `spm`) | no separate denoising; task + drift + confounds fit **simultaneously** in the GLM | the confound betas are estimated *controlling for* the task, so nuisance regression cannot eat the task |
| task melodic (`KUL_fmriproc_conn.sh`) | `--no-confounds --lp 0` — high-pass and smoothing only | separating task from artifact is what the ICA is *for* |
| rsfMRI networks (`-N`) | `--task-signal remove` on task runs | strips the evoked response so task runs can be pooled with genuine rest |

**Why melodic no longer confound-regresses.** The nuisance regressors are routinely
collinear with a block paradigm. On a 30 s on/off language run the 6 motion parameters alone
carried 27.4 % of the design's variance and aCompCor another 23.2 % — together 52.7 %. So
regressing them out before ICA removed most of the task with them, and melodic's best
component correlated only 0.27 with the paradigm. With regression skipped it finds the task
at 0.63–0.77, and that component's spatial map correlates 0.57 with the subject's own GLM
z-map. This is not about *how much* the patient moved — it is about *when*: small task-locked
movement is maximally damaging precisely because it is collinear.

`--task-signal preserve` also exists (joint model, confound betas estimated alongside the
design) but **nothing calls it**: it necessarily preserves task-correlated motion too, which
lights up a quarter to a third of the brain.

Resting-state runs are handled only by `KUL_run_rsfMRI_networks.sh`; `KUL_fmriproc_conn.sh`
skips them, since running melodic on the same data in both places was the same decomposition
computed twice.

### Lesion mask

When a lesion mask exists — `sub-{participant}_lesion_and_cavity.nii.gz` from
`KUL_anat_segment_tumor.sh` for types 1/2, or the manual `lesion.nii.gz` for type 3 — it is:

- binarised into `RESULTS/.../Anat/sub-{participant}_lesion.nii.gz`, next to the other
  T1w-space volumes that feed figures and PACS;
- exported as a Brainlab label (`Karawun/sub-{participant}/labels/Lesion.nii.gz`,
  palette colour 50), so the tumour appears in the same scene as the tracts;
- rendered as a PACS series in `PACS/Clinical_*` under `-R`.

### FAT1w QA volume for Brainlab

`KUL_karawun_prepare.sh` now generates `FAT1w.nii.gz` (= √FA · T1w, after Goedemans et al.,
*Imaging Neuroscience* 2024) via `KUL_FAT1w.py`, and loads it into Brainlab as a second
anatomical alongside the T1w. It makes the FA-to-T1w registration directly inspectable —
if that registration has slipped, every tract in the scene is displaced the same way and
nothing else in the export would reveal it.

This had regressed: nothing called `KUL_FAT1w.py` any more, so the file
`KUL_karawun_prepare.sh` looked for never existed and the QA volume was silently omitted
from every export. It needs `dwiprep/sub-X/sub-X/qa/fa_reg2T1w.nii.gz`, i.e.
`KUL_dwiprep_anat.sh` must have run; without it the step is skipped with a message.

## Options

```
Required:
  -p   participant name (BIDS name, no underscores)

Optional:
  -t   processing type (1-7, see table above; default 1)
  -d   dicom zip file (or directory)
  -s   scaffold a default DICOM and study_config
  -B   delete regenerable intermediates, then archive what remains into a
         password-protected ../Finished_<date>_sub-<p>_type<t>.7z and exit
         (see "Backup and cleanup" below). Destructive and one-way.
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
  -W   skip DSC perfusion processing (KUL_dsc_perfusion.sh). It otherwise runs
       automatically whenever a DSC series is present in BIDS/*/perf/
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
- `RESULTS/.../Anat/sub-{participant}_lesion.nii.gz` — the lesion mask (binarised) alongside the other T1w-space volumes, whichever processing type produced it
- `RESULTS/.../SPM/` — the task-fMRI GLM maps that are actually reviewed and exported. When both GLM variants run, the with-confounds (`_wc`) one is preferred; the plain one is used when `_wc` produced nothing
- `RESULTS/.../SPM_all/` — the complete record: every GLM variant × threshold combination
- `RESULTS/.../Tracto/` — final bundles as MRtrix `.tck`, plus tract density maps regridded to T1w
- `RESULTS/.../TRK/` — the same bundles as TrackVis `.trk`, for **freeview** (which does not read `.tck`). Written by `KUL_FWT_make_TCKs.sh` next to each `.tck` and copied here; the conversion reference is the FA/tracking space, not the FreeSurfer parcellation — see the `.trk` comment in that script for why that distinction matters when the dMRI and T1w come from different sessions
- `RESULTS/.../Perfusion/` — DSC perfusion maps and lesion/NAWM ratios, when DSC data is present (see [KUL_dsc_perfusion](/docs/KUL_dsc_perfusion/KUL_dsc_perfusion.md))
- `REPORT/sub-{participant}_06_Tract_QQ/` — with `-Q`, one self-contained interactive HTML per bundle (drag-to-rotate metric tower, bundle geometry, along-tract profiles for every metric, and an endpoint-connectivity table + heatmap naming the parcels the bundle joins), plus `sub-{participant}_FWT_report.html`: a single page with every bundle's screenshots, switchable between the four renderings, filterable by bundle name, each linking to its own detail page
- `*_figures_*/` — PNG screenshots for review (Tracto, SPM/fMRI and Clinical), per underlay and orientation
- `RESULTS/.../PACS/` — DICOM series for PACS, created only with `-R` (via `KUL_nii2dcm.py`, using a donor DICOM for correct study/series linkage). `PACS/Clinical_*` holds the lesion and DSC perfusion series; `PACS/fMRI_*` and `PACS/Tracto_*` hold the activation and tract series.
- `Karawun/sub-{participant}/` — Brainlab-compatible export, including a `FAT1w.nii.gz` QA volume, a `labels/Lesion.nii.gz` label when a lesion mask exists, and the fMRI activation labels (colours 51-63). See [KUL_karawun_prepare](/docs/KUL_karawun_prepare/KUL_karawun_prepare.md) for the label-colour convention — a label's voxel value *is* its Brainlab colour, and the ranges are non-overlapping by design

To push the PACS DICOMs to an Orthanc/PACS node, see `tools/send_2_orthanc.sh`.

## Dependencies

This pipeline ties together most of KUL_NIS, so it needs the full software stack — see the **Requirements** section of the main [README](/README.md). The key external tools it invokes are: dcm2bids/dcm2niix, ANTs, FSL, FreeSurfer (8.2.0, or FastSurfer with `-X`), fmriprep, MRtrix3 (**3.0.8-2097-g99963980**, `dev` branch, built with CMake+Ninja), SPM12 (MATLAB), synb0-disco, hd-glio-auto, resseg, MSBP, **KUL_VBG**, **KUL_FWT** (via the `scilpy` conda env, see "Conda environments" above), Karawun, and `xvfb-run` (**required** for headless `mrview` screenshots — see below). `KUL_nii2dcm.py` additionally needs python3 with SimpleITK, Pillow and numpy.

### Headless `mrview` screenshots (`xvfb-run`)

Every figure/PACS screenshot in this pipeline is captured by running `mrview` off-screen via `xvfb-run`. `xvfb-run` must be installed and on the `PATH` (`sudo apt install -y xvfb libgl1-mesa-dri`); the script checks for it at the point screenshots are generated and prints a warning if it is missing, but will otherwise fail when it tries to render figures. On systems where MATLAB/MCR (used by SPM12) is also installed, the `mrview` calls are additionally run through a filtered environment that strips the MCR's bundled Qt5 from `LD_LIBRARY_PATH` and forces software rendering (`LIBGL_ALWAYS_SOFTWARE=1`), otherwise `mrview` can abort with `Could not find the Qt platform plugin xcb/offscreen`. Override the `mrview` binary with `MRVIEW_BIN=/path/to/mrview` if `PATH` resolution is ambiguous.
