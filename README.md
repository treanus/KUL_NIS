# KULeuven Neuro Imaging Suite (KUL_NIS)

KUL_NIS provides tools:
- for conversion of dicom data to BIDS format of 
	- T1w (with and without Gd)
	- T2w 
	- FLAIR
	- PDw
	- FGATIR
	- DIR (double inversion recovery, in development)
	- ASL (not fully BIDS compatible yet)
	- DSC perfusion (not fully BIDS compatible yet)
	- MTR
	- fMRI
	- dMRI
	- QSM (in development)

- for analysis of 
	- fMRI data
         - using fmriprep
         - resting stat analysis using FSL melodic
         - task based fMRI

	- diffusion MRI data
         - preprocessing using mrtrix3
         - coregistration to T1w data using synb0-disco and/or Ants
         - normalisation to MNI space using fmriprep
         - group fixel based analysis in group template space
	- myelin markers such as:
		 - MRT (magnetisation transfer ratio)
         - T1/T2 and T1/FLAIR ratio's
	- structural data
         - coregistration of data to the T1w without Gd
         - Freesufer parcellation
         - FastSurfer parcellation
		 - Tumor segmentation

- to work easily with KUL_VBG
	- a pipeline that enables to run Freesurfer and FastSurfer in patients with (large) brain lesions (tumour/stroke)

- to work easily with KUL_FWT
	- an automated csd probabilistic tractography pipeline

- for conversion of the output of pipelines back to DICOM for use in:
	- Brainlab Neurosurgery format using Karawun
	- conversion of results for import into a general pacs using a MevisLab interface

All scripts (only) work with BIDS data.

Note that some output data of the analysis pipelines can be converted back to dicom for import into a PACS system (see above).

The pipelines used should work fine with healthy volunteer data, but are being implemented for use with clinical data (tumors, stroke, MS, PD) but have not yet been tested fully with clinical data.

**Any use in a clinical environment is off-label, not FDA aproved, not CE-labeled or approved. Also see the license file please.**


## Requirements

KUL_NIS is a set of bash/python wrappers around established neuroimaging software. All scripts assume a Linux environment and expect the dependencies below to be installed and on the `PATH`. The easiest way to install most of them is [KUL_Linux_Installation](https://github.com/treanus/KUL_Linux_Installation).

### Core dependencies (used by most tools)

| Software | Version (tested) | Used for |
|---|---|---|
| [MRtrix3](https://www.mrtrix.org/) | **3.0.8-2097-g99963980** (`dev` branch, built with CMake+Ninja — see SOFTWARE_ROOT_SETUP.md) | dwiprep, mrconvert/mrview, tractography, fixel-based analysis, figures |
| [FSL](https://fsl.fmrib.ox.ac.uk/) | 6.x | topup/eddy distortion correction, melodic resting-state, fslmaths |
| [ANTs](https://github.com/ANTsX/ANTs) | latest | registration (rigid/SyN), N4 bias correction |
| [dcm2bids](https://github.com/UNFmontreal/Dcm2Bids) | **≥ 3.0** | dicom → BIDS (v2 supported only via `KUL_dcm2bids_v2bkup.sh`) |
| [dcm2niix](https://github.com/rordenlab/dcm2niix) | latest | dicom → NIfTI backend for dcm2bids |
| python3 | ≥ 3.8 | helper scripts (`KUL_nii2dcm.py`, `KUL_EDs_b2masks.py`, etc.) |
| `xvfb` / `xvfb-run` | **required** | headless `mrview` screenshots in `KUL_clinical_fmridti.sh` — every `mrview` call is wrapped in `xvfb-run`; without it, figure/PACS screenshot generation fails (the script checks and warns at runtime if `xvfb-run` is not on `PATH`) |
| [p7zip](https://www.7-zip.org/) (`7z`) | — | encrypted backup archives (`-B` option) |
| [dcmtk](https://dicom.offis.de/dcmtk) (`dcmsend`) | — | push DICOMs to PACS / Orthanc (`tools/send_2_orthanc.sh`) |

#### Python packages (for KUL_nii2dcm.py and other .py scripts)

| Package | Used for |
|---|---|
| [SimpleITK](https://simpleitk.org/) | reading donor DICOM metadata/geometry and NIfTI/TIFF input, and writing the output DICOM series in `KUL_nii2dcm.py` (no separate DICOM library needed) |
| [Pillow](https://python-pillow.org/) (PIL) | PNG loading, cropping and resizing in `KUL_nii2dcm.py`; version-agnostic (works with both pre- and post-9.1 Pillow resampling APIs). Also required by `KUL_EDs_b2masks.py --bg-image` snapshots (mrview quadrant crop and VTK render encoding) |
| [numpy](https://numpy.org/) | array math throughout |
| [nibabel](https://nipy.org/nibabel/) | NIfTI I/O in Python scripts |
| [scipy](https://scipy.org/) | distance transforms / morphology in `KUL_EDs_b2masks.py` |
| [matplotlib](https://matplotlib.org/) | **required** by `KUL_EDs_b2masks.py` — renders the default QC snapshot (Agg backend) whenever `--bg-image` is not given, and is the fallback if the `mrview`/VTK render fails |
| [vtk](https://vtk.org/) | optional; `KUL_EDs_b2masks.py --bg-image --volume-render` only (off-screen glass-brain snapshot, no display/`xvfb-run` needed) |

> See [KUL_nii2dcm](/docs/KUL_nii2dcm/KUL_nii2dcm.md) for the full option/dependency list, including `python3` (the script requires `python3` explicitly, not `python`).

### Pipeline / analysis dependencies (tool-specific)

| Software | Version (tested) | Used by |
|---|---|---|
| [fmriprep](https://fmriprep.org/) | **25.1.4** (all clinical configs) | `KUL_preproc_all`, clinical pipeline; run via Docker or Singularity |
| [mriqc](https://mriqc.readthedocs.io/) | latest | quality control; run via Docker or Singularity |
| [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) | **8.2.0** | recon-all, cortical parcellation (Lausanne2018, Glasser HCP-MMP1), subregion segmentation |
| [FastSurfer](https://github.com/Deep-MI/FastSurfer) | **v3+** (required for FreeSurfer 8 compatibility) | faster surface reconstruction alternative to recon-all |
| [SPM12](https://www.fil.ion.ucl.ac.uk/spm/) + MATLAB | SPM12 | task-fMRI GLM statistics (`KUL_fmriproc_spm`) |
| [synb0-disco](https://github.com/MASILab/Synb0-DISCO) | **v3.0** (`leonyichencai/synb0-disco:v3.0`) | susceptibility distortion correction when no reverse-PE acquisition exists |
| [HD-BET](https://github.com/MIC-DKFZ/HD-BET) | latest | brain extraction (`KUL_dwiprep -m 1`, `KUL_anat_register`) |
| [hd-glio-auto](https://github.com/NeuroAI-HD/HD-GLIO-AUTO) | latest | AI glioma segmentation (clinical types 1–3) |
| [resseg](https://github.com/fepegar/resseg) | latest | resection-cavity segmentation |
| [LoRE](https://github.com/TissueVisionMics/lore) (`lore_dwi2decomposition`, `lore_decomposition2contrast`) | latest | low-rank DWI decomposition / microstructure contrasts (`-D run_dwiprep_lore_sd.txt`) |
| [scilpy](https://github.com/scilus/scilpy) | **2.3.0** | tractography post-processing (`KUL_tracts_ocd`, `KUL_FWT`); installed by `setup_environment.sh` as a conda env named `scilpy` — `KUL_clinical_fmridti.sh` finds it automatically; override with `-f <env_name>` only if you used a different name |
| pyfMRI (nilearn/nibabel/numpy/scipy/pandas/matplotlib/pyyaml) | — | rsfMRI network mapping (`-N`) and the nilearn task-fMRI GLM engine (`-E nilearn`); installed by `setup_environment.sh`'s `env-pyfmri` section as a conda env named `pyfMRI` — found automatically; override with `-y <env_name>` only if you used a different name |
| [qsiprep](https://qsiprep.readthedocs.io/) | latest | alternative dMRI preprocessing (`KUL_qsiprep`) |

### Sibling KUL repositories

These are separate repos that the clinical pipeline calls and must be installed alongside KUL_NIS:

- [**KUL_VBG**](https://github.com/KUL-Radneuron/KUL_VBG) — Virtual Brain Grafting: enables FreeSurfer/FastSurfer in patients with large lesions. Brain extraction uses `mri_synthstrip` (FreeSurfer built-in; `-B 1` in the clinical pipeline).
- [**KUL_FWT**](https://github.com/KUL-Radneuron/KUL_FWT) — automated CSD probabilistic tractography pipeline. `KUL_clinical_fmridti.sh` no longer auto-prepends a `../KUL_FWT` sibling folder to `PATH`; make sure `KUL_FWT_make_VOIs.sh` / `KUL_FWT_make_TCKs.sh` are already resolvable on `PATH` before running the clinical pipeline. Tractography post-processing also requires a scilpy conda environment — `setup_environment.sh` creates it as `scilpy` and `KUL_clinical_fmridti.sh` finds it automatically; override with `-f <env_name>` only if you used a different name.
- [**KUL_DTI_ALPS**](KUL_DTI_ALPS/) — DTI-ALPS index calculation using MNI-space ROIs (bundled as a subdirectory of KUL_NIS_unified).
- [Karawun](https://github.com/DevelopmentalImagingMCRI/karawun) — convert tractography/segmentation results to Brainlab Neurosurgery format.

> Versions marked "latest" are not pinned by the code and track current releases. Versions in **bold** are explicitly set in study configs or the code itself and represent the values KUL_NIS is currently validated against.


## User guide

### [Tumour work-up with perfusion, Brainlab and PACS](/docs/KUL_tumour_workup/KUL_tumour_workup.md)

Start-to-finish walkthrough of a glioma case: DICOM in, perfusion maps, tract
labels and PACS series out. Covers the DSC perfusion step, the lesion as a
Brainlab label, the FAT1w registration-QA volume, reading the perfusion numbers,
and troubleshooting.


## Clinical pipeline (one command, dicom → figures/PACS)

### [KUL_clinical_fmridti](/docs/KUL_clinical_fmridti/KUL_clinical_fmridti.md)
The clinical batch pipeline. From a single DICOM input it runs the full presurgical / DBS fMRI–dMRI work-up — dcm2bids, tumor segmentation, KUL_VBG, fmriprep, SPM/melodic activation maps, KUL_dwiprep and KUL_FWT tractography — and produces review figures and (optionally, via `-R`) DICOMs for PACS and Brainlab. Seven processing types cover intra-/extra-axial glioma, manual-mask lesions, non-glioma cases, DBS (DRT / CSHD) and DTI-ALPS. Click the link in the header for the full guide.


## Tools for BIDS data conversion

### [KUL_dcm2bids](/docs/KUL_dcm2bids/KUL_dcm2bids.md)
Converts dicom data to BIDS format for multipule MRI vendors. 
Provides slice-timing, total-readout-time, phase encoding direction & other data for Philips scanners (which is automatically defined for Siemens and GE scanners).
Click the link in the header above for more information.
 
### KUL_multisubjects_dcm2bids
KUL_dcm2bids to convert multiple datasets at once.

### KUL_bids_summary
Provides output of multiple parameters of a BIDS dataset, including acquisition date, scanner software verion, voxel spacing (x/y/z), etc... written to `BIDS_info.tsv`, readable in google sheets, excel, etc...



## Tools for structural analysis

### KUL_anat_biascorrect
This script performs Ants N4 bias correction on the structural images of a participant in the BIDS folder.
Run KUL_anat_biascorrect.sh to see information.

### KUL_anat_register_rigid
This script rigidly registers structural images in the BIDS folder to the T1w (without contrast) in the BIDS folder.

### [KUL_anat_segment_tumor](/docs/KUL_anat_segment_tumor/KUL_anat_segment_tumor.md)
This script segments pre- or post-operative brain tumor lesions and/or resection cavities using a combination of AI tools hd-glio-auto, resseg and fastsurfer.
 


## Tools for fMRI analysis

### KUL_preproc_all
This script allows to start an fmriprep analysis with a config file.
See 

### KUL_fmriproc_spm

This script will run an automated analysis using SPM12 of a standard blocked design active fMRI data with a paradigm using 30 seconds BASELINE followed by 30 seconds TASK epochs, after these have been preprocessed with fmriprep.

### KUL_fmriproc_conn

This script will run an automated analysis using FSL melodic on active and resting-state fMRI data, after these have been preprocessed with fmriprep.
 
## Tools for perfusion analysis

### [KUL_dsc_perfusion](/docs/KUL_dsc_perfusion/KUL_dsc_perfusion.md)

Processes DSC (dynamic susceptibility contrast) perfusion data of tumour
patients: denoising, motion and EPI distortion correction, leakage-corrected
rCBV/rCBF/MTT/TTP/TT0 estimation, and normalisation against contralesional
normal-appearing white matter. All maps are delivered in the participant's T1w
space, so they can be sent to PACS and Karawun alongside the fMRI and
tractography results. Runs automatically from `KUL_clinical_fmridti.sh`
whenever a DSC series is present, or standalone via `KUL_dsc_perfusion.sh`.


## Tools for dMRI analysis  

### KUL_dwiprep
Preprocesses diffusion MRI data with MRtrix3: denoising, Gibbs unringing, motion/distortion correction via `dwifslpreproc` (topup/eddy or synb0-disco when no reverse phase-encoding is available), bias correction and brain masking (including a `mri_synthstrip` option). Run `KUL_dwiprep.sh` for options. For old studies whose headers lack phase-encoding info, set `export KUL_dwiprep_custom_dwifslpreproc="..."` to pass explicit parameters.

### KUL_dwiprep_anat
Coregisters the preprocessed dMRI to the subject's T1w (without Gd) and brings anatomical/parcellation information into diffusion space.

### KUL_dwiprep_MNI
Normalises diffusion-space results to MNI space (using fmriprep transforms) for group analysis.

## Tools for mask comparison

### [KUL_EDs_b2masks](/docs/KUL_EDs_b2masks/KUL_EDs_b2masks.md)
Swiss-knife distance tool for binary NIfTI masks: minimum Euclidean, Hausdorff, 95th-percentile Hausdorff and ASSD metrics, one-to-one/one-to-many/all-pairs comparison modes, multi-class label support, parallel workers, and PNG/HTML/CSV reporting per pair.

## Tools for the study of myelin 

### Magnetisation Transfer Ratio and T1/T2 and T1/FLAIR ratio as a myelin marker 

[KUL_T1T2FLAIRMTR_ratio](/docs/KUL_T1T2FLAIRMTR_ratio/KUL_T1T2FLAIRMTR_ratio.md)


## Tools for importing results back into dicom and transfer to PACS/BrainLab

Pipeline results (tracts, activation maps) can be converted back to DICOM so they can be reviewed in a clinical PACS or imported into Brainlab for neuronavigation. In the clinical pipeline this is triggered with `KUL_clinical_fmridti.sh -R <underlay>`, **after** reviewing the figures.

### [KUL_nii2dcm](/docs/KUL_nii2dcm/KUL_nii2dcm.md)
Wraps rendered PNG screenshots into a DICOM series, using a donor DICOM from the same study/session so the result links correctly in PACS and supports multi-planar reconstruction (`KUL_nii2dcm.py`). Click the link in the header for options and dependencies.

### [KUL_karawun_prepare](/docs/KUL_karawun_prepare/KUL_karawun_prepare.md)
Converts tractography/segmentation results into the Brainlab Neurosurgery format (`KUL_karawun_prepare.sh`, `KUL_karawun2brainlab.sh`), see [Karawun](https://github.com/DevelopmentalImagingMCRI/karawun). Writes the tract bundles, the FAT1w registration-QA volume, the lesion (tumour cases) and the thalamic VIM / STN target VOIs (DBS cases).

**Read the doc before changing any label colour.** A Karawun label's voxel value *is* its Brainlab colour, so two labels sharing a value are indistinguishable in the scene. The values are allocated in fixed non-overlapping ranges — tracts 1-41, VIM 16/30, STN 23/24, auto-assigned tracts 42-49, lesion 50, fMRI 51-63 — and stock Karawun clamps everything above index 30 to a single colour, so the KU Leuven extended-palette fork must be pinned for `KarawunDev`.

### send_2_orthanc
Pushes a directory of generated DICOMs to an Orthanc/PACS node with `dcmsend` (`tools/send_2_orthanc.sh`).

### MevisLab
A MevisLab interface is available to convert results for import into a general PACS.


## Other (under dev)

  
## Who are we
Dr. Ahmed Radwan - KUL - ahmed.radwan@kuleuven.be

Prof. Dr. Stefan Sunaert - UZ/KUL - stefan.sunaert@kuleuven.be
