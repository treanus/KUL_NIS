# KULeuven Neuro Imaging Suite (KUL_NIS)

KUL_NIS provides tools:
- for conversion of dicom data to BIDS format of 
	- T1w (with and without Gd)
	- T2w 
	- FLAIR
	- DIR (double inversion recovery)
	- ASL (not fully BIDS compatible yet)
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
- working with KUL_VBG
        - a pipeline to Freesurfer and FastSurfer in patients with brain lesions (tumour/stroke)
- working with KUL_FWT
        - an automated csd probabilistic tractography pipeline
- for conversion of the output of pipelines back to DICOM for use in:
        - Brainlab Neurosurgery format using Karawun
        - covertion for a basic pacs system using a MevisLab interface

All scripts (only) work with BIDS data.

Note that some output data of the analysis pipelines can be converted back to dicom for import into a PACS system (see above).

The pipelines used should work fine with healthy volunteer data, but are being implemented for use with clinical data (tumors, stroke, MS, PD) but have not yet been tested fully with clinical data.

**Any use in a clinical environment is off-label, not FDA aproved, not CE-labeled or approved. Also see the license file please.**


## Tools for BIDS data conversion

### KUL_dcm2bids
Converts dicom data to BIDS format for multipule MRI vendors. 
Provides slice-timing, total-readout-time, phase encoding direction & other data for Philips scanners (which is automatically defined for Siemens and GE scanners).
 
### KUL_multisubjects_dcm2bids
KUL_dcm2bids to convert multiple datasets.

### KUL_bids_summary
Provides output of multiple parameters of a BIDS dataset, including acquisition date, scanner software verion, etc... readable in google sheets, excel, etc...
  
 
## Tools for fMRI analysis

### KUL_preproc_all
This script allows to start an fmriprep analysis with a config file.
See 

### KUL_fmriproc_spm

This script will analyse a standard block based active fMRI data with a paradigm using 30 seconds BASELIBE followed by 30 seconds TASK epochs.

### KUL_fmriproc_conn
 (under development)
 
## Tools for dMRI analysis  

### KUL_dwiprep

### KUL_dwiprep_anat

### KUL_dwiprep_MNI

## Tools for the study of myelin 

### Magnetisation Transfer Ratio

### T1/T2 and T1/FLAIR ratio as a myelin marker


## Tools for importing results back into dicom and transfer to PACS/BrainLab

### Karawun 

### MevisLab


## Other (under dev)

  
## Who are we
Dr. Ahmed Radwan - KUL - ahmed.radwan@kuleuven.be

Prof. Dr. Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be