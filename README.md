# KULeuven Neuro Imaging Suite (KUL_NIS)

KUL_NIS provides tools to analyse (resting/task-based) fMRI, diffusion & structural MRI data.
Pipelines use many open-source packages.


It mainly provides automated pipelines to

1/ convert dicom files to the BIDS format
        - for Siemens, GE and Philips dicom data
        - for Philips data it includes basic computation of slice-timing and EES/ERT in BIDS format

2/ perform mriqc on structural and functional data

3/ perform fmriprep on structural and functional data

4/ perform freesurfer on the structural data organised in the BIDS format

5/ perform fastsurfer as an alternative to freesurfer (again in BIDS format)

6/ perform MRtix3 processing on dMRI data with our own developed 'dwiprep' pipeline (using BIDS input)

7/ perform qsiprep processing on dMRI data

8/ perform Virtual Brain Grafting to run freesurfer/fastsurfer on brains with large lesions, such as tumours or CVA lesions
        see: https://github.com/KUL-Radneuron/KUL_VBG

9/ perform MTC (magnetisation transfer contrast or ratio), T1w/T2w and T1/FLAIR ratio computation
        - using calibration
        - additionally compute MS lesion load using samseg

10/ perform automated fiber tractography using a combination of freesurfer/fastsurfer, KUL_VBG, and ifod2 MRtrix

11/ run AI tools such as HD-BET and HD-GLIO for automated brain extraction and tumour segmentation

12/ provide an automated pipeline to analyse clincal fMRI/dMRI data using:
        - GLM based fMRI SPM analysis
        - Melodic based ICA analysis
        - brain extraction using HD-BET
        - tumour delinaetion using HD-GLIO
        - brain parcellation using freesurfer/fastsurfer
        - automated tractography of all major tracts (own pipeline)



Depencies:
        Runs on linux or Osx/MacOs (the latter with some limitations).
        Needs an NVidia GPU (cuda) for some tools.
                - see https://github.com/treanus/KUL_Linux_Installation



Project PI's: Stefan Sunaert
Contributors: A. Radwan

Requires Mrtrix3, FSL, ants, dcm2niix, docker, fmriprep and mriqc

@ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be

### Changes made by AR:
1- Added KUL_radsyndisco.sh to allow a local run of Synb0DISCO from the fork -> https://github.com/Rad-dude/Synb0-DISCO.git
2- Added UKBB_fMRI_mod templates (modified from UKBB_fMRI template by warping to MNI space) under /atlasses/Temp_4_KUL_dwiprep
3- Modified KUL_dwifslpreproc to use local version of synb0disco via KUL_radsyndisco.sh
