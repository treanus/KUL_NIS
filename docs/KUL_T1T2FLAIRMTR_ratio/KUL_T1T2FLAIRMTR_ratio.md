# KUL_T1T2FLAIRMTR_ratio

## Purpose

Computes a T1w/T2w, T1w/FLAIR and MTC (magnetisation transfer contrast) ratio, using BIDS organised data.
The full methodology is described in "T1w/FLAIR ratio standardization as a myelin marker in MS patients,
        by S Cappelle, D Pareto, S Sunaert, I Smets, A Laenen, B Dubois, Ph Demaerel" - [Article link](https://pubmed.ncbi.nlm.nih.gov/36451354/)

## Required structural input data

Required are data organised in BIDS [format](../KUL_dcm2bids/KUL_dcm2bids.md):
- a T1w image without contrast
- a FLAIR and/or a T2w image
- a MTI image pair (optional)

## Usage

A typical command for running is:  

`KUL_T1T2FLAIRMTR_ratio.sh -a -n 10 -f 1 -m`

You need to be in **main** directory to run this command. This is also where the BIDS folder resides. 

For info about the command just run:

`KUL_T1T2FLAIRMTR_ratio.sh`


## Output in the folder T1T2FLAIRMTR_ratio

Images starting with sub-xxx-ses-yyy_:
- _T1w.nii.gz: the input T1w, converted into 1x1x1 mm, biascorrected
- _T2w_reg2T1w.nii.gz: the input T2w, converted into 1x1x1 mm, biascorrected and rigidly registered the the T1w
- _FLAIR_reg2T1w.nii.gz: the input FLAIR, converted into 1x1x1 mm, biascorrected and rigidly registered the the T1w

Further output is in 1x1x1 mm voxel dimensions.


**T1w/T2w, T1w/FLAIR ratio-maps** starting with sub-xxx-ses-yyy_:
- calib-none.nii.gz: ratio without calibration (raw T1w/T2w or T1w/FLAIR ratio)
- calib-lin.nii.gz: LINEAR histogram matching using eye/muscle tissue (Ganzetti et al. 2014 alike)
- calib-nonlin.nii.gz: NONLINEAR histogram matching using eye/muscle tissue (Ganzetti et al. 2014 alike)
- calb-nonlin2.nii.gz: NONLINEAR histogram matching according to Cappelle et al. 2022
- calb-nonlin3.nii.gz: NONLINEAR histogram matching on brain tissue, excluding white matter lesions (from samseg)

MTR ratio maps:
- ratio-MTC.nii.gz

ROIs:
- rois/sub-x_ses-x_MSLesion: a binary image containing the MS lesions as identified by Freesurfer Samseg

mask:
- in folder masks one can inspect the eye/muscle, brain, background, etc masks.

Note:
- outputs are also given warped to the ICBM 2009a Nonlinear Symmetric MNI space (mni_icbm152_t1_tal_nlin_sym_09a.nii)

## Further processing

Stats could be further generated with KUL_MS_lesion_stats.sh


## Depedencies

A cuda compatible GPU with at least 6GB memory.

Internally KUL_T1T2FLAIRMTR_ratio uses:
- [MRtix3](https://www.mrtrix.org/)
- [ANTs](http://stnava.github.io/ANTs/)
- [HD-BET](https://github.com/MIC-DKFZ/HD-BET)
- optionally [FastSurfer](https://deep-mi.org/research/fastsurfer/)
- optionally [Freesurfer Samseg](https://surfer.nmr.mgh.harvard.edu/fswiki/Samseg)

These could be installed using [KUL_Linux_Installation](https://github.com/treanus/KUL_Linux_Installation)
