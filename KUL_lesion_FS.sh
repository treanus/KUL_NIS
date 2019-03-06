#!/bin/bash -e
# @ Ahmed Radwan & Stefan Sunaert ahmed.radwan@kuleuven.be
#
# v 0.1 - dd 02/02/2019 - dev

v="0.1 - dd 02/02/2019"

# This script is meant for allowing a decent recon-all output in the presence of a large brain lesion
# It is not a final end-all solution but a rather crude and simplistic work around
# The main idea is to replace the lesion with a hole and fill the hole with information from normal hemisphere
# maintains subject specificity and diseased hemisphere information but replaces lesion tissue with sham brain
# should be followed by a loop calculating overlap between fake labels resulting from sham brain with actual lesion
#
# Description:
# We start after fmripep is finished (should also run the T2 through fmriprep)
# We need to always have a T1 and a T2
# Use both native T1&T2 and MNIT1&T2
# 		- use itk-snap for lesion masking on both the native and MNI space T2/T1 WIs
#		- binarize, smooth and invert the lesion mask from itk (giving lmask_bin and lmask_binv)
# 		- make RL flipped struct. images with FSL
#		- make a hole in orig. struct images with fslmaths
#		- derive sham lesioned brain from flipped images
#		- create sham filled struct images using fslmaths
#		- Warp back to native space
#		- same workflow for native space images
# 		- run recon-all with T1 + T2 for both sets of images as timepoints
#		- use fslcc with lobe specific masks to estimate overlap of lesion with lobe specific labels
#		- need to derive lobe specific labels from FS to determine which labels to calculate overlap and stats for
#
# mri_annotation2label --subject ${subj} --sd ${cwd}/freesurfer/sub-${subj} --hemi lh --lobesStrict lobes
# mri_annotation2label --subject ${subj} --sd ${cwd}/freesurfer/sub-${subj} --hemi rh --lobesStrict lobes
# mri_aparc2aseg --s ${subj} --sd ${cwd}/freesurfer/sub-${subj}  --labelwm --hypo-as-wm --rip-unknown \
#   --volmask --o ${cwd}/freesurfer/sub-${subj}/${subj}/wmparc.lobes.mgz --ctxseg aparc+aseg.mgz \
#   --annot lobes --base-offset 200