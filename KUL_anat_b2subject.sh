#!/bin/bash
# Bash shell script to take 1 image in mni space and bring it to native space using the inverse of the warp field from the normalization step

participant=$1
input=$2
output=$(basename $input)_subject_space.nii.gz

cwd=$(pwd)
transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
reference=${cwd}/RESULTS/sub-${participant}/Anat/T1w.nii.gz

antsApplyTransforms -d 3 -i $input -r $reference -t $transform -o $output -n Linear
