#!/bin/bash
# Sarah Cappelle & Stefan Sunaert
# 22/12/2020
# This script is the first part of Sarah's Study1
# 
# This scripts follows the rationale of D. Pareto et al. AJNR 2020
# Starting from 3D-T1w, 3D-FLAIR and 2D-T2w scans we compute:
#  create masked brain images using HD-BET
#  ANTs rigid coregister and reslice all images to the 3D-T1w (in isotropic 1 mm space)
#  bias correct the images using N4biascorrect fron ANTs
#  compute a T1FLAIR_ratio and a T1T2_ratio

# here we give the data
T1w="BIDS/sub-P001/ses-20120707/anat/sub-P001_ses-20120707_T1w.nii.gz"
T2w="BIDS/sub-P001/ses-20120707/anat/sub-P001_ses-20120707_T2w.nii.gz"
FLAIR="BIDS/sub-P001/ses-20120707/anat/sub-P001_ses-20120707_FLAIR.nii.gz"


# cofigure and make the results fodler
images="$T1w $T2w $FLAIR"
out=("T1w" "T2w" "FLAIR")
mkdir -p T1T2FLAIR_ratio/compute

# Brain extract the T1w, T2w and FLAIR
i=0
# perform HD-BET
for nii in $images;do
    echo $nii
    fslreorient2std $nii T1T2FLAIR_ratio/compute/std_${out[i]}
    if [ ! "${out[i]}" == "T2w" ]; then
        mrgrid T1T2FLAIR_ratio/compute/std_${out[i]}.nii.gz crop -axis 0 24,24 -axis 2 48,0 \
            T1T2FLAIR_ratio/compute/std_cropped_${out[i]}.nii.gz
    else
        cp "T1T2FLAIR_ratio/compute/std_${out[i]}.nii.gz" "T1T2FLAIR_ratio/compute/std_cropped_${out[i]}.nii.gz"
    fi
    hd-bet -i T1T2FLAIR_ratio/compute/std_cropped_${out[i]}.nii.gz -o T1T2FLAIR_ratio/compute/${out[i]} #-device cpu -mode fast -tta 0
    i=$((i+1))
done


# Rigidly register the T2w to the T1w
function KUL_rigid_register {
antsRegistration --verbose 1 --dimensionality 3 \
        --output [T1T2FLAIR_ratio/compute/${ants_type},T1T2FLAIR_ratio/compute/${newname}.nii.gz] \
        --interpolation BSpline \
        --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [T1T2FLAIR_ratio/compute/$ants_template,T1T2FLAIR_ratio/compute/$ants_source,1] \
        --transform Rigid[0.1] \
        --metric MI[T1T2FLAIR_ratio/compute/$ants_template,T1T2FLAIR_ratio/compute/$ants_source,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
}

echo " registering the the T2w to the betted T1w image (rigid)..."
ants_type="rigid_T2w_reg2t1"
newname="T2w_reg2t1"
ants_template="T1w.nii.gz"
ants_source="T2w.nii.gz"
KUL_rigid_register

echo " registering the the FLAIR to the betted T1w image (rigid)..."
ants_type="rigid_FLAIR_reg2t1"
newname="FLAIR_reg2t1"
ants_template="T1w.nii.gz"
ants_source="FLAIR.nii.gz"
KUL_rigid_register

# Reslice all images into 1mm isotropic space
mrgrid T1T2FLAIR_ratio/compute/T1w.nii.gz regrid -voxel 1 T1T2FLAIR_ratio/compute/T1w_iso.nii.gz
mrgrid T1T2FLAIR_ratio/compute/T2w_reg2t1.nii.gz regrid -voxel 1 T1T2FLAIR_ratio/compute/T2w_iso.nii.gz
mrgrid T1T2FLAIR_ratio/compute/FLAIR_reg2t1.nii.gz regrid -voxel 1 T1T2FLAIR_ratio/compute/FLAIR_iso.nii.gz

# compute the T1/T2 and T1/FLAIR ratio
mrcalc T1T2FLAIR_ratio/compute/T1w_iso.nii.gz T1T2FLAIR_ratio/compute/T2w_iso.nii.gz -divide T1T2FLAIR_ratio/T1T2ratio.nii.gz
mrcalc T1T2FLAIR_ratio/compute/T1w_iso.nii.gz T1T2FLAIR_ratio/compute/FLAIR_iso.nii.gz -divide T1T2FLAIR_ratio/T1FLAIRratio.nii.gz

