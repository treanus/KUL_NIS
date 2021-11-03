#main_folder=/Users/xm52195/data/Laura
#subject=STR_008
dcm_folder=DICOM
T1_ref=cT1_${1}.nii

#cd $main_folder/${subject}

mrresize -force $T1_ref T1_4mm.nii.gz -voxel 4 

bet $T1_ref T1_mask  -f 0.1 -R -m 

# convert to nii
#mrconvert -strides $T1_ref $dcm_folder anat.nii.gz -force
echo "SELECT THE (unregistered) PWI"
mrconvert -strides $T1_ref $dcm_folder pwi.nii.gz -force
echo "SELECT THE K2 map"
mrconvert -strides $T1_ref $dcm_folder K2map.nii.gz -force

# make mean of registered pwi
fslmaths pwi.nii.gz -Tmean mean_pwi

# register mean pwi to betted T1w 
ants_anat=T1_mask.nii.gz
ant_pwi=mean_pwi.nii.gz
mkdir -p affine
ants_type=affine/affine

antsRegistration --verbose 1 --dimensionality 3 \
    --output [${ants_type}_out,mean_pwi_reg2T1.nii.gz,${ants_type}_outInverseWarped.nii.gz] \
    --interpolation Linear \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$ants_anat,$ant_pwi,1] \
    --transform Rigid[0.1] \
    --metric MI[$ants_anat,$ant_pwi,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox


antsApplyTransforms -d 3 -i K2map.nii.gz -r mean_pwi_reg2T1.nii.gz -o K2_reg2T1.nii.gz -t affine/affine_out0GenericAffine.mat 
antsApplyTransforms -e 3 -i pwi.nii.gz -r T1_4mm.nii.gz -o pwi_reg2T1.nii.gz -t affine/affine_out0GenericAffine.mat 

mrfilter K2_reg2T1.nii.gz smooth -fwhm 8 smK2_reg2T1.nii.gz -force

fslmaths smK2_reg2T1.nii.gz -mas T1_mask_mask smK2_reg2T1_masked
fslmaths K2_reg2T1.nii.gz -mas T1_mask_mask K2_reg2T1_masked

mrconvert pwi_reg2T1.nii.gz -coord 3 5:95 pwi_reg2T1_without_dummy.nii.gz -force

#cd $main_folder