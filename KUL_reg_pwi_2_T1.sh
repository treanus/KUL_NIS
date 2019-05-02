dcm_folder=STR_008
T1_ref=cT1_STR_008.nii

mrresize cT1_STR_008.nii cT1_4mm.nii.gz -voxel 4

# convert to nii
mrconvert -strides $T1_ref $dcm_folder anat.nii.gz -force
mrconvert -strides $T1_ref $dcm_folder pwi.nii.gz -force
mrconvert -strides $T1_ref $dcm_folder K2map.nii.gz -force

# make mean of registered pwi
fslmaths pwi.nii.gz -Tmean MEAN_pwi

# register mean pwi to betted T1w 
ants_anat=$T1_ref
ant_pwi=MEAN_pwi.nii.gz
ants_type=rigid

if [ ! -f rigid_outWarped.nii.gz ]; then

    kul_e2cl " registering the the mean pwi to the betted T1w image (rigid)..." ${log}
    antsRegistration --verbose 1 --dimensionality 3 \
        --output [${ants_type}_out,${ants_type}_outWarped.nii.gz,${ants_type}_outInverseWarped.nii.gz] \
        --interpolation Linear \
        --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [$ants_anat,$ant_pwi,1] \
        --transform Rigid[0.1] \
        --metric MI[$ants_anat,$ant_pwi,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox

else

    echo " registering the T1w image (rigid) already done, skipping..."

fi

antsApplyTransforms -d 3 -i K2map.nii.gz -r ${ants_type}_outWarped.nii.gz -o ${ants_type}_K2.nii.gz -t rigid_out0GenericAffine.mat 
antsApplyTransforms -e 3 -i pwi.nii.gz -r cT1_4mm.nii.gz -o ${ants_type}_pwi.nii.gz -t rigid_out0GenericAffine.mat 