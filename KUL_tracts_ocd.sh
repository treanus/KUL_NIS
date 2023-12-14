#!/bin/bash

function KUL_antsApply_Transform {

    echo "interpolation_type: $interpolation_type"
    antsApplyTransforms -d 3 --float 1 \
        --verbose 1 \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n $interpolation_type
}

function KUL_make_fs_roi {

    n_rois=${#r2m[@]}
    echo $n_rois
    n_roi=1
 
    for r2 in ${r2m[@]}; do

        echo $r2
        mrcalc $atlas_fs $r2 -eq $outdir/tmp/tmp${n_roi}.nii -force
        ((n_roi=n_roi+1))
        
    done
    
    mrmath $outdir/tmp/tmp*.nii sum $outdir/$roi_name -force
    rm -f $outdir/tmp/tmp*.nii

}


# general sttings
participant=StefanSunaert
thr_ATR_L=1.0
thr_ATR_R=1.0
thr_slMFM_L=1.0
thr_slMFB_R=1.0
ncpu=120
source=dwiprep/sub-${participant}/sub-${participant}/response/dhollander_wmfod_reg2T1w.mif 
select_slMFP=20000
outdir=tracto_ocd
mkdir -p $outdir/tmp/ABGT
interpolation_type="NearestNeighbor"
reference=fmriprep/sub-${participant}/anat/sub-${participant}_desc-preproc_T1w.nii.gz 
transform=fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 
atlas_CIT="/usr/local/KUL_apps/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/CIT168_Reinf_Learn (Pauli 2017) - imported from CIT168/mixed"
atlas_ABGT_orig="/usr/local/KUL_apps/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/Atlas of the Basal Ganglia and Thalamus (He 2020)"
cp -R "$atlas_ABGT_orig"/* $outdir/tmp/ABGT
atlas_ABGT=$outdir/tmp/ABGT
atlas_JHU="/usr/local/fsl/data/atlases/JHU"
atlas_fs=BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz
Hemi_L=BIDS/derivatives/KUL_compute/sub-${participant}/FWT/sub-${participant}_VOIs/custom_VOIs/Left_hemir_custom.nii.gz
Hemi_R=BIDS/derivatives/KUL_compute/sub-${participant}/FWT/sub-${participant}_VOIs/custom_VOIs/Right_hemir_custom.nii.gz

# slMFB

# make the VTA include
# atlas = CIT168
mrcalc "$atlas_CIT/VTA.nii.gz" 0.1 -gt $outdir/CIT168toMNI152-2009c_VTA.nii.gz
input=$outdir/CIT168toMNI152-2009c_VTA.nii.gz 
output=$outdir/CIT168toSubject-2009c_VTA.nii.gz
KUL_antsApply_Transform
# make both left and right
mrgrid $Hemi_L regrid -template $reference $outdir/Hemi_L.nii.gz
mrgrid $Hemi_R regrid -template $reference $outdir/Hemi_R.nii.gz
mrcalc $outdir/CIT168toSubject-2009c_VTA.nii.gz $outdir/Hemi_L.nii.gz -mul $outdir/CIT168toSubject-2009c_VTA_L.nii.gz
mrcalc $outdir/CIT168toSubject-2009c_VTA.nii.gz $outdir/Hemi_R.nii.gz -mul $outdir/CIT168toSubject-2009c_VTA_R.nii.gz
rm -f $outdir/CIT168toMNI152-2009c_VTA.nii.gz $outdir/CIT168toSubject-2009c_VTA.nii.gz


# make the ALIC includes
# atlas = JHU WM
# left
mrcalc $atlas_JHU/JHU-ICBM-labels-1mm.nii.gz 18 -eq - | maskfilter -npass 2 - dilate $outdir/JHU-ICMB-alic_L.nii.gz
input=$outdir/JHU-ICMB-alic_L.nii.gz
output=$outdir/JHU-subjet-alic_L.nii.gz
KUL_antsApply_Transform 
#right
mrcalc $atlas_JHU/JHU-ICBM-labels-1mm.nii.gz 17 -eq - | maskfilter -npass 2 - dilate $outdir/JHU-ICMB-alic_R.nii.gz
input=$outdir/JHU-ICMB-alic_R.nii.gz
output=$outdir/JHU-subjet-alic_R.nii.gz
KUL_antsApply_Transform 
rm -f $outdir/JHU-ICMB-alic_L.nii.gz $outdir/JHU-ICMB-alic_R.nii.gz


# make the prefrontal includes
# atlas = FS
# left
r2m=(1012 1014 1019 1032)
roi_name=FS-subject-PFC_L.nii.gz
KUL_make_fs_roi

#right
r2m=(2012 2014 2019 2032)
roi_name=FS-subject-PFC_R.nii.gz
KUL_make_fs_roi

# make BrainStem exclude
mrcalc $atlas_fs 16 -eq $outdir/FS-subject_exclude_BStmp.nii.gz
maskfilter -npass 2 $outdir/FS-subject_exclude_BStmp.nii.gz erode $outdir/FS-subject_exclude_BS.nii.gz
rm -f $outdir/FS-subject_exclude_BStmp.nii.gz

# make extra excludes
r2m=(1003 17 18 250 1004 1006 1024 3024 1016 3016 1035 3035)
roi_name=FS-subject_excludes4slMFB_L.nii.gz
KUL_make_fs_roi
r2m=(2003 53 54 250 2004 2006 2024 4024 2016 4016 2035 4035)
roi_name=FS-subject_excludes4slMFB_R.nii.gz
KUL_make_fs_roi


# Tckgen
# left
incl1=$outdir/CIT168toSubject-2009c_VTA_L.nii.gz
incl2=$outdir/JHU-subjet-alic_L.nii.gz
incl3=$outdir/FS-subject-PFC_L.nii.gz
tract=$outdir/slMFB_L.tck
tckgen -nthreads $ncpu $source -seed_image $incl1 -include $incl2 -include $incl3 \
    -exclude $outdir/FS-subject_exclude_BS.nii.gz -exclude $Hemi_R \
    -exclude $outdir/FS-subject_excludes4slMFB_L.nii.gz $tract \
    -select $select_slMFP
tcksift2 $tract $source $outdir/slMFB_L_weights.txt -force



# right
incl1=$outdir/CIT168toSubject-2009c_VTA_R.nii.gz
incl2=$outdir/JHU-subjet-alic_R.nii.gz
incl3=$outdir/FS-subject-PFC_R.nii.gz
tract=$outdir/slMFB_R.tck
tckgen -nthreads $ncpu $source -seed_image $incl1 -include $incl2 -include $incl3 \
    -exclude $outdir/FS-subject_exclude_BS.nii.gz -exclude $Hemi_L \
    -exclude $outdir/FS-subject_excludes4slMFB_R.nii.gz $tract \
    -select $select_slMFP
tcksift2 $tract $source $outdir/slMFB_R_weights.txt -force

# make maps
tract=$outdir/slMFB_L.tck
tract_center=$outdir/slMFB_center_L.tck
num_tck=$(tckstats -quiet -output count $outdir/slMFB_L.tck)
echo "slMFB_L has $num_tck streamlines"
tract_nii_full=$outdir/slMFB_L_full.nii.gz
tract_nii_thr=$outdir/slMFB_L_center.nii.gz
tract_weight=$outdir/slMFB_L_weights.txt 
tract_color=1
threshold=1
tract_corr_threshold=1
tract_threshold=$(( num_tck*threshold/tract_corr_threshold/10 ))
tract_threshold=$thr_slMFM_L
echo "The compute threshold is: $tract_threshold (with a correction of $tract_corr_threshold)"
tckmap -template $reference -tck_weights_in $tract_weight $tract $tract_nii_full -force
mrcalc $tract_nii_full $tract_threshold -gt ${tract_color} -mul $tract_nii_thr -force
#rm -f $tract_nii_full
tckedit -tck_weights_in $tract_weight $tract $tract_center


tract=$outdir/slMFB_R.tck
tract_center=$outdir/slMFB_center_R.tck
num_tck=$(tckstats -quiet -output count $outdir/slMFB_R.tck)
echo "slMFB_R has $num_tck streamlines"
tract_nii_full=$outdir/slMFB_R_full.nii.gz
tract_nii_thr=$outdir/slMFB_R_center.nii.gz
tract_weight=$outdir/slMFB_R_weights.txt 
tract_color=2
threshold=1
tract_corr_threshold=1
tract_threshold=$(( num_tck*threshold/tract_corr_threshold/10 ))
tract_threshold=$thr_slMFB_R
echo "The compute threshold is: $tract_threshold (with a correction of $tract_corr_threshold)"
tckmap -template $reference -tck_weights_in $tract_weight $tract $tract_nii_full -force
mrcalc $tract_nii_full $tract_threshold -gt ${tract_color} -mul $tract_nii_thr -force
#rm -f $tract_nii_full
tckedit -tck_weights_in $tract_weight $tract $tract_center


#ATR
Tha_ant_L="${atlas_ABGT}/lh/Tha-a.nii.gz"
Tha_ld_L="${atlas_ABGT}/lh/Tha-ld.nii.gz"
Tha_ant_R="${atlas_ABGT}/rh/Tha-a.nii.gz"
Tha_ld_R="${atlas_ABGT}/rh/Tha-ld.nii.gz"

input="$Tha_ant_L"
output=$outdir/Tha_ant_L.nii.gz
KUL_antsApply_Transform
input=$Tha_ant_R
output=$outdir/Tha_ant_R.nii.gz
KUL_antsApply_Transform
input=$Tha_ld_L
output=$outdir/Tha_ld_L.nii.gz
KUL_antsApply_Transform
input=$Tha_ld_R
output=$outdir/Tha_ld_R.nii.gz
KUL_antsApply_Transform

# Thalamus - Tha-a and Th-ld
mrcalc $atlas_fs 10 -eq $outdir/Thalamus_full_L.nii.gz
mrcalc $atlas_fs 49 -eq $outdir/Thalamus_full_R.nii.gz
mrgrid $outdir/Thalamus_full_L.nii.gz regrid -template $reference - | \
    mrcalc - $outdir/Tha_ant_L.nii.gz -subtract $outdir/Tha_ld_L.nii.gz -subtract - | \
    maskfilter -npass 2 - clean $outdir/Thalamus_excl_L.nii.gz
mrgrid $outdir/Thalamus_full_R.nii.gz regrid -template $reference - | \
    mrcalc - $outdir/Tha_ant_R.nii.gz -subtract $outdir/Tha_ld_R.nii.gz -subtract - | \
    maskfilter -npass 2 - clean $outdir/Thalamus_excl_R.nii.gz


#extra excludes
r2m=(17 18 28 250 1002 1003 1004 1006 1016 1024 1035 3003 3024 3016 3035)
roi_name=FS-subject_excludes4ATR_L.nii.gz
KUL_make_fs_roi
r2m=(53 54 60 250 2002 2003 2004 2006 2016 2024 2035 4003 4024 4016 4035)
roi_name=FS-subject_excludes4ATR_R.nii.gz
KUL_make_fs_roi


# Tckgen
source=dwiprep/sub-${participant}/sub-${participant}/response/dhollander_wmfod_reg2T1w.mif 
# left
incl1=$outdir/Tha_ant_L.nii.gz
incl2=$outdir/JHU-subjet-alic_L.nii.gz
incl3=$outdir/FS-subject-PFC_L.nii.gz
tract=$outdir/ATR_L.tck
tract_weight=$outdir/ATR_L_weights.txt 
tckgen -nthreads $ncpu $source -seed_image $incl1 -include $incl2 -include $incl3 \
    -exclude $outdir/Thalamus_excl_L.nii.gz \
    -exclude $Hemi_R -exclude $outdir/FS-subject_excludes4ATR_L.nii.gz $tract
tcksift2 $tract $source $tract_weight -force

# right
incl1=$outdir/Tha_ant_R.nii.gz
incl2=$outdir/JHU-subjet-alic_R.nii.gz
incl3=$outdir/FS-subject-PFC_R.nii.gz
tract=$outdir/ATR_R.tck
tract_weight=$outdir/ATR_R_weights.txt 
tckgen -nthreads $ncpu $source -seed_image $incl1 -include $incl2 -include $incl3 \
    -exclude $outdir/Thalamus_excl_R.nii.gz \
    -exclude $Hemi_L -exclude $outdir/FS-subject_excludes4ATR_R.nii.gz $tract
tcksift2 $tract $source $tract_weight -force

# make maps
tract=$outdir/ATR_L.tck
num_tck=$(tckstats -quiet -output count $tract)
echo "$tract has $num_tck streamlines"
tract_nii_full=$outdir/ATR_L_full.nii.gz
tract_nii_thr=$outdir/ATR_L_center.nii.gz
tract_weight=$outdir/ATR_L_weights.txt 
tract_color=3
threshold=1
tract_corr_threshold=1
tract_threshold=$(( num_tck*threshold/tract_corr_threshold/10 ))
tract_threshold=$thr_ATR_L
echo "The compute threshold is: $tract_threshold (with a correction of $tract_corr_threshold)"
tckmap -template $reference -tck_weights_in $tract_weight $tract $tract_nii_full -force
mrcalc $tract_nii_full $tract_threshold -gt ${tract_color} -mul $tract_nii_thr -force
#rm -f $tract_nii_full

tract=$outdir/ATR_R.tck
num_tck=$(tckstats -quiet -output count $tract)
echo "$tract has $num_tck streamlines"
tract_nii_full=$outdir/ATR_R_full.nii.gz
tract_nii_thr=$outdir/ATR_R_center.nii.gz
tract_weight=$outdir/ATR_R_weights.txt 
tract_color=4
threshold=1
tract_corr_threshold=1
tract_threshold=$(( num_tck*threshold/tract_corr_threshold/10 ))
tract_threshold=$thr_ATR_R
echo "The compute threshold is: $tract_threshold (with a correction of $tract_corr_threshold)"
tckmap -template $reference -tck_weights_in $tract_weight $tract $tract_nii_full -force
mrcalc $tract_nii_full $tract_threshold -gt ${tract_color} -mul $tract_nii_thr -force
#rm -f $tract_nii_full