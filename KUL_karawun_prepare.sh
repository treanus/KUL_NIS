#!/bin/bash
# Bash shell script to prepare fMRI/DTI results for Brainlab Elements Server
#
# Requires Mrtrix3, Karawun
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/12/2020
version="0.1"

participant=sub-Vandervelpen


function KUL_karawun_get_tract {
    cp BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/${tract_name_orig}_output/${tract_name_orig}_fin_BT_iFOD2.tck \
        Karawun/${participant}/tck/${tract_name_final}.tck
    mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/${tract_name_orig}_output/${tract_name_orig}_fin_map_BT_iFOD2.nii.gz \
        regrid -template Karawun/${participant}/T1w.nii.gz \
        - | mrcalc - ${tract_threshold} -gt ${tract_color} -mul \
        Karawun/${participant}/labels/${tract_name_final}_center.nii.gz -force
}

function KUL_karawun_get_voi {
    mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_VOIs/${tract_name_orig}_VOIs/${tract_name_orig}_incs1/${tract_name_orig}_incs1_map.nii.gz \
        regrid -template Karawun/${participant}/T1w.nii.gz \
        - | mrcalc - ${voi_threshold} -gt ${voi_color} -mul \
        Karawun/${participant}/labels/${voi_name_final}.nii.gz -force
}

mkdir -p Karawun/${participant}/labels
mkdir -p Karawun/${participant}/tck
mkdir -p Karawun/${participant}/DICOM

mrcalc RESULTS/${participant}/Anat/T1w.nii 100 -div Karawun/${participant}/T1w.nii.gz -force

tract_name_orig="CSHP_LT"
voi_name_final="DISTAL_STN_MOTOR_Left"
voi_color=3
voi_threshold=0.1
KUL_karawun_get_voi

tract_name_orig="CSHP_RT"
voi_name_final="DISTAL_STN_MOTOR_Right"
voi_color=3
voi_threshold=0.1
KUL_karawun_get_voi

tract_name_orig="CSHP_LT"
tract_name_final="CSHDP_Left"
tract_color=2
tract_threshold=20
KUL_karawun_get_tract

tract_name_orig="CSHP_RT"
tract_name_final="CSHDP_Right"
tract_color=2
tract_threshold=20
KUL_karawun_get_tract

tract_name_orig="CST_LT"
tract_name_final="CST_Left"
tract_color=1
tract_threshold=20
KUL_karawun_get_tract

tract_name_orig="CST_RT"
tract_name_final="CST_Right"
tract_color=1
tract_threshold=20
KUL_karawun_get_tract

exit


conda activate KarawunEnv

importTractography -d DICOM/IM-0001-0001.dcm -o ${participant}_for_elements -n T1w.nii.gz -t tck/*.tck -l labels/*.gz