#!/bin/bash
# Bash shell script to prepare fMRI/DTI results for Brainlab Elements Server
#
# Requires Mrtrix3, Karawun
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/12/2020
version="0.1"

participant=sub-Marinissen

mkdir -p Karawun/${participant}/labels
mkdir -p Karawun/${participant}/tck
mkdir -p Karawun/${participant}/DICOM

mrcalc RESULTS/${participant}/Anat/T1w.nii 100 -div Karawun/${participant}/T1w.nii.gz -force

mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_VOIs/CSHP_LT_VOIs/CSHP_LT_incs1/CSHP_LT_incs1_map.nii.gz \
    regrid -template Karawun/${participant}/T1w.nii.gz \
    - | mrcalc - 0.1 -gt 3 -mul \
    Karawun/${participant}/labels/DISTAL_STN_MOTOR_LT.nii.gz -force
mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_VOIs/CSHP_RT_VOIs/CSHP_RT_incs1/CSHP_RT_incs1_map.nii.gz \
    regrid -template Karawun/${participant}/T1w.nii.gz \
    - | mrcalc - 0.1 -gt 3 -mul \
    Karawun/${participant}/labels/DISTAL_STN_MOTOR_RT.nii.gz -force
cp BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CSHP_LT_output/CSHP_LT_fin_BT_iFOD2.tck \
    Karawun/${participant}/tck/CSHDP_LT.tck
mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CSHP_LT_output/CSHP_LT_fin_map_BT_iFOD2.nii.gz \
    regrid -template Karawun/${participant}/T1w.nii.gz \
    - | mrcalc - 20 -gt 2 -mul \
    Karawun/${participant}/labels/CSHDP_LT_center.nii.gz -force
cp BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CSHP_RT_output/CSHP_RT_fin_BT_iFOD2.tck \
    Karawun/${participant}/tck/CSHDP_RT.tck
mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CSHP_RT_output/CSHP_RT_fin_map_BT_iFOD2.nii.gz \
    regrid -template Karawun/${participant}/T1w.nii.gz \
    - | mrcalc - 20 -gt 2 -mul \
    Karawun/${participant}/labels/CSHDP_RT_center.nii.gz -force
cp BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CST_RT_output/CST_RT_fin_BT_iFOD2.tck \
    Karawun/${participant}/tck/CST_RT.tck
cp BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CST_LT_output/CST_LT_fin_BT_iFOD2.tck \
    Karawun/${participant}/tck/CST_LT.tck
mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CST_LT_output/CST_LT_fin_map_BT_iFOD2.nii.gz \
    regrid -template Karawun/${participant}/T1w.nii.gz \
    - | mrcalc - 20 -gt 1 -mul \
    Karawun/${participant}/labels/CST_LT_center.nii.gz -force
mrgrid BIDS/derivatives/KUL_compute/${participant}/FWT/${participant}_TCKs_output/CST_RT_output/CST_RT_fin_map_BT_iFOD2.nii.gz \
    regrid -template Karawun/${participant}/T1w.nii.gz \
    - | mrcalc - 20 -gt 1 -mul \
    Karawun/${participant}/labels/CST_RT_center.nii.gz -force


conda activate KarawunEnv

importTractography -d DICOM/IM-0001-0001.dcm -o ${participant}_for_elements -n T1w.nii.gz -t tck/*.tck -l labels/*.gz