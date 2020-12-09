#!/bin/bash
# Bash shell script to prepare fMRI/DTI results for iPlan Elements Server
#
# Requires Mrtrix3, Karawun
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/12/2020

t="_flipped"




# Tracts should be in folder '16bit'
mkdir -p for_elements/dti

#ls 16bit/*rBO${t}.hdr
mrconvert 16bit/*rBO${t}.img for_elements/dti/B0.nii

function kul_convert_tracts {
    file1=16bit/*r_$1_L${t}.img
    if [ -f $file1 ]; then
        mrcalc $file1 1 -gt $2 -mult for_elements/dti/$1_L.nii
    fi
    file1=16bit/*r_$1_R${t}.img
    if [ -f $file1 ]; then
        mrcalc $file1 1 -gt $3 -mult for_elements/dti/$1_R.nii
    fi
}

kul_convert_tracts "Cing" 1 2
kul_convert_tracts "CST" 3 4
kul_convert_tracts "FA" 5 6
kul_convert_tracts "FAT" 7 8
kul_convert_tracts "ILFIFOF" 9 10
kul_convert_tracts "SMACST" 11 12
kul_convert_tracts "UNC" 13 14

# fMRI Act should be in folder '16bit'
mkdir -p for_elements/fmri
mrconvert 16bit/Anat/*anat${t}.img for_elements/fmri/T1.nii

function kul_convert_fmri {
    file1=16bit/*_act_$1${t}.img
    echo $file1
    ls $file1
    if [ -f $file1 ]; then
        mrcalc $file1 2 -gt $2 -mult for_elements/fmri/$1_L.nii
    fi
}

kul_convert_fmri "Hand" 15
kul_convert_fmri "Lip" 16
kul_convert_fmri "Voet" 17
kul_convert_fmri "Taal" 18

