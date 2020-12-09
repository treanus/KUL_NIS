#!/bin/bash
# Bash shell script to prepare fMRI/DTI results for iPlan Elements Server
#
# Requires Mrtrix3, Karawun
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/12/2020

#if 0; then
# Prepare the nii files for the tracts
mkdir -p for_elements/dti

# Convert the B0
mrconvert 16bit_Elements/*rBO.img for_elements/B0.nii

# convert each tract
cd 16bit_Elements
files=( $(find -name "*r_*.img" -type f -printf '%P\n') )


for t in "${files[@]}"; do 
#for t in CST_L CST_R SLF_L SLF_R IFOF_L IFOF_R ILF_L ILF_R SMAPT_L SMAPT_R FAT_L FAT_R CINGc_L CINGc_R CINGt_L CINGt_R; do
    
    echo $t
    PS3="File ${t} corresponds to tract: "

    select opt in Other Ignore CST_L CST_R SLF_L SLF_R IFOF_L IFOF_R \
            ILF_L ILF_R SMAPT_L SMAPT_R FAT_L FAT_R CINGc_L \
            CINGc_R CINGt_L CINGt_R UNC_L UNC_R; do 
        #echo "you have selected $REPLY"
        #echo "this is $opt"
        color=$REPLY
        if [ $REPLY -eq 1 ]; then
            read -p "Describe the tract: " tract
        elif [ $REPLY -eq 2 ]; then
            break
        else
            tract=$opt
        fi    
        mrcalc $t 1 -gt $color -mult ../for_elements/dti/${tract}.nii
        break
    done
done

cd ..
#fi

# Prepare the nii for the fMRI act
mkdir -p for_elements/fmri

# Convert the T1w anat
mrconvert 16bit_Anat/*anat.img for_elements/T1.nii

cd 16bit_Elements
files=( $(find -name "*_act_*.img" -type f -printf '%P\n') )

for fmri in "${files[@]}"; do 
    echo $fmri
    PS3="File ${fmri} corresponds to "

    select opt in Other Ignore HAND LIP FOOT LANGUAGE; do
        if [ $REPLY -eq 1 ]; then
            read -p "Describe the task: " task
        elif [ $REPLY -eq 2 ]; then
            break
        else
            task=$opt
        fi
        color=$REPLY    
        echo $task
        echo $color
        mrcalc $fmri 2 -gt $color -mult ../for_elements/fmri/${task}.nii
        break
    done
done
cd ..

# Call Karawun
echo "\n\n\n"
echo "Now preparing for Karawun DICOM conversion"

importTractography -d DICOM_classic/DICOM/IM_0002 -o for_elements/Upload_dti \
    -n for_elements/B0.nii -l for_elements/dti/*.nii


importTractography -d DICOM_classic/DICOM/IM_0002 -o for_elements/Upload_fmri \
    -n for_elements/T1.nii -l for_elements/fmri/*.nii
