#!/bin/bash -e
#
# Creates a summary of information of the BIDS directory
# Information gathered is:
#  - subjects
#  - sessions
#  - available data (T1w, T2w, FLAIR, func, dwi)
#
# Requires Mrtrix3
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 16/02/2019 - alpha version
v="v0.1 - dd 16/02/2019"

# Source KUL_main_functions
# KUL_main_functions will:
#  - Say Welcome
#  - check wether all necessary software is installed (and exit if needed)
#  - provide some general functions like logging
kul_main_dir=$(dirname "$0")
#source $kul_main_dir/KUL_main_functions.sh
script=$(basename "$0")
cwd=$(pwd)


#---------- MAIN -------------------------------------------------
bids_dir=BIDS
output=BIDS_info.tsv

# find all images in the bids directory
search_mri=($(find ${bids_dir} -type f | grep nii.gz))
num_mri=${#search_mri[@]}

echo "Number of nifti data in the BIDS folder: $num_mri"
#echo ${search_mri[@]}

echo -e "MRI-scan \t Subject \t Session \t Type \t Scan \t Site \t Manufacturer \t Model \t TE \t TR \t DIM \t Dim_x \t Dim_y \t Dim_z \t Dynamics" > $output 

for i in `seq 0 $(($num_mri-1))`; do
    
    mri=${search_mri[$i]}
    echo "MRI-scan: $mri"

    sub=$(echo $mri | cut -d/ -f2)
    echo "Subject: $sub"

    ses=$(echo $mri | cut -d/ -f3)
    echo "Session: $ses"

    type=$(echo $mri | cut -d/ -f4)
    echo "Type: $type"

    scan=$(echo $mri | cut -d_ -f 3 | cut -d. -f 1)
    Echo "Scan: $scan"

    json=${mri%%.*}.json
    #echo $json

    site=$(grep StationName  $json | cut -d: -f2 | cut -d, -f 1)
    echo "Site: $site"

    manufacturer=$(grep \"Manufacturer\"  $json | cut -d: -f2 | cut -d, -f 1)
    echo "Manufacturer: $manufacturer"

    model=$(grep \"ManufacturersModelName\"  $json | cut -d: -f2 | cut -d, -f 1)
    echo "Model: $model"

    TE=$(grep EchoTime  $json | cut -d: -f2 | cut -d, -f 1)
    echo "TE: $TE"

    TR=$(grep RepetitionTime  $json | cut -d: -f2 | cut -d, -f 1)
    echo "TR: $TR"

    dim=$(mrinfo $mri -ndim)
    echo "Dimensions: $dim"

    dim_x=$(mrinfo $mri -size | cut -d" " -f 1)
    dim_y=$(mrinfo $mri -size | cut -d" " -f 2)
    dim_z=$(mrinfo $mri -size | cut -d" " -f 3)
    dynamics=$(mrinfo $mri -size | cut -d" " -f 4)

    echo -e "$mri \t $sub \t $ses \t $type \t $scan \t $site \t $manufacturer \
        \t $model \t $TE \t $TR \t $dim \t $dim_x \t $dim_y \t $dim_z \t $dynamics" >> $output 

done

