#!/bin/bash
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
search_mri=($(find ${bids_dir} -type f | grep nii.gz | sort))
num_mri=${#search_mri[@]}

echo "Number of nifti data in the BIDS folder: $num_mri"
#echo ${search_mri[@]}

echo -e "MRI-scan, Subject, Session, Type, Scan, Site, Manufacturer, Model, Software, Coil, MagneticFieldStrength, SeriesDescription, SeriesNumber, AcquisitionType, TE, TR, DIM, Dim_x, Dim_y, Dim_z, Dynamics, ETL" > $output 

for i in `seq 0 $(($num_mri-1))`; do
    
    mri=${search_mri[$i]}
    echo "MRI-scan: $mri"

    sub=$(echo $mri | cut -d/ -f2)
    echo "Subject: $sub"

    ses=$(echo $mri | cut -d/ -f3)
    echo "Session: $ses"

    type=$(echo $mri | cut -d/ -f4)
    echo "Type: $type"

    scan=$(echo $mri | awk -F_ '{print $NF}' | cut -d. -f 1)
    echo "Scan: $scan"

    json=${mri%%.*}.json
    #echo $json

    site=$(grep StationName  $json | cut -d: -f2 | cut -d, -f 1)
    echo "Site: $site"

    manufacturer=$(grep \"Manufacturer\"  $json | cut -d: -f2 | cut -d, -f 1 | tr -d '"')
    echo "Manufacturer: $manufacturer"

    model=$(grep \"ManufacturersModelName\"  $json | cut -d: -f2 | cut -d, -f 1 | tr -d '"')
    echo "Model: $model"
    
    soft=$(grep \"SoftwareVersions\"  $json | cut -d: -f2 | cut -d, -f 1 | tr -d '"')
    echo "Software: $soft"

    coil=$(grep \"CoilString\"  $json | cut -d: -f2 | cut -d, -f 1 | tr -d '"')
    echo "Coil: $coil"

    MagneticFieldStrength=$(grep \"MagneticFieldStrength\"  $json | cut -d: -f2 | cut -d, -f 1)
    echo "MagneticFieldStrength: $MagneticFieldStrength"

    SeriesDescription=$(grep \"SeriesDescription\"  $json | cut -d: -f2 | cut -d, -f 1 | tr -d '"')
    echo "SeriesDescription: $SeriesDescription"

    SeriesNumber=$(grep \"SeriesNumber\"  $json | cut -d: -f2 | cut -d, -f 1)
    echo "SeriesNumber: $SeriesNumber"

    AcquisitionType=$(grep \"MRAcquisitionType\"  $json | cut -d: -f2 | cut -d, -f 1 | tr -d '"')
    echo "AcquisitionType: $AcquisitionType"

    TE=$(grep EchoTime  $json | cut -d: -f2 | cut -d, -f 1)
    echo "TE: $TE"

    TR=$(grep -w RepetitionTime  $json | cut -d: -f2 | cut -d, -f 1)
    echo "TR: $TR"

    dim=$(mrinfo $mri -ndim)
    echo "Dimensions: $dim"

    dim_x=$(mrinfo $mri -size | cut -d" " -f 1)
    dim_y=$(mrinfo $mri -size | cut -d" " -f 2)
    dim_z=$(mrinfo $mri -size | cut -d" " -f 3)
    dynamics=$(mrinfo $mri -size | cut -d" " -f 4)

    ETL=$(grep EchoTrainLength $json | cut -d: -f2 | cut -d, -f 1)
    echo "ETL: $ETL"

    echo -e "$mri, $sub, $ses, $type, $scan, $site, $manufacturer \
       , $model, $soft, $coil, $MagneticFieldStrength, $SeriesDescription, $SeriesNumber \
       , $AcquisitionType, $TE, $TR, $dim, $dim_x, $dim_y, $dim_z, $dynamics, $ETL" >> $output 

done

