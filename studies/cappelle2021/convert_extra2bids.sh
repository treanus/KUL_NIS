#!/bin/bash
# Sarah Cappelle & Stefan Sunaert
# 15/12/2020
# This script is the first part of Sarah's Study1
# We have data from IDE, but this is cluttered
# We sort the data into named DICOM folders
# We delete unwanted derived data (MPRs mostly)
# We convert these dicoms to BIDS format

# STEP 1 - sort the dicoms
# we use dicomsort.py from https://github.com/pieper/dicomsort
# this will warn for some double files, probably beacuse IDE export both dicom classic and enhanced
# The DICOMs from Ide are in /DATA/Data_MRI_Ide/DATA_RAW
#dicomsort.py -k /Users/xm52195/data/Cappelle_MS_DATA_ide/DICOM_extra \
# DICOM_extra_sorted/%PatientName/%StudyDate/%SeriesDescription-%SeriesNumber/%InstanceNumber.dcm

# STEP 2 - clean the DICOM_sorted
# IDE also exported derived images, such as MPRs from the 3D T1w
# these can be found in the series ending with a number higher than 1
#  e.g. series 301 is the 3D T1w
#       series 302 is the COR mpr of 301
#rm -rf DICOM_extra_sorted/*/*/*[2-9]

# STEP 3 - convert to bids
# We loop over all subjects
#  but only P001, P002, etc... (Note there are also non-anonymised data!)
# then we invoke KUL_dcm2bids 
#cd DICOM_sorted
#files=($( find -type d -maxdepth 2 -printf '%P\n' | grep P | grep -w '\w\{1,4\}'))
#cd ..
file="lijst.txt"
name=($(cut -d ',' -f1 $file ))
pat=($(cut -d ',' -f2 $file ))
n=${#name[@]}
n=$(($n-1))

for i in $(seq 0 $n); do
    p=${pat[$i]}
    pn=${name[$i]}
    f="DICOM_extra_sorted/${name[$i]}*"
    ses=($( find $f -type d -maxdepth 1))
    for d in "${ses[@]:1}"; do
        s=$(echo $d | cut -d '/' -f3)
        echo "Converting patient $pn participant $p session $s in directory $d"
        KUL_dcm2bids_new.sh -p $p -s $s -d $d -c study_config/sequences.txt
    done    
done

#rm -rf BIDS/tmp_dcm2bids
