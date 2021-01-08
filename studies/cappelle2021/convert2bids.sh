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
dicomsort.py -k /DATA/Data_MRI_Ide/DATA_RAW \
 DICOM_sorted/%PatientName/%StudyDate/%SeriesDescription-%SeriesNumber/%InstanceNumber.dcm

# STEP 2 - clean the DICOM_sorted
# IDE also exported derived images, such as MPRs from the 3D T1w
# these can be found in the series ending with a number higher than 1
#  e.g. series 301 is the 3D T1w
#       series 302 is the COR mpr of 301
rm -rf DICOM_sorted/P*/*/*[2-9]

# STEP 3 - convert to bids
# We loop over all subjects
#  but only P001, P002, etc... (Note there are also non-anonymised data!)
# then we invoke KUL_dcm2bids 
cd DICOM_sorted
files=($( find -type d -maxdepth 2 -printf '%P\n' | grep P | grep -w '\w\{1,4\}'))
cd ..
for f
 in "${files[@]}"; do
    #echo $f
    s=${f#*/}
    p=${f%/*}
    #echo "p: $p"
    #echo "s: $s"
    if ! [ "$p" = "$s" ]; then
        echo "Converting participant $p session $s"
        KUL_dcm2bids_linux.sh -p $p -s $s -d DICOM_sorted/$f -c study_config/sequences.txt
    fi
done
#rm -rf BIDS/tmp_dcm2bids
