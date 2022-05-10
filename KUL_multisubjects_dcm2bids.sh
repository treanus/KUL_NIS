#!/bin/bash
# Bash shell script wrapper for calling KUL_dcm2bids for multiple subjects
#
# Requires KUL_dcm2bids
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 26/10/2018 - alpha version
version="v0.1 - dd 26/10/2018"



# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl from KUL_main_functions (for logging)

# source general functions
kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd, mrtrix3new & $log_dir is made in main_functions

# BEGIN LOCAL FUNCTIONS --------------

# --- function Usage ---
function Usage {

cat <<USAGE

`basename $0` wrapper for calling KUL_dcm2bids for multiple subjects

Usage:

  `basename $0` -d dicom_dir -c config_file -o bids_dir <OPT_ARGS>

  Depends on a config file that defines parameters with subjects, dicom_zip_files, study_sequence_descriptions e.g.

    BIDS_participant;EAD;dicom_zip;config_file;session
    P001;70386396;John_Lennon.zip;study_sequences_2parts.txt;tp01
    P002;60423284;Curt_Cobain.zip;study_sequences_3parts.txt;tp01
    P003;77503712;Stefan_Sunaert.zip;study_sequences_3parts.txt;tp01

  explains that we scanned 3 subjects in this study. 
    BIDS_participant will be the anonymised participant name in the BIDS subfolder
    EAD  = an internal reference number
    dicom_zip_file = name of the zipfile containing all dicoms
    config_file = description of the subjects (see above)
    session = e.g. tp01 (timepoint 1) in a longitudinal study with multiple timepoints

Example:

  `basename $0` -d DICOM -c Study_config/Subjects_for_dcm2bids_conversion.csv -o BIDS

Required arguments:

     -d:  directory where all dicom_zip_file (the zip or tar.gz containing all your dicoms) are stored
     -c:  config_file
     -o:  bids directory

Optional arguments:

     -t:  temporary directory (default = /tmp)
     -e:  copy task-*_events.tsv from config to BIDS dir
     -v:  verbose

USAGE

    exit 1
}

# END LOCAL FUNCTIONS --------------




# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
silent=1
tmp=/tmp

# Set flags
dcmdir_flag=0
conf_flag=0
bids_flag=0
tmp_flag=0
events_flag=0

if [ "$#" -lt 3 ]; then
    Usage >&2
    exit 1

else

    while getopts "d:c:o:veth" OPT; do

        case $OPT in
        d) #dicom_dir
            dcmdir_flag=1
            dcmdir=$OPTARG
        ;;
        c) #config_file
            conf_flag=1
            conf=$OPTARG
        ;;
        o) #bids output directory
            bids_flag=1
            bids_output=$OPTARG
        ;;
        v) #verbose
            silent=0
        ;;
        e) #events for task based fmri
            events_flag=1
        ;;
        t) #temporary directory
            tmp_flag=1
            tmp=$OPTARG
        ;;
        h) #help
            Usage >&2
            exit 0
        ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo
            Usage >&2
            exit 1
        ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            echo
            Usage >&2
            exit 1
        ;;
        esac

    done

fi

# check for required options
if [ $dcmdir_flag -eq 0 ] ; then 
    echo 
    echo "Option -d is required: give the directory with your dicoms zip files" >&2
    echo
    exit 2 
fi 

if [ $conf_flag -eq 0 ] ; then 
    echo 
    echo "Option -c is required: give the path to the file that describes the subjects" >&2
    echo
    exit 2 

elif [ ! -f $conf ] ; then
    echo 
    echo "The config file $conf does not exist"
    echo
    exit 2
fi 

if [ $bids_flag -eq 0 ] ; then 
    echo 
    echo "Option -o is required: give the path to the BIDS directory" >&2
    echo
    exit 2 
fi 

# INITIATE ---



# ----------- SAY HELLO ----------------------------------------------------------------------------------

if [ $silent -eq 0 ]; then
    echo "  The script you are running has basename `basename "$0"`, located in dirname $kul_main_dir"
    echo "  The present working directory is `pwd`"
fi

echo " reading $conf"

# we read the config file (and it may be csv, tsv or ;-seperated)
while IFS=$'\t,;' read -r BIDS_participant EAD dicom_zip config_file session comment; do
    
    
    if [ "$dicom_zip" = "dicom_zip" ]; then
        
        echo "first line" > /dev/null 2>&1

    else
    
        if [ $silent -eq 0 ]; then

            echo " Information about the config file:"
            echo "  BIDS_participant = $BIDS_participant"
            echo "  EAD = $EAD"
            echo "  dicom_zip_file = $dicom_zip"
            echo "  sequences_config_file = $config_file"
            echo "  session = $session"
            echo "  bids_output = $bids_output"
        fi
        
        # check if already converted
        bids_dir_to_check=${bids_output}/sub-${BIDS_participant}/ses-${session}
        
        #echo $BIDS_participant
        #echo $bids_dir_to_check

        #if [ ! -d $bids_dir_to_check ]; then

            kul_e2cl "Performing KUL_dcm2bids.sh -d $dcmdir/$dicom_zip -p $BIDS_participant -c $config_file -o $bids_output -s "${session}" " $log
            KUL_dcm2bids.sh -d $dcmdir/$dicom_zip -p $BIDS_participant -c $config_file -o $bids_output -s $session
        
        #else
        
        #    echo " BIDS conversion of participant $BIDS_participant already done, skipping..."

        #fi

    fi

done < $conf

# copying task based events.tsv to BIDS directory
if [ $events_flag -eq 1 ]; then
    kul_e2cl " Copying task based events.tsv to BIDS directory" $log
    cp Study_config/task-*_events.tsv $bids_output &
fi

# make a full .cvs file with final dicom tag info of all subjects
kul_e2cl " Making a full .cvs file with final dicom tag info of all subjects" $log
csv_all_subjects=$log_dir/ALL_subjects_dicom_info.csv
csv_single_subject=${preproc}/log/KUL_dcm2bids.sh/\*_final_dicom_info.csv
#echo "cat $csv_single_subject > $csv_all_subjects"
cat $csv_single_subject > $csv_all_subjects &

kul_e2cl "Finished $script" $log
