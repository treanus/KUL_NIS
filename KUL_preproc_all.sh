#!/bin/bash -e
# Bash shell script to preprocces an entire study
#
# Requires dcm2bids (jooh fork), dcm2niix, Mrtrix3
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 06/11/2018 - alpha version
v="v0.1 - dd 06/11/2018"

# TODO
#  - make it work for multiple vendors
#  - wrap around for multiple subjects

ncpu=6
mem_mb=15000

# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl from KUL_main_functions (for logging)
#  - kul_dcmtags (for reading specific parameters from dicom header)

# source general functions
kul_main_dir=`dirname "$0"`
script=`basename "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# BEGIN LOCAL FUNCTIONS --------------

# --- function Usage ---
function Usage {

cat <<USAGE

`basename $0` preproccesses an entire study

Usage:

  `basename $0` -c config_file -o bids_dir

Example:

  `basename $0` -c definitions_of_sequences.txt -o BIDS -t /scratch

Required arguments:

     -c:  description of the subjects (see KUL_multisubjects_dcm2bids)
     -o:  bids directory


Optional arguments:
    
    -t:  temporary directory (default = /tmp)

USAGE

    exit 1
}

# END LOCAL FUNCTIONS --------------



# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
silent=1
tmp=/tmp
#ncpu=18
#mem_mb=32768

# Set flags
conf_flag=0
bids_flag=0
tmp_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "c:o:" OPT; do

        case $OPT in
        c) #config_file
            conf_flag=1
            conf=$OPTARG
        ;;
        o) #bids_dir
            bids_flag=1
            bids_dir=$OPTARG
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


#docker system prune -a


# ----------- STEP 1 ---
kul_e2cl "Performing KUL_multisubjects_dcm2bids... " $log
KUL_multisubjects_dcm2bids.sh -d DICOM -c $conf -o $bids_dir -e -v




# ----------- STEP 2 ---
kul_e2cl "Performing multisubject mriqc... " $log
# we read the config file (and it may be csv, tsv or ;-seperated)
while IFS=$'\t,;' read -r BIDS_participant EAD dicom_zip config_file session comment; do
    
    
    if [ "$dicom_zip" = "dicom_zip" ]; then
        
        echo "first line" > /dev/null 2>&1

    else
        
        # check if already performed mriqc
        mriqc_dir_to_check=mriqc/sub-${BIDS_participant}
        
        #echo $mriqc_dir_to_check

        if [ ! -d $mriqc_dir_to_check ]; then

            kul_e2cl " Performing mriqc on participant $BIDS_participant" $log
            docker run --read-only --tmpfs /run --tmpfs /tmp --rm \
            -v ${cwd}/${bids_dir}:/data:ro -v ${cwd}/mriqc:/out \
            poldracklab/mriqc:latest \
            --participant_label $BIDS_participant \
            --n_procs $ncpu --ants-nthreads $ncpu --no-sub \
            /data /out participant 
        
        else
        
            kul_e2cl " mriqc of participant $BIDS_participant already done, skipping..." $log

        fi

    fi

done < $conf



kul_e2cl "Performing mriqc group summary" $log

# check if already performed group mriqc
mriqc_file_to_check=mriqc/group_bold.html
        
#echo $mriqc_file_to_check

if [ ! -f $mriqc_file_to_check ]; then
    
    docker run --read-only --tmpfs /run --tmpfs /tmp --rm \
            -v ${cwd}/${bids_dir}:/data:ro -v ${cwd}/mriqc:/out \
            poldracklab/mriqc:latest \
            /data /out group

else

    kul_e2cl " group mriqc already done, skipping..." $log

fi



# ----------- STEP 3 ---
kul_e2cl "Performing multisubject fmriprep... " $log
# we read the config file (and it may be csv, tsv or ;-seperated)
while IFS=$'\t,;' read -r BIDS_participant EAD dicom_zip config_file session comment; do
    
    
    if [ "$dicom_zip" = "dicom_zip" ]; then
        
        echo "first line" > /dev/null 2>&1

    else
        
        # check if already performed fmriprep
        fmriprep_dir_to_check=fmriprep/sub-${BIDS_participant}
        
        #echo $fmriprep_dir_to_check

        if [ ! -d $fmriprep_dir_to_check ]; then

            docker run --rm \
                -v ${cwd}/${bids_dir}:/data:delegated \
                -v ${cwd}/fmriprep:/out:delegated \
                -v ${cwd}/fmriprep_work:/scratch:delegated \
                -v /Users/xm52195/apps/Freesurfer_License/license.txt:/opt/freesurfer/license.txt \
                poldracklab/fmriprep:latest \
                --participant_label ${BIDS_participant} \
                -w /scratch \
                --nthreads $ncpu --omp-nthreads $ncpu \
                --mem_mb $mem_mb \
                --notrack \
                --stop-on-first-crash \
                --anat-only \
                /data /out \
                participant

        rm -fr -v ${cwd}/fmriprep_work
        
        else
        
            kul_e2cl "fmriprep of participant $BIDS_participant already done, skipping..." $log

        fi

    fi

done < $conf


STOP 


docker run --rm \
    -v /Users/xm52195/Dropbox/MIND2_T1/BIDS:/data:delegated \
    -v /Users/xm52195/Dropbox/MIND2_T1/fmriprep:/out:delegated \
    -v /Users/xm52195/Dropbox/MIND2_T1/fmriprep_work:/scratch:delegated \
    -v /Users/xm52195/apps/Freesurfer_License/license.txt:/opt/freesurfer/license.txt \
    poldracklab/fmriprep:latest \
    --participant_label P001 \
    -w /scratch \
    --nthreads 12 --omp-nthreads 6 \
    --mem_mb 32768 \
    --fs-no-reconall --notrack \
    --stop-on-first-crash \
    /data /out \
    participant

#docker run -it --read-only --tmpfs /run --tmpfs /tmp --rm -v /Users/xm52195/Dropbox/MIND2_T1/BIDS:/data:ro -v /Users/xm52195/Dropbox/MIND2_T1/mriqc:/out poldracklab/mriqc:latest --participant_label P002 P003 P004 P005 P006 P007 --no-sub /data /out participant 

#--n_procs 28 --ants-nthreads 28 --mem_gb 46 
