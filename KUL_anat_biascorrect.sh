#!/bin/bash
# Bash shell script to register BIDS data to the T1w
#
# Requires ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 20/02/2022
version="0.1"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` corrects the signal bias in structural BIDS data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe

Required arguments:

     -p:  participant name

Optional arguments:

     -s:  session
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ants_verbose=1
verbose_level=1

# Set required options
p_flag=0
d_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:v:s:" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        s) #session
			session=$OPTARG
            s_flag=1
        ;;
        t) #type
			type=$OPTARG
		;;
        v) #verbose
            verbose_level=$OPTARG
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
if [ $p_flag -eq 0 ] ; then
	echo
	echo "Option -p is required: give the BIDS name of the participant." >&2
	echo
	exit 2
fi

if [ $s_flag -eq 1 ] ; then

    echo
    echo "Option -s $session is specified"
    echo
    participant="${participant}/ses-${session}"

fi

KUL_LOG_DIR="KUL_LOG/${script}/sub-${participant}"
mkdir -p $KUL_LOG_DIR

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
    ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
    ants_verbose=1
fi


# Functions

function KUL_check_data {
    
    echo -e "\n\nAn overview of the BIDS data:"
    bidsdir="BIDS/sub-$participant"
    T1w=($(find $bidsdir -name "*T1w.nii.gz" ! -name "*gadolinium*" -type f ))
    nT1w=${#T1w[@]}
    echo "  number of non-contrast T1w: $nT1w"
    cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
    ncT1w=${#cT1w[@]}
    echo "  number of contrast enhanced T1w: $ncT1w"
    FLAIR=($(find $bidsdir -name "*FLAIR.nii.gz" -type f ))
    nFLAIR=${#FLAIR[@]}
    echo "  number of FLAIR: $nFLAIR"
    FGATIR=($(find $bidsdir -name "*FGATIR.nii.gz" -type f ))
    nFGATIR=${#FGATIR[@]}
    echo "  number of FGATIR: $nFGATIR"
    T2w=($(find $bidsdir -name "*T2w.nii.gz" -type f ))
    nT2w=${#T2w[@]}
    echo "  number of T2w: $nT2w"
    SWI=($(find $bidsdir -name "*run-01_SWI.nii.gz" -type f ))
    nSWI=${#SWI[@]}
    SWIp=($(find $bidsdir -name "*run-02_SWI.nii.gz" -type f ))
    nSWIp=${#SWIp[@]}
    echo "  number of SWI magnitude: $nSWI"
    echo "  number of SWI phase: $nSWIp"

    echo -e "\n\n"

}

function KUL_biascorrect {
    bias_output=$kulderivativesdir/${source_mri_label}_bc.nii.gz
    if [ ! -f $bias_output ]; then 
        echo "  doing biascorrection on the ${td}"        
        N4BiasFieldCorrection --verbose $ants_verbose \
        -d 3 \
        -i $source_mri \
        -o $bias_output
    else
        echo "  biascorrection of the ${td} already done"
    fi 
}

function KUL_biascorrect_anat_images {

    if [ $ncT1w -gt 0 ];then
        source_mri_label="cT1w"
        source_mri=$cT1w
        task_in="KUL_biascorrect"
        KUL_task_exec $verbose_level "Bias correcting the $source_mri_label" "anat_biascorrect"
    fi
    if [ $nT1w -gt 0 ];then
        source_mri_label="T1w"
        source_mri=$T1w
        task_in="KUL_biascorrect"
        KUL_task_exec $verbose_level "Bias correcting the $source_mri_label" "anat_biascorrect"
    fi
    if [ $nT2w -gt 0 ];then
        source_mri_label="T2w"
        source_mri=$T2w
        task_in="KUL_biascorrect"
        KUL_task_exec $verbose_level "Bias correcting the $source_mri_label" "anat_biascorrect"
    fi
    if [ $nFLAIR -gt 0 ];then
        source_mri_label="FLAIR"
        source_mri=$FLAIR
        task_in="KUL_biascorrect"
        KUL_task_exec $verbose_level "Bias correcting the $source_mri_label" "anat_biascorrect"
    fi
    if [ $nFGATIR -gt 0 ];then
        source_mri_label="FGATIR"
        source_mri=$FGATIR
        task_in="KUL_biascorrect"
        KUL_task_exec $verbose_level "Bias correcting the $source_mri_label" "anat_biascorrect"
    fi
    if [ $nSWI -gt 0 ];then
        source_mri_label="SWI"
        source_mri=$SWI
        task_in="KUL_biascorrect"
        KUL_task_exec $verbose_level "Bias correcting the $source_mri_label" "anat_biascorrect"
    fi
    
}


# --- MAIN ---

# STEP 1 - SETUP

KUL_check_data

kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute/sub-$participant/KUL_anat_biascorrect
mkdir -p $kulderivativesdir


# STEP 2 - biascorrect all anatomical images
KUL_biascorrect_anat_images 

echo "Finished"
