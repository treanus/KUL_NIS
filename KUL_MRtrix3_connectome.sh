#!/bin/bash -e
# Bash shell script to run mrtrix_connectome
#
# Requires docker
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/09/2021
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` runs MRtrix3_connectome tuned for KUL/UZLeuven data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -l 1 

   this will perform preprocessing on all participants/sessions in the BIDS directory

Required arguments:

     -l:  level (1=preproc, 2=participant, 3=group) 

Optional arguments:
     
     -p:  participant name
     -s:  session
     -g:  use gpu (does not work an MacOs)
     -n:  number of cpu to use (default 15)
     -v:  show output from commands

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ncpu=15
gpu=0
level=0
session=""

# Set required options
p_flag=0
l_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:n:l:gv" OPT; do

		case $OPT in
		p) #participant
			participant="$OPTARG"
            p_flag=1
		;;
		l) #level
			level=$OPTARG
            l_flag=1
		;;
		s) #session
			session=$OPTARG
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
		g) #gpu
			gpu=1
		;;
        v) #verbose
			silent=0
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
if [ $l_flag -eq 0 ] ; then
	echo
	echo "Option -l is required: give the processing level." >&2
	echo
	exit 2
fi


# ----  MAIN  ------------------------

#echo $participant
if [ -z $session ];then
	sessuf=""
	mrtrix_session_label=""
else
	sessuf="/ses-"
	mrtrix_session_label=" --session_label "$session" "
fi

if [ $gpu -eq 1 ]; then
	gpu_cmd1="--gpus all"
else
	gpu_cmd1=""
fi

outputdir="$cwd/MRtrix3_connectome"
scratchdir="$cwd/MRtrix3_connectome_sub-${participant}"

if [ $level -eq 1 ]; then
	
	test_file="$cwd/MRtrix3_connectome/MRtrix3_connectome-preproc/sub-${participant}/dwi/sub-${participant}_desc-preproc_dwi.nii.gz"
	
	if [ ! -f $test_file ];then
		my_cmd="docker run -i --rm \
			-v $cwd/BIDS:/bids_dataset \
			-v $outputdir:/output \
			$gpu_cmd1 \
			treanus/mrtrix3_connectome \
			/bids_dataset /output preproc \
			--participant_label "$participant" \
			$mrtrix_session_label \
			--topup_prefix synb0 \
			--output_verbosity 4 \
			--n_cpus $ncpu"
	else
		my_cmd="echo Already preprocessed"
	fi

elif [ $level -eq 2 ]; then


	my_cmd="docker run -i --rm \
			-v $cwd/BIDS:/bids_dataset \
			-v $outputdir:/output \
			$gpu_cmd1 \
			bids/mrtrix3_connectome \
			/bids_dataset /output participant \
			--participant_label "$participant" \
			$mrtrix_session_label \
			--output_verbosity 2 \
			--template_reg ants \
			--parcellation desikan
			--n_cpus $ncpu"

elif [ $level -eq 3 ]; then

	my_cmd="docker run -i --rm \
			-v $cwd/BIDS:/bids_dataset \
			-v $outputdir:/output \
			$gpu_cmd1 \
			bids/mrtrix3_connectome \
			/bids_dataset /output group \
			$mrtrix_session_label \
			--output_verbosity 4 \
			--n_cpus $ncpu"
			
fi

echo $my_cmd

mkdir -p $outputdir

eval $my_cmd


