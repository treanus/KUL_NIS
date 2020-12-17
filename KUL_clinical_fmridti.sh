#!/bin/bash
# Bash shell script to analyse clinical fMRI/DTI
#
# Requires matlab fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 16/12/2020
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` is a batch analysis of clinical fMRI/DTI data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -d DICOM/JohnDoe.zip

Required arguments:

     -p:  participant name
     -d:  dicom zip file (or directory)

Optional arguments:

     -v:  show output from commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given

# Set required options
p_flag=0
d_flag=0 

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:d:v" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        d) #dicomzip
			dicomzip=$OPTARG
            d_flag=1
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
if [ $d_flag -eq 0 ] ; then
	echo
	echo "Option -d is required: give the DICOM location." >&2
	echo
	exit 2
fi

# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then

	export MRTRIX_QUIET=1

fi

# --- MAIN ---

# convert the DICOM to BIDS
KUL_dcm2bids.sh -d $dicomzip -p $participant -c study_config/sequences.txt -e

# run fmriprep
cp study_config/run_fmriprep.txt KUL_LOG/run_fmriprep_$participant.txt
sed -i "s/BIDS_participants: /BIDS_participants: $participant/" KUL_LOG/run_fmriprep_$participant.txt
KUL_preproc_all.sh -e -c KUL_LOG/run_fmriprep_$participant.txt 

exit

# run SPM12
tcf="/DATA/test/study_config/test.m" #template config file
tjf="/DATA/test/study_config/test_job.m" #template job file
pcf="/DATA/test/study_config/test_$participant.m" #participant config file
pjf="/DATA/test/study_config/test_job_$participant.m" #participant job file
cp $tcf $pcf
cp $tjf $pjf
sed -i "s|###JOBFILE###|$pjf|" $pcf

fmridir="/DATA/test/SPM"
fmrifile="HAND_bold.nii"
fmriresults="/DATA/test/SPM/RESULTS"
sed -i "s|###FMRIDIR###|$fmridir|" $pjf
sed -i "s|###FMRIFILE###|$fmrifile|" $pjf
sed -i "s|###FMRIRESULTS###|$fmriresults|" $pjf

/usr/local/MATLAB/R2018a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('$pcf');exit ; "

