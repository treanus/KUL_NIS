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

Optional arguments:

     -d:  dicom zip file (or directory)
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

# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then

	export MRTRIX_QUIET=1

fi

# --- MAIN ---

# convert the DICOM to BIDS
if [ ! -d "BIDS/sub-$participant" ];then
    KUL_dcm2bids.sh -d $dicomzip -p $participant -c study_config/sequences.txt -e
else
    echo "BIDS conversion already done"
fi

# run fmriprep
if [ ! -f fmriprep/sub-$participant.html ]; then
    cp study_config/run_fmriprep.txt KUL_LOG/run_fmriprep_$participant.txt
    sed -i "s/BIDS_participants: /BIDS_participants: $participant/" KUL_LOG/run_fmriprep_$participant.txt
    KUL_preproc_all.sh -e -c KUL_LOG/run_fmriprep_$participant.txt 
    rm -fr fmriprep_work_$participant
else
    echo "fmriprep already done"
fi


# run SPM12
echo "Preparing for SPM"
fmridatadir="$cwd/SPM/sub-$participant/fmridata"
mkdir -p $fmridatadir
fmriprepdir="fmriprep/sub-$participant/func"
searchtask="smoothAROMAnonaggr_bold.nii"
tasks=( $(find $fmriprepdir -name "*${searchtask}.gz" -type f -printf '%P\n') )
tcf="$kul_main_dir/share/spm12/spm12_fmri_stats_1run.m" #template config file
tjf="$kul_main_dir/share/spm12/spm12_fmri_stats_1run_job.m" #template job file
#echo ${tasks[@]}
matlab_exe=$(which matlab)

for task in ${tasks[@]}; do
    d1=${task#*_task-}
    shorttask=${d1%_space*}
    if [ ! "$shorttask" = "rest" ]; then
        #echo "$task -- $shorttask"
        echo "Analysing task $shorttask"
        fmrifile="${shorttask}_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii"
        cp $fmriprepdir/*$fmrifile.gz $fmridatadir
        gunzip $fmridatadir/*$fmrifile.gz
        fmriresults="$cwd/SPM/sub-$participant/stats_$shorttask"
        mkdir -p $fmriresults
        pcf="${fmriresults}.m" #participant config file
        pjf="${fmriresults}_job.m" #participant job file
        cp $tcf $pcf
        cp $tjf $pjf
        sed -i "s|###JOBFILE###|$pjf|" $pcf
        sed -i "s|###FMRIDIR###|$fmridatadir|" $pjf
        sed -i "s|###FMRIFILE###|$fmrifile|" $pjf
        sed -i "s|###FMRIRESULTS###|$fmriresults|" $pjf
        $matlab_exe -nodisplay -nosplash -nodesktop -r "run('$pcf');exit;"
        cp $fmriresults/spmT_0001.nii $cwd/SPM/sub-$participant/${shorttask}.nii
    fi
done

