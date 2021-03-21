#!/bin/bash -e
# Bash shell script to crop T1w
#
# Requires MRtrix3
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 19/03/2021
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` crops T1w data in the z-axis, to remove the neck

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe 

Arguments:
     -a:  run on all BIDS participants
     -p:  participant name
     -n:  number of cpu to use (default 15)
     -v:  show output from commands

Documentation:
    USE AT OWN RISK!!!
    This is to be used with extreme caution!
    This will crop T1w images given a number of slices to delete from the bottom of the image, this is to delete extracranial tissue, such as the neck.
    It can be needed for qsiprep, since it cannot deal well with large fov acquired sagittal images.
    It will change you BIDS.
    It does not check anything, does not take into account any anatomical feature, just deletes data!
    USE AT OWN RISK!!!

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ncpu=15

# Set required options
p_flag=0
auto=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:n:av" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
        a) #verbose
			auto=1
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

# verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
    ants_verbose=0
    fs_silent=" > /dev/null 2>&1" 
fi

# --- MAIN ---
printf "\n\n\n"

# here we give the data
if [ $auto -eq 0 ]; then
    if [ -z "$session" ]; then
        fullsession1=""
        fullsession2=""
    else
        fullsession1="ses-${session}/"
        fullsession2="ses-${session}_"
    fi
    datadir="$cwd/BIDS/sub-${participant}/${fullsession1}anat"
    T1w=("$datadir/sub-${participant}_${fullsession2}T1w.nii.gz")
    #T2w=("$datadir/sub-${participant}_${fullsession2}T2w.nii.gz")
    #FLAIR=("$datadir/sub-${participant}_${fullsession2}FLAIR.nii.gz")
    #MTI=("$datadir/sub-${participant}_${fullsession2}MTI.nii.gz")
else
    T1w=($(find $cwd/BIDS -type f -name "*T1w.nii.gz" | sort ))
fi


for test_T1w in ${T1w[@]}; do

    my_cmd="mrgrid $test_T1w crop -axis 2 80,0 $test_T1w -force"
    echo $my_cmd
    eval $my_cmd

done
