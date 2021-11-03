#!/bin/bash -e
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@uzleuven.be
#
# v0.1 - dd 19/01/2019 - jurassic version
v="v0.1 - dd 19/01/2019"

# To Do
#  - make this thing work!
#  - This script should be able to do the following:
# Decide what kind of fMRI data you have in your BIDS dir
# carry out a single subject analysis of all this BOLD data automatically
# for rs-fMRI ICA, time series clustering (SLIC/NCUT), atlas based roi to roi analysis and connectome generation should be automated
# we will probably need to use some matlab and python functions and commands here and there
# For tb-fMRI we can do FEAT +/- all the previous
# Include an option for high resolution rendering in the three orthogonal planes and mosaic generation in TIFF/JPEG (slicer)
# For mTE-rsfMRI, after MEICA, run same as single echo rs-fMRI, check out tedana also

# assuming most of the preproc is taken care of by fmriprep


# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl (for logging)
#
# this script uses "preprocessing control", i.e. if some steps are already processed it will skip these

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs dMRI preprocessing.

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001 -n 6 

Required arguments:

     -p:  praticipant (BIDS name of the subject)


Optional arguments:
     
     -s:  session (BIDS session)
     -n:  number of cpu for parallelisation
     -a:  options to pass to atlas based parcellation
     -i:  options to pass to ICA
     -c:  options to pass to clustering
	 -t:  options to pass to FEAT
	 -r:  high_res rendering in three orthogonal planes (WIP: only for FEAT)
     -v:  verbose


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
ncpu=6 # default if option -n is not given
silent=1 # default if option -v is not given
# Specify additional options for FSL eddy


# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:n:d:t:e:v" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            subj=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        n) #ncpu
            ncpu=$OPTARG
        ;;
        a) #atlas_r2r
            atlas_opts=$OPTARG
        ;;
        i) #ica_opts
            ica_opts=$OPTARG
        ;;
        c) #cluster_opts
            clust_opts=$OPTARG
		;;
	    t) #feat_opts
	       	feat_opts=$OPTARG
        ;;
		r) # High res render
			HR_rend=$OPTARG
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