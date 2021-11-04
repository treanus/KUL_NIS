#!/bin/bash
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 22/01/2019 - alpha version
v="v0.1 - dd 22/01/2019"

# Warps FA (and other) maps to MNI space
#   - NNI warp is done by fmrirep on tht T1w anatomy
#   - Since the dwi_preproced_reg2T1w is in the same space FA map can be warped

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

`basename $0` warp FA (and ohter data) to MNI space.

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001 -n 6 

Required arguments:

     -p:  participant (anonymised name of the subject)

Optional arguments:

     -s:  session (of the participant)
     -n:  number of cpu for parallelisation
     -v:  show output from mrtrix commands


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
ncpu=6
silent=1

# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:n:s:vh" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            subj=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        n) #parallel
            ncpu=$OPTARG
        ;;
        v) #verbose
            silent=0
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
if [ $p_flag -eq 0 ] ; then 
    echo 
    echo "Option -s is required: give the anonymised name of a subject (this will create a directory subject_preproc with results)." >&2
    echo
    exit 2 
fi 

# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then 

    export MRTRIX_QUIET=1

fi

# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Some parallelisation
FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=log/log_${d}.txt

# --- MAIN ----------------

bids_subj=BIDS/sub-${subj}

# Either a session is given on the command line
# If not the session(s) need to be determined.
if [ $s_flag -eq 1 ]; then

    # session is given on the command line
    search_sessions=BIDS/sub-${subj}/ses-${ses}

else

    # search if any sessions exist
    search_sessions=($(find BIDS/sub-${subj} -type d | grep dwi))

fi    
 
num_sessions=${#search_sessions[@]}
    
echo "  Number of BIDS sessions: $num_sessions"
echo "    notably: ${search_sessions[@]}"


# ---- BIG LOOP for processing each session
for i in `seq 0 $(($num_sessions-1))`; do

    # set up directories 
    cd $cwd
    long_bids_subj=${search_sessions[$i]}
    #echo $long_bids_subj
    bids_subj=${long_bids_subj%dwi}


    kul_e2cl " Start processing $bids_subj" ${preproc}/${log}

    # CD to the directory of the preprocessed data
    preproc=dwiprep/sub-${subj}/$(basename $bids_subj) 
    #echo $preproc
    cd ${preproc}

    kul_e2cl "Welcome to KUL_dwiprep_anat $v - $d" ${log}


    # STEP 1 - Anatomical Processing ---------------------------------------------
    mkdir -p MNI

    function KUL_antsApply_Transform {

        antsApplyTransforms -d 3 --float 1 \
        --verbose 1 \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n Linear
    }



    fmriprep_subj=fmriprep/"sub-${subj}"
    fmriprep_anat="${cwd}/${fmriprep_subj}/anat/sub-${subj}_desc-preproc_T1w.nii.gz"
    fmriprep_anat_mask="${cwd}/${fmriprep_subj}/anat/sub-${subj}_desc-brain_mask.nii.gz"

    # transform the T1w into MNI space using fmriprep data
    input=T1w/T1w_BrainExtractionBrain.nii.gz
    output=MNI/sub-${subj}_T1w_space-MNI152NLin2009cAsym.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
    reference=${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz
    KUL_antsApply_Transform

    # transform the FA into MNI space using fmriprep data
    input=qa/fa_reg2T1w.nii.gz
    output=MNI/sub-${subj}_FA_space-MNI152NLin2009cAsym.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
    reference=${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz
    KUL_antsApply_Transform

    # transform the ADC into MNI space using fmriprep data
    input=qa/adc_reg2T1w.nii.gz
    output=MNI/sub-${subj}_ADC_space-MNI152NLin2009cAsym.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
    reference=${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz
    KUL_antsApply_Transform

    echo " Finished processing $bids_subj" 
    # ---- END of the BIG loop over sessions

done

# write a file to indicate that dwiprep_anat runned succesfully
#   this file will be checked by KUL_preproc_all

echo "done" > ../dwiprep_MNI_is_done.log

kul_e2cl "Finished " ${log}

