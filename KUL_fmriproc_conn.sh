#!/bin/bash
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@uzleuven.be
#
# v0.1 - dd 19/01/2019 - jurassic version
version="v0.2 - dd 22/12/2021"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions


# Decide what kind of fMRI data you have in your BIDS dir
# carry out a single subject analysis of all this BOLD data automatically

# assuming most of the preproc is taken care of by fmriprep


# -----------------------------------  MAIN  ---------------------------------------------


# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs an automated melodic based fMRI analysis assuming:
    - having run fmriprep with aroma

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001

Required arguments:

     -p:  participant


Optional arguments:
     
     -s:  session
     -v:  verbose (0=silent, 1=normal, 2=verbose; default=1)


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
verbose_level=1


# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:v:" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            participant=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
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

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
	export MRTRIX_QUIET=1
    ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
    ants_verbose=1
fi

# Functions --------------------------------------------------------
function KUL_antsApply_Transform {
    if [ $KUL_DEBUG -gt 0 ]; then
        echo "input=$input"
        echo "output=$output"
        echo "transform=$transform"
        echo "reference=$reference"
    fi
    antsApplyTransforms -d 3 --float 1 \
        --verbose $ants_verbose \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n Linear
}

function KUL_compute_melodic {
    # run FSL Melodic
    echo "Computing Melodic"

    #if [ $silent -eq 1 ] ; then
    #    str_silent_melodic=" >> KUL_LOG/sub-${participant}_melodic.log"
    #fi

    # NOTE: searchtask has no res- token (fMRIPrep inserts "_res-2_" between
    # _space-... and _desc-... when a :res-2 output-space modifier is used) —
    # the * wildcards on either side of searchtask absorb that, same style as
    # KUL_fmri_denoise.sh's own find pattern.
    tasks=( $(find $fmriprepdir -name "*${searchtask}*.gz" -type f) )
    if [ ${#tasks[@]} -eq 0 ]; then
        echo "  ERROR: no BOLD files matched '*${searchtask}*.gz' under $fmriprepdir — melodic.done will NOT be created"
        return 1
    fi
    denoiseddir="$computedir/denoised"
    # we loop over the found tasks
    for task in ${tasks[@]}; do
        taskbase=$(basename $task)
        d1=${taskbase#*_task-}
        shorttask=${d1%%_space-*}
        #echo "$task -- $shorttask"
        echo " Analysing task $shorttask"

        # Denoise (confound regression + bandpass, then SUSAN smoothing) via
        # KUL_fmri_denoise.sh instead of the old raw-copy + ad-hoc fslmaths -s 3
        # smoothing block — avoids double-smoothing since SUSAN already runs as
        # the last step of KUL_fmri_denoise.sh.
        task_in="$kul_main_dir/KUL_fmri_denoise.sh --fmriprep ${cwd}/fmriprep --out $denoiseddir --method nilearn --space MNI152NLin2009cAsym --sub sub-${participant} --task $shorttask"
        KUL_task_exec $verbose_level "Denoising task $shorttask" "0_denoise_$shorttask"

        run_label="${taskbase%_desc-preproc_bold.nii.gz}"
        melodic_in_1="$denoiseddir/sub-${participant}/func/${run_label}_postproc_nilearn/${run_label}_desc-denoised_bold.nii.gz"
        if [ ! -f "$melodic_in_1" ]; then
            echo "  ERROR: expected denoised output not found, skipping task $shorttask: $melodic_in_1"
            continue
        fi

        fmriresults="$computedir/stats_$shorttask"
        mkdir -p $fmriresults
        # find the TR
        # edited by AR 04/11/2022
        tr=$(mrinfo $melodic_in_1 -spacing | cut -d " " -f 4)
        # make model and contrast
        # edited by AR 04/11/2022
        dyn=$(mrinfo $melodic_in_1 -size | cut -d " " -f 4)
        t_glm_con="$kul_main_dir/share/FSL/fsl_glm.con"
        t_glm_mat="$kul_main_dir/share/FSL/fsl_glm_${dyn}dyn.mat"
        # set dimensionality and model for rs-/a-fMRI
        if [[ $shorttask == *"rest"* ]]; then
            dim="--dim=15"
            model=""
        else
            dim=""
            model="--Tdes=$t_glm_mat --Tcon=$t_glm_con"
        fi

        task_in="melodic -i $melodic_in_1 -o $fmriresults --report --tr=$tr --Oall $model $dim"
        KUL_task_exec $verbose_level "Running melodic on $melodic_in_1" "1_melodic"

        # edited by AR 04/11/2022
        # now we compare to known networks
        # KUL_NIT_networks.nii.gz is on the classic FSL-standard MNI152 2mm
        # grid (91x109x91); melodic_IC.nii.gz is on fMRIPrep's actual
        # MNI152NLin2009cAsym:res-2 grid (97x115x97) -- different template,
        # fslcc requires matching grids and errors out ("Mismatch in image
        # dimensions") otherwise. Use the pre-resampled NLin2009cAsym-space
        # copy instead (see KUL_NIT_networks_space-MNI152NLin2009cAsym_res-2.nii.gz's
        # generation -- mrgrid regrid -interp nearest, a one-off, not
        # per-subject, same as the shared-atlas-resample pattern
        # step0_synthseg.sh already uses).
        mkdir -p $fmriresults/kul
        task_in="fslcc --noabs -p 3 -t .204 $kul_main_dir/atlases/Local/Sunaert2021/KUL_NIT_networks_space-MNI152NLin2009cAsym_res-2.nii.gz \
            $fmriresults/melodic_IC.nii.gz > $fmriresults/kul/kul_networks.txt"
        KUL_task_exec $verbose_level "Running fslcc for $melodic_in_1" "2_fslcc"


        while IFS=$' ' read network ic stat; do
            #echo $network
            network_name=$(sed "${network}q;d" $kul_main_dir/atlases/Local/Sunaert2021/KUL_NIT_networks.txt)
            #echo $network_name
            icfile="$fmriresults/stats/thresh_zstat${ic}.nii.gz"
            network_file="$fmriresults/kul/melodic_${network_name}_ic${ic}.nii.gz"
            #echo $icfile
            #echo $network_file
            mrcalc $icfile 2 -gt $icfile -mul $network_file -force

            # since Melodic analysis was in MNI space, we transform back in native space
            input=$network_file
            output=$globalresultsdir/rsfMRI_${shorttask}_${network_name}_ic${ic}.nii
            transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
            find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
            reference=${find_T1w[0]}
            KUL_antsApply_Transform $str_silent_melodic
        done < $fmriresults/kul/kul_networks.txt
    done
    echo "Done computing Melodic"
    touch KUL_LOG/sub-${participant}_melodic.done

}

# MAIN --------------------------------------------------------------
KUL_check_participant


#  setup variables
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
computedir="$kulderivativesdir/sub-$participant/FSL_melodic"
scriptsdir="$computedir/scripts"
fmriprepdir="fmriprep/sub-$participant/func"
globalresultsdir="$cwd/RESULTS/sub-$participant/Melodic"
searchtask="_space-MNI152NLin2009cAsym*_desc-preproc_bold.nii"

if [ $KUL_DEBUG -gt 0 ]; then
    echo "kulderivativesdir: $kulderivativesdir"
    echo "fmriprepdir: $fmriprepdir"
    echo "globalresultsdir: $globalresultsdir"
fi

mkdir -p $scriptsdir
mkdir -p $computedir/RESULTS
mkdir -p $globalresultsdir


if [ $verbose_level -lt 2 ] ; then
    str_silent_SPM=" >> KUL_LOG/$script/sub-${participant}_spm12.log"
fi
#echo $str_silent_SPM


if [ ! -f KUL_LOG/sub-${participant}_melodic.done ]; then

    # KUL_compute_melodic's denoise step below calls KUL_fmri_denoise.sh
    # --method nilearn -- unlike KUL_run_rsfMRI_networks.sh/
    # KUL_fmriproc_nilearn_new.sh, this script never activated the conda
    # env nilearn actually lives in, so it silently ran under whatever
    # environment KUL_clinical_fmridti.sh itself happened to be in and
    # failed with "nilearn not installed" (masked at the KUL_task_exec
    # level by an unrelated near-instant-failure duration-formatting
    # error, which is what actually shows up in the console).
    if ! conda env list | awk '{print $1}' | grep -qx "$KUL_PYFMRI_ENV"; then
        echo "ERROR: conda env '$KUL_PYFMRI_ENV' was not found (checked 'conda env list')." >&2
        echo "  Run the KUL_NIS installer's env-pyfmri section." >&2
        exit 1
    fi
    KUL_activate_conda_env "$KUL_PYFMRI_ENV"

    KUL_compute_melodic

else
    echo "Melodic analysis already done"
fi