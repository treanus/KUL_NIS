#!/bin/bash
# Bash shell script to register BIDS data to the T1w
#
# Requires ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 18/02/2022
version="0.1"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` register all structural BIDS data to the T1w (without Gd)

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe

Required arguments:

     -p:  participant name

    OR

     -t:  target image
     -s:  source image

Optional arguments:

     -s:  session
     -c:  use the bias corrected images as input (from KUL_anat_biascorrect.sh)
     -i:  interpolation type (1=BSpline, 2=NearestNeighbor; default=1)
     -o:  apply the transformation to other images
     -n:  number of cpu to use (default 15)
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ants_verbose=1
ncpu=15
verbose_level=1
bc_in=0
interpolation=1
target=""
source=""
other=""

# Set required options
p_flag=0
d_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:s:i:o:n:v:c" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        t) #target
			target=$OPTARG
		;;
        s) #source
			source=$OPTARG
		;;
        i) #interpolation
			interpolation=$OPTARG
		;;
        o) #other
			other=($OPTARG)
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
        v) #verbose
            verbose_level=$OPTARG
		;;
        c) #bc-input
            bc_in=1
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
if [ ! "$source" == "" ]; then

    direct=1

else

    direct=0

    if [ $p_flag -eq 0 ] ; then
        echo
        echo "Option -p is required: give the BIDS name of the participant." >&2
        echo
        exit 2
    fi

fi

if [ $interpolation -eq 1 ]; then
    interpolation_type="BSpline"
elif [ $interpolation -eq 1 ]; then
    interpolation_type="NearestNeighbor"
fi

KUL_LOG_DIR="KUL_LOG/${script}/sub-${participant}"
mkdir -p $KUL_LOG_DIR

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
	export MRTRIX_QUIET=1
    silent=1
    str_silent=" > /dev/null 2>&1" 
    ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
    silent=0
    str_silent="" 
    ants_verbose=1
fi

function KUL_antsApply_Transform {
    antsApplyTransforms -d 3 --float 1 \
        --verbose 1 \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n $interpolation_type
}


function KUL_check_data {
    
    if [ $bc_in -eq 0 ]; then
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

        # check the T1w
        if [ $nT1w -eq 0 ]; then
            echo "No T1w (without Gd) found. The script will not run since this is the registration target."
            exit 1
        fi 

    elif [ $bc_in -eq 1 ]; then

        echo -e "\n\nAn overview of the bias corrected derivatives data:"
        bidsdir="BIDS/derivatives/KUL_compute/sub-$participant/KUL_biascorrect"
        T1w=($(find $bidsdir -name "T1w_bc.nii.gz" -type f ))
        nT1w=${#T1w[@]}
        echo "  number of non-contrast T1w: $nT1w"
        cT1w=($(find $bidsdir -name "cT1w_bc.nii.gz" -type f ))
        ncT1w=${#cT1w[@]}
        echo "  number of contrast enhanced T1w: $ncT1w"
        FLAIR=($(find $bidsdir -name "FLAIR_bc.nii.gz" -type f ))
        nFLAIR=${#FLAIR[@]}
        echo "  number of FLAIR: $nFLAIR"
        FGATIR=($(find $bidsdir -name "FGATIR_bc.nii.gz" -type f ))
        nFGATIR=${#FGATIR[@]}
        echo "  number of FGATIR: $nFGATIR"
        T2w=($(find $bidsdir -name "T2w_bc.nii.gz" -type f ))
        nT2w=${#T2w[@]}
        echo "  number of T2w: $nT2w"
        SWI=($(find $bidsdir -name "SWI_bc.nii.gz" -type f ))
        nSWI=${#SWI[@]}
        SWIp=($(find BIDS/sub-$participant -name "*run-02_SWI.nii.gz" -type f ))
        nSWIp=${#SWIp[@]}
        echo "  number of SWI magnitude: $nSWI"
        echo "  number of SWI phase: $nSWIp"

    else
        echo "oeps no input found. Exitting"
        exit 1
    fi

    echo -e "\n\n"

}

function KUL_rigid_register {
    warp_field="${registeroutputdir}/${source_mri_label}_reg2_${target_mri_label}"
    output_mri="${kulderivativesdir}/${source_mri_label}_reg2_${target_mri_label}.nii.gz"
    #echo "Rigidly registering $source_mri to $target_mri"
    antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$warp_field,$output_mri] \
    --interpolation $interpolation_type \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$target_mri,$source_mri,1] \
    --transform Rigid[0.1] \
    --metric MI[$target_mri,$source_mri,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    #echo "Done rigidly registering $source_mri to $target_mri"
}

function KUL_register_anatomical_images {

    target_mri=$T1w
    target_mri_label="T1w"
    registeroutputdir="$kulderivativesdir/antsregister"
    mkdir -p $registeroutputdir
    ln -sf $cwd/$T1w $kulderivativesdir/T1w.nii.gz

    if [ $ncT1w -gt 0 ];then
        source_mri_label="cT1w"
        source_mri=$cT1w
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nT2w -gt 0 ];then
        source_mri_label="T2w"
        source_mri=$T2w
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nFLAIR -gt 0 ];then
        source_mri_label="FLAIR"
        source_mri=$FLAIR
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nFGATIR -gt 0 ];then
        source_mri_label="FGATIR"
        source_mri=$FGATIR
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nSWI -gt 0 ];then
        source_mri_label="SWI"
        source_mri=$SWI
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"

        input=$SWIp
        transform="${registeroutputdir}/${source_mri_label}_reg2_T1w0GenericAffine.mat"
        output="${kulderivativesdir}/${source_mri_label}_phase_reg2_T1w.nii.gz"
        reference=$target_mri
        task_in="KUL_antsApply_Transform"
        KUL_task_exec $verbose_level "Applying the rigid registration of SWIm to SWIp too" "anat_register_rigid"
    fi

}


# --- MAIN ---

# STEP 1 - SETUP

if [ $direct -eq 0 ]; then

    KUL_check_data

    kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute/sub-$participant/KUL_anat_register_rigid
    mkdir -p $kulderivativesdir


    # STEP 2 - register all anatomical other data to the T1w without contrast
    KUL_register_anatomical_images

else

    source_mri_label_tmp=$(basename $source)
    source_mri_label=${source_mri_label_tmp%%.*}
    #echo $source_mri_label
    source_mri=$source
    target_mri_label_tmp=$(basename $target)
    target_mri_label=${target_mri_label_tmp%%.*}
    #echo $target_mri_label
    target_mri=$target

    echo "Rigidly registering $source_mri_label to $target_mri_label (interpolation=$interpolation_type)"

    kulderivativesdir=$(pwd)
    registeroutputdir=$(pwd)

    KUL_rigid_register
    
    if [ ${#other[@]} -gt 0 ]; then

        for other_image in ${other[@]}; do

            input=$other_image
            output_tmp=$(basename $other_image)
            output="${output_tmp%%.*}_reg2_${target_mri_label}.nii.gz"
            transform="${warp_field}0GenericAffine.mat"
            echo  $transform
            reference=$target_mri
            KUL_antsApply_Transform

        done

    fi

fi

echo "Finished"
