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
        (in this case a rigid registration is performed of all structural images in the BIDS folder to the T1w)

    OR

     -t:  target image
     -s:  source image
        (in this case just 2 images are registered - rigidly or non-rigidly)


Optional arguments:

     -s:  session
     -c:  use the bias corrected images as input (from KUL_anat_biascorrect.sh)
     -m:  mask the inputs (1=source & target, 2=source only, 3=target only)
     -w:  register non-rigidly using antsRegistrationSyN
     -r:  type of registration
            1: rigid (default)
            2: affine
            3: elastic using antsRegistrationSyn (same as -w option, kept for compatibility)
     -d:  output directory
     -i:  interpolation type (1=BSpline, 2=NearestNeighbor, 3=Linear; default=1)
     -o:  apply the transformation to other images (put these between "")
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
warp2mni=0
output_dir=""
od=0
mask=0
reg_type=1

# Set required options
p_flag=0
d_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:s:i:o:n:v:d:m:r:cw" OPT; do

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
        d) #output_dir
			output_dir=$OPTARG
            od=1
		;;
        i) #interpolation
			interpolation=$OPTARG
		;;
        o) #other
			other=($OPTARG)
		;;
        m) #mask
			mask=$OPTARG
		;;
        r) #reg_type
			reg_type=$OPTARG
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
        w) #warp2mni
            warp2mni=1
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
elif [ $interpolation -eq 2 ]; then
    interpolation_type="NearestNeighbor"
elif [ $interpolation -eq 3 ]; then
    interpolation_type="Linear"
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

    echo "interpolation_type: $interpolation_type"
    antsApplyTransforms -d 3 --float 1 \
        --verbose 1 \
        -i $input \
        -o $output \
        -r $reference \
        $transform \
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
        bidsdir="BIDS/derivatives/KUL_compute/sub-$participant/KUL_anat_biascorrect"
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

    echo "source_mri2: $source_mri2"
    echo "target_mri2: $target_mri2"
    echo "outputwarp: $outputwarp"

    warp_field="${registeroutputdir}/${source_mri_label}_reg2_${target_mri_label}"
    output_mri="${kulderivativesdir}/${source_mri_label}_reg2_${target_mri_label}.nii.gz"
    #echo "Rigidly registering $source_mri to $target_mri"
    antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$warp_field,$output_mri] \
    --interpolation $interpolation_type \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$target_mri2,$source_mri2,1] \
    --transform Rigid[0.1] \
    --metric MI[$target_mri2,$source_mri2,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    #echo "Done rigidly registering $source_mri to $target_mri"
    ConvertTransformFile 3 ${registeroutputdir}/${source_mri_label}_reg2_${target_mri_label}0GenericAffine.mat \
        ${registeroutputdir}/${source_mri_label}_reg2_${target_mri_label}0GenericAffine.txt
}

function KUL_affine_register {

    warp_field="${registeroutputdir}/${source_mri_label}_reg2_${target_mri_label}"
    output_mri="${kulderivativesdir}/${source_mri_label}_reg2_${target_mri_label}.nii.gz"
    #echo "Rigidly registering $source_mri to $target_mri"
    antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$warp_field,$output_mri] \
    --interpolation $interpolation_type \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$target_mri2,$source_mri2,1] \
    --transform Affine[0.1] \
    --metric MI[$target_mri2,$source_mri2,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    #echo "Done rigidly registering $source_mri to $target_mri"
    ConvertTransformFile 3 ${registeroutputdir}/${source_mri_label}_reg2_${target_mri_label}0GenericAffine.mat \
        ${registeroutputdir}/${source_mri_label}_reg2_${target_mri_label}0GenericAffine.txt
}

function KUL_warp2MNI {

    echo "source_mri2: $source_mri2"
    echo "target_mri2: $target_mri2"
    echo "outputwarp: $outputwarp"

    local dim=3
    local its=10000x111110x11110
    local percentage=0.3
    local syn="100x100x50,-0.01,5"
    #local nm=${D}${nm1}_fixed_${nm2}_moving_setting_is_${mysetting}   # construct output prefix
    local nm=$outputwarp

    my_cmd="antsRegistration -d $dim -r [ $f, $m ,1 ]  \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t translation[ 0.1 ] \
                         -c [ $its,1.e-8,20 ]  \
                        -s 4x2x1vox  \
                        -f 6x4x2 -l 1 \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t rigid[ 0.1 ] \
                         -c [ $its,1.e-8,20 ]  \
                        -s 4x2x1vox  \
                        -f 3x2x1 -l 1 \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t affine[ 0.1 ] \
                         -c [ $its,1.e-8,20 ]  \
                        -s 4x2x1vox  \
                        -f 3x2x1 -l 1 \
                        -m mattes[  $f, $m , 0.5 , 32 ] \
                        -m cc[  $f, $m , 0.5 , 4 ] \
                         -t SyN[ .20, 3, 0 ] \
                         -c [ $syn ]  \
                        -s 1x0.5x0vox  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                       -o [ ${nm},${nm}_diff.nii.gz,${nm}_inv.nii.gz]"

    #antsApplyTransforms -d $dim -i $m -r $f -n linear -t ${nm}1Warp.nii.gz -t ${nm}0GenericAffine.mat -o ${nm}_warped.nii.gz

    my_cmd="antsRegistrationSyN.sh -d 3 -f ${target_mri2} -m ${source_mri2} \
        -o ${outputwarp} -n ${ncpu} -j 1 -t s $fs_silent"
    eval $my_cmd
    
    ConvertTransformFile 3 ${outputwarp}0GenericAffine.mat \
        ${outputwarp}0GenericAffine.txt

    antsApplyTransforms -d $dim -i ${source_mri2} -r ${target_mri2} -n BSpline  \
        -t ${outputwarp}1Warp.nii.gz -t ${outputwarp}0GenericAffine.mat -o ${outputwarp}.nii.gz

}

function KUL_register_anatomical_images {

    target_mri2=$T1w
    target_mri_label="T1w"
    registeroutputdir="$kulderivativesdir/antsregister"
    mkdir -p $registeroutputdir
    ln -sf $cwd/$T1w $kulderivativesdir/T1w.nii.gz

    if [ $ncT1w -gt 0 ];then
        source_mri_label="cT1w"
        source_mri2=$cT1w
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nT2w -gt 0 ];then
        source_mri_label="T2w"
        source_mri2=$T2w
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nFLAIR -gt 0 ];then
        source_mri_label="FLAIR"
        source_mri2=$FLAIR
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nFGATIR -gt 0 ];then
        source_mri_label="FGATIR"
        source_mri2=$FGATIR
        task_in="KUL_rigid_register"
        KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "anat_register_rigid"
    fi
    if [ $nSWI -gt 0 ];then
        source_mri_label="SWI"
        source_mri2=$SWI
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

    kulderivativesdir=$(pwd)
    if [ $od -eq 1 ]; then
        kulderivativesdir=$output_dir
        registeroutputdir=$output_dir
    else
        kulderivativesdir=$(pwd)
        registeroutputdir=$(pwd)
    fi

    run_hdbet=0
    if [ $warp2mni -eq 1 ] || [ $reg_type -eq 3 ]; then

        outputwarp=${registeroutputdir}/${source_mri_label}_warp2_${target_mri_label}
        outputwarp_test="${outputwarp}Warped.nii.gz"

        if [ ! -f $outputwarp_test ]; then
            run_hdbet=1
        fi
    
    else

        run_hdbet=1

    fi

    if [ $run_hdbet -eq 1 ]; then
        if [ $mask -eq 0 ]; then
            echo "no masking"
            source_mri2=${source_mri}
            target_mri2=${target_mri}
        elif [ $mask -eq 1 ]; then
            hd-bet -i $source  -o /tmp/${source_mri_label}_brain.nii.gz
            hd-bet -i $target  -o /tmp/${target_mri_label}_brain.nii.gz
            source_mri2=/tmp/${source_mri_label}_brain.nii.gz
            target_mri2=/tmp/${target_mri_label}_brain.nii.gz
        elif [ $mask -eq 2 ]; then
            hd-bet -i $source  -o /tmp/${source_mri_label}_brain.nii.gz
            source_mri2=/tmp/${source_mri_label}_brain.nii.gz
            target_mri2=${target_mri}
        elif [ $mask -eq 3 ]; then
            hd-bet -i $target  -o /tmp/${target_mri_label}_brain.nii.gz
            source_mri2=${source_mri}
            target_mri2=/tmp/${target_mri_label}_brain.nii.gz
        fi
    fi
    


    # KIND OF REGISTRATION - rigid or non-rigid?
    echo "reg_type: $reg_type"
    if [ $warp2mni -eq 1 ] || [ $reg_type -eq 3 ]; then

        #outputwarp=${registeroutputdir}/${source_mri_label}_warp2_${target_mri_label}
        #outputwarp_test="${outputwarp}Warped.nii.gz"
        #echo $outputwarp_test
        if [ ! -f $outputwarp_test ]; then

            echo "Warping non-rigidly $source_mri_label to $target_mri_label (interpolation=$interpolation_type)"
            mkdir -p $registeroutputdir
        
            KUL_warp2MNI

        else
            echo "Warping already done, skipping"
        fi

    elif [ $reg_type -eq 2 ]; then

        echo "Affinely registering $source_mri_label to $target_mri_label (interpolation=$interpolation_type)"
        KUL_affine_register
    
    else 

        echo "Rigidly registering $source_mri_label to $target_mri_label (interpolation=$interpolation_type)"
        KUL_rigid_register

    fi

    if [ ${#other[@]} -gt 0 ]; then

        for other_image in ${other[@]}; do

            input=$other_image
            output_tmp=$(basename $other_image)
            output="$registeroutputdir/${output_tmp%%.*}_reg2_${target_mri_label}.nii.gz"
            
            if [ $warp2mni -eq 1 ] || [ $reg_type -eq 3 ]; then

                transform="-t ${outputwarp}1Warp.nii.gz -t ${outputwarp}0GenericAffine.mat"
            
            else
                transform="-t ${warp_field}0GenericAffine.mat "
            fi
            #echo  $transform
            reference=$target_mri
            KUL_antsApply_Transform

        done

    fi

fi

echo "Finished"
