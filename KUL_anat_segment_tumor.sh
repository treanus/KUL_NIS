#!/bin/bash
# Bash shell script to segment a tumour and/or resection cavity
#
# Requires HD-GLIO-AUTO, HD-BET, resseg
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 14/02/2022
version="0.1"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` segments a tumor and/or resection cavity using AI tools

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe

Required arguments:

     -p:  participant name

Optional arguments:

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

# Set required options
p_flag=0
d_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:n:v:" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        n) #ncpu
			ncpu=$OPTARG
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



# --- functions ---
function KUL_antsApply_Transform {
    antsApplyTransforms -d 3 --float 1 \
        --verbose 1 \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n Linear
}

function KUL_check_data {
    
    echo -e "\n\nAn overview of the bias corrected derivatives data:"
    bidsdir="BIDS/derivatives/KUL_compute/sub-$participant/KUL_register_rigid"
    if [ ! -d $bidsdir ]; then
        echo "No suitable data found in the derivatives folder"
        echo "Run KUL_anat_biascorrect and KUL_anat_register_rigid first"
        exit
    fi
    T1w=($(find -L $bidsdir -name "T1w.nii.gz" -type f ))
    nT1w=${#T1w[@]}
    echo "  number of non-contrast T1w: $nT1w"
    cT1w=($(find $bidsdir -name "cT1w_reg2_T1w.nii.gz" -type f ))
    ncT1w=${#cT1w[@]}
    echo "  number of contrast enhanced T1w: $ncT1w"
    FLAIR=($(find $bidsdir -name "FLAIR_reg2_T1w.nii.gz" -type f ))
    nFLAIR=${#FLAIR[@]}
    echo "  number of FLAIR: $nFLAIR"
    T2w=($(find $bidsdir -name "T2w_reg2_T1w.nii.gz" -type f ))
    nT2w=${#T2w[@]}
    echo "  number of T2w: $nT2w"

    # check hd-glio-auto requirements

    if [ $nT1w -lt 1 ] || [ $ncT1w -lt 1 ] || [ $nT2w -lt 1 ] || [ $nT1w -lt 1 ]; then
        echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
        echo " At least one is missing. Check the derivatives folder"
        
        read -p "Are you sure you want to continue? (y/n)? " answ
        if [[ "$answ" == "n" ]]; then
            exit 1
        fi
    fi


    echo -e "\n\n"

}


function KUL_hd_glio_auto {
    
    # Segmentation of the tumor using HD-GLIO-AUTO
    
    # check if it needs to be performed

    hdglioinputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/hdglio/input"
    hdgliooutputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/hdglio"
    
    # only run if not yet done
    if [ ! -f "$hdgliooutputdir/lesion_solid_tumour.nii.gz" ]; then

        # prepare the inputs
        mkdir -p $hdglioinputdir
        mkdir -p $hdgliooutputdir/output
        ln -s $cwd/$T1w $hdglioinputdir/T1.nii.gz
        ln -s $cwd/$cT1w $hdglioinputdir/CT1.nii.gz
        ln -s $cwd/$FLAIR $hdglioinputdir/FLAIR.nii.gz
        ln -s $cwd/$T2w $hdglioinputdir/T2.nii.gz
        
        # run HD-GLIO-AUTO using docker
        if [ ! -f /usr/local/KUL_apps/HD-GLIO-AUTO/scripts/run.py ]; then
            task_in="docker run --gpus all --mount type=bind,source=$hdglioinputdir,target=/input \
                --mount type=bind,source=$hdgliooutputdir/output,target=/output \
                jenspetersen/hd-glio-auto"
            hdglio_type="docker"
        else
            task_in="python /usr/local/KUL_apps/HD-GLIO-AUTO/scripts/run.py -i $hdglioinputdir -o $hdgliooutputdir/output"
            hdglio_type="local install"
        fi
        KUL_task_exec $verbose_level "HD-GLIO-AUTO using $hdglio_type" "hdglioauto"

        # compute some additional output
        task_in="maskfilter $hdgliooutputdir/output/segmentation.nii.gz dilate $hdgliooutputdir/output/lesion_dil5.nii.gz -npass 5 -nthreads $ncpu -force; \
            maskfilter $hdgliooutputdir/output/lesion_dil5.nii.gz fill $hdgliooutputdir/output/lesion_dil5_fill.nii.gz -nthreads $ncpu -force; \
            maskfilter $hdgliooutputdir/output/lesion_dil5_fill.nii.gz erode $hdgliooutputdir/lesion_total.nii.gz -npass 5 -nthreads $ncpu -force; \
            mrcalc $hdgliooutputdir/output/segmentation.nii.gz 1 -eq $hdgliooutputdir/lesion_perilesional_tissue.nii.gz -force; \
            mrcalc $hdgliooutputdir/output/segmentation.nii.gz 2 -eq $hdgliooutputdir/lesion_solid_tumour.nii.gz -force" 
        KUL_task_exec $verbose_level "compute lesion, perilesion zone & solid parts" "hdglioauto"
        
    else
        echo "Already done HD-GLIO-AUTO"
    fi

}

function KUL_resseg {
    
    # Segmentation of the tumor resection cavity using resseg
    
    resseginputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/resseg/input"
    ressegoutputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/resseg/output"
    
    resseginput="T1"

    # only run if not yet done
    if [ ! -f "$ressegoutputdir/${resseginput}_cavity.nii.gz" ]; then
        
        echo "Running resseg"
        
        # prepare the inputs
        mkdir -p $resseginputdir
        mkdir -p $ressegoutputdir

        cp $T1w $resseginputdir/${resseginput}.nii.gz

        # run resseg
        eval "$(conda shell.bash hook)"
        conda activate resseg
        resseg-mni -t $ressegoutputdir/${resseginput}_reg2mni.tfm \
            -r $ressegoutputdir/${resseginput}_reg2mni.nii.gz \
            $resseginputdir/${resseginput}.nii.gz
        resseg -t $ressegoutputdir/${resseginput}_reg2mni.tfm \
            -o $ressegoutputdir/${resseginput}_cavity.nii.gz \
            $resseginputdir/${resseginput}.nii.gz
        conda deactivate
    
    else
        echo "Already done resseg"
    fi

}

function KUL_fast {
    
    # Segmentation of the image using FSL FAST
    
    echo "Running FSL FAST"
    fastinputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/fast/input"
    fastoutputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/fast/output"
    hdgliooutputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/hdglio/output"

    # only run if not yet done
    if [ ! -f "$fastoutputdir/fast_seg.nii.gz" ]; then

        # prepare the inputs
        mkdir -p $fastinputdir
        mkdir -p $fastoutputdir
        
        ln -s $hdgliooutputdir/T1_r2s_bet_reg.nii.gz $fastinputdir/T1w.nii.gz
        ln -s $hdgliooutputdir/CT1_r2s_bet_reg.nii.gz $fastinputdir/cT1w.nii.gz
        ln -s $hdgliooutputdir/T2_r2s_bet_reg.nii.gz $fastinputdir/T2w.nii.gz
        ln -s $hdgliooutputdir/FLAIR_r2s_bet_reg.nii.gz $fastinputdir/FLAIR.nii.gz

        fast -S 4 -n 4 -H 0.1 -I 4 -l 20.0 -g \
            -o $fastoutputdir/fast \
            $fastinputdir/cT1w.nii.gz \
            $fastinputdir/T1w.nii.gz \
            $fastinputdir/T2w.nii.gz \
            $fastinputdir/FLAIR.nii.gz

    fi    

}

function KUL_fastsurfer {

    fastsurferoutputdir="$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/fastsurfer"
    if [ ! -f $fastsurferoutputdir/$participant/mri/aparc.DKTatlas+aseg.deep.mgz ]; then
        echo "Running segmentation-only fastsufer"
        my_cmd="$FASTSURFER_HOME/run_fastsurfer.sh \
            --sid $participant --sd $fastsurferoutputdir \
            --t1 $cwd/$T1w \
            --seg_only --py python --ignore_fs_version"
        eval $my_cmd
    else
        echo "Already run Fastsurfer"
    fi
}

# --- MAIN ---

# STEP 1 - Check to input data
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
mkdir -p $kulderivativesdir

# Check if fMRI and/or dwi data are present and/or to redo some processing
KUL_check_data

# STEP 2 - run HD-GLIO-AUTO
KUL_hd_glio_auto

# STEP 3 - run resseg
KUL_resseg

# STEP 4 - run Fastsurfer
KUL_fastsurfer

# STEP 5 - make final segmentations
#  HD-GLIO-AUTO nicely segments the tumor and perilesion tissue, but misses any surgical resection cavity
#  resseg find the surgical resection cavity, but overestimates and mislabels ventricles as cavity
#  Fastsurfer identifies the ventricles

echo "Running final segmentations"

input_hdglioauto=$hdgliooutputdir/lesion_total.nii.gz
input_resseg=$ressegoutputdir/T1_cavity.nii.gz
input_fastsurfer=$fastsurferoutputdir/$participant/mri/aparc.DKTatlas+aseg.deep.mgz
outputdir=$kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor

# get the CSF from fastsurfer
mrgrid $input_fastsurfer regrid -template $input_hdglioauto $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/fastsurfer.nii.gz -interp nearest -force
mrcalc $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/fastsurfer.nii.gz \
    4 -eq $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle1.nii.gz -force
maskfilter $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle1.nii.gz connect \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle1c.nii.gz -nthreads $ncpu -force
mrcalc $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/fastsurfer.nii.gz \
    43 -eq $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle2.nii.gz -force
maskfilter $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle2.nii.gz connect \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle2c.nii.gz -nthreads $ncpu -force
mrcalc $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle1c.nii.gz \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricle2c.nii.gz -add \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricles.nii.gz -force

# subtract hdglio-lesion-total and csf from resseg 
mrcalc $input_resseg $input_hdglioauto -subtract 0.9 -gt $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/resseg_subtracted.nii.gz -force
mrcalc $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/resseg_subtracted.nii.gz \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/ventricles.nii.gz -subtract 0.9 -gt \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/resseg_subtracted_no_ventricles.nii.gz -force

maskfilter $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/resseg_subtracted_no_ventricles.nii.gz clean \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/sub-participant_lesion_cavity.nii.gz -nthreads $ncpu -force

# compute the whole lesion 
mrcalc $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/sub-participant_lesion_cavity.nii.gz \
    $input_hdglioauto -add \
    $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/sub-participant_lesion_full.nii.gz -force
cp $input_hdglioauto $kulderivativesdir/sub-${participant}/KUL_anat_segment_tumor/sub-participant_lesion_solid.nii.gz

echo "Finished"
