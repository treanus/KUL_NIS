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

     -R:  open mrview with results
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
result=0

# Set required options
p_flag=0
d_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:v:R" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        R) #results
			result=1
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
    bidsdir="BIDS/derivatives/KUL_compute/sub-$participant/KUL_anat_register_rigid"
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
        exit 1
    fi


    echo -e "\n\n"

}


function KUL_hd_glio_auto {
    
    # Segmentation of the tumor using HD-GLIO-AUTO
    
    # check if it needs to be performed

    hdglioinputdir="$kulderivativesdir/hdglio/input"
    hdgliooutputdir="$kulderivativesdir/hdglio"
    hdglio_output1=$hdgliooutputdir/output/lesion_perilesional_tissue
    hdglio_output2=$hdgliooutputdir/output/lesion_solid_tissue
    hdglio_output3=$hdgliooutputdir/output/lesion_total
    
    # only run if not yet done
    if [ ! -f ${hdglio_output3}.nii.gz ]; then

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

        hdglio_segmentation=$hdgliooutputdir/output/segmentation.nii.gz

        lesion_type_found=$(mrstats -output max $hdglio_segmentation)
        #echo $lesion_type_found
        
        if [ $lesion_type_found -ea 0 ];then

            kul_echo "hd-glio-auto did not find a lesion"
            # output an empty lesion mask
            mrcalc $hdglio_segmentation 1 -eq ${hdglio_output3}.nii.gz

        elif [ $lesion_type_found -ge 1 ];then

            kul_echo "hd-glio-auto found $hdglio_output1"
            mrcalc $hdglio_segmentation 1 -eq - | maskfilter - dilate -npass 5 -nthreads $ncpu ${hdglio_output1}_dil5.nii.gz -force
            maskfilter ${hdglio_output1}_dil5.nii.gz fill ${hdglio_output1}_dil5_fill.nii.gz -nthreads $ncpu -force
            maskfilter ${hdglio_output1}_dil5_fill.nii.gz erode ${hdglio_output1}.nii.gz -npass 5 -nthreads $ncpu -force
            cp ${hdglio_output1}.nii.gz ${hdglio_output3}.nii.gz

        elif [ $lesion_type_found -eq 2 ];then

            kul_echo "hd-glio-auto found $hdglio_output2"
            mrcalc $hdglio_segmentation 2 -eq - | maskfilter - dilate -npass 5 -nthreads $ncpu ${hdglio_output2}_dil5.nii.gz -force
            maskfilter ${hdglio_output2}_dil5.nii.gz fill ${hdglio_output2}_dil5_fill.nii.gz -nthreads $ncpu -force
            maskfilter ${hdglio_output2}_dil5_fill.nii.gz erode ${hdglio_output2}.nii.gz -npass 5 -nthreads $ncpu -force
            
            mrcalc $hdglio_segmentation 2 -le - | maskfilter - dilate -npass 5 -nthreads $ncpu ${hdglio_output3}_dil5.nii.gz -force
            maskfilter ${hdglio_output2}_dil5.nii.gz fill ${hdglio_output3}_dil5_fill.nii.gz -nthreads $ncpu -force
            maskfilter ${hdglio_output2}_dil5_fill.nii.gz erode ${hdglio_output3}.nii.gz -npass 5 -nthreads $ncpu -force
        
        fi
        
        #KUL_task_exec $verbose_level "compute lesion, perilesion zone & solid parts" "hdglioauto"
        
    else
        echo "Already done HD-GLIO-AUTO"
    fi

}

function KUL_resseg {
    
    # Segmentation of the tumor resection cavity using resseg
    
    resseginputdir1="$kulderivativesdir/resseg/input1"
    resseginputdir2="$kulderivativesdir/resseg/input2"
    ressegoutputdir="$kulderivativesdir/resseg/output"
    
    resseginput="T1"

    # only run if not yet done
    if [ ! -f "$ressegoutputdir/${resseginput}_cavity2.nii.gz" ]; then
        
        echo "Running resseg"
        
        # prepare the inputs
        mkdir -p $resseginputdir1
        mkdir -p $resseginputdir2
        mkdir -p $ressegoutputdir

        cp $T1w $resseginputdir1/${resseginput}.nii.gz

        # run resseg 1st time
        eval "$(conda shell.bash hook)"
        conda activate resseg
        resseg-mni -t $ressegoutputdir/${resseginput}_reg2mni.tfm \
            -r $ressegoutputdir/${resseginput}_reg2mni.nii.gz \
            $resseginputdir1/${resseginput}.nii.gz
        resseg -a 3 -t $ressegoutputdir/${resseginput}_reg2mni.tfm \
            -o $ressegoutputdir/${resseginput}_cavity1.nii.gz \
            $resseginputdir1/${resseginput}.nii.gz
        conda deactivate

        # run resseg 2nd time
        hdgliooutputdir="$kulderivativesdir/hdglio"
        # make a betted T1 as input
        maskfilter $hdgliooutputdir/output/mask.nii.gz dilate - -npass 15 | mrcalc $resseginputdir1/${resseginput}.nii.gz - -mul $resseginputdir2/${resseginput}.nii.gz -force
        eval "$(conda shell.bash hook)"
        conda activate resseg
        resseg -a 3 -t $ressegoutputdir/${resseginput}_reg2mni.tfm \
            -o $ressegoutputdir/${resseginput}_cavity2.nii.gz \
            $resseginputdir2/${resseginput}.nii.gz
        conda deactivate


    else
        echo "Already done resseg"
    fi

}

function KUL_fast {
    
    # Segmentation of the image using FSL FAST
    
    echo "Running FSL FAST"
    fastinputdir="$kulderivativesdir/fast/input"
    fastoutputdir="$kulderivativesdir/fast/output"
    hdgliooutputdir="$kulderivativesdir/hdglio/output"

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
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_anat_segment_tumor
mkdir -p $kulderivativesdir
globalresultsdir=$cwd/RESULTS/sub-$participant
mkdir -p $globalresultsdir/Lesion
mkdir -p $globalresultsdir/Anat

# setup
fastsurferoutputdir="$kulderivativesdir/fastsurfer"
input_hdglioauto1=$kulderivativesdir/hdglio/output/lesion_total.nii.gz
input_hdglioauto2=$kulderivativesdir/hdglio/output/lesion_perilesional_tissue.nii.gz
input_hdglioauto3=$kulderivativesdir/hdglio/output/lesion_solid_tissue.nii.gz
input_resseg1=$kulderivativesdir/resseg/output/T1_cavity1.nii.gz
input_resseg2=$kulderivativesdir/resseg/output/T1_cavity2.nii.gz
input_fastsurfer=$fastsurferoutputdir/$participant/mri/aparc.DKTatlas+aseg.deep.mgz

local_output_hdglio1=$kulderivativesdir/sub-${participant}_hdglio_lesion_total.nii.gz
local_output_hdglio2=$kulderivativesdir/sub-${participant}_hdglio_lesion_perilesional_tissue.nii.gz
local_output_hdglio3=$kulderivativesdir/sub-${participant}_hdglio_lesion_solid_tissue.nii.gz
global_output_hdglio1=$globalresultsdir/Lesion/sub-${participant}_hdglio_lesion_total.nii.gz
global_output_hdglio2=$globalresultsdir/Lesion/sub-${participant}_hdglio_lesion_perilesional_tissue.nii.gz
global_output_hdglio3=$globalresultsdir/Lesion/sub-${participant}_hdglio_lesion_solid_tissue.nii.gz

resseg_output=$kulderivativesdir/sub-${participant}_resseg_cavity_clean.nii.gz
fastsurferoutput=$kulderivativesdir/sub-${participant}_fastsurfer_ventricles.nii.gz

mrview_hdglio2=0
mrview_hdglio3=0

if [ -f $globalresultsdir/Lesion/sub-${participant}_tumor_segment.png ];then
    echo "Already done."
    exit
fi

# Check if fMRI and/or dwi data are present and/or to redo some processing
KUL_check_data

if [ $result -eq 0 ]; then
    # Get the data from KUL_anat_register_rigid
    ln -sf $cwd/BIDS/derivatives/KUL_compute/sub-$participant/KUL_anat_register_rigid/*.gz $globalresultsdir/Anat

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

    # link outputs
    ln -sf $input_hdglioauto1 $local_output_hdglio1
    ln -sf $input_hdglioauto1 $global_output_hdglio1

    if [ -f $input_hdglioauto2 ]; then
        ln -sf $input_hdglioauto2 $local_output_hdglio2
        ln -sf $input_hdglioauto2 $global_output_hdglio2
        mrview_hdglio2=1
    fi

    if [ -f $input_hdglioauto3 ]; then
        ln -sf $input_hdglioauto3 $local_output_hdglio3
        ln -sf $input_hdglioauto3 $global_output_hdglio3
        mrview_hdglio3=1
    fi


    # get the CSF from fastsurfer 

    # regrid fastsurfer to T1w space, get L and R ventricle and add them
    mrgrid $input_fastsurfer regrid -template $input_hdglioauto1 \
        $kulderivativesdir/sub-${participant}_fastsurfer_labels.nii.gz -interp nearest -force
    mrcalc $kulderivativesdir/sub-${participant}_fastsurfer_labels.nii.gz \
        4 -eq $kulderivativesdir/ventricle1.nii.gz -force
    mrcalc $kulderivativesdir/sub-${participant}_fastsurfer_labels.nii.gz \
        43 -eq $kulderivativesdir/ventricle2.nii.gz -force
    mrcalc $kulderivativesdir/ventricle1.nii.gz \
        $kulderivativesdir/ventricle2.nii.gz -add \
        $fastsurferoutput -force
    rm -rf $kulderivativesdir/ventricle1.nii.gz \
        $kulderivativesdir/ventricle2.nii.gz
    ln -sf $fastsurferoutput $globalresultsdir/Lesion/

    # Correct RESSEG
    # 1st try to determine if the 2 resseg runs overlap
    # subtract hdglio-lesion-total and csf from resseg 
    #   and clean

    # calculate the overlap between the 2 resseg runs
    mrcalc $input_resseg1 $input_resseg2 -add $kulderivativesdir/resseg/T1_cavity_combined.nii.gz -force
    cav_overlap=$(mrstats $kulderivativesdir/resseg/T1_cavity_combined.nii.gz -output max)

    if [ $cav_overlap -lt 2 ]; then
        # there is no overlap, we discard resseg output
        resseg_keep=""
        resseg_use=0
    else
        # there seems to be overlap, let's keep #1
        resseg_keep=$input_resseg1
        resseg_use=1

        ln -sf $resseg_keep $kulderivativesdir/sub-${participant}_resseg_cavity1_noclean.nii.gz

        # now compute resseg - hd_glio - csf
        mrcalc $resseg_keep $input_hdglioauto1 -subtract 0.9 -gt \
            $kulderivativesdir/resseg/resseg_subtracted.nii.gz -force
        mrcalc $kulderivativesdir/resseg/resseg_subtracted.nii.gz \
            $fastsurferoutput -subtract 0.9 -gt \
            $kulderivativesdir/resseg/resseg_subtracted_no_ventricles.nii.gz -force
        maskfilter $kulderivativesdir/resseg/resseg_subtracted_no_ventricles.nii.gz clean \
            $resseg_output -nthreads $ncpu -force
        ln -s $resseg_output $globalresultsdir/Lesion/sub-${participant}_resseg_cavity.nii.gz
    fi 


    # compute a refined whole lesion + cavity
    if [ $resseg_use -eq 1 ]; then
        mrcalc $resseg_output \
            $input_hdglioauto1 -add \
            $kulderivativesdir/tmp_lesion_full.nii.gz -force

    else

        ln -sf $input_hdglioauto1 \
            $kulderivativesdir/tmp_lesion_full.nii.gz

    fi

    maskfilter $kulderivativesdir/tmp_lesion_full.nii.gz dilate \
        $kulderivativesdir/tmp_lesion_full_dil.nii.gz -npass 5 -force
    maskfilter $kulderivativesdir/tmp_lesion_full_dil.nii.gz connect \
        $kulderivativesdir/tmp_lesion_full_dil_connect.nii.gz -force
    maskfilter $kulderivativesdir/tmp_lesion_full_dil_connect.nii.gz erode \
        $kulderivativesdir/sub-${participant}_lesion_full.nii.gz  -force

    #rm -rf $kulderivativesdir/tmp_*.nii.gz
    ln -sf $kulderivativesdir/sub-${participant}_lesion_full.nii.gz \
        $globalresultsdir/Lesion/sub-${participant}_lesion_full.nii.gz

else 

    if [ -f $resseg_output ];then
        resseg_use=1
    else
        resseg_use=0
    fi
    mrview_exit=""

fi

# create a figure
rm -f $globalresultsdir/Lesion/tmp*.png

underlay=$globalresultsdir/Anat/FLAIR_reg2_T1w.nii.gz
underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF)}')

mrview_hdglio2_overlay=""
mrview_hdglio3_overlay=""
if [ $mrview_hdglio2 -eq 1 ]; then
    mrview_hdglio2_overlay="-overlay.load $local_output_hdglio2 -overlay.opacity 0.5 -overlay.colour 0,1,1"
fi

if [ $mrview_hdglio3 -eq 1 ]; then
    mrview_hdglio3_overlay="-overlay.load $local_output_hdglio3 -overlay.opacity 0.5 -overlay.colour 0,0,1"
fi

if [ $resseg_use -eq 1 ]; then
    mrview_resseg_overlay="-overlay.load $resseg_output -overlay.opacity 0.5 -overlay.colour 0,255,0"
fi


if [ $result -eq 0 ]; then
    i=0
    voxel_index="-capture.folder $globalresultsdir/Lesion -capture.prefix tmp -noannotations "
    while [ $i -lt $underlay_slices ]
    do
        #echo Number: $i
        voxel_index="$voxel_index -voxel 0,0,$i -capture.grab"
        let "i+=7" 
    done
    mode_plane="-mode 1 -plane 2"
    mrview_exit="-exit"
else
    voxel_index=""
    mode_plane="-mode 2"
    mrview_exit=""
fi
#echo $voxel_index

cmd="mrview -load $underlay 
    $mode_plane \
    -overlay.load $globalresultsdir/Lesion/sub-${participant}_lesion_full.nii.gz -overlay.opacity 0.5 -overlay.colour 255,255,0 \
    -overlay.load $local_output_hdglio1 -overlay.opacity 0.5 -overlay.colour 85,0,255 \
    $mrview_hdglio2_overlay \
    $mrview_hdglio3_overlay \
    -overlay.load $fastsurferoutput -overlay.opacity 0.5 -overlay.colour 0,85,127 \
    $mrview_resseg_overlay \
    $voxel_index \
    -force \
    $mrview_exit"
#echo $cmd
eval $cmd

if [ $result -eq 0 ];then

    montage $globalresultsdir/Lesion/tmp*.png -mode Concatenate $globalresultsdir/Lesion/sub-${participant}_tumor_segment.png
    rm -f $globalresultsdir/Lesion/tmp*.png
fi

echo "Finished"
