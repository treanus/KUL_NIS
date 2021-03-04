#!/bin/bash -e
# Sarah Cappelle & Stefan Sunaert
# 22/12/2020 - v1.0
# 18/02/2021 - v1.1 (adding calibration)
# This script is the first part of Sarah's Study1
# This script computes a T1/T2, T1/FLAIR and MTC (magnetisation transfer contrast) ratio
# 
# This scripts follows the rationales of 
#   - D. Pareto et al. AJNR 2020
#   - Ganzetti et al. Frontiers in human neuroscience 2014
#
# Starting from 3D-T1w, 3D-FLAIR and 2D-T2w scans we compute:
#  bias correct the images using N4biascorrect from ANTs
#  ANTs rigid coregister and reslice all images to the 3D-T1w (in isotropic 1 mm space)
#  create masked brain images using HD-BET
#  calibrate the images according to Ganzetti
#  compute a T1FLAIR_ratio, a T1T2_ratio and a MTR
v="1.1"

kul_main_dir=`dirname "$0"`
script=$0
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` computes a T1/T2 and a T1/FLAIR ratio image.

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -f 2 -m -v 

Required arguments:

     -p:  participant name


Optional arguments:

     -s:  session of the participant
     -a:  automatic mode (just work on all images in the BIDS folder)
     -n:  number of cpu to use (default 15)
     -m:  also run MS lesion segmentation using Freesurfer7 SamSeg
     -f:  also run fastsurfer (1=full, 3, segmentation with CC, 3=segmentation without CC)
     -v:  show output from commands

Documentation:
This script computes a T1/T2, T1/FLAIR and MTC (magnetisation transfer contrast) ratio, using BIDS organised data. 
A T1w is mandatory. 
A T1w/T2w ratio is computed if a T2w is present.
A T1w/FLAIR ratio is computed if a FLAIR is present.
A MTR is computed if an MTI pair is available.
It follows the rationale of Ganzetti et al. Frontiers in human neuroscience 2014 and D. Pareto et al. AJNR 2020.
It also calculates lesion using Freesurfer Samseg if T1w and FLAIR are present and option -m is given.
It also calculates a FastSurfer parcellation if option -f is used.

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
auto=0 # default if option -s is not given
silent=1 # default if option -v is not given
outputdir="$cwd/T1T2FLAIRMTR_ratio"
ms=0
ncpu=15
fastsurf=0

# Set required options
#p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:f:n:amv" OPT; do

		case $OPT in
		a) #automatic mode
			auto=1
		;;
		p) #participant
			participant=$OPTARG
		;;
        s) #session
			session=$OPTARG
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
        f) #fastsurfer
			fastsurf=$OPTARG
		;;
        m) #MS lesion segmentation
			ms=1
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
#if [ $p_flag -eq 0 ] ; then
#	echo
#	echo "Option -p is required: give the BIDS name of the participant." >&2
#	echo
#	exit 2
#fi

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ncpu}
ants_verbose=1
fs_silent=""
# verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
    ants_verbose=0
    fs_silent=" > /dev/null 2>&1" 
fi

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=log/log_${d}.txt

# --- FUNCTIONS ---

function KUL_antsApply_Transform {
    antsApplyTransforms -d 3 \
    --verbose $ants_verbose \
    -i $input \
    -o $output \
    -r $reference \
    -t $transform \
    -n Linear
}
function KUL_antsApply_Transform_MNI {
    antsApplyTransforms -d 3 \
        --verbose $ants_verbose \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform1 -t $transform2 \
        -n $interpolation_type
}

function KUL_iso_biascorrect {
    bias_output=$outdir/compute/${output}_iso_biascorrected.nii.gz
    if [ ! -f $bias_output ]; then 
        mrgrid $input regrid -voxel 1 $outdir/compute/${output}_iso.nii.gz -force
        bias_input=$outdir/compute/${output}_iso.nii.gz        
        N4BiasFieldCorrection --verbose $ants_verbose \
        -d 3 \
        -i $bias_input \
        -o $bias_output
    fi 
}

function KUL_MTI_reorient_crop_hdbet_iso {
    #fslreorient2std $input $outdir/compute/${output}
    #cp $input $outdir/compute/${output}.nii.gz
    #mrgrid $outdir/compute/${output}.nii.gz crop -axis 0 $crop_x,$crop_x -axis 2 $crop_z,0 \
    #    $outdir/compute/${output}_cropped.nii.gz -nthreads $ncpu -force
    #mrmath $outdir/compute/${output}.nii.gz mean $outdir/compute/${output}_mean.nii.gz \
    # -axis 3 -nthreads $ncpu -force
    #echo "  doing hd-bet on ${output}_mean.nii.gz"
    #if [ ! -f $outdir/compute/${output}_mean_brain.nii.gz ]; then
    #    my_cmd="hd-bet -i $outdir/compute/${output}_mean.nii.gz \
    #     -o $outdir/compute/${output}_mean_brain.nii.gz $fs_silent"
    #   eval my_cmd
    #fi
    #mrcalc $outdir/compute/${output}.nii.gz $outdir/compute/${output}_mean_brain_mask.nii.gz \
    #    -mul $outdir/compute/${output}_brain.nii.gz -force
    iso_output=$outdir/compute/${output}_iso.nii.gz
    mrgrid $input regrid -voxel 1 $iso_output -force
}

# Rigidly register the input to the T1w
function KUL_rigid_register {
antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$outdir/compute/${ants_type},$outdir/compute/${newname}] \
    --interpolation BSpline \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$outdir/compute/$ants_template,$outdir/compute/$ants_source,1] \
    --transform Rigid[0.1] \
    --metric MI[$outdir/compute/$ants_template,$outdir/compute/$ants_source,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
}

# Register and compute the ratio
function KUL_register_computeratio {
    base0=${test_T1w##*/}
    base=${base0%_T1w*}
    ants_type="${base}_rigid_${td}_reg2t1_"
    ants_template="${base}_T1w_iso_biascorrected.nii.gz"
    ants_source="${base}_${td}_iso_biascorrected.nii.gz"
    newname="${base}_${td}_iso_biascorrected_reg2T1w.nii.gz"
    finalname="${base}_${td}_reg2T1w.nii.gz"
    KUL_rigid_register

    # Calibrate
    echo "  performing (non)linear histogram matching"
    mrhistmatch \
        -mask_input $outdir/compute/${base}_eye_and_muscle.nii.gz \
        -mask_target /tmp/mni_eye_and_muscle.nii.gz \
        linear \
        $outdir/compute/$newname \
        $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t2_tal_nlin_sym_09a.nii \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-lin_reg2T1w.nii.gz -force
    mrhistogram -bin 100 -ignorezero  \
        $outdir/compute/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz \
        $outdir/compute/${base}_${td}_iso_biascorrected_reg2T1w_histogram.csv -force
    mrhistogram -bin 100 -ignorezero  \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-lin_reg2T1w.nii.gz \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-lin_reg2T1w_histogram.csv -force

    # do the NON-LINEAR histogram matching using non-brain tissue
    mrhistmatch \
        -mask_input $outdir/compute/${base}_MNI2subj_brainmask_mni_dilated_inverse.nii.gz \
        -mask_target $kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated_inverse.nii \
        nonlinear \
        $outdir/compute/$newname \
        $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t2_tal_nlin_sym_09a.nii \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz -force
    mrhistogram -bin 100 -ignorezero  \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin_reg2T1w_histogram.csv -force

    # Methods of Cappelle & Sunaert
    echo "  performing second (Cappelle) nonlinear histogram matching"

    # make a mask
    echo "  doing hd-bet on ${base}_${td}_iso_biascorrected_reg2T1w.nii.gz"
    if [ ! -f $outdir/compute/${base}_${td}_iso_biascorrected_brain_reg2T1w.nii.gz ]; then
        my_cmd="hd-bet -i $outdir/compute/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz \
         -o $outdir/compute/${base}_${td}_iso_biascorrected_brain_reg2T1w.nii.gz $fs_silent"
        eval $my_cmd
    fi 
    #maskfilter ${output} erode $outdir/compute/${base}_${td}_mask_eroded.nii.gz -nthreads $ncpu -force
    
    mask1="$outdir/compute/${base}_${td}_iso_biascorrected_brain_reg2T1w_mask.nii.gz"
    mask2="$outdir/compute/${base}_${td}_iso_biascorrected_brain_reg2T1w_inverted_mask.nii.gz"
    mrcalc $mask1 0.1 -lt $mask2 -force

    # find the ventricles by thresholding the T2w
    max_T2w=$(mrstats $outdir/compute/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz -output max)
    echo "  max signal of the T2w is $max_T2w"
    ventricles="$outdir/compute/${base}_T2w_iso_biascorrected_reg2T1w_ventricules.nii.gz"
    mrcalc $outdir/compute/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz $max_T2w 0.75 -mul -gt \
        $ventricles -force
    skull_and_ventricules="$outdir/compute/${base}_T2w_iso_biascorrected_reg2T1w_skull_and_ventricules.nii.gz"
    mrcalc $mask2 $ventricles -add $skull_and_ventricules -force
    reference_histo_mask="$kul_main_dir/atlasses/Local/Cappelle2021/T2wFLAIR_template_skull_and_ventricles_mask.nii.gz"
    reference_histo_image="$kul_main_dir/atlasses/Local/Cappelle2021/${td}_template.nii.gz"
    # do the scond NON-LINEAR histogram matching using non-brain tissue
    mrhistmatch \
        -mask_input $skull_and_ventricules \
        -mask_target $reference_histo_mask \
        nonlinear \
        $outdir/compute/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz \
        $reference_histo_image \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz -force
    mrhistogram -bin 100 -ignorezero  \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz \
        $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w_histogram.csv -force

    # Calculate the ratio
    mrcalc \
        $outdir/compute/${base}_T1w_iso_biascorrected.nii.gz \
         $outdir/compute/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz -divide \
        $outdir/compute/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
        $outdir/${base}_ratio-T1${td}_calib-none.nii.gz -nthreads $ncpu -force
    mrcalc \
        $outdir/compute/${base}_T1w_iso_biascorrected_calib-lin.nii.gz \
         $outdir/compute/${base}_${td}_iso_biascorrected_calib-lin_reg2T1w.nii.gz -divide \
        $outdir/compute/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
        $outdir/${base}_ratio-T1${td}_calib-lin.nii.gz -nthreads $ncpu -force
    mrcalc \
        $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin.nii.gz \
         $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz -divide \
        $outdir/compute/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
        $outdir/${base}_ratio-T1${td}_calib-nonlin.nii.gz -nthreads $ncpu -force
    mrcalc \
        $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin2.nii.gz \
         $outdir/compute/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz -divide \
        $outdir/compute/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
        $outdir/${base}_ratio-T1${td}_calib-nonlin2.nii.gz -nthreads $ncpu -force

    # Also warp to MNI space
    reference=$fix_im
    transform1="$outdir/warp2mni/${base}_T1w2MNI_1Warp.nii.gz"
    transform2="[$outdir/warp2mni/${base}_T1w2MNI_0GenericAffine.mat,0]"
    interpolation_type="Linear"
    input="$outdir/${base}_ratio-T1${td}_calib-none.nii.gz"
    output="$outdir/${base}_ratio-T1${td}_calib-none_MNI.nii.gz"
    KUL_antsApply_Transform_MNI
    input="$outdir/${base}_ratio-T1${td}_calib-lin.nii.gz"
    output="$outdir/${base}_ratio-T1${td}_calib-lin_MNI.nii.gz"
    KUL_antsApply_Transform_MNI
    input="$outdir/${base}_ratio-T1${td}_calib-nonlin.nii.gz"
    output="$outdir/${base}_ratio-T1${td}_calib-nonlin_MNI.nii.gz"
    KUL_antsApply_Transform_MNI
    input="$outdir/${base}_ratio-T1${td}_calib-nonlin2.nii.gz"
    output="$outdir/${base}_ratio-T1${td}_calib-nonlin2_MNI.nii.gz"
    KUL_antsApply_Transform_MNI
    reference=$ref_im
    input="$outdir/compute/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz"
    output="$outdir/${base}_${td}_reg2T1w_MNI.nii.gz"
    KUL_antsApply_Transform_MNI

    # save a final result
    cp $outdir/compute/$newname $outdir/$finalname
}

function KUL_MTI_register_computeratio {
    base0=${test_T1w##*/}
    base=${base0%_T1w*}
    # convert the 4D MTI to single 3Ds
    input="$outdir/compute/${base}_${td}_iso.nii.gz"
    S0="$outdir/compute/${base}_${td}_iso_S0.nii.gz"
    Smt="$outdir/compute/${base}_${td}_iso_Smt.nii.gz"
    mrconvert $input -coord 3 0 $S0 -force
    mrconvert $input -coord 3 1 $Smt -force
    # determine the registration
    ants_type="${base}_rigid_${td}_reg2t1_"
    ants_template="${base}_T1w_iso_biascorrected.nii.gz"
    ants_source="${base}_${td}_iso_Smt.nii.gz"
    newname="${base}_${td}_iso_Smt_reg2T1w.nii.gz"
    finalname="${base}_${td}_reg2T1w.nii.gz"
    KUL_rigid_register
    Smt="$outdir/compute/$newname"
    # Now apply the coregistration to the 4D MTI 
    input=$S0
    output="$outdir/compute/${base}_${td}_iso_S0_reg2T1w.nii.gz"
    S0=$output
    transform="$outdir/compute/${base}_rigid_${td}_reg2t1_0GenericAffine.mat"
    reference="$outdir/compute/${base}_T1w_iso_biascorrected.nii.gz"
    KUL_antsApply_Transform
    # make a better mask
    mask=$outdir/compute/${base}_T1w_iso_biascorrected_brain_mask.nii.gz
    # MTR formula: (S0 - Smt)/S0
    mrcalc $S0 $Smt -subtract $S0 -divide $mask -multiply \
     $outdir/${base}_ratio-MTC.nii.gz -nthreads $ncpu -force

    # Also warp to MNI space
    reference=$ref_im
    transform1="$outdir/warp2mni/${base}_T1w2MNI_1Warp.nii.gz"
    transform2="[$outdir/warp2mni/${base}_T1w2MNI_0GenericAffine.mat,0]"
    interpolation_type="Linear"
    input="$outdir/${base}_ratio-MTC.nii.gz"
    output="$outdir/${base}_ratio-MTC_MNI.nii.gz"
    KUL_antsApply_Transform_MNI

    cp $outdir/compute/$newname $outdir/$finalname
}

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
    T2w=("$datadir/sub-${participant}_${fullsession2}T2w.nii.gz")
    FLAIR=("$datadir/sub-${participant}_${fullsession2}FLAIR.nii.gz")
    MTI=("$datadir/sub-${participant}_${fullsession2}MTI.nii.gz")
else
    T1w=($(find BIDS -type f -name "*T1w.nii.gz" | sort ))
fi

d=0
t2=0
flair=0
mti=0
for test_T1w in ${T1w[@]}; do

    base0=${test_T1w##*/};base=${base0%_T1w*}
    local_participant=${base%_ses*}
    local_session="ses-${base##*ses-}"
    outdir=$outputdir/$local_participant/$local_session
    check_done="$outdir/compute/${base}.done"
    check_done2="$outdir/${base}_T1w.nii.gz"

    if [ ! -f $check_done2 ];then

        # Test whether T2 and/or FLAIR also exist
        test_T2w="${test_T1w%_T1w*}_T2w.nii.gz"
        if [ -f $test_T2w ];then
            #echo "The T2 exists"
            d=$((d+1))
            t2=1
        fi
        test_FLAIR="${test_T1w%_T1w*}_FLAIR.nii.gz"
        if [ -f $test_FLAIR ];then
            #echo "The FLAIR exists"
            d=$((d+1))
            flair=1
        fi
        test_MTI="${test_T1w%_T1w*}_MTI.nii.gz"
        if [ -f $test_MTI ];then
            #echo "The MTI exists"
            d=$((d+1))
            mti=1
        fi

        # If a T1w and a T2w and/or a FLAIR exists
        if [ $d -gt 0 ]; then
            mkdir -p $outdir/compute
            mkdir -p $outputdir/log
            kul_e2cl "KUL_T1T2FLAIR_ratio is processing $local_participant and session $local_session" ${outputdir}/${log}
            
            # for the T1w
            input=$test_T1w
            output=${test_T1w##*/}
            output=${output%%.*}
            echo "  doing biascorrection on image $output"
            KUL_iso_biascorrect
            
            # run fastsurfer       
            #echo "fastsurf: $fastsurf"     
            if [ $fastsurf -gt 0 ]; then
                mkdir -p $outdir/fs
                if [ $fastsurf -eq 1 ]; then
                    echo "  running full fastsufer"
                    my_cmd="$FASTSURFER_HOME/run_fastsurfer.sh \
                        --sid $base --sd $outdir/fs \
                        --t1 $cwd/$test_T1w \
                        --fs_license $FREESURFER_HOME/license.txt \
                        --vol_segstats --py python --parallel --threads $ncpu $fs_silent"
                elif [ $fastsurf -eq 2 ]; then
                    echo "  running segmentation-with-CC & stats fastsufer"
                    my_cmd="$FASTSURFER_HOME/run_fastsurfer.sh \
                        --sid $base --sd $outdir/fs \
                        --t1 $cwd/$test_T1w \
                        --seg_with_cc_only --vol_segstats --py python --ignore_fs_version --threads $ncpu $fs_silent"
                elif [ $fastsurf -eq 3 ]; then
                    echo "  running segmentation-only fastsufer"
                    my_cmd="$FASTSURFER_HOME/run_fastsurfer.sh \
                        --sid $base --sd $outdir/fs \
                        --t1 $cwd/$test_T1w \
                        --seg_only --py python --ignore_fs_version --threads $ncpu $fs_silent"
                fi
                eval $my_cmd
            fi

            # hd-bet brain extraction of the T1w
            echo "  doing hd-bet on ${output}_iso_biascorrected.nii.gz"
            if [ ! -f $outdir/compute/${output}_iso_biascorrected_brain.nii.gz ]; then 
                my_cmd="hd-bet -i $outdir/compute/${output}_iso_biascorrected.nii.gz \
                 -o $outdir/compute/${output}_iso_biascorrected_brain.nii.gz $fs_silent"
                eval $my_cmd
            fi

            # Spatially normalise (warp) the T1w to the MNI atlas
            fix_im="$kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a.nii"
            ref_im="$kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a_with_neck.nii"
            mov_im=$bias_output
            output="$outdir/warp2mni/${base}_T1w2MNI_"
            mkdir -p $outdir/warp2mni
            if [ ! -f $outdir/warp2mni/${base}_T1w2MNI_Warped.nii.gz ]; then 
                echo "  warping to MNI (takes about 20 minutes)"
                my_cmd="antsRegistrationSyN.sh -d 3 -f ${fix_im} -m ${mov_im} \
                 -o ${output} -n ${ncpu} -j 1 -t s $fs_silent"
                eval $my_cmd
            else
                echo "  skipping MNI spatial normalisation, since it exists already"
            fi
            
            # Warp the eye and muscle back to subject space
            reference=$mov_im
            transform1="$outdir/warp2mni/${base}_T1w2MNI_1InverseWarp.nii.gz"
            transform2="[$outdir/warp2mni/${base}_T1w2MNI_0GenericAffine.mat,1]"
            interpolation_type="NearestNeighbor"

            input="$kul_main_dir/atlasses/Ganzetti2014/eyemask.nii"
            output="$outdir/compute/${base}_eye.nii.gz"
            KUL_antsApply_Transform_MNI

            input="$kul_main_dir/atlasses/Ganzetti2014/tempmask.nii"
            output="$outdir/compute/${base}_tempmuscle.nii.gz"
            KUL_antsApply_Transform_MNI
            
            # Warp the brain_mask and its inverse to subject space
            input="$kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated.nii"
            output="$outdir/compute/${base}_MNI2subj_brainmask_mni_dilated.nii.gz"
            KUL_antsApply_Transform_MNI
            input="$kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated_inverse.nii"
            output="$outdir/compute/${base}_MNI2subj_brainmask_mni_dilated_inverse.nii.gz"
            KUL_antsApply_Transform_MNI

            # sum the masks
            echo "  performing linear histogram matching"
            mrcalc $outdir/compute/${base}_eye.nii.gz $outdir/compute/${base}_tempmuscle.nii.gz -add \
             $outdir/compute/${base}_eye_and_muscle.nii.gz -force
            mrcalc $kul_main_dir/atlasses/Ganzetti2014/eyemask.nii $kul_main_dir/atlasses/Ganzetti2014/tempmask.nii -add \
             /tmp/mni_eye_and_muscle.nii.gz -force
            
            # do the LINEAR histogram matching using eye/muscle tissue
            mrhistmatch \
             -mask_input $outdir/compute/${base}_eye_and_muscle.nii.gz \
             -mask_target /tmp/mni_eye_and_muscle.nii.gz \
             linear \
             $outdir/compute/${base}_T1w_iso_biascorrected.nii.gz \
             $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a.nii \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-lin.nii.gz -force

            mrhistogram -bin 100 -ignorezero \
             $outdir/compute/${base}_T1w_iso_biascorrected.nii.gz \
             $outdir/compute/${base}_T1w_iso_biascorrected_histogram.csv -force
            mrhistogram -bin 100 -ignorezero \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-lin.nii.gz \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-lin_histogram.csv -force

            # do the NON-LINEAR histogram matching using non-brain tissue
            echo "  performing nonlinear histogram matching"
            mrhistmatch \
             -mask_input $outdir/compute/${base}_MNI2subj_brainmask_mni_dilated_inverse.nii.gz \
             -mask_target $kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated_inverse.nii \
             nonlinear \
             $outdir/compute/${base}_T1w_iso_biascorrected.nii.gz \
             $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a.nii \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin.nii.gz -force

            mrhistogram -bin 100 -ignorezero  \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin.nii.gz \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin_histogram.csv -force
             

            # Methods of Cappelle & Sunaert
            echo "  performing second (Cappelle) nonlinear histogram matching"
            mask1="$outdir/compute/${base}_T1w_iso_biascorrected_brain_mask.nii.gz"
            mask2="$outdir/compute/${base}_T1w_iso_biascorrected_brain_inverted_mask.nii.gz"
            mrcalc $mask1 0.1 -lt $mask2 -force

            reference_histo_mask="$kul_main_dir/atlasses/Local/Cappelle2021/T1w_template_brain_mask_inverse.nii.gz"
            reference_histo_image="$kul_main_dir/atlasses/Local/Cappelle2021/T1w_template.nii.gz"

            # do the scond NON-LINEAR histogram matching using non-brain tissue
            mrhistmatch \
             -mask_input $mask2 \
             -mask_target $reference_histo_mask \
             nonlinear \
             $outdir/compute/${base}_T1w_iso_biascorrected.nii.gz \
             $reference_histo_image \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin2.nii.gz -force

            mrhistogram -bin 100 -ignorezero  \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin2.nii.gz \
             $outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin2_histogram.csv -force

            # Warp into MNI space too
            echo "  warping T1w to MNI"
            reference=$ref_im
            transform1="$outdir/warp2mni/${base}_T1w2MNI_1Warp.nii.gz"
            transform2="[$outdir/warp2mni/${base}_T1w2MNI_0GenericAffine.mat,0]"
            interpolation_type="LanczosWindowedSinc"

            input="$outdir/compute/${base}_T1w_iso_biascorrected_calib-lin.nii.gz"
            output="$outdir/compute/${base}_T1w_iso_biascorrected_calib-lin_MNI.nii.gz"
            KUL_antsApply_Transform_MNI

            input="$outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin.nii.gz"
            output="$outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin_MNI.nii.gz"
            KUL_antsApply_Transform_MNI

            input="$outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin2.nii.gz"
            output="$outdir/compute/${base}_T1w_iso_biascorrected_calib-nonlin2_MNI.nii.gz"
            KUL_antsApply_Transform_MNI

            input="$outdir/compute/${base}_T1w_iso_biascorrected.nii.gz"
            output="$outdir/compute/${base}_T1w_iso_biascorrected_MNI.nii.gz"

            KUL_antsApply_Transform_MNI
            # TODO - add mask here

            cp $outdir/compute/${base}_T1w_iso_biascorrected.nii.gz $outdir/${base}_T1w.nii.gz
            cp $outdir/compute/${base}_T1w_iso_biascorrected_MNI.nii.gz $outdir/${base}_T1w_MNI.nii.gz


            if [ $t2 -eq 1 ];then
                input=$test_T2w
                output=${test_T2w##*/}
                output=${output%%.*}
                crop_x=0
                crop_z=0
                echo "  doing biascorrection of the T2w"
                KUL_iso_biascorrect

                td="T2w"
                echo "  coregistering T2 to T1 and computing the ratio"
                KUL_register_computeratio
            fi

            if [ $flair -eq 1 ];then
                input=$test_FLAIR
                output=${test_FLAIR##*/}
                output=${output%%.*}
                crop_x=0
                crop_z=0
                echo "  doing biascorrection of the FLAIR"
                KUL_iso_biascorrect
                
                td="FLAIR"
                echo "  coregistering FLAIR to T1 and computing the ratio"
                KUL_register_computeratio
            fi

            # if MS lesion segmentation
            if [ $ms -eq 1 ];then
                if [ $flair -eq 1 ];then
                    MSlesion="$outdir/${base}_MSLesion.nii.gz"
                    if [ ! -f $MSlesion ];then
                        echo "  running samseg (takes about 20 minutes)"
                        T1w_iso="$outdir/${base}_T1w.nii.gz"
                        FLAIR_reg2T1w="$outdir/${base}_FLAIR_reg2T1w.nii.gz"
                        my_cmd="run_samseg --input $T1w_iso $FLAIR_reg2T1w --pallidum-separate \
                        --lesion --lesion-mask-pattern 0 1 --output $outdir/samseg \
                        --threads $ncpu $fs_silent"
                        eval $my_cmd
                        SamSeg="$outdir/samseg/seg.mgz"
                        mrcalc $SamSeg 99 -eq $MSlesion -force -nthreads $ncpu
                    else
                        echo "  already ran samseg"
                    fi
                else
                    echo "  Warning! No Flair available to do lesion MS segmentation"
                fi        
            fi

            if [ $mti -eq 1 ];then
                input=$test_MTI
                output=${test_MTI##*/}
                output=${output%%.*}
                crop_x=0
                crop_z=0
                echo "  doing hd-bet of the MTI"
                KUL_MTI_reorient_crop_hdbet_iso
                td="MTI"
                echo "  coregistering MTI to T1 and computing the MTC ratio"
                KUL_MTI_register_computeratio
            fi

            #rm -fr $outdir/compute/${base}*.gz
            #cd T1T2FLAIRMTR_ratio/compute/
            #rm -fr *std.nii.gz *iso.nii.gz *ted.nii.gz *reg2T1w.nii.gz *mask.nii.gz *eye.nii.gz *muscle.nii.gz *.csv *MTI* *MNI.nii.gz *verse.nii.gz *inv.nii.gz *brain.nii.gz
            #rm -fr *std.nii.gz *iso.nii.gz *ted.nii.gz *reg2T1w.nii.gz *eye.nii.gz *muscle.nii.gz *.csv *MTI* *MNI.nii.gz *verse.nii.gz *inv.nii.gz *brain.nii.gz
            #cd ../..

            touch $check_done

            echo " done"

        else
            echo " Nothing to do here"
        fi

    else
        echo "  $base already done"
    fi

done
