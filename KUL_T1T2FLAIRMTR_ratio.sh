#!/bin/bash -e 
# Sarah Cappelle & Stefan Sunaert
# 22/12/2020
# This script is the first part of Sarah's Study1
# This script computes a T1/T2, T1/FLAIR and MTC (magnetisation transfer contrast) ratio
# 
# This scripts follows the rationale of D. Pareto et al. AJNR 2020
# Starting from 3D-T1w, 3D-FLAIR and 2D-T2w scans we compute:
#  create masked brain images using HD-BET
#  bias correct the images using N4biascorrect from ANTs
#  ANTs rigid coregister and reslice all images to the 3D-T1w (in isotropic 1 mm space)
#  compute a T1FLAIR_ratio, a T1T2_ratio and a MTR

kul_main_dir=`dirname "$0"`
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

  `basename $0` -p JohnDoe -v 

Required arguments:

	 -p:  participant name


Optional arguments:

     -s:  session of the participant
     -a:  automatic mode (just work on all images in the BIDS folder)
     -v:  show output from commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
auto=0 # default if option -s is not given
silent=1 # default if option -v is not given
outputdir="T1T2FLAIRMTR_ratio"

# Set required options
#p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:av" OPT; do

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

ants_verbose=1
# verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
    ants_verbose=0
fi

# --- FUNCTIONS ---

function KUL_antsApply_Transform {
    antsApplyTransforms -d 3 --float 1 \
    --verbose $ants_verbose \
    -i $input \
    -o $output \
    -r $reference \
    -t $transform \
    -n Linear
}

function KUL_reorient_crop_hdbet_biascorrect_iso {
    fslreorient2std $input $outputdir/compute/${output}_std
    mrgrid $outputdir/compute/${output}_std.nii.gz crop -axis 0 $crop_x,$crop_x -axis 2 $crop_z,0 \
            $outputdir/compute/${output}_std_cropped.nii.gz -force
    result=$(hd-bet -i $outputdir/compute/${output}_std_cropped.nii.gz -o $outputdir/compute/${output}_std_cropped_brain.nii.gz 2>&1)
    if [ $silent -eq 0 ]; then
        echo $result
    fi 
    bias_input=$outputdir/compute/${output}_std_cropped_brain.nii.gz
    mask=$outputdir/compute/${output}_std_cropped_brain_mask.nii.gz
    bias_output=$outputdir/compute/${output}_std_cropped_brain_biascorrected.nii.gz
    N4BiasFieldCorrection --verbose $ants_verbose \
     -d 3 \
     -i $bias_input \
     -o $bias_output
    iso_output=$outputdir/compute/${output}_std_cropped_brain_biascorrected_iso.nii.gz
    mrgrid $bias_output regrid -voxel 1 $iso_output -force
    iso_output2=$outputdir/compute/${output}_std_cropped_brain_mask_iso.nii.gz
    mrgrid $mask regrid -voxel 1 $iso_output2 -force
    mv $iso_output $outputdir
}

function KUL_MTI_reorient_crop_hdbet_iso {
    fslreorient2std $input $outputdir/compute/${output}_std
    mrgrid $outputdir/compute/${output}_std.nii.gz crop -axis 0 $crop_x,$crop_x -axis 2 $crop_z,0 $outputdir/compute/${output}_std_cropped.nii.gz -force
    mrmath $outputdir/compute/${output}_std_cropped.nii.gz mean $outputdir/compute/${output}_mean_std_cropped.nii.gz -axis 3
    result=$(hd-bet -i $outputdir/compute/${output}_mean_std_cropped.nii.gz -o $outputdir/compute/${output}_mean_std_cropped_brain.nii.gz 2>&1)
    if [ $silent -eq 0 ]; then
        echo $result
    fi 
    mrcalc $outputdir/compute/${output}_std_cropped.nii.gz $outputdir/compute/${output}_mean_std_cropped_brain_mask.nii.gz \
        -mul $outputdir/compute/${output}_std_cropped_brain.nii.gz -force
    iso_output=$outputdir/compute/${output}_std_cropped_brain_iso.nii.gz
    mrgrid $outputdir/compute/${output}_std_cropped_brain.nii.gz regrid -voxel 1 $iso_output -force
}

# Rigidly register the input to the T1w
function KUL_rigid_register {
antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$outputdir/compute/${ants_type},$outputdir/compute/${newname}] \
    --interpolation BSpline \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$outputdir/compute/$ants_template,$outputdir/compute/$ants_source,1] \
    --transform Rigid[0.1] \
    --metric MI[$outputdir/compute/$ants_template,$outputdir/compute/$ants_source,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox

    # also apply the registration to the mask
    input=$mask
    output=${mask##*/}
    output=${output%%.*}
    output=$outputdir/compute/${output}_reg2T1w.nii.gz
    transform=$outputdir/compute/${ants_type}0GenericAffine.mat
    reference=$outputdir/compute/$ants_template
    #echo "input $input"
    #echo "output $output"
    #echo "transform $transform"
    #echo "reference $reference"
    KUL_antsApply_Transform
}

# Register and compute the ratio
function KUL_register_computeratio {
    base0=${test_T1w##*/}
    base=${base0%_T1w*}
    ants_type="${base}_rigid_${td}_reg2t1_"
    ants_template="${base}_T1w_std_cropped_brain_biascorrected_iso.nii.gz"
    ants_source="${base}_${td}_std_cropped_brain_biascorrected_iso.nii.gz"
    newname="${base}_${td}_std_cropped_brain_biascorrected_iso_reg2T1w.nii.gz"
    KUL_rigid_register
    # make a better mask
    maskfilter ${output} erode $outputdir/compute/${base}_${td}_mask_eroded.nii.gz -force
    #mrcalc $outputdir/compute/$ants_template $outputdir/compute/$newname -divide \
    #    $outputdir/${base}_T1${td}_ratio_a.nii.gz
    mrcalc $outputdir/compute/$ants_template $outputdir/compute/$newname -divide \
        $outputdir/compute/${base}_${td}_mask_eroded.nii.gz -multiply $outputdir/${base}_T1${td}_ratio.nii.gz -force
    mv $newname $outputdir
}

function KUL_MTI_register_computeratio {
    base0=${test_T1w##*/}
    base=${base0%_T1w*}
    # convert the 4D MTI to single 3Ds
    input="$outputdir/compute/${base}_${td}_std_cropped_brain_iso.nii.gz"
    S0="$outputdir/compute/${base}_${td}_std_cropped_brain_iso_S0.nii.gz"
    Smt="$outputdir/compute/${base}_${td}_std_cropped_brain_iso_Smt.nii.gz"
    mrconvert $input -coord 3 0 $S0 -force
    mrconvert $input -coord 3 1 $Smt -force
    # determine the registration
    ants_type="${base}_rigid_${td}_reg2t1_"
    ants_template="${base}_T1w_std_cropped_brain_biascorrected_iso.nii.gz"
    ants_source="${base}_${td}_std_cropped_brain_iso_Smt.nii.gz"
    newname="${base}_${td}_std_cropped_brain_iso_Smt_reg2T1w.nii.gz"
    KUL_rigid_register
    Smt="$outputdir/compute/$newname"
    # Now apply the coregistration to the 4D MTI 
    input=$S0
    output="$outputdir/compute/${base}_${td}_std_cropped_brain_iso_S0_reg2T1w.nii.gz"
    S0=$output
    transform="$outputdir/compute/${base}_rigid_${td}_reg2t1_0GenericAffine.mat"
    reference="$outputdir/compute/${base}_T1w_std_cropped_brain_biascorrected_iso.nii.gz"
    KUL_antsApply_Transform
    # make a better mask
    mask=$outputdir/compute/${base}_T1w_std_cropped_brain_mask_iso.nii.gz
    # MTR formula: (S0 - Smt)/S0
    mrcalc $S0 $Smt -subtract $S0 -divide $mask -multiply \
     $outputdir/${base}_MTC_ratio.nii.gz -force
}

# --- MAIN ---
printf "\n\n\n"

# here we give the data
if [ $auto -eq 0 ]; then
    datadir="$cwd/BIDS/sub-${participant}/ses-$session/anat"
    T1w=("$datadir/sub-${participant}_ses-${session}_T1w.nii.gz")
    T2w=("$datadir/sub-${participant}_ses-${session}_T2w.nii.gz")
    FLAIR=("$datadir/sub-${participant}_ses-${session}_FLAIR.nii.gz")
    MTI=("$datadir/sub-${participant}_ses-${session}_MTI.nii.gz")
else
    T1w=($(find BIDS -type f -name "*T1w.nii.gz" | sort ))
fi

d=0
t2=0
flair=0
mti=0
for test_T1w in ${T1w[@]}; do

    base0=${test_T1w##*/};base=${base0%_T1w*}
    check_done="$outputdir/compute/${base}.done"

    if [ ! -f $check_done ];then

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

        # If a T2 and/or a FLAIR exists
        if [ $d -gt 0 ]; then
            echo "KUL_T1T2FLAIR_ratio is starting"
            mkdir -p $outputdir/compute

            # for the T1w
            input=$test_T1w
            output=${test_T1w##*/}
            output=${output%%.*}
            echo " doing hd-bet and biascorrection on image $output"
            crop_x=0
            crop_z=0
            KUL_reorient_crop_hdbet_biascorrect_iso
            mask_T1W=$mask
            cp $iso_output $outputdir/${base}_T1w.nii.gz

            if [ $t2 -eq 1 ];then
                input=$test_T2w
                output=${test_T2w##*/}
                output=${output%%.*}
                crop_x=0
                crop_z=0
                echo " doing hd-bet and biascorrection of the T2w"
                KUL_reorient_crop_hdbet_biascorrect_iso
                td="T2w"
                echo " coregistering T2 to T1 and computing the ratio"
                KUL_register_computeratio
            fi

            if [ $flair -eq 1 ];then
                input=$test_FLAIR
                output=${test_FLAIR##*/}
                output=${output%%.*}
                crop_x=0
                crop_z=0
                echo " doing hd-bet and biascorrection of the FLAIR"
                KUL_reorient_crop_hdbet_biascorrect_iso
                td="FLAIR"
                echo " coregistering FLAIR to T1 and computing the ratio"
                KUL_register_computeratio
            fi

            if [ $mti -eq 1 ];then
                input=$test_MTI
                output=${test_MTI##*/}
                output=${output%%.*}
                crop_x=0
                crop_z=0
                echo " doing hd-bet of the MTI"
                KUL_MTI_reorient_crop_hdbet_iso
                td="MTI"
                echo " coregistering MTI to T1 and computing the MTC ratio"
                KUL_MTI_register_computeratio
            fi

            #rm -fr $outputdir/compute/${base}*.gz
            touch $check_done

            echo " done"

        else
            echo " Nothing to do here"
        fi

    else
        echo " $base already done"
    fi

done
