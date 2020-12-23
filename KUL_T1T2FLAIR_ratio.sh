#!/bin/bash -e
# Sarah Cappelle & Stefan Sunaert
# 22/12/2020
# This script is the first part of Sarah's Study1
# 
# This scripts follows the rationale of D. Pareto et al. AJNR 2020
# Starting from 3D-T1w, 3D-FLAIR and 2D-T2w scans we compute:
#  create masked brain images using HD-BET
#  bias correct the images using N4biascorrect from ANTs
#  ANTs rigid coregister and reslice all images to the 3D-T1w (in isotropic 1 mm space)
#  compute a T1FLAIR_ratio and a T1T2_ratio

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

	 -a:  automatic mode (just work on all images in the BIDS folder)
	 -v:  show output from mrtrix commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
auto=0 # default if option -s is not given
silent=1 # default if option -v is not given
outputdir="T1T2FLAIR_ratio"

# Set required options
#p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:av" OPT; do

		case $OPT in
		a) #automatic mode
			auto=1
		;;
		p) #participant
			participant=$OPTARG
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

ants_verbose=0
# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then

	export MRTRIX_QUIET=1
    ants_verbose=1

fi

# --- FUNCTIONS ---
 

# Function to do biascorrection
function KUL_ants_biascorrect {
N4BiasFieldCorrection --verbose $ants_verbose \
    -d 3 \
    -i $bias_input \
    -o $bias_output
}

# Function to do step 1
function KUL_reorient_crop_hdbet_biascorrect {
    fslreorient2std $input $outputdir/compute/std_$output
    mrgrid $outputdir/compute/std_$output crop -axis 0 $crop_x,$crop_x -axis 2 $crop_z,0 \
            $outputdir/compute/cropped_std_$output
    hd-bet -i $outputdir/compute/cropped_std_$output -o $outputdir/compute/brain_cropped_std_$output
    bias_input=$outputdir/compute/brain_cropped_std_$output
    bias_output=$outputdir/compute/biascorrected_brain_cropped_std_$output
    iso_output=$outputdir/compute/iso_biascorrected_brain_cropped_std_$output
    KUL_ants_biascorrect
    mrgrid $bias_output regrid -voxel 1 $iso_output
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
}


# --- MAIN ---
printf "\n\n\n"

# here we give the data
if [ $auto -eq 0 ]; then
    datadir="$cwd/BIDS/sub-${participant}/anat"
    T1w=("$datadir/sub-${participant}_T1w.nii.gz")
    T2w=("$datadir/sub-${participant}_T2w.nii.gz")
    FLAIR=("$datadir/sub-${participant}_FLAIR.nii.gz")
else
    T1w=($(find BIDS -type f -name "*T1w.nii.gz"))
fi

#n_T1w=${#T1w[@]}
#echo "there are $n_T1w T1w images"
d=0
t2=0
flair=0
for test_T1w in ${T1w[@]}; do

    # Test whether T2 and/or FLAIR also exist
    #echo "notably: $test_T1w"
    test_T2w="${test_T1w%_T1w*}_T2w.nii.gz"
    #echo $test_T2w
    if [ -f $test_T2w ];then
        #echo "The T2 exists"
        d=$((d+1))
        t2=1
    fi
    test_FLAIR="${test_T1w%_T1w*}_FLAIR.nii.gz"
    #echo $test_FLAIR
    if [ -f $test_FLAIR ];then
        #echo "The FLAIR exists"
        d=$((d+1))
        flair=1
    fi


    # If a T2 and/or a FLAIR exists
    if [ $d -gt 0 ]; then
        #echo "something to do here"
        mkdir -p $outputdir/compute

        # for the T1w
        input=$test_T1w
        output=${test_T1w##*/}
        crop_x=24
        crop_z=48
        KUL_reorient_crop_hdbet_biascorrect

        if [ $t2 -eq 1 ];then
            input=$test_T2w
            output=${test_T2w##*/}
            crop_x=0
            crop_z=0
            KUL_reorient_crop_hdbet_biascorrect

            base0=${test_T1w##*/}
            base=${base0%_T1w*}
            ants_type="rigid_T2w_reg2t1"
            newname="1iso_biascorrected_brain_cropped_std_${base}_T2w_reg2T1w.nii.gz"
            ants_template="iso_biascorrected_brain_cropped_std_${base}_T1w.nii.gz"
            ants_source="iso_biascorrected_brain_cropped_std_${base}_T2w.nii.gz"
            KUL_rigid_register

            mrcalc $outputdir/compute/$ants_template $outputdir/compute/$newname -divide  $outputdir/${base}-T1T2ratio.nii.gz

        fi

        if [ $flair -eq 1 ];then
            input=$test_FLAIR
            output=${test_FLAIR##*/}
            crop_x=0
            crop_z=0
            KUL_reorient_crop_hdbet_biascorrect

            base0=${test_T1w##*/}
            base=${base0%_T1w*}
            ants_type="rigid_FLAIR_reg2t1"
            newname="1iso_biascorrected_brain_cropped_std_${base}_FLAIR_reg2T1w.nii.gz"
            ants_template="iso_biascorrected_brain_cropped_std_${base}_T1w.nii.gz"
            ants_source="iso_biascorrected_brain_cropped_std_${base}_FLAIR.nii.gz"
            KUL_rigid_register

            mrcalc $outputdir/compute/$ants_template $outputdir/compute/$newname -divide  $outputdir/${base}-T1FLAIRratio.nii.gz

        fi

    else
        echo "Nothing to do here"
    fi

done
