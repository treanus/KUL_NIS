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
     -f:  also run fastsurfer (1=full, 2=segmentation with CC & stats, 3=fast segmentation without CC)
     -d:  only perform part of the workflow
            1 = bias-correction
            2 = rigid registration of the T2w and/or FLAIR to the T1w
            3 = spatially normalise (warp) the T1w to the MNI atlas
            4 = run fastsufer (req option -f)
            5 = run samseg (req option -m)
            6 = hd-bet the images
            7 = calibration
            8 = compute the T1w/T2w and/or T1w/FLAIR ratio
            9 = compute the MTR
            10 = warp the results to MNI
     -c:  do not run, but output files for the VSC HPC
     -v:  show output from commands

Documentation:

    This script computes a T1/T2, T1/FLAIR and MTC (magnetisation transfer contrast) ratio, using BIDS organised data. 
    A T1w is mandatory. 
    A T1w/T2w ratio is computed if a T2w is present.
    A T1w/FLAIR ratio is computed if a FLAIR is present.
    A MTR is computed if an MTI pair is available.
    It follows the rationale of Ganzetti et al. Frontiers in human neuroscience 2014 and D. Pareto et al. AJNR 2020.
    It also segments (MS or T1w-hypo/FLAIR-hyper) lesions using Freesurfer Samseg if T1w and FLAIR are present and option -m is given.
    It also calculates a FastSurfer parcellation if option -f is used.

References:
    @ Sarah Cappelle & Stefan Sunaert

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
hpc=0
deel=10

# Set required options
#p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:f:n:d:acmv" OPT; do

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
        d) #new workflow
			deel=$OPTARG
		;;
        m) #MS lesion segmentation
			ms=1
		;;
        c) #VSC HPC
			hpc=1
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
    local input=$1
    local output
    output=${input##*/}
    output=${output%%.*}
    bias_output=$outdir/tmp/${output}_iso_biascorrected.nii.gz
    if [ ! -f $bias_output ]; then 
        echo "  doing biascorrection on the ${td}"
        mrgrid $input regrid -voxel 1 $outdir/tmp/${output}_iso.nii.gz -force
        local bias_input=$outdir/tmp/${output}_iso.nii.gz        
        N4BiasFieldCorrection --verbose $ants_verbose \
        -d 3 \
        -i $bias_input \
        -o $bias_output
    else
        echo "  biascorrection of the ${td} already done"
    fi 
}

function KUL_MTI_reorient_crop_hdbet_iso {
    iso_output=$outdir/tmp/${output}_iso.nii.gz
    mrgrid $input regrid -voxel 1 $iso_output -force
}

# Rigidly register the input to the T1w
function KUL_rigid_register {
antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$outdir/tmp/${ants_type},$outdir/tmp/${newname}] \
    --interpolation BSpline \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$outdir/tmp/$ants_template,$outdir/tmp/$ants_source,1] \
    --transform Rigid[0.1] \
    --metric MI[$outdir/tmp/$ants_template,$outdir/tmp/$ants_source,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
}

# Register and compute the ratio
function KUL_reg2t1 {
    local base0=${test_T1w##*/}
    local base=${base0%_T1w*}
    ants_type="${base}_rigid_${td}_reg2t1_"
    ants_template="${base}_T1w_iso_biascorrected.nii.gz"
    ants_source="${base}_${td}_iso_biascorrected.nii.gz"
    newname="${base}_${td}_iso_biascorrected_reg2T1w.nii.gz"
    #finalname="${base}_${td}_reg2T1w.nii.gz"
    if [ ! -f $outdir/tmp/${newname} ]; then
        echo "  rigid registration of the ${td} to the T1w"
        KUL_rigid_register
    else
        echo "  skipping rigid registration of the ${td} to the T1w - already done"
    fi
}

function KUL_computeratio {
    if [ ! -f $outdir/${base}_ratio-T1${td}_calib-nonlin3.nii.gz ];then
        echo "  computing the ratio T1w/${td}"
        mrcalc \
            $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz \
            $outdir/tmp/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz -divide \
            $outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
            $outdir/${base}_ratio-T1${td}_calib-none.nii.gz -nthreads $ncpu -force
        mrcalc \
            $outdir/tmp/${base}_T1w_iso_biascorrected_calib-lin.nii.gz \
            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-lin_reg2T1w.nii.gz -divide \
            $outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
            $outdir/${base}_ratio-T1${td}_calib-lin.nii.gz -nthreads $ncpu -force
        mrcalc \
            $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin.nii.gz \
            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz -divide \
            $outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
            $outdir/${base}_ratio-T1${td}_calib-nonlin.nii.gz -nthreads $ncpu -force
        mrcalc \
            $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin2.nii.gz \
            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz -divide \
            $outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
            $outdir/${base}_ratio-T1${td}_calib-nonlin2.nii.gz -nthreads $ncpu -force
        mrcalc \
            $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin3.nii.gz \
            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-nonlin3_reg2T1w.nii.gz -divide \
            $outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz -multiply \
            $outdir/${base}_ratio-T1${td}_calib-nonlin3.nii.gz -nthreads $ncpu -force
    else
        echo "  the ratio T1w/${td} is already computed"
    fi
}

function KUL_apply_warp2mni {
    echo "  warping ${td} results to MNI"

    reference=$fix_im
    transform1="$outdir/warp2mni/${base}_T1w2MNI_1Warp.nii.gz"
    transform2="[$outdir/warp2mni/${base}_T1w2MNI_0GenericAffine.mat,0]"
    interpolation_type="LanczosWindowedSinc"

    if [[ ${td} == "T1w" ]];then
        local ext=""
    elif [[ ! ${td} == "MTI" ]];then
        local ext="_reg2T1w"
        input="$outdir/${base}_ratio-T1${td}_calib-none.nii.gz"
        output="$outdir/${base}_space-MNI_ratio-T1${td}_calib-none.nii.gz"
        KUL_antsApply_Transform_MNI
        input="$outdir/${base}_ratio-T1${td}_calib-lin.nii.gz"
        output="$outdir/${base}_space-MNI_ratio-T1${td}_calib-lin.nii.gz"
        KUL_antsApply_Transform_MNI
        input="$outdir/${base}_ratio-T1${td}_calib-nonlin.nii.gz"
        output="$outdir/${base}_space-MNI_ratio-T1${td}_calib-nonlin.nii.gz"
        KUL_antsApply_Transform_MNI
        input="$outdir/${base}_ratio-T1${td}_calib-nonlin2.nii.gz"
        output="$outdir/${base}_space-MNI_ratio-T1${td}_calib-nonlin2.nii.gz"
        KUL_antsApply_Transform_MNI
        input="$outdir/${base}_ratio-T1${td}_calib-nonlin3.nii.gz"
        output="$outdir/${base}_space-MNI_ratio-T1${td}_calib-nonlin3.nii.gz"
        KUL_antsApply_Transform_MNI
    elif [[ ${td} == "MTI" ]];then
        local ext="_reg2T1w"
        input="$outdir/${base}_ratio-MTC.nii.gz"
        output="$outdir/${base}_space-MNI_ratio-MTC.nii.gz"
        KUL_antsApply_Transform_MNI
    fi

    input="$outdir/${base}_${td}${ext}.nii.gz"
    output="$outdir/${base}_space-MNI_${td}${ext}.nii.gz"
    KUL_antsApply_Transform_MNI
}

function KUL_MTI_register_computeratio {
    base0=${test_T1w##*/}
    base=${base0%_T1w*}
    # convert the 4D MTI to single 3Ds
    input="$outdir/tmp/${base}_${td}_iso.nii.gz"
    S0="$outdir/tmp/${base}_${td}_iso_S0.nii.gz"
    Smt="$outdir/tmp/${base}_${td}_iso_Smt.nii.gz"
    mrconvert $input -coord 3 0 $S0 -force
    mrconvert $input -coord 3 1 $Smt -force
    # determine the registration
    ants_type="${base}_rigid_${td}_reg2t1_"
    ants_template="${base}_T1w_iso_biascorrected.nii.gz"
    ants_source="${base}_${td}_iso_Smt.nii.gz"
    newname="${base}_${td}_iso_Smt_reg2T1w.nii.gz"
    finalname="${base}_${td}_reg2T1w.nii.gz"
    KUL_rigid_register
    Smt="$outdir/tmp/$newname"
    # Now apply the coregistration to the 4D MTI 
    input=$S0
    output="$outdir/tmp/${base}_${td}_iso_S0_reg2T1w.nii.gz"
    S0=$output
    transform="$outdir/tmp/${base}_rigid_${td}_reg2t1_0GenericAffine.mat"
    reference="$outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz"
    KUL_antsApply_Transform
    # make a better mask
    mask=$outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz
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

    cp $outdir/tmp/$newname $outdir/$finalname
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
    T1w=($(find $cwd/BIDS -type f -name "*T1w.nii.gz" | sort ))
fi

hpc_counter=0
hpc_task=1

for test_T1w in ${T1w[@]}; do

    # defining the basename from the array elements
    base0=${test_T1w##*/};base=${base0%_T1w*}
    local_participant=${base%_ses*}
    local_session="ses-${base##*ses-}"
    outdir=$outputdir/$local_participant/$local_session
    check_done="$outdir/${base}.done"
    # resetting found files    
    d=0
    t2=0
    flair=0
    mti=0
    # read the level of completed workflow
    if [ ! -f $check_done ];then
        level_done=0
    else
        level_done=$(cat $check_done)
        if [ -z $level_done ];then
            echo "empty"
            level_done=10
        fi
    fi
    echo "level_done: $level_done"

    # only execute if workflow level in not yet reached
    if [ $level_done -lt $deel ];then

        # write files for the VSC
        if [ $hpc -eq 1 ];then
            
            mkdir -p VSC
            if [ ! -f VSC/VSC_commands.sh ];then
                echo "#!/bin/bash -e" > VSC/VSC_commands.sh
                echo "participant, session" > VSC/pbs_data_${hpc_task}.csv
            fi
            vsc_participant=${local_participant##*sub-}
            vsc_session=${local_session##*ses-}
            vsc_cmd="KUL_T1T2FLAIRMTR_ratio_new.sh -p $vsc_participant -s $vsc_session -m -f 2 -d 3"
            echo $vsc_cmd >> VSC/VSC_commands.sh
            echo "${vsc_participant},${vsc_session}" >> VSC/pbs_data_${hpc_task}.csv
            
            hpc_counter=$((hpc_counter+1))
            if [ $hpc_counter -ge 34 ];then
                hpc_task=$((hpc_task+1))
                hpc_counter=0
            fi

        else
        
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
                mkdir -p $outdir/tmp
                mkdir -p $outdir/histograms
                mkdir -p $outdir/masks
                mkdir -p $outdir/rois
                mkdir -p $outputdir/log
                kul_e2cl "KUL_T1T2FLAIR_ratio is processing $local_participant and session $local_session" ${outputdir}/${log}
                
                # PART 1 - make all images conform and bias-correct
                if [ $deel -ge 1 ]; then
                    # for the T1w
                    td="T1w"
                    KUL_iso_biascorrect $test_T1w
                    cp $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz $outdir/${base}_T1w.nii.gz

                    if [ $t2 -eq 1 ];then
                        td="T2w"
                        KUL_iso_biascorrect $test_T2w
                    fi

                    if [ $flair -eq 1 ];then
                        td="FLAIR"
                        KUL_iso_biascorrect $test_FLAIR
                    fi
                fi

                # PART 2 - register the T2w and/or FLAIR to the T1w
                if [ $deel -ge 2 ]; then
                    if [ $t2 -eq 1 ];then
                        td="T2w"
                        KUL_reg2t1
                        cp $outdir/tmp/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz $outdir/${base}_${td}_reg2T1w.nii.gz
                    fi

                    if [ $flair -eq 1 ];then
                        td="FLAIR"
                        KUL_reg2t1
                        cp $outdir/tmp/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz $outdir/${base}_${td}_reg2T1w.nii.gz
                    fi                    
                fi
                
                # PART 3 - Spatially normalise (warp) the T1w to the MNI atlas
                if [ $deel -ge 3 ]; then
                    fix_im="$kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a.nii"
                    ref_im="$kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a_with_neck.nii"
                    #input=${test_T1w##*/}
                    #input=${input%%.*}
                    #mov_im=$outdir/tmp/${input}_iso_biascorrected.nii.gz
                    mov_im=$outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz
                    output="$outdir/warp2mni/${base}_T1w2MNI_"
                    mkdir -p $outdir/warp2mni
                    if [ ! -f $outdir/warp2mni/${base}_T1w2MNI_Warped.nii.gz ]; then 
                        echo "  computing the spatial normalisation (warp) of the T1w to MNI (takes about 20 minutes)"
                        my_cmd="antsRegistrationSyN.sh -d 3 -f ${fix_im} -m ${mov_im} \
                        -o ${output} -n ${ncpu} -j 1 -t s $fs_silent"
                        eval $my_cmd
                    else
                        echo "  skipping MNI spatial normalisation, since it's already computed"
                    fi
                fi

                # PART 4 - run fastsurfer (optional)      
                if [ $deel -ge 4 ] && [ $fastsurf -gt 0 ]; then
                    mkdir -p $outdir/fs
                    if [ ! -f $outdir/fs/$base/mri/aparc.DKTatlas+aseg.deep.mgz ]; then
                        if [ $fastsurf -eq 1 ]; then
                            echo "  running full fastsufer"
                            my_cmd="$FASTSURFER_HOME/run_fastsurfer.sh \
                                --sid $base --sd $outdir/fs \
                                --t1 $test_T1w \
                                --fs_license $FS_LICENSE \
                                --vol_segstats --py python --parallel --threads $ncpu $fs_silent"
                        elif [ $fastsurf -eq 2 ]; then
                            echo "  running segmentation-with-CC & stats fastsufer"
                            my_cmd="$FASTSURFER_HOME/run_fastsurfer.sh \
                                --sid $base --sd $outdir/fs \
                                --t1 $test_T1w \
                                --fs_license $FS_LICENSE \
                                --seg_with_cc_only --vol_segstats --py python --ignore_fs_version --threads $ncpu $fs_silent"
                        elif [ $fastsurf -eq 3 ]; then
                            echo "  running segmentation-only fastsufer"
                            my_cmd="$FASTSURFER_HOME/run_fastsurfer.sh \
                                --sid $base --sd $outdir/fs \
                                --t1 $test_T1w \
                                --seg_only --py python --ignore_fs_version --threads $ncpu $fs_silent"
                        fi
                        eval $my_cmd
                    else
                        echo "  fastsurfer seems to have run already"
                    fi
                fi

                # PART 5 - MS lesion segmentation (optional)
                if [ $deel -ge 5 ] && [ $ms -eq 1 ];then
                    if [ $flair -eq 1 ];then
                        MSlesion="$outdir/rois/${base}_MSLesion.nii.gz"
                        if [ ! -f $MSlesion ];then
                            echo "  running samseg (takes about 20 minutes)"
                            T1w_iso="$outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz"
                            FLAIR_reg2T1w="$outdir/tmp/${base}_FLAIR_iso_biascorrected_reg2T1w.nii.gz"
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

                # PART 6 - hd-bet and background search
                if [ $deel -ge 6 ];then
                    # hd-bet brain extraction of the T1w
                    
                    if [ ! -f $outdir/tmp/${base}_T1w_iso_biascorrected_brain.nii.gz ]; then 
                        echo "  doing hd-bet on ${base}_T1w_iso_biascorrected.nii.gz"
                        my_cmd="hd-bet -i $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz \
                        -o $outdir/tmp/${base}_T1w_iso_biascorrected_brain.nii.gz $fs_silent"
                        eval $my_cmd
                        mv $outdir/tmp/${base}_T1w_iso_biascorrected_brain_mask.nii.gz $outdir/masks
                    else
                        echo "  skipping hd-bet on ${base}_T1w_iso_biascorrected.nii.gz - already done"
                    fi

                    # find the background by thresholding the T1w
                    background="$outdir/masks/${base}_T1w_iso_biascorrected_background.nii.gz"
                    if [ ! -f $background ]; then
                        max_T1w=$(mrstats $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz -output max)
                        echo "  thresholding background using 0.5% of max signal of the T1w - $max_T1w"
                        mrcalc $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz $max_T1w 0.005 -mul -lt \
                            $background -force 
                    else
                        echo "  thresholding background already done"
                    fi
                fi

                # PART 7 - Calibration
                check7="$outdir/histograms/${base}_T1w_iso_biascorrected_calib-nonlin3_histogram.csv"
                if [ $deel -ge 7 ] && [ ! -f $check7 ];then
                    # METHOD 1 - LINEAR histogram matching using eye/muscle tissue
                    echo "  performing linear histogram matching"
                    
                    # Warp the eye and muscle back to subject space
                    echo "  warping eye/muscle and skull/air back to subject space"
                    reference=$mov_im
                    transform1="$outdir/warp2mni/${base}_T1w2MNI_1InverseWarp.nii.gz"
                    transform2="[$outdir/warp2mni/${base}_T1w2MNI_0GenericAffine.mat,1]"
                    interpolation_type="NearestNeighbor"

                    input="$kul_main_dir/atlasses/Ganzetti2014/eyemask.nii"
                    output="$outdir/masks/${base}_eye.nii.gz"
                    KUL_antsApply_Transform_MNI

                    input="$kul_main_dir/atlasses/Ganzetti2014/tempmask.nii"
                    output="$outdir/masks/${base}_tempmuscle.nii.gz"
                    KUL_antsApply_Transform_MNI

                    # sum the masks
                    mrcalc $outdir/masks/${base}_eye.nii.gz $outdir/masks/${base}_tempmuscle.nii.gz -add \
                    $outdir/masks/${base}_eye_and_muscle.nii.gz -force
                    mrcalc $kul_main_dir/atlasses/Ganzetti2014/eyemask.nii $kul_main_dir/atlasses/Ganzetti2014/tempmask.nii -add \
                    $outdir/masks/mni_eye_and_muscle.nii.gz -force
                    
                    # do the LINEAR histogram matching using eye/muscle tissue
                    mrhistmatch \
                    -mask_input $outdir/masks/${base}_eye_and_muscle.nii.gz \
                    -mask_target $outdir/masks/mni_eye_and_muscle.nii.gz \
                    linear \
                    $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz \
                    $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a.nii \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-lin.nii.gz -force

                    mrhistogram -bin 100 -ignorezero \
                    $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz \
                    $outdir/histograms/${base}_T1w_iso_biascorrected_histogram.csv -force
                    mrhistogram -bin 100 -ignorezero \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-lin.nii.gz \
                    $outdir/histograms/${base}_T1w_iso_biascorrected_calib-lin_histogram.csv -force

                    if [ $t2 -eq 1 ];then
                        mrhistmatch \
                        -mask_input $outdir/masks/${base}_eye_and_muscle.nii.gz \
                        -mask_target $outdir/masks/mni_eye_and_muscle.nii.gz \
                        linear \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz \
                        $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t2_tal_nlin_sym_09a.nii \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_calib-lin_reg2T1w.nii.gz -force

                        mrhistogram -bin 100 -ignorezero \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_T2w_iso_biascorrected_reg2T1w_histogram.csv -force
                        mrhistogram -bin 100 -ignorezero \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_calib-lin_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_T2w_iso_biascorrected_calib-lin_reg2T1w_histogram.csv -force
                    fi

                    if [ $flair -eq 1 ];then
                        # Note: this is probably not a good calibration, since Ganzetti did not have a MNI-FLAIR 
                        mrhistmatch \
                        -mask_input $outdir/masks/${base}_eye_and_muscle.nii.gz \
                        -mask_target $outdir/masks/mni_eye_and_muscle.nii.gz \
                        linear \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_reg2T1w.nii.gz \
                        $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t2_tal_nlin_sym_09a.nii \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_calib-lin_reg2T1w.nii.gz -force

                        mrhistogram -bin 100 -ignorezero \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_FLAIR_iso_biascorrected_reg2T1w_histogram.csv -force
                        mrhistogram -bin 100 -ignorezero \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_calib-lin_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_FLAIR_iso_biascorrected_calib-lin_reg2T1w_histogram.csv -force
                    fi
                    
                    # METHOD 2 - the NON-LINEAR histogram matching using non-brain tissue
                    echo "  performing nonlinear histogram matching"
                    # Warp the brain_mask and its inverse to subject space
                    input="$kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated.nii"
                    output="$outdir/masks/${base}_MNI2subj_brainmask_mni_dilated.nii.gz"
                    KUL_antsApply_Transform_MNI
                    input="$kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated_inverse.nii"
                    output="$outdir/masks/${base}_MNI2subj_brainmask_mni_dilated_inverse.nii.gz"
                    KUL_antsApply_Transform_MNI

                    mrhistmatch \
                    -mask_input $outdir/masks/${base}_MNI2subj_brainmask_mni_dilated_inverse.nii.gz \
                    -mask_target $kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated_inverse.nii \
                    nonlinear \
                    $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz \
                    $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a.nii \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin.nii.gz -force

                    mrhistogram -bin 100 -ignorezero  \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin.nii.gz \
                    $outdir/histograms/${base}_T1w_iso_biascorrected_calib-nonlin_histogram.csv -force
                    
                    if [ $t2 -eq 1 ];then
                        mrhistmatch \
                        -mask_input $outdir/masks/${base}_MNI2subj_brainmask_mni_dilated_inverse.nii.gz \
                        -mask_target $kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated_inverse.nii \
                        nonlinear \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz \
                        $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t2_tal_nlin_sym_09a.nii \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz -force

                        mrhistogram -bin 100 -ignorezero  \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_T2w_iso_biascorrected_calib-nonlin_reg2T1w_histogram.csv -force
                    fi
                    
                    if [ $flair -eq 1 ];then
                        mrhistmatch \
                        -mask_input $outdir/masks/${base}_MNI2subj_brainmask_mni_dilated_inverse.nii.gz \
                        -mask_target $kul_main_dir/atlasses/Ganzetti2014/brainmask_mni_dilated_inverse.nii \
                        nonlinear \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_reg2T1w.nii.gz \
                        $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t2_tal_nlin_sym_09a.nii \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz -force

                        mrhistogram -bin 100 -ignorezero  \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_calib-nonlin_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_FLAIR_iso_biascorrected_calib-nonlin_reg2T1w_histogram.csv -force
                    fi


                    # Method 3 - Cappelle & Sunaert
                    echo "  performing second (Cappelle) nonlinear histogram matching"
                    mask1="$outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz"
                    mask2="$outdir/masks/${base}_T1w_iso_biascorrected_brain_inverted_mask.nii.gz"
                    mrcalc $mask1 0.1 -lt $mask2 -force

                    mask3="$outdir/masks/${base}_T1w_iso_biascorrected_brain_inverted_mask_nobackground.nii.gz"
                    mrcalc $mask1 0.1 -lt $background -sub 0.1 -gt $mask3 -force
                    
                    reference_histo_mask="$kul_main_dir/atlasses/Local/Cappelle2021/T1w_template_brain_mask_inverse_nobackground.nii.gz"
                    reference_histo_image="$kul_main_dir/atlasses/Local/Cappelle2021/T1w_template.nii.gz"

                    # do the scond NON-LINEAR histogram matching using non-brain tissue
                    mrhistmatch \
                    -mask_input $mask3 \
                    -mask_target $reference_histo_mask \
                    nonlinear \
                    $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz \
                    $reference_histo_image \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin2.nii.gz -force

                    mrhistogram -bin 100 -ignorezero  \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin2.nii.gz \
                    $outdir/histograms/${base}_T1w_iso_biascorrected_calib-nonlin2_histogram.csv -force
                    
                    if [ $t2 -eq 1 ];then
                        td="T2w"
                        # find the ventricles by thresholding the T2w
                        max_T2w=$(mrstats $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz -output max)
                        echo "  max signal of the T2w is $max_T2w"
                        ventricles="$outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w_ventricules.nii.gz"
                        mrcalc $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz $max_T2w 0.75 -mul -gt \
                            $ventricles -force 
                        skull_and_ventricules="$outdir/masks/${base}_T2w_iso_biascorrected_reg2T1w_skull_and_ventricules.nii.gz"
                        mrcalc $mask3 $ventricles -add 0.1 -gt $skull_and_ventricules -force

                        reference_histo_image="$kul_main_dir/atlasses/Local/Cappelle2021/${td}_template.nii.gz"                        
                        reference_histo_mask_nobackground="$kul_main_dir/atlasses/Local/Cappelle2021/T2wFLAIR_template_skull_and_ventricles_nobackground_mask.nii.gz"
                        # do the scond NON-LINEAR histogram matching using non-brain tissue
                        mrhistmatch \
                            -mask_input $skull_and_ventricules \
                            -mask_target $reference_histo_mask_nobackground \
                            nonlinear \
                            $outdir/tmp/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz \
                            $reference_histo_image \
                            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz -force
                        mrhistogram -bin 100 -ignorezero  \
                            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz \
                            $outdir/histograms/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w_histogram.csv -force
                    fi
                    
                    if [ $flair -eq 1 ];then
                        td="FLAIR"
                        # find the ventricles by thresholding the T2w
                        max_T2w=$(mrstats $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz -output max)
                        echo "  max signal of the T2w is $max_T2w"
                        ventricles="$outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w_ventricules.nii.gz"
                        mrcalc $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz $max_T2w 0.75 -mul -gt \
                            $ventricles -force 
                        skull_and_ventricules="$outdir/masks/${base}_T2w_iso_biascorrected_reg2T1w_skull_and_ventricules.nii.gz"
                        mrcalc $mask3 $ventricles -add 0.1 -gt $skull_and_ventricules -force

                        reference_histo_image="$kul_main_dir/atlasses/Local/Cappelle2021/${td}_template.nii.gz"                        
                        reference_histo_mask_nobackground="$kul_main_dir/atlasses/Local/Cappelle2021/T2wFLAIR_template_skull_and_ventricles_nobackground_mask.nii.gz"
                        # do the scond NON-LINEAR histogram matching using non-brain tissue
                        mrhistmatch \
                            -mask_input $skull_and_ventricules \
                            -mask_target $reference_histo_mask_nobackground \
                            nonlinear \
                            $outdir/tmp/${base}_${td}_iso_biascorrected_reg2T1w.nii.gz \
                            $reference_histo_image \
                            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz -force
                        mrhistogram -bin 100 -ignorezero  \
                            $outdir/tmp/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w.nii.gz \
                            $outdir/histograms/${base}_${td}_iso_biascorrected_calib-nonlin2_reg2T1w_histogram.csv -force
                    fi
                    
                    # Last method of calibration
                    # do it in the whole image except for the MS lesions
                    mask1="$outdir/masks/${base}_T1w_iso_biascorrected_brain_mask.nii.gz"
                    mask_subj="$outdir/masks/${base}_brain_mask_without_lesions.nii.gz"
                    MSlesion_dilated="$outdir/masks/${base}_MSlesion_dilated.nii.gz"
                    maskfilter $MSlesion dilate $MSlesion_dilated -force
                    mrcalc $mask1 $MSlesion_dilated -subtract $mask_subj -force
                    mask_target=$kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a_mask.nii

                    mrhistmatch \
                    -mask_input $mask_subj \
                    -mask_target $mask_target \
                    nonlinear \
                    $outdir/tmp/${base}_T1w_iso_biascorrected.nii.gz \
                    $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t1_tal_nlin_sym_09a.nii \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin3.nii.gz -force

                    mrhistogram -bin 100 -ignorezero  \
                    $outdir/tmp/${base}_T1w_iso_biascorrected_calib-nonlin3.nii.gz \
                    $outdir/histograms/${base}_T1w_iso_biascorrected_calib-nonlin3_histogram.csv -force

                    if [ $t2 -eq 1 ];then
                        mrhistmatch \
                        -mask_input $mask_subj \
                        -mask_target $mask_target \
                        nonlinear \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_reg2T1w.nii.gz \
                        $kul_main_dir/atlasses/Ganzetti2014/mni_icbm152_t2_tal_nlin_sym_09a.nii \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_calib-nonlin3_reg2T1w.nii.gz -force

                        mrhistogram -bin 100 -ignorezero  \
                        $outdir/tmp/${base}_T2w_iso_biascorrected_calib-nonlin3_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_T2w_iso_biascorrected_calib-nonlin3_histogram.csv -force
                    fi

                    if [ $flair -eq 1 ];then
                        mrhistmatch \
                        -mask_input $mask_subj \
                        -mask_target $mask_target \
                        nonlinear \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_reg2T1w.nii.gz \
                        $kul_main_dir/atlasses/Local/Cappelle2021/Winkler2009_GG-366-FLAIR_1.0mm_adapted.nii.gz \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_calib-nonlin3_reg2T1w.nii.gz -force

                        mrhistogram -bin 100 -ignorezero  \
                        $outdir/tmp/${base}_FLAIR_iso_biascorrected_calib-nonlin3_reg2T1w.nii.gz \
                        $outdir/histograms/${base}_FLAIR_iso_biascorrected_calib-nonlin3_histogram.csv -force
                    fi


                fi


                # PART 8 - compute the ratio
                if [ $deel -ge 8 ];then
                    if [ $t2 -eq 1 ];then
                        td="T2w"
                        KUL_computeratio
                    fi

                    if [ $flair -eq 1 ];then
                        td="FLAIR"
                        KUL_computeratio
                    fi
                fi

                # PART 9 - MTI and MTR
                test_done_MTI="$outdir/${base}_ratio-MTC_MNI.nii.gz"
                if [ $deel -ge 9 ] && [ $mti -eq 1 ] && [ ! -f $test_done_MTI ];then
                    input=$test_MTI
                    output=${test_MTI##*/}
                    output=${output%%.*}
                    echo "  doing hd-bet of the MTI"
                    KUL_MTI_reorient_crop_hdbet_iso
                    td="MTI"
                    echo "  coregistering MTI to T1 and computing the MTC ratio"
                    KUL_MTI_register_computeratio
                fi

                # PART 10  - warping results to mni too
                if [ $deel -ge 10 ];then
                    td="T1w"
                    KUL_apply_warp2mni
                    if [ $t2 -eq 1 ];then
                        td="T2w"
                        KUL_apply_warp2mni
                    fi

                    if [ $flair -eq 1 ];then
                        td="FLAIR"
                        KUL_apply_warp2mni
                    fi
                    
                    if [ $mti -eq 1 ];then
                        td="MTI"
                        KUL_apply_warp2mni
                    fi
                fi              

                #rm -fr $outdir/tmp/${base}*.gz
                #cd T1T2FLAIRMTR_ratio/tmp/
                #rm -fr *std.nii.gz *iso.nii.gz *ted.nii.gz *reg2T1w.nii.gz *mask.nii.gz *eye.nii.gz *muscle.nii.gz *.csv *MTI* *MNI.nii.gz *verse.nii.gz *inv.nii.gz *brain.nii.gz
                #rm -fr *std.nii.gz *iso.nii.gz *ted.nii.gz *reg2T1w.nii.gz *eye.nii.gz *muscle.nii.gz *.csv *MTI* *MNI.nii.gz *verse.nii.gz *inv.nii.gz *brain.nii.gz
                #cd ../..

                echo ${deel} > $check_done

                echo "  done"
             
            else
                echo "  Nothing to do here"
            fi
        fi 

    else
        echo "  $base already done"
    fi

done
