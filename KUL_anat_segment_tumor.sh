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

# Run resseg on a brain-mask-restricted input, with a sanity check against
# hallucination: resseg's model was never trained on hard-zero-masked
# images, and can occasionally predict a "cavity" entirely inside the
# zeroed-out (masked-away) region -- i.e. pure hallucination on background
# it has never seen anything like, unrelated to any real anatomy. Detected
# by checking what fraction of the predicted cavity falls outside the mask
# that built its own input; if that's the majority of the prediction, retry
# once with a larger dilation margin (more real tissue context around the
# mask boundary tends to avoid this), and just keep whichever result is
# sane. Globals used: $resseg_intermediate1, $resseginput, $ressegoutputdir,
# $ncpu, $resseg_seed, $verbose_level, $resseg_mask_dilation_npass,
# $resseg_mask_dilation_npass_fallback.
function KUL_resseg_masked_run {
    local base_mask="$1"     # undilated, native-space binary mask
    local masked_input="$2"  # output: masked T1 to feed resseg
    local cavity_out="$3"    # output: resseg's cavity segmentation
    local label="$4"         # for logging / log-file naming
    local dil_mask="$ressegoutputdir/${label}_maskdil.nii.gz"
    local npass n_found n_outside outside_frac

    for npass in $resseg_mask_dilation_npass $resseg_mask_dilation_npass_fallback; do
        maskfilter $base_mask dilate -npass $npass -nthreads $ncpu $dil_mask -force
        mrcalc $resseg_intermediate1/${resseginput}.nii.gz $dil_mask -mult $masked_input -force

        KUL_activate_conda_env resseg
            task_in="resseg -a 3 -s $resseg_seed -t $ressegoutputdir/${resseginput}_reg2mni.tfm \
                -o $cavity_out $masked_input"
            KUL_task_exec $verbose_level "resseg $label (npass=$npass)" "resseg_${label}"
        conda deactivate

        n_found=$(mrstats $cavity_out -output count -ignorezero)
        if [ "$n_found" -gt 0 ]; then
            n_outside=$(mrcalc $masked_input 0 -eq $cavity_out -mult - -quiet | mrstats - -output count -ignorezero -quiet)
            outside_frac=$(echo "scale=3; $n_outside / $n_found" | bc -l)
            kul_echo "  resseg $label (npass=$npass): $n_found voxels found, $n_outside outside the masked input (fraction: $outside_frac)"
            if (( $(echo "$outside_frac < 0.5" | bc -l) )); then
                break
            fi
            kul_echo "  resseg $label: hallucinated on masked-out background at npass=$npass"
        else
            kul_echo "  resseg $label (npass=$npass): nothing found"
        fi
        [ "$npass" -eq "$resseg_mask_dilation_npass_fallback" ] || kul_echo "  retrying resseg $label with a larger dilation margin"
    done
}

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
        kul_echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
        kul_echo " At least one is missing. Check the derivatives folder"
        #exit 1
    fi


    echo -e "\n\n"

}


function KUL_hd_glio_auto {
    
    # Segmentation of the tumor using HD-GLIO-AUTO
    
    # only run if not yet done

    if [ ! -f ${hdgliooutputdir}/output/segmentation.nii.gz ]; then

        # prepare the inputs
        mkdir -p $hdglioinputdir
        mkdir -p $hdgliooutputdir/output
        cp -f $cwd/$T1w $hdglioinputdir/T1.nii.gz
        cp -f $cwd/$cT1w $hdglioinputdir/CT1.nii.gz
        cp -f $cwd/$FLAIR $hdglioinputdir/FLAIR.nii.gz
        cp -f $cwd/$T2w $hdglioinputdir/T2.nii.gz
        
        # Run HD-GLIO-AUTO from the native (non-docker) install created by the
        # installer's env-hdglio section (KUL_Linux_setup).
        #
        # run.py invokes 'hd-bet' and 'hd_glio_predict' as subprocesses, by name,
        # so the hdglio env's bin must be on PATH -- and it must be *that* env's
        # hd-bet: HD-GLIO-AUTO calls 'hd-bet -device 0', which HD-BET 2.x (the
        # version in hd-bet-env, used elsewhere in KUL_NIS) rejects. The two are
        # installed separately and deliberately. run.py also needs FSL on PATH
        # for fslreorient2std/flirt/fslmaths.
        #
        # This used to fall back to bare 'hd-bet'/'hd_glio_predict' when the
        # local install was absent. Those are not on PATH in any standard
        # KUL_NIS setup, and the calls were not wrapped in KUL_task_exec, so the
        # step failed with no log, no warning and no segmentation -- and the run
        # continued to produce a lesion mask with the tumour missing from it.
        # Failing loudly here is the point of this block.
        # Search candidate software roots rather than trusting one variable:
        # SOFTWARE_ROOT is exported by the installer's bashrc block, which is
        # not guaranteed to have been sourced in the shell this runs from, and
        # /usr/local/KUL_apps is where installs lived before that move.
        local _hdglio_run="" _hdglio_bin="" _sw
        for _sw in "${SOFTWARE_ROOT:-}" /opt/kul_software /usr/local/KUL_apps; do
            [ -z "$_sw" ] && continue
            if [ -f "$_sw/src/HD-GLIO-AUTO/scripts/run.py" ] && [ -x "$_sw/miniforge3/envs/hdglio/bin/python" ]; then
                _hdglio_run="$_sw/src/HD-GLIO-AUTO/scripts/run.py"
                _hdglio_bin="$_sw/miniforge3/envs/hdglio/bin"
                break
            fi
        done

        if [ -z "$_hdglio_run" ]; then
            kul_echo "ERROR: HD-GLIO-AUTO is not installed."
            kul_echo "  looked for <root>/src/HD-GLIO-AUTO/scripts/run.py together with"
            kul_echo "  <root>/miniforge3/envs/hdglio/bin/python, under:"
            kul_echo "    ${SOFTWARE_ROOT:-(SOFTWARE_ROOT unset)}, /opt/kul_software, /usr/local/KUL_apps"
            kul_echo "  install it with the KUL_Linux_setup installer:"
            kul_echo "    https://github.com/Rad-dude/KUL_Linux_setup"
            kul_echo "    ./setup_environment.sh --only env-hdglio"
            kul_echo "  (a jenspetersen/hd-glio-auto docker image also exists, but this"
            kul_echo "   script does not drive it -- the native install is what is used.)"
            exit 1
        fi

        hdglio_type="native install ($_hdglio_run)"
        local _path_before="$PATH"
        export PATH="$_hdglio_bin:$PATH"
        task_in="$_hdglio_bin/python $_hdglio_run -i $hdglioinputdir -o $hdgliooutputdir/output -v -np"
        KUL_task_exec $verbose_level "HD-GLIO-AUTO using $hdglio_type" "hdglioauto"
        export PATH="$_path_before"

        # KUL_task_exec reports a non-zero exit, but the interesting failure is a
        # missing output: run.py writes segmentation.nii.gz several minutes before
        # it finishes, so a late crash can leave a usable-looking directory. Check
        # the file the rest of this script actually consumes.
        if [ ! -f ${hdgliooutputdir}/output/segmentation.nii.gz ]; then
            kul_echo "ERROR: HD-GLIO-AUTO produced no segmentation.nii.gz."
            kul_echo "  see ${KUL_LOG_DIR}/hdglioauto.error.log for what it did"
            exit 1
        fi

    else
        kul_echo "Already done HD-GLIO-AUTO"
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
        
        kul_echo "Running resseg"
        
        # prepare the inputs
        mkdir -p $resseginputdir1
        mkdir -p $resseginputdir2
        mkdir -p $ressegoutputdir

        # everything that isn't itself fed to resseg lives in intermediate1/
        # -- resseginputdir1's own top level holds only the final masked
        # image resseg run 1 actually takes as input
        resseg_intermediate1="$resseginputdir1/intermediate"
        mkdir -p $resseg_intermediate1

        cp $T1w $resseg_intermediate1/${resseginput}.nii.gz

        # fixed seed so the two resseg runs (and repeat invocations of this
        # script) are reproducible -- resseg's TTA (-a 3) is stochastic
        # without one, and the two runs' arbitration (STEP 5C, below) is only
        # meaningful if each run's result is stable from one invocation to
        # the next
        resseg_seed=42

        # both masks start at the same (tighter) margin, so the only
        # difference between the two runs is which mask (subject-specific
        # vs normal-population template) -- not how generously either is
        # dilated -- unless a run hallucinates, in which case
        # KUL_resseg_masked_run retries it at the larger fallback margin.
        resseg_mask_dilation_npass=5
        resseg_mask_dilation_npass_fallback=15

        # resseg-mni still runs once, up front -- its transform is used
        # internally by resseg itself for its own MNI-space inference
        # (unrelated to the brain-to-brain registration below), and is
        # shared by both runs
        KUL_activate_conda_env resseg
            task_in="resseg-mni -t $ressegoutputdir/${resseginput}_reg2mni.tfm \
                -r $ressegoutputdir/${resseginput}_reg2mni.nii.gz \
                $resseg_intermediate1/${resseginput}.nii.gz"
            KUL_task_exec $verbose_level "resseg running mni" "resseg_mni"
        conda deactivate

        # Build the NORMAL-BRAIN (Colin27) mask in native space for run 1.
        # A brain-to-brain rigid+affine registration (not resseg-mni's own
        # full-head-to-full-head one) gives a more accurate mapping in the
        # region that matters. Unlike run 2's patient-specific BET mask, a
        # template mask reflects where a normal brain would be, so it cannot
        # exclude a resection cavity the way this patient's own (possibly
        # post-surgical) brain extraction might.
        colin_t1="$HOME/.cache/torchio/mni_colin27_1998_nifti/colin27_t1_tal_lin.nii.gz"
        colin_mask="$HOME/.cache/torchio/mni_colin27_1998_nifti/colin27_t1_tal_lin_mask.nii.gz"
        mrcalc $colin_t1 $colin_mask -mult $ressegoutputdir/colin27_brain.nii.gz -force
        mrcalc $resseg_intermediate1/${resseginput}.nii.gz $kulderivativesdir/hdglio/output/mask.nii.gz -mult \
            $resseg_intermediate1/${resseginput}_brain.nii.gz -force

        antsRegistrationSyN.sh -d 3 \
            -f $ressegoutputdir/colin27_brain.nii.gz \
            -m $resseg_intermediate1/${resseginput}_brain.nii.gz \
            -t a -n $ncpu \
            -o $ressegoutputdir/native2colin_

        antsApplyTransforms -d 3 --float 1 \
            -i $colin_mask \
            -o $resseg_intermediate1/colin_mask_native.nii.gz \
            -r $resseg_intermediate1/${resseginput}.nii.gz \
            -t [$ressegoutputdir/native2colin_0GenericAffine.mat,1] -n NearestNeighbor

        # run resseg 1st time (Colin27-derived mask), with hallucination
        # retry
        KUL_resseg_masked_run \
            "$resseg_intermediate1/colin_mask_native.nii.gz" \
            "$resseginputdir1/${resseginput}_mnimasked.nii.gz" \
            "$ressegoutputdir/${resseginput}_cavity1.nii.gz" \
            "run1"

        # run resseg 2nd time (subject-specific BET mask), with
        # hallucination retry
        KUL_resseg_masked_run \
            "$kulderivativesdir/hdglio/output/mask.nii.gz" \
            "$resseginputdir2/${resseginput}.nii.gz" \
            "$ressegoutputdir/${resseginput}_cavity2.nii.gz" \
            "run2"


    else
        echo "Already done resseg"
    fi

}

function KUL_fast {
    
    # Segmentation of the image using FSL FAST
    
    kul_echo "Running FSL FAST"
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
        kul_echo "Running segmentation-only fastsufer"
        KUL_activate_conda_env fastsurfer_gpu
        task_in="$FASTSURFER_HOME/run_fastsurfer.sh \
            --sid $participant --sd $fastsurferoutputdir \
            --t1 $cwd/$T1w \
            --seg_only --py python --ignore_fs_version"
        KUL_task_exec $verbose_level "Running FastSurfer" "Fastsurfer"
        conda deactivate
    else
        kul_echo "Already run Fastsurfer"
    fi
}

# --- MAIN ---

# STEP 1 - Setup & Check to input data
# setup
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_anat_segment_tumor
globalresultsdir=$cwd/RESULTS/sub-$participant

hdglioinputdir="$kulderivativesdir/hdglio/input"
hdgliooutputdir="$kulderivativesdir/hdglio"
hdglio_segmentation=$hdgliooutputdir/output/segmentation.nii.gz
hdglio_output0=$hdgliooutputdir/hdglio_lesion_empty.nii.gz
hdglio_output1=$hdgliooutputdir/hdglio_lesion_perilesional_tissue.nii.gz
hdglio_output2=$hdgliooutputdir/hdglio_lesion_solid_tissue.nii.gz
hdglio_output3=$hdgliooutputdir/hdglio_lesion_total.nii.gz
local_output_hdglio1=$kulderivativesdir/sub-${participant}_hdglio_lesion_perilesional_tissue.nii.gz
local_output_hdglio2=$kulderivativesdir/sub-${participant}_hdglio_lesion_solid_tissue.nii.gz
local_output_hdglio3=$kulderivativesdir/sub-${participant}_hdglio_lesion_total.nii.gz
global_output_hdglio1=$globalresultsdir/Lesion/sub-${participant}_hdglio_lesion_perilesional_tissue.nii.gz
global_output_hdglio2=$globalresultsdir/Lesion/sub-${participant}_hdglio_lesion_solid_tissue.nii.gz
global_output_hdglio3=$globalresultsdir/Lesion/sub-${participant}_hdglio_lesion_total.nii.gz


input_resseg1=$kulderivativesdir/resseg/output/T1_cavity1.nii.gz
input_resseg2=$kulderivativesdir/resseg/output/T1_cavity2.nii.gz
local_output_resseg=$kulderivativesdir/sub-${participant}_resseg_cavity_only.nii.gz
global_output_resseg=$globalresultsdir/Lesion/sub-${participant}_resseg_cavity_only.nii.gz

fastsurferoutputdir="$kulderivativesdir/fastsurfer"
input_fastsurfer=$fastsurferoutputdir/$participant/mri/aparc.DKTatlas+aseg.deep.mgz
fastsurferoutput=$kulderivativesdir/sub-${participant}_fastsurfer_ventricles.nii.gz


global_output_full=$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz


if [ -f $globalresultsdir/Lesion/sub-${participant}_tumor_segment.png ] && [ $result -eq 0 ];then
    echo "Already done."
    exit
fi

# Check if fMRI and/or dwi data are present and/or to redo some processing
KUL_check_data


mkdir -p $kulderivativesdir
globalresultsdir=$cwd/RESULTS/sub-$participant
mkdir -p $globalresultsdir/Lesion
mkdir -p $globalresultsdir/Anat


if [ $result -eq 0 ]; then
    # Get the data from KUL_anat_register_rigid
    cp $cwd/BIDS/derivatives/KUL_compute/sub-$participant/KUL_anat_register_rigid/*.gz $globalresultsdir/Anat

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

    kul_echo "Running final segmentations"

    # compute some additional output

    # STEP 5A - HD-GLIO-AUTO
    hdglio_type_found=$(mrstats -output max $hdglio_segmentation)
    #echo $hdglio_type_found
    
    if [ $hdglio_type_found -eq 0 ];then

        kul_echo "hd-glio-auto did not find a lesion"
        # output an empty lesion mask
        mrcalc $hdglio_segmentation 1 -eq ${hdglio_output0}.nii.gz -force

    fi   

    if [ $hdglio_type_found -le 2 ];then

        kul_echo "hd-glio-auto found $hdglio_output1"
        #echo $hdglio_type_found
        cmd="mrcalc $hdglio_segmentation 1 -eq - | maskfilter - dilate -npass 10 -nthreads $ncpu - | \
        maskfilter - fill - -nthreads $ncpu | \
        maskfilter - erode ${hdglio_output1} -npass 10 -nthreads $ncpu -force"
        #echo $cmd
        eval $cmd 

        cp ${hdglio_output1} ${hdglio_output3}
        
        ln -sf $hdglio_output1 $local_output_hdglio1
        ln -sf $hdglio_output1 $local_output_hdglio3
        ln -sf $hdglio_output1 $global_output_hdglio1
        ln -sf $hdglio_output1 $global_output_hdglio3

    fi

    if [ $hdglio_type_found -eq 2 ];then

        kul_echo "hd-glio-auto found $hdglio_output2"
        #echo $hdglio_type_found
        cmd="mrcalc $hdglio_segmentation 2 -eq - | maskfilter - dilate -npass 10 -nthreads $ncpu - | \
        maskfilter - fill - -nthreads $ncpu | \
        maskfilter - erode ${hdglio_output2} -npass 10 -nthreads $ncpu -force"
        #echo $cmd
        eval $cmd

        mrcalc ${hdglio_output1} ${hdglio_output2} -add 0.9 -gt ${hdglio_output3} -force
        mrcalc ${hdglio_output3} ${hdglio_output2} -subtract 0.9 -gt ${hdglio_output1} -force
        
        ln -sf $hdglio_output1 $local_output_hdglio1
        ln -sf $hdglio_output2 $local_output_hdglio2
        ln -sf $hdglio_output3 $local_output_hdglio3
        ln -sf $hdglio_output1 $global_output_hdglio1
        ln -sf $hdglio_output2 $global_output_hdglio2
        ln -sf $hdglio_output3 $global_output_hdglio3

    fi


    # STEP 5B - FASTSURFER
    # get the CSF from fastsurfer 

    # regrid fastsurfer to T1w space, get L and R ventricle and add them
    mrgrid $input_fastsurfer regrid -template $T1w \
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
    cp $fastsurferoutput $globalresultsdir/Lesion/


    # STEP 5C - Correct RESSEG
    # resseg is run twice (run 1: raw full head; run 2: dilated-BET-masked) as
    # a lightweight consistency check. The two runs can legitimately disagree
    # -- run 1 in particular, having no brain mask, can hallucinate a "cavity"
    # on extracranial structures (e.g. the maxillary sinus, which is a dark,
    # fluid-filled shape not unlike a resection cavity). Arbitrate cheapest
    # signal first, only escalating when inconclusive:
    #   1. spatial overlap between the two runs
    #   2. distance to the HD-GLIO lesion (a genuine cavity sits directly
    #      adjacent to/overlapping the resection target; a stray extracranial
    #      detection does not) -- only when HD-GLIO actually found a lesion
    # If neither resolves it, we discard resseg (as before), but now say so
    # loudly instead of silently -- this mask feeds VBG for the rest of the
    # pipeline, so a silent wrong/missing cavity is worse than a visible one.

    # calculate the overlap between the 2 resseg runs
    mrcalc $input_resseg1 $input_resseg2 -add $kulderivativesdir/resseg/T1_cavity_combined.nii.gz -force
    cav_overlap=$(mrstats $kulderivativesdir/resseg/T1_cavity_combined.nii.gz -output max)

    resseg_keep=""
    resseg_use=0
    resseg_decision="none"

    if [ $cav_overlap -ge 2 ]; then
        # the 2 runs agree spatially -- keep #1, as before
        resseg_keep=$input_resseg1
        resseg_use=1
        resseg_decision="overlap"

    elif [ $hdglio_type_found -ne 0 ]; then
        kul_echo "resseg runs disagree (no overlap); arbitrating by distance to the HD-GLIO lesion"
        n1=$(mrstats $input_resseg1 -output count -ignorezero)
        n2=$(mrstats $input_resseg2 -output count -ignorezero)

        if [ "$n1" -gt 0 ] && [ "$n2" -gt 0 ]; then
            ImageMath 3 $kulderivativesdir/resseg/lesion_dist.nii.gz MaurerDistance $hdglio_output3
            dist1=$(mrstats $kulderivativesdir/resseg/lesion_dist.nii.gz -mask $input_resseg1 -output min)
            dist2=$(mrstats $kulderivativesdir/resseg/lesion_dist.nii.gz -mask $input_resseg2 -output min)
            kul_echo "  run 1 nearest distance to lesion: ${dist1}mm; run 2: ${dist2}mm"
            if (( $(echo "$dist1 < $dist2" | bc -l) )); then
                resseg_keep=$input_resseg1
            else
                resseg_keep=$input_resseg2
            fi
            resseg_use=1
            resseg_decision="distance"
        elif [ "$n1" -gt 0 ]; then
            resseg_keep=$input_resseg1; resseg_use=1; resseg_decision="distance-only-run1"
        elif [ "$n2" -gt 0 ]; then
            resseg_keep=$input_resseg2; resseg_use=1; resseg_decision="distance-only-run2"
        fi
    fi

    if [ $resseg_use -eq 0 ]; then
        kul_echo "WARNING: could not confidently arbitrate the resseg cavity candidates for sub-${participant} (decision path: $resseg_decision)."
        kul_echo "  Proceeding WITHOUT a resseg cavity contribution -- MANUAL REVIEW RECOMMENDED for the resection cavity."
    fi
    
    # STEP 5D - Depending on output compute
    if [ $hdglio_type_found -eq 0 ]; then
        # HD-GLIO did not find anything
        # did resseg find anything?
        if [ $resseg_use -eq 0 ]; then
            # Oeps, resseg did not find anything
            # we keep nothing
            kul_echo "Sorry, nothing found, we exit here. Perform a manual segmentation please."
            exit

        else 
            # now we keep resseg - ventricles
            mrcalc $resseg_keep $fastsurferoutput -sub 0.9 -gt \
                $kulderivativesdir/tmp_sub-${participant}_cavity_only.nii.gz -force
        
        fi    

    else
        # HD-GLIO did find perilesion and/or solid tissue 
        # did resseg find anything?
        if [ $resseg_use -eq 1 ]; then

            # now we keep hd-glio - resseg - ventricles
            mrcalc $resseg_keep ${hdglio_output3} -subtract $fastsurferoutput -subtract 0.9 -gt \
                $kulderivativesdir/tmp_sub-${participant}_cavity_only.nii.gz -force

        fi

    fi
    

    # clean and fill the cavity
    maskfilter -nthreads $ncpu -npass 5 $kulderivativesdir/tmp_sub-${participant}_cavity_only.nii.gz dilate - | \
        maskfilter -nthreads $ncpu - fill - | \
        maskfilter -nthreads $ncpu -npass 5 - erode $local_output_resseg -force
    rm -rf $kulderivativesdir/tmp_*.gz
    cp $local_output_resseg $global_output_resseg


    # compute a refined whole lesion + cavity
    if [ $resseg_use -eq 1 ]; then
        mrcalc $local_output_resseg \
            $hdglio_output3 -add \
            $kulderivativesdir/tmp_lesion_full.nii.gz -force

    else

        ln -sf $hdglio_output3 \
            $kulderivativesdir/tmp_lesion_full.nii.gz

    fi

    maskfilter $kulderivativesdir/tmp_lesion_full.nii.gz dilate - -npass 5 -nthreads $ncpu | \
    maskfilter - connect - -nthreads $ncpu | \
    maskfilter - erode -npass 5 $kulderivativesdir/sub-${participant}_lesion_and_cavity.nii.gz  -nthreads $ncpu -force
    rm -rf $kulderivativesdir/tmp_*.nii.gz

    cp $kulderivativesdir/sub-${participant}_lesion_and_cavity.nii.gz \
        $global_output_full

fi


# create a figure
rm -f $globalresultsdir/Lesion/tmp*.png

underlay=$globalresultsdir/Anat/FLAIR_reg2_T1w.nii.gz
underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF)}')


mrview_global_output_full=""
mrview_hdglio1_overlay=""
mrview_hdglio2_overlay=""
mrview_ventricles_overlay=""
if [ -f $global_output_full ]; then
    mrview_global_output_full="-overlay.load $global_output_full -overlay.opacity 0.4 -overlay.colour 255,255,0 -overlay.threshold_min 0.1"
fi
if [ -f $global_output_hdglio1 ]; then
    mrview_hdglio1_overlay="-overlay.load $global_output_hdglio1 -overlay.opacity 0.4 -overlay.colour 85,0,255 -overlay.threshold_min 0.1"
fi
if [ -f $global_output_hdglio2 ]; then
    mrview_hdglio2_overlay="-overlay.load $global_output_hdglio2 -overlay.opacity 0.4 -overlay.colour 255,0,0 -overlay.threshold_min 0.1"
fi
if [ -f $global_output_resseg ]; then
    mrview_resseg_overlay="-overlay.load $global_output_resseg -overlay.opacity 0.4 -overlay.colour 170,85,0 -overlay.threshold_min 0.1"
fi
if [ -f $fastsurferoutput ]; then
    mrview_ventricles_overlay="-overlay.load $fastsurferoutput -overlay.opacity 0.4 -overlay.colour 0,85,127 -overlay.threshold_min 0.1"
fi


if [ $result -eq 0 ]; then
    fig2f="-t 2"
else
    fig2f=""
fi

config_mrview=study_config/mrview_overlay_segment_tumor.txt
overlays="$mrview_hdglio1_overlay $mrview_hdglio2_overlay $mrview_ventricles_overlay $mrview_resseg_overlay"
echo $overlays > $config_mrview
KUL_mrview_figure.sh -p ${participant} \
    -u $underlay -o $config_mrview \
    -d $globalresultsdir/Lesion \
    -f tumor_segment \
    $fig2f \
    -v $verbose_level


kul_echo "Finished"
