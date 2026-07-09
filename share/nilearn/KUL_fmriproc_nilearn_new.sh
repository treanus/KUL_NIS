#!/bin/bash
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@uzleuven.be
#
# v0.1 - dd 19/01/2019 - jurassic version
version="v2.0-nilearn - dd 08/07/2026"
# Nilearn port: first-level GLM done in Python/Nilearn instead of MATLAB/SPM12.
# Requires: python3 with nilearn, nibabel, numpy, pandas (no MATLAB, no SPM).

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions


# To Do
#  - This script should be able to do the following:
# Decide what kind of fMRI data you have in your BIDS dir
# carry out a single subject analysis of all this BOLD data automatically
# for rs-fMRI ICA, time series clustering (SLIC/NCUT), atlas based roi to roi analysis and connectome generation should be automated
# we will probably need to use some matlab and python functions and commands here and there
# For tb-fMRI we can do FEAT +/- all the previous
# Include an option for high resolution rendering in the three orthogonal planes and mosaic generation in TIFF/JPEG (slicer)
# For mTE-rsfMRI, after MEICA, run same as single echo rs-fMRI, check out tedana also

# assuming most of the preproc is taken care of by fmriprep


# -----------------------------------  MAIN  ---------------------------------------------


# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs an automated task based fMRI Nilearn GLM analysis assuming:
    - a 30 seconds REST followed by 30 seconds TASK epochs
    - having run fmriprep with aroma

Note: requires python3 with nilearn, nibabel, numpy, pandas (no matlab/spm)

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001

Required arguments:

     -p:  participant


Optional arguments:

     -s:  session
     -S:  smoothing FWHM in mm for SUSAN (default: adaptive = mean voxel size)
            e.g. -S 5 for 5mm FWHM, -S 6 for 6mm FWHM
     -P:  FWE-corrected p-value for Bizzi thresholding (default: 0.01)
            e.g. -P 0.01, -P 0.005, -P 0.001
     -v:  verbose (0=silent, 1=normal, 2=verbose; default=1)


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
verbose_level=1
smooth_fwhm=0  # 0 = adaptive (mean voxel size)
pfwe=0.01


# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:S:P:v:" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            participant=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        S) #smoothing FWHM
            smooth_fwhm=$OPTARG
        ;;
        P) #FWE p-value threshold
            pfwe=$OPTARG
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

function KUL_tsv_filter {
    awk -F'\t' -v cols="trans_x,trans_y,trans_z,rot_x,rot_y,rot_z,a_comp_cor_00,a_comp_cor_01,a_comp_cor_02,a_comp_cor_03,a_comp_cor_04" '
    BEGIN {
        split(cols, col_arr, ",")
        for (i in col_arr) col_set[col_arr[i]] = 1
    }
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i in col_set) {
                col_idx[i] = 1
            }
        }
    }
    NR>1 {
        out = ""
        for (i=1; i<=NF; i++) {
            if (i in col_idx) {
                if (out == "") out = $i
                else out = out "\t" $i
            }
        }
        print out
    }' $filter_input > $filter_output
}

function KUL_compute_nilearn {
    # Nilearn first-level GLM (replaces KUL_compute_SPM_matlab + KUL_threshold_SPM_Bizzi).
    # Caller must set: nl_bolds[] nl_confounds_list[] nl_tr nl_mask nl_taskname fmrifile spm_type
    if [ $spm_type -eq 1 ]; then
        echo "   Nilearn analysis without confounds"
        fmriresults="$computedir/RESULTS/stats_$fmrifile"
        nl_hrf_deriv=0
        use_confounds=0
        global_result=${globalresultsdir}/afMRI_${fmrifile}.nii
        wc_suffix=""
    else
        echo "   Nilearn analysis with confounds"
        fmriresults="$computedir/RESULTS/stats_${fmrifile}_wc"
        nl_hrf_deriv=1
        use_confounds=1
        global_result=${globalresultsdir}/afMRI_${fmrifile}_wc.nii
        wc_suffix="_wc"
    fi

    # clean a possible old result
    rm -rf "$fmriresults"
    mkdir -p "$fmriresults"

    # assemble the (possibly multi-run) --bold and --confounds argument lists
    local bold_args=""
    for b in "${nl_bolds[@]}"; do bold_args="$bold_args \"$b\""; done
    local conf_args=""
    if [ $use_confounds -eq 1 ]; then
        conf_args="--confounds"
        for c in "${nl_confounds_list[@]}"; do conf_args="$conf_args \"$c\""; done
    fi

    cmd="$python_exe \"$nilearn_glm_script\" \
        --bold $bold_args \
        --mask \"$nl_mask\" \
        --tr $nl_tr \
        --task-name \"$nl_taskname\" \
        --output-dir \"$fmriresults\" \
        --rest-dur 30 --task-dur 30 --first rest \
        --hrf-derivative $nl_hrf_deriv \
        --pfwe $pfwe --cluster-k 50 \
        $conf_args"
    echo "  $cmd"
    eval $cmd $str_silent_SPM

    # Nilearn output is in MNI space; warp back to T1w space for display.
    local mni_to_t1w=$(find_first_match "${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_*from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5" "MNI-to-T1w transform")
    local find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
    reference=${find_T1w[0]}
    transform="$mni_to_t1w"

    # unthresholded stat map
    if [ -f "$fmriresults/spmT_0001.nii" ]; then
        input="$fmriresults/spmT_0001.nii"
        output="$global_result"
        KUL_antsApply_Transform
    fi

    # thresholded maps (same tags as the old Bizzi step, so downstream is unchanged)
    pfwe_tag=$(echo "$pfwe" | sed 's/0\.//' | sed 's/0*$//')
    fwe_tag="FWE${pfwe_tag}_k50"
    for thresh_tag in "p001unc_k50" "${fwe_tag}"; do
        thresh_nii="$fmriresults/spmT_0001_${thresh_tag}.nii"
        if [ -f "$thresh_nii" ]; then
            input="$thresh_nii"
            output="${globalresultsdir}/afMRI_${fmrifile}${wc_suffix}_${thresh_tag}.nii"
            transform="$mni_to_t1w"
            KUL_antsApply_Transform
        fi
    done
}

# MAIN --------------------------------------------------------------
python_exe=$(which python3)
# the Nilearn GLM worker lives next to this script (share/nilearn/)
nilearn_glm_script="$kul_main_dir/share/nilearn/KUL_nilearn_glm.py"

if [ $KUL_DEBUG -gt 0 ]; then
    echo "python3 lives at $python_exe"
    echo "nilearn glm script: $nilearn_glm_script"
fi

if [[ -z "$python_exe" ]]; then
    echo "python3 is required but not found on path. Exitting"
    exit 1
fi

if [ ! -f "$nilearn_glm_script" ]; then
    echo "Nilearn GLM worker not found at $nilearn_glm_script. Exitting"
    exit 1
fi

if ! "$python_exe" -c "import nilearn, nibabel, numpy, pandas" 2>/dev/null ; then
    echo "python3 is missing required packages. Install with:" >&2
    echo "    pip install nilearn nibabel numpy pandas" >&2
    exit 1
fi


KUL_check_participant


#  setup variables
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
computedir="$kulderivativesdir/sub-$participant/SPM"
fmridatadir="$computedir/fmridata"
scriptsdir="$computedir/scripts"
confoundsdir="$computedir/confounds"
fmriprepdir="${cwd}/fmriprep/sub-$participant"
globalresultsdir=$cwd/RESULTS/sub-$participant/SPM


if [ $KUL_DEBUG -gt 0 ]; then 
    echo "kulderivativesdir: $kulderivativesdir"
    echo "fmridatadir: $fmridatadir"
    echo "fmriprepdir: $fmriprepdir"
    echo "globalresultsdir: $globalresultsdir"
fi

mkdir -p $fmridatadir
mkdir -p $scriptsdir
mkdir -p $confoundsdir
mkdir -p $computedir/RESULTS
mkdir -p $globalresultsdir

# Use MNI-space fmriprep output so SPM stats are in MNI and the MNI→T1w
# warp-back at the end produces correctly aligned native-space results.
fmriprep_output_type="_space-MNI152NLin2009cAsym"


if [ $verbose_level -lt 2 ] ; then
    mkdir -p "KUL_LOG/$script"
    str_silent_SPM=" >> KUL_LOG/$script/sub-${participant}_nilearn.log"
fi

if [ ! -f KUL_LOG/sub-${participant}_SPM.done ]; then

    # find the output of fmriprep
    fmriprep_match=($(find $fmriprepdir/func \
        -name "*${fmriprep_output_type}*_desc-preproc_bold.nii.gz" \
        -type f))
    echo " fmriprep found the following files: "
    for match in ${fmriprep_match[@]}; do
        echo "  $match"
    done    

    # Extract unique task names from fmriprep_match
    tasks=()
    for match in "${fmriprep_match[@]}"; do
        task=$(basename "$match" | sed -E 's/.*_task-([A-Za-z0-9]+).*_desc-preproc_bold\.nii\.gz/\1/')
        tasks+=("$task")
    done
    unique_tasks=($(printf "%s\n" "${tasks[@]}" | sort -u))
    echo "Unique tasks: ${unique_tasks[@]}"

    for task in "${unique_tasks[@]}"; do

        if [[ ! "$task" = *"rest"* ]]; then

            task_files=()
            for match in "${fmriprep_match[@]}"; do
                # Extract the task name from the filename
                match_task=$(basename "$match" | sed -E 's/.*_task-([A-Za-z0-9]+).*/\1/')
                if [[ "$match_task" == "$task" ]]; then
                    task_files+=("$match")
                fi
            done
            echo "Files for task $task:"
            printf "  %s\n" "${task_files[@]}"
            # Now you can use "${task_files[@]}" for further processing

            tf_taskname=()
            tf_TR=()
            tf_mean_voxel_size=()
            tf_bt=()
            tf_smooth=()
            tf_confounds=()
            tf_mask=()
            for task_file in "${task_files[@]}"; do     
                echo " Processing file: $task_file"
                
                # get the task name
                taskname=$(basename "$task_file" | sed -E 's/.*_task-([A-Za-z0-9]+(_run-[0-9]+)?).*_desc-preproc_bold\.nii\.gz/\1/')
                tf_taskname+=("$taskname")
                echo "  task name: $taskname"
                
                # determine the TR
                TR=($(mrinfo $task_file -spacing | awk '{print $(NF)}'))
                echo "  repetition time (TR): $TR"
                tf_TR+=("$TR")

                # determine the voxel size
                spacing=($(mrinfo $task_file -spacing))
                echo "  the voxel size: ${spacing[0]} + ${spacing[1]} + ${spacing[2]}"
                mean=$(echo "(${spacing[0]} + ${spacing[1]} + ${spacing[2]}) / 3" | bc -l)
                tf_mean_voxel_size+=("$mean")
                echo "  mean voxel size: $mean"

                # determine SUSAN sigma: fixed FWHM if -S given, else adaptive (mean voxel size)
                if (( $(echo "$smooth_fwhm > 0" | bc -l) )); then
                    sigma=$(echo "$smooth_fwhm / 2.3548" | bc -l)
                    echo "  SUSAN sigma: ${sigma}mm (fixed FWHM=${smooth_fwhm}mm)"
                else
                    sigma=$mean
                    echo "  SUSAN sigma: ${sigma}mm (adaptive = mean voxel size)"
                fi

                # determine the brightness threshold for susan smoothing
                mask=$(dirname ${task_file})/$(basename ${task_file} "-preproc_bold.nii.gz")-brain_mask.nii.gz
                tf_mask+=("$mask")   # keep per-run mask for the Nilearn GLM

                # ---- brain-mask the preproc BOLD BEFORE smoothing --------------------
                # fmriprep (>=23.2) leaves the MNI/func output UNMASKED. With a tight
                # FOV the sinc/Lanczos resampling extrapolates a bright halo past the
                # superior brain edge; SUSAN would smear that halo into superior cortex
                # and it would end up in the SPM betas. Zero it out first (mask -> smooth
                # -> mask), and mask the SUSAN reference too so edge voxels aren't pulled
                # toward the (now zero) exterior.
                masked_file="$fmridatadir/$(basename ${task_file} ".nii.gz")_masked.nii.gz"
                if [ ! -f "$masked_file" ]; then
                    echo "  Masking $task_file with $mask -> $masked_file"
                    fslmaths "$task_file" -mas "$mask" "$masked_file"
                fi
                # ----------------------------------------------------------------------

                p50=$(fslstats ${masked_file} -k $mask -p 50)
                bt=$(echo "$p50 * 0.66666" | bc -l)
                tf_bt+=("$bt")
                echo "  brightness threshold (0.6666 * median): $bt"

                # run susan smoothing if not done yet
                smooth_file="$fmridatadir/$(basename ${task_file} "-preproc_bold.nii.gz")_smooth.nii"
                tf_smooth+=("$smooth_file")
                # find the boldref file (MNI-space boldref has no desc- prefix)
                boldref=$(dirname ${task_file})/$(basename ${task_file} "_desc-preproc_bold.nii.gz")_boldref.nii.gz
                echo "  boldref: $boldref"
                # mask the boldref (USAN reference) with the same brain mask
                boldref_masked="$fmridatadir/$(basename ${boldref} ".nii.gz")_masked.nii.gz"
                if [ ! -f "$boldref_masked" ]; then
                    fslmaths "$boldref" -mas "$mask" "$boldref_masked"
                fi
                if [ ! -f $smooth_file ]; then
                    echo "  Smoothing (masked) $masked_file to $smooth_file"
                    cmd="susan ${masked_file} $bt ${sigma} 3 1 1 \
                        $boldref_masked $bt \
                        $smooth_file"
                    echo "  $cmd"
                    eval $cmd
                    # clip any smoothing that bled just past the brain edge
                    echo "  Re-masking smoothed output with $mask"
                    fslmaths "$smooth_file" -mas "$mask" "$smooth_file"
                    gunzip -f $smooth_file.gz
                else
                    echo "  Smoothed file already exists: $smooth_file"
                fi

                # filter confounds
                # Confounds TSV never carries a _space- entity — find it by task+run key
                _task_run_key=$(basename "$task_file" | grep -oE 'task-[A-Za-z0-9]+(_run-[0-9]+)?')
                filter_input=$(find $(dirname ${task_file}) -name "*${_task_run_key}_desc-confounds_timeseries.tsv" | head -1)
                confounds_file="$confoundsdir/${_task_run_key}_confounds.txt"
                tf_confounds+=("$confounds_file")
                run_number=$(basename "$task_file" | sed -E 's/.*_run-0*([0-9]+).*/\1/')
                if [ -n "$run_number" ]; then
                    fmrifile="${task}_run-${run_number}"
                else
                    fmrifile="${task}"
                fi
                filter_output=$confounds_file
                echo "  filtering confounds from $filter_input to $filter_output"
                KUL_tsv_filter

            done

            unique_TRs=($(printf "%s\n" "${tf_TR[@]}" | sort -u))
            echo "  Unique TRs for task $task: ${unique_TRs[@]}"

            if [ ${#unique_TRs[@]} -gt 1 ]; then
                echo "  Warning: Multiple TRs found for task $task. Please check the data."
            fi 

            # Run Nilearn GLM for each run separately
            for i in "${!task_files[@]}"; do
                run_taskname="${tf_taskname[$i]}"
                fmrifile=$run_taskname
                nl_bolds=("${tf_smooth[$i]}")
                nl_confounds_list=("${tf_confounds[$i]}")
                nl_tr="${tf_TR[$i]}"
                nl_mask="${tf_mask[$i]}"
                nl_taskname="$run_taskname"
                echo "  Running Nilearn GLM for $run_taskname with TR=$nl_tr"

                # Run with and without confounds (threshold + warp-back happen inside)
                spm_types=(1 2) # 1: without confounds, 2: with confounds
                for spm_type in "${spm_types[@]}"; do
                    KUL_compute_nilearn
                done
            done

            # Aggregate all runs of this task together (fixed-effects in Nilearn).
            # Nilearn handles any number of runs, so no per-count templates needed.
            if [ ${#task_files[@]} -gt 1 ]; then
                fmrifile=$task
                nl_bolds=("${tf_smooth[@]}")
                nl_confounds_list=("${tf_confounds[@]}")
                nl_tr="${tf_TR[0]}"       # assumes a common TR across runs (as before)
                nl_mask="${tf_mask[0]}"
                nl_taskname="$task"
                echo "  Running aggregate Nilearn GLM for task $task with ${#task_files[@]} runs"
                spm_types=(1 2)
                for spm_type in "${spm_types[@]}"; do
                    KUL_compute_nilearn
                done
            fi
        fi
    done

    # cleanup
    #rm -rf $fmridatadir

    touch KUL_LOG/sub-${participant}_SPM.done
    echo "Done computing Nilearn GLM"
else
    echo "Nilearn analysis already done"
fi
