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

Requires the conda env named by \$KUL_PYFMRI_ENV (default 'pyfMRI'; nilearn,
nibabel, numpy, pandas — created by the KUL_NIS installer's env-pyfmri
section, see KUL_main_functions.sh). Exits if that env isn't found.

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
     -c:  MAX CORES the script may use (auto-scheduling). RECOMMENDED.
            Give this one number (e.g. -c 48, matching the parent pipeline's
            thread budget) and the script decides -j/-J/-n/-T per task:
              prep : prep_jobs = min(#runs, cores)
              GLM  : jobs = min(#analyses, cores); n_jobs = cores/jobs; BLAS=1
            so cores are always packed - many analyses -> more parallel jobs;
            few analyses -> each GLM uses more cores. Overrides -j/-J/-n/-T.
     -j:  (manual) number of GLM analyses in parallel (default: 1 = serial)
     -J:  (manual) number of preprocessing (mask+SUSAN) jobs in parallel
            (default: same as -j; SUSAN is light so this can be higher)
     -n:  (manual) Nilearn n_jobs = cores INSIDE a single GLM fit (default: 1)
     -T:  (manual) BLAS/OMP threads per GLM job (default: 1)
     -C:  conda env to use instead of \$KUL_PYFMRI_ENV (default 'pyfMRI'; you
            shouldn't normally need this — KUL_clinical_fmridti.sh's -y flag
            sets it for you if you do)
     -v:  verbose (0=silent, 1=normal, 2=verbose; default=1)

   Use -c for hands-off scheduling, or the -j/-J/-n/-T knobs for manual control.
   Manual CPU budget (GLM phase): roughly  -j x -n x -T  <=  physical cores.


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
verbose_level=1
smooth_fwhm=0  # 0 = adaptive (mean voxel size)
pfwe=0.01
max_cores=0       # -c : total core budget for auto-scheduling (0 = manual mode)
parallel_jobs=1   # -j : concurrent GLM analyses (1 = serial, current behaviour)
prep_jobs=0       # -J : concurrent preprocessing jobs (0 -> default to -j)
nl_njobs=1        # -n : Nilearn n_jobs inside a single GLM fit
nl_threads=1      # -T : BLAS/OMP threads per parallel job


# Set required options
p_flag=0
s_flag=0
pyfmri_env="$KUL_PYFMRI_ENV"

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:S:P:c:j:J:n:T:C:v:" OPT; do

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
        c) #max cores (auto-scheduling)
            max_cores=$OPTARG
        ;;
        j) #parallel GLM analyses
            parallel_jobs=$OPTARG
        ;;
        J) #parallel preprocessing jobs
            prep_jobs=$OPTARG
        ;;
        n) #nilearn n_jobs inside one GLM
            nl_njobs=$OPTARG
        ;;
        T) #BLAS threads per job
            nl_threads=$OPTARG
        ;;
        C) #conda env override (defaults to $KUL_PYFMRI_ENV)
            pyfmri_env=$OPTARG
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

# preprocessing concurrency defaults to the GLM concurrency unless -J is given
if [ "$prep_jobs" -le 0 ] 2>/dev/null; then
    prep_jobs=$parallel_jobs
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

function KUL_throttle {
    # Block until fewer than $1 background jobs are running.
    # Portable: uses `wait -n` when available (bash>=4.3), else polls.
    local max="$1"
    [ "$max" -lt 1 ] && max=1
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$max" ]; do
        wait -n 2>/dev/null || sleep 0.5
    done
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

function KUL_prep_run {
    # Heavy per-run preprocessing: brain-mask -> SUSAN smooth -> re-mask, plus
    # confound filtering. Fully independent per run (writes only to this run's
    # deterministic filenames), so it is safe to run in parallel. Reads the
    # GLOBAL g_* arrays by index ($1) - these span ALL runs of ALL tasks, so the
    # prep pool parallelizes across tasks, not just within one.
    local i="$1"
    local task_file="${g_bold[$i]}"
    local mask="${g_mask[$i]}"
    local sigma="${g_sigma[$i]}"
    local smooth_file="${g_smooth[$i]}"
    local boldref="${g_boldref[$i]}"
    local filter_input="${g_filterinput[$i]}"
    local filter_output="${g_confounds[$i]}"

    echo "  [prep] run $i ($(basename $task_file))"

    # ---- brain-mask the preproc BOLD BEFORE smoothing --------------------
    # fmriprep (>=23.2) leaves the MNI/func output UNMASKED. With a tight FOV
    # the sinc/Lanczos resampling extrapolates a bright halo past the superior
    # brain edge; SUSAN would smear it into cortex. Zero it out first
    # (mask -> smooth -> mask), and mask the SUSAN reference too.
    local masked_file="$fmridatadir/$(basename ${task_file} ".nii.gz")_masked.nii.gz"
    if [ ! -f "$masked_file" ]; then
        echo "  [prep] masking $task_file -> $masked_file"
        fslmaths "$task_file" -mas "$mask" "$masked_file"
    fi

    # brightness threshold = 0.6666 * in-brain median
    local p50=$(fslstats "$masked_file" -k "$mask" -p 50)
    local bt=$(echo "$p50 * 0.66666" | bc -l)
    echo "  [prep] brightness threshold: $bt"

    # mask the SUSAN reference (boldref) with the same brain mask
    local boldref_masked="$fmridatadir/$(basename ${boldref} ".nii.gz")_masked.nii.gz"
    if [ ! -f "$boldref_masked" ]; then
        fslmaths "$boldref" -mas "$mask" "$boldref_masked"
    fi

    # SUSAN smoothing (masked input + masked USAN, then clip the output)
    if [ ! -f "$smooth_file" ]; then
        echo "  [prep] SUSAN smoothing -> $smooth_file"
        susan "$masked_file" $bt $sigma 3 1 1 "$boldref_masked" $bt "$smooth_file"
        fslmaths "$smooth_file" -mas "$mask" "$smooth_file"
        gunzip -f "$smooth_file.gz"
    else
        echo "  [prep] smoothed file already exists: $smooth_file"
    fi

    # filter confounds (KUL_tsv_filter reads $filter_input / $filter_output)
    echo "  [prep] filtering confounds $filter_input -> $filter_output"
    KUL_tsv_filter
}

function KUL_compute_nilearn {
    # Nilearn first-level GLM (replaces KUL_compute_SPM_matlab + KUL_threshold_SPM_Bizzi).
    # Caller must set: nl_bolds[] nl_confounds_list[] nl_tr nl_mask nl_taskname fmrifile spm_type
    if [ $spm_type -eq 1 ]; then
        echo "   Nilearn analysis without confounds"
        fmriresults="$computedir/RESULTS/stats_$fmrifile"
        nl_hrf_deriv=0
        use_confounds=0
        global_result=${globalresultsdir_all}/afMRI_${fmrifile}.nii
        wc_suffix=""
    else
        echo "   Nilearn analysis with confounds"
        fmriresults="$computedir/RESULTS/stats_${fmrifile}_wc"
        nl_hrf_deriv=1
        use_confounds=1
        global_result=${globalresultsdir_all}/afMRI_${fmrifile}_wc.nii
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
        --n-jobs $nl_njobs \
        $conf_args"
    echo "  $cmd"
    # Cap BLAS/OMP threads so parallel jobs (-j) don't oversubscribe the CPU,
    # and give each job its own log so parallel output doesn't interleave.
    local thr="${nl_threads:-1}"
    local env_pfx="OMP_NUM_THREADS=$thr OPENBLAS_NUM_THREADS=$thr MKL_NUM_THREADS=$thr NUMEXPR_NUM_THREADS=$thr"
    if [ $verbose_level -lt 2 ]; then
        local job_log="KUL_LOG/$script/sub-${participant}_${fmrifile}${wc_suffix}.log"
        mkdir -p "$(dirname "$job_log")"
        eval "$env_pfx $cmd" >> "$job_log" 2>&1
    else
        eval "$env_pfx $cmd"
    fi

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
            # Hardwired downstream selection: only wc + p001unc_k50 lands in
            # globalresultsdir (RESULTS/.../SPM); everything else goes to
            # globalresultsdir_all (RESULTS/.../SPM_all).
            if [ "$wc_suffix" == "_wc" ] && [ "$thresh_tag" == "p001unc_k50" ]; then
                output="${globalresultsdir}/afMRI_${fmrifile}${wc_suffix}_${thresh_tag}.nii"
            else
                output="${globalresultsdir_all}/afMRI_${fmrifile}${wc_suffix}_${thresh_tag}.nii"
            fi
            transform="$mni_to_t1w"
            KUL_antsApply_Transform
        fi
    done
}

# MAIN --------------------------------------------------------------
# conda env existence check (mirrors KUL_run_rsfMRI_networks.sh)
if ! conda env list | awk '{print $1}' | grep -qx "$pyfmri_env"; then
    echo "ERROR: conda env '$pyfmri_env' was not found (checked 'conda env list')." >&2
    echo "  Run the KUL_NIS installer's env-pyfmri section, or pass -C <env> to use a different one." >&2
    exit 1
fi
KUL_activate_conda_env "$pyfmri_env"

python_exe=$(which python3)
# the Nilearn GLM worker lives next to this script (share/nilearn/)
nilearn_glm_script="$kul_main_dir/share/nilearn/KUL_nilearn_glm.py"

if [ $KUL_DEBUG -gt 0 ]; then
    echo "python3 lives at $python_exe"
    echo "nilearn glm script: $nilearn_glm_script"
fi

if [ ! -f "$nilearn_glm_script" ]; then
    echo "Nilearn GLM worker not found at $nilearn_glm_script. Exitting"
    exit 1
fi

if ! "$python_exe" -c "import nilearn, nibabel, numpy, pandas" 2>/dev/null ; then
    echo "ERROR: conda env '$pyfmri_env' is missing required packages (nilearn, nibabel, numpy, pandas)." >&2
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
# globalresultsdir holds ONLY the hardwired downstream-selected outcome (wc,
# p001unc_k50 thresholded) per task; every other combination (nc, FWE, raw
# unthresholded) goes to globalresultsdir_all instead.
globalresultsdir=$cwd/RESULTS/sub-$participant/SPM
globalresultsdir_all=$cwd/RESULTS/sub-$participant/SPM_all


if [ $KUL_DEBUG -gt 0 ]; then
    echo "kulderivativesdir: $kulderivativesdir"
    echo "fmridatadir: $fmridatadir"
    echo "fmriprepdir: $fmriprepdir"
    echo "globalresultsdir: $globalresultsdir"
    echo "globalresultsdir_all: $globalresultsdir_all"
fi

mkdir -p $fmridatadir
mkdir -p $scriptsdir
mkdir -p $confoundsdir
mkdir -p $computedir/RESULTS
mkdir -p $globalresultsdir
mkdir -p $globalresultsdir_all

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

    # =====================================================================
    # GLOBAL PHASE 1 (serial, cheap): resolve deterministic paths/values for
    # EVERY non-rest run of EVERY task, filling the global g_* arrays. g_task
    # records each run's (coarse) task so the GLM phase can regroup them.
    # =====================================================================
    g_bold=(); g_task=(); g_taskname=(); g_TR=()
    g_sigma=(); g_smooth=(); g_confounds=(); g_mask=(); g_boldref=(); g_filterinput=()
    for match in "${fmriprep_match[@]}"; do
        _coarse_task=$(basename "$match" | sed -E 's/.*_task-([A-Za-z0-9]+).*_desc-preproc_bold\.nii\.gz/\1/')
        [[ "$_coarse_task" == *"rest"* ]] && continue
        task_file="$match"

        taskname=$(basename "$task_file" | sed -E 's/.*_task-([A-Za-z0-9]+(_run-[0-9]+)?).*_desc-preproc_bold\.nii\.gz/\1/')
        TR=($(mrinfo $task_file -spacing | awk '{print $(NF)}'))
        spacing=($(mrinfo $task_file -spacing))
        mean=$(echo "(${spacing[0]} + ${spacing[1]} + ${spacing[2]}) / 3" | bc -l)
        if (( $(echo "$smooth_fwhm > 0" | bc -l) )); then
            sigma=$(echo "$smooth_fwhm / 2.3548" | bc -l)
        else
            sigma=$mean
        fi
        mask=$(dirname ${task_file})/$(basename ${task_file} "-preproc_bold.nii.gz")-brain_mask.nii.gz
        smooth_file="$fmridatadir/$(basename ${task_file} "-preproc_bold.nii.gz")_smooth.nii"
        boldref=$(dirname ${task_file})/$(basename ${task_file} "_desc-preproc_bold.nii.gz")_boldref.nii.gz
        _task_run_key=$(basename "$task_file" | grep -oE 'task-[A-Za-z0-9]+(_run-[0-9]+)?')
        filter_input=$(find $(dirname ${task_file}) -name "*${_task_run_key}_desc-confounds_timeseries.tsv" | head -1)
        confounds_file="$confoundsdir/${_task_run_key}_confounds.txt"

        g_bold+=("$task_file");    g_task+=("$_coarse_task"); g_taskname+=("$taskname"); g_TR+=("$TR")
        g_sigma+=("$sigma");       g_smooth+=("$smooth_file"); g_confounds+=("$confounds_file")
        g_mask+=("$mask");         g_boldref+=("$boldref");    g_filterinput+=("$filter_input")
        echo "  [fill] $taskname  TR=$TR  sigma=$sigma"
    done

    n_runs_total=${#g_bold[@]}
    if [ "$n_runs_total" -eq 0 ]; then
        echo "No non-rest task runs found. Nothing to do."
    else

    # =====================================================================
    # GLOBAL PHASE 2 (parallel, heavy): mask + SUSAN + confound filtering for
    # ALL runs of ALL tasks in ONE pool (so a 1-run task preprocesses next to
    # another task's runs, not after it). One barrier before any GLM.
    # =====================================================================
    if [ "$max_cores" -gt 0 ]; then
        prep_jobs=$(( n_runs_total < max_cores ? n_runs_total : max_cores ))
    fi
    echo "Preprocessing (mask+SUSAN+confounds) for ALL $n_runs_total runs in parallel (prep_jobs=$prep_jobs)"
    for i in "${!g_bold[@]}"; do
        KUL_throttle "$prep_jobs"
        if [ $verbose_level -lt 2 ]; then
            prep_log="KUL_LOG/$script/sub-${participant}_prep_${g_taskname[$i]}.log"
            mkdir -p "$(dirname "$prep_log")"
            KUL_prep_run "$i" >> "$prep_log" 2>&1 &
        else
            KUL_prep_run "$i" &
        fi
    done
    wait   # barrier: ALL preprocessing done before any GLM starts

    # =====================================================================
    # PHASE 3 (parallel, heavy): GLM analyses for ALL tasks in ONE pool, so
    # Lip and Taal analyses run together within the core budget (not task by
    # task). Total analyses A = 2 per run (nc/wc) over every run, + 2 per
    # multi-run task (aggregate nc/wc). Auto-schedule from that total, dispatch
    # everything, then a single barrier.
    # =====================================================================
    # count total analyses across all tasks
    A_total=$(( 2 * n_runs_total ))
    for task in "${unique_tasks[@]}"; do
        [[ "$task" == *"rest"* ]] && continue
        _c=0; for k in "${!g_task[@]}"; do [ "${g_task[$k]}" == "$task" ] && _c=$((_c+1)); done
        [ "$_c" -gt 1 ] && A_total=$(( A_total + 2 ))
    done

    if [ "$max_cores" -gt 0 ]; then
        parallel_jobs=$(( A_total < max_cores ? A_total : max_cores )); [ "$parallel_jobs" -lt 1 ] && parallel_jobs=1
        nl_njobs=$(( max_cores / parallel_jobs )); [ "$nl_njobs" -lt 1 ] && nl_njobs=1
        nl_threads=1
        echo "[auto] GLM: cores=$max_cores total_analyses=$A_total -> -j$parallel_jobs -n$nl_njobs -T1"
    fi
    echo "Dispatching ALL GLM analyses across tasks (total=$A_total; parallel_jobs=$parallel_jobs, n_jobs=$nl_njobs)"

    # --- per-run analyses for every run of every task ---
    for k in "${!g_bold[@]}"; do
        fmrifile="${g_taskname[$k]}"
        nl_bolds=("${g_smooth[$k]}")
        nl_confounds_list=("${g_confounds[$k]}")
        nl_tr="${g_TR[$k]}"
        nl_mask="${g_mask[$k]}"
        nl_taskname="${g_taskname[$k]}"
        for spm_type in 1 2; do
            KUL_throttle "$parallel_jobs"
            KUL_compute_nilearn &
        done
    done

    # --- aggregate (fixed-effects) analyses for each multi-run task ---
    for task in "${unique_tasks[@]}"; do
        [[ "$task" == *"rest"* ]] && continue
        idxs=(); for k in "${!g_task[@]}"; do [ "${g_task[$k]}" == "$task" ] && idxs+=("$k"); done
        [ ${#idxs[@]} -gt 1 ] || continue

        # TR sanity check (aggregate assumes a common TR across runs)
        _trs=(); for k in "${idxs[@]}"; do _trs+=("${g_TR[$k]}"); done
        [ $(printf "%s\n" "${_trs[@]}" | sort -u | wc -l) -gt 1 ] && echo "  Warning: Multiple TRs for task $task."

        fmrifile="$task"
        nl_bolds=(); nl_confounds_list=()
        for k in "${idxs[@]}"; do nl_bolds+=("${g_smooth[$k]}"); nl_confounds_list+=("${g_confounds[$k]}"); done
        nl_tr="${g_TR[${idxs[0]}]}"
        nl_mask="${g_mask[${idxs[0]}]}"
        nl_taskname="$task"
        for spm_type in 1 2; do
            KUL_throttle "$parallel_jobs"
            KUL_compute_nilearn &
        done
    done

    wait   # barrier: all GLM analyses done

    fi   # n_runs_total > 0

    # cleanup
    #rm -rf $fmridatadir

    touch KUL_LOG/sub-${participant}_SPM.done
    echo "Done computing Nilearn GLM"
else
    echo "Nilearn analysis already done"
fi
