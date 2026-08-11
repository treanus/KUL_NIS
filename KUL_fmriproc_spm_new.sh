#!/bin/bash
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@uzleuven.be
#
# v0.1 - dd 19/01/2019 - jurassic version
version="v2.0 - dd 03/07/2027"

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

`basename $0` performs an automated task based fMRI spm12 GLM analysis assuming:
    - a 30 seconds REST followed by 30 seconds TASK epochs
    - having run fmriprep with aroma

Note: requires matlab and spm12 installed

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
            One number (e.g. -c 48, matching the parent pipeline budget) and
            the script sets -j/-J/-T per task:
              prep : prep_jobs = min(#runs, cores)
              GLM  : jobs = min(#analyses, cores); MATLAB threads = cores/jobs
            NOTE each parallel GLM is a full MATLAB instance (RAM + possibly a
            license seat each), so keep -c within your license/RAM limits.
     -j:  (manual) number of SPM/MATLAB analyses in parallel (default: 1)
     -J:  (manual) number of preprocessing (mask+SUSAN) jobs in parallel
            (default: same as -j; SUSAN is light so this can be higher)
     -T:  (manual) MATLAB computational threads per instance
            (default: 0 = unrestricted, i.e. one MATLAB uses all cores)
     -v:  verbose (0=silent, 1=normal, 2=verbose; default=1)

   Use -c for hands-off scheduling, or the -j/-J/-T knobs for manual control.
   SPM has no within-GLM n_jobs like Nilearn, so leftover cores are given to
   MATLAB's internal multithreading (maxNumCompThreads) rather than a voxel-loop
   split - a weaker speed-up per analysis than the Nilearn port.


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
parallel_jobs=1   # -j : concurrent SPM/MATLAB analyses (1 = serial)
prep_jobs=0       # -J : concurrent preprocessing jobs (0 -> default to -j)
matlab_threads=0  # -T : MATLAB comp threads per instance (0 = unrestricted)


# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:S:P:c:j:J:T:v:" OPT; do

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
        j) #parallel SPM analyses
            parallel_jobs=$OPTARG
        ;;
        J) #parallel preprocessing jobs
            prep_jobs=$OPTARG
        ;;
        T) #MATLAB comp threads per instance
            matlab_threads=$OPTARG
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

function KUL_threshold_SPM_Bizzi {
    if [ ! -f "$fmriresults/SPM.mat" ]; then
        echo " Bizzi thresholding skipped: SPM.mat not found in $fmriresults"
        return
    fi
    # Detect _wc suffix from fmriresults path so output names don't collide
    local wc_suffix=""
    [[ "$fmriresults" == *_wc ]] && wc_suffix="_wc"
    spm_bizzi_script="${scriptsdir}/thresh_Bizzi_${fmrifile}${wc_suffix}.m"
    spm_bizzi_script=${spm_bizzi_script/run-/run}
    cp "$kul_main_dir/share/spm12/spm12_threshold_Bizzi.m" "$spm_bizzi_script"
    sed -i.bck "s|###FMRIRESULTS###|$fmriresults|" "$spm_bizzi_script"
    sed -i.bck "s|###PFWE###|$pfwe|g" "$spm_bizzi_script"
    rm -f "${spm_bizzi_script}.bck"
    local mct=""
    [ "${matlab_threads:-0}" -gt 0 ] && mct="maxNumCompThreads(${matlab_threads});"
    cmd="$matlab_exe -nodisplay -nosplash -nodesktop -r \"${mct}run('$spm_bizzi_script');exit;\""
    if [ $verbose_level -lt 2 ]; then
        local job_log="KUL_LOG/$script/sub-${participant}_${fmrifile}${wc_suffix}.log"
        mkdir -p "$(dirname "$job_log")"
        eval "$cmd" >> "$job_log" 2>&1
    else
        eval "$cmd"
    fi
    pfwe_tag=$(echo "$pfwe" | sed 's/0\.//' | sed 's/0*$//')
    fwe_tag="FWE${pfwe_tag}_k50"
    mni_to_t1w=$(find_first_match "${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_*from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5" "MNI-to-T1w transform")
    find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
    reference=${find_T1w[0]}
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
    # confound filtering. Independent per run (writes only this run's files), so
    # safe to run in parallel. Reads the GLOBAL g_* arrays by index ($1) - these
    # span ALL runs of ALL tasks, so the prep pool parallelizes across tasks.
    # IMPORTANT: the _masked intermediates go to $prep_tmpdir (a SUBDIR of
    # fmridatadir). SPM selects its input by a non-recursive FPList filter on the
    # task name inside fmridatadir, so only the final _smooth.nii must live there
    # - otherwise the masked/boldref files would be picked up as extra scans.
    local i="$1"
    local task_file="${g_bold[$i]}"
    local mask="${g_mask[$i]}"
    local sigma="${g_sigma[$i]}"
    local smooth_file="${g_smooth[$i]}"
    local boldref="${g_boldref[$i]}"
    local filter_input="${g_filterinput[$i]}"
    local filter_output="${g_confounds[$i]}"

    echo "  [prep] run $i ($(basename $task_file))"

    # brain-mask the preproc BOLD BEFORE smoothing (halo removal; see header)
    local masked_file="$prep_tmpdir/$(basename ${task_file} ".nii.gz")_masked.nii.gz"
    if [ ! -f "$masked_file" ]; then
        fslmaths "$task_file" -mas "$mask" "$masked_file"
    fi

    # brightness threshold = 0.6666 * in-brain median
    local p50=$(fslstats "$masked_file" -k "$mask" -p 50)
    local bt=$(echo "$p50 * 0.66666" | bc -l)

    # mask the SUSAN reference (boldref) with the same brain mask
    local boldref_masked="$prep_tmpdir/$(basename ${boldref} ".nii.gz")_masked.nii.gz"
    if [ ! -f "$boldref_masked" ]; then
        fslmaths "$boldref" -mas "$mask" "$boldref_masked"
    fi

    # SUSAN smoothing (masked input + masked USAN, then clip the output)
    if [ ! -f "$smooth_file" ]; then
        susan "$masked_file" $bt $sigma 3 1 1 "$boldref_masked" $bt "$smooth_file"
        fslmaths "$smooth_file" -mas "$mask" "$smooth_file"
        gunzip -f "$smooth_file.gz"
    fi

    # filter confounds (KUL_tsv_filter reads $filter_input / $filter_output)
    KUL_tsv_filter
}

function KUL_compute_SPM_matlab {

    # which type of SPM analysis
    if [ $spm_type -eq 1 ]; then
        echo "   SPM analysis without confounds"
        fmriresults="$computedir/RESULTS/stats_$fmrifile"
        # prepare the job and config files
        spm_participant_config_file="${scriptsdir}/stats_${fmrifile}.m" #participant config file
        spm_participant_job_file="${scriptsdir}/stats_${fmrifile}_job.m" #participant job file
        spm_hrfderivs="[0 0]"
        global_result=${globalresultsdir_all}/afMRI_${fmrifile}.nii
    else
        echo "   SPM analysis with confounds"
        fmriresults="$computedir/RESULTS/stats_${fmrifile}_wc"
        # prepare the job and config files
        spm_participant_config_file="${scriptsdir}/stats_${fmrifile}_wc.m" #participant config file
        spm_participant_job_file="${scriptsdir}/stats_${fmrifile}_job_wc.m" #participant job file
        spm_hrfderivs="[1 0]"
        global_result=${globalresultsdir_all}/afMRI_${fmrifile}_wc.nii
    fi

    # clean a possible old result
    rm -rf $fmriresults
    mkdir -p $fmriresults

    # get rid of - in filename, since this breaks -r in matlab
    spm_participant_config_file=${spm_participant_config_file/run-/run}
    spm_participant_job_file=${spm_participant_job_file/run-/run}
    cp $spm_template_config_file $spm_participant_config_file
    cp $spm_template_job_file $spm_participant_job_file
    sed -i.bck "s|###JOBFILE###|$spm_participant_job_file|" $spm_participant_config_file
    sed -i.bck "s|###FMRIDIR###|$fmridatadir|" $spm_participant_job_file
    #sed -i.bck "s|###FMRIFILE###|$fmrifile|" $spm_participant_job_file
    sed -i.bck "s|###FMRIRESULTS###|$fmriresults|" $spm_participant_job_file
    sed -i.bck "s|###HRFDERIVS###|$spm_hrfderivs|" $spm_participant_job_file
    #sed -i.bck "s|###TR###|$TR|" $spm_participant_job_file
    #sed -i.bck "s|###CONFOUNDSFILE###|$run_confounds|" $spm_participant_job_file

    
    if [ $spm_nruns -eq 1 ]; then
        # Here $i is the index of the current seperate run analyis
        spm_taskname=${tf_taskname[$i]}
        spm_filter=${tf_filter[$i]}
        spm_TR=${tf_TR[$i]}
        if [ $spm_type -eq 1 ]; then
            spm_confounds_file=""
        else
            spm_confounds_file="${tf_confounds[$i]}"
        fi
        sed -i.bck "s|###FMRIFILE###|$spm_filter|" $spm_participant_job_file
        sed -i.bck "s|###CONFOUNDSFILE###|$spm_confounds_file|" $spm_participant_job_file
        sed -i.bck "s|###TR###|$spm_TR|" $spm_participant_job_file
    
    else
        # Here multiple runs together - we assume all runs have the same TR
        spm_TR=${tf_TR[0]}
        sed -i.bck "s|###TR###|$spm_TR|" $spm_participant_job_file

        for ((j=0; j<spm_nruns; j++)); do
            spm_taskname="${tf_taskname[$j]}"
            spm_filter="${tf_filter[$j]}"
            if [ $spm_type -eq 1 ]; then
                spm_confounds_file=""
            else
                spm_confounds_file="${tf_confounds[$j]}"
            fi
            sed -i.bck "s|###FMRIFILE$j###|$spm_filter|" $spm_participant_job_file
            cmd="sed -i.bck \"s|###CONFOUNDSFILE$j###|$spm_confounds_file|\" $spm_participant_job_file"
            eval $cmd
            
        done
    fi

    rm -f "${spm_participant_config_file}.bck"
    rm -f "${spm_participant_job_file}.bck"

    # call matlab and execute
    # cap MATLAB computational threads when parallelizing (0 = unrestricted);
    # write each analysis to its own log so parallel output does not interleave.
    local mct=""
    [ "${matlab_threads:-0}" -gt 0 ] && mct="maxNumCompThreads(${matlab_threads});"
    local job_tag="$fmrifile"; [ $spm_type -eq 2 ] && job_tag="${fmrifile}_wc"
    cmd="$matlab_exe -nodisplay -nosplash -nodesktop -r \"${mct}run('$spm_participant_config_file');exit;\""
    echo "$cmd"
    if [ $verbose_level -lt 2 ]; then
        local job_log="KUL_LOG/$script/sub-${participant}_${job_tag}.log"
        mkdir -p "$(dirname "$job_log")"
        eval "$cmd" >> "$job_log" 2>&1
    else
        eval "$cmd"
    fi

    # SPM output is in MNI space; warp back to T1w space for display
    input=$fmriresults/spmT_0001.nii
    output=$global_result
    transform=$(find_first_match "${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_*from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5" "MNI-to-T1w transform")
    find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
    reference=${find_T1w[0]}
    KUL_antsApply_Transform

} 

# MAIN --------------------------------------------------------------
matlab_exe=$(which matlab)

if [ $KUL_DEBUG -gt 0 ]; then 
    echo "matlab lives at $matlab_exe"
fi

if [[ -z "$matlab_exe" ]]; then
    echo "Matlab is required but not found on path. Exitting"
    exit 1
fi


KUL_check_participant


#  setup variables
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
computedir="$kulderivativesdir/sub-$participant/SPM"
fmridatadir="$computedir/fmridata"
prep_tmpdir="$fmridatadir/intermediate"   # masked intermediates (kept out of FPList)
scriptsdir="$computedir/scripts"
confoundsdir="$computedir/confounds"
fmriprepdir="${cwd}/fmriprep/sub-$participant"
# globalresultsdir holds ONLY the hardwired downstream-selected outcome (wc,
# p001unc_k50 Bizzi-thresholded) per task; every other combination (nc, FWE,
# raw unthresholded) goes to globalresultsdir_all instead.
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
mkdir -p $prep_tmpdir
mkdir -p $scriptsdir
mkdir -p $confoundsdir
mkdir -p $computedir/RESULTS
mkdir -p $globalresultsdir
mkdir -p $globalresultsdir_all

# Use MNI-space fmriprep output so SPM stats are in MNI and the MNI→T1w
# warp-back at the end produces correctly aligned native-space results.
fmriprep_output_type="_space-MNI152NLin2009cAsym"


if [ $verbose_level -lt 2 ] ; then
    str_silent_SPM=" >> KUL_LOG/$script/sub-${participant}_spm12.log"
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
    g_bold=(); g_task=(); g_taskname=(); g_TR=(); g_filter=()
    g_sigma=(); g_smooth=(); g_confounds=(); g_mask=(); g_boldref=(); g_filterinput=()
    for match in "${fmriprep_match[@]}"; do
        _coarse_task=$(basename "$match" | sed -E 's/.*_task-([A-Za-z0-9]+).*_desc-preproc_bold\.nii\.gz/\1/')
        # Exclude resting-state runs from this task/activation GLM pipeline.
        # "*rest*" alone misses task labels like "rsfMRI" (no literal "rest"
        # substring), which then got processed as an activation task and
        # ended up named afMRI_rsfMRI_... — case-insensitive match on both
        # "rest" and "rsfmri" so any of task-rest/task-RestingState/task-rsfMRI
        # naming is excluded.
        _coarse_task_lc=$(echo "$_coarse_task" | tr '[:upper:]' '[:lower:]')
        [[ "$_coarse_task_lc" == *"rest"* || "$_coarse_task_lc" == *"rsfmri"* ]] && continue
        task_file="$match"

        # Match task- and run- as SEPARATE entities -- see the identical fix and
        # rationale in KUL_fmriproc_nilearn_new.sh. taskname becomes $fmrifile,
        # i.e. stats_<name>/ and afMRI_<name>.nii, so a bare 'TAAL' for every run
        # made all per-run GLMs and the aggregate write the same two output dirs
        # concurrently.
        #
        # In THIS engine it also corrupted the input selection, not just the
        # names: taskname is substituted into the SPM job as ###FMRIFILE###,
        # which is cfg_basicio's file_fplist.filter -- a REGEXP matched against
        # everything in $fmridatadir. A filter of 'TAAL' matches the smoothed
        # files of run-01 AND run-02, so each supposedly single-run GLM was fed
        # both runs, and each session of the multi-run job got the same doubled
        # set. g_filter below is therefore built to match exactly one run.
        _task_lbl=$(basename "$task_file" | grep -oE 'task-[A-Za-z0-9]+' | head -1)
        _run_lbl=$(basename "$task_file" | grep -oE '_run-[0-9]+' | head -1)
        _task_run_key="${_task_lbl}${_run_lbl}"
        taskname="${_task_lbl#task-}${_run_lbl}"
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
        # (_task_lbl/_run_lbl/_task_run_key are derived above, alongside taskname,
        # so the confounds key and the output name are built from the same two
        # entities -- both were previously wrong in the same way.)
        #
        # SPM file-selector filter for exactly this run. Anchored on the smoothed
        # file's own basename, so it cannot match a sibling run no matter which
        # entities sit between task- and run-. '.' is left unescaped: as a regexp
        # it matches the literal dot too, and escaping it here would have to
        # survive the sed substitution below.
        spm_filter="^$(basename "$smooth_file")$"
        # Derive the confounds sidecar from the BOLD name rather than globbing:
        # everything from _space- onward is fMRIPrep's output-space decoration
        # (including any _res-<N> token), and what remains is the confounds stem.
        _bold_stem=$(basename "$task_file"); _bold_stem="${_bold_stem%%_space-*}"
        filter_input="$(dirname ${task_file})/${_bold_stem}_desc-confounds_timeseries.tsv"
        if [ ! -s "$filter_input" ]; then
            echo "  ERROR: confounds TSV missing or empty for ${_task_run_key}: $filter_input" >&2
        fi
        confounds_file="$confoundsdir/${_task_run_key}_confounds.txt"

        g_bold+=("$task_file");    g_task+=("$_coarse_task"); g_taskname+=("$taskname"); g_TR+=("$TR")
        g_filter+=("$spm_filter")
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
    # PHASE 3 (parallel, heavy): SPM/MATLAB GLM analyses for ALL tasks in ONE
    # pool. Total analyses A = 2 per run (nc/wc) + 2 per multi-run task
    # (aggregate nc/wc). Auto-schedule MATLAB threads from that total. Each
    # analysis is a separate MATLAB instance writing to its own stats_*/afMRI_*
    # outputs; backgrounding snapshots the per-analysis session arrays
    # (tf_taskname/tf_confounds/tf_TR), spm_nruns, templates, i, fmrifile,
    # spm_type. compute+threshold are grouped (threshold needs the SPM.mat).
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
        matlab_threads=$(( max_cores / parallel_jobs )); [ "$matlab_threads" -lt 1 ] && matlab_threads=1
        echo "[auto] GLM: cores=$max_cores total_analyses=$A_total -> -j$parallel_jobs MATLAB threads/job=$matlab_threads"
    fi
    echo "Dispatching ALL SPM GLM analyses across tasks (total=$A_total; parallel_jobs=$parallel_jobs, matlab_threads=$matlab_threads)"

    # --- per-run analyses for every run of every task ---
    spm_nruns=1
    for k in "${!g_bold[@]}"; do
        spm_template_config_file="$kul_main_dir/share/spm12/spm12_new_fmri_stats_1run.m"
        spm_template_job_file="$kul_main_dir/share/spm12/spm12_new_fmri_stats_1run_job.m"
        # single-element session arrays consumed by KUL_compute_SPM_matlab ($i=0)
        tf_taskname=("${g_taskname[$k]}"); tf_confounds=("${g_confounds[$k]}"); tf_TR=("${g_TR[$k]}")
        tf_filter=("${g_filter[$k]}")
        i=0
        fmrifile="${g_taskname[$k]}"
        for spm_type in 1 2; do
            KUL_throttle "$parallel_jobs"
            { KUL_compute_SPM_matlab; KUL_threshold_SPM_Bizzi; } &
        done
    done

    # --- aggregate (multi-run) analyses per task ---
    for task in "${unique_tasks[@]}"; do
        [[ "$task" == *"rest"* ]] && continue
        idxs=(); for k in "${!g_task[@]}"; do [ "${g_task[$k]}" == "$task" ] && idxs+=("$k"); done
        spm_nruns=${#idxs[@]}
        [ "$spm_nruns" -gt 1 ] || continue

        if [ "$spm_nruns" -eq 2 ]; then
            spm_template_config_file="$kul_main_dir/share/spm12/spm12_new_fmri_stats_2runs.m"
            spm_template_job_file="$kul_main_dir/share/spm12/spm12_new_fmri_stats_2runs_job.m"
        elif [ "$spm_nruns" -eq 3 ]; then
            spm_template_config_file="$kul_main_dir/share/spm12/spm12_new_fmri_stats_3runs.m"
            spm_template_job_file="$kul_main_dir/share/spm12/spm12_new_fmri_stats_3runs_job.m"
        else
            echo "  Warning: aggregate not defined for >3 runs; skipping task $task aggregate."
            continue
        fi

        # session arrays in run order for this task (consumed via j-loop + tf_TR[0])
        tf_taskname=(); tf_confounds=(); tf_TR=(); tf_filter=()
        for k in "${idxs[@]}"; do
            tf_taskname+=("${g_taskname[$k]}"); tf_confounds+=("${g_confounds[$k]}"); tf_TR+=("${g_TR[$k]}")
            tf_filter+=("${g_filter[$k]}")
        done
        fmrifile="$task"
        for spm_type in 1 2; do
            KUL_throttle "$parallel_jobs"
            { KUL_compute_SPM_matlab; KUL_threshold_SPM_Bizzi; } &
        done
    done

    wait   # barrier: all GLM analyses done

    fi   # n_runs_total > 0

    # cleanup
    #rm -rf $fmridatadir

    touch KUL_LOG/sub-${participant}_SPM.done
    echo "Done computing SPM"
else
    echo "SPM analysis already done"
fi
