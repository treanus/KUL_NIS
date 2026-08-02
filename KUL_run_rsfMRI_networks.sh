#!/bin/bash
# Bash shell script for presurgical / eloquent-cortex rsfMRI network mapping
#
# Denoises EVERY functional BOLD run for a subject (resting-state AND any
# task-fMRI runs, including runs where the task paradigm could not be
# reliably performed), then runs seed-based + RSN-level functional
# connectivity (share/rsfmri_pipeline/) and generates a patient-own-data
# PDF report. Opt-in, standalone step — does not replace KUL_fmriproc_conn.sh.
#
# @ Ahmed Radwan - KU Leuven, Translational MRI - ahmed.radwan@kuleuven.be
version="v1.0"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs presurgical/eloquent-cortex rsfMRI network mapping:
    denoise (all BOLD runs, rest + task) -> SynthSeg -> seed/RSN FC -> PDF report

Usage:

  `basename $0` -p subject -c conda_env <OPT_ARGS>

Example:

  `basename $0` -p pat001 -c rsfmri_env

Required arguments:

     -p:  participant
     -c:  name of the conda env with the rsfmri_pipeline dependencies installed
          (numpy, nibabel, nilearn, pandas, scipy, matplotlib, pyyaml)

Optional arguments:

     -P:  condition profile from share/rsfmri_pipeline/config/profiles.yaml (default: Presurgical)
     -I:  also generate masked-ICA thumbnail pages (opt-in; adds runtime)
     -v:  verbose (0=silent, 1=normal, 2=verbose; default=1)

USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
verbose_level=1
rsfmri_profile="Presurgical"
include_ica=0

# Set required options
p_flag=0
c_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:c:P:Iv:" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            participant=$OPTARG
        ;;
        c) #conda env
            c_flag=1
            rsfmri_env=$OPTARG
        ;;
        P) #condition profile
            rsfmri_profile=$OPTARG
        ;;
        I) #opt-in masked-ICA thumbnail pages
            include_ica=1
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
if [ $c_flag -eq 0 ] ; then
    echo
    echo "Option -c is required: give the name of the conda env with rsfmri_pipeline dependencies." >&2
    echo
    exit 2
fi

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
    export MRTRIX_QUIET=1
fi

KUL_check_participant

# conda env existence check (mirrors KUL_dwiprep.sh's lore_sd pattern)
if ! conda env list | awk '{print $1}' | grep -qx "$rsfmri_env"; then
    echo "ERROR: conda env '$rsfmri_env' was not found (checked 'conda env list')." >&2
    exit 1
fi
KUL_activate_conda_env "$rsfmri_env"

# MAIN --------------------------------------------------------------
pipeline_dir="$kul_main_dir/share/rsfmri_pipeline"

#  setup variables — single shared (not per-participant-nested) base dirs
#  under $kulderivativesdir: the pipeline's own step code already nests by
#  sub-<ID> internally, and this script only ever handles one participant.
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
globalresultsdir="$cwd/RESULTS/sub-$participant/rsfMRI_Networks"

export RSFMRI_FMRIPREP_DIR="$cwd/fmriprep"
export RSFMRI_DENOISED_DIR="$kulderivativesdir/rsfMRI_networks/denoised"
export RSFMRI_ANALYSIS_DIR="$kulderivativesdir/rsfMRI_networks/analysis"
export RSFMRI_SYNTHSEG_DIR="$RSFMRI_ANALYSIS_DIR/synthseg"

if [ $KUL_DEBUG -gt 0 ]; then
    echo "pipeline_dir: $pipeline_dir"
    echo "RSFMRI_FMRIPREP_DIR: $RSFMRI_FMRIPREP_DIR"
    echo "RSFMRI_DENOISED_DIR: $RSFMRI_DENOISED_DIR"
    echo "RSFMRI_ANALYSIS_DIR: $RSFMRI_ANALYSIS_DIR"
    echo "RSFMRI_SYNTHSEG_DIR: $RSFMRI_SYNTHSEG_DIR"
    echo "globalresultsdir: $globalresultsdir"
fi

mkdir -p "$RSFMRI_DENOISED_DIR" "$RSFMRI_ANALYSIS_DIR" "$globalresultsdir"

done_flag="${cwd}/KUL_LOG/sub-${participant}_rsfMRI_networks.done"

if [ -f "$done_flag" ]; then
    echo "rsfMRI networks already done"
    exit 0
fi

# 1) Denoise EVERY functional run for this subject — deliberately no --task
#    filter, so resting-state AND all task-fMRI runs (HAND/FOOT/LIP/TAAL/etc,
#    regardless of whether the task GLM found activation) get pooled into the
#    same connectivity estimate. get_runs()/bold_path() in utils.py are
#    already task-label-agnostic, so no pipeline code changes are needed for
#    this beyond omitting --task here.
task_in="$kul_main_dir/KUL_fmri_denoise.sh --fmriprep $RSFMRI_FMRIPREP_DIR --out $RSFMRI_DENOISED_DIR --method nilearn --space MNI152NLin2009cAsym --sub sub-${participant}"
KUL_task_exec $verbose_level "rsfMRI networks: denoise all BOLD runs" "rsfmri_1_denoise" || { echo "ERROR: rsfMRI denoising failed for sub-${participant} — aborting (rsfMRI_networks.done will not be created)" >&2; exit 1; }

# 2) Subject-specific SynthSeg segmentation + subject-specific atlases
task_in="bash $pipeline_dir/src/step0_synthseg.sh --subjects $participant"
KUL_task_exec $verbose_level "rsfMRI networks: SynthSeg" "rsfmri_2_synthseg" || { echo "ERROR: SynthSeg failed for sub-${participant} — aborting (rsfMRI_networks.done will not be created)" >&2; exit 1; }

# 3) Seed-based FC (step1) + RSN-level FC/coupling (step2) + masked ICA (step3)
#    --profile is mandatory in run_pipeline.py's own argparse, hence the
#    Presurgical default rather than leaving it optional.
task_in="python3 $pipeline_dir/src/run_pipeline.py --profile $rsfmri_profile --subjects $participant --steps 1 2 3"
KUL_task_exec $verbose_level "rsfMRI networks: seed/RSN pipeline" "rsfmri_3_pipeline" || { echo "ERROR: rsfmri_pipeline run failed for sub-${participant} — aborting (rsfMRI_networks.done will not be created)" >&2; exit 1; }

# 4) Patient-own-data PDF report (no normative/HV comparison — step4 is out
#    of scope for this integration). A report failure is non-fatal: the
#    seed/RSN/ICA analysis outputs from step 3 already exist and are the
#    valuable artifact regardless of whether the PDF renders.
_ica_opt=""
[ $include_ica -eq 1 ] && _ica_opt="--include-ica"
task_in="python3 $pipeline_dir/src/step5_lite_report.py --profile $rsfmri_profile --rsn --patients $participant $_ica_opt"
KUL_task_exec $verbose_level "rsfMRI networks: PDF report" "rsfmri_4_report" || echo "WARNING: rsfMRI networks report generation failed for sub-${participant} — analysis outputs exist regardless (check rsfmri_4_report.error.log)"

cp -f "$RSFMRI_ANALYSIS_DIR/reports/sub-${participant}_rsfmri_networks_report.pdf" "$globalresultsdir/" 2>/dev/null

touch "$done_flag"
echo "Done: rsfMRI networks for sub-${participant}"
