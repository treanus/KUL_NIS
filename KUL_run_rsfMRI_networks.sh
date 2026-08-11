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

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001

Requires the conda env named by \$KUL_PYFMRI_ENV (default 'pyfMRI'; numpy,
nibabel, nilearn, pandas, scipy, matplotlib, pyyaml — created by the KUL_NIS
installer's env-pyfmri section, see KUL_main_functions.sh). Exits if that env
isn't found.

Required arguments:

     -p:  participant

Optional arguments:

     -c:  conda env to use instead of \$KUL_PYFMRI_ENV (default 'pyfMRI'; you
          shouldn't normally need this — KUL_clinical_fmridti.sh's -y flag
          sets it for you if you do)
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
pyfmri_env="$KUL_PYFMRI_ENV"

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
        c) #conda env override (defaults to $KUL_PYFMRI_ENV)
            pyfmri_env=$OPTARG
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
# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
    export MRTRIX_QUIET=1
fi

KUL_check_participant

# conda env existence check (mirrors KUL_dwiprep.sh's lore_sd pattern)
if ! conda env list | awk '{print $1}' | grep -qx "$pyfmri_env"; then
    echo "ERROR: conda env '$pyfmri_env' was not found (checked 'conda env list')." >&2
    echo "  Run the KUL_NIS installer's env-pyfmri section, or pass -c <env> to use a different one." >&2
    exit 1
fi
KUL_activate_conda_env "$pyfmri_env"

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
# step0_synthseg.sh derives its default paths from its own location, which after
# relocation into share/ resolves to the KUL_NIS repo rather than the study dir.
# The three above were already exported; these two were not, so the FreeSurfer
# and VBG lookups silently pointed at <KUL_NIS>/BIDS/derivatives/freesurfer --
# a path that cannot exist. Every Lausanne lookup therefore failed regardless of
# whether the subject actually had one, and reported "run KUL_FS_multiparc.sh
# first" even when multiparc had produced all five scales.
export RSFMRI_FS_DIR="$cwd/BIDS/derivatives/freesurfer"
export RSFMRI_VBG_DIR="$cwd/KUL_VBG"

if [ $KUL_DEBUG -gt 0 ]; then
    echo "pipeline_dir: $pipeline_dir"
    echo "RSFMRI_FMRIPREP_DIR: $RSFMRI_FMRIPREP_DIR"
    echo "RSFMRI_DENOISED_DIR: $RSFMRI_DENOISED_DIR"
    echo "RSFMRI_ANALYSIS_DIR: $RSFMRI_ANALYSIS_DIR"
    echo "RSFMRI_SYNTHSEG_DIR: $RSFMRI_SYNTHSEG_DIR"
    echo "RSFMRI_FS_DIR: $RSFMRI_FS_DIR"
    echo "RSFMRI_VBG_DIR: $RSFMRI_VBG_DIR"
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
#    Task runs are additionally stripped of their evoked response
#    (--task-signal remove) so that pooling them with genuine rest is
#    defensible: what is wanted here is the ongoing/background coupling, and
#    leaving a block paradigm in would drive spurious co-activation between
#    every region the task engages. Relying on the nuisance regressors to
#    remove it incidentally is not enough -- they only took ~58% of it on this
#    study's language run, leaving a substantial task-locked residual.
#    Rest runs have no design, so this is a no-op for them.
#    Set KUL_RSFMRI_TASK_SIGNAL=ignore to pool the task runs as-is instead.
_rsfmri_task_signal="${KUL_RSFMRI_TASK_SIGNAL:-remove}"
# Variant-keyed output: this product genuinely differs from the one the task
# melodic path needs, so they must not share a directory. Where a recipe does
# coincide the denoise script's own skip-if-exists check reuses it.
_rsfmri_denoised_variant="$RSFMRI_DENOISED_DIR/task-${_rsfmri_task_signal}"
export RSFMRI_DENOISED_DIR="$_rsfmri_denoised_variant"
task_in="$kul_main_dir/KUL_fmri_denoise.sh --fmriprep $RSFMRI_FMRIPREP_DIR --out $_rsfmri_denoised_variant --method nilearn --space MNI152NLin2009cAsym --sub sub-${participant} --task-signal $_rsfmri_task_signal"
KUL_task_exec $verbose_level "rsfMRI networks: denoise all BOLD runs (task-signal=$_rsfmri_task_signal)" "rsfmri_1_denoise" || { echo "ERROR: rsfMRI denoising failed for sub-${participant} — aborting (rsfMRI_networks.done will not be created)" >&2; exit 1; }

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

# Also link the actual per-patient statmaps (seed-based and RSN-FC), not just
# the PDF summary — these are the clinically relevant zscore/pval/loosig maps
# consumed by step5_report.py, and previously only lived under
# BIDS/derivatives/KUL_compute/rsfMRI_networks/analysis/{sba,rsn_fc}/, never
# copied into RESULTS/ like every other modality's output.
#
# The pipeline computes these on the MNI152NLin2009cAsym-space denoised BOLD
# (see the --space above), so — same as KUL_fmriproc_nilearn_new.sh does for
# its own MNI-space GLM output — warp back to native T1w before landing in
# RESULTS/, so they register onto the same T1w grid as SPM, Tracto and the
# Karawun labels instead of sitting in a different space unannounced.
mkdir -p "$globalresultsdir/sba" "$globalresultsdir/rsn_fc"
_rsfmri_mni2t1w=$(compgen -G "${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_*from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5" | head -1)
_rsfmri_t1w_ref=$(find "${cwd}/BIDS/sub-${participant}/anat/" -name "*_T1w.nii.gz" ! -name "*gadolinium*" 2>/dev/null | head -1)
# Recursive, not a flat glob: step1_sba.py writes Lausanne-based seeds (the
# Lip/Hand/Foot somatotopic set) into an sba/sub-<ID>/lausanne_scale3_sba/
# subdirectory -- see seed_subdir() in share/rsfmri_pipeline/src/config.py. A
# flat '*.nii.gz' would silently skip exactly those maps. The subdirectory is
# mirrored under RESULTS/ so the split survives into the deliverable, and the
# path is derived per-file rather than assumed.
_sba_root="$RSFMRI_ANALYSIS_DIR/sba/sub-${participant}"
_rsn_root="$RSFMRI_ANALYSIS_DIR/rsn_fc/sub-${participant}"
if [ -n "$_rsfmri_mni2t1w" ] && [ -n "$_rsfmri_t1w_ref" ]; then
    while IFS= read -r _statmap; do
        [ -f "$_statmap" ] || continue
        case "$_statmap" in
            */sba/*) _outdir="$globalresultsdir/sba/$(dirname "${_statmap#$_sba_root/}")" ;;
            *)       _outdir="$globalresultsdir/rsn_fc/$(dirname "${_statmap#$_rsn_root/}")" ;;
        esac
        _outdir="${_outdir%/.}"     # dirname gives '.' for files at the root
        mkdir -p "$_outdir"
        antsApplyTransforms -d 3 --float 1 \
            -i "$_statmap" -o "$_outdir/$(basename "$_statmap")" \
            -r "$_rsfmri_t1w_ref" -t "$_rsfmri_mni2t1w" -n Linear
    done < <(find "$_sba_root" "$_rsn_root" -name '*.nii.gz' 2>/dev/null)
else
    echo "WARNING: MNI-to-T1w transform or native T1w not found for sub-${participant} — copying rsfMRI statmaps as-is (still in MNI space, not warped)"
    # -a preserves the subdirectory split; '/.' copies contents, not the dir itself
    cp -a "$_sba_root/." "$globalresultsdir/sba/" 2>/dev/null
    cp -a "$_rsn_root/." "$globalresultsdir/rsn_fc/" 2>/dev/null
fi

# Only claim success if the run actually produced statmaps. Every step above
# can "succeed" while writing nothing: run_pipeline.py logs per-seed ERRORs and
# still exits 0, and step5_lite_report.py happily renders a report whose pages
# were all skipped for want of maps. That combination once wrote a .done on a
# run with 0 seeds and 0 RSNs, so the next invocation skipped straight past it
# ("rsfMRI networks already done") and the empty result looked deliberate.
# find is already recursive, so this counts the lausanne_scale3_sba/ subdirectory
# too -- important, since a profile whose seeds are ALL Lausanne-based would
# otherwise look like a zero-statmap failure and abort a successful run.
_n_statmaps=$(find "$_sba_root" "$_rsn_root" -name '*.nii.gz' 2>/dev/null | wc -l)
if [ "$_n_statmaps" -eq 0 ]; then
    echo "ERROR: rsfMRI networks produced no statmaps for sub-${participant} — rsfMRI_networks.done will NOT be created." >&2
    echo "       Check rsfmri_3_pipeline.log for per-seed errors (a missing denoised BOLD is the usual cause)." >&2
    exit 1
fi

touch "$done_flag"
echo "Done: rsfMRI networks for sub-${participant} ($_n_statmaps statmaps)"
