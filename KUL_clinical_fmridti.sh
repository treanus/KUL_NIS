#!/bin/bash
# Bash shell script to analyse clinical fMRI/DTI
#

# Requires matlab fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - KU Leuven, Translational MRI - ahmed.radwan@kuleuven.be
# 07/07/2026
version="2.0"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` is a batch analysis of clinical fMRI/DTI data

Usage:

  `basename $0` <OPT_ARGS>

Examples (paths/names fictitious — see docs/KUL_clinical_fmridti/KUL_clinical_fmridti.md for more):

  # 0. first time for this patient: scaffold a DICOM/ + study_config/ folder, then exit
  `basename $0` -p JaneDoe -s -t 1

  # 1. basic glioma work-up (type 1, default SPM engine), from a zip archive
  `basename $0` -p JaneDoe -d DICOM/JaneDoe.zip -n 32

  # 2. DBS case (type 3, manual mask) + nilearn GLM + lore_sd FOD + rsfMRI networks + tractometry
  `basename $0` -p M0012 -d ./DICOM/M0012 -t 3 \\
      -E nilearn -D run_dwiprep_lore_sd.txt -U -N -Q -n 32 -v 2

Required arguments:

     -p:  participant name

Optional arguments:

     -t:  processing type
        1: (DEFAULT) do automatic tumor intra-axial segmentation and vbg (tumor with T1w, cT1w, T2w and FLAIR)
        2: do automatic tumor extra-axial segmentation and vbg (tumor with T1w, cT1w, T2w and FLAIR)
        3: do vbg with manual mask (tumor but missing one of T1w, cT1w, T2w and FLAIR; 
                    put lesion.nii.gz in RESULTS/sub-{participant}/Lesion)
        4: dMRI/fMRI without glioma (cavernoma, epilepsy, etc... cT1w)
        5: dMRI for DBS of essential tremor (DRT tract)
        6: dMRI for DBS of Parkinson's disease (CSHD pathway)
        7: processing stream for DTI_ALPS
     -d:  dicom zip file (or directory)
     -s:  scaffold (make a default DICOM and study_config)
     -B:  delete the regenerable intermediates, then archive everything that is
          left into a password-protected ../Finished_<date>_sub-<p>_type<t>.7z
          and exit. Removed: fmriprep_work, BIDS/tmp_dcm2bids, dwiprep raw/dwi
          working dirs, the denoised BOLD variants and the SUSAN-smoothed GLM
          inputs (all rebuilt from fmriprep output + confounds, and several GB
          per subject). Kept: BIDS, fmriprep, dwiprep derivatives, RESULTS,
          REPORT, Karawun and KUL_LOG - plus dwi/geomcorr.mif, the topup+eddy
          corrected DWI, so that step never has to be re-run.
          Note this is destructive and one-way - re-running afterwards has to
          redo the denoise/smoothing steps.
     -r:  redo certain steps (program will ask)
     -R:  generate DICOMs for PACS and Karawun (run this AFTER reviewing figures)
        1: use cT1w as underlay
        2: use FLAIR as underlay
        3: use SWI as underlay
        4: use T1w as underlay
        5: use FGATIR as underlay
        6: use DIR as underlay
        7: use MP2RAGE (INV2) as underlay
     -O:  orientations to render, comma-separated (default: TRA,SAG,COR)
     -e:  add edge outline to SPM/Melodic overlays (dark blue contour at threshold boundary)
     -T:  fixed threshold for ALL SPM/Melodic overlays (default: auto = max/3 per map)
          if not given and running interactively, you will be prompted to enter
          one threshold per map (space-separated, matching the listed order)
     -a:  opacity of SPM/Melodic (fMRI) activation overlays (0=transparent, 1=opaque; default 0.7)
          lower values let underlying anatomy show through on figures and PACS DICOMs
     -D:  dwiprep config file to use from study_config/ (default: run_dwiprep.txt)
          use e.g. -D run_dwiprep_lore_sd.txt to run lore-sd based FOD estimation
     -S:  fMRI SUSAN smoothing FWHM in mm (default: adaptive = mean voxel size)
          e.g. -S 6 for 6mm FWHM, -S 8 for 8mm FWHM
     -P:  FWE-corrected p-value for Bizzi fMRI thresholding (default: 0.01)
          e.g. -P 0.01, -P 0.005, -P 0.001
     -n:  number of threads to use (default 48)
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)
     -X:  use FastSurfer instead of plain recon-all for the reconstruction step
          in types 4, 5, 6 (faster, requires GPU; default is FreeSurfer 8.2.0 recon-all)
     -f:  conda env to use instead of \$KUL_SCILPY_ENV (default 'scilpy') for
          FWT (automated tractography). You shouldn't normally need this —
          the KUL_NIS installer creates 'scilpy' with what FWT needs.
     -E:  fMRI GLM engine to use: spm or nilearn (default: spm)
          spm    : KUL_fmriproc_spm_new.sh    (MATLAB/SPM12, requires a MATLAB license)
          nilearn: KUL_fmriproc_nilearn_new.sh (python3 nilearn/nibabel/numpy/pandas, no MATLAB)
          both engines are auto-scheduled across $ncpu cores via their -c option
     -U:  EXPERIMENTAL opt-in: pass -U through to KUL_FWT_make_TCKs.sh so it prefers the
          rfa-modulated lore_sd FOD (rfa_modulated_fod_reg2T1w.mif, from KUL_dwiprep.sh),
          if found, over the plain lore_sd ODF. Without -U, tractography is unchanged.
     -Q:  opt-in: run KUL_FWT's per-bundle tractometry (-Q, along-tract scalar profiles).
          Off by default — adds substantial runtime (real per-bundle work across ~50
          bundle/hemisphere combinations, processed sequentially, no parallelism yet).
          Tractography/tracts themselves are generated either way.
     -W:  skip DSC perfusion processing (KUL_dsc_perfusion.sh). By default it
          runs automatically whenever a DSC series is present in BIDS/*/perf/,
          the same way fMRI and dMRI data are picked up.
     -N:  opt-in: run presurgical/eloquent-cortex rsfMRI network mapping
          (KUL_run_rsfMRI_networks.sh). Off by default.
     -C:  condition profile for -N, from share/rsfmri_pipeline/config/profiles.yaml
          (default: Presurgical)
     -y:  conda env to use instead of \$KUL_PYFMRI_ENV (default 'pyfMRI') for
          -N (rsfMRI network mapping) and -E nilearn (nilearn task-fMRI GLM).
          You shouldn't normally need this — the KUL_NIS installer creates
          'pyfMRI' with everything both steps need.
     -m:  conda env to use instead of \$KUL_DICOM_ENV (default 'KUL_dicom') for
          the DICOM generation in -R/-F (KUL_nii2dcm.py; needs SimpleITK,
          Pillow, numpy). If the env doesn't exist the step falls back to
          plain 'python3' with a warning, so this is only needed if you
          named the env differently. Create it with:
            mamba env create -f \$kul_main_dir/share/envs/KUL_dicom.yml

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ants_verbose=1
ncpu=48
bc=0
type=1
dwiprep_config_file="run_dwiprep.txt"
redo=0
results=0
make_dcm=0
verbose_level=1
dbs=0
scaffold=0
alps=0
msbp=0
multiparc=0
fwt=1
scilpy="$KUL_SCILPY_ENV"
orientations="TRA,SAG,COR"
spm_edge=0
spm_thresh_override=""
spm_opacity=0.7
smooth_fwhm=5
pfwe=0.01
use_fastsurfer=0
fmri_engine="spm"
use_rfa_mod_fod=0
run_fwt_tractometry=0
rsfmri_networks=0
skip_dsc=0
rsfmri_profile="Presurgical"
pyfmri_env_override=""
dicom_env_override=""
declare -A spm_thresh_map=()

# Set required options
p_flag=0
d_flag=0


if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:d:n:v:R:F:O:a:f:T:D:S:P:E:NC:y:m:XBrseUQW" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        t) #type
			type=$OPTARG
		;;
        W) #skip DSC perfusion
			skip_dsc=1
		;;
        d) #dicomzip
			dicomzip=$OPTARG
            d_flag=1
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
		B) #backup&clean
			bc=1
		;;
        r) #redo
			redo=1
		;;
        s) #scaffold
			scaffold=1
		;;
        R) #generate DICOMs for PACS and Karawun (user-explicit, run after reviewing figures)
			results=$OPTARG
            make_dcm=1
		;;
        F) #internal: generate figures/screenshots only (no DICOMs)
            results=$OPTARG
        ;;
        O) #orientations
            orientations=$OPTARG
        ;;
        e) #spm edge outline
            spm_edge=1
        ;;
        T) #spm threshold override
            spm_thresh_override=$OPTARG
        ;;
        a) #spm/melodic overlay opacity
            spm_opacity=$OPTARG
        ;;
        D) #dwiprep config file
            dwiprep_config_file=$OPTARG
        ;;
        S) #fMRI smoothing FWHM
            smooth_fwhm=$OPTARG
        ;;
        P) #FWE p-value for Bizzi thresholding
            pfwe=$OPTARG
        ;;
        v) #verbose
            verbose_level=$OPTARG
		;;
        X) #use FastSurfer instead of plain recon-all for types 4/5/6
            use_fastsurfer=1
        ;;
        f) # conda env override for FWT (default: $KUL_SCILPY_ENV)
            scilpy=$OPTARG
        ;;
        E) # fMRI GLM engine: spm or nilearn
            fmri_engine=$OPTARG
        ;;
        U) # EXPERIMENTAL opt-in: prefer the rfa-modulated lore_sd FOD in KUL_FWT, if found
            use_rfa_mod_fod=1
        ;;
        Q) # opt-in: run KUL_FWT's per-bundle tractometry (-Q). Off by default —
           # adds substantial runtime (real per-bundle work, ~50 bundle/hemisphere
           # combinations processed sequentially with no bundle-level parallelism yet).
            run_fwt_tractometry=1
        ;;
        N) # opt-in: run presurgical/eloquent-cortex rsfMRI network mapping
           # (KUL_run_rsfMRI_networks.sh). Off by default.
            rsfmri_networks=1
        ;;
        C) # condition profile for -N (default: Presurgical)
            rsfmri_profile=$OPTARG
        ;;
        y) # conda env override for -N and -E nilearn (default: $KUL_PYFMRI_ENV)
            pyfmri_env_override=$OPTARG
        ;;
        m) # conda env override for the -R/-F DICOM step (default: $KUL_DICOM_ENV)
            dicom_env_override=$OPTARG
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
# participant name is compulsory
if [ $p_flag -eq 0 ] ; then
	echo
	echo "Option -p is required: give the BIDS name of the participant." >&2
	echo
	exit 2
fi

if [ "$fmri_engine" != "spm" ] && [ "$fmri_engine" != "nilearn" ]; then
	echo
	echo "Option -E must be 'spm' or 'nilearn' (got '$fmri_engine')." >&2
	echo
	exit 2
fi
echo " fMRI GLM engine set to: $fmri_engine"

function KUL_scaffold {

    echo "Making scaffold for clinical_sub-${participant}_type${type}"
    mkdir -p $cwd/clinical_sub-${participant}_type${type}/DICOM
    mkdir -p $cwd/clinical_sub-${participant}_type${type}/study_config
    rm -fr $cwd/clinical_sub-${participant}_type${type}/study_config/*
    if [ $type -lt 5 ]; then
        echo "Setting up for a tumor/epilepsy/... patient (type: $type)"
        cp ${kul_main_dir}/study_config/clinical_fmri_dmri/* $cwd/clinical_sub-${participant}_type${type}/study_config
    elif [ $type -eq 5 ]; then
        echo "Setting up for a DBS patient (type: $type)"
        cp ${kul_main_dir}/study_config/clinical_dmri_dbs_drt/* $cwd/clinical_sub-${participant}_type${type}/study_config
    elif [ $type -eq 6 ]; then
        echo "Setting up for a DBS patient (type: $type)"
        cp ${kul_main_dir}/study_config/clinical_dmri_dbs_hdp/* $cwd/clinical_sub-${participant}_type${type}/study_config
    elif [ $type -eq 7 ]; then
        echo "Setting up for DTI-ALPS processing (type: $type)"
        cp ${kul_main_dir}/study_config/DTI_ALPS_proc/* $cwd/clinical_sub-${participant}_type${type}/study_config
    fi

    exit 0

}

# -R and -F render and export what is already in RESULTS/; they run no
# preprocessing at all. Every pre-flight check below guards a preprocessing step
# those runs never reach, and each one either scaffolds a study directory that
# was not asked for or exits 2 outright -- so a tree holding nothing but
# RESULTS/sub-<p>/ could not be exported to PACS without first being dressed up
# as a full study. Skip them for export-only runs, so the PACS workflow can be
# driven on its own the way KUL_karawun_prepare.sh can.
#
# Deliberately narrow: this only skips checks, never any processing decision,
# and a normal run (results=0) still fails fast on all of them exactly as before.
if [ $results -gt 0 ]; then
    export_only=1
else
    export_only=0
fi

# Scaffold: explicit (-s) or automatic (first run for this patient, no study_config/ yet).
# Must run before the pre-flight -D check below -- that check unconditionally looks for
# study_config/${dwiprep_config_file} (default "run_dwiprep.txt") and exits with an error
# if study_config/ doesn't exist yet, which would otherwise always pre-empt the automatic
# scaffold on a genuinely fresh patient folder.
if [ $scaffold -eq 1 ]; then
    KUL_scaffold
fi
if [ $export_only -eq 0 ] && [ ! -d $cwd/study_config ]; then
    KUL_scaffold
fi


# Pre-flight check of the -D dwiprep config: catch a lore_sd env
# misconfiguration here, before fmriprep/dwiprep are launched, instead of
# failing deep inside KUL_preproc_all.sh after other pipeline steps (and
# their downstream dependents: VBG, dwiprep_MNI, FWT) have already run.
if [ $export_only -eq 0 ]; then
    if [ ! -f study_config/${dwiprep_config_file} ]; then
        echo
        echo "ERROR: dwiprep config file study_config/${dwiprep_config_file} (given with -D) does not exist." >&2
        echo
        exit 2
    fi
    dwiprep_do_check=$(grep -E "^do_dwiprep:" study_config/${dwiprep_config_file} | grep -v \# | cut -d':' -f2- | tr -d '\r' | sed 's/^ *//;s/ *$//')
    dwiprep_options_check=$(grep -E "^dwiprep_options:" study_config/${dwiprep_config_file} | grep -v \# | cut -d':' -f2-)
    loresd_env_check=$(grep -E "^loresd_env:" study_config/${dwiprep_config_file} | grep -v \# | cut -d':' -f2- | tr -d '\r' | sed 's/^ *//;s/ *$//')
    # blank in the config -> defaults to $KUL_LORESD_ENV (mirrors KUL_preproc_all.sh)
    [ -z "$loresd_env_check" ] && loresd_env_check="$KUL_LORESD_ENV"
    if [ "$dwiprep_do_check" = "1" ] && [[ "$dwiprep_options_check" == *"lore_sd"* ]] && ! conda env list | awk '{print $1}' | grep -qx "$loresd_env_check"; then
        echo
        echo "ERROR: study_config/${dwiprep_config_file} requests 'lore_sd', but conda env '$loresd_env_check' was not found (checked 'conda env list')." >&2
        echo "  Run the KUL_NIS installer's env-lore-sd section, or set loresd_env: <env> in the config to use a different one." >&2
        echo "  Available conda envs:" >&2
        conda env list 2>/dev/null | tail -n +3 | awk '{print "    "$1}' >&2
        echo
        exit 2
    fi
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

# Determine what to process depending on patient lesion type
if [ $type -eq 1 ]; then
    hdglio=1; vbg=1; multiparc=0; fwt=1; alps=0
elif [ $type -eq 2 ]; then
    hdglio=1; vbg=2; multiparc=0; fwt=1; alps=0
elif [ $type -eq 3 ]; then
    hdglio=0; vbg=3; multiparc=0; fwt=1; alps=0
elif [ $type -eq 4 ]; then
    hdglio=0; vbg=0; multiparc=1; fwt=1; alps=0
elif [ $type -eq 5 ]; then
    hdglio=0; vbg=0; dbs=1; multiparc=1; fwt=1; alps=0
elif [ $type -eq 6 ]; then
    hdglio=0; vbg=0; dbs=2; multiparc=1; fwt=1; alps=0
elif [ $type -eq 7 ]; then
    hdglio=0; vbg=0; dbs=0; multiparc=0; fwt=0; alps=1
fi

figs=0

# fwt (automated tractography) needs the scilpy conda env; defaults to
# $KUL_SCILPY_ENV (the fixed name the installer creates it under), override
# with -f if you need a differently-named env
if [ $export_only -eq 0 ] && [ $fwt -eq 1 ] && ! conda env list | awk '{print $1}' | grep -qx "$scilpy" ; then
	echo
	echo "ERROR: conda env '$scilpy' (for FWT) was not found (checked 'conda env list')." >&2
	echo "  Run the KUL_NIS installer's env-scilpy section, or pass -f <env> to use a different one." >&2
	echo
	exit 2
fi

if [ -n "$scilpy" ] ; then
    echo " Scilpy conda env name for FWT set to $scilpy"
fi

# GLOBAL defs
globalresultsdir=$cwd/RESULTS/sub-$participant
derivativesdir=${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}
# KUL_VBG writes very deep subject/session-nested paths internally
# (output_VBG/sub-X/sub-X_FS_output/sub-X/...). Keep it flat and top-level,
# a sibling of dwiprep/fmriprep/BIDS, matching how those tools are already
# laid out (they're not nested under BIDS/derivatives/ either) and giving
# the most path-length headroom against FreeSurfer 8.2.0's mris_register
# buffer overflow on long SUBJECTS_DIR paths (see
# VBG_FS820_mris_register_bugreport.md — the actual fix is the symlink
# workaround in KUL_VBG.sh itself; this is extra margin, not the fix).
vbg_dir=${cwd}/KUL_VBG

# NOTE: these three live here, above the `if [ $results -gt 0 ]` block below,
# rather than with the other functions further down. That block runs at top
# level and ends in `exit`, so any function defined after it has simply not been
# executed yet when it runs -- bash registers a function when its definition
# runs, not when the file is parsed. Calling KUL_resolve_lesion from the PACS
# code inside that block while it was defined further down failed silently:
# "command not found" on stderr, the `if` took the else branch, exit status 0,
# and the lesion was just missing from the export. `bash -n` does not catch it.

# Resolve the lesion mask this participant actually has, whatever produced it.
# Types 1/2 get it from KUL_anat_segment_tumor; type 3 is the manual mask the
# user drops in themselves. Echoes the path, or nothing if there is none.
function KUL_resolve_lesion {
    local _c
    for _c in "$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" \
              "$globalresultsdir/Lesion/lesion.nii.gz"; do
        if [ -f "$_c" ]; then
            echo "$_c"
            return 0
        fi
    done
    return 1
}

# Put the lesion alongside the other T1w-space volumes in RESULTS/*/Anat, so it
# sits with the underlays that feed the figures and the PACS export rather than
# only in its own Lesion/ folder.
function KUL_copy_lesion_to_anat {
    local _lesion
    _lesion=$(KUL_resolve_lesion) || { echo "No lesion mask to copy to Anat"; return 0; }

    local _dst="$globalresultsdir/Anat/sub-${participant}_lesion.nii.gz"
    if [ ! -f "$_dst" ] || [ "$_lesion" -nt "$_dst" ]; then
        # binarised on the way over: the hd-glio output carries per-class values,
        # and anything overlaying this wants a clean 0/1 mask
        mrcalc "$_lesion" 0 -gt "$_dst" -force -quiet && \
            echo "Copied $(basename "$_lesion") to Anat/ as $(basename "$_dst")"
    fi
}

# Create RESULTS/sub-*/PACS_input/{overlays,series_quantitative} and the README
# that explains them, and set _dropdir/_drop_thr/_drop_quant for the caller.
#
# Called from two places, which is why it is a function rather than inline:
#   - the end of a normal pipeline run, so the folders exist for the operator to
#     review and copy into BEFORE the first -R. They used to be created only
#     inside the -R block, i.e. only once it was already too late to put
#     anything in them;
#   - the -R/-F block itself, so a standalone export still works against a tree
#     that never saw a full pipeline run.
function KUL_make_pacs_dropdirs {
    _dropdir="$globalresultsdir/PACS_input"
    _drop_thr="$_dropdir/overlays"
    _drop_quant="$_dropdir/series_quantitative"
    mkdir -p "$_drop_thr" "$_drop_quant"

    # Echo a README once so the folders explain themselves on the filesystem.
    if [ ! -f "$_dropdir/README.txt" ]; then
        cat > "$_dropdir/README.txt" <<'DROPREADME'
Drop NIfTI files here to export them to PACS with KUL_clinical_fmridti.sh -R.

  overlays/              colour overlay fused on the anatomical, for review.

                         Colour windowing is always computed from the map's
                         own robust (2-98%) range, so a perfusion map with
                         values in the thousands and an fMRI t-map with a
                         narrow range both render legibly. Nothing to set.

                         Thresholding, i.e. which voxels are drawn at all:
                           default          automatic, max/3 (suits stat maps)
                           <name>.thresh    a file holding one number, to pin
                                            a threshold (e.g. 1.75 for nrCBV)
                           <name>.thresh    holding the word "none" to draw
                                            the whole map with a colourbar --
                                            use this for rCBV, ALFF, ReHo
                           -T <value>       overrides everything, all maps

  series_quantitative/   a standalone measurable DICOM series (no overlay, no
                         rendering, no threshold). Values carry
                         RescaleSlope/Intercept so an ROI on PACS reads true
                         units. Name a file <name>.label.nii.gz for an integer
                         label map, which is written with no rescaling so each
                         label keeps its value.

While either folder holds a file, the automatic discovery of
RESULTS/sub-*/{SPM,Melodic,Perfusion,Lesion} is skipped. Empty them to get
the automatic behaviour back.

Note that automatic overlay discovery only looks at SPM/, Melodic/, the lesion,
and *only* nrCBV_corrected/nrCBF from Perfusion/. Anything else -- a raw
rCBV_corrected, an ADC map, a map from outside this pipeline -- has to be
dropped in here to reach PACS.

Files must already be registered to the anatomical (same scanner coordinate
frame as the underlay) — they inherit the donor's Frame of Reference.
DROPREADME
    fi
}




# The BACKUP and clean option
if [ $bc -eq 1 ]; then
    # clean some stuff
    clean_dwiprep="./dwiprep/sub-${participant}/sub-${participant}/*dwifsl*tmp* \
        ./dwiprep/sub-${participant}/sub-${participant}/raw \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi_orig*"
    clean_other="./fmriprep_work* \
        ./BIDS/tmp_dcm2bids"

    # dwi/ is scratch (degibbs, biascorr, noiselevel, nonbzeros...) with one
    # exception: geomcorr.mif is what dwifslpreproc produced, i.e. the output of
    # topup+eddy -- by far the most expensive step here, and whose own scratch
    # (*dwifsl*tmp*, ~19 GB) is deleted just above. Keeping it means bias
    # correction and everything downstream can be redone without paying for eddy
    # again. dwi_preproced.mif one level up is the POST-bias-correction volume,
    # so it is not a substitute. geomcorr_grad_checked.b is the gradient table
    # dwigradcheck corrected post-eddy and has to travel with it.
    clean_dwi_dir="./dwiprep/sub-${participant}/sub-${participant}/dwi"

    # Denoised BOLD copies and the SUSAN-smoothed GLM inputs. Both are whole 4D
    # series regenerated from the fmriprep output and the confounds TSVs, so they
    # are pure intermediates -- but they are large, and now multiplied: the task
    # melodic, resting-state and pseudo-rest paths each keep their own variant
    # (~6 GB combined on a 3-run subject, plus ~2 GB of smoothed copies). Left in
    # place they were silently archived into the .7z.
    clean_fmri="./BIDS/derivatives/KUL_compute/sub-${participant}/FSL_melodic/denoised \
        ./BIDS/derivatives/KUL_compute/rsfMRI_networks/denoised \
        ./BIDS/derivatives/KUL_compute/sub-${participant}/SPM/fmridata"

    echo "Removing regenerable intermediates before archiving (keeping geomcorr.mif):"
    { du -shc $clean_dwiprep $clean_other $clean_fmri 2>/dev/null | tail -1
      find $clean_dwi_dir -type f -not -name 'geomcorr*' -print0 2>/dev/null \
        | du -shc --files0-from=- 2>/dev/null | tail -1
    } | awk '{print "  " $0}'
    rm -fr $clean_dwiprep $clean_other $clean_fmri
    find $clean_dwi_dir -type f -not -name 'geomcorr*' -delete 2>/dev/null

    while true; do
        read -s -p "Give a password to encrypt the backup with: " password
        echo
        read -s -p "Password (again): " password2
        echo
        [ "$password" = "$password2" ] && break
        echo "Please try again"
    done    
    today=$(date +%Y_%m_%d)
    7z a -bd -y -p${password} -mhe=on -mx=7 ../Finished_${today}_sub-${participant}_type${type}.7z * 

    exit 0
fi

# The make RESULTS option
if [ $results -gt 0 ];then

    mrview_tracts[0]="Tract-csd_CST_LT"
    mrview_rgb[0]="0.678,0.847,0.902"
    mrview_tracts[1]="Tract-csd_CST_RT"
    mrview_rgb[1]="0,0,1"
    mrview_tracts[2]="Tract-csd_AF_all_LT"
    mrview_rgb[2]="1,0,0"
    mrview_tracts[3]="Tract-csd_AF_all_RT"
    mrview_rgb[3]="0,1,0"
    mrview_tracts[4]="Tract-csd_CCing_LT"
    mrview_rgb[4]="1,1,0"
    mrview_tracts[5]="Tract-csd_CCing_RT"
    mrview_rgb[5]="1,0.647,0"
    mrview_tracts[6]="Tract-csd_TCing_LT"
    mrview_rgb[6]="1,1,0"
    mrview_tracts[7]="Tract-csd_TCing_RT"
    mrview_rgb[7]="1,0.646,0"
    mrview_tracts[8]="Tract-csd_FAT_LT"
    mrview_rgb[8]="1,0.647,0"
    mrview_tracts[9]="Tract-csd_FAT_RT"
    mrview_rgb[9]="1,1,0"
    mrview_tracts[10]="Tract-csd_ILF_LT"
    mrview_rgb[10]="0,0,1"
    mrview_tracts[11]="Tract-csd_ILF_RT"
    mrview_rgb[11]="0.678,0.847,0.784"
    mrview_tracts[12]="Tract-csd_IFOF_LT"
    mrview_rgb[12]="0.75,0.25,0.75"
    mrview_tracts[13]="Tract-csd_IFOF_RT"
    mrview_rgb[13]="1,0.75,0.8"
    mrview_tracts[14]="Tract-csd_UF_LT"
    mrview_rgb[14]="0,0.784,0"
    mrview_tracts[15]="Tract-csd_UF_RT"
    mrview_rgb[15]="0.784,0,0"
    mrview_tracts[16]="Tract-csd_OR_occlobe_LT"
    mrview_rgb[16]="0.2,0.784,0.4"
    mrview_tracts[17]="Tract-csd_OR_occlobe_RT"
    mrview_rgb[17]="0.784,0.2,0.4"
    mrview_tracts[18]="Tract-csd_MdLF_LT"
    mrview_rgb[18]="0.6,0.784,0.04"
    mrview_tracts[19]="Tract-csd_MdLF_RT"
    mrview_rgb[19]="0.784,0.6,0.04"
    mrview_tracts[20]="Tract-csd_ML_LT"
    mrview_rgb[20]="0,0.5,0.04"
    mrview_tracts[21]="Tract-csd_ML_RT"
    mrview_rgb[21]="0,0.5,0.5"

    mrview_tracts[22]="Tract-csd_DRT_LT"
    mrview_rgb[22]="1,0,0.23"
    mrview_tracts[23]="Tract-csd_DRT_RT"
    mrview_rgb[23]="0.23,1,0"

    mrview_tracts[24]="Tract-csd_CSHDP_LT"
    mrview_rgb[24]="0.23,0.12,0"
    mrview_tracts[25]="Tract-csd_CSHDP_RT"
    mrview_rgb[25]="1,0.12,0.20"

    mrview_tracts[26]="Tract-csd_SLF_all_LT"
    mrview_rgb[26]="0.1,0.6,1"
    mrview_tracts[27]="Tract-csd_SLF_all_RT"
    mrview_rgb[27]="1,0.6,0.1"
    mrview_tracts[28]="Tract-csd_SLF_I_LT"
    mrview_rgb[28]="0.5,0.9,0.9"
    mrview_tracts[29]="Tract-csd_SLF_I_RT"
    mrview_rgb[29]="0.9,0.5,0.9"
    mrview_tracts[30]="Tract-csd_SLF_II_LT"
    mrview_rgb[30]="0.3,0.3,0.9"
    mrview_tracts[31]="Tract-csd_SLF_II_RT"
    mrview_rgb[31]="0.9,0.3,0.3"
    mrview_tracts[32]="Tract-csd_SLF_III_LT"
    mrview_rgb[32]="0.1,0.8,0.3"
    mrview_tracts[33]="Tract-csd_SLF_III_RT"
    mrview_rgb[33]="0.8,0.1,0.5"

    mrview_tracts[34]="Tract-csd_Ant_Comm"
    mrview_rgb[34]="0.5,0,0.5"
    mrview_tracts[35]="Tract-csd_Post_Comm"
    mrview_rgb[35]="0.8,0.4,0.8"
    mrview_tracts[36]="Tract-csd_CC_PreF_Comm"
    mrview_rgb[36]="0.2,0.8,0.8"
    mrview_tracts[37]="Tract-csd_CC_Motor_Comm"
    mrview_rgb[37]="0.8,0.4,0"
    mrview_tracts[38]="Tract-csd_CC_PMandSM_Comm"
    mrview_rgb[38]="1,0.5,0"
    mrview_tracts[39]="Tract-csd_CC_Sensory_Comm"
    mrview_rgb[39]="0.4,0.8,0.4"
    mrview_tracts[40]="Tract-csd_CC_Parietal_Comm"
    mrview_rgb[40]="0.4,0,0.8"
    mrview_tracts[41]="Tract-csd_CC_Temporal_Comm"
    mrview_rgb[41]="0.8,0,0.4"
    mrview_tracts[42]="Tract-csd_CC_Occipital_Comm"
    mrview_rgb[42]="0,0.6,0.8"

    ntracts_paired=34   # indices 0-33: lateralized LT/RT pairs, step 2
    ntracts_total=42    # index of last commissural tract (indices 34-42), step 1

    # Drop folders first, before anything that can bail out. A normal pipeline
    # run creates these at its end, but an export-only run (a tree holding just
    # RESULTS/) has never had that step, and the operator needs somewhere to copy
    # maps into. Creating them here means even a run that dies at the mrview
    # pre-flight below leaves a usable, self-describing tree behind, rather than
    # forcing a throwaway -R just to make the folders appear.
    KUL_make_pacs_dropdirs

    # Sync FWT output to RESULTS/Tracto before generating screenshots
    # Tract maps are resampled to T1w resolution so they can be compared
    # directly with anatomical images and lesion masks.
    mkdir -p $globalresultsdir/Tracto
    rm -fr $globalresultsdir/Tracto/*
    _t1w_ref="$globalresultsdir/Anat/T1w.nii.gz"
    [[ ! -f "$_t1w_ref" ]] && _t1w_ref="fmriprep/sub-${participant}/anat/sub-${participant}_desc-preproc_T1w.nii.gz"
    for tck_outdir in "$derivativesdir/FWT/sub-${participant}_TCKs_output"/*_output; do
        [ -d "$tck_outdir" ] || continue
        tract_name=$(basename "$tck_outdir" _output)
        fin_tck="${tck_outdir}/${tract_name}_fin_BT_iFOD2.tck"
        fin_map="${tck_outdir}/${tract_name}_fin_map_BT_iFOD2.nii.gz"
        _out_map="$globalresultsdir/Tracto/Tract-csd_${tract_name}.nii.gz"
        if [ -f "$fin_tck" ]; then
            cp "$fin_tck" "$globalresultsdir/Tracto/Tract-csd_${tract_name}.tck"
            [ -f "$fin_map" ] && mrgrid "$fin_map" regrid -template "$_t1w_ref" -interp linear "$_out_map" -force -quiet
        else
            use_tck=$(ls "${tck_outdir}/${tract_name}_filt"*"_BT_iFOD2.tck" 2>/dev/null | grep -v "_inMNI" | sort -V | tail -1)
            use_map=$(ls "${tck_outdir}/${tract_name}_filt"*"_map_BT_iFOD2.nii.gz" 2>/dev/null | grep -v "_inMNI" | sort -V | tail -1)
            if [ -n "$use_tck" ]; then
                echo "  ${tract_name}: fin not found, using $(basename $use_tck)"
                cp "$use_tck" "$globalresultsdir/Tracto/Tract-csd_${tract_name}.tck"
                [ -n "$use_map" ] && mrgrid "$use_map" regrid -template "$_t1w_ref" -interp linear "$_out_map" -force -quiet
            fi
        fi
    done

    #echo "ntracts: $ntracts"
    result_type=0

    if [ $results -eq 1 ]; then

        underlay=$globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_cT1w"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_cT1w"

    elif [ $results -eq 2 ]; then

        underlay=$globalresultsdir/Anat/FLAIR_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_FLAIR"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_FLAIR"

    elif [ $results -eq 3 ]; then

        underlay=$globalresultsdir/Anat/SWI_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_SWI"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_SWI"

    elif [ $results -eq 5 ]; then

        underlay=$globalresultsdir/Anat/FGATIR_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_FGATIR"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_FGATIR"

    elif [ $results -eq 6 ]; then

        underlay=$globalresultsdir/Anat/DIR_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_DIR"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_DIR"

    elif [ $results -eq 7 ]; then

        underlay=$globalresultsdir/Anat/MP2RAGE_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_MP2RAGE"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_MP2RAGE"

    else

        underlay=$globalresultsdir/Anat/T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_T1w"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_T1w"

    fi

    mrview_resolution=512

    # --- Python interpreter for KUL_nii2dcm.py ------------------------------
    # Resolved once here and reused by every conversion below. We call the env's
    # python by ABSOLUTE PATH rather than 'conda activate' on purpose: this block
    # also runs mrview/mrinfo/mrstats/mrcalc from PATH, and activating an env
    # mid-script reorders PATH for all of them. Absolute-path invocation is just
    # as automatic for the user and touches nothing else.
    # Falls back to plain python3 (the historical behaviour) if the env is
    # missing, so hosts that never created it keep working.
    _dicom_env="${dicom_env_override:-$KUL_DICOM_ENV}"
    _nii2dcm_py=""
    if command -v conda >/dev/null 2>&1; then
        _conda_base=$(conda info --base 2>/dev/null)
        if [ -n "$_conda_base" ] && [ -x "$_conda_base/envs/$_dicom_env/bin/python" ]; then
            _nii2dcm_py="$_conda_base/envs/$_dicom_env/bin/python"
        fi
    fi
    # SOFTWARE_ROOT may not be exported; probe the usual roots too (same idea as
    # the HD-GLIO lookup in KUL_anat_segment_tumor.sh).
    if [ -z "$_nii2dcm_py" ]; then
        for _sw in "${SOFTWARE_ROOT:-}" /opt/kul_software /usr/local/KUL_apps; do
            [ -z "$_sw" ] && continue
            if [ -x "$_sw/miniforge3/envs/$_dicom_env/bin/python" ]; then
                _nii2dcm_py="$_sw/miniforge3/envs/$_dicom_env/bin/python"
                break
            fi
        done
    fi
    if [ -n "$_nii2dcm_py" ]; then
        echo "Using conda env '$_dicom_env' for DICOM generation: $_nii2dcm_py"
    else
        _nii2dcm_py="python3"
        echo "WARNING: conda env '$_dicom_env' not found — falling back to plain 'python3'."
        echo "         Create it with: mamba env create -f $kul_main_dir/share/envs/KUL_dicom.yml"
        echo "         or point at an existing env with -m <envname>."
    fi
    # Verify the interpreter can actually import what KUL_nii2dcm.py needs, so a
    # broken env is reported here rather than as N identical failures later.
    if ! "$_nii2dcm_py" -c "import SimpleITK, PIL, numpy" >/dev/null 2>&1; then
        echo "ERROR: '$_nii2dcm_py' cannot import SimpleITK/Pillow/numpy."
        echo "       DICOM generation will fail. Fix the env (see -m in the usage) and re-run."
    fi
    _nii2dcm_script="$kul_main_dir/KUL_nii2dcm.py"
    if [ ! -f "$_nii2dcm_script" ]; then
        echo "ERROR: $_nii2dcm_script not found — DICOM generation will fail."
    fi

    # Donor DICOM: used to copy patient/study metadata into PACS DICOMs.
    # Search Karawun first, then RESULTS/sub-*/DICOM (case-insensitive, .dcm and .ima).
    donor_dcm=$(find "Karawun/sub-${participant}/DICOM" \( -iname "*.dcm" -o -iname "*.ima" \) -type f 2>/dev/null | sort | head -1)
    if [ -z "$donor_dcm" ]; then
        donor_dcm=$(find "$globalresultsdir/DICOM" \( -iname "*.dcm" -o -iname "*.ima" \) -type f 2>/dev/null | sort | head -1)
    fi
    # Philips exports are often extensionless, so fall back to any file in
    # these known DICOM-only directories (same trust model already used for
    # the SmartBrain/Localizer donor copy above), excluding non-DICOM
    # housekeeping files that sometimes ship alongside a DICOM export.
    if [ -z "$donor_dcm" ]; then
        donor_dcm=$(find "Karawun/sub-${participant}/DICOM" -type f -not -name ".*" -not -iname "DICOMDIR" 2>/dev/null | sort | head -1)
    fi
    if [ -z "$donor_dcm" ]; then
        donor_dcm=$(find "$globalresultsdir/DICOM" -type f -not -name ".*" -not -iname "DICOMDIR" 2>/dev/null | sort | head -1)
    fi
    if [ $make_dcm -eq 1 ] && [ -z "$donor_dcm" ]; then
        echo ""
        echo "ERROR: No donor DICOM found in Karawun/sub-${participant}/DICOM/ or $globalresultsdir/DICOM/"
        echo "       Copy one (a single file is enough) into either folder before running -R again:"
        echo "         one slice from a high-resolution anatomical series (T1w, FLAIR, T2, ...)"
        echo "         same donor is used for both PACS and Karawun/Brainlab output"
        echo "       DICOM output will be skipped for this run."
        echo ""
    elif [ -n "$donor_dcm" ]; then
        echo "Using donor DICOM: $donor_dcm"
    fi

    # --- SeriesNumber base, anchored to the donor ---------------------------
    # Exported series get  D*100 + map_index*10 + orientation_index , where D is
    # the donor's own SeriesNumber. That keeps our series adjacent to the study
    # they belong to and ordered, while landing well above the range scanners
    # actually use (Philips 101/201/.../1501, Siemens 1-30), so we never
    # renumber over a real acquisition. Previously SeriesNumber was left empty
    # on every series, which most PACS sort last or arbitrarily.
    _donor_series_base=90
    if [ -n "$donor_dcm" ]; then
        _dsn=$("$_nii2dcm_py" -c "
import sys, SimpleITK as sitk
try:
    r = sitk.ImageFileReader(); r.SetFileName(sys.argv[1]); r.ReadImageInformation()
    for k in ('0020|0011', '0020|0011 '):
        if r.HasMetaDataKey(k.strip()):
            print(int(str(r.GetMetaData(k.strip())).strip())); break
except Exception:
    pass
" "$donor_dcm" 2>/dev/null)
        if [[ "$_dsn" =~ ^[0-9]+$ ]] && [ "$_dsn" -gt 0 ]; then
            _donor_series_base=$_dsn
        else
            echo "Note: donor DICOM has no usable SeriesNumber (0020,0011); using base ${_donor_series_base} for exported series."
        fi
    fi
    echo "SeriesNumber base (from donor): ${_donor_series_base} -> series numbered from $(( _donor_series_base * 100 + 10 ))"

    # Fixed orientation->digit table. Deliberately NOT the position in $orientations:
    # -O SAG,TRA would otherwise silently renumber everything between runs.
    declare -A _ori_index=( [TRA]=0 [COR]=1 [SAG]=2 )

    # map_index comes from a registry file written ONLY by the parent shell.
    # The render functions get backgrounded, so if a child could append here two
    # concurrent renders would race and could claim the same index. Registration
    # therefore happens in the parent right before a name is dispatched, and the
    # children only ever look up.
    #
    # The registry PERSISTS per subject, and is never rewritten -- only appended.
    # A per-run temp file would restart map_index at 1 every time, so exporting
    # fMRI today and perfusion tomorrow, or re-running against a second underlay
    # (-R 4 then -R 2), would issue SeriesNumbers already used by series sitting
    # in the same PACS study. Delete this file only if you want the whole
    # subject renumbered from scratch.
    mkdir -p "$globalresultsdir/PACS"
    _pacs_index_file="$globalresultsdir/PACS/.series_index"
    [ -f "$_pacs_index_file" ] || : > "$_pacs_index_file"
    # Parent-only: give $1 a stable 1-based index for the rest of the run.
    _pacs_register() {
        grep -q -x -F "$1" "$_pacs_index_file" 2>/dev/null || printf '%s\n' "$1" >> "$_pacs_index_file"
    }
    # Read-only: SeriesNumber for a registered name + orientation. Returns empty
    # for an unregistered name, in which case the caller simply omits -n and gets
    # the historical (empty SeriesNumber) behaviour rather than a colliding one.
    _pacs_series_number() {
        local _idx _oi
        _idx=$(grep -n -x -F "$1" "$_pacs_index_file" 2>/dev/null | head -1 | cut -d: -f1)
        if [ -z "$_idx" ]; then
            echo "WARNING: '$1' was not registered for SeriesNumber allocation" >&2
            return 0
        fi
        _oi=${_ori_index[$2]:-0}
        printf '%s' "$(( _donor_series_base * 100 + _idx * 10 + _oi ))"
    }

    # Compute correct PixelSpacing for each orientation from the underlay geometry.
    # mrview -size N,N renders the full image FOV into exactly N pixels per side;
    # the larger in-plane physical dimension spans mrview_resolution pixels.
    _dims=($(mrinfo $underlay -size))
    _vox=($(mrinfo $underlay -spacing))
    px_tra=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[1]}*${_vox[1]})/$mrview_resolution)")
    px_sag=$(python3 -c "print(max(${_dims[1]}*${_vox[1]},${_dims[2]}*${_vox[2]})/$mrview_resolution)")
    px_cor=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[2]}*${_vox[2]})/$mrview_resolution)")

    # How many mrview instances run at once, and Mesa llvmpipe threads for each.
    # Total cores used = mrview_threads * _render_par. This now governs the
    # fMRI/Melodic/Clinical renders as well as the tract renders; those used to
    # run strictly serially while still being sized as if 3 were in flight, so
    # they used a third of the box.
    # Override with KUL_RENDER_PAR=1 to render serially. Worth doing on hosts
    # where concurrent mrview instances contend on a global SysV semaphore:
    # one instance runs and the others block on it with zero CPU until the
    # timeout kills them, so 3-way parallelism ends up slower than serial and
    # can look like a hang. Symptom: `cat /proc/<pid>/wchan` on the stalled
    # mrview reads do_semtimedop, and `ipcs -s` shows waiters (ncount > 0).
    _render_par=${KUL_RENDER_PAR:-3}
    [ "$_render_par" -lt 1 ] 2>/dev/null && _render_par=1
    mrview_threads=$(( ncpu / _render_par ))
    [ $mrview_threads -lt 1 ] && mrview_threads=1
    echo "Render parallelism: $_render_par concurrent mrview, $mrview_threads Mesa thread(s) each"

    # --- Portable headless-mrview environment (Linux Mint / Ubuntu, GPU or not) ---
    # Every mrview capture below is prefixed with $_mrview_env. This:
    #   * strips the MATLAB MCR's Qt from LD_LIBRARY_PATH. That Qt ships no platform
    #     plugins, so if it shadows the system Qt you get the classic abort:
    #       "Could not find the Qt platform plugin xcb/offscreen in ''".
    #     All MCR lib dirs live under .../glnxa64, so removing that single token is a
    #     precise filter that leaves CUDA (and everything else) intact. It is applied
    #     ONLY to the mrview call, so any MCR-based step elsewhere keeps its runtime.
    #     /snap/ is filtered for the same reason -- a snap's bundled libs are built
    #     against a different glibc than the host's.
    #   * forces Mesa llvmpipe (LIBGL_ALWAYS_SOFTWARE=1) so software GL is used
    #     deterministically even on machines that DO have a GPU but no display.
    #   * drops any leaked QT_QPA_PLATFORM=offscreen / QT_PLUGIN_PATH from the shell,
    #     which would otherwise override the xcb plugin xvfb-run provides.
    #   * drops the GTK/GLib module search vars. This one is not theoretical: when
    #     the pipeline is launched from a terminal inside snap-packaged VS Code,
    #     GTK_PATH points at /snap/code/<rev>/usr/lib/x86_64-linux-gnu/gtk-3.0.
    #     Qt loads libcanberra-gtk-module.so from there, and that module's RPATH
    #     drags /snap/core20/current/lib/x86_64-linux-gnu/libpthread.so.0 in
    #     against the host glibc, so mrview dies before it renders anything:
    #       symbol lookup error: .../libpthread.so.0: undefined symbol:
    #       __libc_pthread_init, version GLIBC_PRIVATE
    #     GTK_PATH alone is enough to trigger it and enough to fix it; the rest
    #     are unset alongside because they leak from the same snap wrapper.
    #     This is why the step could work on one machine and fail on another --
    #     it depends on how the terminal was launched, not on the host.
    # Override the binary with MRVIEW_BIN=/path/to/mrview if PATH is ambiguous.
    _mrview_ld=$(printf '%s' "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | grep -v 'glnxa64' | grep -v '^/snap/' | paste -sd:)
    _mrview_bin="${MRVIEW_BIN:-$(command -v mrview)}"
    _mrview_env="env -u QT_QPA_PLATFORM -u QT_PLUGIN_PATH \
        -u GTK_PATH -u GTK_MODULES -u GTK_EXE_PREFIX -u GTK_IM_MODULE_FILE \
        -u GIO_MODULE_DIR -u GSETTINGS_SCHEMA_DIR -u LOCPATH \
        -u GDK_PIXBUF_MODULE_FILE -u GDK_PIXBUF_MODULEDIR \
        LD_LIBRARY_PATH=$_mrview_ld LIBGL_ALWAYS_SOFTWARE=1"
    if [ -z "$_mrview_bin" ]; then
        echo "ERROR: mrview not found on PATH. Set MRVIEW_BIN=/path/to/mrview or fix PATH."
        _mrview_bin="mrview"   # fall through; the failure will be explicit
    fi
    if ! command -v xvfb-run >/dev/null 2>&1; then
        echo "WARNING: xvfb-run not found — mrview screenshots will fail on this host."
        echo "         Install with: sudo apt install -y xvfb libgl1-mesa-dri"
    fi

    # Base X display number for this run. xvfb-run -a searches upward from here
    # for a free server, so a stale /tmp/.X<n>-lock left by a killed run (or an
    # ssh -X session, which allocates from :10) no longer wedges the render --
    # that used to fail with zero PNGs and no error. The base is randomised so
    # two concurrent pipeline runs don't both start probing the same number.
    # Defined here because the preflight render below already needs it.
    _xvfb_base=$(( 50 + RANDOM % 40 ))

    # Preflight A: warn about mrview processes wedged from an earlier run.
    #
    # A hung mrview holds a SysV semaphore (seen with `ipcs -s`, wchan
    # do_semtimedop) that it never releases, and EVERY later mrview on the host
    # then blocks on it forever. One wedged run therefore poisons every
    # subsequent one — renders that produce no PNGs, no error, and eventually
    # just hit the timeout. Detected and reported rather than killed
    # automatically: on a shared machine those processes may belong to someone
    # else's live run, and that call is not this script's to make.
    _stale_mrview=$(pgrep -x mrview 2>/dev/null | wc -l)
    if [ "$_stale_mrview" -gt 0 ]; then
        echo ""
        echo "WARNING: ${_stale_mrview} mrview process(es) are already running on this host."
        echo "         If they are wedged from an earlier run they hold a semaphore that will"
        echo "         block every render below indefinitely. Check their age with:"
        echo "           ps -o pid,etime,cmd -p \$(pgrep -dx, mrview) | cut -c1-100"
        echo "         If they are stale, clear them (and any orphaned Xvfb) with:"
        echo "           pkill -x mrview; pkill -f 'Xvfb :'"
        echo "           ipcs -s | awk '/^0x/ {print \$2}' | xargs -r -n1 ipcrm -s"
        echo ""
    fi

    # Preflight B: render one real slice before committing to hundreds.
    #
    # Deliberately a capture and not `mrview --version`: --version never opens a
    # window or touches OpenGL, so it succeeds on a host where every actual
    # render hangs. This catches both failure modes we have actually seen — the
    # snap/GTK library clash (mrview dies immediately) and the wedged-semaphore
    # cascade (mrview starts and blocks forever) — in a few seconds, instead of
    # after a 10-minute timeout on each of ~130 series.
    _pf_dir=$(mktemp -d "${TMPDIR:-/tmp}/kul_mrview_preflight_XXXXXX")
    _pf_tmp=$(mktemp -d "${TMPDIR:-/tmp}/kul_mrview_pftmp_XXXXXX")
    _mrview_check=$(eval "$_mrview_env TMPDIR=$_pf_tmp timeout -k 10 120 \
        xvfb-run -a -n $_xvfb_base -e /dev/stderr \
        --server-args=\"-screen 0 ${mrview_resolution}x${mrview_resolution}x24\" \
        $_mrview_bin -size $mrview_resolution,$mrview_resolution \
        -load $underlay -mode 1 -plane 2 -noannotations \
        -capture.folder $_pf_dir -capture.prefix pf -voxel 0,0,0 -capture.grab \
        -force -exit" 2>&1)
    _pf_rc=$?
    _pf_n=$(ls "$_pf_dir"/*.png 2>/dev/null | wc -l)
    rm -rf "$_pf_dir" "$_pf_tmp"
    if [ $_pf_rc -ne 0 ] || [ "$_pf_n" -eq 0 ]; then
        echo ""
        echo "ERROR: mrview cannot render. Screenshots (and therefore all PACS/Karawun"
        echo "       DICOM output) would fail. Test render exit=${_pf_rc}, PNGs=${_pf_n}."
        printf '         %s\n' "$(printf '%s' "$_mrview_check" | grep -viE 'keysym|xkbcomp|^>|Errors from' | head -5)"
        echo ""
        if [ $_pf_rc -eq 124 ] || [ $_pf_rc -eq 137 ]; then
            echo "       It HUNG rather than failed. That is almost always a wedged mrview"
            echo "       from an earlier run holding a semaphore every later one waits on:"
            echo "         ps -o pid,etime,cmd -p \$(pgrep -dx, mrview) | cut -c1-100"
            echo "         pkill -x mrview; pkill -f 'Xvfb :'"
            echo "         ipcs -s | awk '/^0x/ {print \$2}' | xargs -r -n1 ipcrm -s"
        else
            echo "       If the output mentions /snap/, the environment is leaking a snap's"
            echo "       libraries into mrview. This block already unsets the usual culprits"
            echo "       (GTK_PATH etc.); check for others with:  env | grep -i snap"
            echo "       Running from a plain terminal rather than a snap-packaged editor's"
            echo "       integrated terminal is the quickest workaround."
        fi
        echo ""
        exit 1
    fi
    echo "mrview preflight OK (rendered a test slice)"

    mkdir -p $resultsdir_png
    mkdir -p $resultsdir_dcm

    # Derive underlay suffix for SPM figure directories (matches Tracto naming)
    _ulsuffix="${resultsdir_png##*_figures_}"
    spm_resultsdir_png="$globalresultsdir/SPM_figures_${_ulsuffix}"
    spm_resultsdir_dcm="$globalresultsdir/PACS/fMRI_${_ulsuffix}"
    mkdir -p "$spm_resultsdir_png"
    mkdir -p "$spm_resultsdir_dcm"

    # --- Failure collection + per-series logs -------------------------------
    # Render/convert failures used to be swallowed ("|| echo WARNING") and the
    # run always exited 0. Children are backgrounded, so failures are collected
    # in a file rather than a bash array (a child's array write never reaches
    # the parent). Single short appends to the same file are atomic enough.
    _pacs_logdir="$cwd/KUL_LOG/sub-${participant}_PACS"
    mkdir -p "$_pacs_logdir"
    _pacs_fail_file=$(mktemp "${TMPDIR:-/tmp}/kul_pacs_fail_XXXXXX")
    : > "$_pacs_fail_file"
    _pacs_note_fail() {
        # kind, label, detail
        printf '%-8s %-48s %s\n' "$1" "$2" "$3" >> "$_pacs_fail_file"
    }

    # Run one mrview capture and verify it produced what it should.
    #   $1 slot   $2 png_dir   $3 expected PNG count   $4 label   $5 mrview args
    # Returns 0 only if mrview exited cleanly AND the full slice count landed.
    # A short render is retried once from scratch; previously any single PNG
    # counted as "done", so a run killed by the timeout left a truncated series
    # that every later run skipped.
    _mrview_capture() {
        local slot="$1" png_dir="$2" expected="$3" label="$4" args="$5"
        local attempt rc got logf mtmp
        logf="$_pacs_logdir/${label}.log"
        for attempt in 1 2; do
            mkdir -p "$png_dir"
            # Private TMPDIR per invocation. Qt builds a QSystemSemaphore whose
            # SysV key is ftok()'d from a backing file it creates in $TMPDIR
            # (.../qipc_systemsem_<hash>). Every mrview on the host hashes to
            # the same name, so they all share ONE semaphore and serialise on
            # it -- and if one is SIGKILLed while holding it (which the timeout
            # below does), every later mrview blocks in semtimedop forever, with
            # zero CPU, including ones started by hand from another terminal.
            # That is the hang that wedged this host for 8.5 h. A per-process
            # TMPDIR gives each render its own key, so they cannot collide and
            # nothing is left behind for the next run to trip over.
            mtmp=$(mktemp -d "${TMPDIR:-/tmp}/kul_mrview_XXXXXX")
            eval "$_mrview_env TMPDIR=$mtmp LP_NUM_THREADS=$mrview_threads timeout -k 30 600 \
                xvfb-run -a -n $(( _xvfb_base + slot )) -e /dev/stderr \
                --server-args=\"-screen 0 ${mrview_resolution}x${mrview_resolution}x24\" \
                $_mrview_bin $args" >>"$logf" 2>&1
            rc=$?
            rm -rf "$mtmp"
            got=$(ls "$png_dir"/*.png 2>/dev/null | wc -l)
            if [ $rc -eq 0 ] && [ "$got" -eq "$expected" ]; then
                return 0
            fi
            echo "WARNING: ${label} rendered ${got}/${expected} PNGs (mrview exit ${rc}), attempt ${attempt}/2"
            if [ $attempt -lt 2 ]; then
                rm -rf "$png_dir"
            fi
        done
        _pacs_note_fail "render" "$label" "${got}/${expected} PNGs, mrview exit ${rc} — see ${logf}"
        return 1
    }

    # True when a screenshot folder already holds a COMPLETE set of slices.
    _png_dir_complete() {
        local d="$1" expected="$2" got
        got=$(ls "$d"/*.png 2>/dev/null | wc -l)
        [ "$got" -eq "$expected" ] && [ "$expected" -gt 0 ]
    }

    # Number of slices mrview will capture for an orientation, from the underlay.
    # Computed once here (it is the same for every map/bundle) rather than
    # re-running mrinfo per orientation per render as before.
    _ul_slices_SAG=${_dims[0]}
    _ul_slices_COR=${_dims[1]}
    _ul_slices_TRA=${_dims[2]}
    _expected_slices() {
        case "$1" in
            TRA) printf '%s' "$_ul_slices_TRA" ;;
            SAG) printf '%s' "$_ul_slices_SAG" ;;
            COR) printf '%s' "$_ul_slices_COR" ;;
            *)   printf '0' ;;
        esac
    }
    # mrview -plane index for an orientation (0=sagittal, 1=coronal, 2=axial).
    _plane_of() {
        case "$1" in TRA) printf '2' ;; SAG) printf '0' ;; COR) printf '1' ;; *) printf '2' ;; esac
    }

    # Robust display range for an overlay: 2nd-98th percentile over the
    # non-zero, finite voxels (optionally restricted to voxels at or above a
    # threshold, so a thresholded map's colours span what is actually shown).
    #
    # Every overlay gets one. mrview's -overlay.threshold_min only decides which
    # voxels are drawn, not how values map to colours: without an explicit
    # -overlay.intensity, mrview windows the colourmap over the volume's full
    # range, so a map whose values run into the thousands (rCBV, ALFF) renders
    # as a saturated all-white blob, while a narrow-range fMRI t-map happens to
    # look fine. That is why fusing perfusion maps on an anatomical never
    # worked. A single hot vessel voxel would also compress the whole brain
    # into the bottom of the colourmap, hence percentiles rather than min/max.
    # max / p99 over non-zero finite voxels: how outlier-driven the maximum is.
    # A statistical map sits near 1-3; a physiological map with vessel voxels
    # (rCBV, ALFF) runs far higher. Used to spot maps for which the max/3 auto
    # threshold is meaningless -- see _render_one_spm.
    _map_tail_ratio() {
        "$_nii2dcm_py" -c "
import sys, numpy as np, SimpleITK as sitk
a = sitk.GetArrayFromImage(sitk.ReadImage(sys.argv[1])).astype('float64')
a = a[np.isfinite(a)]; a = a[a != 0]
if a.size == 0: sys.exit(1)
mx = float(a.max()); p99 = float(np.percentile(a, 99))
print(f'{(mx/p99) if p99 > 0 else 0:.3f}')
" "$1" 2>/dev/null
    }

    # Reference tissue for anchoring the display window of a continuous overlay.
    # Grey matter, i.e. the cortex, from the anatomical segmentation.
    _gm_ref="$globalresultsdir/Anat/T1w_GM.nii.gz"
    # Upper end of the window, as a multiple of the map's median in that tissue.
    # 4x puts normal cortex at a quarter scale and saturates only the hottest
    # ~3-4% of the brain, which for rCBV is tumour and vessel.
    _win_ref_mult=4

    _map_intensity_range() {
        "$_nii2dcm_py" -c "
import sys, numpy as np, SimpleITK as sitk
a = sitk.GetArrayFromImage(sitk.ReadImage(sys.argv[1])).astype('float64')
thr = float(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] else None
ref = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] else None
mult = float(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else 4.0

# Upper bound, preferred: a multiple of the map's median inside a reference
# tissue. A percentile of the map itself drifts with how much tumour and vessel
# happen to sit in the field of view, so the same physiology gets a different
# colour in different patients. Anchoring to cortex makes the scale mean
# 'x cortical value' and stay comparable between studies and scanners.
hi_ref = None
if ref:
    try:
        m = sitk.GetArrayFromImage(sitk.ReadImage(ref)).astype('float64')
        if m.shape == a.shape:
            sel = (m > 0.5) & np.isfinite(a) & (a != 0)
            if sel.sum() > 1000:
                med = float(np.median(a[sel]))
                if np.isfinite(med) and med > 0:
                    hi_ref = mult * med
    except Exception:
        pass

a = a[np.isfinite(a)]
a = a[a >= thr] if thr is not None else a[a != 0]
if a.size == 0:
    sys.exit(1)
lo = float(np.percentile(a, 2))
hi = hi_ref if hi_ref is not None else float(np.percentile(a, 98))
if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
    lo, hi = float(a.min()), float(a.max())
if hi <= lo:
    hi = lo + 1e-6

# Round the ends so the colourbar reads cleanly: whole numbers once the values
# are >= 1, otherwise at most 3 decimals (0.01, 0.005). Reverted if rounding
# would collapse or invert the range, which it can for very small-valued maps.
def nice(v):
    return float(round(v)) if abs(v) >= 1 else round(v, 3)
r_lo, r_hi = nice(lo), nice(hi)
if r_hi > r_lo:
    lo, hi = r_lo, r_hi

print(f'{lo:g},{hi:g}')
" "$1" "${2:-}" "${3:-}" "$_win_ref_mult" 2>/dev/null
    }

    # Wrapper around KUL_nii2dcm.py: uses the resolved interpreter, logs, and
    # records failures instead of discarding the exit status.
    #   $1 label (also the log name), rest = arguments to KUL_nii2dcm.py
    _run_nii2dcm() {
        local label="$1"; shift
        local logf="$_pacs_logdir/${label}.dcm.log"
        local rc
        # Status captured directly, not after an `if`: a false `if` with no
        # `else` yields 0, so reading $? there reported every failure as "exit 0".
        "$_nii2dcm_py" "$_nii2dcm_script" "$@" >>"$logf" 2>&1
        rc=$?
        [ $rc -eq 0 ] && return 0
        echo "WARNING: KUL_nii2dcm.py failed for ${label} (exit ${rc}) — see ${logf}"
        _pacs_note_fail "convert" "$label" "KUL_nii2dcm.py exit ${rc} — see ${logf}"
        return $rc
    }

    # Render one SPM/Melodic map (all orientations).
    # Threshold is computed per-map as max/3 (same as the report logic).
    _render_one_spm() {
        local spmfile="$1" spmname="$2" slot="${3:-0}" render_mode="${4:-thresholded}"
        local ori
        IFS=',' read -ra ori <<< "$orientations"

        # Resolve the threshold. "none" (from a <name>.thresh sidecar) means show
        # the whole range -- right for physiological maps like rCBV/ALFF/ReHo,
        # where thresholding would hide most of what the reader wants to see.
        local _thresh="" _max_T
        if [ "$render_mode" = "continuous" ]; then
            _thresh=""
        elif [ -n "${spm_thresh_map[$spmname]+x}" ]; then
            _thresh=${spm_thresh_map[$spmname]}
            if [ "$_thresh" = "none" ]; then
                _thresh=""
                echo "Map ${spmname}: no threshold (full range)"
            else
                echo "Map ${spmname}: using per-map threshold=${_thresh}"
            fi
        elif [ -n "$spm_thresh_override" ]; then
            _thresh=$spm_thresh_override
            echo "Map ${spmname}: using fixed threshold=${_thresh}"
        else
            # Automatic threshold is max/3, which assumes the maximum is a
            # meaningful peak -- true for a statistical map, false for a
            # physiological one whose maximum is a vessel voxel. On real data
            # rCBV had max=12111 against a p99 of 1543, so max/3 kept 704 of
            # 2,083,162 voxels (0.03%) and the overlay was effectively blank.
            # When the maximum is that far out in the tail, threshold the map at
            # all is the wrong question: these are brain-masked physiological
            # maps meant to be read brain-wide, so render them continuous
            # (full range, robust window, colourbar) instead.
            # An explicit threshold (-T, sidecar, or the fixed clinical cutoffs
            # for nrCBV/nrCBF) always wins -- this only affects the auto case.
            local _ratio
            _ratio=$(_map_tail_ratio "$spmfile")
            if [ -n "$_ratio" ] && awk "BEGIN {exit !($_ratio > 3)}"; then
                echo "Map ${spmname}: max/p99 = ${_ratio} (outlier-driven maximum)"
                echo "  -> rendering brain-wide with a colourbar instead of an auto threshold"
                _thresh=""
            else
                _max_T=$(mrstats -output max "$spmfile")
                _thresh=$(awk "BEGIN {print $_max_T/3}")
                echo "Map ${spmname}: max=${_max_T}, auto threshold=${_thresh}"
            fi
        fi

        # Colour windowing, ALWAYS. mrview's threshold decides which voxels are
        # drawn; the intensity range decides how their values map to colours.
        # Without this a map whose values run into the thousands (rCBV, ALFF)
        # renders as a saturated all-white blob, while a narrow-range fMRI t-map
        # happens to look right -- the reason perfusion overlays never worked.
        # For a thresholded map the range is taken over the suprathreshold
        # voxels, so the colourmap spans exactly what is on screen.
        local _overlay_opts _colourbar=0 _range
        if [ -z "$_thresh" ]; then
            # Continuous (physiological) map: anchor the top of the window to
            # cortical values so the colours mean the same thing across studies.
            _range=$(_map_intensity_range "$spmfile" "" "$_gm_ref")
        else
            # Thresholded (statistical) map: the cortex median of a t-map is not
            # a meaningful anchor, so window over the suprathreshold voxels.
            _range=$(_map_intensity_range "$spmfile" "$_thresh")
        fi
        if [ -z "$_range" ]; then
            echo "WARNING: no usable intensity range for ${spmname} (empty, all-NaN, or nothing above threshold)"
            _pacs_note_fail "render" "$spmname" "no usable intensity range"
            return 1
        fi
        echo "Map ${spmname}: overlay intensity range ${_range}"
        _overlay_opts="-overlay.intensity $_range"
        if [ -n "$_thresh" ]; then
            _overlay_opts="$_overlay_opts -overlay.threshold_min $_thresh"
        else
            # No threshold: show the full map and give the reader a scale.
            _overlay_opts="$_overlay_opts -overlay.no_threshold_min -overlay.no_threshold_max"
            _colourbar=1
        fi

        # Optionally compute a 1-voxel edge mask for the dark-blue outline
        # (only meaningful for a thresholded map -- there is no edge without one)
        local _edge_overlay="" _tmp_mask="" _tmp_eroded="" _tmp_edge=""
        if [ $spm_edge -eq 1 ] && [ -n "$_thresh" ]; then
            _tmp_mask=$(mktemp /tmp/spm_mask_XXXXXX.nii.gz)
            _tmp_eroded=$(mktemp /tmp/spm_eroded_XXXXXX.nii.gz)
            _tmp_edge=$(mktemp /tmp/spm_edge_XXXXXX.nii.gz)
            mrcalc "$spmfile" $_thresh -ge "$_tmp_mask" -force -quiet
            maskfilter "$_tmp_mask" erode "$_tmp_eroded" -npass 1 -force -quiet
            mrcalc "$_tmp_mask" "$_tmp_eroded" -sub "$_tmp_edge" -force -quiet
            _edge_overlay="-overlay.load $_tmp_edge -overlay.opacity 1.0 -overlay.colour 0,0,0.8 -overlay.threshold_min 0.5"
        fi

        local _px_tra=$px_tra _px_sag=$px_sag _px_cor=$px_cor _sample_png=""
        for orient in "${ori[@]}"; do
            local plane underlay_slices
            underlay_slices=$(_expected_slices "$orient")
            plane=$(_plane_of "$orient")

            local png_dir="$spm_resultsdir_png/${spmname}_${orient}"
            if _png_dir_complete "$png_dir" "$underlay_slices"; then
                echo "Skipping screenshots for ${spmname}_${orient} (complete: ${underlay_slices} slices)"
            else
                mkdir -p "$png_dir"
                local voxel_index="" i=0
                while [ $i -lt $underlay_slices ]; do
                    if [[ "$orient" == "TRA" ]]; then
                        voxel_index="$voxel_index -voxel 0,0,$i -capture.grab"
                    elif [[ "$orient" == "SAG" ]]; then
                        voxel_index="$voxel_index -voxel $i,0,0 -capture.grab"
                    else
                        voxel_index="$voxel_index -voxel 0,$i,0 -capture.grab"
                    fi
                    let "i+=1"
                done
                echo "Making ${spmname}_${orient} SPM on $(basename $underlay)"
                _mrview_capture "$slot" "$png_dir" "$underlay_slices" "${spmname}_${orient}" \
                    "-size $mrview_resolution,$mrview_resolution \
                    -load $underlay -mode 1 -plane $plane \
                    -overlay.load $spmfile -overlay.opacity $spm_opacity -overlay.colourmap 1 \
                        $_overlay_opts \
                    $_edge_overlay \
                    -noannotations -orientlabel 0 -voxelinfo 0 -colourbar $_colourbar \
                    -capture.folder $png_dir -capture.prefix ${spmname}_${orient} \
                    $voxel_index -force -exit"
            fi
            [ -z "$_sample_png" ] && _sample_png=$(ls "$png_dir"/*.png 2>/dev/null | head -1)
        done

        [ -n "$_tmp_mask" ] && rm -f "$_tmp_mask" "$_tmp_eroded" "$_tmp_edge"

        # Recompute PixelSpacing from first screenshot (same as tract logic)
        if [ -n "$_sample_png" ]; then
            local _png_max
            _png_max=$("$_nii2dcm_py" -c "from PIL import Image; w,h=Image.open('$_sample_png').size; print(max(w,h))")
            _px_tra=$("$_nii2dcm_py" -c "print(max(${_dims[0]}*${_vox[0]},${_dims[1]}*${_vox[1]})/$_png_max)")
            _px_sag=$("$_nii2dcm_py" -c "print(max(${_dims[1]}*${_vox[1]},${_dims[2]}*${_vox[2]})/$_png_max)")
            _px_cor=$("$_nii2dcm_py" -c "print(max(${_dims[0]}*${_vox[0]},${_dims[2]}*${_vox[2]})/$_png_max)")
        fi

        # DICOM conversion (only when explicitly requested via -R)
        # Skip if this map was not selected by the user (_dcm_spm_set empty = all selected)
        if [ $make_dcm -eq 1 ] && [ -n "$donor_dcm" ] && \
           { [ ${#_dcm_spm_set[@]} -eq 0 ] || [ -n "${_dcm_spm_set[$spmname]+x}" ]; }; then
            local dcm_label="${spmname}_on_${_ulsuffix}"
            for orient in "${ori[@]}"; do
                local ps snum _nopt
                case $orient in TRA) ps=$_px_tra ;; SAG) ps=$_px_sag ;; COR) ps=$_px_cor ;; esac
                local png_dir="$spm_resultsdir_png/${spmname}_${orient}"
                local dcmdir="$spm_resultsdir_dcm/${dcm_label}_${orient}"
                # Only convert a complete screenshot set — a short render would
                # otherwise be silently wrapped into a truncated series.
                if ! _png_dir_complete "$png_dir" "$(_expected_slices "$orient")"; then
                    echo "Skipping DICOMs for ${dcm_label}_${orient} (screenshots incomplete)"
                    continue
                fi
                snum=$(_pacs_series_number "map:${_ulsuffix}:$spmname" "$orient")
                _nopt=""; [ -n "$snum" ] && _nopt="-n $snum"
                if [ -n "$(ls "$dcmdir"/*.dcm 2>/dev/null | head -1)" ]; then
                    echo "Skipping DICOMs for ${dcm_label}_${orient} (already exist)"
                else
                    mkdir -p "$dcmdir"
                    if [[ "$orient" == "SAG" ]]; then
                        echo "Making dicoms in $dcmdir (donor-match mode, SeriesNumber=${snum:-<none>})"
                        _run_nii2dcm "${dcm_label}_${orient}" -s "${dcm_label}_${orient}" $_nopt \
                            -u "$underlay" -o "$orient" -M \
                            "$png_dir" "$donor_dcm" "$dcmdir"
                    else
                        echo "Making dicoms in $dcmdir (PixelSpacing=${ps}mm, SeriesNumber=${snum:-<none>})"
                        _run_nii2dcm "${dcm_label}_${orient}" -s "${dcm_label}_${orient}" -p $ps $_nopt \
                            -u "$underlay" -o "$orient" \
                            "$png_dir" "$donor_dcm" "$dcmdir"
                    fi
                fi
            done
        fi
    }

    # Render one bundle (all orientations sequentially) on its own X display
    _render_one_bundle() {
        local slot="$1" tractname="$2" mrview_tck="$3"
        local ori
        IFS=',' read -ra ori <<< "$orientations"

        for orient in "${ori[@]}"; do
            local underlay_slices plane png_dir
            underlay_slices=$(_expected_slices "$orient")
            plane=$(_plane_of "$orient")
            png_dir="$resultsdir_png/${tractname}_${orient}"
            if _png_dir_complete "$png_dir" "$underlay_slices"; then
                echo "Skipping screenshots for ${tractname}_${orient} (complete: ${underlay_slices} slices)"
            else
                echo "Making ${tractname}_${orient} on $(basename $underlay)"
                mkdir -p "$png_dir"
                local voxel_index="-capture.folder $png_dir -capture.prefix ${tractname}_${orient}"
                local i=0
                while [ $i -lt $underlay_slices ]; do
                    if [[ "$orient" == "TRA" ]]; then
                        voxel_index="$voxel_index -voxel 0,0,$i -capture.grab"
                    elif [[ "$orient" == "SAG" ]]; then
                        voxel_index="$voxel_index -voxel $i,0,0 -capture.grab"
                    else
                        voxel_index="$voxel_index -voxel 0,$i,0 -capture.grab"
                    fi
                    let "i+=1"
                done
                _mrview_capture "$slot" "$png_dir" "$underlay_slices" "${tractname}_${orient}" \
                    "-size $mrview_resolution,$mrview_resolution \
                    -load $underlay -mode 1 -plane $plane \
                    -tractography.lighting 1 -tractography.slab 1.5 -tractography.thickness 0.3 \
                    -noannotations -orientlabel 0 -voxelinfo 0 -colourbar 0 \
                    $mrview_tck $voxel_index -force -exit"
            fi
        done

        # Recompute PixelSpacing from first available screenshot
        local _px_tra=$px_tra _px_sag=$px_sag _px_cor=$px_cor _sample_png=""
        for orient in "${ori[@]}"; do
            _sample_png=$(ls "$resultsdir_png/${tractname}_${orient}"/*.png 2>/dev/null | head -1)
            [ -n "$_sample_png" ] && break
        done
        if [ -n "$_sample_png" ]; then
            local _png_max
            _png_max=$("$_nii2dcm_py" -c "from PIL import Image; w,h=Image.open('$_sample_png').size; print(max(w,h))")
            _px_tra=$("$_nii2dcm_py" -c "print(max(${_dims[0]}*${_vox[0]},${_dims[1]}*${_vox[1]})/$_png_max)")
            _px_sag=$("$_nii2dcm_py" -c "print(max(${_dims[1]}*${_vox[1]},${_dims[2]}*${_vox[2]})/$_png_max)")
            _px_cor=$("$_nii2dcm_py" -c "print(max(${_dims[0]}*${_vox[0]},${_dims[2]}*${_vox[2]})/$_png_max)")
            echo "Actual PNG max dim: ${_png_max}px → PixelSpacing TRA=${_px_tra}mm SAG=${_px_sag}mm COR=${_px_cor}mm"
        fi

        # DICOM conversion (only when explicitly requested via -R)
        if [ $make_dcm -eq 1 ] && [ -n "$donor_dcm" ]; then
            for orient in "${ori[@]}"; do
                local ps snum _nopt
                case $orient in TRA) ps=$_px_tra ;; SAG) ps=$_px_sag ;; COR) ps=$_px_cor ;; esac
                local png_dir="$resultsdir_png/${tractname}_${orient}"
                local dcmdir="$resultsdir_dcm/${tractname}_${orient}"
                # Only convert a complete screenshot set — a short render would
                # otherwise be silently wrapped into a truncated series.
                if ! _png_dir_complete "$png_dir" "$(_expected_slices "$orient")"; then
                    echo "Skipping DICOMs for ${tractname}_${orient} (screenshots incomplete)"
                    continue
                fi
                snum=$(_pacs_series_number "tract:${_ulsuffix}:$tractname" "$orient")
                _nopt=""; [ -n "$snum" ] && _nopt="-n $snum"
                if [ -n "$(ls "$dcmdir"/*.dcm 2>/dev/null | head -1)" ]; then
                    echo "Skipping DICOMs for ${tractname}_${orient} (already exist)"
                else
                    mkdir -p "$dcmdir"
                    if [[ "$orient" == "SAG" ]]; then
                        echo "Making dicoms in $dcmdir (donor-match mode, SeriesNumber=${snum:-<none>})"
                        _run_nii2dcm "FT_${tractname}_${orient}" \
                            -s "FT_${tractname}_${orient}_${_ulsuffix}" $_nopt \
                            -u "$underlay" -o "$orient" -M \
                            "$png_dir" "$donor_dcm" "$dcmdir"
                    else
                        echo "Making dicoms in $dcmdir (PixelSpacing=${ps}mm, SeriesNumber=${snum:-<none>})"
                        _run_nii2dcm "FT_${tractname}_${orient}" \
                            -s "FT_${tractname}_${orient}_${_ulsuffix}" -p $ps $_nopt \
                            -u "$underlay" -o "$orient" \
                            "$png_dir" "$donor_dcm" "$dcmdir"
                    fi
                fi
            done
        fi
    }

    # Launch current batch of bundles in parallel, then wait.
    # SeriesNumbers are claimed HERE, in the parent, before anything is
    # backgrounded — two children racing to claim an index would otherwise be
    # able to land on the same number.
    _flush_bundle_batch() {
        for slot in "${!_btractnames[@]}"; do
            _pacs_register "tract:${_ulsuffix}:${_btractnames[$slot]}"
        done
        for slot in "${!_btractnames[@]}"; do
            _render_one_bundle "$slot" "${_btractnames[$slot]}" "${_btcks[$slot]}" &
        done
        wait
        _btractnames=()
        _btcks=()
    }

    # Same batching for the map renders (fMRI/Melodic/Clinical/drop-folder).
    # These were serial, and could not simply be backgrounded before, because
    # every instance rendered on the same hardcoded display :20; each now gets
    # its own slot off $_xvfb_base.
    _bspmfiles=(); _bspmnames=(); _bspmmodes=()
    _flush_spm_batch() {
        [ ${#_bspmnames[@]} -eq 0 ] && return 0
        local slot
        for slot in "${!_bspmnames[@]}"; do
            _pacs_register "map:${_ulsuffix}:${_bspmnames[$slot]}"
        done
        for slot in "${!_bspmnames[@]}"; do
            _render_one_spm "${_bspmfiles[$slot]}" "${_bspmnames[$slot]}" "$slot" "${_bspmmodes[$slot]}" &
        done
        wait
        _bspmfiles=(); _bspmnames=(); _bspmmodes=()
    }
    # Queue one map; flushes automatically once a full batch has accumulated.
    _queue_spm() {
        _bspmfiles+=("$1"); _bspmnames+=("$2"); _bspmmodes+=("${3:-thresholded}")
        [ ${#_bspmnames[@]} -ge $_render_par ] && _flush_spm_batch
        return 0
    }

    # ── Manual drop folders ─────────────────────────────────────────────────
    # Auto-discovery below globs SPM/, Melodic/, Perfusion/ and Lesion/, which
    # only works for maps this pipeline produced under the names it expects.
    # These folders let you export anything by copying the file in — no
    # renaming, no processing. Two categories, by what you want out:
    #
    #   overlays/             colour overlay fused on the anatomical, for
    #                         review. Colour windowing always adapts to the
    #                         map's own range, so a wide-range map (rCBV) and a
    #                         narrow-range one (fMRI t-map) both render legibly.
    #   series_quantitative/  a real measurable DICOM series, no rendering.
    #                         Pixel values carry RescaleSlope/Intercept so an
    #                         ROI drawn on PACS reads true units — the same
    #                         thing the scanner's own ADC maps give you.
    #
    # The same file may be copied (or symlinked) into both to get a fusion
    # overlay and a measurable series from one map.
    #
    # Discovery is CONDITIONAL: if both are empty or absent, the automatic
    # SPM/Melodic/Perfusion/Lesion discovery runs exactly as it always has. As
    # soon as either holds a file, only the drop folders are used.
    # _dropdir/_drop_thr/_drop_quant were set by KUL_make_pacs_dropdirs, called
    # near the top of this block (and at the end of a normal pipeline run).

    # The two folders gate INDEPENDENTLY. Counting them together meant a single
    # file dropped in series_quantitative/ switched the overlay side to "manual"
    # as well, and since overlays/ was empty that silently produced no overlays
    # at all -- the automatic SPM/Melodic/Perfusion/Lesion discovery having been
    # turned off by a file that had nothing to do with it.
    _drop_count=$(find "$_drop_thr" \( -name '*.nii' -o -name '*.nii.gz' \) -type f 2>/dev/null | wc -l)
    _drop_quant_count=$(find "$_drop_quant" \( -name '*.nii' -o -name '*.nii.gz' \) -type f 2>/dev/null | wc -l)
    if [ "$_drop_count" -gt 0 ]; then
        _use_dropdir=1
        echo ""
        echo "Overlays: using $_drop_count file(s) from $_drop_thr"
        echo "          automatic SPM/Melodic/Perfusion/Lesion discovery skipped for this run."
        echo ""
    else
        _use_dropdir=0
    fi
    if [ "$_drop_quant_count" -gt 0 ]; then
        echo "Quantitative: using $_drop_quant_count file(s) from $_drop_quant"
    fi

    # ── fMRI (SPM & Melodic) → PACS ─────────────────────────────────────────
    # Deliberately BEFORE the tract renders below. This section opens with two
    # interactive prompts (per-map thresholds, then which maps to export as
    # DICOM); when it sat after the tract loop the operator had to wait out
    # every bundle render before they could answer, and could not leave the
    # terminal. Tract rendering needs no input, so putting fMRI first lets the
    # prompts be answered in the first seconds and the rest run unattended.
    # Nothing here depends on the tract renders: _render_one_spm, _ulsuffix and
    # spm_resultsdir_{png,dcm} are all defined further up.
    # If interactive and no global -T override, ask for a per-map threshold list.
    _all_spm_names=()
    if [ $_use_dropdir -eq 1 ]; then
        for _spm in "$_drop_thr"/*.nii.gz "$_drop_thr"/*.nii; do
            [ -f "$_spm" ] || continue
            _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
            _all_spm_names+=("$_spmname")
        done
    else
        for _spm in "$globalresultsdir/SPM/"*.nii.gz "$globalresultsdir/SPM/"*.nii \
                    "$globalresultsdir/Melodic/"*.nii.gz "$globalresultsdir/Melodic/"*.nii; do
            [ -f "$_spm" ] || continue
            _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
            _all_spm_names+=("$_spmname")
        done
    fi

    # A <name>.thresh sidecar pins a threshold without the interactive prompt,
    # so drop-folder runs can be fully unattended.
    for _idx in "${!_all_spm_names[@]}"; do
        _thrfile="$_drop_thr/${_all_spm_names[$_idx]}.thresh"
        if [ $_use_dropdir -eq 1 ] && [ -f "$_thrfile" ]; then
            _thrval=$(head -1 "$_thrfile" | tr -d '[:space:]')
            if [ -n "$_thrval" ]; then
                spm_thresh_map["${_all_spm_names[$_idx]}"]="$_thrval"
                echo "Threshold for ${_all_spm_names[$_idx]} pinned to $_thrval by $(basename "$_thrfile")"
            fi
        fi
    done

    if [ ${#_all_spm_names[@]} -gt 0 ] && [ -z "$spm_thresh_override" ] && [ -t 0 ]; then
        echo "fMRI/Melodic maps found for PACS conversion:"
        for _idx in "${!_all_spm_names[@]}"; do
            echo "  $((_idx+1)). ${_all_spm_names[$_idx]}"
        done
        read -r -p "Enter threshold values, space-separated, in the same order as above (Enter = auto max/3 for all): " -a _thresh_input
        if [ ${#_thresh_input[@]} -gt 0 ]; then
            if [ ${#_thresh_input[@]} -ne ${#_all_spm_names[@]} ]; then
                echo "Warning: got ${#_thresh_input[@]} value(s) for ${#_all_spm_names[@]} map(s) — ignoring, using auto threshold for all"
            else
                for _idx in "${!_all_spm_names[@]}"; do
                    spm_thresh_map["${_all_spm_names[$_idx]}"]="${_thresh_input[$_idx]}"
                done
            fi
        fi
    fi

    # When generating DICOMs (-R), let the user select which maps to export.
    # Empty _dcm_spm_set means "all maps".
    declare -A _dcm_spm_set=()
    if [ $make_dcm -eq 1 ] && [ ${#_all_spm_names[@]} -gt 0 ] && [ -t 0 ] && [ -z "$spm_thresh_override" ]; then
        echo ""
        echo "Select fMRI/Melodic maps to export as PACS DICOMs (underlay: ${_ulsuffix}):"
        for _idx in "${!_all_spm_names[@]}"; do
            echo "  $((_idx+1)). ${_all_spm_names[$_idx]}"
        done
        read -r -p "Enter numbers to include (space-separated), or Enter for all: " -a _sel_input
        if [ ${#_sel_input[@]} -gt 0 ]; then
            for _num in "${_sel_input[@]}"; do
                _sel_idx=$(( _num - 1 ))
                if [ $_sel_idx -ge 0 ] && [ $_sel_idx -lt ${#_all_spm_names[@]} ]; then
                    _dcm_spm_set["${_all_spm_names[$_sel_idx]}"]=1
                else
                    echo "Warning: ignoring out-of-range selection '$_num'"
                fi
            done
            echo "Selected for DICOM export: ${!_dcm_spm_set[*]}"
        else
            echo "All maps selected for DICOM export."
        fi
    fi

    if [ $_use_dropdir -eq 1 ]; then
        for _spm in "$_drop_thr"/*.nii.gz "$_drop_thr"/*.nii; do
            [ -f "$_spm" ] || continue
            _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
            _queue_spm "$_spm" "$_spmname" thresholded
        done
        _flush_spm_batch
    else
        for _spm in "$globalresultsdir/SPM/"*.nii.gz "$globalresultsdir/SPM/"*.nii; do
            [ -f "$_spm" ] || continue
            _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
            _queue_spm "$_spm" "$_spmname" thresholded
        done
        _flush_spm_batch
    fi


    _btractnames=()
    _btcks=()

    tract_set_i=0
    while true; do

        if [ $tract_set_i -lt $ntracts_paired ]; then
            # Paired mode: lateralized LT/RT bundles, step 2
            # tractname = first bundle name without Tract-csd_ prefix (e.g. CST_LT)
            tract_set=(${mrview_tracts[@]:$tract_set_i:2})
            tractname="${tract_set[0]#Tract-csd_}"; tractname="${tractname%_LT}"
            tract_i=$tract_set_i
            tract_set_i=$(($tract_set_i + 2))
        elif [ $tract_set_i -le $ntracts_total ]; then
            # Solo mode: commissural bundles, step 1
            tract_set=("${mrview_tracts[$tract_set_i]}")
            tractname="${mrview_tracts[$tract_set_i]#Tract-csd_}"
            tract_i=$tract_set_i
            tract_set_i=$(($tract_set_i + 1))
        else
            break
        fi

        mrview_tck=""
        tracts_found=0
        for tract in ${tract_set[@]}; do
            if [ -f $globalresultsdir/Tracto/${tract}.tck ]; then
                mrview_tck="$mrview_tck -tractography.load $globalresultsdir/Tracto/${tract}.tck \
                    -tractography.colour ${mrview_rgb[$tract_i]}"
                tracts_found=$(($tracts_found+1))
            fi
            tract_i=$(($tract_i+1))
        done

        echo "tract_i: $tract_i"
        echo "tracts_found: $tracts_found"

        if [ $tracts_found -gt 0 ]; then
            _btractnames+=("$tractname")
            _btcks+=("$mrview_tck")
            if [ ${#_btractnames[@]} -ge 3 ]; then
                _flush_bundle_batch
            fi
        else
            echo "No ${tractname} found"
        fi

    done
    [ ${#_btractnames[@]} -gt 0 ] && _flush_bundle_batch

    # ── Group renders: Projection / Associative / Commissural (3-way parallel) ─
    _btractnames=()
    _btcks=()

    _proj_tck=""
    for _gi in 0 1 16 17 20 21 22 23 24 25; do
        [ -f "$globalresultsdir/Tracto/${mrview_tracts[$_gi]}.tck" ] && \
            _proj_tck="$_proj_tck -tractography.load $globalresultsdir/Tracto/${mrview_tracts[$_gi]}.tck -tractography.colour ${mrview_rgb[$_gi]}"
    done
    if [ -n "$_proj_tck" ]; then
        _btractnames+=("Projection"); _btcks+=("$_proj_tck")
    else
        echo "No Projection tracts found"
    fi

    # SLF (26-33: SLF_all/I/II/III, LT/RT) is an association tract, not
    # commissural — it was previously grouped with Commissural below (a
    # leftover from before SLF was added to the tract array without
    # updating these ranges), so it never appeared in the Associative
    # DICOM output despite being one.
    _assoc_tck=""
    for _gi in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 18 19 26 27 28 29 30 31 32 33; do
        [ -f "$globalresultsdir/Tracto/${mrview_tracts[$_gi]}.tck" ] && \
            _assoc_tck="$_assoc_tck -tractography.load $globalresultsdir/Tracto/${mrview_tracts[$_gi]}.tck -tractography.colour ${mrview_rgb[$_gi]}"
    done
    if [ -n "$_assoc_tck" ]; then
        _btractnames+=("Associative"); _btcks+=("$_assoc_tck")
    else
        echo "No Associative tracts found"
    fi

    # Indices 34-42 are the true commissural tracts (Ant_Comm, Post_Comm,
    # and the 7 CC_*_Comm segments — see ntracts_paired/ntracts_total above).
    # The old range (26-34) both misclassified SLF as commissural and missed
    # 8 of these 9 genuine commissural tracts (only Ant_Comm/34 was caught).
    _comm_tck=""
    for _gi in 34 35 36 37 38 39 40 41 42; do
        [ -f "$globalresultsdir/Tracto/${mrview_tracts[$_gi]}.tck" ] && \
            _comm_tck="$_comm_tck -tractography.load $globalresultsdir/Tracto/${mrview_tracts[$_gi]}.tck -tractography.colour ${mrview_rgb[$_gi]}"
    done
    if [ -n "$_comm_tck" ]; then
        _btractnames+=("Commissural"); _btcks+=("$_comm_tck")
    else
        echo "No Commissural tracts found"
    fi

    [ ${#_btractnames[@]} -gt 0 ] && _flush_bundle_batch

    # ── Per-task scaled fMRI labels for Karawun/Brainlab ────────────────────
    # RESULTS/sub-.../SPM now holds exactly one (hardwired wc_p001unc_k50)
    # map per task. Binarize each at its resolved threshold and multiply by a
    # stable per-task integer (1..N, tasks sorted alphabetically) so each task
    # gets a distinct value — written as SEPARATE files, matching the existing
    # Karawun/Brainlab convention (see KUL_karawun2brainlab.sh: each tract/task
    # is its own scaled .nii, and `importTractography -l <dir>/*` already
    # takes multiple separate label files directly — no combining needed).
    # Only runs once Karawun prep has produced its own T1w.nii.gz (the grid
    # every Karawun label is regridded onto).
    if [ -f "Karawun/sub-${participant}/T1w.nii.gz" ]; then
        _spm_task_names=()
        for _spm in "$globalresultsdir/SPM/"*.nii.gz "$globalresultsdir/SPM/"*.nii; do
            [ -f "$_spm" ] || continue
            _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
            _spm_task_names+=("$_spmname")
        done
        if [ ${#_spm_task_names[@]} -gt 0 ]; then
            IFS=$'\n' _spm_task_names_sorted=($(sort <<<"${_spm_task_names[*]}")); unset IFS
            mkdir -p "Karawun/sub-${participant}/labels"

            # Palette indices reserved for fMRI activation labels.
            #
            # These used to be numbered 1..N, which put them straight on top of
            # the tract colours: task 1 got colour 1 = Arcuate Fasciculus, so on
            # a language case the activation and the language tract rendered
            # identically -- the one pair you most need to tell apart. Tasks 3,
            # 4 and 5 landed on CST_LT, CST_RT and Cingulum.
            #
            # Palette budget (see KUL_karawun_prepare.sh for the full map):
            #   1-41 known tracts, 42-49 auto tracts, 50 lesion, 51-63 fMRI.
            #
            # REQUIRES the extended-palette karawun fork. Every value here is
            # above 30, and stock karawun clamps anything above 30 to its last
            # entry -- so without the fork every activation renders in the same
            # colour as every other one. That is the deliberate trade for
            # keeping the whole low block available to tracts.
            _fmri_label_colors=(51 52 53 54 55 56 57 58 59 60 61 62 63)

            for _idx in "${!_spm_task_names_sorted[@]}"; do
                _spmname="${_spm_task_names_sorted[$_idx]}"
                if [ $_idx -lt ${#_fmri_label_colors[@]} ]; then
                    _label_int=${_fmri_label_colors[$_idx]}
                else
                    # more tasks than reserved colours; keep going past the end
                    # of the list rather than silently reusing one
                    _label_int=$((63 + _idx - ${#_fmri_label_colors[@]} + 1))
                    echo "WARNING: more fMRI maps than reserved label colours; using ${_label_int}"
                fi
                echo "Karawun fMRI label: task '${_spmname}' scaled to value ${_label_int}"
                _spmfile="$globalresultsdir/SPM/${_spmname}.nii"
                [ -f "$_spmfile" ] || _spmfile="$globalresultsdir/SPM/${_spmname}.nii.gz"
                [ -f "$_spmfile" ] || continue

                if [ -n "${spm_thresh_map[$_spmname]+x}" ]; then
                    _label_thresh=${spm_thresh_map[$_spmname]}
                elif [ -n "$spm_thresh_override" ]; then
                    _label_thresh=$spm_thresh_override
                else
                    _label_max_T=$(mrstats -output max "$_spmfile")
                    _label_thresh=$(awk "BEGIN {print $_label_max_T/3}")
                fi

                # SPM/*.nii(.gz) filenames already come out of
                # KUL_fmriproc_spm_new.sh prefixed with "afMRI_" — strip it
                # before re-adding once here, so this doesn't produce
                # afMRI_afMRI_... regardless of what the upstream naming does.
                _label_basename="${_spmname#afMRI_}"
                mrgrid "$_spmfile" regrid -template "Karawun/sub-${participant}/T1w.nii.gz" -interp linear - -quiet | \
                    mrcalc - $_label_thresh -ge $_label_int -mult \
                    "Karawun/sub-${participant}/labels/afMRI_${_label_basename}.nii.gz" -force -quiet
            done
        fi
    else
        # Karawun prep runs at the end of every pipeline run now, so this should
        # not normally fire. It still can on an export-only tree (no FWT output
        # to prepare from), or if prep failed and left no marker.
        echo "Karawun/sub-${participant}/T1w.nii.gz not found — skipping Karawun fMRI labels."
        echo "  Karawun prep normally runs at the end of a full pipeline run. To do it now:"
        echo "    KUL_karawun_prepare.sh -p ${participant} -t <1=tumor|2=DBS ET|3=DBS Parkinson>"
        echo "  then re-run this -R to add the fMRI activation labels."
    fi

    if [ $_use_dropdir -eq 0 ]; then
        for _spm in "$globalresultsdir/Melodic/"*.nii.gz "$globalresultsdir/Melodic/"*.nii; do
            [ -f "$_spm" ] || continue
            _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
            _queue_spm "$_spm" "$_spmname" thresholded
        done
        _flush_spm_batch
    fi

    # ── Lesion & DSC perfusion → PACS ───────────────────────────────────────
    # Same renderer as the fMRI maps, pointed at a separate output folder so
    # these don't get mixed in with the SPM/Melodic series on PACS.
    # Skipped entirely when drop folders are in use — the operator has said
    # explicitly what to export.
    _extra_names=(); _extra_files=(); _extra_thresh=()

    if [ $_use_dropdir -eq 0 ] && _lesion_pacs=$(KUL_resolve_lesion); then
        # a mask, so any threshold in (0,1) selects it; 0.5 is the obvious one
        _extra_names+=("Lesion"); _extra_files+=("$_lesion_pacs"); _extra_thresh+=("0.5")
    fi

    # Fixed thresholds, not the auto max/3 used for the fMRI maps. These are
    # NAWM-normalised ratios, so a threshold has a fixed clinical meaning that
    # max/3 would throw away: 1.75 is the conventional high-grade glioma rCBV
    # cutoff, and 1.0 on nrCBF is simply "above contralesional normal WM".
    # -T still overrides both.
    _perf_dir="$globalresultsdir/Perfusion"
    if [ $_use_dropdir -eq 0 ]; then
        for _pm in nrCBV_corrected:1.75 nrCBF:1.0; do
            _pf="$_perf_dir/sub-${participant}_${_pm%%:*}.nii.gz"
            [ -f "$_pf" ] || continue
            _extra_names+=("${_pm%%:*}"); _extra_files+=("$_pf"); _extra_thresh+=("${_pm##*:}")
        done
    fi

    if [ ${#_extra_names[@]} -gt 0 ]; then
        _saved_png="$spm_resultsdir_png"; _saved_dcm="$spm_resultsdir_dcm"
        # _dcm_spm_set carries the user's fMRI map selection. These maps are not
        # in it, and a non-empty set means "only these", so it would silently
        # exclude them -- clear it for this pass and put it back afterwards.
        _saved_sel=("${!_dcm_spm_set[@]}")
        unset _dcm_spm_set; declare -A _dcm_spm_set=()

        spm_resultsdir_png="$globalresultsdir/Clinical_figures_${_ulsuffix}"
        spm_resultsdir_dcm="$globalresultsdir/PACS/Clinical_${_ulsuffix}"
        mkdir -p "$spm_resultsdir_png" "$spm_resultsdir_dcm"

        # Thresholds are set in the PARENT here, before anything is queued:
        # _render_one_spm only reads spm_thresh_map, and a backgrounded child's
        # write would never reach the parent anyway.
        for _i in "${!_extra_names[@]}"; do
            if [ -n "${_extra_thresh[$_i]}" ] && [ -z "$spm_thresh_override" ] && \
               [ -z "${spm_thresh_map[${_extra_names[$_i]}]+x}" ]; then
                spm_thresh_map["${_extra_names[$_i]}"]="${_extra_thresh[$_i]}"
            fi
            echo "PACS extra: ${_extra_names[$_i]} <- $(basename "${_extra_files[$_i]}")"
            _queue_spm "${_extra_files[$_i]}" "${_extra_names[$_i]}" thresholded
        done
        _flush_spm_batch

        spm_resultsdir_png="$_saved_png"; spm_resultsdir_dcm="$_saved_dcm"
        unset _dcm_spm_set; declare -A _dcm_spm_set=()
        for _k in "${_saved_sel[@]}"; do _dcm_spm_set["$_k"]=1; done
    fi

    # ── Quantitative series → PACS ──────────────────────────────────────────
    # Unlike everything above, these are not screenshots. The NIfTI's own voxels
    # are written as a 16-bit DICOM series carrying RescaleSlope/Intercept, so
    # an ROI drawn on PACS reads real units (rCBV, ALFF, ReHo, FA, ...) rather
    # than display colours. No mrview, no xvfb, no threshold, no underlay.
    if [ $make_dcm -eq 1 ] && [ -n "$donor_dcm" ]; then
        _quant_dcm_root="$globalresultsdir/PACS/Quant"

        # Explicit wins: whatever is in series_quantitative/ is what gets
        # exported. With nothing dropped there, fall back to every DSC map in
        # Perfusion/ -- those are the ones anyone actually measures, and the
        # scanner gives you ADC the same way. (PCASL deliberately excluded:
        # not processed yet.)
        _quant_files=()
        for _qf in "$_drop_quant"/*.nii.gz "$_drop_quant"/*.nii; do
            [ -f "$_qf" ] && _quant_files+=("$_qf")
        done
        if [ ${#_quant_files[@]} -eq 0 ]; then
            for _qf in "$globalresultsdir/Perfusion"/*.nii.gz "$globalresultsdir/Perfusion"/*.nii; do
                [ -f "$_qf" ] || continue
                case "$(basename "$_qf")" in
                    *_mask.nii*|*NAWM*) continue ;;   # masks, not measurements
                esac
                _quant_files+=("$_qf")
            done
            [ ${#_quant_files[@]} -gt 0 ] && \
                echo "No files in series_quantitative/ — exporting ${#_quant_files[@]} DSC map(s) from Perfusion/"
        fi

        for _qf in "${_quant_files[@]}"; do
            [ -f "$_qf" ] || continue
            mkdir -p "$_quant_dcm_root"
            _qname=$(basename "$_qf"); _qname=${_qname%.nii.gz}; _qname=${_qname%.nii}
            # "<name>.label.nii.gz" marks an integer label/segmentation map,
            # which must be written without any rescaling.
            _qmode="-q"
            case "$_qname" in
                *.label) _qmode="--label"; _qname="${_qname%.label}" ;;
            esac
            # These inherit the donor's Frame of Reference, so a map in a
            # different space would land misaligned on PACS with nothing to
            # indicate it. Warn rather than silently export.
            if ! mrinfo "$_qf" -quiet >/dev/null 2>&1; then
                echo "WARNING: $_qf is not readable as an image — skipping"
                _pacs_note_fail "convert" "$_qname" "unreadable input"
                continue
            fi
            if [ "$(mrinfo "$_qf" -size)" != "$(mrinfo "$underlay" -size)" ] || \
               [ "$(mrinfo "$_qf" -spacing)" != "$(mrinfo "$underlay" -spacing)" ]; then
                echo "NOTE: $(basename "$_qf") is not on the underlay grid"
                echo "      ($(mrinfo "$_qf" -size) @ $(mrinfo "$_qf" -spacing) vs $(mrinfo "$underlay" -size) @ $(mrinfo "$underlay" -spacing))."
                echo "      That is fine if it shares the scanner coordinate frame — geometry is taken"
                echo "      from the file itself — but it will be misaligned on PACS if it is not registered."
            fi
            _pacs_register "quant:$_qname"
            _qsnum=$(_pacs_series_number "quant:$_qname" "NONE")
            _qnopt=""; [ -n "$_qsnum" ] && _qnopt="-n $_qsnum"
            _qdir="$_quant_dcm_root/$_qname"
            if [ -n "$(ls "$_qdir"/*.dcm 2>/dev/null | head -1)" ]; then
                echo "Skipping quantitative DICOMs for ${_qname} (already exist)"
                continue
            fi
            mkdir -p "$_qdir"
            # Extra display headroom for CBV maps only. Their bright tail is
            # choroid plexus -- real, very vascular anatomy that saturated as
            # white blobs at a plain p99.5. Deliberately not applied to the
            # other DSC maps: K2 and MTT have heavier tails still, but widening
            # their window only flattens them, so this cannot be inferred from
            # the data and is keyed on the map instead.
            _qhead=""
            case "$_qname" in
                *CBV*|*cbv*) _qhead="--window-headroom 1.3" ;;
            esac
            echo "Making quantitative dicoms in $_qdir (${_qmode}${_qhead:+, +30% window headroom}, SeriesNumber=${_qsnum:-<none>})"
            _run_nii2dcm "quant_${_qname}" -s "$_qname" $_qnopt $_qmode $_qhead \
                "$_qf" "$donor_dcm" "$_qdir"
        done
    fi

    # ── Run summary ─────────────────────────────────────────────────────────
    # Failures used to be printed as passing WARNING lines and the run always
    # exited 0, so a half-finished PACS export looked identical to a good one.
    if [ -s "$_pacs_fail_file" ]; then
        echo ""
        echo "=============================================================="
        echo " PACS/figure generation finished WITH FAILURES"
        echo "=============================================================="
        printf '%-8s %-48s %s\n' "STAGE" "SERIES" "DETAIL"
        cat "$_pacs_fail_file"
        echo "--------------------------------------------------------------"
        echo " $(grep -c '' "$_pacs_fail_file") failure(s). Logs: $_pacs_logdir"
        echo " Re-running with -R retries only what is missing or incomplete."
        echo "=============================================================="
        rm -f "$_pacs_fail_file"
        exit 1
    fi
    echo ""
    echo "PACS/figure generation completed with no failures. Logs: $_pacs_logdir"
    rm -f "$_pacs_fail_file"

    exit

fi

# --- functions ---
function KUL_check_redo {
    if [ $redo -eq 1 ];then

        # Questions follow the same order as the main processing steps.

        # Step 3 — tumor segmentation (types 1 and 2 only)
        if [ $type -lt 3 ]; then
            read -p "Redo: tumor segmentation? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -rf ${cwd}/KUL_LOG/sub-${participant}_anat_*.done >/dev/null 2>&1
                rm -rf $derivativesdir/KUL_anat_biascorrect >/dev/null 2>&1
                rm -rf $derivativesdir/KUL_anat_register_rigid >/dev/null 2>&1
                rm -rf $derivativesdir/KUL_anat_segment_tumor >/dev/null 2>&1
                rm -rf $globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz >/dev/null 2>&1
            fi
        fi

        # Step 2 — anatomical registration (all types)
        # Clears the derivative as well, so this is a genuine re-registration.
        # To only refresh the copies in RESULTS/Anat without paying for ANTs
        # again, use the "refresh RESULTS" question at the end instead.
        read -p "Redo: anatomical registration? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -f ${cwd}/KUL_LOG/sub-${participant}_anat_reg.done >/dev/null 2>&1
            rm -rf $derivativesdir/KUL_anat_register_rigid >/dev/null 2>&1
        fi

        # Step 4 — fmriprep (all types)
        read -p "Redo: fmriprep? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -rf ${cwd}/fmriprep/sub-${participant} >/dev/null 2>&1
            rm -rf ${cwd}/fmriprep_work >/dev/null 2>&1
            rm -f ${cwd}/fmriprep/sub-${participant}.html >/dev/null 2>&1
        fi

        # Step 5 — dwiprep (all types)
        read -p "Redo: KUL_dwiprep? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -f ${cwd}/KUL_LOG/sub-${participant}_run_dwiprep.txt >/dev/null 2>&1
            rm -rf ${cwd}/dwiprep/sub-${participant} >/dev/null 2>&1
            rm -rf $derivativesdir/synb0 >/dev/null 2>&1
        fi

        # Step 6 — Gd contrast T1w subtraction (only when a contrast T1w exists)
        # Pure mrcalc on two files already in RESULTS/Anat, so the marker is all
        # there is to clear.
        if [ ${ncT1w:-0} -gt 0 ] || [ ${ncT1w:-0} -eq -1 ]; then
            read -p "Redo: Gd contrast subtraction? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -f ${cwd}/KUL_LOG/sub-${participant}_cT1w_subtraction.done >/dev/null 2>&1
            fi
        fi

        # Step 8 — SPM and Melodic (only if fMRI data present)
        if [ $n_fMRI -gt 0 ]; then
            read -p "Redo: SPM? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -f ${cwd}/KUL_LOG/sub-${participant}_SPM.done >/dev/null 2>&1
                rm -fr $derivativesdir/SPM/* >/dev/null 2>&1
                rm -fr ${cwd}/RESULTS/sub-${participant}/SPM/* >/dev/null 2>&1
                rm -fr ${cwd}/RESULTS/sub-${participant}/SPM_all/* >/dev/null 2>&1
            fi
            read -p "Redo: Melodic? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -f ${cwd}/KUL_LOG/sub-${participant}_melodic.done >/dev/null 2>&1
                rm -fr $derivativesdir/FSL_melodic/* >/dev/null 2>&1
                rm -fr ${cwd}/RESULTS/sub-${participant}/Melodic/* >/dev/null 2>&1
            fi
        fi

        # Step 9 — VBG (types 1, 2, 3 only)
        if [ $vbg -gt 0 ]; then
            read -p "Redo: KUL_VBG? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -f ${cwd}/KUL_LOG/sub-${participant}_VBG.log >/dev/null 2>&1
                rm -fr $vbg_dir/* >/dev/null 2>&1
            fi
        fi

        # Step 9c — DSC perfusion (only if perfusion data present)
        #
        # Two levels, because they cost very different things. KUL_dsc_fit
        # returns early when fit/rCBV_corrected.nii.gz exists, so clearing only
        # the marker re-exports the existing maps into RESULTS/Perfusion in
        # seconds. Deleting the derivative throws that away and refits from the
        # raw DSC series, which is the expensive part -- worth doing when the fit
        # itself is what you are unhappy with, and a waste otherwise.
        if [ ${n_dsc:-0} -gt 0 ]; then
            read -p "Redo: DSC perfusion? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -f ${cwd}/KUL_LOG/sub-${participant}_DSC.done >/dev/null 2>&1
                rm -fr ${cwd}/RESULTS/sub-${participant}/Perfusion/* >/dev/null 2>&1
                read -p "  Also refit from the raw DSC series? Slow; say n to just re-export the existing fit. (y/n) " answ
                if [[ "$answ" == "y" ]]; then
                    rm -fr $derivativesdir/KUL_dsc_perfusion >/dev/null 2>&1
                fi
            fi
        fi

        # Step 11 — dwiprep_anat (all types)
        read -p "Redo: KUL_dwiprep_anat? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -f ${cwd}/KUL_LOG/sub-${participant}_dwiprep_anat.done >/dev/null 2>&1
        fi

        # Step 12 — dwiprep_MNI (all types)
        read -p "Redo: KUL_dwiprep_MNI? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -f ${cwd}/KUL_LOG/sub-${participant}_dwiprep_MNI.done >/dev/null 2>&1
        fi

        # Step 13 — DTI-ALPS (type 7 only)
        if [ $alps -eq 1 ]; then
            read -p "Redo: KUL_calc_DTI_ALPS? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -rf ${cwd}/KUL_dwiprep/sub-${participant}/sub-${participant}/DTI_ALPS >/dev/null 2>&1
                rm -f ${cwd}/KUL_LOG/sub-${participant}_dti_ALPS.done >/dev/null 2>&1
            fi
        fi

        # Step 14 — FWT (types 1–6)
        if [ $fwt -eq 1 ]; then
            read -p "Redo: KUL_FWT? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                rm -f ${cwd}/KUL_LOG/sub-${participant}_FWT.done >/dev/null 2>&1
                # KUL_FWT writes to $derivativesdir/FWT, not $derivativesdir/KUL_FWT. The old
                # path matched nothing, so "redo" cleared the .done gate and left every
                # <bundle>_VOIs.done marker and .tck in place -- make_VOIs then skipped every
                # bundle as "already generated" and the rerun reproduced the previous result
                # exactly, including after a recipe change.
                rm -fr $derivativesdir/FWT/* >/dev/null 2>&1
            fi
        fi

        # Step 15 — Karawun/Brainlab folder
        # Marker only. Deleting Karawun/sub-*/ itself would take the donor DICOM
        # in Karawun/sub-*/DICOM/ with it, which nothing regenerates -- and prep
        # overwrites its own outputs with -force anyway.
        read -p "Redo: Karawun/Brainlab folder? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -f ${cwd}/KUL_LOG/sub-${participant}_karawun_prepare.done >/dev/null 2>&1
        fi

        # Step 17 — figures
        read -p "Redo: figures? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -f ${cwd}/KUL_LOG/sub-${participant}_figures.done >/dev/null 2>&1
        fi

        # Redo the results-production block: repopulate RESULTS/ and Karawun/
        # from analyses that have already run.
        #
        # Asked last, and separately, because it answers a different question.
        # The steps above are "I was unhappy with this analysis, run it again".
        # This one is "the analyses are fine, but RESULTS/ and Karawun/ do not
        # reflect them" -- cleaned out by hand, or simply never being sure when
        # they get repopulated.
        #
        # It only ever clears a marker when that step's *input* still exists in
        # BIDS/derivatives, which is what makes it cheap: every one of these
        # steps re-uses an existing result rather than recomputing it (the GLM
        # re-uses its spmT maps, melodic its decomposition, DSC its fit), so
        # clearing the marker re-exports into RESULTS instead of re-analysing.
        # Where the input is gone, the marker is left alone and said so: that is
        # a genuine re-analysis and belongs to its own question above.
        #
        #   marker | dependency that makes re-export cheap | what it repopulates
        local _refresh=(
            "anat_reg|$derivativesdir/KUL_anat_register_rigid/T1w.nii.gz|Anat/ (T1w, *reg2_T1w)"
            "cT1w_subtraction|$globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz|Anat/cT1w_T1w_subtracted"
            "SPM|$(ls -d $derivativesdir/SPM/RESULTS/stats_* 2>/dev/null | head -1)|SPM/ and SPM_all/"
            "melodic|$(ls $derivativesdir/FSL_melodic/stats_*/stats/thresh_zstat1.nii.gz 2>/dev/null | head -1)|Melodic/"
            "DSC|$derivativesdir/KUL_dsc_perfusion/fit/rCBV_corrected.nii.gz|Perfusion/"
            "rsfMRI_networks|$kulderivativesdir/rsfMRI_networks/analysis|rsfMRI_Networks/"
            "karawun_prepare|$derivativesdir/FWT/sub-${participant}_TCKs_output|Karawun/ (T1w, tck, labels)"
        )
        local _rspec _rm _rdep _rwhat _ready=() _notready=()
        for _rspec in "${_refresh[@]}"; do
            IFS='|' read -r _rm _rdep _rwhat <<< "$_rspec"
            if [ -n "$_rdep" ] && [ -e "$_rdep" ]; then
                _ready+=("$_rm|$_rwhat")
            else
                _notready+=("$_rm|$_rwhat")
            fi
        done

        echo ""
        echo "  Repopulate RESULTS/ and Karawun/ from analyses that already ran?"
        echo "  This re-exports; it re-analyses nothing."
        if [ ${#_ready[@]} -gt 0 ]; then
            echo "    will rebuild:"
            for _rspec in "${_ready[@]}"; do
                IFS='|' read -r _rm _rwhat <<< "$_rspec"
                printf '      %-28s (from %s)\n' "$_rwhat" "$_rm"
            done
        fi
        if [ ${#_notready[@]} -gt 0 ]; then
            echo "    cannot rebuild - the analysis output is gone, so these would be"
            echo "    real re-runs; use their own questions above:"
            for _rspec in "${_notready[@]}"; do
                IFS='|' read -r _rm _rwhat <<< "$_rspec"
                printf '      %-28s (%s)\n' "$_rwhat" "$_rm"
            done
        fi
        echo "    never touched: Lesion/, PACS_input/, DICOM/ - your own files"
        echo "                   Tracto/ and TRK/ re-sync from FWT every run anyway"
        echo "    REPORT/        rebuilt every run from whatever sources exist"
        echo "                   (fmriprep html, eddy QC, FA overlays, tract summary);"
        echo "                   the -F/-R figures come back with -F/-R"
        if [ ${#_ready[@]} -gt 0 ]; then
            read -p "  Repopulate RESULTS and Karawun? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                for _rspec in "${_ready[@]}"; do
                    IFS='|' read -r _rm _rwhat <<< "$_rspec"
                    # Marker only, never the derivative: that is what keeps this a
                    # re-export. "Redo: <step>?" above is the one that pays again.
                    rm -f ${cwd}/KUL_LOG/sub-${participant}_${_rm}.done >/dev/null 2>&1
                done
                echo "    ${#_ready[@]} step(s) will re-export during this run"
            fi
        fi

    fi
}

# Report RESULTS/Karawun output that a completed step should have produced but
# which is no longer on disk -- the state you get by deleting a results folder to
# "start fresh" while KUL_LOG still holds the .done markers. Every step gates on
# its marker alone and never checks its own output, so that combination otherwise
# produces a run that skips everything, reports success, and leaves RESULTS empty.
#
# The marker is the discriminator between "deleted" and "never generated", and it
# is a reliable one: each step puts its applicability test OUTSIDE the marker
# check (KUL_fmriproc only reaches SPM.done when n_fMRI>0, the DSC block only when
# n_dsc>0) and only touches the marker on success. So a marker existing means
# applicable AND ran AND succeeded:
#
#   marker absent                -> never generated. Normal; nothing is said.
#   marker present, output there -> fine.
#   marker present, output gone  -> it existed once and does not now.  <- reported
#
# Repair is deliberately nothing more than removing the marker. This runs before
# every processing step, so the pipeline's own code regenerates the output moments
# later in the same invocation, applying the same transforms in the same order
# (Tracto maps regrid onto Anat/T1w.nii.gz, the GM mask is recomputed from
# fmriprep's dseg, ...). Nothing here re-implements a processing step, so nothing
# here can drift away from one.
function KUL_verify_results {

    local _spec _m _p _kind _label _marker _missing _repairable=() _n=0 answ

    # Type 3's lesion mask is not a pipeline product: it is hand-drawn and copied
    # in by the user ("put lesion.nii.gz in RESULTS/sub-{participant}/Lesion").
    # No marker, no derivative copy, nothing that can regenerate it -- so this is
    # a stop, not a repair offer. Continuing gives a silently lesion-free run: no
    # VBG lesion, no Karawun lesion label, no PACS lesion overlay, no error.
    if [ $type -eq 3 ] && [ ! -f "$globalresultsdir/Lesion/lesion.nii.gz" ]; then
        echo "" >&2
        echo "ERROR: -t 3 requires a manual lesion mask and it is not present:" >&2
        echo "         ${globalresultsdir#$cwd/}/Lesion/lesion.nii.gz" >&2
        echo "       Nothing in this pipeline can regenerate it - it is the mask you" >&2
        echo "       drew and copied in. Restore it from wherever it came from." >&2
        echo "       Continuing would silently produce a lesion-free run." >&2
        echo "" >&2
        exit 2
    fi

    # marker | sentinel path | file|dir | description
    #
    # A directory sentinel is "exists and is non-empty": these steps name their
    # outputs after tasks/networks/bundles, so there is no fixed filename to test.
    local _checks=(
        "anat_reg|$globalresultsdir/Anat/T1w.nii.gz|file|anatomical registration (Anat/)"
        "cT1w_subtraction|$globalresultsdir/Anat/cT1w_T1w_subtracted.nii.gz|file|Gd contrast subtraction"
        "SPM|$globalresultsdir/SPM|dir|fMRI GLM (SPM/)"
        "melodic|$globalresultsdir/Melodic|dir|melodic (Melodic/)"
        "DSC|$globalresultsdir/Perfusion|dir|DSC perfusion (Perfusion/)"
        "karawun_prepare|$cwd/Karawun/sub-${participant}/T1w.nii.gz|file|Karawun/Brainlab folder"
    )

    for _spec in "${_checks[@]}"; do
        IFS='|' read -r _m _p _kind _label <<< "$_spec"
        _marker="${cwd}/KUL_LOG/sub-${participant}_${_m}.done"
        [ -f "$_marker" ] || continue          # never generated -- nothing to report
        _missing=0
        if [ "$_kind" = "dir" ]; then
            { [ -d "$_p" ] && [ -n "$(ls -A "$_p" 2>/dev/null)" ]; } || _missing=1
        else
            [ -f "$_p" ] || _missing=1
        fi
        [ $_missing -eq 0 ] && continue

        if [ $_n -eq 0 ]; then
            echo ""
            echo "=============================================================="
            echo " RESULTS/Karawun integrity check"
            echo "=============================================================="
        fi
        _n=$((_n + 1))
        echo "  ${_label}"
        echo "    KUL_LOG/sub-${participant}_${_m}.done records this step as done,"
        echo "    but ${_p#$cwd/} is missing or empty."
        _repairable+=("${_m}|${_label}")
    done

    # T1w_GM has no marker of its own -- it belongs to the fmriprep step, which
    # gates on fmriprep/sub-*.html. That step now rebuilds it (and the fmriprep
    # report symlink, and Anat/T1w_fmriprep.nii.gz) whenever the source exists
    # rather than only when fmriprep itself runs, so this is a note, not a repair.
    if [ -f "$cwd/fmriprep/sub-${participant}.html" ] && \
       [ ! -f "$globalresultsdir/Anat/T1w_GM.nii.gz" ]; then
        if [ $_n -eq 0 ]; then
            echo ""
            echo "=============================================================="
            echo " RESULTS/Karawun integrity check"
            echo "=============================================================="
        fi
        _n=$((_n + 1))
        echo "  cortical GM mask (Anat/T1w_GM.nii.gz)"
        echo "    Missing - it will be rebuilt from fmriprep's dseg later in this run."
        echo "    It is what continuous PACS overlays window against; without it they"
        echo "    fall back to a p98 top."
    fi

    [ $_n -eq 0 ] && return 0

    if [ ${#_repairable[@]} -gt 0 ]; then
        echo ""
        echo " Clearing a marker makes this run regenerate that output using the"
        echo " pipeline's own code. Anatomical and Karawun steps are quick; the GLM,"
        echo " melodic and DSC are full re-runs."
        echo ""
        if [ ${redo:-0} -eq 1 ]; then
            # -r is the explicit "ask me about redoing things" mode and now has a
            # question for every check above, including a "Refresh RESULTS?" one
            # that covers the cheap copy/derive steps in a single answer. Prompting
            # here too would just ask everything twice.
            echo " You passed -r, so the questions that follow will cover these."
        elif [ -t 0 ]; then
            for _spec in "${_repairable[@]}"; do
                IFS='|' read -r _m _label <<< "$_spec"
                read -p "  Regenerate ${_label}? (y/n) " answ
                if [[ "$answ" == "y" ]]; then
                    rm -f "${cwd}/KUL_LOG/sub-${participant}_${_m}.done" >/dev/null 2>&1
                    echo "    marker cleared - will be regenerated during this run"
                else
                    echo "    left as is"
                fi
            done
        else
            # Batch/nohup runs have no tty: `read` would take the next line of the
            # script's stdin or block. Report and change nothing.
            echo " Not an interactive terminal, so nothing was changed. Re-run from a"
            echo " terminal to be asked, or clear the markers yourself:"
            for _spec in "${_repairable[@]}"; do
                IFS='|' read -r _m _label <<< "$_spec"
                echo "   rm -f KUL_LOG/sub-${participant}_${_m}.done   # ${_label}"
            done
        fi
    fi
    echo "=============================================================="
    echo ""
    return 0
}

# Verify the fMRI python env is actually usable, before anything expensive runs.
#
# Deliberately not modelled on the lore_sd/scilpy checks above, which test
# `conda env list | grep -qx <name>`. That test would have PASSED the failure
# this exists for: the env was present and had every package. What broke was
# activation -- conda's bin/ was not on PATH in that shell -- and the fMRI step
# reported it as "missing required packages", which sends you to inspect the env
# instead of the shell. So this does what the step does: activate, then import.
#
# Runs in a subshell so a successful activation does not leak into the pipeline's
# own environment.
function KUL_check_pyfmri_env {
    local _env="${pyfmri_env_override:-$KUL_PYFMRI_ENV}"
    local _err _rc

    # Bootstrap in this shell as well, not only inside the subshell below, or the
    # "does the env exist" branch cannot run `conda env list` and misreports a
    # broken env as a missing one. Sourcing conda.sh only defines the function --
    # unlike `conda shell.bash hook` it activates nothing, so nothing leaks.
    KUL_conda_bootstrap || true

    _err=$( ( KUL_conda_bootstrap && conda activate "$_env" && \
              python -c "import nilearn, nibabel, numpy, pandas" ) 2>&1 )
    _rc=$?
    if [ $_rc -eq 0 ]; then
        echo " fMRI python env '$_env': usable"
        return 0
    fi

    echo "" >&2
    echo "ERROR: the fMRI python env '$_env' is not usable, and this run needs it:" >&2
    echo "       the task-fMRI GLM (-E nilearn), melodic, and the rsfMRI networks (-N)." >&2
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$_env"; then
        echo "  The env EXISTS, so this is an activation or package problem rather than a" >&2
        echo "  missing env. If conda is not initialised in this shell, run KUL_Linux_setup's" >&2
        echo "  bashrc section (./setup_environment.sh --only bashrc), or export" >&2
        echo "  KUL_CONDA_BASE=<conda root>." >&2
    else
        echo "  The env does NOT exist. Create it with:" >&2
        echo "    ./setup_environment.sh --only env-pyfmri     (in KUL_Linux_setup)" >&2
    fi
    [ -n "$_err" ] && echo "  Detail: $(echo "$_err" | grep -v '^$' | tail -2 | tr '\n' ' ')" >&2
    echo "  Use -y <env> if you keep it under a different name." >&2
    echo "" >&2
    return 1
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

function KUL_pad_anat {
    local pad_check=${cwd}/KUL_LOG/sub-${participant}_pad_anat.done
    if [ -f $pad_check ]; then
        echo "Anatomical FOV padding already done"
        return
    fi
    local pad_voxels=10
    local anat_dir="${cwd}/BIDS/sub-${participant}/anat"
    local padded_any=0
    for nii in "${anat_dir}"/*.nii.gz; do
        [ -f "$nii" ] || continue
        # get image dimensions and check if any brain voxel is within pad_voxels of any edge
        read -r nx ny nz <<< $(mrinfo "$nii" -size | awk '{print $1, $2, $3}')
        # create a quick brain mask and check proximity to FOV edge
        local tmpdir=$(mktemp -d)
        mrthreshold "$nii" -abs 0 "${tmpdir}/mask.mif" -force -quiet 2>/dev/null
        # check bounding box of non-zero voxels
        read -r x0 x1 y0 y1 z0 z1 <<< $(mrstats "${tmpdir}/mask.mif" -mask "${tmpdir}/mask.mif" -output bbox 2>/dev/null | awk '{print $1,$2,$3,$4,$5,$6}')
        rm -rf "$tmpdir"
        local needs_pad=0
        for val in $x0 $y0 $z0; do
            [ -n "$val" ] && (( $(echo "$val < $pad_voxels" | bc -l) )) && needs_pad=1
        done
        [ -n "$x1" ] && (( $(echo "$nx - $x1 < $pad_voxels" | bc -l) )) && needs_pad=1
        [ -n "$y1" ] && (( $(echo "$ny - $y1 < $pad_voxels" | bc -l) )) && needs_pad=1
        [ -n "$z1" ] && (( $(echo "$nz - $z1 < $pad_voxels" | bc -l) )) && needs_pad=1
        if [ $needs_pad -eq 1 ]; then
            echo "  Padding ${nii} (brain near FOV edge, adding ${pad_voxels} voxels)"
            mrgrid "$nii" pad -uniform $pad_voxels "${nii%.nii.gz}_padded.nii.gz" -force -quiet
            mv "${nii%.nii.gz}_padded.nii.gz" "$nii"
            padded_any=1
        fi
    done
    [ $padded_any -eq 0 ] && echo "  No anatomical images required FOV padding"
    touch $pad_check
}

function KUL_convert2bids {
    # convert the DICOM to BIDS
    if [ ! -d "BIDS/sub-${participant}" ];then
        if [ $d_flag -eq 1 ]; then
            KUL_dcm2bids.sh -d $dicomzip -p ${participant} -c study_config/sequences.txt -e -v
        else
            echo "Error: no dicom zip file given."
        fi
    else
        echo "BIDS conversion already done"
    fi
}

function KUL_check_data {
    
    mkdir -p $globalresultsdir
    echo -e "\n\nAn overview of the BIDS data:"
    bidsdir="BIDS/sub-$participant"
    T1w=($(find $bidsdir -name "*T1w.nii.gz" ! -name "*gadolinium*" -type f ))
    nT1w=${#T1w[@]}
    echo "  number of non-contrast T1w: $nT1w"

    cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
    ncT1w=${#cT1w[@]}
    # check if it has been moved away
    if [ $ncT1w -eq 0 ]; then
        if [ -f "$kulderivativesdir/sub-${participant}/cT1w/sub-${participant}_ce-gadolinium_T1w.nii.gz" ]; then
            ncT1w=-1
        fi
    fi
    echo "  number of contrast enhanced T1w: $ncT1w"


    FLAIR=($(find $bidsdir -name "*FLAIR.nii.gz" -type f ))
    nFLAIR=${#FLAIR[@]}
    echo "  number of FLAIR: $nFLAIR"
    FGATIR=($(find $bidsdir -name "*FGATIR.nii.gz" -type f ))
    nFGATIR=${#FGATIR[@]}
    echo "  number of FGATIR: $nFGATIR"
    T2w=($(find $bidsdir -name "*T2w.nii.gz" -type f ))
    nT2w=${#T2w[@]}
    echo "  number of T2w: $nT2w"
    SWI=($(find $bidsdir -name "*_SWI.nii.gz" ! -name "*SWIp*" -type f | sort))
    nSWI=${#SWI[@]}
    echo "  number of SWI magnitude: $nSWI"

    _swip_all=($(find $bidsdir -name "*_SWIp.nii.gz" -type f | sort))
    nSWIp=${#_swip_all[@]}
    if [ $nSWIp -eq 1 ]; then
        SWIp=${_swip_all[0]}
    elif [ $nSWIp -gt 1 ]; then
        if [ -t 0 ]; then
            echo "Multiple SWIp volumes found:"
            for _i in "${!_swip_all[@]}"; do
                echo "  $((_i+1)). $(basename ${_swip_all[$_i]})"
            done
            read -r -p "Select SWIp volume to use (Enter = first): " _sel
            if [[ "$_sel" =~ ^[0-9]+$ ]] && [ "$_sel" -ge 1 ] && [ "$_sel" -le $nSWIp ]; then
                SWIp=${_swip_all[$((_sel-1))]}
            else
                SWIp=${_swip_all[0]}
            fi
        else
            SWIp=${_swip_all[0]}
        fi
        nSWIp=1
    fi
    echo "  number of SWI phase: $nSWIp"

    DIR=($(find $bidsdir -name "*DIR.nii.gz" -type f ))
    nDIR=${#DIR[@]}
    echo "  number of DIR: $nDIR"

    # MP2RAGE: user selects which volume (TI) to use; auto-selects highest TI if non-interactive
    MP2RAGE=""
    nMP2RAGE=0
    _mp2rage_files=($(find $bidsdir -name "*MP2RAGE*.nii.gz" -type f | sort))
    if [ ${#_mp2rage_files[@]} -eq 1 ]; then
        nMP2RAGE=1
        MP2RAGE=${_mp2rage_files[0]}
    elif [ ${#_mp2rage_files[@]} -gt 1 ]; then
        nMP2RAGE=1
        if [ -t 0 ]; then
            echo "Multiple MP2RAGE volumes found:"
            for _i in "${!_mp2rage_files[@]}"; do
                _json="${_mp2rage_files[$_i]%.nii.gz}.json"
                _ti=$(python3 -c "import json; d=json.load(open('$_json')); print(d.get('TriggerDelayTime','?'))" 2>/dev/null || echo "?")
                echo "  $((_i+1)). $(basename ${_mp2rage_files[$_i]}) (TI=${_ti}ms)"
            done
            read -r -p "Select MP2RAGE volume to use (Enter = auto highest TI / INV2): " _sel
            if [[ "$_sel" =~ ^[0-9]+$ ]] && [ "$_sel" -ge 1 ] && [ "$_sel" -le ${#_mp2rage_files[@]} ]; then
                MP2RAGE=${_mp2rage_files[$((_sel-1))]}
            fi
        fi
        if [ -z "$MP2RAGE" ]; then
            _best_ti=0
            for _f in "${_mp2rage_files[@]}"; do
                _json="${_f%.nii.gz}.json"
                _ti=$(python3 -c "import json; d=json.load(open('$_json')); print(d.get('TriggerDelayTime',0))" 2>/dev/null || echo 0)
                if (( $(echo "$_ti > $_best_ti" | bc -l) )); then _best_ti=$_ti; MP2RAGE=$_f; fi
            done
        fi
    fi
    echo "  MP2RAGE selected: ${MP2RAGE:-(none)}"

    # check the T1w
    if [ $nT1w -eq 0 ]; then
        echo "No T1w (without Gd) found. Fmriprep will not run."
        echo " Is the BIDS dataset correct?"
        read -p "Are you sure you want to continue? (y/n)? " answ
        if [[ "$answ" == "n" ]]; then
            exit 1
        fi
    fi 

    # check the cT1w
    if [ $ncT1w -eq -1 ]; then
        echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
        echo " At least one is missing."
        echo " The contrast T1w has been moved to the derivatives folder, due to previous processing."
        read -p "Do you want to restore it? (y/n)? " answ
        if [[ "$answ" == "y" ]]; then
            cp -f $kulderivativesdir/sub-${participant}/cT1w/sub-${participant}_ce-gadolinium_T1w.* $bidsdir/anat
            cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
            ncT1w=${#cT1w[@]}
        fi
    fi

    # check hd-glio-auto requirements
    if [ $hdglio -eq 1 ]; then
        if [ $nT1w -lt 1 ] || [ $ncT1w -lt 1 ] || [ $nT2w -lt 1 ] || [ $nT1w -lt 1 ]; then
            echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
            echo " At least one is missing. Is the BIDS dataset correct?"
            read -p "Are you sure you want to continue? (y/n)? " answ
            if [[ "$answ" == "n" ]]; then
                exit 1
            fi
        fi
    fi 

    # check the BIDS
    find_fmri=($(find ${cwd}/BIDS/sub-${participant} -name "*_bold.nii.gz"))
    n_fMRI=${#find_fmri[@]}
    if [ $n_fMRI -eq 0 ]; then
        echo "WARNING: no fMRI data"
    fi

    find_dwi=($(find ${cwd}/BIDS/sub-${participant} -name "*_dwi.nii.gz"))
    n_dwi=${#find_dwi[@]}
    if [ $n_dwi -eq 0 ]; then
        echo "WARNING: no dwi data"
    fi

    find_dsc=($(find ${cwd}/BIDS/sub-${participant} -name "*_dsc.nii.gz"))
    n_dsc=${#find_dsc[@]}
    echo "  number of DSC perfusion series: $n_dsc"
    echo -e "\n\n"

}

function KUL_rigid_register {
    warp_field="${registeroutputdir}/${source_mri_label}_reg2_T1w"
    output_mri="${globalresultsdir}/Anat/${source_mri_label}_reg2_T1w.nii.gz"
    #echo "Rigidly registering $source_mri to $target_mri"
    antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$warp_field,$output_mri] \
    --interpolation BSpline \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$target_mri,$source_mri,1] \
    --transform Rigid[0.1] \
    --metric MI[$target_mri,$source_mri,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    #echo "Done rigidly registering $source_mri to $target_mri"
}

function KUL_run_fmriprep {
    if [ ! -f fmriprep/sub-${participant}.html ]; then
        
        # preparing for fmriprep
        cp study_config/run_fmriprep.txt KUL_LOG/sub-${participant}_run_fmriprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        if [ $n_fMRI -gt 0 ]; then
            #fmriprep_options="--fs-no-reconall --use-aroma --use-syn-sdc "
            fmriprep_options="--fs-no-reconall "
        else
            fmriprep_options="--fs-no-reconall --anat-only "
        fi
        sed -i.bck "s/fmriprep_options: /fmriprep_options: ${fmriprep_options}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        
        # running fmriprep
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_fmriprep.txt 
        
        # cleaning the working directory
        rm -fr fmriprep_work_${participant}
        
    else
        echo "Fmriprep already done"
    fi

    # Export whenever fmriprep's output exists, not only when it just ran.
    #
    # These three used to live inside the branch above, so once
    # fmriprep/sub-*.html existed they were unreachable: deleting T1w_GM.nii.gz or
    # the report symlink meant a full fmriprep re-run to get a file that takes a
    # second to make. Everything here reads from fmriprep/ and writes into
    # RESULTS/ or REPORT/, so it is safe and cheap to redo on every run.
    local _fp_anat="fmriprep/sub-$participant/anat"
    if [ -f "$_fp_anat/sub-${participant}_desc-preproc_T1w.nii.gz" ]; then
        cp -f "$_fp_anat/sub-${participant}_desc-preproc_T1w.nii.gz" \
            $globalresultsdir/Anat/T1w_fmriprep.nii.gz
    fi
    # GM mask: only rebuilt when absent, since it is the one that costs anything
    # (two maskfilter passes) and nothing upstream of it changes between runs.
    if [ -f "$_fp_anat/sub-${participant}_dseg.nii.gz" ] && \
       [ ! -f "$globalresultsdir/Anat/T1w_GM.nii.gz" ]; then
        echo "Rebuilding Anat/T1w_GM.nii.gz from fmriprep's dseg"
        mrcalc "$_fp_anat/sub-${participant}_dseg.nii.gz" 1 -eq \
            "$_fp_anat/sub-${participant}_dseg.nii.gz" -mul - | \
            maskfilter - median - | \
            maskfilter - dilate $globalresultsdir/Anat/T1w_GM.nii.gz
    fi
    # A symlink, so re-point it rather than leaving a stale or missing one.
    if [ -f fmriprep/sub-${participant}.html ]; then
        ln -sfn ${cwd}/fmriprep/sub-${participant}.html \
            ${cwd}/REPORT/sub-${participant}_03_fmriprep.html
    fi
}

function KUL_run_dwiprep {
    if [ $n_dwi -gt 0 ];then
        if [ ! -f dwiprep/sub-${participant}/dwiprep_is_done.log ]; then
            cp study_config/${dwiprep_config_file} KUL_LOG/sub-${participant}_run_dwiprep.txt
            sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_dwiprep.txt
            rm -f KUL_LOG/sub-${participant}_run_dwiprep.txt.bck
            
            KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_dwiprep.txt 
            
        else
            echo "Dwiprep already done"
        fi
    fi
}

function KUL_run_freesurfer {
    if [ ! -f BIDS/derivatives/freesurfer/${participant}_freesurfer_is.done ]; then
        cp study_config/run_freesurfer.txt KUL_LOG/sub-${participant}_run_freesurfer.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_freesurfer.txt
        rm -f KUL_LOG/sub-${participant}_run_freesurfer.txt.bck
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_freesurfer.txt 
    else
        echo "Freesurfer already done"
    fi
}


function KUL_run_fastsurfer {
if [ ! -f KUL_LOG/sub-${participant}_FastSurfer.done ]; then
    # make your log file
    kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"

    fs_output="${cwd}/BIDS/derivatives/freesurfer"
    fasu_output="$derivativesdir/FastSurfer"
    if [ $vbg -eq 1 ];then
        T1_4_parc=$vbg_dir/output_VBG/sub-${participant}/sub-${participant}_T1_nat_4parc.mgz
    else
        T1_4_parc="${cwd}/$T1w"
    fi
    mkdir -p ${fs_output} >/dev/null 2>&1
    mkdir -p ${fasu_output} >/dev/null 2>&1

    # Run recon-all and convert the T1 to .mgz for display
    # running with -noskulltrip and using brain only inputs
    # for recon-all
    # if we can run up to skull strip, break, fix with hd-bet result then continue it would be much better
    # if we can switch to fast-surf, would be great also
    # another possiblity is using recon-all -skullstrip -clean-bm -gcut -subjid <subject name>
    
    task_in="recon-all -i ${T1_4_parc} -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -autorecon1 -no-isrunning"
    KUL_task_exec $verbose_level "FastSurfer part 1: recon-all stage 1" "$KUL_LOG_DIR/FastSurfer"

    FaSu_loc=$(which run_fastsurfer.sh)

    # --- FastSurfer / FreeSurfer version compatibility check ---
    fs_version_str=$(recon-all --version 2>/dev/null)
    fs_major=$(echo "$fs_version_str" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 | cut -d. -f1)
    fasu_version_str=$(run_fastsurfer.sh --version 2>/dev/null || echo "")
    fasu_major=$(echo "$fasu_version_str" | grep -oE '[0-9]+\.[0-9]+' | head -1 | cut -d. -f1)
    kul_echo "FreeSurfer version detected: ${fs_version_str:-unknown}"
    kul_echo "FastSurfer version detected: ${fasu_version_str:-unknown}"
    if [[ -n "$fs_major" && -n "$fasu_major" ]]; then
        if [[ "$fs_major" -ge 8 && "$fasu_major" -lt 3 ]]; then
            kul_echo "WARNING: FreeSurfer v${fs_major}.x is installed but FastSurfer v${fasu_major}.x was detected."
            kul_echo "WARNING: FastSurfer v3+ is required for FreeSurfer v8 compatibility."
            kul_echo "WARNING: Please update FastSurfer: https://github.com/Deep-MI/FastSurfer/releases"
            kul_echo "WARNING: Proceeding with --ignore_fs_version but results may be unreliable."
        fi
    elif [[ -z "$fasu_version_str" && -n "$fs_major" && "$fs_major" -ge 8 ]]; then
        kul_echo "WARNING: FreeSurfer v${fs_major}.x detected but FastSurfer version could not be determined."
        kul_echo "WARNING: Please ensure FastSurfer v3+ is installed for FreeSurfer v8 compatibility."
        kul_echo "WARNING: Get FastSurfer v3+: https://github.com/Deep-MI/FastSurfer/releases"
    fi
    # --- end version check ---

    user_id_str=$(id -u $(whoami))
    T1_4_FaSu=$(basename ${T1_4_parc})
    nvram=$(echo $(nvidia-smi --query-gpu=memory.free --format=csv) | rev | cut -d " " -f2 | rev)
    if [[ ! -z ${nvram} ]]; then
        if [[ ${nvram} -lt 6000 ]]; then
            batch_fasu="4"
        elif [[ ${nvram} -gt 6500 ]] && [[ ${nvram} -lt 7000 ]]; then
            batch_fasu="6"
        elif [[ ${nvram} -gt 7000 ]]; then
            batch_fasu="8"
        fi
    else
        batch_fasu="2"
    fi


    if [[ ! -z ${FaSu_loc} ]]; then

        if [ ${nvram} -lt 5500 ]; then
            FaSu_cpu=" --no_cuda "
            FaSu_mode="cpu"
            #echo " Running FastSurfer without CUDA " | tee -a ${prep_log}
        else
            FaSu_cpu=""
            FaSu_mode="cuda-gpu"
            #echo " Running FastSurfer with CUDA " | tee -a ${prep_log}
        fi

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        task_in="run_fastsurfer.sh --t1 ${T1_4_parc} \
        --sid sub-${participant} --sd ${fasu_output} --fsaparc --parallel --threads ${ncpu} \
        --fs_license $FS_LICENSE --py python ${FaSu_cpu} --ignore_fs_version --batch ${batch_fasu}"
        kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"
        KUL_task_exec $verbose_level "FastSurfer part 2: Fastsurfer itself (script & $FaSu_mode mode)" "$KUL_LOG_DIR/FastSurfer"

    else

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        echo "Local FastSurfer not found, switching to Docker version" | tee -a ${prep_log}
        T1_4_FaSu=$(basename ${T1_4_parc})
        dir_4_FaSu=$(dirname ${T1_4_parc})

        if [ ${nvram} -lt 5500 ]; then
            FaSu_v="cpu"
        else
            FaSu_v="gpu"
        fi

        task_in="docker run -v ${dir_4_FaSu}:/data -v ${fasu_output}:/output \
        -v $FREESURFER_HOME:/fs60 --rm --user ${user_id_str} fastsurfer:${FaSu_v} \
        --fs_license /fs60/$(basename $FS_LICENSE) --sid sub-${participant} \
        --sd /output/ --t1 /data/${T1_4_FaSu} \
        --parallel --threads ${ncpu}"
        kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"
        KUL_task_exec $verbose_level "FastSurfer part 2: Fastsurfer itself (docker & $FaSu_mode mode)" "$KUL_LOG_DIR/FastSurfer"

    fi

    #fs_output="${cwd}/BIDS/derivatives/freesurfer"
    #fasu_output="$derivativesdir/FastSurfer"
    

    # time to copy the surfaces and labels from FaSu to FS dirtask_exec
    # here we run FastSurfer first and 

    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/surf ${output_d}/${participant}_FS_output/${participant}/
    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/label ${output_d}/${participant}_FS_output/${participant}/
    
    #cp -rf ${output_d}/sub-${participant}/surf/* $fs_output/sub-${participant}/surf/
    #cp -rf ${output_d}/sub-${participant}/label/* $fs_output/sub-${participant}/surf/label/ 

    rsync -azv ${fasu_output}/sub-${participant}/ ${fs_output}/sub-${participant}/

    #task_in="recon-all -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -all -noskullstrip"
    #task_exec

    task_in="recon-all -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} \
        -parallel -no-isrunning -make all"
    KUL_task_exec $verbose_level "FastSurfer part 3: recon-all -make-all" "FastSurfer" || { kul_echo "FastSurfer part 3 failed — NOT writing FastSurfer.done"; return 1; }

    #fs_parc_mgz="${fs_output}/${participant}/mri/aparc+aseg.mgz"
    touch KUL_LOG/sub-${participant}_FastSurfer.done

else
    echo "Already done FastSurfer"
fi
}


function KUL_segment_tumor {
    
    # Segmentation of the tumor
    
    # check if it needs to be performed
    if [ $hdglio -eq 1 ];then

        if [ ! -f "$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" ]; then
            KUL_anat_segment_tumor.sh -p $participant -v $verbose_level
            lesion_png=RESULTS/sub-${participant}/Lesion/sub-${participant}_tumor_segment.png
            if [ -f $lesion_png ]; then
                cp -f $lesion_png REPORT/sub-${participant}_01_tumor_segment.png
            fi
            if [ $vbg -gt 0 ]; then
                if [ $vbg -eq 1 ]; then
                    vbg_lesion="$derivativesdir/KUL_anat_segment_tumor/sub-${participant}_lesion_and_cavity.nii.gz"
                elif [ $vbg -eq 2 ]; then
                    vbg_lesion="$derivativesdir/KUL_anat_segment_tumor/sub-${participant}_lesion_and_cavity.nii.gz"
                elif [ $vbg -eq 3 ]; then
                    vbg_lesion="${cwd}/RESULTS/sub-${participant}/Lesion/lesion.nii.gz"
                fi
                cp $vbg_lesion BIDS/sub-${participant}/anat/sub-${participant}_T1w_label-lesion_roi.nii.gz
            fi
        fi
    fi

}

function KUL_run_VBG {
    if [ $vbg -gt 0 ]; then
        if [ $vbg -eq 1 ]; then
            vbg_extra_axial=""
            vbg_lesion="$derivativesdir/KUL_anat_segment_tumor/sub-${participant}_lesion_and_cavity.nii.gz"
        elif [ $vbg -eq 2 ]; then
            vbg_extra_axial="-E"
            vbg_lesion="$derivativesdir/KUL_anat_segment_tumor/sub-${participant}_lesion_and_cavity.nii.gz"
        elif [ $vbg -eq 3 ]; then
            vbg_extra_axial=""
            vbg_lesion="${cwd}/RESULTS/sub-${participant}/Lesion/lesion.nii.gz"
        fi

        vbg_test="$vbg_dir/output_VBG/sub-${participant}/sub-${participant}_T1_nat_filled.nii.gz"
        if [[ ! -f $vbg_test ]]; then
            echo "Computing KUL_VBG"
            mkdir -p ${cwd}/BIDS/derivatives/freesurfer/sub-${participant}
            mkdir -p $vbg_dir
            
            # Use whatever FreeSurfer is already configured on PATH (exported
            # by the KUL_Linux_setup installer) instead of a separately hardcoded path
            # here, which can silently drift out of sync with the actual
            # install location (as it did: this used to point at a directory
            # that no longer exists).
            if [ -z "$FREESURFER_HOME" ]; then
                kul_echo "ERROR: FREESURFER_HOME is not set. Source the KUL environment block (see KUL_Linux_setup) before running this pipeline."
                return 1
            fi
            if [ ! -f "$FREESURFER_HOME/license.txt" ] && [ ! -f "$FREESURFER_HOME/.license" ]; then
                kul_echo "ERROR: no license.txt/.license found in FREESURFER_HOME ($FREESURFER_HOME)."
                return 1
            fi
            export SUBJECTS_DIR=$FREESURFER_HOME/subjects
            export FS_LICENSE=$FREESURFER_HOME/license.txt
            source $FREESURFER_HOME/SetUpFreeSurfer.sh
            fs_v=$(recon-all --version)
            kul_echo "Using $fs_v for VBG"

            task_in="KUL_VBG.sh -S ${participant} \
                -l $vbg_lesion \
                -o $vbg_dir \
                -m $vbg_dir \
                $vbg_extra_axial \
                -z T1 -b -B 1 -t -P 1 -M -O -H -n $ncpu"
            KUL_task_exec $verbose_level "KUL_VBG" "7_VBG" || { kul_echo "KUL_VBG failed — not copying possibly incomplete output to freesurfer derivatives"; return 1; }

            # copy the output of VBG to the derivatives freesurfer directory
            cp -r $vbg_dir/output_VBG/sub-${participant}_FS_output/sub-${participant} \
                BIDS/derivatives/freesurfer/

            # add to the report
            KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz \
                -t 2 -d REPORT -f 04_VBG_input
            KUL_mrview_figure.sh -p ${participant} -u $vbg_dir/output_VBG/sub-${participant}/sub-${participant}_T1_nat_filled.nii.gz \
                -t 2 -d REPORT -f 04_VBG_output
            #convert label:"Results of VBG:\nTop: Original T1w\nBottom: Filled T1w" REPORT/VBG_temp_caption.png
            montage REPORT/sub-${participant}_04_VBG_input.png REPORT/sub-${participant}_04_VBG_output.png -tile 1x2 -geometry +0+0 \
                -mode Concatenate REPORT/sub-${participant}_04_VBG.png   
            rm REPORT/sub-${participant}_04_VBG_input.png REPORT/sub-${participant}_04_VBG_output.png  

            echo "Done computing KUL_VBG"

        else
            echo "KUL_VBG has already run"
        fi

    fi
}


function KUL_run_multiparc {
    # Types 4, 5, 6: no VBG, but FWT still needs aparc+aseg and multiparc output.
    # Calls KUL_FS_multiparc.sh which handles recon-all (or FastSurfer with -X)
    # then adds Lausanne2018 x5, Glasser, thalamic, and brainstem parcellations.
    if [ $vbg -gt 0 ] || [ $fwt -eq 0 ]; then
        return 0
    fi

    if [ $multiparc -gt 0 ] && [ ! -f KUL_LOG/sub-${participant}_multiparc.done ]; then

        _fs_multiparc="${kul_main_dir}/KUL_FS_multiparc.sh"
        _fastsurfer_flag=""
        [ $use_fastsurfer -eq 1 ] && _fastsurfer_flag="-X"

        task_in="${_fs_multiparc} \
            -s sub-${participant} \
            -f ${cwd}/BIDS/derivatives/freesurfer \
            -i ${T1w[0]} \
            -n ${ncpu} ${_fastsurfer_flag}"
        KUL_task_exec $verbose_level "KUL_FS_multiparc (recon-all + parcellation)" "09_multiparc" || { kul_echo "KUL_FS_multiparc failed — NOT writing multiparc.done"; return 1; }

        touch KUL_LOG/sub-${participant}_multiparc.done

    else
        echo "multiparc already done"
    fi
}


function KUL_run_FWT {
    if [ $fwt -ne 1 ]; then
        echo "FWT not required for this analysis type"
        return 0
    fi
    # # Ensure updated KUL_FWT scripts (sibling repo) take priority over any older installation in PATH
    # _kul_nis_dir="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
    # export PATH="${_kul_nis_dir}/../KUL_FWT:$PATH"
    if [ $n_dwi -gt 0 ];then
        config="tracks_list.txt"
        if [ ! -f KUL_LOG/sub-${participant}_FWT.done ]; then

            # Resolve FS aparc+aseg path: prefer VBG FS output, fall back to standard freesurfer derivatives
            _vbg_fs_apas="$vbg_dir/output_VBG/sub-${participant}/sub-${participant}_FS_output/sub-${participant}/mri/aparc+aseg.mgz"
            _std_fs_apas="$cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz"
            if [ -f "$_vbg_fs_apas" ]; then
                _fs_apas="$_vbg_fs_apas"
            else
                _fs_apas="$_std_fs_apas"
            fi
            _fs_scale3="$(dirname $_fs_apas)/lausanne2018.scale3+aseg.mgz"

            task_in="KUL_FWT_make_VOIs.sh -p ${participant} \
            -F ${_fs_apas} \
            -c $cwd/study_config/${config} \
            -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
            -o $kulderivativesdir/sub-${participant}/FWT \
            -n $ncpu"
            KUL_task_exec $verbose_level "KUL_FWT voi generation" "12_FWTvoi" || { kul_echo "FWT VOI generation failed — NOT writing FWT.done"; return 1; }

            KUL_activate_conda_env ${scilpy}
            _fwt_rfa_opt=""
            [ $use_rfa_mod_fod -eq 1 ] && _fwt_rfa_opt="-U"
            _fwt_qq_opt=""
            [ $run_fwt_tractometry -eq 1 ] && _fwt_qq_opt="-Q"
            task_in="KUL_FWT_make_TCKs.sh -p ${participant} \
            -F ${_fs_apas} \
            -c $cwd/study_config/${config} \
            -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
            -o $kulderivativesdir/sub-${participant}/FWT \
            -T 1 -a iFOD2 \
            -f 1 \
            -S ${_fwt_rfa_opt} ${_fwt_qq_opt} \
            -n $ncpu"
            KUL_task_exec $verbose_level "KUL_FWT tract generation" "12_FWTtck" || { kul_echo "FWT tract generation failed — NOT writing FWT.done"; return 1; }

            touch KUL_LOG/sub-${participant}_FWT.done

        else
            echo "FWT already done"
        fi

        # Always sync FWT output to RESULTS/Tracto so re-runs and new bundles are picked up.
        # Tract maps are resampled to T1w resolution for direct comparison with anatomy/lesion masks.
        mkdir -p $globalresultsdir/Tracto
        rm -fr $globalresultsdir/Tracto/*
        _t1w_ref="$globalresultsdir/Anat/T1w.nii.gz"
        [[ ! -f "$_t1w_ref" ]] && _t1w_ref="fmriprep/sub-${participant}/anat/sub-${participant}_desc-preproc_T1w.nii.gz"
        for tck_outdir in "$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output"/*_output; do
            tract_name=$(basename "$tck_outdir" _output)
            fin_tck="${tck_outdir}/${tract_name}_fin_BT_iFOD2.tck"
            fin_map="${tck_outdir}/${tract_name}_fin_map_BT_iFOD2.nii.gz"
            _out_map="$globalresultsdir/Tracto/Tract-csd_${tract_name}.nii.gz"
            if [ -f "$fin_tck" ]; then
                cp "$fin_tck" "$globalresultsdir/Tracto/Tract-csd_${tract_name}.tck"
                [ -f "$fin_map" ] && mrgrid "$fin_map" regrid -template "$_t1w_ref" -interp linear "$_out_map" -force -quiet
            else
                use_tck=$(ls "${tck_outdir}/${tract_name}_filt"*"_BT_iFOD2.tck" 2>/dev/null | grep -v "_inMNI" | sort -V | tail -1)
                use_map=$(ls "${tck_outdir}/${tract_name}_filt"*"_map_BT_iFOD2.nii.gz" "${tck_outdir}/${tract_name}_filt"*"_BT_iFOD2_map.nii.gz" 2>/dev/null | grep -v "_inMNI" | sort -V | tail -1)
                if [ -n "$use_tck" ]; then
                    echo "  ${tract_name}: fin not found, using $(basename $use_tck)"
                    cp "$use_tck" "$globalresultsdir/Tracto/Tract-csd_${tract_name}.tck"
                    [ -n "$use_map" ] && mrgrid "$use_map" regrid -template "$_t1w_ref" -interp linear "$_out_map" -force -quiet
                fi
            fi
        done
        pdfunite $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*_output/Screenshots/*fin_BT_iFOD2_inMNI_screenshot2_niGB.pdf $globalresultsdir/Tracto/Tracts_Summary.pdf 2>/dev/null || true
        cp -f $globalresultsdir/Tracto/Tracts_Summary.pdf REPORT/sub-${participant}_06_Tract_Summary.pdf 2>/dev/null || true

        # Surface the -Q tractometry output, which otherwise never leaves KUL_compute.
        # KUL_FWT writes the per-bundle spider report into <bundle>_output/QQ/, and
        # nothing above copies it out: the loop takes only .tck and _fin_map, and
        # RESULTS/Tracto is wiped at the top of this block, so hand-placed copies
        # would not survive a re-run either.
        #
        # The .html alone, deliberately: it already renders every along-tract
        # profile that the sibling *_scores_*_plot.pdf files show, the spider PDF's
        # content as its interactive tower, and (since the connectivity table was
        # added to it) the endpoint parcel pairs that the connectivity PNG/CSV
        # encoded. It is fully self-contained -- inline SVG and JS, no external
        # assets -- so it survives being copied away from its QQ directory.
        #
        # These appear only once the whole subject's FWT run is done
        # (KUL_FWT_bundle_spider_plot.py runs after every bundle, since it
        # normalizes each metric across bundles), so a partial run copies nothing
        # and a later re-run picks them all up.
        _qq_src="$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output"
        _qq_dst="REPORT/sub-${participant}_06_Tract_QQ"
        # The screenshot contact sheet (KUL_FWT_bundle_report.py) links to each
        # bundle's spider page by bare filename, so both have to land in the same
        # directory for "metrics ->" to resolve. Hence the flat copy here rather
        # than per-bundle subdirectories.
        if compgen -G "${_qq_src}/*_output/QQ/*spider3d*.html" > /dev/null 2>&1 || \
           compgen -G "${_qq_src}/sub-${participant}*_FWT_report.html" > /dev/null 2>&1; then
            mkdir -p "$_qq_dst"
            cp -f "${_qq_src}"/*_output/QQ/*spider3d*.html "$_qq_dst/" 2>/dev/null || true
            cp -f "${_qq_src}"/sub-${participant}*_FWT_report.html "$_qq_dst/" 2>/dev/null || true
            echo "  tractometry (QQ) reports copied to $_qq_dst"
        fi

        # Collect the .trk copies KUL_FWT_make_TCKs.sh writes next to each .tck
        # (for freeview, which does not read MRtrix .tck). Conversion lives in FWT
        # rather than here so it happens where the tracking reference is known --
        # see the .trk comment there for why that reference is subj_FA and not the
        # FS parcellation. Own directory, so Tracto stays .tck-only.
        if compgen -G "$_qq_src/*_output/*_fin_*.trk" > /dev/null 2>&1; then
            mkdir -p $globalresultsdir/TRK
            rm -f $globalresultsdir/TRK/*.trk
            _n_trk=0
            for _src_trk in "$_qq_src"/*_output/*_fin_*.trk; do
                _trk_name=$(basename "$(dirname "$_src_trk")" _output)
                if cp -f "$_src_trk" "$globalresultsdir/TRK/Tract-csd_${_trk_name}.trk"; then
                    _n_trk=$((_n_trk + 1))
                fi
            done
            echo "  copied $_n_trk .trk bundles to $globalresultsdir/TRK"
        fi
    fi
}

function KUL_anatomical_biascorrect {
    check="KUL_LOG/sub-${participant}_anat_biascorrect.done"
    if [ ! -f $check ]; then

        KUL_anat_biascorrect.sh -p $participant -v $verbose_level

        touch $check
    else 
        echo "Anatomical bias correction already done"
    fi
}

function KUL_register_anatomical_images {
    check="KUL_LOG/sub-${participant}_anat_reg.done"
    if [ ! -f $check ]; then

        export KUL_MP2RAGE_FILE="$MP2RAGE"
        export KUL_SWIP_FILE="$SWIp"
        KUL_anat_register.sh -p $participant -c -v $verbose_level
        unset KUL_MP2RAGE_FILE KUL_SWIP_FILE
        cp  $cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_anat_register_rigid/*reg2_T1w.nii.gz $globalresultsdir/Anat/
        cp  $cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_anat_register_rigid/T1w.nii.gz $globalresultsdir/Anat/
        touch $check
    
    else 
        echo "Anatomical registration already done"
    fi
}

function KUL_clear_cT1w {
    
    # a funtion to remove the cT1w (gadolinium enhanced T1w) away since it conflicts during msbp
    clear_cT1w_outputdir="$kulderivativesdir/sub-${participant}/cT1w"
    mkdir -p $clear_cT1w_outputdir

    if [ $ncT1w -gt 0 ]; then
        source_mri="${cT1w%*.nii.gz}*"
        if [ -f $cT1w ]; then
            mv $source_mri $clear_cT1w_outputdir 
        fi
    fi

}

function KUL_fmriproc {

    if [ $n_fMRI -gt 0 ];then

        if [ ! -f ${cwd}/KUL_LOG/sub-${participant}_SPM.done ]; then
            if [ "$fmri_engine" == "nilearn" ]; then
                _pyfmri_C_opt=""
                [ -n "$pyfmri_env_override" ] && _pyfmri_C_opt="-C $pyfmri_env_override"
                task_in="KUL_fmriproc_nilearn_new.sh -p $participant -S $smooth_fwhm -P $pfwe -c $ncpu $_pyfmri_C_opt"
                KUL_task_exec $verbose_level "KUL_fmriproc_nilearn_new" "7_fmriproc_nilearn" || kul_echo "KUL_fmriproc_nilearn_new failed for sub-${participant} — SPM.done will not be created (check 7_fmriproc_nilearn.error.log)"
            else
                task_in="KUL_fmriproc_spm_new.sh -p $participant -S $smooth_fwhm -P $pfwe -c $ncpu"
                KUL_task_exec $verbose_level "KUL_fmriproc_spm_new" "7_fmriproc_spm" || kul_echo "KUL_fmriproc_spm_new failed for sub-${participant} — SPM.done will not be created (check 7_fmriproc_spm.error.log)"
            fi

        fi

        if [ ! -f ${cwd}/KUL_LOG/sub-${participant}_melodic.done ]; then
            task_in="KUL_fmriproc_conn.sh -p $participant"
            KUL_task_exec $verbose_level "KUL_fmriproc_conn" "8_fmriproc_conn" || kul_echo "KUL_fmriproc_conn failed for sub-${participant} — melodic.done will not be created (check 8_fmriproc_conn.error.log)"
        fi

        # ── fMRI figures and reports for REPORT/ ────────────────────────────
        #
        # Deliberately outside the two branches above. These used to be produced
        # only in the run that computed the maps, so once the .done marker
        # existed they could never come back: the maps were still there, the
        # markers said "done", and nothing redrew the pictures. Deleting REPORT/
        # lost every fMRI visual permanently.
        #
        # Each item is produced only when it is missing, so a healthy tree
        # renders nothing and this costs nothing per run.

        # GLM activations: the hardwired wc (with-confounds) Bizzi-thresholded
        # maps only (p<0.001 unc, k>=50) — matches the RESULTS/.../SPM selection.
        for bizzi_map in $derivativesdir/SPM/*_wc/spmT_0001_p001unc_k50.nii; do
            [ -f "$bizzi_map" ] || continue
            task=$(basename $(dirname $bizzi_map))
            task=${task%_wc}
            [ -f "REPORT/sub-${participant}_05_afMRI_${task}_p001unc_k50.png" ] && continue
            KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz \
                -o "$bizzi_map" -t 2 -d REPORT -f 05_afMRI_${task}_p001unc_k50
        done

        # Melodic networks. The threshold is baked into the filename, so match on
        # the stem rather than trying to predict it.
        for spm in RESULTS/sub-${participant}/Melodic/*.nii; do
            [ -f "$spm" ] || continue
            task=$(basename $spm)
            [ -n "$(ls REPORT/sub-${participant}_05_rsfMRI_${task}_Thr_*.png 2>/dev/null)" ] && continue
            max_T=$(mrstats -output max $spm)
            thresh=$(awk "BEGIN {print $max_T/3}")
            mrcalc $spm $thresh -gt REPORT/spm_tmp_$task -force -quiet
            KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz -o REPORT/spm_tmp_$task \
                -t 2 -d REPORT -f 05_rsfMRI_${task}_Thr_${thresh}
            rm -f REPORT/spm_tmp_$task
        done

        # melodic's own per-run HTML report, which otherwise never leaves the
        # derivative. Symlinked, so it stays in step with the decomposition and
        # costs nothing.
        for _mrep in $derivativesdir/FSL_melodic/stats_*/report/00index.html; do
            [ -f "$_mrep" ] || continue
            _mtask=$(basename $(dirname $(dirname "$_mrep")))
            ln -sfn "$_mrep" "REPORT/sub-${participant}_05_melodic_${_mtask#stats_}.html"
        done

        # The rsfMRI-networks step (-N) writes its report into its own derivative
        # and into RESULTS/rsfMRI_Networks, but never into REPORT -- so the one
        # place a clinician looks did not have it.
        _rsn_rep="$kulderivativesdir/rsfMRI_networks/analysis/reports/sub-${participant}_rsfmri_networks_report"
        for _ext in pdf html; do
            [ -f "${_rsn_rep}.${_ext}" ] || continue
            cp -f "${_rsn_rep}.${_ext}" \
                "REPORT/sub-${participant}_05_rsfMRI_networks_report.${_ext}" 2>/dev/null
        done
    fi

}

function KUL_run_rsfMRI_networks {
    if [ $rsfmri_networks -eq 1 ]; then
        rsfmri_check=${cwd}/KUL_LOG/sub-${participant}_rsfMRI_networks.done
        if [ ! -f $rsfmri_check ]; then
            _pyfmri_c_opt=""
            [ -n "$pyfmri_env_override" ] && _pyfmri_c_opt="-c $pyfmri_env_override"
            task_in="KUL_run_rsfMRI_networks.sh -p $participant $_pyfmri_c_opt -P $rsfmri_profile -v $verbose_level"
            KUL_task_exec $verbose_level "KUL_run_rsfMRI_networks" "9_rsfmri_networks" || kul_echo "KUL_run_rsfMRI_networks failed for sub-${participant} — rsfMRI_networks.done will not be created (check 9_rsfmri_networks.error.log)"
        else
            echo "rsfMRI networks already done"
        fi
    fi
}

function KUL_run_dsc {

    if [ $skip_dsc -eq 1 ]; then
        echo "DSC perfusion processing skipped (-W)"
        return 0
    fi
    if [ ${n_dsc:-0} -eq 0 ]; then
        return 0
    fi

    if [ -f KUL_LOG/sub-${participant}_DSC.done ]; then
        echo "DSC perfusion already done"
        return 0
    fi

    # Runs after VBG/multiparc on purpose: the contralesional NAWM reference
    # needs a FreeSurfer aseg, which only exists once one of those has run.
    # KUL_dsc_perfusion.sh still writes the parametric maps if the aseg or the
    # lesion mask is missing, and can be re-run later to add the ratios.
    _dsc_anat="$globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz"
    [ -f "$_dsc_anat" ] || _dsc_anat="$globalresultsdir/Anat/T1w.nii.gz"

    _dsc_y_opt=""
    [ -n "$pyfmri_env_override" ] && _dsc_y_opt="-y $pyfmri_env_override"

    task_in="${kul_main_dir}/KUL_dsc_perfusion.sh -p ${participant} \
        -a ${_dsc_anat} ${_dsc_y_opt} \
        -n ${ncpu} -v ${verbose_level}"
    KUL_task_exec $verbose_level "KUL_dsc_perfusion (DSC perfusion maps)" "13_DSC" || \
        { kul_echo "KUL_dsc_perfusion failed - NOT writing DSC.done"; return 1; }

    dsc_png=REPORT/sub-${participant}_06_DSC_rCBV.png
    if [ -f $dsc_png ]; then
        echo "  DSC QC figure: $dsc_png"
    fi

    touch KUL_LOG/sub-${participant}_DSC.done

}


function KUL_run_dwiprep_anat {

    dwi_anat_check=${cwd}/KUL_LOG/sub-${participant}_dwiprep_anat.done
    if [ ! -f $dwi_anat_check ]; then
        task_in="KUL_dwiprep_anat.sh -p $participant -n $ncpu"
        KUL_task_exec $verbose_level "KUL_dwiprep_anat" "11_dwiprep_anat" || { kul_echo "KUL_dwiprep_anat failed — NOT writing dwiprep_anat.done"; return 1; }

        touch $dwi_anat_check
    fi

    # Report copies live outside the branch above, so deleting REPORT/ does not
    # require re-running dwiprep_anat to get three files back that are pure
    # copies of QA output sitting in dwiprep/.
    local _qa="dwiprep/sub-${participant}/sub-${participant}/qa"
    if [ -f "$_qa/sub-${participant}_T1w_with_fa.png" ]; then
        cp -f "$_qa/sub-${participant}_T1w_with_fa.png" \
            REPORT/sub-${participant}_02_T1w_with_fa.png
        cp -f "$_qa/sub-${participant}_T1w_brain_with_fa.png" \
            REPORT/sub-${participant}_02_T1w_brain_with_fa.png 2>/dev/null
    fi
    if [ -f "dwiprep/sub-${participant}/sub-${participant}/eddy_qc/quad/qc.pdf" ]; then
        cp -f dwiprep/sub-${participant}/sub-${participant}/eddy_qc/quad/qc.pdf \
            REPORT/sub-${participant}_02_eddy_qc.pdf
    fi

}

function KUL_run_dwiprep_MNI {

    dwi_MNI_check=${cwd}/KUL_LOG/sub-${participant}_dwiprep_MNI.done
    if [ ! -f $dwi_MNI_check ]; then
        task_in="KUL_dwiprep_MNI.sh -p $participant -n $ncpu"
        KUL_task_exec $verbose_level "KUL_dwiprep_MNI" "12_dwiprep_MNI" || { kul_echo "KUL_dwiprep_MNI failed — NOT writing dwiprep_MNI.done"; return 1; }

        touch $dwi_MNI_check
    fi

}

function KUL_calc_DTI_ALPS {

    if [ $alps -eq 1 ]; then
        dti_ALPS_check=${cwd}/KUL_LOG/sub-${participant}_dti_ALPS.done
        if [ ! -f $dti_ALPS_check ]; then
            task_in="${kul_main_dir}/KUL_DTI_ALPS/KUL_calc_DTIALPS.sh -p $participant -n $ncpu"
            KUL_task_exec $verbose_level "KUL_calc_DTIALPS" "13_dti_ALPS" || { kul_echo "KUL_calc_DTIALPS failed — NOT writing dti_ALPS.done"; return 1; }

            touch $dti_ALPS_check
        fi
    fi

}

function KUL_run_cT1w_subtraction {

    cT1w_subt_check=${cwd}/KUL_LOG/sub-${participant}_cT1w_subtraction.done
    if [ ! -f $cT1w_subt_check ]; then
        if [ $nT1w -gt 0 ] && [ $ncT1w -gt 0 ]; then
            echo "Creating T1w subtraction image"
            mkdir -p $derivativesdir/cT1w_subtraction
            # get the mask
            cp -f fmriprep/sub-$participant/anat/sub-${participant}_desc-brain_mask.nii.gz $derivativesdir/cT1w_subtraction/brain_mask.nii.gz
            # get the lesion & exclude the tumor if any
            #echo $hdglio
            if [ $hdglio -eq 1 ];then
                if [  -f "$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" ]; then
                    mrcalc $derivativesdir/cT1w_subtraction/brain_mask.nii.gz $globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz \
                        -sub $derivativesdir/cT1w_subtraction/brain_final_mask.nii.gz -force
                    tomatch=$derivativesdir/cT1w_subtraction/brain_final_mask.nii.gz
                else
                    tomatch=$derivativesdir/cT1w_subtraction/brain_mask.nii.gz
                fi
            else
                tomatch=$derivativesdir/cT1w_subtraction/brain_mask.nii.gz
            fi
            # historgam match the cT1w and T1w
            mrhistmatch linear \
                -mask_input $tomatch \
                -mask_target $tomatch \
                $globalresultsdir/Anat/T1w.nii.gz $globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz \
                $globalresultsdir/Anat/T1w_matched2_cT1w.nii.gz -force
            # subtract them    
            mrcalc $globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz $globalresultsdir/Anat/T1w_matched2_cT1w.nii.gz -sub \
                $globalresultsdir/Anat/cT1w_T1w_subtracted.nii.gz -force
            touch $cT1w_subt_check
        fi
    else
        echo "T1w subtraction already done"    
    fi

}


# --- MAIN ---

# STEP 1 - BIDS conversion
KUL_convert2bids

KUL_check_participant

kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
mkdir -p $kulderivativesdir
mkdir -p $globalresultsdir/Anat
mkdir -p $globalresultsdir/SPM
mkdir -p $globalresultsdir/SPM_all
mkdir -p $globalresultsdir/Melodic
mkdir -p $globalresultsdir/Tracto
mkdir -p $globalresultsdir/PACS/fMRI
mkdir -p $cwd/RESULTS
mkdir -p $cwd/REPORT
# Donor-DICOM drop points, created up front so -R needs no manual mkdir. Both are
# searched (Karawun first, then RESULTS) and a single file in either is enough,
# so the user drops one donor in one place rather than making a folder and
# guessing which one matters.
mkdir -p $cwd/Karawun/sub-${participant}/DICOM
mkdir -p $globalresultsdir/DICOM

if [ $KUL_DEBUG -gt 0 ]; then 
    echo "kulderivativesdir: $kulderivativesdir"
    echo "globalresultsdir: $globalresultsdir"
fi

# Run BIDS validation
check_in=${KUL_LOG_DIR}/1_bidscheck.done
if [ ! -f $check_in ]; then

    docker run -ti --rm -v ${cwd}/BIDS:/data:ro bids/validator /data

    read -p "Are you happy? (y/n) " answ
    if [[ ! "$answ" == "y" ]]; then
        exit 1
    else
        touch $check_in
    fi
fi


# Check if fMRI and/or dwi data are present and/or to redo some processing
echo "Starting KUL_clinical_fmridti"
KUL_check_data

# Only meaningful once KUL_check_data has counted the fMRI runs, and pointless
# for an export-only run, which touches no python. Fails the run rather than
# warning: with fMRI data present, every fMRI step would otherwise get most of
# the way through the pipeline and then die one at a time.
if [ $export_only -eq 0 ] && [ ${n_fMRI:-0} -gt 0 ]; then
    KUL_check_pyfmri_env || exit 2
fi

KUL_verify_results
KUL_check_redo

# STEP 1 - run bias correction
KUL_anatomical_biascorrect

# STEP 2 - register all anatomical other data to the T1w without contrast
KUL_register_anatomical_images

# STEP 3 - run tumor segmentation 
KUL_segment_tumor

# STEP 3b - make the lesion available alongside the other T1w-space volumes
KUL_copy_lesion_to_anat
    
# STEP 4 - run fmriprep and continue
KUL_run_fmriprep &

# STEP 5 - run dwiprep and continue
KUL_run_dwiprep &
wait

# STEP 6 - generate Gd contrast T1w subtraction
KUL_run_cT1w_subtraction

# STEP 7 - get rid of the Gadolinium T1w image since it may conflict with msbp
KUL_clear_cT1w

# STEP 8 - run SPM & melodic
KUL_fmriproc

# STEP 9 - run VBG
KUL_run_VBG
wait

# STEP 9b - FastSurfer + KUL_multiparc (types 4, 5, 6 only — no VBG)
# this should only run if VBG is not used!!
KUL_run_multiparc
wait

# STEP 9b-bis - run rsfMRI network analysis (opt-in, -N)
# After VBG/multiparc, not before: the pipeline's step0 warps
# lausanne2018.scale3+aseg.mgz into its analysis space, and the seeds in the
# Presurgical_Somatotopic profile (Lip/Hand/Foot L/R) are defined on it. That
# file is written by KUL_VBG.sh -M on types 1/2/3 and by KUL_FS_multiparc.sh on
# types 4/5/6 -- both of which used to run *after* this step, so step0 logged
# "lausanne2018.scale3+aseg.mgz not found ... run KUL_FS_multiparc.sh first"
# and silently dropped those seeds on every type.
KUL_run_rsfMRI_networks

# STEP 9c - DSC perfusion (runs whenever perf data exists, unless -W)
# after VBG/multiparc, so the FreeSurfer aseg the NAWM reference needs exists
KUL_run_dsc

# STEP 10 - run SPM/melodic/msbp
# this is not okay!
# KUL_run_msbp
# wait

# STEP 10 run dwiprep_anat
KUL_run_dwiprep_anat

# STEP 11 run dwiprep_MNI (needed for DTI-ALPS; skipped for other types)
KUL_run_dwiprep_MNI

# STEP 12a run DTI_ALPS calculation (type 7 only)
KUL_calc_DTI_ALPS
wait

# STEP 12b - run Fun With Tracts
KUL_run_FWT

# STEP 15 - Prepare the Karawun folder, and the PACS drop folders beside it.
#
# Runs at the end of every pipeline run, NOT behind -R. It used to be gated on
# `make_dcm -eq 1`, and -R is the only flag that sets make_dcm -- but -R also
# sets results>0, which enters the export block near the top of this file and
# exits there, ~1200 lines before this point. So this step was unreachable on
# every realistic invocation, and the per-task fMRI activation labels that the
# -R block writes (only if Karawun/sub-*/T1w.nii.gz already exists) were
# therefore never produced: the Brainlab export shipped tracts and no fMRI.
#
# Ordering is the reason it lives exactly here: it needs the FWT output from
# STEP 12b above, and the -R block needs the T1w.nii.gz it writes. It depends on
# nothing that -R produces, so it belongs in the main flow.
#
# The importTractography command is printed at the end — run it manually, after
# copying a donor DICOM into Karawun/sub-${participant}/DICOM/ (or
# RESULTS/sub-${participant}/DICOM/) and reviewing the labels.
_fwt_tck_root="${derivativesdir}/FWT/sub-${participant}_TCKs_output"
karawun_prepare_check=${cwd}/KUL_LOG/sub-${participant}_karawun_prepare.done
if [ ! -d "$_fwt_tck_root" ]; then
    # Now that this is unconditional, a subject whose FWT produced nothing would
    # otherwise fail prep at the end of every single run. Skip cleanly and leave
    # no marker, so it retries once FWT has actually run.
    echo "No FWT tract output at $_fwt_tck_root"
    echo "  Skipping Karawun prep - it will retry on the next run."
elif [ ! -f $karawun_prepare_check ]; then
    _karawun_rc=0
    _karawun_type=""
    if [ $type -lt 5 ]; then
        _karawun_type=1; _karawun_thr=3;  _karawun_what="Preparing Karawun folder"
    elif [ $type -eq 5 ]; then
        _karawun_type=2; _karawun_thr=10; _karawun_what="Preparing Karawun folder (DBS ET)"
    elif [ $type -eq 6 ]; then
        _karawun_type=3; _karawun_thr=10; _karawun_what="Preparing Karawun folder (DBS Parkinson)"
    fi
    if [ -z "$_karawun_type" ]; then
        # type 7 (DTI-ALPS), and anything added later. This used to fall through
        # every branch in silence: _karawun_rc stayed 0, the marker was touched,
        # and the run reported a successful prep having done nothing at all.
        echo "No Karawun mapping for type ${type} - skipping Karawun prep."
        echo "  Mapped: types 1-4 (tumor), 5 (DBS ET), 6 (DBS Parkinson)."
        echo "  Run 'KUL_karawun_prepare.sh -p ${participant} -t <1|2|3>' by hand if you want it."
    else
        # RESULTS/Anat/T1w.nii.gz is the volume prep rescales into
        # Karawun/sub-*/T1w.nii.gz, and every label is then regridded onto that.
        # Without it prep still exits 0: it prints "T1w voxel sizes  -" from empty
        # mrinfo output, copies the .tck files, and produces no T1w and no labels
        # at all. Check the input up front rather than discovering it afterwards.
        if [ ! -f "$globalresultsdir/Anat/T1w.nii.gz" ]; then
            echo "ERROR: cannot prepare Karawun - ${globalresultsdir#$cwd/}/Anat/T1w.nii.gz is missing." >&2
            echo "       Every Karawun label is regridded onto it, so prep would produce" >&2
            echo "       a folder with tracts but no T1w and no labels." >&2
            echo "       Restore it (it is a copy of" >&2
            echo "       BIDS/derivatives/KUL_compute/sub-${participant}/KUL_anat_register_rigid/T1w.nii.gz)" >&2
            echo "       and re-run; no marker is written, so this will retry." >&2
            _karawun_rc=1
        else
            kul_echo "$_karawun_what"
            KUL_karawun_prepare.sh -p ${participant} -t $_karawun_type -r $_karawun_thr || _karawun_rc=$?
        fi

        # Only claim success if it actually produced something. The exit code
        # alone is not enough: KUL_karawun_prepare.sh returns 0 even when its
        # inputs were missing and it wrote no T1w and no labels, so trusting it
        # wrote the .done marker over a broken folder -- which then made every
        # later run print "already prepared" and skip it. Verify the two outputs
        # everything downstream depends on instead.
        if [ $_karawun_rc -eq 0 ] && \
           [ -f "$cwd/Karawun/sub-${participant}/T1w.nii.gz" ] && \
           [ -n "$(ls -A "$cwd/Karawun/sub-${participant}/labels" 2>/dev/null)" ]; then
            touch $karawun_prepare_check
        elif [ $_karawun_rc -eq 0 ]; then
            echo "ERROR: Karawun prep reported success but produced no T1w.nii.gz and/or no" >&2
            echo "       labels in Karawun/sub-${participant}/. NOT writing the .done marker," >&2
            echo "       so this retries rather than silently skipping from now on." >&2
        else
            echo "ERROR: Karawun prep failed (exit $_karawun_rc) — NOT writing $karawun_prepare_check" >&2
            echo "       fix the cause and re-run; it will retry rather than skip." >&2
        fi
    fi
else
    # Deliberate: the Karawun folder (tracts, T1w, labels) does not depend on the
    # -R underlay, so repeating -R with a different underlay regenerates the PACS
    # DICOMs but correctly leaves Karawun alone. Delete the marker to force a
    # rebuild.
    echo "Karawun folder already prepared (delete KUL_LOG/sub-${participant}_karawun_prepare.done to rebuild)"
fi

# The PACS drop folders, so the operator has somewhere to copy the maps they
# want exported BEFORE the first -R run. Same function the -R block calls; these
# used to be created only inside that block, i.e. only once it was already too
# late to put anything in them.
KUL_make_pacs_dropdirs


# STEP 17 - figure generation is no longer automatic.
# Review results first, then run with -F <underlay> to generate screenshots, or -R <underlay> to also push DICOMs.
if [ $figs -eq 1 ]; then
fig_check=${cwd}/KUL_LOG/sub-${participant}_figures.done
if [ ! -f $fig_check ]; then
    if [ $ncT1w -gt 0 ] || [ $ncT1w -eq -1 ]; then
        KUL_clinical_fmridti.sh -p $participant -t $type -F 1 -O "$orientations"
    fi
    if [ $nFLAIR -gt 0 ]; then
        KUL_clinical_fmridti.sh -p $participant -t $type -F 2 -O "$orientations"
    fi
    if [ $nSWI -gt 0 ]; then
        KUL_clinical_fmridti.sh -p $participant -t $type -F 3 -O "$orientations"
    fi
    if [ $nT1w -gt 0 ] && [ $ncT1w -lt 1 ] && [ $nFLAIR -eq 0 ]; then
        KUL_clinical_fmridti.sh -p $participant -t $type -F 4 -O "$orientations"
    fi
    if [ $nFGATIR -gt 0 ]; then
        KUL_clinical_fmridti.sh -p $participant -t $type -F 5 -O "$orientations"
    fi
    if [ $nDIR -gt 0 ]; then
        KUL_clinical_fmridti.sh -p $participant -t $type -F 6 -O "$orientations"
    fi
    if [ $nMP2RAGE -gt 0 ]; then
        KUL_clinical_fmridti.sh -p $participant -t $type -F 7 -O "$orientations"
    fi
    touch $fig_check
else
    echo "Figures already done"
fi
else
    echo "Tractography figure generation not required"
fi



# STEP 18 - DICOM generation is no longer automatic.
# Review figures in RESULTS/sub-${participant}/, then run:
#   KUL_clinical_fmridti.sh -p ${participant} -t ${type} -R <1|2|3|4> [-O orientations]
# to generate PACS DICOMs.

if [ $results -eq 0 ]; then
    echo ""
    echo "================================================================"
    echo " Finished processing sub-${participant}."
    echo ""
    echo " Nothing has been sent to PACS/Karawun yet - that step needs an"
    echo " explicit -R run on purpose, so you look at the output first:"
    echo ""
    echo "   1. Review the output maps and tractograms in"
    echo "        $globalresultsdir/"
    echo "          Anat/  SPM/  Melodic/  Tracto/ (.tck)  TRK/ (.trk, freeview)"
    echo "        REPORT/"
    echo "          sub-${participant}_06_Tract_QQ/    one HTML per bundle, plus"
    echo "          sub-${participant}_FWT_report.html  all bundle screenshots in one page"
    echo "      The fMRI/tract thresholds used for PACS export are only"
    echo "      computed when you run -F or -R below (auto = max/3 per"
    echo "      map, or interactively/-T if you override) - this review"
    echo "      is your chance to catch a threshold that looks wrong"
    echo "      before it is baked into the exported DICOMs."
    echo ""
    if [ -f "Karawun/sub-${participant}/T1w.nii.gz" ]; then
        echo "      The Karawun folder for Brainlab has already been prepared:"
        echo "        Karawun/sub-${participant}/   T1w, FAT1w, tck/, labels/"
        echo "      Review it too. It is gated by"
        echo "      KUL_LOG/sub-${participant}_karawun_prepare.done - delete that"
        echo "      marker to have the next run rebuild it."
    else
        echo "      NOTE: the Karawun/Brainlab folder was NOT prepared (see the"
        echo "      Karawun message earlier in this run). Without it, -R cannot"
        echo "      write the fMRI activation labels for Brainlab."
    fi
    echo ""
    echo "   2. Generate screenshots only, no PACS/Karawun push yet:"
    echo "        KUL_clinical_fmridti.sh -p ${participant} -t ${type} -F <1-7> [-O orientations]"
    echo ""
    echo "   3. Choose what goes to PACS by copying maps into:"
    echo "        $globalresultsdir/PACS_input/"
    echo "          overlays/              colour overlay fused on the anatomical"
    echo "          series_quantitative/   measurable series, ROI reads true units"
    echo "      Both folders already exist and carry a README.txt. Leave them"
    echo "      empty to get the automatic discovery instead - but note that"
    echo "      auto-discovery only picks up SPM/, Melodic/, the lesion, and"
    echo "      *only* nrCBV_corrected/nrCBF from Perfusion/. Anything else"
    echo "      (a raw rCBV_corrected, an ADC map, anything from outside this"
    echo "      pipeline) has to be copied in here to reach PACS."
    echo ""
    echo "   4. Drop ONE donor DICOM into either of these:"
    echo "        Karawun/sub-${participant}/DICOM/     <- searched first"
    echo "        $globalresultsdir/DICOM/"
    echo "      One file is enough - a single slice from a high-resolution"
    echo "      anatomical series (T1w, FLAIR, T2, ...). The same donor is"
    echo "      used for both PACS and Karawun/Brainlab, so one is all you"
    echo "      need; Karawun/ is the better place since the import command"
    echo "      in step 6 reads from there."
    echo ""
    echo "   5. Once happy with the review and the donor DICOM is in place:"
    echo "        KUL_clinical_fmridti.sh -p ${participant} -t ${type} -R <1-7> [-O orientations]"
    echo "      This writes the PACS DICOMs, and also adds the fMRI activation"
    echo "      labels to Karawun/sub-${participant}/labels/ for Brainlab."
    echo "      Re-running -R with a different underlay (e.g. -R 4 then -R 2)"
    echo "      regenerates the PACS DICOMs for that underlay; the Karawun"
    echo "      folder does not depend on the underlay and is left alone."
    echo ""
    echo "   6. -R does not push to Brainlab by itself. Karawun prep printed"
    echo "      the importTractography command to run manually (conda"
    echo "      activate KarawunDev first); it writes"
    echo "      Karawun/sub-${participant}/sub-${participant}_for_elements."
    echo "================================================================"
    echo ""
elif [ $make_dcm -eq 0 ]; then
    echo ""
    echo "Screenshots done - review them in $globalresultsdir/ before running -R to push to PACS/Karawun."
    echo ""
fi

echo "Finished"
