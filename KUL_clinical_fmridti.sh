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
     -B:  make a backup and cleanup 
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
declare -A spm_thresh_map=()

# Set required options
p_flag=0
d_flag=0


if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:d:n:v:R:F:O:a:f:T:D:S:P:E:NC:y:XBrseUQW" OPT; do

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

# Scaffold: explicit (-s) or automatic (first run for this patient, no study_config/ yet).
# Must run before the pre-flight -D check below -- that check unconditionally looks for
# study_config/${dwiprep_config_file} (default "run_dwiprep.txt") and exits with an error
# if study_config/ doesn't exist yet, which would otherwise always pre-empt the automatic
# scaffold on a genuinely fresh patient folder.
if [ $scaffold -eq 1 ]; then
    KUL_scaffold
fi
if [ ! -d $cwd/study_config ]; then
    KUL_scaffold
fi


# Pre-flight check of the -D dwiprep config: catch a lore_sd env
# misconfiguration here, before fmriprep/dwiprep are launched, instead of
# failing deep inside KUL_preproc_all.sh after other pipeline steps (and
# their downstream dependents: VBG, dwiprep_MNI, FWT) have already run.
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
if [ $fwt -eq 1 ] && ! conda env list | awk '{print $1}' | grep -qx "$scilpy" ; then
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

# NOTE: these two live here, above the `if [ $results -gt 0 ]` block below,
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




# The BACKUP and clean option
if [ $bc -eq 1 ]; then
    # clean some stuff
    clean_dwiprep="./dwiprep/sub-${participant}/sub-${participant}/*dwifsl*tmp* \
        ./dwiprep/sub-${participant}/sub-${participant}/raw \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi_orig*"
    clean_other="./fmriprep_work* \
        ./BIDS/tmp_dcm2bids"

    rm -fr $clean_dwiprep $clean_other

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

    # Compute correct PixelSpacing for each orientation from the underlay geometry.
    # mrview -size N,N renders the full image FOV into exactly N pixels per side;
    # the larger in-plane physical dimension spans mrview_resolution pixels.
    _dims=($(mrinfo $underlay -size))
    _vox=($(mrinfo $underlay -spacing))
    px_tra=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[1]}*${_vox[1]})/$mrview_resolution)")
    px_sag=$(python3 -c "print(max(${_dims[1]}*${_vox[1]},${_dims[2]}*${_vox[2]})/$mrview_resolution)")
    px_cor=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[2]}*${_vox[2]})/$mrview_resolution)")

    # Mesa software renderer threads per mrview instance.
    # 3 bundles run in parallel, so total cores = mrview_threads * 3.
    mrview_threads=$(( ncpu / 3 ))
    [ $mrview_threads -lt 1 ] && mrview_threads=1

    # --- Portable headless-mrview environment (Linux Mint / Ubuntu, GPU or not) ---
    # Every mrview capture below is prefixed with $_mrview_env. This:
    #   * strips the MATLAB MCR's Qt5 from LD_LIBRARY_PATH. That Qt5 ships no platform
    #     plugins, so if it shadows the system Qt you get the classic abort:
    #       "Could not find the Qt platform plugin xcb/offscreen in ''".
    #     All MCR lib dirs live under .../glnxa64, so removing that single token is a
    #     precise filter that leaves CUDA (and everything else) intact. It is applied
    #     ONLY to the mrview call, so any MCR-based step elsewhere keeps its runtime.
    #   * forces Mesa llvmpipe (LIBGL_ALWAYS_SOFTWARE=1) so software GL is used
    #     deterministically even on machines that DO have a GPU but no display.
    #   * drops any leaked QT_QPA_PLATFORM=offscreen / QT_PLUGIN_PATH from the shell,
    #     which would otherwise override the xcb plugin xvfb-run provides.
    # Override the binary with MRVIEW_BIN=/path/to/mrview if PATH is ambiguous.
    _mrview_ld=$(printf '%s' "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | grep -v 'glnxa64' | paste -sd:)
    _mrview_bin="${MRVIEW_BIN:-$(command -v mrview)}"
    _mrview_env="env -u QT_QPA_PLATFORM -u QT_PLUGIN_PATH LD_LIBRARY_PATH=$_mrview_ld LIBGL_ALWAYS_SOFTWARE=1"
    if [ -z "$_mrview_bin" ]; then
        echo "ERROR: mrview not found on PATH. Set MRVIEW_BIN=/path/to/mrview or fix PATH."
        _mrview_bin="mrview"   # fall through; the failure will be explicit
    fi
    if ! command -v xvfb-run >/dev/null 2>&1; then
        echo "WARNING: xvfb-run not found — mrview screenshots will fail on this host."
        echo "         Install with: sudo apt install -y xvfb libgl1-mesa-dri"
    fi

    mkdir -p $resultsdir_png
    mkdir -p $resultsdir_dcm

    # Derive underlay suffix for SPM figure directories (matches Tracto naming)
    _ulsuffix="${resultsdir_png##*_figures_}"
    spm_resultsdir_png="$globalresultsdir/SPM_figures_${_ulsuffix}"
    spm_resultsdir_dcm="$globalresultsdir/PACS/fMRI_${_ulsuffix}"
    mkdir -p "$spm_resultsdir_png"
    mkdir -p "$spm_resultsdir_dcm"

    # Render one SPM/Melodic map (all orientations) on display :20.
    # Threshold is computed per-map as max/3 (same as the report logic).
    _render_one_spm() {
        local spmfile="$1" spmname="$2"
        local ori
        IFS=',' read -ra ori <<< "$orientations"

        local _max_T _thresh
        if [ -n "${spm_thresh_map[$spmname]+x}" ]; then
            _thresh=${spm_thresh_map[$spmname]}
            echo "SPM ${spmname}: using per-map threshold=${_thresh}"
        elif [ -n "$spm_thresh_override" ]; then
            _thresh=$spm_thresh_override
            echo "SPM ${spmname}: using fixed threshold=${_thresh}"
        else
            _max_T=$(mrstats -output max "$spmfile")
            _thresh=$(awk "BEGIN {print $_max_T/3}")
            echo "SPM ${spmname}: max=${_max_T}, auto threshold=${_thresh}"
        fi

        # Optionally compute a 1-voxel edge mask for the dark-blue outline
        local _edge_overlay="" _tmp_mask="" _tmp_eroded="" _tmp_edge=""
        if [ $spm_edge -eq 1 ]; then
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
            if [[ "$orient" == "TRA" ]]; then
                underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF)}'); plane=2
            elif [[ "$orient" == "SAG" ]]; then
                underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF-2)}'); plane=0
            else
                underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF-1)}'); plane=1
            fi

            local png_dir="$spm_resultsdir_png/${spmname}_${orient}"
            if [ -n "$(ls "$png_dir"/*.png 2>/dev/null | head -1)" ]; then
                echo "Skipping screenshots for ${spmname}_${orient} (already exist)"
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
                eval "$_mrview_env LP_NUM_THREADS=$mrview_threads timeout 600 xvfb-run -n 20 --server-args=\"-screen 0 ${mrview_resolution}x${mrview_resolution}x24\" $_mrview_bin -size $mrview_resolution,$mrview_resolution \
                    -load $underlay -mode 1 -plane $plane \
                    -overlay.load $spmfile -overlay.opacity $spm_opacity -overlay.colourmap 1 \
                        -overlay.threshold_min $_thresh \
                    $_edge_overlay \
                    -noannotations -orientlabel 0 -voxelinfo 0 -colourbar 0 \
                    -capture.folder $png_dir -capture.prefix ${spmname}_${orient} \
                    $voxel_index -force -exit"
                [ -z "$_sample_png" ] && _sample_png=$(ls "$png_dir"/*.png 2>/dev/null | head -1)
            fi
        done

        [ -n "$_tmp_mask" ] && rm -f "$_tmp_mask" "$_tmp_eroded" "$_tmp_edge"

        # Recompute PixelSpacing from first screenshot (same as tract logic)
        if [ -n "$_sample_png" ]; then
            local _png_max
            _png_max=$(python3 -c "from PIL import Image; w,h=Image.open('$_sample_png').size; print(max(w,h))")
            _px_tra=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[1]}*${_vox[1]})/$_png_max)")
            _px_sag=$(python3 -c "print(max(${_dims[1]}*${_vox[1]},${_dims[2]}*${_vox[2]})/$_png_max)")
            _px_cor=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[2]}*${_vox[2]})/$_png_max)")
        fi

        # DICOM conversion (only when explicitly requested via -R)
        # Skip if this map was not selected by the user (_dcm_spm_set empty = all selected)
        if [ $make_dcm -eq 1 ] && [ -n "$donor_dcm" ] && \
           { [ ${#_dcm_spm_set[@]} -eq 0 ] || [ -n "${_dcm_spm_set[$spmname]+x}" ]; }; then
            local dcm_label="${spmname}_on_${_ulsuffix}"
            for orient in "${ori[@]}"; do
                local ps
                case $orient in TRA) ps=$_px_tra ;; SAG) ps=$_px_sag ;; COR) ps=$_px_cor ;; esac
                local png_dir="$spm_resultsdir_png/${spmname}_${orient}"
                local dcmdir="$spm_resultsdir_dcm/${dcm_label}_${orient}"
                if [ -n "$(ls "$dcmdir"/*.dcm 2>/dev/null | head -1)" ]; then
                    echo "Skipping DICOMs for ${dcm_label}_${orient} (already exist)"
                else
                    mkdir -p "$dcmdir"
                    if [[ "$orient" == "SAG" ]]; then
                        echo "Making dicoms in $dcmdir (donor-match mode)"
                        KUL_nii2dcm.py -s "${dcm_label}_${orient}" \
                            -u "$underlay" -o "$orient" -M \
                            "$png_dir" "$donor_dcm" "$dcmdir" \
                            || echo "WARNING: KUL_nii2dcm.py failed for ${dcm_label}_${orient}"
                    else
                        echo "Making dicoms in $dcmdir (PixelSpacing=${ps}mm)"
                        KUL_nii2dcm.py -s "${dcm_label}_${orient}" -p $ps \
                            -u "$underlay" -o "$orient" \
                            "$png_dir" "$donor_dcm" "$dcmdir" \
                            || echo "WARNING: KUL_nii2dcm.py failed for ${dcm_label}_${orient}"
                    fi
                fi
            done
        fi
    }

    # Render one bundle (all orientations sequentially) on display :$((10+slot))
    _render_one_bundle() {
        local slot="$1" tractname="$2" mrview_tck="$3"
        local ori
        IFS=',' read -ra ori <<< "$orientations"

        for orient in "${ori[@]}"; do
            if [ -n "$(ls "$resultsdir_png/${tractname}_${orient}"/*.png 2>/dev/null | head -1)" ]; then
                echo "Skipping screenshots for ${tractname}_${orient} (already exist)"
            else
                local underlay_slices plane
                if [[ "$orient" == "TRA" ]]; then
                    underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF)}'); plane=2
                elif [[ "$orient" == "SAG" ]]; then
                    underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF-2)}'); plane=0
                else
                    underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF-1)}'); plane=1
                fi
                echo "Making ${tractname}_${orient} on $(basename $underlay)"
                mkdir -p "$resultsdir_png/${tractname}_${orient}"
                local voxel_index="-capture.folder $resultsdir_png/${tractname}_${orient} -capture.prefix ${tractname}_${orient}"
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
                eval "$_mrview_env LP_NUM_THREADS=$mrview_threads timeout 600 xvfb-run -n $((10 + slot)) --server-args=\"-screen 0 ${mrview_resolution}x${mrview_resolution}x24\" $_mrview_bin -size $mrview_resolution,$mrview_resolution \
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
            _png_max=$(python3 -c "from PIL import Image; w,h=Image.open('$_sample_png').size; print(max(w,h))")
            _px_tra=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[1]}*${_vox[1]})/$_png_max)")
            _px_sag=$(python3 -c "print(max(${_dims[1]}*${_vox[1]},${_dims[2]}*${_vox[2]})/$_png_max)")
            _px_cor=$(python3 -c "print(max(${_dims[0]}*${_vox[0]},${_dims[2]}*${_vox[2]})/$_png_max)")
            echo "Actual PNG max dim: ${_png_max}px → PixelSpacing TRA=${_px_tra}mm SAG=${_px_sag}mm COR=${_px_cor}mm"
        fi

        # DICOM conversion (only when explicitly requested via -R)
        if [ $make_dcm -eq 1 ] && [ -n "$donor_dcm" ]; then
            for orient in "${ori[@]}"; do
                local ps
                case $orient in TRA) ps=$_px_tra ;; SAG) ps=$_px_sag ;; COR) ps=$_px_cor ;; esac
                local dcmdir="$resultsdir_dcm/${tractname}_${orient}"
                if [ -n "$(ls "$dcmdir"/*.dcm 2>/dev/null | head -1)" ]; then
                    echo "Skipping DICOMs for ${tractname}_${orient} (already exist)"
                else
                    mkdir -p "$dcmdir"
                    if [[ "$orient" == "SAG" ]]; then
                        echo "Making dicoms in $dcmdir (donor-match mode)"
                        KUL_nii2dcm.py -s "FT_${tractname}_${orient}_${_ulsuffix}" \
                            -u "$underlay" -o "$orient" -M \
                            "$resultsdir_png/${tractname}_${orient}" "$donor_dcm" "$dcmdir" \
                            || echo "WARNING: KUL_nii2dcm.py failed for ${tractname}_${orient}"
                    else
                        echo "Making dicoms in $dcmdir (PixelSpacing=${ps}mm)"
                        KUL_nii2dcm.py -s "FT_${tractname}_${orient}_${_ulsuffix}" -p $ps \
                            -u "$underlay" -o "$orient" \
                            "$resultsdir_png/${tractname}_${orient}" "$donor_dcm" "$dcmdir" \
                            || echo "WARNING: KUL_nii2dcm.py failed for ${tractname}_${orient}"
                    fi
                fi
            done
        fi
    }

    # Launch current batch of bundles in parallel (max 3), then wait
    _flush_bundle_batch() {
        for slot in "${!_btractnames[@]}"; do
            _render_one_bundle "$slot" "${_btractnames[$slot]}" "${_btcks[$slot]}" &
        done
        wait
        _btractnames=()
        _btcks=()
    }

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

    # ── SPM & Melodic → PACS ────────────────────────────────────────────────
    # If interactive and no global -T override, ask for a per-map threshold list.
    _all_spm_names=()
    for _spm in "$globalresultsdir/SPM/"*.nii.gz "$globalresultsdir/SPM/"*.nii \
                "$globalresultsdir/Melodic/"*.nii.gz "$globalresultsdir/Melodic/"*.nii; do
        [ -f "$_spm" ] || continue
        _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
        _all_spm_names+=("$_spmname")
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

    for _spm in "$globalresultsdir/SPM/"*.nii.gz "$globalresultsdir/SPM/"*.nii; do
        [ -f "$_spm" ] || continue
        _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
        _render_one_spm "$_spm" "$_spmname"
    done

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
            # fe3b800 merged left/right colours for well-separated tracts to
            # free budget "for fMRI labels without colliding", but the labels
            # themselves were never moved into the freed range. This finishes
            # that: the low values here are exactly the indices the tract table,
            # the DBS VOIs (23, 24) and the lesion (16) leave unused, so it
            # works on stock karawun. Anything past those spills into 50+,
            # which needs the extended-palette fork (values >30 clamp to the
            # last entry on stock karawun).
            #
            # The low seven are ordered by measured worst-case CIEDE2000
            # distance to every tract/lesion colour and to each other, best
            # first, so the slots actually used on a typical 1-3 task case get
            # the most distinguishable colours:
            #   2 -> 15.3   6 -> 12.7   12 -> 9.9   8 -> 7.2
            #   14 -> 6.9   10 -> 6.1   30 -> 4.1
            # For reference, fMRI colour 2 sits dE 77.7 from the Arcuate, which
            # is the comparison that motivated all of this.
            #
            # The fork's 50+ colours are better separated (~13) than low slots
            # 4-7, but they are listed after, not before: on stock karawun
            # anything >30 clamps to a single entry, so leading with them would
            # make every task past the third render identically. Degrading to a
            # merely-mediocre colour is a better failure than degrading to no
            # distinction at all.
            _fmri_label_colors=(2 6 12 8 14 10 30 50 51 52 53 54 55 56 57 58)

            for _idx in "${!_spm_task_names_sorted[@]}"; do
                _spmname="${_spm_task_names_sorted[$_idx]}"
                if [ $_idx -lt ${#_fmri_label_colors[@]} ]; then
                    _label_int=${_fmri_label_colors[$_idx]}
                else
                    # more tasks than reserved colours; keep going past the end
                    # of the list rather than silently reusing one
                    _label_int=$((59 + _idx - ${#_fmri_label_colors[@]}))
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
        echo "Karawun/sub-${participant}/T1w.nii.gz not found — skipping Karawun fMRI labels (run Karawun prep first)"
    fi

    for _spm in "$globalresultsdir/Melodic/"*.nii.gz "$globalresultsdir/Melodic/"*.nii; do
        [ -f "$_spm" ] || continue
        _spmname=$(basename "$_spm"); _spmname=${_spmname%.nii.gz}; _spmname=${_spmname%.nii}
        _render_one_spm "$_spm" "$_spmname"
    done

    # ── Lesion & DSC perfusion → PACS ───────────────────────────────────────
    # Same renderer as the fMRI maps, pointed at a separate output folder so
    # these don't get mixed in with the SPM/Melodic series on PACS.
    _extra_names=(); _extra_files=(); _extra_thresh=()

    if _lesion_pacs=$(KUL_resolve_lesion); then
        # a mask, so any threshold in (0,1) selects it; 0.5 is the obvious one
        _extra_names+=("Lesion"); _extra_files+=("$_lesion_pacs"); _extra_thresh+=("0.5")
    fi

    # Fixed thresholds, not the auto max/3 used for the fMRI maps. These are
    # NAWM-normalised ratios, so a threshold has a fixed clinical meaning that
    # max/3 would throw away: 1.75 is the conventional high-grade glioma rCBV
    # cutoff, and 1.0 on nrCBF is simply "above contralesional normal WM".
    # -T still overrides both.
    _perf_dir="$globalresultsdir/Perfusion"
    for _pm in nrCBV_corrected:1.75 nrCBF:1.0; do
        _pf="$_perf_dir/sub-${participant}_${_pm%%:*}.nii.gz"
        [ -f "$_pf" ] || continue
        _extra_names+=("${_pm%%:*}"); _extra_files+=("$_pf"); _extra_thresh+=("${_pm##*:}")
    done

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

        for _i in "${!_extra_names[@]}"; do
            if [ -n "${_extra_thresh[$_i]}" ] && [ -z "$spm_thresh_override" ] && \
               [ -z "${spm_thresh_map[${_extra_names[$_i]}]+x}" ]; then
                spm_thresh_map["${_extra_names[$_i]}"]="${_extra_thresh[$_i]}"
            fi
            echo "PACS extra: ${_extra_names[$_i]} <- $(basename "${_extra_files[$_i]}")"
            _render_one_spm "${_extra_files[$_i]}" "${_extra_names[$_i]}"
        done

        spm_resultsdir_png="$_saved_png"; spm_resultsdir_dcm="$_saved_dcm"
        unset _dcm_spm_set; declare -A _dcm_spm_set=()
        for _k in "${_saved_sel[@]}"; do _dcm_spm_set["$_k"]=1; done
    fi

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
                rm -fr $derivativesdir/KUL_FWT/* >/dev/null 2>&1
            fi
        fi

        # Step 15 — figures
        read -p "Redo: figures? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            rm -f ${cwd}/KUL_LOG/sub-${participant}_figures.done >/dev/null 2>&1
        fi

    fi
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
        
        # copying the result to the global results dir
        cp -f fmriprep/sub-$participant/anat/sub-${participant}_desc-preproc_T1w.nii.gz $globalresultsdir/Anat/T1w_fmriprep.nii.gz
        #gunzip -f $globalresultsdir/Anat/T1w.nii.gz
        
        # create a GM mask in the global results dir
        mrcalc fmriprep/sub-$participant/anat/sub-${participant}_dseg.nii.gz 1 -eq \
            fmriprep/sub-$participant/anat/sub-${participant}_dseg.nii.gz -mul - | \
            maskfilter - median - | \
            maskfilter - dilate $globalresultsdir/Anat/T1w_GM.nii.gz

        # add to the report
        if [ -f fmriprep/sub-${participant}.html ]; then
            ln -s ${cwd}/fmriprep/sub-${participant}.html ${cwd}/REPORT/sub-${participant}_03_fmriprep.html
        fi

    else
        echo "Fmriprep already done"
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
            # by setup_environment.sh) instead of a separately hardcoded path
            # here, which can silently drift out of sync with the actual
            # install location (as it did: this used to point at a directory
            # that no longer exists).
            if [ -z "$FREESURFER_HOME" ]; then
                kul_echo "ERROR: FREESURFER_HOME is not set. Source setup_environment.sh before running this pipeline."
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

            # add to report using the hardwired wc (with-confounds) Bizzi-thresholded
            # maps only (p<0.001 unc, k>=50) — matches the RESULTS/.../SPM selection
            for bizzi_map in $derivativesdir/SPM/*_wc/spmT_0001_p001unc_k50.nii; do
                [ -f "$bizzi_map" ] || continue
                task=$(basename $(dirname $bizzi_map))
                task=${task%_wc}
                KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz \
                    -o "$bizzi_map" -t 2 -d REPORT -f 05_afMRI_${task}_p001unc_k50
            done
        fi

        if [ ! -f ${cwd}/KUL_LOG/sub-${participant}_melodic.done ]; then
            task_in="KUL_fmriproc_conn.sh -p $participant"
            KUL_task_exec $verbose_level "KUL_fmriproc_conn" "8_fmriproc_conn" || kul_echo "KUL_fmriproc_conn failed for sub-${participant} — melodic.done will not be created (check 8_fmriproc_conn.error.log)"

            # add to report
            for spm in RESULTS/sub-${participant}/Melodic/*.nii; do
                #echo $spm
                max_T=$(mrstats -output max $spm)
                #echo $max_T
                thresh=$(awk "BEGIN {print $max_T/3}")
                task=$(basename $spm)
                mrcalc $spm $thresh -gt REPORT/spm_tmp_$task
                KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz -o REPORT/spm_tmp_$task \
                    -t 2 -d REPORT -f 05_rsfMRI_${task}_Thr_${thresh}
                rm -f REPORT/spm_tmp_$task
            done
        fi
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

        if [ -f dwiprep/sub-${participant}/sub-${participant}/qa/sub-${participant}_T1w_with_fa.png ]; then
            cp -f dwiprep/sub-${participant}/sub-${participant}/qa/sub-${participant}_T1w_with_fa.png \
                REPORT/sub-${participant}_02_T1w_with_fa.png
            cp -f dwiprep/sub-${participant}/sub-${participant}/qa/sub-${participant}_T1w_brain_with_fa.png \
                REPORT/sub-${participant}_02_T1w_brain_with_fa.png
            cp -f dwiprep/sub-${participant}/sub-${participant}/eddy_qc/quad/qc.pdf REPORT/sub-${participant}_02_eddy_qc.pdf
        fi

        touch $dwi_anat_check
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

# STEP 8b - run rsfMRI network analysis (opt-in, -N)
KUL_run_rsfMRI_networks

# STEP 9 - run VBG
KUL_run_VBG 
wait

# STEP 9b - FastSurfer + KUL_multiparc (types 4, 5, 6 only — no VBG)
# this should only run if VBG is not used!!
KUL_run_multiparc
wait

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

# STEP 15 - Prepare Karawun folder. No longer automatic: like the PACS
# DICOM export, this only runs when explicitly requested via -R (review
# the FWT tract/VOI output first). The importTractography command is
# printed at the end — run it manually after copying a donor DICOM into
# Karawun/sub-${participant}/DICOM/ (or RESULTS/sub-${participant}/DICOM/)
# and curating the fMRI label maps.
if [ $make_dcm -eq 1 ]; then
    karawun_prepare_check=${cwd}/KUL_LOG/sub-${participant}_karawun_prepare.done
    if [ ! -f $karawun_prepare_check ]; then
        if [ $type -lt 5 ]; then
            kul_echo "Preparing Karawun folder"
            KUL_karawun_prepare.sh -p ${participant} -t 1 -r 3
        elif [ $type -eq 5 ]; then
            kul_echo "Preparing Karawun folder (DBS ET)"
            KUL_karawun_prepare.sh -p ${participant} -t 2 -r 10
        elif [ $type -eq 6 ]; then
            kul_echo "Preparing Karawun folder (DBS Parkinson)"
            KUL_karawun_prepare.sh -p ${participant} -t 3 -r 10
        fi
        touch $karawun_prepare_check
    else
        echo "Karawun folder already prepared"
    fi
else
    echo "Karawun folder prep skipped (run with -R to prepare it)"
fi


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
    echo "        BIDS/derivatives/KUL_compute/sub-${participant}/FWT/"
    echo "      The fMRI/tract thresholds used for PACS export are only"
    echo "      computed when you run -F or -R below (auto = max/3 per"
    echo "      map, or interactively/-T if you override) - this review"
    echo "      is your chance to catch a threshold that looks wrong"
    echo "      before it is baked into the exported DICOMs."
    echo ""
    echo "   2. Generate screenshots only, no PACS/Karawun push yet:"
    echo "        KUL_clinical_fmridti.sh -p ${participant} -t ${type} -F <1-7> [-O orientations]"
    echo ""
    echo "   3. Before running -R, create RESULTS/sub-${participant}/DICOM/"
    echo "      (it is not made for you) and drop a donor DICOM into it -"
    echo "      a single file is enough, one slice from a high-resolution"
    echo "      anatomical series (T1w, FLAIR, T2, ...). The same donor is"
    echo "      used for both the PACS and the Karawun/Brainlab output, so"
    echo "      keep it to one file/series for consistency."
    echo ""
    echo "   4. -R also triggers Karawun prep (it no longer runs on its"
    echo "      own either). Once happy with the review and the donor"
    echo "      DICOM is in place, run:"
    echo "        KUL_clinical_fmridti.sh -p ${participant} -t ${type} -R <1-7> [-O orientations]"
    echo "================================================================"
    echo ""
elif [ $make_dcm -eq 0 ]; then
    echo ""
    echo "Screenshots done - review them in $globalresultsdir/ before running -R to push to PACS/Karawun."
    echo ""
fi

echo "Finished"
