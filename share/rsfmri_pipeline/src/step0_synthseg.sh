#!/usr/bin/env bash
# step0_synthseg.sh
# Run mri_synthseg, 5ttgen (FreeSurfer mode), MNI warping, union tissue maps,
# and subject-specific Yeo17 atlas for every subject's fMRIPrep T1w (native space).
#
# Usage:
#   bash step0_synthseg.sh --subjects <ID> [<ID> ...] [--force | --force2]
#   bash step0_synthseg.sh -h
#
# Options:
#   --subjects  one or more subject IDs (required, e.g. HV01 HV02 PT01)
#   --force     reprocess ALL steps, including mri_synthseg (slow)
#   --force2    reprocess steps 1-5 only; skip mri_synthseg if output exists
#   -h, --help  show this help message and exit
#
# Output per subject (analysis/synthseg/sub-<ID>/):
#   sub-<ID>_synthseg.nii.gz                 parcellation (native)
#   sub-<ID>_synthseg_vols.csv               volume table
#   sub-<ID>_synthseg_qc.csv                 QC metrics
#   sub-<ID>_5tt.nii.gz                      5-tissue-type image (native)
#   sub-<ID>_synthseg_MNI.nii.gz             parcellation (MNI)
#   sub-<ID>_5tt_MNI.nii.gz                  5TT image (MNI)
#   sub-<ID>_union_GM.nii.gz                 union GM mask (native)
#   sub-<ID>_union_WM.nii.gz                 union WM mask (native)
#   sub-<ID>_union_CSF.nii.gz                union CSF mask (native)
#   sub-<ID>_Yeo17_subject_specific_MNI.nii.gz  Yeo17 propagated through subject GM (MNI)

set -euo pipefail
# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"   # fallback for standalone use (resolves to KUL_NIS/ after relocation)
PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
FMRIPREP_DIR="${RSFMRI_FMRIPREP_DIR:-${BASE_DIR}/fmriprep}"
SYNTHSEG_DIR="${RSFMRI_SYNTHSEG_DIR:-${BASE_DIR}/analysis/synthseg}"
LOG_DIR="${SYNTHSEG_DIR}/logs"
# Standard freesurfer derivatives location; per-subject VBG output (preferred
# when present, same precedence KUL_run_FWT already uses) is checked directly
# in the per-subject loop below since it needs the subject ID either way.
FS_DIR="${RSFMRI_FS_DIR:-${BASE_DIR}/BIDS/derivatives/freesurfer}"
THREADS=4

YEO_DIR="${PIPELINE_DIR}/Yeo_JNeurophysiol11_MNI152"
YEOBUCKNER17="${YEO_DIR}/Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz"
# YEOBUCKNER17RS is set per-subject inside the main loop below (its correct
# target grid is subject-specific -- see the per-subject declaration for why).

# FSL atlases for subject-specific SBA propagation
# Use $FSLDIR if already set (standard when FSL is sourced); otherwise derive
# from the location of flirt on PATH.
if [[ -z "${FSLDIR:-}" ]]; then
    _flirt=$(command -v flirt 2>/dev/null || true)
    if [[ -z "${_flirt}" ]]; then
        echo "ERROR: FSLDIR is not set and flirt not found on PATH." >&2
        echo "       Source your FSL setup script (e.g. source \$FSLDIR/etc/fslconf/fsl.sh)." >&2
        exit 1
    fi
    FSLDIR=$(dirname "$(dirname "$(realpath "${_flirt}")")")
    echo "FSLDIR resolved from PATH: ${FSLDIR}"
fi
FSL_ATLAS_DIR="${FSLDIR}/data/atlases"
HO_CORT_ATLAS="${FSL_ATLAS_DIR}/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz"
FSL_STRIATUM_ATLAS="${FSL_ATLAS_DIR}/Striatum/striatum-con-label-thr25-7sub-2mm.nii.gz"
SHARED_ATLAS_DIR="${SYNTHSEG_DIR}/atlases_2mm"
HO_CORT_RS="${SHARED_ATLAS_DIR}/ho_cort_rs2mm.nii.gz"
FSL_STRIATUM_RS="${SHARED_ATLAS_DIR}/fsl_striatum_rs2mm.nii.gz"

mkdir -p "${LOG_DIR}" "${SHARED_ATLAS_DIR}"

# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------
usage() {
    sed -n '2,/^$/p' "$0" | grep '^#' | sed 's/^# \?//'
    exit 0
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
SUBJECTS=()
FORCE=0
FORCE2=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        --subjects)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                SUBJECTS+=("$1")
                shift
            done
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --force2)
            FORCE2=1
            shift
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# FORCE_DS (downstream): true when either flag is set — used for steps 1-5
FORCE_DS=$(( FORCE || FORCE2 ))

if [[ ${#SUBJECTS[@]} -eq 0 ]]; then
    echo "ERROR: --subjects is required." >&2
    echo "Run with -h for usage." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Check required tools and static inputs
# ---------------------------------------------------------------------------
for tool in mri_synthseg 5ttgen antsApplyTransforms ImageMath mrconvert mrcalc mrinfo mri_convert; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "ERROR: ${tool} not found on PATH." >&2
        exit 1
    fi
done

if [[ ! -f "${YEOBUCKNER17}" ]]; then
    echo "ERROR: Yeo17 atlas not found at ${YEOBUCKNER17}" >&2
    exit 1
fi
echo "Yeo17 atlas  : ${YEOBUCKNER17}"

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
ts() { date '+%H:%M:%S'; }

img_info() {
    local F="$1"
    [[ -f "${F}" ]] || return 0
    local SIZE SPACING
    SIZE=$(mrinfo    "${F}" -size    2>/dev/null || echo "?")
    SPACING=$(mrinfo "${F}" -spacing 2>/dev/null || echo "?")
    echo "    dim: ${SIZE}   vox: ${SPACING} mm"
}

log_tail() {
    local LOG="$1"
    [[ -f "${LOG}" ]] || return 0
    echo "    --- last lines of ${LOG##*/} ---" >&2
    tail -10 "${LOG}" | sed 's/^/    | /' >&2
}

step_header() {
    local N="$1" LABEL="$2"
    echo "  -- [$(ts)] Step ${N}: ${LABEL}"
}

# ---------------------------------------------------------------------------
# Helper: warp a label image to MNI with GenericLabel interpolation
# Usage: warp_to_mni <input> <output> <xfm> <log> [etype]
# ---------------------------------------------------------------------------
warp_to_mni() {
    local INPUT="$1" OUTPUT="$2" XFM="$3" LOG="$4" ETYPE="${5:-0}"
    if [[ ! -f "${XFM}" ]]; then
        echo "    SKIP: transform not found at ${XFM}" >&2
        return 1
    fi
    echo "    input : ${INPUT##*/}"
    echo "    output: ${OUTPUT##*/}"
    echo "    xfm   : ${XFM##*/}"
    if antsApplyTransforms \
            -d 3 \
            -e "${ETYPE}" \
            -i "${INPUT}" \
            -r "${MNI_REF}" \
            -t "${XFM}" \
            -n GenericLabel \
            -o "${OUTPUT}" \
            > "${LOG}" 2>&1; then
        echo "    DONE ($(ts))"
        img_info "${OUTPUT}"
    else
        echo "    FAILED — see ${LOG}" >&2
        log_tail "${LOG}"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Process each subject
# ---------------------------------------------------------------------------
N_DONE=0
N_SKIP=0
N_FAIL=0
N_TOTAL=${#SUBJECTS[@]}
N_SUBJ=0

for SUBJ in "${SUBJECTS[@]}"; do

    (( N_SUBJ++ )) || true

    T1W="${FMRIPREP_DIR}/sub-${SUBJ}/anat/sub-${SUBJ}_desc-preproc_T1w.nii.gz"
    # Prefer explicit-resolution (res-2) fmriprep anat output if present; fall
    # back to whatever resolution fmriprep actually produced (not every
    # fmriprep run is invoked with --output-spaces ...:res-2). Resolved once
    # and reused for T1W/GM/WM/CSF below so they all share one voxel grid --
    # they're combined voxelwise later and must not come from mixed resolutions.
    if [[ -f "${FMRIPREP_DIR}/sub-${SUBJ}/anat/sub-${SUBJ}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz" ]]; then
        RES_SUFFIX="_res-2"
    else
        RES_SUFFIX=""
    fi
    MNI_REF="${FMRIPREP_DIR}/sub-${SUBJ}/anat/sub-${SUBJ}_space-MNI152NLin2009cAsym${RES_SUFFIX}_desc-preproc_T1w.nii.gz"
    OUT_DIR="${SYNTHSEG_DIR}/sub-${SUBJ}"

    OUT_SEG="${OUT_DIR}/sub-${SUBJ}_synthseg.nii.gz"
    OUT_VOL="${OUT_DIR}/sub-${SUBJ}_synthseg_vols.csv"
    OUT_QC="${OUT_DIR}/sub-${SUBJ}_synthseg_qc.csv"
    LOG_SEG="${LOG_DIR}/sub-${SUBJ}_synthseg.log"

    OUT_5TT="${OUT_DIR}/sub-${SUBJ}_5tt.nii.gz"
    LOG_5TT="${LOG_DIR}/sub-${SUBJ}_5ttgen.log"

    OUT_MNI_SEG="${OUT_DIR}/sub-${SUBJ}_synthseg_MNI.nii.gz"
    OUT_MNI_5TT="${OUT_DIR}/sub-${SUBJ}_5tt_MNI.nii.gz"
    XFM="${FMRIPREP_DIR}/sub-${SUBJ}/anat/sub-${SUBJ}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5"

    FP_GM="${FMRIPREP_DIR}/sub-${SUBJ}/anat/sub-${SUBJ}_space-MNI152NLin2009cAsym${RES_SUFFIX}_label-GM_probseg.nii.gz"
    FP_WM="${FMRIPREP_DIR}/sub-${SUBJ}/anat/sub-${SUBJ}_space-MNI152NLin2009cAsym${RES_SUFFIX}_label-WM_probseg.nii.gz"
    FP_CSF="${FMRIPREP_DIR}/sub-${SUBJ}/anat/sub-${SUBJ}_space-MNI152NLin2009cAsym${RES_SUFFIX}_label-CSF_probseg.nii.gz"
    OUT_UNION_GM="${OUT_DIR}/sub-${SUBJ}_union_GM.nii.gz"
    OUT_UNION_WM="${OUT_DIR}/sub-${SUBJ}_union_WM.nii.gz"
    OUT_UNION_CSF="${OUT_DIR}/sub-${SUBJ}_union_CSF.nii.gz"
    LOG_UNION="${LOG_DIR}/sub-${SUBJ}_union.log"

    OUT_YEO_SUBJ_MNI="${OUT_DIR}/sub-${SUBJ}_Yeo17_subject_specific_MNI.nii.gz"
    # Per-subject, NOT shared: union_GM's grid is derived from this subject's
    # own native T1w acquisition matrix and differs subject to subject, so a
    # single cached-once resample target (the old behaviour) silently mismatches
    # every subject's union_GM except whichever one populated the cache first --
    # ImageMath PropagateLabelsThroughMask doesn't validate matching grids, so
    # that mismatch corrupted (partially or fully zeroed) networks silently.
    YEOBUCKNER17RS="${OUT_DIR}/sub-${SUBJ}_Yeo17_rs2mm.nii.gz"
    LOG_YEO="${LOG_DIR}/sub-${SUBJ}_yeo17.log"

    OUT_HO_CORT_SUBJ="${OUT_DIR}/sub-${SUBJ}_ho_cort_subject_MNI.nii.gz"
    OUT_FSL_STRIATUM_SUBJ="${OUT_DIR}/sub-${SUBJ}_fsl_striatum_subject_MNI.nii.gz"
    LOG_ATLAS="${LOG_DIR}/sub-${SUBJ}_atlas_propagation.log"

    # VBG's own FreeSurfer output is preferred when present, same precedence
    # KUL_run_FWT already uses for aparc+aseg -- fall back to the standard
    # freesurfer derivatives location otherwise.
    LAUSANNE_VBG="${BASE_DIR}/KUL_VBG/output_VBG/sub-${SUBJ}_FS_output/sub-${SUBJ}/mri/lausanne2018.scale3+aseg.mgz"
    LAUSANNE_STD="${FS_DIR}/sub-${SUBJ}/mri/lausanne2018.scale3+aseg.mgz"
    OUT_LAUSANNE_T1W="${OUT_DIR}/sub-${SUBJ}_lausanne_scale3_T1w.nii.gz"
    OUT_LAUSANNE_MNI="${OUT_DIR}/sub-${SUBJ}_lausanne_scale3_MNI.nii.gz"
    LOG_LAUSANNE="${LOG_DIR}/sub-${SUBJ}_lausanne_scale3.log"

    echo "============================================================"
    echo "Subject ${N_SUBJ}/${N_TOTAL}: sub-${SUBJ}  [$(ts)]"
    echo "============================================================"

    # Check inputs exist
    if [[ ! -f "${T1W}" ]]; then
        echo "  ERROR: T1w not found: ${T1W}" >&2
        (( N_FAIL++ )) || true
        continue
    fi
    if [[ ! -f "${MNI_REF}" ]]; then
        echo "  ERROR: fMRIPrep MNI T1w not found: ${MNI_REF}" >&2
        (( N_FAIL++ )) || true
        continue
    fi
    echo "  T1w    : ${T1W}"
    echo "  MNI ref: ${MNI_REF}"

    mkdir -p "${OUT_DIR}"

    # -----------------------------------------------------------------------
    # Step 0: mri_synthseg
    # --robust is important for patients who may have brain lesions;
    # it has negligible cost on healthy volunteers.
    # -----------------------------------------------------------------------
    step_header "0/5" "mri_synthseg"
    if [[ ! -f "${OUT_SEG}" || "${FORCE}" -eq 1 ]]; then
        echo "    input : ${T1W}"
        echo "    output: ${OUT_SEG}"
        echo "    log   : ${LOG_SEG}"
        if mri_synthseg \
                --i  "${T1W}"     \
                --o  "${OUT_SEG}" \
                --parc            \
                --robust          \
                --vol "${OUT_VOL}" \
                --qc  "${OUT_QC}" \
                --threads "${THREADS}" \
                > "${LOG_SEG}" 2>&1; then
            echo "    DONE ($(ts))"
            img_info "${OUT_SEG}"
            (( N_DONE++ )) || true
        else
            echo "    FAILED — see ${LOG_SEG}" >&2
            log_tail "${LOG_SEG}"
            (( N_FAIL++ )) || true
            continue
        fi
    else
        echo "    SKIP: output exists (use --force to reprocess)"
        img_info "${OUT_SEG}"
        (( N_SKIP++ )) || true

    fi

    # -----------------------------------------------------------------------
    # Step 1: 5ttgen freesurfer (driven by SynthSeg parcellation)
    # -----------------------------------------------------------------------
    step_header "1/5" "5ttgen freesurfer"
    if [[ ! -f "${OUT_5TT}" || "${FORCE_DS}" -eq 1 ]]; then
        echo "    input : ${OUT_SEG##*/}"
        echo "    output: ${OUT_5TT##*/}"
        echo "    log   : ${LOG_5TT##*/}"
        if 5ttgen freesurfer "${OUT_SEG}" "${OUT_5TT}" -force \
                > "${LOG_5TT}" 2>&1; then
            echo "    DONE ($(ts))"
            img_info "${OUT_5TT}"
        else
            echo "    FAILED — see ${LOG_5TT}" >&2
            log_tail "${LOG_5TT}"
            (( N_FAIL++ )) || true
        fi
    else
        echo "    SKIP: output exists (use --force or --force2 to reprocess)"
        img_info "${OUT_5TT}"
    fi

    # -----------------------------------------------------------------------
    # Step 2: Warp SynthSeg parcellation and 5TT to MNI space
    # -----------------------------------------------------------------------
    step_header "2/5" "warp to MNI"

    if [[ ! -f "${OUT_MNI_SEG}" || "${FORCE_DS}" -eq 1 ]]; then
        echo "  [2a] SynthSeg parcellation → MNI"
        warp_to_mni "${OUT_SEG}" "${OUT_MNI_SEG}" "${XFM}" \
            "${LOG_DIR}/sub-${SUBJ}_synthseg_mni.log" || (( N_FAIL++ )) || true
    else
        echo "  [2a] SynthSeg MNI SKIP: output exists"
        img_info "${OUT_MNI_SEG}"
    fi

    if [[ -f "${OUT_5TT}" ]]; then
        if [[ ! -f "${OUT_MNI_5TT}" || "${FORCE_DS}" -eq 1 ]]; then
            echo "  [2b] 5TT → MNI"
            warp_to_mni "${OUT_5TT}" "${OUT_MNI_5TT}" "${XFM}" \
                "${LOG_DIR}/sub-${SUBJ}_5tt_mni.log" 3 || (( N_FAIL++ )) || true
        else
            echo "  [2b] 5TT MNI SKIP: output exists"
            img_info "${OUT_MNI_5TT}"
        fi
    else
        echo "  [2b] 5TT MNI SKIP: ${OUT_5TT##*/} not found" >&2
    fi

    # -----------------------------------------------------------------------
    # Step 3: Union tissue maps (5ttgen volumes + fMRIPrep probseg)
    # union_GM  = (5tt cortex vol-0 + BGT vol-1 > 0) | (fMRIPrep GM  > 0.25)
    # union_WM  = (5tt WM  vol-2 > 0)                | (fMRIPrep WM  > 0.25)
    # union_CSF = (5tt CSF vol-3 > 0)                | (fMRIPrep CSF > 0.25)
    # -----------------------------------------------------------------------
    step_header "3/5" "union tissue maps"

    MISS_DEP=0
    [[ -f "${OUT_MNI_5TT}" ]] || { echo "    MISSING: ${OUT_MNI_5TT##*/}" >&2; MISS_DEP=1; }
    [[ -f "${FP_GM}"       ]] || { echo "    MISSING: ${FP_GM##*/}"       >&2; MISS_DEP=1; }
    [[ -f "${FP_WM}"       ]] || { echo "    MISSING: ${FP_WM##*/}"       >&2; MISS_DEP=1; }
    [[ -f "${FP_CSF}"      ]] || { echo "    MISSING: ${FP_CSF##*/}"      >&2; MISS_DEP=1; }

    if [[ "${MISS_DEP}" -eq 0 ]]; then
        if [[ ! -f "${OUT_UNION_GM}" || "${FORCE_DS}" -eq 1 ]]; then
            echo "    5TT MNI : ${OUT_MNI_5TT##*/}"
            echo "    fMRIPrep: ${FP_GM##*/} (and WM, CSF)"
            TMP_CTX="${OUT_DIR}/tmp_cortical_gm.nii.gz"
            TMP_BGT="${OUT_DIR}/tmp_bgt.nii.gz"
            TMP_WM5="${OUT_DIR}/tmp_wm5tt.nii.gz"
            TMP_CSF5="${OUT_DIR}/tmp_csf5tt.nii.gz"
            if {
                mrconvert "${OUT_MNI_5TT}" -coord 3 0 -axes 0,1,2 "${TMP_CTX}"  -force &&
                mrconvert "${OUT_MNI_5TT}" -coord 3 1 -axes 0,1,2 "${TMP_BGT}"  -force &&
                mrconvert "${OUT_MNI_5TT}" -coord 3 2 -axes 0,1,2 "${TMP_WM5}"  -force &&
                mrconvert "${OUT_MNI_5TT}" -coord 3 3 -axes 0,1,2 "${TMP_CSF5}" -force &&
                mrcalc "${TMP_CTX}" "${TMP_BGT}" -add 0 -gt \
                       "${FP_GM}"  0.25 -gt -add 0 -gt "${OUT_UNION_GM}"  -force &&
                mrcalc "${TMP_WM5}"  0 -gt \
                       "${FP_WM}"  0.25 -gt -add 0 -gt "${OUT_UNION_WM}"  -force &&
                mrcalc "${TMP_CSF5}" 0 -gt \
                       "${FP_CSF}" 0.25 -gt -add 0 -gt "${OUT_UNION_CSF}" -force
                #rm -f "${TMP_CTX}" "${TMP_BGT}" "${TMP_WM5}" "${TMP_CSF5}"
            } >> "${LOG_UNION}" 2>&1; then
                echo "    DONE ($(ts))"
                img_info "${OUT_UNION_GM}"
            else
                echo "    FAILED — see ${LOG_UNION}" >&2
                log_tail "${LOG_UNION}"
                # rm -f "${TMP_CTX}" "${TMP_BGT}" "${TMP_WM5}" "${TMP_CSF5}"
                (( N_FAIL++ )) || true
            fi
        else
            echo "    SKIP: output exists (use --force or --force2 to reprocess)"
            img_info "${OUT_UNION_GM}"
        fi
    else
        echo "    SKIP: missing dependencies (see above)" >&2
    fi

    # -----------------------------------------------------------------------
    # Step 4: Subject-specific Yeo17 atlas (MNI space)
    # union_GM is already in MNI — propagate Yeo17 labels through it directly.
    # 4a resamples the atlas onto THIS subject's own union_GM grid (must be
    # per-subject, not shared -- union_GM's grid varies by subject); 4b then
    # propagates through that same subject's mask, so both inputs to
    # PropagateLabelsThroughMask are guaranteed to be on one matching grid.
    # -----------------------------------------------------------------------
    step_header "4/5" "subject-specific Yeo17 (MNI)"

    if [[ ! -f "${OUT_UNION_GM}" ]]; then
        echo "    SKIP: ${OUT_UNION_GM##*/} not found" >&2
    else
        # [4a] Resample Yeo17 onto this subject's union_GM grid
        if [[ ! -f "${YEOBUCKNER17RS}" || "${FORCE_DS}" -eq 1 ]]; then
            echo "  [4a] Resample Yeo17 → 2mm MNI"
            echo "    input : ${YEOBUCKNER17##*/}"
            echo "    ref   : ${OUT_UNION_GM##*/}"
            echo "    output: ${YEOBUCKNER17RS##*/}"
            if antsApplyTransforms \
                    -d 3 \
                    -i "${YEOBUCKNER17}" \
                    -r "${OUT_UNION_GM}" \
                    -n GenericLabel \
                    -o "${YEOBUCKNER17RS}" \
                    >> "${LOG_YEO}" 2>&1; then
                echo "    DONE ($(ts))"
                img_info "${YEOBUCKNER17RS}"
            else
                echo "    FAILED — see ${LOG_YEO}" >&2
                log_tail "${LOG_YEO}"
                (( N_FAIL++ )) || true
            fi
        else
            echo "  [4a] Yeo17 2mm SKIP: output exists"
            img_info "${YEOBUCKNER17RS}"
        fi

        # [4b] Propagate labels through subject GM mask
        if [[ ! -f "${YEOBUCKNER17RS}" ]]; then
            echo "  [4b] SKIP: resampled Yeo17 atlas not available" >&2
        elif [[ ! -f "${OUT_YEO_SUBJ_MNI}" || "${FORCE_DS}" -eq 1 ]]; then
            echo "  [4b] PropagateLabelsThroughMask"
            echo "    mask  : ${OUT_UNION_GM##*/}"
            echo "    labels: ${YEOBUCKNER17RS##*/}"
            echo "    output: ${OUT_YEO_SUBJ_MNI##*/}"
            if ImageMath 3 "${OUT_YEO_SUBJ_MNI}" PropagateLabelsThroughMask \
                    "${OUT_UNION_GM}" "${YEOBUCKNER17RS}" 2 \
                    >> "${LOG_YEO}" 2>&1; then
                echo "    DONE ($(ts))"
                img_info "${OUT_YEO_SUBJ_MNI}"
            else
                echo "    FAILED — see ${LOG_YEO}" >&2
                log_tail "${LOG_YEO}"
                (( N_FAIL++ )) || true
            fi
        else
            echo "  [4b] Yeo17 propagation SKIP: output exists"
            img_info "${OUT_YEO_SUBJ_MNI}"
        fi
    fi

    # -----------------------------------------------------------------------
    # Step 5: Subject-specific SBA atlas propagation (ho_cort, fsl_striatum)
    # Propagates each atlas through the subject's union GM mask so that
    # SBA seeds are constrained to the individual's grey matter.
    # 5a/5c: resample atlases to MNI152NLin2009cAsym 2mm grid (once, shared).
    # 5b/5d: PropagateLabelsThroughMask per subject.
    # -----------------------------------------------------------------------
    step_header "5/5" "subject-specific SBA atlas propagation"

    if [[ ! -f "${OUT_UNION_GM}" ]]; then
        echo "    SKIP: ${OUT_UNION_GM##*/} not found — run Step 3 first" >&2
    else
        # ---- ho_cort --------------------------------------------------------
        if [[ ! -f "${HO_CORT_ATLAS}" ]]; then
            echo "    SKIP ho_cort: atlas not found at ${HO_CORT_ATLAS}" >&2
        else
            # [5a] Resample to MNI152NLin2009cAsym 2mm — shared, created once
            if [[ ! -f "${HO_CORT_RS}" || "${FORCE_DS}" -eq 1 ]]; then
                echo "  [5a] Resample ho_cort → 2mm MNI"
                if antsApplyTransforms \
                        -d 3 \
                        -i "${HO_CORT_ATLAS}" \
                        -r "${OUT_UNION_GM}" \
                        -n GenericLabel \
                        -o "${HO_CORT_RS}" \
                        >> "${LOG_ATLAS}" 2>&1; then
                    echo "    DONE ($(ts))"; img_info "${HO_CORT_RS}"
                else
                    echo "    FAILED — see ${LOG_ATLAS}" >&2; log_tail "${LOG_ATLAS}"
                fi
            else
                echo "  [5a] ho_cort 2mm SKIP: output exists"; img_info "${HO_CORT_RS}"
            fi

            # [5b] Propagate through subject GM
            if [[ ! -f "${HO_CORT_RS}" ]]; then
                echo "  [5b] SKIP: resampled ho_cort not available" >&2
            elif [[ ! -f "${OUT_HO_CORT_SUBJ}" || "${FORCE_DS}" -eq 1 ]]; then
                echo "  [5b] PropagateLabelsThroughMask: ho_cort"
                if ImageMath 3 "${OUT_HO_CORT_SUBJ}" PropagateLabelsThroughMask \
                        "${OUT_UNION_GM}" "${HO_CORT_RS}" 2 \
                        >> "${LOG_ATLAS}" 2>&1; then
                    echo "    DONE ($(ts))"; img_info "${OUT_HO_CORT_SUBJ}"
                else
                    echo "    FAILED — see ${LOG_ATLAS}" >&2; log_tail "${LOG_ATLAS}"
                    (( N_FAIL++ )) || true
                fi
            else
                echo "  [5b] ho_cort propagation SKIP: output exists"
                img_info "${OUT_HO_CORT_SUBJ}"
            fi
        fi

        # ---- fsl_striatum ---------------------------------------------------
        if [[ ! -f "${FSL_STRIATUM_ATLAS}" ]]; then
            echo "    SKIP fsl_striatum: atlas not found at ${FSL_STRIATUM_ATLAS}" >&2
        else
            # [5c] Resample to MNI152NLin2009cAsym 2mm — shared, created once
            if [[ ! -f "${FSL_STRIATUM_RS}" || "${FORCE_DS}" -eq 1 ]]; then
                echo "  [5c] Resample fsl_striatum → 2mm MNI"
                if antsApplyTransforms \
                        -d 3 \
                        -i "${FSL_STRIATUM_ATLAS}" \
                        -r "${OUT_UNION_GM}" \
                        -n GenericLabel \
                        -o "${FSL_STRIATUM_RS}" \
                        >> "${LOG_ATLAS}" 2>&1; then
                    echo "    DONE ($(ts))"; img_info "${FSL_STRIATUM_RS}"
                else
                    echo "    FAILED — see ${LOG_ATLAS}" >&2; log_tail "${LOG_ATLAS}"
                fi
            else
                echo "  [5c] fsl_striatum 2mm SKIP: output exists"; img_info "${FSL_STRIATUM_RS}"
            fi

            # [5d] Propagate through subject GM
            if [[ ! -f "${FSL_STRIATUM_RS}" ]]; then
                echo "  [5d] SKIP: resampled fsl_striatum not available" >&2
            elif [[ ! -f "${OUT_FSL_STRIATUM_SUBJ}" || "${FORCE_DS}" -eq 1 ]]; then
                echo "  [5d] PropagateLabelsThroughMask: fsl_striatum"
                if ImageMath 3 "${OUT_FSL_STRIATUM_SUBJ}" PropagateLabelsThroughMask \
                        "${OUT_UNION_GM}" "${FSL_STRIATUM_RS}" 2 \
                        >> "${LOG_ATLAS}" 2>&1; then
                    echo "    DONE ($(ts))"; img_info "${OUT_FSL_STRIATUM_SUBJ}"
                else
                    echo "    FAILED — see ${LOG_ATLAS}" >&2; log_tail "${LOG_ATLAS}"
                    (( N_FAIL++ )) || true
                fi
            else
                echo "  [5d] fsl_striatum propagation SKIP: output exists"
                img_info "${OUT_FSL_STRIATUM_SUBJ}"
            fi
        fi
    fi

    # -----------------------------------------------------------------------
    # Step 6: Lausanne2018 scale3 subject-specific atlas (hand/foot/lip SBA
    # seeds -- KUL_FS_multiparc.sh's output, FreeSurfer-conformed space).
    # Two-stage warp, same pattern already used for FS-derived label volumes
    # elsewhere (KUL_dwiprep_anat.sh's aparc+aseg warp): mri_convert with
    # -rl/-rt nearest handles the FS-conformed -> T1w-native resample (header-
    # based, no separate registration needed), then warp_to_mni's existing
    # GenericLabel antsApplyTransforms takes T1w-native -> MNI, same as every
    # other subject-specific atlas in this script.
    # -----------------------------------------------------------------------
    step_header "6/6" "Lausanne2018 scale3 subject-specific atlas"

    if [[ -f "${LAUSANNE_VBG}" ]]; then
        LAUSANNE_SRC="${LAUSANNE_VBG}"
    elif [[ -f "${LAUSANNE_STD}" ]]; then
        LAUSANNE_SRC="${LAUSANNE_STD}"
    else
        LAUSANNE_SRC=""
    fi

    if [[ -z "${LAUSANNE_SRC}" ]]; then
        echo "    SKIP: lausanne2018.scale3+aseg.mgz not found (checked VBG and standard freesurfer derivatives) -- run KUL_FS_multiparc.sh first" >&2
    else
        if [[ ! -f "${OUT_LAUSANNE_T1W}" || "${FORCE_DS}" -eq 1 ]]; then
            echo "  [6a] mri_convert: FS-conformed -> T1w-native"
            if mri_convert -rl "${T1W}" -rt nearest "${LAUSANNE_SRC}" "${OUT_LAUSANNE_T1W}" \
                    > "${LOG_LAUSANNE}" 2>&1; then
                echo "    DONE ($(ts))"; img_info "${OUT_LAUSANNE_T1W}"
            else
                echo "    FAILED — see ${LOG_LAUSANNE}" >&2; log_tail "${LOG_LAUSANNE}"
                (( N_FAIL++ )) || true
            fi
        else
            echo "  [6a] T1w-native SKIP: output exists"; img_info "${OUT_LAUSANNE_T1W}"
        fi

        if [[ ! -f "${OUT_LAUSANNE_T1W}" ]]; then
            echo "  [6b] SKIP: T1w-native lausanne not available" >&2
        elif [[ ! -f "${OUT_LAUSANNE_MNI}" || "${FORCE_DS}" -eq 1 ]]; then
            echo "  [6b] warp T1w-native -> MNI"
            warp_to_mni "${OUT_LAUSANNE_T1W}" "${OUT_LAUSANNE_MNI}" "${XFM}" "${LOG_LAUSANNE}" \
                || (( N_FAIL++ )) || true
        else
            echo "  [6b] MNI lausanne SKIP: output exists"; img_info "${OUT_LAUSANNE_MNI}"
        fi
    fi

done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Done  [$(ts)]"
echo "  SynthSeg: ${N_DONE} completed  |  ${N_SKIP} skipped  |  ${N_FAIL} failed"
echo "  Outputs : ${SYNTHSEG_DIR}"
