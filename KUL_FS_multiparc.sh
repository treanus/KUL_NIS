#!/bin/bash
# KUL_FS_multiparc.sh
# Ahmed M. Radwan, KU Leuven Translational MRI
#
# FreeSurfer recon-all + multi-scale parcellation for non-lesioned subjects
# (types 4, 5, 6 in KUL_clinical_fmridti.sh — no VBG needed).
#
# If a completed recon-all is already present under <fs_subjects_dir>/<subject_id>
# (detected by scripts/recon-all.done) the recon step is skipped and only the
# parcellations are (re-)run.
#
# Parcellations added (same set as KUL_VBG_latest/KUL_VBG_multiparc.sh):
#   - Lausanne2018 scales 1-5
#   - Glasser HCP-MMP1
#   - Thalamic nuclei          (segment_subregions thalamus)
#   - Brainstem substructures  (segment_subregions brainstem)
#   - Hippocampus / amygdala   (segment_subregions hippo-amygdala)
#
# Usage:
#   KUL_FS_multiparc.sh -s <subject_id> -f <fs_subjects_dir> -i <T1w.nii.gz>
#                       [-n <threads>] [-X] [-v] [-h]
#
# Example:
#   KUL_FS_multiparc.sh -s sub-Patient01 \
#       -f /data/study/BIDS/derivatives/freesurfer \
#       -i /data/study/BIDS/sub-Patient01/anat/sub-Patient01_T1w.nii.gz \
#       -n 16

set -e

version="0.2 — 2026-07-07"

# ── Script location — atlases live in sibling KUL_VBG_latest repo ─────────────
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
atlases_dir="${script_dir}/atlases"
lausanne_dir="${atlases_dir}/lausanne2008"
glasser_dir="${atlases_dir}/glasser"
remap_py="${lausanne_dir}/remap_lausanne_to_msbp.py"

# ── Defaults ──────────────────────────────────────────────────────────────────
ncpu=6
verbose=0
use_fastsurfer=0
subj=""
fs_dir=""
t1w_in=""

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
KUL_FS_multiparc.sh v${version}

Usage: $(basename "$0") -s <subject_id> -f <fs_subjects_dir> -i <T1w.nii.gz> [options]

Required:
  -s <subject_id>       Subject ID — must match the folder name to be created
                        under <fs_subjects_dir> (e.g. sub-Patient01)
  -f <fs_subjects_dir>  Path to the FreeSurfer subjects directory.
                        If <fs_subjects_dir>/<subject_id>/scripts/recon-all.done
                        does not exist, recon-all is run automatically.
  -i <T1w.nii.gz>       T1w image used as input for recon-all (only needed when
                        recon-all has not yet been run)

Options:
  -n <n>   Number of threads (default: ${ncpu})
  -X       Use FastSurfer instead of plain recon-all for the reconstruction step
  -v       Verbose output
  -h       Show this help

Outputs (written to <fs_subjects_dir>/<subject_id>/mri/):
  aparc+aseg.mgz                        (from recon-all)
  lausanne2018.scale1-5+aseg.mgz        Lausanne2018 cortical parcellation
  HCPMMP1+aseg.mgz                      Glasser HCP-MMP1
  ThalamicNuclei.v12.T1.FSvoxelSpace.mgz
  brainstemSsLabels.v12.FSvoxelSpace.mgz

EOF
    exit 0
}

# ── Argument parsing ──────────────────────────────────────────────────────────
while getopts ":s:f:i:n:Xvh" opt; do
    case $opt in
        s) subj="$OPTARG" ;;
        f) fs_dir="$OPTARG" ;;
        i) t1w_in="$OPTARG" ;;
        n) ncpu="$OPTARG" ;;
        X) use_fastsurfer=1 ;;
        v) verbose=1 ;;
        h) usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
       \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# ── Preflight ─────────────────────────────────────────────────────────────────
[[ -z "$subj" ]]   && { echo "ERROR: -s <subject_id> required"; exit 1; }
[[ -z "$fs_dir" ]] && { echo "ERROR: -f <fs_subjects_dir> required"; exit 1; }

subj_dir="${fs_dir}/${subj}"
scripts_dir="${subj_dir}/scripts"
mri_dir="${subj_dir}/mri"

# ── Clean up broken stub BEFORE any mkdir (must happen before creating subj_dir,
#    otherwise recon-all refuses -i because it sees an existing subject folder) ──
if [[ -d "${subj_dir}" && ! -f "${scripts_dir}/recon-all.done" && ! -f "${mri_dir}/orig.mgz" ]]; then
    echo "Removing incomplete subject dir (no orig.mgz) and starting fresh..."
    rm -rf "${subj_dir}"
fi

for atlas_file in \
    "${lausanne_dir}/lh.lausanne2018.scale1.annot" \
    "${glasser_dir}/lh.HCPMMP1.annot" \
    "${remap_py}"; do
    [[ ! -f "$atlas_file" ]] && {
        echo "ERROR: atlas file not found: ${atlas_file}"
        echo "       Is KUL_VBG_latest present as a sibling of KUL_NIS_unified?"
        exit 1
    }
done

for py_dep in nibabel numpy; do
    python3 -c "import ${py_dep}" 2>/dev/null || {
        echo "ERROR: Python package '${py_dep}' not found (needed for MSBP remapping)"
        exit 1
    }
done

# ── FreeSurfer environment ────────────────────────────────────────────────────
[[ -z "$FREESURFER_HOME" ]] && { echo "ERROR: FREESURFER_HOME not set — source FreeSurfer setup first"; exit 1; }
echo "Using FreeSurfer at ${FREESURFER_HOME}"

# Log lives in the parent freesurfer/ dir (not inside the subject dir) so that
# recon-all does not see an existing subject folder on a fresh run.
mkdir -p "${fs_dir}"
log_file="${fs_dir}/KUL_FS_multiparc_${subj}_$(date +%Y%m%d_%H%M%S).log"

log() { echo "$*" | tee -a "${log_file}"; }
run() {
    log "  CMD: $*"
    if [[ $verbose -eq 1 ]]; then
        eval "$*" 2>&1 | tee -a "${log_file}"
    else
        eval "$*" >> "${log_file}" 2>&1
    fi
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        log "ERROR: command failed — see log: ${log_file}"
        exit 1
    fi
}

log "========================================================"
log " KUL_FS_multiparc.sh v${version}"
log " Subject:    ${subj}"
log " FS dir:     ${fs_dir}"
log " Threads:    ${ncpu}"
log " Engine:     $([ $use_fastsurfer -eq 1 ] && echo FastSurfer || echo recon-all)"
log " Started:    $(date)"
log "========================================================"

export OMP_NUM_THREADS=${ncpu}
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ncpu}

# ── STEP 1: recon-all (skip if already done) ──────────────────────────────────
if [[ -f "${scripts_dir}/recon-all.done" ]]; then
    log "recon-all.done found — skipping reconstruction"
else
    [[ -z "$t1w_in" ]] && { echo "ERROR: -i <T1w.nii.gz> required — no recon-all.done found in ${scripts_dir}"; exit 1; }
    [[ ! -f "$t1w_in" ]] && { echo "ERROR: T1w image not found: ${t1w_in}"; exit 1; }

    if [[ $use_fastsurfer -eq 1 ]]; then
        # ── FastSurfer path ───────────────────────────────────────────────────
        log ""
        log "Running FastSurfer + recon-all -make-all..."
        fasu_output="${fs_dir}/../FastSurfer"
        mkdir -p "${fasu_output}"

        nvram=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1 || echo "0")
        if [[ ${nvram:-0} -lt 5500 ]]; then
            FaSu_cpu="--no_cuda"
        else
            FaSu_cpu=""
        fi

        run "run_fastsurfer.sh --t1 ${t1w_in} \
            --sid ${subj} --sd ${fasu_output} --fsaparc --parallel --threads ${ncpu} \
            --fs_license ${FS_LICENSE} --py python ${FaSu_cpu} --ignore_fs_version"

        run "rsync -a ${fasu_output}/${subj}/ ${subj_dir}/"

        run "recon-all -s ${subj} -sd ${fs_dir} -openmp ${ncpu} -parallel -no-isrunning -make all"

    else
        # ── Plain recon-all (default) ─────────────────────────────────────────
        log ""
        # If orig.mgz exists a prior run got past the import stage — continue
        # without -i so recon-all picks up where it left off.
        # Otherwise the stub was already cleaned at startup; start fresh.
        if [[ -f "${mri_dir}/orig.mgz" ]]; then
            log "Partial recon-all found (orig.mgz present) — continuing..."
            run "recon-all -s ${subj} -sd ${fs_dir} \
                -openmp ${ncpu} -parallel -all -no-isrunning"
        else
            log "Running recon-all (FreeSurfer)..."
            run "recon-all -i ${t1w_in} -s ${subj} -sd ${fs_dir} \
                -openmp ${ncpu} -parallel -all -no-isrunning"
        fi
    fi
fi

# ── STEP 2: multi-scale parcellation ─────────────────────────────────────────
if [[ -f "${scripts_dir}/multiscale_parc.done" ]]; then
    log "multiscale_parc.done exists — parcellation already complete."
    log "Delete ${scripts_dir}/multiscale_parc.done to force rerun."
    exit 0
fi

[[ ! -d "${fs_dir}/fsaverage" ]] && \
    ln -sf "${FREESURFER_HOME}/subjects/fsaverage" "${fs_dir}/fsaverage"

# Lausanne2018 scales 1-5
log ""
log "Running Lausanne2018 parcellation (scales 1-5)..."
for scale in 1 2 3 4 5; do
    log "  Scale ${scale}..."
    for hemi in lh rh; do
        run "mri_surf2surf \
            --srcsubject fsaverage \
            --trgsubject ${subj} \
            --hemi ${hemi} \
            --sval-annot ${lausanne_dir}/${hemi}.lausanne2018.scale${scale}.annot \
            --tval ${subj_dir}/label/${hemi}.lausanne2018.scale${scale}.annot \
            --sd ${fs_dir}"
    done

    _raw="${mri_dir}/lausanne2018.scale${scale}+aseg_raw.mgz"
    _out="${mri_dir}/lausanne2018.scale${scale}+aseg.mgz"
    _lut="${lausanne_dir}/label-L2018_desc-scale${scale}_atlas_FreeSurferColorLUT.txt"

    run "mri_aparc2aseg \
        --s ${subj} \
        --sd ${fs_dir} \
        --annot lausanne2018.scale${scale} \
        --o ${_raw}"

    if [[ -f "${_lut}" ]]; then
        run "python3 ${remap_py} \
            --input   ${_raw} \
            --lh_annot ${lausanne_dir}/lh.lausanne2018.scale${scale}.annot \
            --rh_annot ${lausanne_dir}/rh.lausanne2018.scale${scale}.annot \
            --lut     ${_lut} \
            --output  ${_out}"
        rm -f "${_raw}"
    else
        log "  WARNING: no MSBP LUT for scale ${scale} — IDs not remapped"
        mv "${_raw}" "${_out}"
    fi
done

# Glasser HCP-MMP1
log ""
log "Running Glasser HCP-MMP1 parcellation..."
for hemi in lh rh; do
    run "mri_surf2surf \
        --srcsubject fsaverage \
        --trgsubject ${subj} \
        --hemi ${hemi} \
        --sval-annot ${glasser_dir}/${hemi}.HCPMMP1.annot \
        --tval ${subj_dir}/label/${hemi}.HCPMMP1.annot \
        --sd ${fs_dir}"
done
run "mri_aparc2aseg \
    --s ${subj} \
    --sd ${fs_dir} \
    --annot HCPMMP1 \
    --o ${mri_dir}/HCPMMP1+aseg.mgz"

# Thalamic nuclei
log ""
log "Running thalamic nuclei segmentation..."
run "segment_subregions thalamus --cross ${subj} --sd ${fs_dir} --threads ${ncpu}"

# Brainstem substructures
log ""
log "Running brainstem substructure segmentation..."
run "segment_subregions brainstem --cross ${subj} --sd ${fs_dir} --threads ${ncpu}"

# Hippocampus / amygdala
log ""
log "Running hippocampus-amygdala segmentation..."
run "segment_subregions hippo-amygdala --cross ${subj} --sd ${fs_dir} --threads ${ncpu}"

# Hypothalamic subunits
log ""
log "Running hypothalamic subunit segmentation..."
run "mri_segment_hypothalamic_subunits --s ${subj} --sd ${fs_dir} --threads ${ncpu}"

# ── Done ──────────────────────────────────────────────────────────────────────
touch "${scripts_dir}/multiscale_parc.done"

log ""
log "========================================================"
log " KUL_FS_multiparc.sh complete."
log " Outputs in ${mri_dir}/"
log " Finished: $(date)"
log "========================================================"
