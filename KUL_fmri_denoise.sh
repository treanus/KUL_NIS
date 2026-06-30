#!/usr/bin/env bash
# KUL_fmri_denoise.sh — Post-fMRIPrep denoising
# Confound regression + bandpass filter, then SUSAN spatial smoothing
#
# --method nilearn  simultaneous regression + bandpass via nilearn.image.clean_img
#                   (methodologically cleaner: single GLM pass)
#                   requires: pip install nilearn
#
# --method fsl      sequential fsl_glm (regression) then fslmaths -bptf (bandpass)
#                   requires only FSL + nibabel (already installed)

set -euo pipefail

# ── Defaults ───────────────────────────────────────────────────────────────────
FMRIPREP_DIR=""
OUT_DIR=""
SPACE="MNI152NLin2009cAsym"
METHOD="nilearn"

SUB_FILTER=""
SES_FILTER=""
TASK_FILTER=""
RUN_FILTER=""
DRYRUN=0

DO_SUSAN=1
FWHM=6
SUSAN_DIM=3
USE_MEDIAN=1
N_USANS=1

HP=0.008
LP=0.09

COLS_STR=""
COLS_FILE=""
ALL_CONFOUNDS=0

# ── Usage ──────────────────────────────────────────────────────────────────────
usage() {
  cat <<EOF
Usage: $(basename "$0") --fmriprep <dir> --out <dir> [options]

Core:
  --method <nilearn|fsl>    Denoising backend (default: nilearn)
  --space <SPACE>           fMRIPrep output space (default: MNI152NLin2009cAsym)
  --sub <sub-XX>
  --ses <01|ses-01>
  --task <taskname>
  --run <01|run-01>
  --dry-run

SUSAN:
  --no-susan
  --fwhm <mm>               (default: 6)
  --susan-dim <2|3>         (default: 3)

Denoising:
  --hp <Hz>                 High-pass cutoff (default: 0.008)
  --lp <Hz>                 Low-pass cutoff  (default: 0.09)
  --all-confounds           Use 24 motion regressors (adds squares + squared derivatives
                            to the default 12); still combined with aCompCor
  --cols "<c1 c2 ...>"      Fully override confound columns
  --cols-file <file.txt>    Read confound columns from file (one per line, # comments ok)

Methods:
  nilearn   Simultaneous confound regression + bandpass via nilearn.image.clean_img.
            Temporal mean restored after cleaning. SUSAN follows denoising.
            Install: pip install nilearn

  fsl       fsl_glm confound regression, then fslmaths -bptf bandpass.
            Temporal mean restored after filtering. SUSAN follows denoising.
            Requires only FSL + nibabel (no extra install needed).
EOF
}

# ── Helpers ────────────────────────────────────────────────────────────────────
norm_ses() { [[ -z "${1:-}" ]] && echo "" && return; [[ "$1" == ses-* ]] && echo "$1" || echo "ses-$1"; }
norm_run() { [[ -z "${1:-}" ]] && echo "" && return; [[ "$1" == run-* ]] && echo "$1" || echo "run-$1"; }

run_cmd() {
  echo "+ $*"
  if [[ "$DRYRUN" -eq 0 ]]; then
    eval "$@"
  fi
}

conf_from_bold() {
  local base="$1"
  local stripped
  stripped="$(echo "$base" | sed -E 's/_space-[^_]+//g; s/_res-[^_]+//g')"
  echo "${stripped/_desc-preproc_bold.nii.gz/_desc-confounds_timeseries.tsv}"
}

get_outlier_cols() {
  python3 - "$1" <<'PY'
import sys
hdr = open(sys.argv[1], encoding="utf-8", errors="ignore").readline().rstrip("\n")
cols = hdr.split("\t")
keep = [c for c in cols if c.startswith(("motion_outlier", "non_steady_state_outlier", "outlier"))]
print(" ".join(keep))
PY
}

find_boldref() {
  local dir="$1" bold_base="$2" outdir="$3"
  local prefix="${bold_base%_desc-preproc_bold.nii.gz}"
  for cand in \
      "${dir}/${prefix}_boldref.nii.gz" \
      "${dir}/${prefix}_desc-coreg_boldref.nii.gz" \
      "${dir}/${prefix}_desc-boldref.nii.gz"; do
    [[ -f "$cand" ]] && echo "$cand" && return 0
  done
  local fallback="${outdir}/${prefix}_MEAN_boldref.nii.gz"
  [[ -f "$fallback" ]] || run_cmd "fslmaths \"${dir}/${bold_base}\" -Tmean \"$fallback\""
  echo "$fallback"
}

get_tr() {
  python3 - "$1" <<'PY'
import sys, json
try:
    data = json.load(open(sys.argv[1]))
    print(data["RepetitionTime"])
except Exception as e:
    print(f"ERROR reading TR: {e}", file=sys.stderr)
    sys.exit(1)
PY
}

# ── Argument parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fmriprep)   FMRIPREP_DIR="$2"; shift 2;;
    --out)        OUT_DIR="$2"; shift 2;;
    --space)      SPACE="$2"; shift 2;;
    --method)     METHOD="$2"; shift 2;;
    --sub)        SUB_FILTER="$2"; shift 2;;
    --ses)        SES_FILTER="$(norm_ses "$2")"; shift 2;;
    --task)       TASK_FILTER="$2"; shift 2;;
    --run)        RUN_FILTER="$(norm_run "$2")"; shift 2;;
    --dry-run)    DRYRUN=1; shift;;
    --no-susan)   DO_SUSAN=0; shift;;
    --fwhm)       FWHM="$2"; shift 2;;
    --susan-dim)  SUSAN_DIM="$2"; shift 2;;
    --hp)         HP="$2"; shift 2;;
    --lp)         LP="$2"; shift 2;;
    --all-confounds) ALL_CONFOUNDS=1; shift;;
    --cols)       COLS_STR="$2"; shift 2;;
    --cols-file)  COLS_FILE="$2"; shift 2;;
    -h|--help)    usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -n "$FMRIPREP_DIR" && -n "$OUT_DIR" ]] || { usage; exit 1; }
[[ "$METHOD" == "nilearn" || "$METHOD" == "fsl" ]] || {
  echo "ERROR: --method must be 'nilearn' or 'fsl'"; exit 1
}

# ── Dependency checks ──────────────────────────────────────────────────────────
for cmd in fslmaths fslstats susan; do
  command -v "$cmd" >/dev/null || { echo "ERROR: $cmd not in PATH"; exit 1; }
done
command -v python3 >/dev/null || { echo "ERROR: python3 not in PATH"; exit 1; }
python3 -c "import nibabel" 2>/dev/null || { echo "ERROR: nibabel not installed (pip install nibabel)"; exit 1; }

if [[ "$METHOD" == "nilearn" ]]; then
  python3 -c "import nilearn" 2>/dev/null || {
    echo "ERROR: nilearn not installed. Run:  pip install nilearn"
    exit 1
  }
elif [[ "$METHOD" == "fsl" ]]; then
  command -v fsl_glm >/dev/null || { echo "ERROR: fsl_glm not in PATH"; exit 1; }
fi

# ── Confound columns ───────────────────────────────────────────────────────────
# Default: 12 motion parameters (6 + derivatives) + 5 WM aCompCor + 5 CSF aCompCor
# --all-confounds expands to 24 motion params (adds squares + squared derivatives)
# Outlier spike regressors (motion_outlier*, non_steady_state_outlier*) are
# added automatically per-run from the confounds TSV header.
DEFAULT_COLS=(
  trans_x trans_y trans_z rot_x rot_y rot_z
  trans_x_derivative1 trans_y_derivative1 trans_z_derivative1
  rot_x_derivative1 rot_y_derivative1 rot_z_derivative1
  w_comp_cor_00 w_comp_cor_01 w_comp_cor_02 w_comp_cor_03 w_comp_cor_04
  c_comp_cor_00 c_comp_cor_01 c_comp_cor_02 c_comp_cor_03 c_comp_cor_04
)

ALL_COLS=(
  trans_x trans_y trans_z rot_x rot_y rot_z
  trans_x_derivative1 trans_y_derivative1 trans_z_derivative1
  rot_x_derivative1 rot_y_derivative1 rot_z_derivative1
  trans_x_power2 trans_y_power2 trans_z_power2
  rot_x_power2 rot_y_power2 rot_z_power2
  trans_x_derivative1_power2 trans_y_derivative1_power2 trans_z_derivative1_power2
  rot_x_derivative1_power2 rot_y_derivative1_power2 rot_z_derivative1_power2
  w_comp_cor_00 w_comp_cor_01 w_comp_cor_02 w_comp_cor_03 w_comp_cor_04
  c_comp_cor_00 c_comp_cor_01 c_comp_cor_02 c_comp_cor_03 c_comp_cor_04
)

COLS=()
if [[ -n "$COLS_FILE" ]]; then
  [[ -f "$COLS_FILE" ]] || { echo "ERROR: --cols-file not found: $COLS_FILE"; exit 1; }
  while IFS= read -r line; do
    line="${line%%#*}"
    line="$(echo "$line" | awk '{$1=$1;print}')"
    [[ -z "$line" ]] && continue
    COLS+=("$line")
  done < "$COLS_FILE"
elif [[ -n "$COLS_STR" ]]; then
  read -r -a COLS <<< "$COLS_STR"
elif [[ "$ALL_CONFOUNDS" -eq 1 ]]; then
  COLS=("${ALL_COLS[@]}")
else
  COLS=("${DEFAULT_COLS[@]}")
fi
[[ "${#COLS[@]}" -gt 0 ]] || { echo "ERROR: no confound columns selected"; exit 1; }

DT=$(echo "$FWHM/2.355" | bc -l)   # SUSAN spatial sigma (mm)

# ── Find BOLD files ────────────────────────────────────────────────────────────
mkdir -p "$OUT_DIR"

mapfile -t BOLDS < <(
  find "$FMRIPREP_DIR" \
    -path "$FMRIPREP_DIR/derivatives" -prune -o \
    -type f -path "*/func/*" \
    -name "*_space-${SPACE}_*_desc-preproc_bold.nii.gz" -print \
  | sort
)

[[ "${#BOLDS[@]}" -gt 0 ]] || {
  echo "No BOLD files found for space=${SPACE} in $FMRIPREP_DIR"
  exit 1
}
echo "Found ${#BOLDS[@]} BOLD run(s)  [method=${METHOD}  space=${SPACE}]"

# ── Per-run loop ───────────────────────────────────────────────────────────────
for bold in "${BOLDS[@]}"; do
  base="$(basename "$bold")"
  dir="$(dirname "$bold")"

  [[ -n "$SUB_FILTER"  && "$bold" != *"/${SUB_FILTER}/"*  ]] && continue
  [[ -n "$SES_FILTER"  && "$bold" != *"/${SES_FILTER}/"*  ]] && continue
  [[ -n "$TASK_FILTER" && "$base" != *"_task-${TASK_FILTER}_"* ]] && continue
  [[ -n "$RUN_FILTER"  && "$base" != *"_${RUN_FILTER}_"*  ]] && continue

  conf="${dir}/$(conf_from_bold "$base")"
  mask="${dir}/${base/_desc-preproc_bold.nii.gz/_desc-brain_mask.nii.gz}"
  json_sidecar="${dir}/${base/_desc-preproc_bold.nii.gz/_desc-preproc_bold.json}"

  if [[ ! -f "$conf" ]]; then
    echo "WARN: missing confounds TSV — skip"; echo "      expected: $conf"; continue
  fi
  if [[ ! -f "$mask" ]]; then
    echo "WARN: missing brain mask — skip"; echo "      expected: $mask"; continue
  fi
  if [[ ! -f "$json_sidecar" ]]; then
    echo "WARN: missing JSON sidecar — skip"; echo "      expected: $json_sidecar"; continue
  fi

  # output mirrors fMRIPrep directory tree
  rel="${bold#${FMRIPREP_DIR}/}"
  out_func_dir="${OUT_DIR}/$(dirname "$rel")"
  mkdir -p "$out_func_dir"

  run_label="${base/_desc-preproc_bold.nii.gz/}"
  run_out="${out_func_dir}/${run_label}_postproc_${METHOD}"
  mkdir -p "$run_out"

  echo "────────────────────────────────────────"
  echo "BOLD : $bold"
  echo "CONF : $conf"
  echo "MASK : $mask"
  echo "OUT  : $run_out"

  TR="$(get_tr "$json_sidecar")"
  echo "TR   : ${TR}s"

  # Auto-detect spike regressors from TSV header
  OUTLIERS_STR="$(get_outlier_cols "$conf")"
  RUN_COLS=("${COLS[@]}")
  OUTLIER_COUNT=0
  if [[ -n "$OUTLIERS_STR" ]]; then
    read -r -a OUTLIERS <<< "$OUTLIERS_STR"
    RUN_COLS+=("${OUTLIERS[@]}")
    OUTLIER_COUNT="${#OUTLIERS[@]}"
  fi
  echo "Regressors: ${#RUN_COLS[@]} (${OUTLIER_COUNT} spike columns)"

  # Final denoised output (both methods write here)
  denoised_out="${run_out}/${run_label}_desc-denoised_bold.nii.gz"
  if [[ -f "$denoised_out" ]]; then
    echo "[SKIP] Denoised output already exists: $denoised_out"
    continue
  fi

  tmean_out="${run_out}/${run_label}_Tmean.nii.gz"

  # If SUSAN follows denoising, write to intermediate; else write final directly
  if [[ "$DO_SUSAN" -eq 1 ]]; then
    denoise_target="${run_out}/${run_label}_desc-denoised_bold_presmooth.nii.gz"
  else
    denoise_target="$denoised_out"
  fi

  # ── Step 1: Denoise — nilearn ──────────────────────────────────────────────
  if [[ "$METHOD" == "nilearn" ]]; then

    if [[ "$DRYRUN" -eq 1 ]]; then
      echo "+ python3 [nilearn clean_img]"
      echo "  in:    $bold"
      echo "  tmean: $tmean_out"
      echo "  out:   $denoise_target"
      echo "  cols:  ${RUN_COLS[*]}"
      [[ "$DO_SUSAN" -eq 1 ]] && echo "+ susan [on denoised] → $denoised_out"
      continue
    fi

    python3 - \
      "$bold" "$conf" "$mask" "$denoise_target" "$tmean_out" \
      "$TR" "$HP" "$LP" \
      "${RUN_COLS[@]}" <<'PY'
import sys, csv, math
import numpy as np
import nibabel as nib
from nilearn.image import clean_img

bold_f   = sys.argv[1]
conf_f   = sys.argv[2]
mask_f   = sys.argv[3]
out_f    = sys.argv[4]
tmean_f  = sys.argv[5]
tr       = float(sys.argv[6])
hp       = float(sys.argv[7])
lp       = float(sys.argv[8])
cols     = sys.argv[9:]

with open(conf_f, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    rows = list(reader)

matrix = []
for row in rows:
    vec = []
    for c in cols:
        val = row.get(c, 'n/a')
        try:
            v = float(val)
            vec.append(0.0 if math.isnan(v) else v)
        except (ValueError, TypeError):
            vec.append(0.0)
    matrix.append(vec)

confounds_array = np.array(matrix, dtype=np.float64)
mask_img = nib.load(mask_f)

# Compute and save temporal mean before cleaning (clean_img removes it via detrending)
bold_img = nib.load(bold_f)
bold_data = bold_img.get_fdata(dtype=np.float32)
tmean_data = bold_data.mean(axis=-1)
nib.save(nib.Nifti1Image(tmean_data, bold_img.affine, bold_img.header), tmean_f)
print(f"  Tmean saved: {tmean_f}")

print(f"  {len(rows)} timepoints  |  {len(cols)} regressors  |  TR={tr}s  |  bandpass {hp}–{lp} Hz")

cleaned = clean_img(
    bold_f,
    confounds=confounds_array,
    t_r=tr,
    high_pass=hp,
    low_pass=lp,
    detrend=True,
    standardize=None,
    mask_img=mask_img,
)

# Restore temporal mean and re-apply mask
cleaned_data = cleaned.get_fdata(dtype=np.float32)
result_data = (cleaned_data + tmean_data[..., np.newaxis]).astype(np.float32)
mask_arr = nib.load(mask_f).get_fdata(dtype=np.float32) > 0
result_data *= mask_arr[..., np.newaxis]
nib.save(nib.Nifti1Image(result_data, cleaned.affine, cleaned.header), out_f)
print(f"  Saved: {out_f}")
PY

  # ── Step 1: Denoise — FSL ──────────────────────────────────────────────────
  elif [[ "$METHOD" == "fsl" ]]; then

    # Convert Hz to sigma in volumes for fslmaths -bptf
    hp_sigma=$(python3 -c "tr,f=float('$TR'),float('$HP'); print(f'{1/(2*tr*f):.4f}')")
    lp_sigma=$(python3 -c "tr,f=float('$TR'),float('$LP'); print(f'{1/(2*tr*f):.4f}')")
    echo "FSL bptf: hp_sigma=${hp_sigma} vols  lp_sigma=${lp_sigma} vols"

    conf_matrix="${run_out}/confounds_matrix.txt"
    residuals="${run_out}/${run_label}_residuals.nii.gz"
    bp_out="${run_out}/${run_label}_bp.nii.gz"

    if [[ "$DRYRUN" -eq 1 ]]; then
      echo "+ python3 [write confound matrix] → $conf_matrix"
      echo "+ fsl_glm --in=$bold --design=$conf_matrix --out_res=$residuals --demean"
      echo "+ fslmaths $bold -Tmean $tmean_out"
      echo "+ fslmaths $residuals -bptf $hp_sigma $lp_sigma $bp_out"
      echo "+ fslmaths $bp_out -add $tmean_out -mas $mask $denoise_target"
      [[ "$DO_SUSAN" -eq 1 ]] && echo "+ susan [on denoised] → $denoised_out"
      continue
    fi

    # Write confound matrix: intercept column + selected regressors, NaN → 0
    python3 - "$conf" "$conf_matrix" "${RUN_COLS[@]}" <<'PY'
import sys, csv, math
conf_f, out_f = sys.argv[1], sys.argv[2]
cols = sys.argv[3:]

with open(conf_f, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    rows = list(reader)

lines = []
for row in rows:
    vals = ['1']  # intercept
    for c in cols:
        val = row.get(c, 'n/a')
        try:
            v = float(val)
            vals.append('0' if math.isnan(v) else f'{v:.8f}')
        except (ValueError, TypeError):
            vals.append('0')
    lines.append(' '.join(vals))

with open(out_f, 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print(f"  Confound matrix: {len(rows)} rows × {len(cols)+1} cols (incl. intercept)")
PY

    # Confound regression (fsl_glm demeaned data → zero-mean residuals)
    run_cmd "fsl_glm --in=\"$bold\" --design=\"$conf_matrix\" --out_res=\"$residuals\" --demean"

    # Save temporal mean to restore after bandpass
    run_cmd "fslmaths \"$bold\" -Tmean \"$tmean_out\""

    # Bandpass filter
    run_cmd "fslmaths \"$residuals\" -bptf \"$hp_sigma\" \"$lp_sigma\" \"$bp_out\""

    # Restore mean, apply brain mask, write denoised output
    run_cmd "fslmaths \"$bp_out\" -add \"$tmean_out\" -mas \"$mask\" \"$denoise_target\""

    # Remove large intermediates
    rm -f "$residuals" "$bp_out"

  fi

  # ── Step 2: SUSAN smoothing ─────────────────────────────────────────────────
  if [[ "$DO_SUSAN" -eq 1 ]]; then
    p50=$(fslstats "$tmean_out" -k "$mask" -p 50)
    bt=$(echo "$p50 * 0.6666667" | bc -l)
    echo "SUSAN: p50=$p50  bt=$bt  dt=$DT"
    run_cmd "susan \"$denoise_target\" $bt $DT $SUSAN_DIM $USE_MEDIAN $N_USANS \"$tmean_out\" $bt \"$denoised_out\""
    rm -f "$denoise_target"
  fi

done

echo "════════════════════════════════════════"
echo "Done."
