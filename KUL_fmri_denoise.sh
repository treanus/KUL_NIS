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
GSR=0
ACOMPCOR_N=5

# How the task design is treated during confound regression. Matters because the
# confound set is routinely collinear with a block paradigm: on a 30s-on/30s-off
# language run, the 22 default regressors explained 52.7% of the design's
# variance (motion alone 27.4%, aCompCor 23.2%) -- so plain regression stripped
# ~58% of the task signal and left melodic with nothing task-locked to find.
#   ignore   regress confounds as-is (historical behaviour, still the default)
#   preserve orthogonalise confounds against the design first, so they can only
#            remove non-task variance -- for task analyses (melodic + GLM)
#   remove   additionally regress the design out, turning a task run into
#            genuine pseudo-rest -- for resting-state-style connectivity
# Rest runs have no design, so every mode is equivalent for them.
TASK_SIGNAL="ignore"
EVENTS_DIR=""

# Skip confound regression altogether (high-pass + smoothing only). This is what
# FSL's own FEAT/MELODIC workflow feeds ICA, and for a decomposition it is the
# right input: separating neural signal from motion and physiological artifact is
# the job ICA is being asked to do, so regressing confounds first destroys the
# structure it would have isolated. Measured on a 30s-on/30s-off language run,
# melodic found no task component at all after the standard 29-regressor denoise
# (best IC |r|=0.27); with regression skipped it recovered one at |r|=0.65 whose
# spatial map correlates 0.57 with the subject's own GLM z-map.
NO_CONFOUNDS=0

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
  --acompcor-n <N>          Number of WM/CSF aCompCor components each (default: 5).
                            More components capture more physiological noise
                            variance without touching gray-matter signal directly.
  --gsr                     Add global_signal + its derivative to the confound
                            regression (fMRIPrep already computes these in the
                            confounds TSV -- no extra computation needed). Composes
                            with any of the options above. Standard fix for
                            widespread whole-brain "hyperconnectivity" from
                            unremoved physiological/global noise, but can introduce
                            artificial negative correlations elsewhere -- a known,
                            debated tradeoff; opt-in rather than default for that
                            reason.
  --cols "<c1 c2 ...>"      Fully override confound columns
  --cols-file <file.txt>    Read confound columns from file (one per line, # comments ok)
  --no-confounds            Skip confound regression entirely (high-pass +
                            smoothing only). The right input for ICA/melodic --
                            see the note at NO_CONFOUNDS in this script.
                            --lp 0 additionally disables the low-pass.
  --task-signal <mode>      How to treat the task design during regression:
                              ignore   (default) regress confounds as-is
                              preserve orthogonalise confounds against the design,
                                       so they cannot remove task variance
                              remove   also regress the design out (pseudo-rest)
                            Needs an events TSV; a run whose task has none falls
                            back to 'ignore' with a warning. No-op for rest runs.
                            NOTE: 'preserve' and 'remove' produce genuinely
                            different data -- give each its own --out directory
                            or one will be mistaken for the other on reuse.
  --events-dir <dir>        Where to look for task-<TASK>_events.tsv
                            (default: ./study_config, then the BIDS root)

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
    --acompcor-n) ACOMPCOR_N="$2"; shift 2;;
    --gsr)        GSR=1; shift;;
    --cols)       COLS_STR="$2"; shift 2;;
    --cols-file)  COLS_FILE="$2"; shift 2;;
    --no-confounds) NO_CONFOUNDS=1; shift;;
    --task-signal) TASK_SIGNAL="$2"; shift 2;;
    --events-dir) EVENTS_DIR="$2"; shift 2;;
    -h|--help)    usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -n "$FMRIPREP_DIR" && -n "$OUT_DIR" ]] || { usage; exit 1; }
[[ "$METHOD" == "nilearn" || "$METHOD" == "fsl" ]] || {
  echo "ERROR: --method must be 'nilearn' or 'fsl'"; exit 1
}
case "$TASK_SIGNAL" in
  ignore|preserve|remove) ;;
  *) echo "ERROR: --task-signal must be 'ignore', 'preserve' or 'remove'"; exit 1;;
esac
# preserve/remove need the design, which only the nilearn path can apply -- the
# fsl path builds its regression matrix with fsl_glm and would silently ignore it
if [[ "$TASK_SIGNAL" != "ignore" && "$METHOD" != "nilearn" ]]; then
  echo "ERROR: --task-signal $TASK_SIGNAL requires --method nilearn"; exit 1
fi

# Resolve task-<TASK>_events.tsv for a run. Events live in study_config/ in the
# KUL layout, but a plain BIDS tree keeps them at the dataset root or beside the
# BOLD, so try all three before giving up.
events_for_task() {
  local task="$1" cand
  for cand in \
    "${EVENTS_DIR}/task-${task}_events.tsv" \
    "$(pwd)/study_config/task-${task}_events.tsv" \
    "${FMRIPREP_DIR}/../BIDS/task-${task}_events.tsv" \
    "${FMRIPREP_DIR}/../task-${task}_events.tsv"; do
    [[ -n "$EVENTS_DIR" || "$cand" != "/task-${task}_events.tsv" ]] || continue
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  return 1
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
# Default: 12 motion parameters (6 + derivatives) + N WM aCompCor + N CSF aCompCor
# (N set by --acompcor-n, default 5). --all-confounds expands to 24 motion
# params (adds squares + squared derivatives). --gsr adds global_signal on top
# of whichever of these (or --cols/--cols-file) is selected. Outlier spike
# regressors (motion_outlier*, non_steady_state_outlier*) are added
# automatically per-run from the confounds TSV header.
ACOMPCOR_W_COLS=()
ACOMPCOR_C_COLS=()
for ((_i=0; _i<ACOMPCOR_N; _i++)); do
  _idx=$(printf "%02d" "$_i")
  ACOMPCOR_W_COLS+=("w_comp_cor_${_idx}")
  ACOMPCOR_C_COLS+=("c_comp_cor_${_idx}")
done

DEFAULT_COLS=(
  trans_x trans_y trans_z rot_x rot_y rot_z
  trans_x_derivative1 trans_y_derivative1 trans_z_derivative1
  rot_x_derivative1 rot_y_derivative1 rot_z_derivative1
  "${ACOMPCOR_W_COLS[@]}" "${ACOMPCOR_C_COLS[@]}"
)

ALL_COLS=(
  trans_x trans_y trans_z rot_x rot_y rot_z
  trans_x_derivative1 trans_y_derivative1 trans_z_derivative1
  rot_x_derivative1 rot_y_derivative1 rot_z_derivative1
  trans_x_power2 trans_y_power2 trans_z_power2
  rot_x_power2 rot_y_power2 rot_z_power2
  trans_x_derivative1_power2 trans_y_derivative1_power2 trans_z_derivative1_power2
  rot_x_derivative1_power2 rot_y_derivative1_power2 rot_z_derivative1_power2
  "${ACOMPCOR_W_COLS[@]}" "${ACOMPCOR_C_COLS[@]}"
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

if [[ "$GSR" -eq 1 ]]; then
  COLS+=(global_signal global_signal_derivative1)
fi

if [[ "$NO_CONFOUNDS" -eq 1 ]]; then
  COLS=()
else
  [[ "${#COLS[@]}" -gt 0 ]] || { echo "ERROR: no confound columns selected"; exit 1; }
fi

DT=$(echo "$FWHM/2.355" | bc -l)   # SUSAN spatial sigma (mm)

# ── Find BOLD files ────────────────────────────────────────────────────────────
mkdir -p "$OUT_DIR"

mapfile -t BOLDS < <(
  find "$FMRIPREP_DIR" \
    -path "$FMRIPREP_DIR/derivatives" -prune -o \
    -type f -path "*/func/*" \
    -name "*_space-${SPACE}*_desc-preproc_bold.nii.gz" -print \
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
  # --no-confounds means no regressors at all, spikes included: they still cost
  # degrees of freedom, and the ICA input this mode exists for is meant to be
  # unregressed
  OUTLIERS_STR=""
  [[ "$NO_CONFOUNDS" -eq 1 ]] || OUTLIERS_STR="$(get_outlier_cols "$conf")"
  RUN_COLS=("${COLS[@]}")
  OUTLIER_COUNT=0
  if [[ -n "$OUTLIERS_STR" ]]; then
    read -r -a OUTLIERS <<< "$OUTLIERS_STR"
    RUN_COLS+=("${OUTLIERS[@]}")
    OUTLIER_COUNT="${#OUTLIERS[@]}"
  fi
  echo "Regressors: ${#RUN_COLS[@]} (${OUTLIER_COUNT} spike columns)"

  # Resolve this run's design. A run with no task entity, or a task with no
  # events file, degrades to 'ignore' rather than failing -- rest runs legitimately
  # have no design, and a missing events TSV is a study-config gap, not a reason
  # to abandon an otherwise valid denoise.
  RUN_TASK_SIGNAL="$TASK_SIGNAL"
  RUN_EVENTS="NONE"
  if [[ "$TASK_SIGNAL" != "ignore" ]]; then
    run_task="$(sed -n 's/.*_task-\([^_]*\)_.*/\1/p' <<< "$base")"
    if [[ -z "$run_task" ]]; then
      echo "  task-signal: no task entity in filename — using 'ignore'"
      RUN_TASK_SIGNAL="ignore"
    elif [[ "$run_task" == "rest" ]]; then
      # expected and uninteresting: rest has no design to preserve or remove,
      # so say so quietly rather than warning about a missing events file
      echo "  task-signal: rest run, nothing to preserve or remove — using 'ignore'"
      RUN_TASK_SIGNAL="ignore"; RUN_EVENTS="NONE"
    elif RUN_EVENTS="$(events_for_task "$run_task")"; then
      echo "  task-signal: ${TASK_SIGNAL}  (design from $(basename "$RUN_EVENTS"))"
    else
      echo "  WARN: task-signal ${TASK_SIGNAL} requested but no events TSV for task-${run_task} — using 'ignore'"
      RUN_TASK_SIGNAL="ignore"; RUN_EVENTS="NONE"
    fi
  fi

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
      "$RUN_TASK_SIGNAL" "$RUN_EVENTS" \
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
# 0 disables a cutoff -- nilearn takes None, not 0, to mean "no filter on this
# side", and a literal 0 Hz high-pass would be a no-op that still costs a
# butterworth pass
hp       = float(sys.argv[7]) or None
lp       = float(sys.argv[8]) or None
task_sig = sys.argv[9]
events_f = sys.argv[10]
cols     = sys.argv[11:]

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

# --no-confounds passes no columns at all; nilearn wants None rather than an
# empty (n_timepoints, 0) array
confounds_array = np.array(matrix, dtype=np.float64) if cols else None

# ── Task design handling ──────────────────────────────────────────────────────
# Built here rather than taken from a canned .mat so it follows this run's actual
# events, TR and volume count instead of assuming a fixed paradigm length.
def build_design(events_file, n_vols, tr):
    import pandas as pd
    from nilearn.glm.first_level import make_first_level_design_matrix
    ev = pd.read_csv(events_file, sep="\t")
    if not {"onset", "duration"} <= set(ev.columns):
        return None, "events TSV lacks onset/duration"
    if "trial_type" not in ev.columns:
        ev = ev.assign(trial_type="task")
    frame_times = np.arange(n_vols) * tr
    dm = make_first_level_design_matrix(
        frame_times, ev, hrf_model="glover", drift_model=None)
    # drop the intercept: it carries no task information and would make the
    # orthogonalisation below strip each confound's mean, which clean_img's own
    # detrending already handles
    task_cols = [c for c in dm.columns if c != "constant"]
    if not task_cols:
        return None, "design has no task columns"
    return dm[task_cols].to_numpy(dtype=np.float64), None

design = None
if task_sig != "ignore" and events_f != "NONE":
    n_vols_hdr = nib.load(bold_f).shape[-1]
    design, why = build_design(events_f, n_vols_hdr, tr)
    if design is None:
        print(f"  WARN: {why} — falling back to task-signal=ignore")
        task_sig = "ignore"

mask_img = nib.load(mask_f)

# Compute and save temporal mean before cleaning (clean_img removes it via detrending)
bold_img = nib.load(bold_f)
bold_data = bold_img.get_fdata(dtype=np.float32)
tmean_data = bold_data.mean(axis=-1)
nib.save(nib.Nifti1Image(tmean_data, bold_img.affine, bold_img.header), tmean_f)
print(f"  Tmean saved: {tmean_f}")

print(f"  {len(rows)} timepoints  |  {len(cols)} regressors  |  TR={tr}s  |  bandpass {hp}–{lp} Hz")

if design is None:
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
    cleaned_data = cleaned.get_fdata(dtype=np.float32)
    affine, header = cleaned.affine, cleaned.header
else:
    # Joint model: fit [design | confounds] together, then subtract only the
    # confound part (preserve) or both parts (remove).
    #
    # Simply orthogonalising the confounds against the design beforehand is NOT
    # equivalent and does not work here: clean_img band-passes the confounds
    # afterwards, which destroys the orthogonality, and with ~29 regressors
    # against the ~40 effective DoF the 0.008-0.09 Hz passband leaves, the fit
    # then reabsorbs the task. Measured on a 30s-on/30s-off language run:
    # orthogonalise-then-filter recovered nothing (0.203 vs 0.222 for plain
    # regression), while this joint model gives 0.729 -- above the 0.523 of the
    # unfiltered input, because the task is protected while the bandpass still
    # removes noise. Estimating the confound betas *controlling for* the design
    # is what makes it exact rather than approximate.
    from nilearn.signal import clean as _clean
    ck = dict(detrend=True, standardize=None, t_r=tr, high_pass=hp, low_pass=lp)

    mask_flat = mask_img.get_fdata(dtype=np.float32).reshape(-1) > 0
    n_t = bold_data.shape[-1]
    Y = bold_data.reshape(-1, n_t)[mask_flat].T.astype(np.float64)

    # every term goes through the identical filter, so the regression happens in
    # one consistent space
    Yf = _clean(Y, **ck)
    Df = _clean(design, **ck)
    # --no-confounds leaves nothing to regress; 'remove' still has work to do
    # (strip the design), 'preserve' becomes a no-op beyond filtering
    Cf = _clean(confounds_array, **ck) if confounds_array is not None else None

    n_d = Df.shape[1]
    n_c = 0 if Cf is None else Cf.shape[1]
    X = np.column_stack([Df] + ([Cf] if Cf is not None else []) + [np.ones(n_t)])
    beta, *_ = np.linalg.lstsq(X, Yf, rcond=None)
    resid = Yf - (Cf @ beta[n_d:n_d + n_c] if Cf is not None else 0.0)
    if task_sig == "remove":
        resid = resid - Df @ beta[:n_d]
    print(f"  task-signal={task_sig}: joint model with {n_d} design + "
          f"{n_c} confound regressors")

    cleaned_data = np.zeros((mask_flat.size, n_t), dtype=np.float32)
    cleaned_data[mask_flat] = resid.T.astype(np.float32)
    cleaned_data = cleaned_data.reshape(bold_data.shape)
    affine, header = bold_img.affine, bold_img.header

# Restore temporal mean and re-apply mask
result_data = (cleaned_data + tmean_data[..., np.newaxis]).astype(np.float32)
mask_arr = nib.load(mask_f).get_fdata(dtype=np.float32) > 0
result_data *= mask_arr[..., np.newaxis]
nib.save(nib.Nifti1Image(result_data, affine, header), out_f)
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
