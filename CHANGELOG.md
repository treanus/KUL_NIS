# Changelog

## Unreleased (working tree, 2026-07-09)

### KUL_clinical_fmridti.sh
- Replaced the `msbp` (Docker MultiScaleBrainParcellator) path with a
  `multiparc` flag; the `KUL_run_msbp` function (Docker-based) was removed
  and its call site disabled — replaced by the native-FS `KUL_run_multiparc`
  path for types 4/5/6.
- `KUL_run_multiparc` now only runs when `multiparc=1`, and its `.done`
  log tag was renumbered `10_multiparc` -> `09_multiparc`.
- FreeSurfer version used by VBG bumped 6.0.0 -> 8.2.0.
- `KUL_run_VBG` now passes `-M -O -H` to `KUL_VBG.sh` in addition to the
  existing flags.
- Pipeline step comments/log tags renumbered after the MSBP removal.
- Added `-E spm|nilearn` (default `spm`) to choose the fMRI GLM engine;
  `KUL_fmriproc` now dispatches to `KUL_fmriproc_spm_new.sh` or
  `KUL_fmriproc_nilearn_new.sh` accordingly, passing `-c $ncpu` to both so
  their auto core-scheduling (see below) is actually used.
- Added scilpy conda-env validation for `-f`/FWT — fixed on this pass:
  the original `[ ! $scilpy ]` / `[ -z $scilpy ]` checks were unquoted and
  inverted (the "env set to X" message printed exactly when the env was
  *not* set); now uses `[ -z "$scilpy" ]` / `[ -n "$scilpy" ]` correctly.

### KUL_dwiprep.sh
- Added `-x <dir>` (custom preproc/output directory) and `-f <conda_env>`
  (conda env with lore_sd installed) options.
- lore_sd now requires `-f`; the script activates/deactivates that conda
  env around the `lore_dwi2decomposition` step via new
  `activate_conda_env`/`deactivate_conda_env` helpers.
  Fixed on this pass: the error message shown when `-f` is missing
  incorrectly told the user to supply `-c <conda_env_name>` — `-c` is an
  unrelated pre-existing flag (`shard_recon`) in this script; corrected to
  reference `-f`.
- `lore_dwi2decomposition` call no longer passes `--mask ${dwi_mask_input}`
  (commented out) — confirm whether this was intentional.
- Reindented `KUL_dwiprep_convert` (whitespace only).

### KUL_dwiprep_anat.sh
- `mrtransform ... dwi_preproced_reg2T1w.mif` now passes
  `-reorient_fod no`. Note: `-force` currently appears twice on that line
  (`-nthreads $ncpu -force -reorient_fod no  -force`) — harmless but looks
  like a stray duplicate from editing.

### KUL_fmriproc_spm_new.sh
- Reworked into a parallel job-pool scheduler: new `-c` (auto-schedule
  from a total core budget) or manual `-j`/`-J`/`-T` knobs; new
  `KUL_throttle` (bash job-pool limiter) and `KUL_prep_run` (per-run
  mask+SUSAN+confounds prep) functions.
- Preprocessing (mask+SUSAN+confounds) for all runs of all tasks now runs
  in one parallel pool before any GLM starts (previously serial, per-task).
- GLM analyses (per-run + per-task aggregates) for all tasks now dispatch
  into one parallel pool instead of a nested per-task/per-run serial loop.
- Version string bumped to "v2.0 - dd 03/07/2027" (note: date looks like
  a typo for 2026).

### KUL_fmriproc_nilearn_new.sh (new, was untracked)
- New Python/Nilearn port of the SPM GLM pipeline — same job/task
  structure and the same `-c`/`-j`/`-J`/`-n`/`-T` scheduling design as the
  reworked `KUL_fmriproc_spm_new.sh`, but running the GLM in
  nilearn instead of MATLAB/SPM12 (no MATLAB/SPM license required).
- Companion `share/nilearn/KUL_nilearn_glm.py` holds the actual GLM code
  called by the wrapper script.
- Note: `share/nilearn/` also contains an older, smaller duplicate copy of
  `KUL_fmriproc_nilearn_new.sh` (different permissions, `600`/no-exec vs.
  `700` for the top-level copy) — looks like a stray leftover from copying
  between workstations rather than intentional; worth deleting once
  confirmed.

### KUL_preproc_all.sh / study_config templates
- Wired the new `KUL_dwiprep.sh -f <loresd_env>` flag through the config
  chain: added a `loresd_env:` key to the three `run_dwiprep_lore_sd.txt`
  study_config templates (`clinical_fmri_dmri`, `clinical_dmri_dbs_drt`,
  `clinical_dmri_dbs_hdp`), parsed in `KUL_preproc_all.sh` and passed as
  `-f $loresd_env` whenever `dwiprep_options` includes `lore_sd`.
  Without this, any pipeline run through those lore_sd configs would now
  hard-fail against KUL_dwiprep.sh's new mandatory `-f` check.
  **Action needed:** fill in `loresd_env:` in each of the three templates
  with the actual conda env name before running lore_sd.

### Known issue not fixed
- `KUL_NIS_old` (sibling directory) has its own unrelated uncommitted
  changes (atlas file deletions, an in-progress `atlasses`->`atlases`
  rename, stray `VSC/`/`atlases/` untracked dirs) — left untouched.
