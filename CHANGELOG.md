# Changelog

## Unreleased (working tree, 2026-07-09, batch 3 — wire -U through; fix Karawun labels)

### KUL_clinical_fmridti.sh
- Added `-U` (EXPERIMENTAL opt-in), passed straight through to
  `KUL_FWT_make_TCKs.sh -U` in `KUL_run_FWT`. Previously the rfa-modulated
  lore_sd FOD (see batch 2) was only reachable by calling
  `KUL_FWT_make_TCKs.sh` directly — `KUL_run_FWT`'s call to it used a fixed
  flag list with no pass-through. Without `-U`, behavior is unchanged.
- **Corrected the Karawun fMRI label output from batch 2.** The earlier
  implementation combined all tasks' thresholded/scaled maps into one
  `afMRI_multilabel.nii.gz` via voxel-wise max. That doesn't match how this
  is actually consumed: the existing `KUL_karawun2brainlab.sh` establishes
  the real convention — each tract/task is its own separately-scaled
  `.nii` file (`mrcalc <mask> <thresh> -gt <color> -mult`), and
  `importTractography -l <dir>/*` already accepts multiple separate label
  files directly. Combining was solving a non-existent problem. Now writes
  one `Karawun/sub-*/labels/afMRI_<taskname>.nii.gz` per task, each scaled
  to its own distinct integer value (still 1..N, tasks sorted
  alphabetically) — no merging, no `mrmath`, no tie-break-on-overlap
  behavior to worry about.

## Unreleased (working tree, 2026-07-09, batch 2 — wishlist items 1/3/5)

### KUL_fmriproc_spm_new.sh / KUL_fmriproc_nilearn_new.sh
- Hardwired the downstream-selected fMRI output: only the `wc` (with
  confounds) + `p001unc_k50` Bizzi-thresholded map per task now lands in
  `RESULTS/sub-*/SPM` (`globalresultsdir`). Every other combination — nc
  variant, FWE-thresholded variant, and the raw unthresholded map — now
  goes to the new `RESULTS/sub-*/SPM_all` (`globalresultsdir_all`) instead.

### KUL_clinical_fmridti.sh
- `SPM_all` added alongside `SPM` in the redo-cleanup and initial mkdir
  blocks.
- The report-figure loop (Bizzi-thresholded map -> REPORT) now only picks
  up the `*_wc` stats folder, matching the hardwired selection.
- New: a combined multi-label fMRI overlay for Karawun. In the interactive
  figure/threshold flow, each task's canonical SPM map is regridded onto
  `Karawun/sub-*/T1w.nii.gz` (continuous data regridded *before*
  thresholding, matching the existing VOI/tract-label idiom in
  `KUL_karawun_prepare.sh` — avoids corrupting a binary mask with
  interpolation), binarized at its resolved threshold, multiplied by a
  stable per-task integer label (tasks sorted alphabetically, 1..N), and
  voxel-wise max-combined into `Karawun/sub-*/labels/afMRI_multilabel.nii.gz`.
  Known tie-break: where two tasks' activations overlap spatially, the
  higher-numbered task's label wins (max-combine). Skipped with a message
  if Karawun prep hasn't produced `T1w.nii.gz` yet.

### KUL_preproc_all.sh / study_config templates
- `use_native_dwi` now defaults to `1` (native resolution) both in
  `KUL_preproc_all.sh`'s fallback and in every checked-in
  `study_config/**/run_dwiprep*.txt` template (previously defaulted to `0`,
  i.e. upsampled), matching `KUL_dwiprep.sh`'s own `-u` default.
- Verified (no code change needed): the lore_sd chain already threads
  through end-to-end regardless of native vs. upsampled — `KUL_dwiprep.sh`
  → `KUL_dwiprep_anat.sh` (FOD/contrast registration to T1w) →
  `KUL_FWT_make_TCKs.sh` (auto-prefers lore_sd ODF over dhollander CSD FOD
  when present).

### KUL_dwiprep.sh / KUL_dwiprep_anat.sh — experimental `rfa_modulated_fod`
- `KUL_dwiprep.sh`: after lore_sd's `rfa.mif` contrast is computed, also
  computes `response/lore_sd/rfa_modulated_fod.mif = odf.mif .* rfa.mif`
  (voxel-wise FOD-amplitude modulation by the rfa contrast) — a first
  experiment aimed at improving tractography specificity through
  pathology by suppressing FOD amplitude in low-rfa tissue. Computed
  unconditionally whenever lore_sd runs (cheap), on the native odf.mif/
  rfa.mif pair (same grid `KUL_dwiprep_anat.sh` already registers odf.mif
  from).
- `KUL_dwiprep_anat.sh`: registers `rfa_modulated_fod.mif` ->
  `rfa_modulated_fod_reg2T1w.mif` the same way as plain `odf.mif`.
- **Opt-in only** — see `KUL_FWT` changelog for the consumer side. Without
  explicitly opting in, tractography behavior is unchanged.

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
