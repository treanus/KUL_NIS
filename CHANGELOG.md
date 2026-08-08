# Changelog

## Unreleased (committed locally, 2026-08-08 — Karawun lesion label, FAT1w QA, lesion/perfusion PACS)

### Fixed

- **`KUL_FAT1w.py` had become orphaned.** Nothing called it, so
  `BIDS/derivatives/KUL_compute/sub-X/KUL_FAT1/FAT1w.nii.gz` was never created,
  and the `if [ -f $FAT1 ]` guard in `KUL_karawun_prepare.sh` silently left
  `FAT1w=""` on every run — the QA volume was quietly missing from every
  Brainlab export. `KUL_karawun_prepare.sh` now generates it itself (no `-s`:
  smoothing blurs the very boundaries the image exists to let you check). It
  needs `dwiprep/sub-X/sub-X/qa/fa_reg2T1w.nii.gz`, so `KUL_dwiprep_anat.sh`
  must have run; otherwise the step is skipped with a message.

  This matters because FAT1w (= √FA · T1w, Goedemans et al., Imaging Neurosci
  2024) is what makes the FA-to-T1w registration inspectable in Brainlab. If
  that registration has slipped, every tract in the scene is displaced the same
  way and nothing else in the export would show it.

### New

- **Lesion as a Brainlab label.** `sub-X_lesion_and_cavity.nii.gz` (types 1/2),
  falling back to the manual `lesion.nii.gz` (type 3), is regridded onto the
  Karawun T1w grid and written to `Karawun/sub-X/labels/Lesion.nii.gz` with
  palette colour 16 — free in the 0-30 range and clear of the fMRI task labels,
  which count up from 1. Resampled with `GenericLabel` (see below), which
  interpolates each label's indicator and takes the argmax, so boundaries follow
  the anatomy rather than the grid and no label is invented that was not in the
  input.
- **Lesion in `RESULTS/sub-X/Anat/`.** Binarised to
  `sub-X_lesion.nii.gz` next to the other T1w-space volumes that feed the
  figures and the PACS export, rather than living only in `Lesion/`.
- **Lesion and DSC perfusion PACS series.** Under `-R`, the lesion mask and the
  NAWM-normalised `nrCBV_corrected` / `nrCBF` maps are rendered through the same
  `_render_one_spm` path as the fMRI maps, into `PACS/Clinical_<underlay>` and
  `Clinical_figures_<underlay>` so they stay separate from the SPM/Melodic
  series. The perfusion maps use fixed thresholds (1.75 on nrCBV, the
  conventional high-grade glioma cutoff; 1.0 on nrCBF, "above contralesional
  normal WM") rather than the auto `max/3` used for activation maps, since these
  are normalised ratios whose thresholds carry a fixed clinical meaning. `-T`
  still overrides. The user's interactive fMRI map selection (`_dcm_spm_set`) is
  saved and cleared around this pass — a non-empty selection means "only these"
  and would otherwise silently exclude the new series.

### Fixed (found while reviewing the above)

- **The lesion PACS export called `KUL_resolve_lesion` before it was defined.**
  `KUL_clinical_fmridti.sh` runs the results/PACS work inside a top-level
  `if [ $results -gt 0 ]` block that ends in `exit`, and the helper was defined
  ~600 lines further down. Bash registers a function when its definition
  *executes*, not when the file is parsed, so on the `-R` path that definition
  was never reached. It failed silently: "command not found" on stderr, the `if`
  took the else branch, exit status 0, and the lesion was simply absent from the
  export. `bash -n` does not catch this. Both helpers now live above that block,
  with a comment explaining why they cannot be moved back.

### Fixed (found by running the Brainlab import on real data)

- **FAT1w overflowed the DICOM `LargestImagePixelValue` tag.** It was copied to
  the Karawun folder unscaled. karawun writes `(0028,0106)`/`(0028,0107)` from
  the image's *original* intensities rather than from the values it rescales the
  pixel data to, and both tags are US — capped at 65535. FAT1w is
  sqrt(FA) * T1w, which reached 97575 on the first real case, so
  `importTractography` aborted with
  `'H' format requires 0 <= number <= 65535 ... (0028,0107) US: 97575`.
  It is now rescaled into 16-bit range before being handed over, the same way
  the T1w beside it already was.

  This path had never executed before — `KUL_FAT1w.py` was orphaned, so the
  volume was never produced and the copy never ran. Restoring the QA volume is
  what exposed it.

### Changed (interpolation)

- `KUL_dsc_perfusion.sh` passes `-p BSpline` to `antsMotionCorr` instead of
  taking its `Linear` default. Motion correction is the first of two resamplings
  the series goes through, and linear interpolation softens the data before the
  sinc-interpolated distortion correction ever sees it — giving away the
  resolution that second step exists to preserve. The blur leaks arterial signal
  (an order of magnitude more dR2* than tissue) into neighbouring voxels and
  biases peritumoral rCBV upward; unlike a per-voxel scale factor, that leak
  does not cancel in the S(t)/S0 ratio.
- The Karawun lesion label is resampled with `antsApplyTransforms -n GenericLabel`
  rather than `mrgrid -interp nearest`. Today this is a no-op — Karawun's
  `T1w.nii.gz` is the RESULTS T1w with only its intensities rescaled, so the
  grids are identical and any interpolator is exact (verified: 2969 -> 2969
  voxels) — but that is an assumption about an upstream step, and nearest is the
  option that degrades worst if it stops holding.

### Karawun palette

`lookup_cie()` clamps any label value above the palette size (31 entries,
indices 0-30) to the last colour, printing "Error - too many labels". The tract
table assigns colours 31-41 to 13 bundles and `next_auto_color` started at 100,
so all of those rendered as the *same* colour in Brainlab.

`next_auto_color` is now 42, one past the table's highest index. Against the
stock palette that still clamps, so it is no worse today; against an extended
palette it does the right thing.

The palette itself is a karawun-side change, committed in the karawun clone at
`/opt/kul_software/src/karawun` on branch `kul-extended-palette` (see
`KUL_PALETTE_NOTES.md` there): 31 -> 64 entries, indices 0-30 byte-identical to
upstream and verified equal at the DICOM value level, so existing scenes are
unaffected. The 33 appended colours are each at least CIEDE2000 11.6 from every
other entry. For context the *existing* palette's worst pair is dE 2.80 (labels
3 and 26), below the just-noticeable threshold. Extending it makes the existing
31-41 table assignments work as intended with no further KUL_NIS change.

That fork must be pinned for the `KarawunDev` install to benefit; `KarawunEnv`
is a separate, non-editable conda-forge install that will not pick it up.

## Unreleased (committed locally, 2026-08-08 — DSC perfusion support)

Brings the standalone DSC workflow in `tools/KUL_DSC_analysis/`
(`DSC_proc_script_WIP3.sh` + `Good_DSCLC_fit5.py`) into KUL_NIS as a first-class
step for tumour patients. The originals are left in place untouched.

### New

- **`KUL_dsc_perfusion.sh`** — preprocessing (denoise, motion correction,
  PE-restricted SyN distortion correction, N4), leakage-corrected
  rCBV/rCBF/MTT/TTP/TT0/K1/K2, and normalisation against contralesional NAWM.
  All maps land in the participant's T1w space, i.e. the same space as the
  fMRI/tractography results, so they can go to PACS and Karawun. Stats go to
  `RESULTS/sub-X/Perfusion/*_perfusion_{stats,summary}.tsv`.
- **`share/dsc/KUL_dsc_fit.py`** — the quantification engine.
- **`docs/KUL_dsc_perfusion/KUL_dsc_perfusion.md`**.

### Changed

- `KUL_dcm2bids.sh`: new `DSC` identifier writing the 4D series to the `perf/`
  datatype directory (sibling of `anat`/`func`/`dwi`) as `*_dsc.nii.gz`, with
  the `pe_dir` column carried into the sidecar as `PhaseEncodingDirection`.
  `**/perf/*dsc*` added to `.bidsignore`, as for ASL — BIDS has no standardised
  DSC suffix yet. Matched on `SeriesDescription` alone: DSC `ImageType` strings
  vary too much between vendors for a tighter criterion to be safe.
- `study_config/sequences.txt`: a `DSC,T2_DSC_Perfusion,-,2,j` entry.
- `KUL_clinical_fmridti.sh`: runs DSC automatically whenever a series is present
  (the way fMRI/dMRI are picked up), placed after VBG/multiparc so the
  FreeSurfer aseg the NAWM reference needs exists. New `-W` skips it.

### Fixed (carried over from the prototype)

- `scipy.integrate.cumtrapz`/`simps` were removed in scipy ≥ 1.14 — the fit
  could not import at all on a current environment.
- **`DenoiseImage -d 4` corrupted the first and last frames of the series.** It
  treats time as a fourth spatial axis, so its neighbourhoods are truncated at
  the boundaries; measured at ~13 % darkening on frame 0 of a 50-frame series,
  as deep as the bolus itself, contaminating rCBV in every voxel. The time axis
  is now mirror-padded before denoising and cropped afterwards.
- **The baseline window was hardcoded to frames 5–10.** A bolus arriving before
  frame 10 put S₀ part-way down the bolus, flipping the sign of ΔR2\* and of
  every downstream map, MTT included, with no warning. Now detected from the
  data, with a hard guard on a non-positive AIF area.
- TE was scraped from the `mrinfo` comments field with `"0.0"` prepended; now
  read from the sidecar `EchoTime`, with a millisecond-vs-second sanity check.
- The AIF candidate pool used an absolute signal-drop threshold that crashed PCA
  when nothing cleared it, and was a no-op when anything did.
- `mri_synthstrip -g` was unconditional, aborting on CPU-only nodes and on GPU
  nodes whose FreeSurfer ships a CPU-only torch. Now probed once.
- Hardcoded `/mnt/DATA1/aradwa0/...` paths and a directory-dependent
  `python ./Good_DSCLC_fit5.py` call.

### Method changes (opt out with `-L`)

- Deconvolution now uses a `dt`-scaled AIF with truncated SVD (0.2·S_max)
  instead of a z-scored AIF with `1/(S+1e-3)` damping, so rCBF is in 1/s and MTT
  in seconds. On a phantom varying only bolus amplitude the legacy path gave an
  rCBF lesion/normal ratio of 1.04 (no contrast); the new one gives 2.94, with
  MTT correctly amplitude-invariant at 1.01.
- Distortion correction is a direct rigid+affine+SyN to the anatomy with the
  deformation restricted to the phase-encoding axis, replacing
  `antsIntermodalityIntrasubject.sh` — which demanded a template and a
  subject-to-template warp solely to emit template-space outputs nothing
  consumed, making a full SyN registration to MNI a prerequisite of every run.
- The per-voxel leakage-fit and deconvolution loops are vectorised.
- The leakage model still regresses on the integral of the AIF rather than on a
  non-enhancing reference tissue curve as textbook Boxerman-Schmainda-Weisskoff
  does. Left as the prototype had it rather than changed silently; flagged in
  the docs under "Known deviations" as worth revisiting.

Committed locally, not pushed. Validated end-to-end on a synthetic DSC phantom
with known ground truth (see `/opt/kul_software/tmp/KUL_DSC_phantom/`), not yet
on clinical data -- the phantom cannot exercise AIF detection realistically, so
the first real case is the actual test. Check the auto-detected `baseline_frames`
in the log on that run.

## Unreleased (working tree, 2026-07-28 — auto-discover bundles in KUL_karawun_prepare.sh)

`KUL_karawun_prepare.sh` unconditionally attempted a hardcoded list of ~40
named tracts every run, regardless of `-t`/processing type -- `-t` never
actually filtered this list, it only changes the threshold formula for type 1
(see `KUL_karawun_get_tract`). In practice this meant any bundle actually
generated by this patient's FWT config that wasn't already in that hardcoded
list (a new bundle added to FWT, or a custom one) was silently never picked
up, no matter what `-t` was set to.

Replaced it with the same auto-discovery pattern `KUL_clinical_fmridti.sh`
already uses for screenshots/PACS export: glob over
`.../FWT/sub-<participant>_TCKs_output/*_output`, deriving the bundle name
from each directory name rather than assuming a fixed list. Known bundles
still get their existing Brainlab display name/color/thresholds from a
lookup table (verified byte-for-byte identical to the previous ~40 hardcoded
values); anything not in the table falls back to sensible defaults (raw FWT
name, auto-assigned color starting at 100) instead of being skipped. The two
CSHDP VOI-derived labels (distal STN-motor) stay explicit, since they come
from a different source path (`_VOIs`, not `_output`) that the auto-discovery
glob doesn't cover.

Not yet committed -- working-tree edit only, pending your review/testing
(per your standing instruction not to commit in KUL_NIS without asking each
time).

## Unreleased (working tree, 2026-07-11 — extend status-reporting pass to flagged scripts)

Follow-up to the 2026-07-10 pass below: a final grep sweep at the end of that pass
turned up more `.done`-touching scripts outside its original scope
(`KUL_preproc_all.sh`/`KUL_clinical_fmridti.sh` only). Applying the same audit here.

### KUL_DRT.sh
- `KUL_run_msbp`/`KUL_run_FWT`: gated their `KUL_task_exec` calls on return value
  before touching `MSBP.done`/`FWT.done`, matching the established pattern.
- The MAIN flow's separate, currently-live MSBP block (a raw `docker run`, not
  routed through `task_exec`/`KUL_task_exec` at all — the `KUL_run_msbp` function
  above it is actually dead code, never called) had zero status checking before
  its own unconditional `touch MSBP.done`. Added an explicit exit-code check.
- STEP 6's per-session `KUL_run_FWT` loop and STEP 5's `KUL_dwiprep_anat` call
  discarded their return values entirely; now warn/abort respectively.

### KUL_dwiprep_anat.sh
- Three blocks (`status.mrtransformNL.done`, `status.freesurfer.done`,
  `status.labelconvert.done`) each run several raw FreeSurfer/MRtrix commands
  (none via `task_exec`/`KUL_task_exec`) followed by an unconditional `touch`.
  Wrapped each block's commands in a `( set -e; ... )` subshell so any command
  failing anywhere in the block aborts it, and gated the `touch` on the
  subshell's exit code instead of always running regardless.

## Unreleased (working tree, 2026-07-10 — logging/status-reporting correctness pass)

Prompted by inconsistent status reporting across real test runs — some steps logged
"success" despite the underlying command having actually failed, because their
`KUL_task_exec` return value was never checked at the call site.

### KUL_main_functions.sh
- `kul_log_file`/`kul_errorlog_file` use a fixed, non-timestamped name per task label
  (e.g. `7_VBG.log`), always written with `>>`/`tee -a`, and were never truncated.
  Traced this to the exact real-world symptom that prompted this whole pass: a VBG run
  that failed fast (bad lesion path, "exitcode 2" in `.error.log`) followed by a retry
  219 minutes later that succeeded — but since neither file was ever cleared, the retry's
  "Success" sat right next to the first attempt's stale failure, reading as a
  contradiction. Now truncates both files right before each `KUL_task_exec` invocation
  writes to them.
- `KUL_task_exec` itself was already correct (captures the real exit code via `wait`).
  Fixed two smaller bugs inside it:
  - `kul_echo`'s `verbose_level == 0` case had no branch at all, so every message was
    silently discarded (not even written to the log file) when running with `-v 0`.
    Changed so level 0 behaves like level 1 (file-only, no terminal echo) — "quiet
    terminal" should never mean "no record exists."
  - The aggregate Success/Fail line at the end of `KUL_task_exec` used
    `${kul_log_file}` unindexed, which defaults to array index `[0]` — when multiple
    tasks are batched in one call, only the first task's log ever got this line. Now
    loops over every task index so each task's own log file gets the line.

### KUL_preproc_all.sh
- 4 batch `KUL_task_exec` call sites (mriqc, fmriprep, freesurfer, dwiprep) discarded
  the return value entirely — a failure for any participant in the batch was
  completely invisible; the script just moved on. Added a `|| kul_echo "WARNING: ..."`
  to each so a batch failure is now at least logged clearly. Not changed to `return`/
  `exit`, since these are fire-and-forget batch launches across a loop of participant
  groups with no per-call `.done` marker to gate — other batches should still run.

### KUL_clinical_fmridti.sh
- Gated the remaining ~7 ungated `KUL_task_exec` call sites, matching the pattern
  already used correctly at 3-4 other sites in this file
  (`KUL_task_exec ... || { kul_echo "...failed..."; return 1; }`):
  - `KUL_run_VBG`: previously, a failed `KUL_VBG.sh` run still had its (possibly
    incomplete/missing) output `cp -r`'d into the freesurfer derivatives directory
    as if nothing had gone wrong. Now returns before the copy on failure.
  - `KUL_run_dwiprep_anat`, `KUL_run_dwiprep_MNI`, `KUL_calc_DTI_ALPS`: each
    unconditionally touched its own `.done` marker regardless of whether the
    preceding step actually succeeded. Now gated the same way.
  - `KUL_fmriproc` (nilearn/spm branch, melodic/conn branch): these two don't touch a
    `.done` marker directly in this file (that happens inside the called sub-script
    itself, out of scope here) so a hard `return 1` isn't appropriate — instead added
    a clear `kul_echo` warning naming which `.done` marker will consequently not be
    created, so a failure here is no longer silent.

## Unreleased (working tree, 2026-07-10 — make FWT tractometry opt-in)

### KUL_clinical_fmridti.sh
- Added `-Q` (opt-in, off by default): runs `KUL_FWT_make_TCKs.sh -Q`
  (per-bundle tractometry). Previously `-Q` was hardcoded on in
  `KUL_run_FWT`. Found via a real timed test run: with ~50 bundle/hemisphere
  combinations processed strictly sequentially (no bundle-level
  parallelism), tractography+tractometry took 4+ hours for one participant.
  Tracts/tractography themselves are still generated either way — only the
  along-tract scalar-profiling pass is now gated behind `-Q`.

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
