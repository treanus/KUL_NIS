# Changelog

## Unreleased (2026-09-05 — resseg cavity segmentation rebuilt; Karawun T1w rescale was silently broken)

### resseg was finding cavities outside the brain, and its arbitration was a coin toss

`KUL_anat_segment_tumor.sh`'s two resseg runs both operated on inputs that let
the tool go wrong in different ways, and the step that reconciled them discarded
one result outright when they disagreed.

- **Run 1 now uses a Colin27 normal-brain mask warped to native space**
  (`antsRegistrationSyN.sh -t a` onto `colin27_t1_tal_lin`, inverse-applied)
  rather than the raw full head. An unmasked head lets resseg "find" a cavity on
  extracranial structures — the maxillary sinus is a dark fluid-filled shape of
  about the right size. A *template* mask is the right tool here precisely
  because, unlike the patient's own BET used by run 2, it cannot exclude a real
  post-surgical cavity.
- **Hallucination detection with one retry** (new `KUL_resseg_masked_run`).
  resseg was never trained on hard-zero-masked images and can place a cavity
  entirely inside the zeroed background. The fraction of predicted voxels landing
  outside the mask that built the input is measured; above 0.5 the run is retried
  once at a larger dilation (`npass` 5 → 15) and the sane result kept. Both runs
  now share the tighter margin, so mask *source* is the only difference between
  them.
- **Fixed seed 42.** resseg's test-time augmentation (`-a 3`) is stochastic, so
  results — and therefore the arbitration below — were not reproducible
  run-to-run.
- **STEP 5C arbitrates instead of discarding.** Previously, no overlap between
  the two runs meant resseg's result was silently dropped. Now: overlap keeps run
  1 as before; otherwise, when HD-GLIO found a lesion, the candidates are scored
  by `ImageMath MaurerDistance` to that lesion and the nearer one kept
  (single-nonempty cases handled explicitly); if still unresolved it prints a
  loud WARNING naming the decision path and recommending manual review. This mask
  feeds VBG for the remainder of the pipeline, so a silent wrong answer here
  propagates everywhere.

### Karawun's T1w rescale factor was computed as zero

`KUL_karawun_prepare.sh` had:

```bash
T1w_factor=$(scale=10; echo "($T1w_max-($T1w_min))/32767" | bc)
```

`scale=10;` runs as a *shell* assignment, not as bc input — so bc used its
default `scale=0` and did integer division. Whenever `max - min < 32767` the
16-bit scale factor truncated to **0**, and the subsequent `mrcalc ... -div`
produced a garbage T1w for Brainlab.

- `scale=10` moved inside the bc expression, plus a guard that copies the T1w
  unscaled with a warning if the factor comes back empty or zero, and logs the
  factor when it doesn't.

### `KUL_nii2dcm.py`

- **`InstanceNumber` (0020,0013) was 0-based** in both PNG-overlay writers, which
  mis-sorts in some viewers. Now `+1`, matching the fix `writeSlices()` already
  carried and which was never propagated to the other two.
- Removed the dead `donor_frames` IPP/SliceLocation table (built and sorted,
  never read).
- The PNG-vs-donor-frame-count mismatch message claimed "using nearest-neighbour
  matching", which never happens — geometry comes wholly from the underlay.
  Downgraded to an informational note.

### Config surface

- **`do_dwiprep_MNI`** added to all three dwiprep configs and gated in
  `KUL_run_dwiprep_MNI`. It previously ran for every type unconditionally, and
  nothing in KUL_NIS — DTI-ALPS included — reads its `MNI/` output.
- **`fwt_cutoff`** is passed to `KUL_FWT_make_TCKs.sh` as `-X` only when
  non-empty; empty means "let FWT decide" (MRtrix default for CSD, 0.05 for
  lore_sd). It was not reachable from this level before.
- **`do_fmri_glm`**, **`do_karawun`**, **`karawun_threshold_tumor`/`_dbs`** added
  as switches. Karawun now reports "switched off" distinctly from "this `-t` has
  no Karawun case".
- CST streamline targets raised to **8000** across all three types
  (`_base/tracks_list.txt` 5000 → 8000; DBS DRT/HDP 3000 → 8000), so every type
  now uses the same target.
- `_base/` `*_ncpu` placeholders lowered 64/32 → **16**, noting that
  `KUL_clinical_fmridti.sh` overwrites them from `-n` — the value only matters to
  a standalone `KUL_preproc_all.sh` run.

## Unreleased (2026-09-05 — nilearn engine over-populated SPM/, Karawun QQ copy broke its own links)

### The nilearn GLM engine put 3 files per task into `RESULTS/.../SPM/`

`SPM/` is contractually the single downstream-selected map per task (with
confounds, `p001unc_k50` thresholded); every other combination belongs in
`SPM_all/`. The SPM12 engine enforces that with an explicit per-file
conditional. `KUL_fmriproc_nilearn_new.sh`'s `KUL_populate_SPM_results` instead
globbed `afMRI_${task}${variant}*.nii`, which also matches `_wc.nii`,
`_wc_p001unc_k50.nii` **and** `_wc_FWE01_k50.nii` — so five real tasks put
fifteen files in `SPM/`.

That is not cosmetic downstream: `KUL_clinical_fmridti.sh` assigns one Karawun
palette slot per file found in `SPM/`, so five tasks consumed slots 51-65
instead of 51-55. Values above 63 exceed even the extended-palette fork, and
`importTractography` aborted with `Error - too many labels`.

- The glob is replaced by an explicit pick of the `_p001unc_k50` file for the
  chosen wc/plain variant, matching the SPM12 engine.

### `REPORT/sub-*_06_Tract_QQ/` copy broke the report's own iframe links

`KUL_FWT_bundle_report.py` embeds each bundle's spider page as
`<iframe src="<bundle>_output/QQ/<file>.html">` — a real relative path, since the
report and the per-bundle pages are not siblings at the original location. The
copy into `REPORT/` flattened every spider page into one directory (a comment
there still described the older bare-filename linking), so the copied report's
iframes pointed at subdirectories that did not exist in the copy. It worked in
`TCKs_output/` and broke once copied.

- The copy now mirrors the same `<bundle>_output/QQ/` nesting under the
  destination.

### Karawun instructions were printed exactly once, ever

`KUL_karawun_prepare.sh` prints the `importTractography` command, but prep only
runs on the first run that creates the folder; every later run took the "already
prepared" branch and said nothing further. So the command was visible only in
whichever run first built the folder — and a later `-F`/`-R` run can still have
added fMRI labels the folder did not have then.

- That branch now reprints the full command (including `FAT1w` when present).
- The end-of-run summary now distinguishes Karawun **prep** (automatic, every
  run) from the **fMRI activation labels** (added only by a subsequent `-F`/`-R`),
  which was the actual source of "when does Karawun run?" confusion.

### DSC perfusion sequence entries

- `study_config/_base/sequences.txt` carries both vendor spellings of the DSC
  series: Philips `T2_DSC_Perfusion` and Siemens `_perf_` (the latter taken from
  a production config; `mb`/`pe_dir` are protocol-specific and should be checked
  against the actual series).

## Unreleased (2026-09-03 — dead study_config files removed, KUL_DRT.sh sed fixed)

- **`KUL_DRT.sh`** had the same append-instead-of-replace bug in both of its
  config stamps (lines 128, 260): `s/BIDS_participants: /...${participant}/`
  appends to whatever the template already holds, so a template shipping
  `BIDS_participants: 001` produced `BIDS_participants: JaneDoe001`. Anchored to
  `^BIDS_participants:.*` so it replaces the value.
- **Nine dead files removed from the top of `study_config/`**: `sequences.txt`,
  `tracks_list.txt`, `run_{dwiprep,fmriprep,freesurfer}.txt`,
  `bids_filter_no_gadolinium.json` and three `task-*_events.tsv`. Only
  `KUL_scaffold` reads the repo-root `study_config/`, and only from `_base/` and
  the type directories, so nothing copied or read these. Every other reference
  in the tree is either usage text describing a path inside the *user's* study
  folder or a script reading the patient's own copy at runtime.

  `sequences.txt` carried uncommitted local edits (the DSC vendor comment and
  the `DSC,_perf_` row); both had already been merged into `_base/`, verified
  row by row before deleting. The only rows not carried over are
  `FLAIR,3D_FLAIR`, deliberately dropped because it collides with `FLAIR,FLAIR`,
  and `cT1w,T1_POST`, covered by the broader `cT1w,POST`.

  **Kept**, each with a live consumer: `sequences_expert.txt`
  (`KUL_dcm2bids_new.sh:901` reads the 7-column expert format),
  `subjects_and_options{,_expert_mode}.csv` (`KUL_preproc_all.sh`, batch and
  `-e` modes), `tracto_{rois,tracts}.csv` (`KUL_dwiprep_fibertract.sh`).
- **Two more per-type duplicates gone**: `clinical_dmri_dbs_{drt,hdp}/run_dwiprep.txt`
  differed from `_base/` in nothing but the hardcoded `dwiprep_ncpu`, which `-n`
  now supplies. `run_dwiprep_lore_sd.txt` went the same way once `--niter=8`
  was added to it, leaving the DBS type directories holding only
  `tracks_list.txt` and `run_fmriprep.txt`.

## Unreleased (2026-09-03 — config reads in KUL_preproc_all.sh are anchored)

Every config read was an unanchored `grep key $conf | grep -v \#`, which matches
any line *containing* the key. Nothing in the shipped configs collided, but only
by luck: adding `do_dwiprep_MNI` next to `do_dwiprep` this session was a near
miss, saved solely because that one grep happened to include a colon. 50 of the
53 distinct patterns had no colon at all.

All 66 reads are now `grep -E "^[[:space:]]*key:" $conf | head -n 1` with an
inline `#` comment stripped. Two silent failure modes go away with it:

| config | key | old | new |
|---|---|---|---|
| `do_dwiprep: 1` + `do_dwiprep_MNI: 0` | `do_dwiprep` | `1\n0` | `1` |
| `fmriprep_options: --skip-bids-validation   # note` | `fmriprep_options` | *empty* | `--skip-bids-validation` |

The second is the nastier one: annotating a setting made `grep -v \#` drop the
whole line, so the value silently reverted to its built-in default.

Verified by evaluating the old and the new expression for all 58 parsed
variables against every shipped config: identical values throughout, so nothing
that works today changes.

### Smaller fixes

- **`VOFc` was unreachable.** `KUL_FWT_tracks_list.txt` listed `VOF_LT`/`VOF_RT`,
  but the recipes are `track_recipes_v2/VOFc_LT.txt`/`VOFc_RT.txt`. A name
  mismatch, not missing anatomy: `KUL_FWT_make_VOIs.sh` logged "no recipe file
  found ... skipping" and the bundle was never built. Renamed; every entry in
  that list now resolves to a recipe.
- **`do_dwiprep_mni` renamed to `do_dwiprep_MNI`**, matching the key
  `KUL_preproc_all.sh` already uses for the same operation. Two keys differing
  only in case, controlling the same script from two orchestrators, was a trap.
- **DTI-ALPS docs corrected.** They claimed `KUL_calc_DTIALPS.sh` is "not part
  of KUL_NIS_unified" and "must be on the PATH". It ships in `KUL_DTI_ALPS/`
  and is invoked by explicit path. Also recorded that it does not depend on
  `KUL_dwiprep_MNI.sh`, which the step ordering wrongly implies.

## Unreleased (2026-09-03 — the dwi distortion-correction scheme is detected, not assumed)

`synbzero_disco_instead_of_topup` and `rev_phase_for_topup_only` describe the
*acquisition*, but the templates pinned them per processing type: types 1-4
shipped `1`/`1`. On a full AP/PA pair — what the current Siemens dMRI protocol
acquires — both are wrong. `-b` means "use Synb0-DISCO **instead of** topup", so
it replaces a real fieldmap with a synthetic one when a real one was acquired,
and `-r` then discards the reverse-phase volumes as data. Half the diffusion
data thrown away, silently, on every type-1-4 Siemens case.

Both keys now default to `auto`, resolved by `KUL_detect_dwi_acq` from the BIDS
dwi series after dcm2bids has run:

| BIDS holds | synb0 | revonly |
|---|---|---|
| one phase-encoding direction | 1 | 0 |
| both, reverse phase b0-only | 0 | 1 |
| both, reverse phase diffusion-weighted | 0 | 0 |

The resolved numbers are written into the `KUL_LOG/` copy, so `KUL_preproc_all.sh`
never sees the word `auto`. An explicit `0`/`1` in the config still wins.

If the scheme cannot be determined — no `PhaseEncodingDirection` in the sidecars,
or no dwi at all — the run **stops** with an error naming the two keys, rather
than guessing. Either guess corrupts the diffusion preprocessing of a clinical
case in a way nothing downstream would flag.

Verified against a real Siemens study (`ep2d_diff_b2500_1.5i_p2_s3_AP` / `_PA`):
detected `REV_FULL` -> `0`/`0`, independently reproducing the setting that study
had been hand-edited to use. Synthetic trees cover the single-direction,
b0-only-reverse, no-dwi and missing-`PhaseEncodingDirection` cases.

### `sequences.txt` covers Siemens as well as Philips

The template only ever matched Philips series names. Merged in the strings from
a working Siemens study, grouped by scanner: `n_t1`, `c_t1`, `flair`, `3D_t2`,
the `fMRI_HAND`/`fMRI_FOOT`/`TAAL_EN`/`TAAL_DE` task runs, and `s3_AP`/`s3_PA`
for the reverse-phase dMRI pair. Later generalized again to cover the DBS and
DTI-ALPS protocols, retiring their per-type copies: 46 rows to 60, nothing
dropped except `FLAIR,3D_FLAIR`, which collides with the broader `FLAIR,FLAIR`.

Siemens rows carry `-` for mb and pe_dir, which is what those columns are worth
there: `KUL_dcm2bids.sh` reads DICOM `(0008,0070)` and, on Siemens, takes both
from the dcm2niix sidecar and never reads the config values.

The dMRI rows match `s3_AP`/`s3_PA` rather than bare `AP`/`PA`. Unanchored, `AP`
also matches the Philips `dMRI_Linear_AP` row, and a series matching two entries
is dropped by dcm2bids ("Several Pairing") rather than placed — the scan would
be lost with only a warning. Simulated the fnmatch of every dwi row against both
scanners' descriptions: each series resolves to exactly one row.

The header comment also documented the wrong script — it described
`KUL_dcm2bids.py`'s `grep`, not `KUL_dcm2bids.sh`'s dcm2bids/fnmatch matching,
and omitted the `^`/`$` anchors that are supported.

## Unreleased (2026-09-03 — every sub-command gets a config file)

`KUL_clinical_fmridti.sh` had grown to ~30 flags, of which 13 existed only to
forward one value to one sub-script — `-A` is named that way because FWT's `-R`
collided with "generate DICOMs" at this level, which is the point where flag
space had clearly run out. Meanwhile `KUL_VBG.sh` and KUL_FWT, the two steps
whose settings matter most, had **no** config surface at all: `-z T1 -b -B 1 -t
-P 1 -M -O -H` and `-T 1 -a iFOD2 -f 1 -S` were literals in the middle of the
script.

The existing `run_dwiprep.txt` / `run_fmriprep.txt` / `run_freesurfer.txt` were
already per-sub-command configs, complete with their own `do_<step>:` toggles.
This extends that pattern to the rest of the pipeline instead of inventing a
new one.

### `study_config/_base/` + per-type overrides

`KUL_scaffold` now copies `_base/` first and the `-t` directory over it, so a
type directory holds only the files it genuinely changes. All four type
directories previously carried a full copy of every config, and they had
drifted: `run_freesurfer.txt` existed in three versions whose only real
difference was a hardcoded `freesurfer_ncpu` (32 vs 24). 40 duplicated files
became 21 in `_base/` plus 7 real overrides, once the later passes in this
same batch removed the copies that differed only in a hardcoded core count.

Overriding is per **file**, not per line — a type directory shipping
`run_fmriprep.txt` replaces the base one wholesale.

### `-n` now reaches every step

`fmriprep_ncpu: 64` / `freesurfer_ncpu: 32` / `dwiprep_ncpu: 64` were baked into
the templates, so `-n 32` on a 32-core box still asked fmriprep for 64. A new
`KUL_prepare_step_config` writes both the participant and `$ncpu` into the copy
placed in `KUL_LOG/`. It also anchors the participant substitution: the old
`s/BIDS_participants: /...` appended to whatever the template already held, so a
template shipping `BIDS_participants: 001` produced `BIDS_participants: P001`.

### Seven new config files

`run_vbg.txt`, `run_fwt.txt`, `run_fmri_glm.txt`, `run_rsfmri_networks.txt`,
`run_dsc.txt`, `run_multiparc.txt`, `run_karawun.txt` — each owning one
sub-command, each key mapping to one flag of the script it drives.

Decisions that follow from `-t` are deliberately **not** keys, since `-t` would
silently win: whether VBG runs and whether it is extra-axial, whether FWT runs,
whether multiparc runs, and which Karawun case applies. The files say so.

Read by `KUL_read_config`, anchored on `^key:` — `KUL_preproc_all.sh` uses an
unanchored `grep key $conf`, which also matches a comment mentioning the key.
A missing file or missing key falls back to the built-in default, so patient
folders scaffolded before this change keep working.

### 13 flags retired, `-x` added

`-S -P -E -N -C -W -X -U -Q -A` became config keys; `-f -y -m` became the
`$KUL_SCILPY_ENV` / `$KUL_PYFMRI_ENV` / `$KUL_DICOM_ENV` variables they already
defaulted to. Passing one now **errors** naming its replacement rather than
being silently ignored — a run that quietly dropped `-S 6` would produce wrongly
smoothed results with nothing to show for it.

`-x key=value` (repeatable) overrides any key for one run without editing the
study's config, and the overrides are echoed at startup.

Surviving flags: `-p -t -d -n -v -s -r -R -F -B -D -O -e -T -a -x`. The 11
self-invocations pass only `-p -t -F -O`, all of which survive, so they were
untouched.

### Verified

- Generated option strings compared against the old hardcoded ones with the
  shipped defaults: `KUL_VBG.sh` and `KUL_FWT_make_TCKs.sh` both come out
  flag-for-flag identical.
- Scaffolding `-t 1`, `-t 5`, `-t 7`: 21 files each, type 5 gets the DBS
  `sequences.txt` and DRT `tracks_list.txt`, type 7 gets `do_freesurfer: 0`.
- Config reads: a value edited in the patient's file is honoured; `-x` beats the
  file; a deleted key falls back to its default without error.
- Retired flags error with their replacement; `-x` without `=` is rejected.
- Not run against real patient data — the next real run is what confirms no
  regression in the steps themselves.

## Unreleased (2026-09-02 — `-t` is validated up front)

`KUL_clinical_fmridti.sh` never checked the value of `-t`, and every consumer of
it is an `if`/`elif` chain with no `else`: the scaffold template dispatch in
`KUL_scaffold`, the lesion-type dispatch that sets `hdglio`/`vbg`/`multiparc`/
`fwt`/`alps`, and the Karawun mapping. An out-of-range type therefore failed
silently and late instead of loudly and early.

- `-t 8` matched no scaffold branch, so `KUL_scaffold` created an **empty**
  `study_config/` and `exit 0`'d as though it had worked. The next run then died
  on the unrelated `-D` pre-flight (`run_dwiprep.txt does not exist`), which
  points at the wrong thing entirely.
- The same out-of-range value left `hdglio`/`vbg`/`alps` unset in the lesion-type
  dispatch, surfacing much later as bash integer errors.
- A non-numeric value (a typo'd flag swallowing the next word) hit
  `[: too many arguments` at the first arithmetic test.

Now validated once, immediately after the `-E` engine check and before
`KUL_scaffold` can run: `-t` must be an integer in 1–7, otherwise the script
prints the value it got plus the list of the seven types and exits 2. Every
downstream chain can keep assuming 1–7, so none of them needed an `else`.

Verified by running the script with `-t 8`, `-t 0`, `-t abc` and `-t ''` (all
exit 2, nothing created) and `-t 3` (scaffolds `clinical_fmri_dmri` as before).

## Unreleased (2026-08-14 — the results layer rebuilds itself, and Karawun prep actually runs)

The theme running through this batch: `RESULTS/`, `Karawun/` and `REPORT/` are
*products*, exported from analyses whose own output lives in `BIDS/derivatives/`.
Several of those exports were trapped inside "did this step just run" branches, so
once a step was marked done its products could never be rebuilt — deleting a
folder meant re-analysing, or in the worst case meant a run that skipped
everything, reported success, and left the folder empty.

The rule now holds everywhere: **clearing a step's marker re-exports; deleting its
derivative re-analyses.**

### `conda activate` worked only in shells that happened to have conda on PATH

`KUL_activate_conda_env` used `source activate <env>` — the pre-4.4 form, which
is a *script* that exists only once conda's `bin/` is on `PATH`. An interactive
login shell carrying the KUL_Linux_setup bashrc block has that; `cron`, `nohup`,
`bash -c` and any non-login shell do not. There it failed with

```
KUL_main_functions.sh: line 330: activate: No such file or directory
```

and the caller's own package check then reported it as **"conda env 'pyfMRI' is
missing required packages"** — pointing at the env, which was present and
complete, instead of at `PATH`. That misdiagnosis cost real debugging time.

`conda activate` needs conda's shell *function*, not its binary, so a new
`KUL_conda_bootstrap` sources `etc/profile.d/conda.sh` whenever that function is
absent, trying `$KUL_CONDA_BASE`, the base reported by any `conda` already on
`PATH`, then the usual install locations. Verified working in a shell started
with `env -i PATH=/usr/bin:/bin` — no conda anywhere. Activation failures are now
reported as activation failures. Shared by all eight scripts that call it.

### The fMRI python env had no pre-flight check

`lore_sd` and `scilpy` are checked before anything expensive runs; `pyfMRI` —
which the GLM, melodic and the `-N` networks all need — was not, so a broken env
surfaced only after fmriprep, dwiprep, VBG and FWT had already been walked.

Both existing checks test `conda env list | grep -qx <name>`, and that test would
have *passed* the failure above: the env existed. So `KUL_check_pyfmri_env` does
what the step does — activate, then `import nilearn, nibabel, numpy, pandas` — in
a subshell so a successful activation does not leak. It runs after
`KUL_check_data` (where `n_fMRI` is known), is skipped for `-R`/`-F`, and
distinguishes the three cases: env missing, env present but unusable, and conda
not initialised at all — each with the command that fixes it.

### A deleted RESULTS folder used to be undetectable

Every processing step gates on its `.done` marker in `KUL_LOG/` and never checks
whether its own output is still present. `KUL_register_anatomical_images` is
representative: marker present → `echo "Anatomical registration already done"`,
and it copies nothing.

Since the markers live in `KUL_LOG/` and the outputs in `RESULTS/`, deleting a
results folder to "start fresh" produced the worst possible state — every step
skipped, the run reported success, and RESULTS stayed empty. It also cascaded:
`KUL_karawun_prepare.sh` reads `RESULTS/sub-*/Anat/T1w.nii.gz`, so Karawun prep
then failed for a reason that pointed at the wrong place.

`KUL_verify_results` now runs alongside `KUL_check_redo`, before any processing.
The marker is what distinguishes "deleted" from "never generated", and it is a
reliable discriminator because every step puts its applicability test *outside*
the marker check (`KUL_fmriproc` only reaches `SPM.done` when `n_fMRI>0`; the DSC
block only when `n_dsc>0`) and only touches the marker on success:

| marker | output | meaning |
|---|---|---|
| absent | absent | never generated — normal, silent |
| present | present | fine |
| present | **absent** | it existed once and does not now — reported |

Because it needs no `-t` flags and no BIDS scanning, it cannot drift out of step
with the pipeline's own decisions about what to run.

Outputs with fixed names are checked directly (`Anat/T1w.nii.gz`,
`Anat/cT1w_T1w_subtracted.nii.gz`); those named after tasks, networks or bundles
are checked as "directory exists and is non-empty" (`SPM/`, `Melodic/`,
`Perfusion/`). `Tracto/` needs no check — it is re-synced from FWT on every run.

**Repair is only ever "remove the marker".** The function runs before the
processing steps, so the pipeline's own code regenerates the output moments later
in the same invocation, applying the same transforms in the same order — the
Tracto maps regrid onto `Anat/T1w.nii.gz`, the GM mask is recomputed from
fmriprep's `dseg`. Nothing here re-implements a processing step, so nothing here
can drift from one. On a healthy tree it prints nothing and asks nothing.

Three deliberate exceptions:

- **Type 3's lesion mask is a hard stop, not a prompt.** `RESULTS/sub-*/Lesion/lesion.nii.gz`
  for `-t 3` is hand-drawn and copied in by the user — no marker, no derivative
  copy, nothing that can regenerate it. RESULTS is therefore not purely an output
  directory; it has one irreplaceable input in it. Missing → `exit 2`, because
  continuing gives a silently lesion-free run: no VBG lesion, no Karawun lesion
  label, no PACS overlay, no error.
- **`Anat/T1w_GM.nii.gz` is reported, not repaired.** It comes from the fmriprep
  step, which gates on `fmriprep/sub-*.html`, so the only lever that would rebuild
  it is re-running fmriprep — hours, for a file that takes seconds. The three-line
  `mrcalc`/`maskfilter` recipe is printed instead.
- **No tty, no prompts.** Under `nohup`/batch, `read` would consume the script's
  stdin or block, so a non-interactive run reports and changes nothing, printing
  the `rm -f` commands to run by hand.

Detection runs unconditionally rather than under `-r`, deliberately: the failure
mode is not knowing anything is wrong, so a check you have to opt into would
never fire when it is needed. The *prompting* defers to `-r` entirely — that mode
is the explicit "ask me about redoing things" conversation, and it now has a
question for every check, so asking here too would ask everything twice.

### `-r` was missing questions for half the pipeline

`KUL_check_redo` had no prompt for anatomical registration, the Gd contrast
subtraction, DSC perfusion or the Karawun folder — so a `-r` run could not redo
them at all, and the integrity check above was their only offer. All four added,
following the existing pattern (clear the marker, and the derivative where one
exists). The Karawun entry deliberately clears only the marker: deleting
`Karawun/sub-*/` would take the donor DICOM in `Karawun/sub-*/DICOM/` with it, and
prep overwrites its own outputs with `-force` anyway.

### REPORT/ had no fMRI content once a subject had been processed

The GLM activation figures (`05_afMRI_*`) and the melodic network figures
(`05_rsfMRI_*`) were generated *inside* the `if [ ! -f …SPM.done ]` and
`if [ ! -f …melodic.done ]` branches — so they existed only in the run that first
computed the maps. Afterwards the maps were still on disk and the markers said
"done", and nothing would ever redraw them: deleting `REPORT/` lost every fMRI
visual permanently, which is how a subject ended up with a REPORT holding only
eddy QC, the fmriprep link and the tract summary.

Both loops moved out of those branches, and each figure is now drawn only when it
is missing — so a healthy tree renders nothing and a wiped `REPORT/` refills on
the next run. The melodic guard matches on the filename stem, since the threshold
is baked into the name.

Two reports that were being computed and then never surfaced are now copied in as
well:

- **melodic's own per-run HTML report** (`stats_*/report/00index.html`), which
  never left the derivative — symlinked, so it stays in step with the
  decomposition;
- **the rsfMRI-networks report** (`-N`), which was written into its derivative and
  into `RESULTS/rsfMRI_Networks/` but not into the one place a clinician looks.

### REPORT/ and the fmriprep exports were trapped behind expensive gates

Same bug class as below, in the two foundational steps. `KUL_run_fmriprep` gates
on `fmriprep/sub-*.html`, and its exports sat *inside* that branch — so once
fmriprep had run, `Anat/T1w_fmriprep.nii.gz`, `Anat/T1w_GM.nii.gz` and the
`REPORT/sub-*_03_fmriprep.html` symlink were unreachable. Restoring a **symlink**
cost a full fmriprep re-run. `KUL_run_dwiprep_anat` had the same shape around its
three QA copies (`02_eddy_qc.pdf`, the two FA overlays).

All of these are pure copies from sources that survive in `fmriprep/` and
`dwiprep/`, so they now run whenever their source exists rather than only when
the step does. The GM mask is the one exception: it is rebuilt only when absent,
since it is the only one that costs anything (two `maskfilter` passes) and
nothing upstream of it changes between runs. The report symlink is re-pointed
with `ln -sfn` rather than left stale.

This also removes the last item the integrity check could only advise on: it no
longer prints a `mrcalc` recipe for `T1w_GM.nii.gz`, because the run rebuilds it.

Deleting `REPORT/` is now recoverable from a normal run for everything except the
VBG montage and the tumour-segmentation PNG, which are genuinely derived, and the
`-F`/`-R` figures, which come back with `-F`/`-R`.

### Analyses now re-export instead of re-running

`RESULTS/` is the last layer of a stack: fmriprep/dwiprep/VBG/FreeSurfer are
foundational, the GLM, melodic, DSC, rsfMRI and FWT sit on top of those, and
`RESULTS/` + `Karawun/` are the clinical products exported from *them*. The
outputs of that middle layer live in `BIDS/derivatives/`, which survives anything
done to `RESULTS/` — but two of them deleted their own work before redoing it, so
there was no way to rebuild the product layer without re-analysing:

- `KUL_fmriproc_nilearn_new.sh` ran `rm -rf "$fmriresults"` before every GLM, so
  the existing `stats_*/spmT_0001*.nii` were destroyed and the warp into
  `SPM_all/` only ever saw freshly computed maps.
- `KUL_fmriproc_conn.sh` re-ran `melodic` unconditionally, though its
  decomposition and `kul_networks.txt` were sitting in the derivative and the
  network matching + warp below are cheap.

Both now re-use an existing result and skip straight to the export, gated on the
file the export itself reads (`spmT_0001.nii`, `stats/thresh_zstat*.nii.gz`).
`KUL_NILEARN_FORCE=1` / `KUL_MELODIC_FORCE=1` recompute anyway. This mirrors what
`KUL_dsc_fit` already did with its own `fit/rCBV_corrected.nii.gz`, so all three
analyses now behave the same way: clearing the marker re-exports, and only
deleting the derivative re-analyses.

### New: "Repopulate RESULTS and Karawun?" as an explicit `-r` step

The other `-r` questions mean "I was unhappy with this analysis, run it again".
There was no way to say the different thing — "the analyses are fine, but
`RESULTS/` does not reflect them" — and no clear answer to when RESULTS actually
gets repopulated, which is genuinely hard to infer: some of it is copied from
derivatives, some computed from other RESULTS files, and some written directly by
a sub-script.

Asked last, it rebuilds the product layer from analyses that already ran, and
decides what it can rebuild by **checking whether each step's input still
exists**. Where the derivative is there, the marker is cleared and the step
re-exports; where it is gone, the marker is left alone and reported, because that
would be a real re-analysis and belongs to its own question. It prints the split
rather than leaving it to be guessed:

```
    will rebuild:
      Anat/ (T1w, *reg2_T1w)       (from anat_reg)
      SPM/ and SPM_all/            (from SPM)
      Perfusion/                   (from DSC)
      Karawun/ (T1w, tck, labels)  (from karawun_prepare)
    cannot rebuild - the analysis output is gone, so these would be
    real re-runs; use their own questions above:
      Melodic/                     (melodic)
    never touched: Lesion/, PACS_input/, DICOM/ - your own files
                   Tracto/ and TRK/ re-sync from FWT every run anyway
```

It never deletes a derivative — that is what keeps it an export rather than a
re-analysis.

The `anat_reg` distinction is the point of it: this clears the marker but leaves
`KUL_anat_register_rigid/` in place, so the step re-copies. "Redo: anatomical
registration?" above deletes the derivative too and pays for ANTs again. Same
marker, two intents, now separately reachable.

### Karawun prep was unreachable code

`KUL_clinical_fmridti.sh` never ran `KUL_karawun_prepare.sh`. The step that was
meant to (STEP 15) was gated on `make_dcm -eq 1`, and `-R` is the only flag that
sets `make_dcm`. But `-R` also sets `results>0`, which enters the render/export
block near the top of the file and ends there in an unconditional `exit` —
~1200 lines before STEP 15 was reached. So:

- a plain run had `make_dcm=0` and printed "Karawun folder prep skipped (run with -R to prepare it)";
- an `-R` run exited long before the step;
- `-F` did both.

The only invocation that could have reached it was `-R 0`, which is not a valid
underlay. The consequence was not just a missing folder: the `-R` block writes
the per-task fMRI activation labels (palette values 51-63) only when
`Karawun/sub-*/T1w.nii.gz` already exists, which is what Karawun prep produces.
On a clean machine that file never existed, so the Brainlab export shipped
tracts and lesion labels but **no fMRI activations**, silently.

Karawun prep now runs at the end of every pipeline run, where the ordering
already works: it needs the FWT output from the step above it, and `-R` needs
the `T1w.nii.gz` it writes. It depends on nothing `-R` produces.

Two failure modes it used to have are fixed at the same time:

- **No FWT output.** Previously irrelevant (the step never ran); now that it is
  unconditional it would fail at the end of every run for a subject whose FWT
  produced nothing. It skips with a message and writes no marker, so it retries.
- **Unmapped processing type.** The `type<5 / =5 / =6` chain had no `else`, so
  type 7 (DTI-ALPS) fell through every branch: the return code stayed 0, the
  `.done` marker was written, and the run reported a successful prep having done
  nothing. Unmapped types now say so and leave no marker.

### PACS drop folders were created too late to use

`RESULTS/sub-*/PACS_input/{overlays,series_quantitative}/` are where an operator
copies the maps they want exported. They were created only inside the `-R`/`-F`
block — that is, only once the export run had already started and it was too
late to put anything in them. First-time users had to run `-R` once to have the
folders appear, then run it again to use them.

Creation moved into `KUL_make_pacs_dropdirs`, called both from the `-R` block and
from the end of a normal pipeline run, so the folders (and their `README.txt`)
are in place during the review step where they are needed. The function is
defined above the `if [ $results -gt 0 ]` block, alongside `KUL_resolve_lesion`
and `KUL_copy_lesion_to_anat`, for the reason documented there: that block is
top-level code, bash registers functions when their definition executes, and a
definition further down the file is simply not visible to it.

The README now also states what automatic discovery does *not* pick up — only
`nrCBV_corrected`/`nrCBF` from `Perfusion/`, so a raw `rCBV_corrected` or an ADC
map has to be dropped in by hand.

### PACS export no longer requires a full study directory

`-R`/`-F` render and export what is already in `RESULTS/`; they run no
preprocessing. But four pre-flight checks stood in front of them, every one
guarding a preprocessing step those runs never reach:

| check | effect on `-R` |
|---|---|
| missing `study_config/` | scaffolded a study directory that was not asked for |
| missing dwiprep config | `exit 2` |
| missing `lore_sd` conda env | `exit 2` |
| missing `scilpy` conda env | `exit 2` |

Exporting a map to PACS therefore required a valid dwiprep config and
preprocessing conda envs. These checks are now skipped when `results>0`, so a
directory holding nothing but `RESULTS/sub-<p>/` can be exported:

```bash
cd /path/to/anything/with/RESULTS
KUL_clinical_fmridti.sh -p <p> -t 1 -R 4 -O SAG
```

Narrow by construction: it skips checks only, never a processing decision, and a
normal run still fails fast on all four exactly as before.

### Closing banner rewritten

The end-of-run instructions said "-R also triggers Karawun prep (it no longer
runs on its own either)" — which was false, and is now the opposite of true. The
banner now reflects the real sequence: review (including the already-prepared
Karawun folder) → copy maps into `PACS_input/` → `-R` → `importTractography`.

## Unreleased (2026-08-12 — measurable DICOMs, and why sagittal export kept failing)

### Sagittal DICOM export died on any multi-frame donor

`KUL_nii2dcm.py`'s donor-match path (`-M`, used for SAG and only SAG) read the
donor with `sitk.ReadImage`. A single **multi-frame** DICOM — what Philips
exports, e.g. a 100-frame SmartBrain localiser in one file — comes back as a
**4-D** image `(cols, rows, frames, 1)`, so `GetDirection()` is a 4×4 matrix. The
code extracted the row direction as `(d[0], d[3], d[6])`, correct for 3×3 but on
a 4×4 giving `(0, 0, 0)`; normalising that raised `ZeroDivisionError`. TRA and
COR take a different path, so this presented as "sagittal DICOMs fail, the rest
are fine". The donor is now forced to exactly three dimensions first, and a
degenerate direction matrix reports what is wrong instead of a bare traceback.

### Every mrview on the host shared one Qt semaphore — the cause of the hangs

Qt creates a `QSystemSemaphore` whose SysV key is `ftok()`'d from a backing file
it writes into `$TMPDIR`:

```
openat("$TMPDIR/qipc_systemsem_bbadddcbefcffa02624b2d14c203d71b03f8027fd9e...",
       O_RDWR|O_CREAT|O_EXCL) = -1 EEXIST
semget(0x5102008e, 1, IPC_CREAT|IPC_EXCL|0600)
semtimedop(..., sem_op=-1, SEM_UNDO)      <- lock
```

Every mrview hashes to the same filename, so **every instance on the machine
shared a single semaphore** and serialised on it. Whenever one stalled at the
head of that queue, all the others blocked in `semtimedop` with **zero CPU** —
including instances started by hand from an unrelated terminal, which is how this
was spotted. Symptom: `cat /proc/<pid>/wchan` reads `do_semtimedop` and `ipcs -s`
shows waiters (`ncount > 0`).

Each render now gets a private `TMPDIR`, so each gets its own key and no two
mrview processes can be coupled. Verified: with a shared TMPDIR two concurrent
renders both used `0x5102008e`/`0x5102008f`; with private ones they used
`0x51028bef`/`0x51028bf1` and `0x51028bee`/`0x51028bf3`. On a real subject run,
concurrent renders now both accumulate CPU (`cpu=00:00:44`, `wchan=0`) where
previously one worked and the rest sat at `cpu=00:00:00`.

Render parallelism is also now overridable with `KUL_RENDER_PAR` (default 3) —
`KUL_RENDER_PAR=1` serialises, which was the workaround before this fix.

**Not fully characterised:** what stalls the process at the *head* of the queue.
A SIGKILL mid-render does not do it — `SEM_UNDO` cleans up correctly, verified.
The coupling is what turned one stalled process into a machine-wide freeze, and
that coupling is gone; a single stalled render can now only affect itself, and is
caught by the existing timeout and retry.

### One wedged mrview poisoned every later render on the host

A hung `mrview` holds a SysV semaphore it never releases, and every subsequent
`mrview` blocks on it indefinitely — renders producing no PNGs, no error, just a
timeout. Observed live: three processes wedged for 8.5 h, after which a render
that takes **0.53 s** on a clean machine was exceeding a 600 s timeout. `-R` now
warns when other `mrview` processes are present and prints the commands to clear
them (it does not kill them automatically — on a shared box they may be someone
else's live run).

The preflight is now a real one-slice capture rather than `mrview --version`,
because `--version` never opens a window or touches OpenGL and so passes happily
on a host where every actual render hangs.

### The automatic threshold blanked physiological maps entirely

`max/3` assumes the maximum is a meaningful peak. That holds for a statistical
map and fails for a physiological one, whose maximum is a vessel voxel. On real
data, rCBV had max 12111 against a p99 of 1543, so the automatic
threshold landed at 4037 and kept **704 of 2,083,162 voxels — 0.03%**. The
overlay was, in the operator's words, "almost entirely clipped away, only a few
voxels here and there".

When the maximum is that far out in the tail, thresholding is the wrong question:
these are brain-masked physiological maps meant to be read brain-wide. Maps whose
`max/p99` exceeds 3 now render continuous — full range, robust window, colourbar
— instead of being auto-thresholded. Measured on the same three slices:

| | coloured-pixel share |
|---|---|
| auto `max/3` = 4037 | 0.007% / 0.045% / 0.019% |
| continuous, windowed 3.2–1272 | 21.2% / 21.4% / 21.3% |

Statistical maps are untouched: the two fMRI maps have ratios 2.6 and 1.5,
below the trigger, and keep `max/3` exactly. Any explicit threshold — `-T`, a
`<name>.thresh` sidecar, or the fixed clinical cutoffs for nrCBV/nrCBF — still
wins outright, so this only affects the automatic case. That is also why it went
unnoticed: answering the interactive threshold prompt bypasses it, and only
unattended runs fall through to `max/3`.

### The default window on quantitative series was computed in the wrong units

Window Center and Width are applied *after* the Modality LUT (DICOM PS3.3
C.11.2), so they live in rescaled real-world units. They were being computed
from the stored int16 values, which made the error scale with each map's slope:

| map | slope | window ÷ data range | how it looked |
|---|---|---|---|
| rCBF | 9.7e-05 | 1522× | flat |
| K2 | 2.1e-05 | 1128× | flat, negative lobe invisible |
| TT0 / TTP | ~3e-03 | 221× / 148× | flat |
| K1 | 7.6e-03 | 6.4× | flat-ish |
| MTT | 0.105 | 0.002× | far too narrow, blown out |
| rCBV_corrected | 0.803 | 0.1× | too narrow, too dark |

Only maps whose slope happened to sit near 1 looked approximately right, which
is why this presented as "sometimes too bright, more often too dark" rather than
as an obvious failure. Now computed in real units.

Two related fixes to the same window:

- **Foreground only.** Percentiles had been taken over the whole volume, but
  these maps are brain-masked and ~87% exact zero, so the window collapsed onto
  the mask rather than the tissue.
- **Top at p99.5, not p98.** On a DSC map the choroid plexus (genuinely very
  vascular) and CSF (deconvolution garbage) reach ~46× the cortical median, and
  a p98 top left the ventricles as blown-out white blobs. p99.5 cuts saturated
  tissue from 2.0% to 0.5%, at the cost of cortex sitting at ~14% of the scale
  rather than ~21% — the right trade for series meant to be re-windowed on PACS.

Label maps span their full range instead, since labels are categorical.

CBV maps additionally get 30% extra headroom above p99.5
(`--window-headroom`, applied by the pipeline to `*CBV*` only), which drops
their saturated tissue from 0.50% to 0.20%. This is keyed on the map rather than
inferred, because no statistic separates the maps that want it: on real DSC data
K2 (max/p99.5 = 47) and MTT (260) have far heavier tails than rCBV (6.5), so any
tail-based rule would widen exactly the wrong windows and flatten them.

### Colourbar ends are rounded

Overlay intensity ranges are rounded so the colourbar reads cleanly: whole
numbers once values reach 1, at most 3 decimals below that. rCBV's window prints
as `3,1054` rather than `3.23908,1053.84`, rCBF as `0.005,0.322`. Reverted
automatically if rounding would collapse the range.

### Overlay windows are anchored to cortex, not to the map's own percentiles

A percentile window drifts with how much tumour and vessel happen to sit in the
field of view, so identical physiology gets a different colour in different
patients. Continuous overlays now take the top of the window from the map's
median inside grey matter (`Anat/T1w_GM.nii.gz`), times 4, which makes the scale
read as "x cortical value" and stay comparable between studies and scanners.
Falls back to the 98th percentile when no GM segmentation is available or it does
not share the map's grid.

On the validation subject, rCBV cortical median 263.5, so the window becomes [3.2, 1053.8] rather
than [3.2, 1272.4]. Rendered brain-wide at ~21% coverage with 0.05% saturating
the top of the scale. Thresholded statistical maps are unaffected — the cortex
median of a t-map is not a meaningful anchor, so those still window over their
suprathreshold voxels.

Tunable via `_win_ref_mult` (default 4) and `_gm_ref` in the `-R` block.

### Quantitative export had no automatic fallback

Measurable series were only produced from the `series_quantitative/` drop folder.
With nothing dropped there, every DSC map in `Perfusion/` is now exported
automatically (masks excluded) — 8 series on the validation subject: rCBV corrected/uncorrected,
rCBF, MTT, TTP, TT0, K1, K2. PCASL is deliberately not included; it is not
processed yet. Anything in the drop folder still wins.

### Overlays washed out whenever the map's range was wide

`-overlay.threshold_min` controls which voxels are drawn, not how values map to
colours. With no explicit `-overlay.intensity`, mrview windows the colourmap over
the volume's full range, so a perfusion map with values in the thousands rendered
as a saturated white blob while a narrow-range fMRI t-map happened to look right.
Every overlay now gets a robust 2–98 % window computed from the map itself
(restricted to suprathreshold voxels when a threshold applies), so wide- and
narrow-range maps are both legible. Applies to the automatic
SPM/Melodic/Perfusion/Lesion discovery too, not just the new drop folders.

### NIfTI → DICOM produced numbers nobody could measure

The NIfTI path scaled values to fill int16 and discarded the scale factor, never
wrote `StudyInstanceUID` (so output landed as an orphan study) or
`SOPInstanceUID`, and hardcoded Secondary Capture — an IOD with no Modality LUT
module. New `-q/--quantitative` writes MR Image Storage with
RescaleSlope/Intercept, so an ROI on PACS reads real units; `--label` preserves
integer labels unscaled. Verified against the DSC phantom's ground truth: lesion
ROI mean **41.32314** on the exported DICOM vs **41.3231** in
`sub-TEST_perfusion_stats.tsv`.

The int16 conversion was also unsafe — numpy wraps modulo 2¹⁶ rather than
clipping, so strongly negative voxels came back as large positives; values are
now rounded and clipped, NaN is handled, and all-negative input no longer
inverts. Direction cosines were emitted with `str()`, which can exceed DS's
16-byte limit; all DS values are now formatted to fit.

### Also

- **Drop folders.** `RESULTS/sub-*/PACS_input/{overlays,series_quantitative}/`
  export anything copied in, no renaming or processing. Automatic
  SPM/Melodic/Perfusion/Lesion discovery runs unchanged while both are empty.
- **SeriesNumber** was empty on every exported series. Now
  `donor_SeriesNumber × 100 + map_index × 10 + orientation`, keeping the series
  tied to the study and ordered, clear of the range scanners use.
- **Failures are no longer swallowed.** Exit statuses are checked, per-series
  logs land in `KUL_LOG/sub-*_PACS/`, and the run ends with a summary and a
  non-zero exit instead of always reporting success.
- **Truncated renders are repairable.** The skip guard compared against "any PNG
  exists", so a run killed mid-render was cached as complete forever; it now
  compares the actual slice count and retries once.
- **`xvfb-run -a`** replaces hardcoded display numbers (`:20`, `:10`–`:12`) —
  the latter collide with `ssh -X`, which allocates from `:10`.
- **fMRI/Melodic/Clinical renders run 3-way parallel**, as the tract renders
  already did; they were serial while still being sized as if 3 were in flight.
- **`GTK_PATH` and friends are stripped** around mrview. Launched from a
  snap-packaged VS Code terminal, Qt loaded the snap's `libcanberra-gtk-module`,
  whose RPATH dragged in a conflicting glibc and killed mrview outright.
- **Dedicated `KUL_dicom` conda env** (`share/envs/KUL_dicom.yml`, `-m` to
  override) instead of whatever `python3` was on PATH.

## Unreleased (2026-08-11 — task-fMRI reaches melodic; three silent BIDS-entity misses)

### melodic was being handed data with the task already removed

`KUL_fmriproc_conn.sh` denoised before running melodic, using the resting-state
confound set and no task term. The nuisance regressors are routinely collinear
with a block paradigm, so this did not merely fail to help — it removed the
signal melodic exists to find. Measured on a 30 s on/off language run:

| confound block | share of the task design's variance |
|---|---|
| 6 motion parameters | 27.4 % |
| 10 aCompCor (WM+CSF) | 23.2 % |
| full 22-regressor set | **52.7 %** |

Task correlation in the top-500 task voxels fell 0.522 → 0.224, and melodic's
best component correlated 0.27 with the paradigm — i.e. nothing task-like
survived. This is not about how much the patient moved: `trans_x` correlated
0.307 with the design, and *task-locked* movement is maximally damaging however
small, because it is collinear.

The task GLM was never affected, and that is the clue: nilearn's
`FirstLevelModel` fits task, drift and confounds in **one** design matrix, so its
confound betas are estimated controlling for the task.

Task melodic now runs on unregressed data (`--no-confounds --lp 0`: high-pass and
smoothing only), which is what FSL's own FEAT/MELODIC feeds ICA — separating
signal from artifact is the job being handed to the ICA. Result: the task
component appears at |r| 0.63 (run-01) and 0.77 (run-02, **variance rank 2**),
and its spatial map correlates 0.57 with the subject's own GLM z-map with every
peak voxel inside the GLM's activation.

`KUL_fmri_denoise.sh` gains:

- `--no-confounds` — skip confound regression entirely (spike regressors
  included; they cost degrees of freedom too). `--lp 0` disables the low-pass.
- `--task-signal {ignore|preserve|remove}` — `remove` regresses the design out so
  a task run becomes genuine pseudo-rest and can be pooled with real rest;
  `ignore` (default) is the historical behaviour. Needs an events TSV; a task
  without one falls back to `ignore` with a warning, and rest runs say so quietly.
- `--events-dir` to point at the events TSVs.

`preserve` is implemented and documented but **deliberately unwired**: protecting
the design necessarily protects task-correlated motion too, which pushes a
quarter to a third of the brain above |r| 0.4. Orthogonalising confounds against
the design beforehand does not work either — nilearn band-passes confounds
*before* regressing (correctly, per Lindquist et al. 2018), which destroys the
orthogonality.

`KUL_run_rsfMRI_networks.sh` now denoises task runs with `--task-signal remove`
(override with `KUL_RSFMRI_TASK_SIGNAL=ignore`), and `KUL_fmriproc_conn.sh` skips
resting-state entirely — the rsFMRI pipeline's masked ICA already decomposes it,
so running melodic there as well was the same work twice. Outputs are
variant-keyed so the two differently denoised products cannot be mistaken for one
another.

### Three BIDS-entity misses, same root cause

fMRIPrep writes `_acq-`/`_dir-`/`_echo-` *between* `task-` and `run-`, so
`task-[A-Za-z0-9]+(_run-[0-9]+)?` matches only `task-TAAL` and silently drops the
run index. In `KUL_fmriproc_nilearn_new.sh` and `KUL_fmriproc_spm_new.sh` that
produced a confounds glob matching nothing; `awk` then fell back to stdin and
wrote a zero-byte confounds file, and the failure surfaced much later as pandas'
`No columns to parse from file` inside the GLM. The with-confounds GLM had
therefore **never produced output**, and because `RESULTS/SPM` was populated only
from that variant, it held nothing but the engine marker while every real map sat
in `SPM_all`. The per-run confounds filename collided between runs for the same
reason.

Both scripts now match the two entities separately and derive the confounds
sidecar from the BOLD name (everything from `_space-` on is fMRIPrep's
output-space decoration, `_res-<N>` included), which is inherently
resolution-agnostic. `KUL_tsv_filter` refuses to write rather than emitting an
empty file.

### `res-2` no longer assumed

Audited every `res-2` reference. `utils.py` and `step0_synthseg.sh` already
preferred `res-2` with a resolution-agnostic fallback; the one hard failure was
`KUL_fmriproc_conn.sh`'s `fslcc` call, which compares melodic's ICs against a
network atlas pre-resampled to the `res-2` grid. `fslcc` requires identical
grids, so a run without `--output-spaces ...:res-2` died with *"Mismatch in image
dimensions"*. The atlas is now resampled to whatever grid melodic actually
produced (nearest-neighbour, cached per grid) and only when they differ.

### `RESULTS/SPM` vs `SPM_all`

Every GLM variant is warped once into `SPM_all` (the complete record).
`RESULTS/SPM` — what is reviewed and exported — is populated after all GLMs
finish, by copying the preferred variant (`_wc`, falling back to plain). Stale
maps from a previous run are cleared first, so a run whose preferred variant
changed cannot leave both behind for the figure/DICOM loop to export. The
`engine_nilearn.txt` marker moved to `BIDS/derivatives/.../SPM/`, out of a folder
that is iterated over.

### KUL_FWT outputs

- **`.trk` for freeview**, written next to each `.tck` by `KUL_FWT_make_TCKs.sh`
  and collected into `RESULTS/.../TRK/`. The conversion reference is `subj_FA`,
  not the `-F` parcellation: KUL_FWT registers FS to FA itself and tracks in FA
  space, so an FS-space reference would write a `.trk` that renders offset
  whenever the dMRI was never aligned to the T1w (a separate session, say).
- **Endpoint connectivity in the per-bundle report** — a ranked parcel-pair table
  and a labelled heatmap restricted to the parcels the bundle actually touches,
  replacing an 89×89 imshow that was ~99.9 % zeros. Validated against known
  anatomy: CST → precentral↔brainstem 74 %, FAT → parsopercularis↔superiorfrontal
  97 %, MdLF → superiorparietal↔superiortemporal 68 %.
- **`KUL_FWT_bundle_report.py`** — one self-contained HTML with every bundle's
  screenshots, switchable between the four renderings and filterable by bundle,
  each linking to its detail page. Source PNGs total ~32 MB per subject; they are
  autocropped, downscaled and JPEG'd to ~14 MB so the page actually opens.
- The per-bundle HTML and the contact sheet are copied into `REPORT/`, flat, so
  the report's `metrics →` links resolve.

### Karawun

`sub-{participant}_karawun_prepare.done` was touched unconditionally, so a failed
prep still gated every later `-R` into printing "already prepared" over an
incomplete folder. It is now written only on success. Both donor-DICOM
directories (`Karawun/sub-*/DICOM/` and `RESULTS/sub-*/DICOM/`) are created up
front — the search already preferred the Karawun one, but the end-of-run text
named only the fallback and told the user to `mkdir` it. That text now also
explains what `-R` does *not* do (it never pushes to Brainlab; `importTractography`
is printed to run manually) and that repeating `-R` with a different underlay
regenerates the PACS DICOMs but deliberately not Karawun.

### SPM12 / MATLAB

`KUL_fmriproc_spm.sh` removed in favour of `KUL_fmriproc_spm_new.sh`. The
`share/spm12/*.m` templates resolve SPM via `$KUL_MATLAB_APPS` (falling back to
the legacy `$KUL_apps_DIR`) and now fail with a clear error when unset or when
`spm.m` is absent, instead of `addpath('/spm12')`-ing silently and dying later at
the first `spm()` call. The installer gains a `spm12` section that lays SPM12
down in `$SOFTWARE_ROOT/src/matlab_apps/spm12` and exports `$KUL_MATLAB_APPS`;
MATLAB itself is commercial and is not installed, and `-E nilearn` remains the
MATLAB-free path.

### `-B` now removes what it used to archive

The denoised BOLD copies and the SUSAN-smoothed GLM inputs are whole 4D series
rebuilt from the fmriprep output plus the confounds TSVs — pure intermediates —
but they were neither cleaned nor excluded, so `-B` packed them into the `.7z`.
The variant split above made it worse: the task-melodic, resting-state and
pseudo-rest paths each keep their own copy (~6 GB on a three-run subject, plus
~2 GB of smoothed inputs). They are now removed alongside `fmriprep_work` and the
dwiprep working dirs; on a real subject the cleanup drops ~43 GB of 140 GB.

One carve-out: dwiprep's `dwi/` is emptied *except* for `geomcorr.mif` and its
`geomcorr_grad_checked.b`. That volume is what `dwifslpreproc` produced — the
output of topup + eddy, the most expensive step in dwiprep — and its own scratch
(`*dwifsl*tmp*`, ~19 GB) is deleted, so without it a corrected DWI can only be
recovered by re-running eddy. `dwi_preproced.mif` one level up is the
post-bias-correction volume and does not substitute for it.

`-B`'s help text was a single line ("make a backup and cleanup") and now states
exactly what is removed, what is kept, and that it is destructive and one-way.

### Known gaps

- Rest runs are still denoised twice (`conn`'s `standard` and the rsFMRI
  pipeline's `task-remove` are identical for them). Real deduplication needs the
  rsfMRI pipeline's path resolution to become variant-aware.
- `share/nilearn/KUL_fmriproc_nilearn_new.sh` is a stale, non-executable
  duplicate of the live root copy (544 lines divergent) and still carries the old
  entity bug. The root copy is the one on `PATH`.
- `KUL_FWT_bundle_report.py` needs Pillow, which the `scilpy` env only has
  transitively (via matplotlib/fury/scikit-image) rather than as an explicit
  dependency.
- The Karawun changes are untested — no `-R` run has exercised them, and no
  `-B` run has exercised the new cleanup list.

## Unreleased (2026-08-09 — `use_native_dwi` now reaches the data it names)

`KUL_dwiprep.sh -u` / `use_native_dwi: 1` promised native-resolution processing
and only half delivered. The upsampling to 1.3 mm ran **unconditionally**; the
flag merely decided whether *estimation* read the result. So a native-mode run
produced FODs at the acquired resolution and everything else at 1.3 mm:
`dwi_mask.nii.gz`, `dwi_preproced_reg2T1w.mif`, the tensor built from it, and
therefore `qa/fa_reg2T1w.nii.gz` — which is the space KUL_FWT builds every VOI
in.

Nothing complained, because the two places that meet those grids both tolerate a
mismatch: `tckgen` interpolates its mask, and `mrtransform -linear` without
`-template` only rewrites the header. `voxel2fixel` does neither, which is why
KUL_FWT's tractometry died on every bundle (see the KUL_FWT changelog) and why
the LoRE DEC needed a resampled copy of the ODF to match a mask that should
never have been upsampled.

`dwi_preproced.mif` remains the canonical output and is still always produced —
what changes with `-u` is its resolution. Both modes are now internally
consistent end to end, so the fixel mismatch is impossible rather than patched at
the point it surfaces.

Also here:

- **A stale-file guard.** A native and an upsampled `dwi_preproced.mif` are
  indistinguishable by name, so re-running a subject with the flag flipped would
  silently reuse the wrong-resolution volume and every downstream grid would
  follow the file rather than the request. The voxel size is now checked, not
  just existence.
- **`odf_resampled.mif` removed.** It existed solely because the LoRE DEC read a
  hardcoded `dwi_mask.nii.gz` while its neighbours used `${dwi_mask_input}`; with
  the mask honouring `-u`, `odf.mif` can be masked directly. The DEC now matches
  the dhollander/tax/tournier calls and cannot break on the removed file.
- `-u`'s help text and the numbered description both claimed unconditional
  upsampling; corrected.

**Not yet exercised.** No dwiprep run has used this. Two known gaps:
`KUL_dwiprep_anat.sh` gates on its own outputs, so on an *existing* subject the
guard rebuilds `dwi_preproced.mif` while `dwi_preproced_reg2T1w`/`fa_reg2T1w`
stay stale — flipping the flag on a processed subject needs those cleared too.
And `KUL_dwiprep_group_fba.sh` still hardcodes `population_template -voxel_size
1.3`, which is a legitimate template-space choice but is the one remaining place
that assumes it.

## Unreleased (2026-08-09 — HD-GLIO-AUTO installed natively; tumour segmentation failed silently)

`KUL_anat_segment_tumor.sh` ran HD-GLIO through a local-install path
(`/usr/local/KUL_apps/HD-GLIO-AUTO/scripts/run.py`) left behind by the move to
`$SOFTWARE_ROOT`, then fell back to bare `hd-bet` / `hd_glio_predict` — neither
on PATH in any standard setup. Neither call was wrapped in `KUL_task_exec`, so
**the tumour segmentation produced nothing, with no log and no warning**, and the
run continued to build a lesion mask without the tumour in it. `resseg` run 2
then failed too: its input is derived from HD-GLIO's brain mask, not from run 1.

`setup_environment.sh` gains `section_env_hdglio`, a native (non-docker) install:
HD-GLIO-AUTO at `6acaad8`, HD-BET **1.0** at `5f63601` in its own env and
checkout (HD-GLIO-AUTO calls `hd-bet -device 0`, which HD-BET 2.x rejects — the
existing `hd-bet-env` is untouched), plus shared model weights under
`$SOFTWARE_ROOT/share/hd_models`.

Transplanting a 2020 image onto a current Python needed eleven fixes, of which
three produce **no install-time signal at all**: `numpy<2` (installs fine, then
`RuntimeError: Numpy is not available` at the first array↔tensor conversion), a
missing shebang on `HD_BET/hd-bet` (`Exec format error`), and an upstream bug in
`HD_BET/data_loading.py` that compares an array against a 3-element shape — which
numpy <1.25 silently evaluated as `False` and numpy ≥1.25 raises on, after
inference has completed. Three of the eleven patch **source files** in git
checkouts, so they are reapplied idempotently rather than assumed.

Verified end to end on a real 4-contrast clinical study: 4m34s,
`segmentation.nii.gz` + `volumes.txt`, volume cross-checked against an
independent count of the label map.

`KUL_anat_segment_tumor.sh` now locates the install by searching
`$SOFTWARE_ROOT`, `/opt/kul_software`, `/usr/local/KUL_apps` (verified to resolve
with `SOFTWARE_ROOT` unset), puts the env's `bin` on PATH only for that call,
wraps it in `KUL_task_exec`, exits with the install command if absent, and checks
`segmentation.nii.gz` exists afterwards — `run.py` writes it minutes before it
finishes, so a late crash otherwise leaves a directory that looks complete.

## Unreleased (2026-08-09 — dcm2bids: BIDS validation, and it now stays quiet)

The validator at the end of every conversion reported two real errors, both ours:

- **`ContrastBolusIngredient: "gadolinium"`** — that field is a controlled enum;
  the lowercase value failed schema validation on every post-contrast T1w this
  has ever produced. Now `GADOLINIUM`.
- **`tmp_dcm2bids/`** — `.bidsignore` covered its contents but not the directory
  itself. Moved under `sourcedata/`, which BIDS recognises and does not validate.
  Moved rather than deleted: it holds every series that did *not* match the
  config, which is exactly what you need when a scan is missing.

Plus the placeholder warnings: `BIDSVersion` had a stray leading `v`, `Name` was
empty, and the README was *appended* to on every conversion (stacking duplicate
paragraphs) rather than written.

`GeneratedBy` was deliberately **not** added despite the validator recommending
it: in BIDS that marks a dataset as a *derivative*, and derivative anatomicals
then require `SkullStripped` in every sidecar — three cosmetic warnings become
ten hard errors on what is correctly raw data.

**The validator no longer runs by default** (`KUL_BIDS_VALIDATE=1` to enable).
What remains on a correct conversion is metadata that simply is not in the
DICOMs, and printing it every time trains you to ignore the output — which is
worse than not running it, since real errors arrive the same way. The invocation
is kept in the file with the reasoning next to it. `docker run -ti` also became
`-i`: `-t` fails with "the input device is not a TTY" whenever this runs
non-interactively, which silently skipped validation inside the pipeline.

## Unreleased (2026-08-09 — melodic network comparison compared across templates)

`KUL_fmriproc_conn.sh` ran `fslcc` between `KUL_NIT_networks.nii.gz` (classic
FSL MNI152 2 mm, 91×109×91) and fMRIPrep's `melodic_IC.nii.gz`
(MNI152NLin2009cAsym res-2, 97×115×97). Different templates, and `fslcc`
requires matching grids, so it errored out. Now uses a pre-resampled
NLin2009cAsym copy of the atlas, generated once (not per subject) the same way
`step0_synthseg.sh` handles shared atlases.

## Unreleased (2026-08-08 — resseg has never had working weights)

`resseg` could not load its pretrained checkpoint, so **every resection-cavity
segmentation this install has ever attempted failed**, reported only as
`resseg run 1 might have failed` in a log while the pipeline carried on.

Two causes, both upstream. `pip install resseg` ships no checkpoint. And
`resseg/model.py` looks for it at `Path(__file__).parent.parent`, which assumed
`model.py` sat one level below the repo root beside the weights — untrue since
upstream moved to a `src/` layout (`src/resseg/model.py` → `parent.parent` is
`src/`, weights are at the repo root) and untrue for a pip install (→
`site-packages/`). It fails whichever copy of the package gets imported.

`setup_environment.sh` now fetches the checkpoint into `resseg/weights/` and
patches the lookup to search there first, falling back to the historical paths,
the torch.hub checkout, and finally a download. Idempotent, and it runs even
when the env already exists so existing broken installs are repaired.

**The checkpoint must never be placed directly in `site-packages/`**, which is
the obvious way to satisfy the upstream path. Python's `site` module parses
every `*.pth` file there as a UTF-8 path-configuration file at interpreter
startup, so a binary checkpoint makes that environment's python refuse to start
(`Fatal Python error: init_import_size: Failed to import the site module`).
Confirmed the hard way; noted in the code so nobody retries it.

The install verifier now instantiates the pretrained model
(`from resseg.model import ressegnet; ressegnet()`) instead of only doing
`import resseg, ants`, which passed happily on a broken install and is why this
went unnoticed. Checked both ways: it WARNs with the weights removed and the
upstream `model.py` restored, and passes once repaired.

Verified end to end on a real subject: the previously failing command produced
`T1_cavity1.nii.gz` (6648 voxels) in 16.5 s.

**Not fixed, and blocking the rest of that step:** `KUL_anat_segment_tumor.sh`
runs HD-GLIO via a local install path that does not exist here, then falls back
to bare `hd-bet` / `hd_glio_predict`, which are not on PATH (`hd-bet` lives only
inside the `hd-bet-env` conda env; `hd_glio_predict` is not installed at all).
Neither call is wrapped in `KUL_task_exec`, so this fails with **no log and no
warning**, leaving no tumour segmentation. `resseg` run 2 then fails too — its
input is built from HD-GLIO's brain mask, not from run 1's output. The
installer already pulls `jenspetersen/hd-glio-auto` (present on this machine),
but the script has no docker branch to use it.

## Unreleased (2026-08-08 — dcm2bids: anchored search-strings, ASL on Siemens, vendor detection)

Four bugs, all found on one Siemens MAGNETOM Cima.X (syngo MR XA61) study, and
all failing quietly enough to reach the end of a pipeline run unnoticed.

### Search-strings can now be anchored

The second config column matched anywhere in the `SeriesDescription`
(`t1_mprage` → `*t1_mprage*`). It can now be anchored: `^str`, `str$`, `^str$`.
Unanchored stays the default, so every existing config is unaffected.

This is not a convenience. A post-contrast series is very often named `c_` +
the pre-contrast name, so one search-string matches both — and when a series
matches two descriptions **dcm2bids places it nowhere**, logging `Several
Pairing` and leaving it in `tmp_dcm2bids/`. On the study at hand that meant no
post-contrast T1w at all, which `-t 1` needs, reported only as a warning inside
dcm2bids' own log.

The same scanner writes its console post-processing as sibling series
(`..._Brainmask`, `..._SS`, `..._SS_N4`, `..._SS_N4_Dn`), all tagged `ORIGINAL`
so no image-type filter separates them. `t2_space` matched all five, giving
`run-01` … `run-05` of which four were derivatives. `^3D_t2_space_sag_cs3_iso$`
now selects the acquisition alone.

The `*…*` wrapping had been hand-written at 19 call sites; it is now built once
per config line. One of those 19 had the leading `*` outside the quotes, where
it was subject to pathname expansion — fixed in passing.

### ASL converts on Siemens, and keeps the derived series

The ASL identifier pinned `ImageType` to `ORIGINAL\PRIMARY\PERFUSION\NONE`, a
Philips spelling. dcm2bids compares `ImageType` element-wise *and* requires
equal length, so a Siemens pCASL (`ORIGINAL\PRIMARY\ASL\NONE\MAGNITUDE`) failed
on the length check before one string was compared — and converted nothing,
without an error, because "no series matched" is not one. Now matched on
`SeriesDescription` alone, as DSC already was.

An ASL protocol also yields several series sharing a `ProtocolName`: raw
label/control plus the console's subtraction and rCBF maps. `acq_label` is now
honoured, so they can be separated (`_acq-raw_asl`, `_acq-deltam_asl`,
`_acq-cbf_asl`). The derived ones were additionally invisible to
`kul_find_relevant_dicom_file`, which filters on `ORIGINAL` — correct for
anatomicals, wrong for perfusion maps that are `DERIVED` by definition. ASL now
falls back to an unfiltered search; no other identifier's behaviour changes.

Still **not** BIDS-valid ASL: no `_aslcontext.tsv`, and `LabelingDuration`,
`BackgroundSuppression` and `M0Type` are absent. `perf/*asl*` stays in
`.bidsignore`. The per-volume `M0_SCAN`/`LABEL`/`CONTROL` labels do sit in
`ImageComments`, so the context file is tractable when someone needs it.

### Two ASL header traps, documented

Both found while scoping that work, both capable of producing quietly wrong
numbers rather than an error. Written up in the ASL section of
[the dcm2bids doc](docs/KUL_dcm2bids/KUL_dcm2bids.md).

**`PostLabelingDelay` in the sidecar may be the labeling duration.** The only
timing in a standard DICOM tag is `(0018,9258) ASLPulseTrainDuration`, which
the standard defines as the *labeling* pulse train duration; dcm2niix writes
`PostLabelingDelay`. On the study to hand both are 1800 ms — confirmed against
the vendor protocol (`sAsl.ulLabelingDuration`, `sAsl.sPostLabelingDelay[0]`) —
so the sidecar is accidentally right and the two candidate sources cannot be
told apart. Where labeling duration ≠ PLD, which is the usual case, the field
could be mislabelled. CBF scales with both. Unresolved pending a protocol where
the values differ; flagged rather than worked around.

**`ASLContext` (0018,9257) is unreliable.** Per-frame, and on this scanner
shifted by one with no `M0` state — it calls the proton-density M0 volume a
LABEL. `ImageComments` is self-consistent. Any future `_aslcontext.tsv` must
read `ImageComments` and validate the pattern, not trust the standard tag.

### Vendor detection recognises the Siemens XA line

`kul_dcmtags` tested `"$manufacturer" = "SIEMENS"`. XA scanners write `Siemens
Healthineers`, so every XA study was treated as Philips and printed `It's NOT
original dicom data (anonymised?)` once per series. Now a case-insensitive
prefix match. Effect was cosmetic — dcm2niix fills the sidecar regardless,
verified — but the message was alarming and wrong, and the Philips ees/trt
calculation was being attempted on data with no Philips tags.

Verified end to end on the Cima.X study: `ce-gadolinium_T1w` present, one `T2w`
instead of five, three `perf/*_asl` series, zero spurious vendor warnings.

## Unreleased (committed locally, 2026-08-08 — shard-recon: multiband factor and MRtrix3 incompatibility)

### `mb` is no longer hardcoded

`KUL_dwiprep.sh` hardcoded `mb=2` for shard-recon's slice-timing model
(`-mb ${mb} -sorder 1,0`), with its own comment noting it "needs to be turned
into a configurable parameter if shard-recon is used".

That matters more than a stale default usually would: **a wrong multiband factor
does not fail.** It tells dwimotioncorrect which slices were acquired
simultaneously, so an incorrect value yields a plausible-looking but wrong motion
correction, with nothing to indicate it. Of the two studies to hand, one is MB2
and the other MB3 — so the hardcoded value was already wrong for half of them.

Now read from the dwi BIDS sidecar's `MultibandAccelerationFactor` (present and
correct in both studies checked), overridable with the new `-M`. If neither
yields a value the run **stops** rather than defaulting, for the reason above.

### shard-recon is incompatible with current MRtrix3 — now said out loud

shard-recon targets the MRtrix3 **3.0.x** python API: its scripts end in
`mrtrix3.execute()`. The dev/CMake branch moved that entry point to
`mrtrix3.app._execute(usage, execute)`, so against the installed
3.0.8-2099-geba5dd55 every shard-recon command dies immediately with

```
AttributeError: module 'mrtrix3' has no attribute 'execute'
```

`shard-recon/bin/mrtrix3.py` is a symlink into the MRtrix3 install, so it picks
up whatever is there. This is not a PATH or build problem — the function it calls
no longer exists.

`KUL_dwiprep.sh -c` now probes `dwimotioncorrect` up front and exits with an
explanation, instead of failing at "part 3" after denoise, degibbs and topup have
already run. The fix it points to is a classic (non-CMake) MRtrix3 3.0.x with
shard-recon rebuilt against it — which shard-recon's own build requires anyway,
having no CMake equivalent.

Documented in `-c`/`-M` usage text.

## Unreleased (committed locally, 2026-08-08 — rsfMRI FreeSurfer lookup never resolved)

`step0_synthseg.sh` derives its default paths from its own location:

```bash
BASE_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
```

After relocation into `share/`, that resolves to the **KUL_NIS repo**, not the
study directory. `KUL_run_rsfMRI_networks.sh` exported `RSFMRI_FMRIPREP_DIR`,
`RSFMRI_DENOISED_DIR`, `RSFMRI_ANALYSIS_DIR` and `RSFMRI_SYNTHSEG_DIR` to
compensate — but **not** `RSFMRI_FS_DIR`, and nothing overrode the VBG path at
all. So the FreeSurfer lookups pointed at
`<KUL_NIS>/BIDS/derivatives/freesurfer/...`, which cannot exist.

Consequences, all silent:

- Step 6 (Lausanne2018 scale3 subject-specific atlas) **never ran for any
  subject**, regardless of whether the parcellation existed. Confirmed on a real
  study whose `BIDS/derivatives/freesurfer/sub-X/mri/` held all five Lausanne
  scales: no `*_lausanne_scale3_*` output had ever been produced.
- The Lausanne-based somatotopic SBA seeding added in `636ae80` therefore could
  never find its input as shipped.
- The message it printed — "run KUL_FS_multiparc.sh first" — actively misled,
  since multiparc had run and produced the files.

Fixed by exporting `RSFMRI_FS_DIR` and a new `RSFMRI_VBG_DIR` from the caller,
and having `step0_synthseg.sh` honour the latter for its VBG lookup. The VBG
path *layout* in step0 was already correct (`output_VBG/sub-X_FS_output/sub-X`,
verified against a real KUL_VBG output tree) — only its root was wrong.

Note `KUL_clinical_fmridti.sh`'s FWT lookup uses a different VBG layout
(`output_VBG/sub-X/sub-X_FS_output/sub-X`) which does not match what KUL_VBG
writes. That one is harmless in practice — `KUL_run_VBG` copies the FS tree into
`BIDS/derivatives/freesurfer/`, which is the fallback it lands on — so it is
noted rather than changed here.

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

### Fixed — VIM FreeSurfer lookup picked the wrong tree

Found on the first real test. The candidate loop selected the first FreeSurfer
directory that *existed*, not the first that actually contained the thalamic
segmentation. A tumour case can carry both a KUL_VBG FreeSurfer output and a
plain one, and `segment_subregions` may have been run on only one — so the loop
locked onto the VBG tree and reported "no thalamic segmentation found" while the
file sat in `BIDS/derivatives/freesurfer/`. It now globs for
`ThalamicNuclei*FSvoxelSpace.mgz` in each candidate and takes the first real hit.

### New — thalamic VIM labels for ET DBS (type 2)

`KUL_karawun_prepare.sh` now exports the thalamic VIM as a Brainlab label for
type 2 (essential tremor / DRT) cases. The DRT tract already reached Brainlab;
the target it is aimed at did not.

Source is FreeSurfer's thalamic subnuclei segmentation (`segment_subregions
thalamus`, run by `KUL_FS_multiparc.sh`), labels **8129 Left-VLp** and **8229
Right-VLp** — VLp is the standard FreeSurfer analogue of the VIM target, there
being no nucleus literally named VIM in that atlas. Resampled onto the Karawun
grid with `mri_vol2vol --regheader --nearest`, written as `VIM_Left.nii.gz` /
`VIM_Right.nii.gz` at colours **16** and **30**. A side with under 10 voxels is
dropped with a warning rather than exported as a sliver.

Verified on real data (sub-VanRooyRosalia, FreeSurfer 8.2): physical volume is
preserved across the resample onto the 0.625 mm³ Karawun grid — Left-VLp
823.0 → 818.7 mm³, Right-VLp 752.0 → 757.5 mm³ (nearest-neighbour rounding,
<1%), at values 16 and 30, with the left/right asymmetry intact.

Worth contrasting with the STN VOIs beside it: those are an atlas region
(DISTAL, *symmetrised*, so left and right are mirror images by construction)
warped into the subject and repurposed from the tractography inclusion VOI. The
VIM label is segmented from the subject's own T1w, so it carries real individual
anatomy and genuine left/right asymmetry.

The filename is globbed, not hardcoded: FS 8.x writes
`ThalamicNuclei.FSvoxelSpace.mgz` while FS 7.x wrote
`ThalamicNuclei.v12.T1.FSvoxelSpace.mgz`. The docstring in
`KUL_FS_multiparc.sh` still named the 7.x file and has been corrected.

Groundwork for a later refinement: intersecting this label with the DRT
streamlines that actually pass through it, for a more specific target.

### Palette budget (this is now an explicit, verified allocation)

With the extended 64-entry fork (indices 0-63, 0 = background):

| range | use |
|---|---|
| 1-41 | known-tract table (fixed; changing these breaks scene continuity) |
| 23, 24 | DBS STN VOIs (pre-existing, in a gap the table leaves free) |
| 16, 30 | thalamic VIM left/right (type 2 / ET), likewise in a free gap |
| 42-49 | auto-assigned tracts (bundles outside the table) |
| 50 | lesion |
| 51-63 | fMRI activation labels |

Tracts get the whole low block, then one slot for the lesion, then fMRI.
`2, 6, 8, 10, 12, 14` are left free for future tract entries.

Left/right bundle colours stay merged: splitting them needs exactly the 8 free
low slots, which would leave nothing for VIM, and would change the colour of
every right-sided bundle relative to scenes already reviewed. Deferred.

**This makes the extended-palette fork a requirement for fMRI labels.** Every
value above is >30, and stock karawun clamps anything >30 to its last entry, so
without the fork every activation renders in the same colour as every other one.
That is the deliberate trade for keeping the entire low block available to
tracts.

Verified disjoint programmatically: no tract, lesion, VOI or fMRI colour can
coincide, and nothing exceeds the ceiling of 63.

Two things forced the ranges rather than a simple offset:

- **Auto-assigned tract colours were unbounded.** They started at 42 and
  incremented per unknown bundle. The shipped `tracks_list.txt` has 71 bundles,
  34 of them outside the known table, so they would have run 42..75 — through
  the fMRI range and past the palette ceiling. They now cycle within 42-55 and
  warn once when they wrap. Two unknown bundles sharing a colour is a mild
  annoyance; a bundle taking an fMRI activation's colour is a misread waiting
  to happen.
- **Total demand exceeds supply** (41 + 34 + 13 + 1 = 89 > 63), so some reuse is
  unavoidable. Confining it to *within* the auto-tract range is the point.

### Fixed (`KUL_make_fMRI_labels.sh`)

This interactive helper — orphaned, nothing calls it — writes a hand-made fMRI
label straight into `Karawun/*/labels/`. Every value it hardcoded collided:
20/21/22/25/26/27 are tract colours and **23/24 are the DBS STN VOIs**, so a
manually added `afMRI_TAAL` came out the same colour as the left STN VOI. It is
the origin of the stale `afMRI_TAAL_thres_5.5.nii.gz` (value 23) found in a real
patient's Karawun folder. Now draws from the same reserved pool, and its header
says explicitly that it is a manual one-off helper, not part of the `-R` flow.

### Fixed (fMRI label colours)

- **fMRI activation labels shared colours with tracts.** They were numbered
  `1..N`, and the tract table also starts at 1, so task 1 got colour 1 =
  Arcuate Fasciculus: on a language case the activation and the language tract
  rendered identically — the one pair you most need to tell apart. Tasks 3, 4
  and 5 landed on CST_LT, CST_RT and Cingulum.

  `fe3b800` merged left/right colours for well-separated tracts specifically to
  free budget "for fMRI labels without colliding", but the labels were never
  moved into the freed range, so indices 2, 6, 8, 10, 12, 14 and 30 sat unused
  while the labels collided. This finishes that job: fMRI labels now draw from
  an explicit reserved list, ordered by measured worst-case CIEDE2000 distance
  to every tract and lesion colour, best first. fMRI colour 2 sits dE 77.7 from
  the Arcuate, versus dE 0.00 before.

  The list stays inside 1-30 for the first seven tasks so it works on stock
  karawun, and only then spills into the fork's 50+ range. Those fork colours
  are better separated (~13) than low slots 4-7, but are listed after rather
  than before: on stock karawun anything >30 clamps to one entry, so leading
  with them would make every task past the third render identically.

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
