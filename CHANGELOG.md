# Changelog

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
