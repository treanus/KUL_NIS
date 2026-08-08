# KUL_dsc_perfusion.sh

DSC (dynamic susceptibility contrast) perfusion processing for tumour patients:
preprocessing of the 4D series, leakage-corrected rCBV/rCBF/MTT estimation, and
normalisation of those maps against contralesional normal-appearing white
matter.

Every output is delivered in the participant's T1w space — the same space as
everything else under `RESULTS/sub-{participant}/` — so the perfusion maps can
go to PACS and Karawun alongside the fMRI activations and the tractography.

Origin: the workflow comes from `tools/KUL_DSC_analysis/` (A. Radwan), which was
a standalone research pipeline. See [Changes from the prototype](#changes-from-the-prototype)
for what moved and why.

---

## Quick start

Inside a KUL_NIS study directory, after `KUL_clinical_fmridti.sh` has run:

```bash
KUL_dsc_perfusion.sh -p JaneDoe -n 32
```

Standalone, outside a study directory:

```bash
KUL_dsc_perfusion.sh -p JaneDoe \
    -d /data/raw/dsc.nii.gz \
    -a /data/raw/T1w.nii.gz \
    -e 0.030 -r 1.5 -n 32
```

Within the clinical pipeline it runs by itself — no flag needed — whenever a DSC
series is present, exactly the way fMRI and dMRI data are picked up:

```bash
KUL_clinical_fmridti.sh -p JaneDoe -d DICOM/JaneDoe.zip -n 32
KUL_clinical_fmridti.sh -p JaneDoe -d DICOM/JaneDoe.zip -n 32 -W   # skip DSC
```

---

## Getting DSC into BIDS

Add a `DSC` line to `study_config/sequences.txt`, giving enough of the series
description to be unambiguous:

```
# Identifier,search-string,task,mb,pe_dir,acq_label
DSC,T2_DSC_Perfusion,-,2,j
```

`KUL_dcm2bids.sh` then writes the series to the `perf/` datatype directory, a
sibling of `anat/`, `func/` and `dwi/`:

```
BIDS/sub-JaneDoe/
├── anat/
├── dwi/
├── func/
└── perf/
    ├── sub-JaneDoe_dsc.nii.gz
    └── sub-JaneDoe_dsc.json
```

Two notes on that config line:

- The **`pe_dir` column matters.** It is written into the sidecar as
  `PhaseEncodingDirection`, and the distortion correction restricts its
  deformation to that axis. Leave it as `-` and the script assumes `j`.
- The DSC series is matched on `SeriesDescription` alone, unlike the ASL block
  which also pins `ImageType`. DSC `ImageType` strings vary a lot between
  vendors and software levels, and an over-tight criterion fails by converting
  nothing at all rather than by converting the wrong thing.

Like `perf/*asl*`, `perf/*dsc*` is added to `.bidsignore`: BIDS has no
standardised DSC suffix yet, so the validator would otherwise reject it.

---

## What it does

| # | Step | Notes |
|---|------|-------|
| 1 | Resolve inputs | DSC series, anatomical reference, TE/TR and PE axis from the BIDS sidecar |
| 2 | Prepare the anatomy | `mri_synthstrip`, then regrid to the DSC voxel size |
| 3 | Denoise | `DenoiseImage -d 4`, with the time axis mirror-padded (see below) |
| 4 | Motion correct | `antsMotionCorr`, affine |
| 5 | Distortion correct | rigid + affine + SyN to the anatomy, deformation restricted to the PE axis |
| 6 | Resample + bias correct | 4D resample onto the anatomical grid, then `N4BiasFieldCorrection -d 4` |
| 7 | Brain mask | `mri_synthstrip` on the temporal median, intersected with the anatomical mask |
| 8 | Fit | `share/dsc/KUL_dsc_fit.py` — ΔR2\*, leakage correction, SVD deconvolution |
| 9 | Deliver | resample the 3D maps to the full-resolution T1w grid |
| 10 | Vendor maps | optional (`-I`), co-registered for side-by-side QC |
| 11 | Normalise | contralesional NAWM reference, ratios and normalised maps |
| 12 | QC figure | normalised rCBV over the anatomy, into `REPORT/` |

The fit runs on the anatomy-regridded DSC grid rather than the full-resolution
T1w grid. Resampling a 4D series to 1 mm would multiply its size by roughly an
order of magnitude without adding information; only the 3D output maps go to
full resolution.

### Outputs

`RESULTS/sub-{participant}/Perfusion/`:

```
sub-X_rCBV_corrected.nii.gz     sub-X_nrCBV_corrected.nii.gz     (NAWM-normalised)
sub-X_rCBV_uncorrected.nii.gz   sub-X_nrCBV_uncorrected.nii.gz
sub-X_rCBF.nii.gz               sub-X_nrCBF.nii.gz
sub-X_MTT.nii.gz    sub-X_TTP.nii.gz    sub-X_TT0.nii.gz
sub-X_K1.nii.gz     sub-X_K2.nii.gz
sub-X_contralesional_NAWM_mask.nii.gz
sub-X_perfusion_stats.tsv        per-map, per-ROI mean/median/std/voxel count
sub-X_perfusion_summary.tsv      lesion median, NAWM median, normalised ratio
sub-X_perfusion_reference.txt    how the reference was chosen, and with what TE/TR
```

Units, with the default deconvolution: `rCBF` in 1/s, `MTT`/`TTP`/`TT0` in
seconds, `rCBV` in arbitrary units (the integral of ΔR2\*, which is why it is
normalised against NAWM before being read).

### The NAWM reference

Normalisation needs a lesion mask and a FreeSurfer `aseg.mgz`, both of which the
clinical pipeline already produces (`KUL_anat_segment_tumor.sh` and
KUL_VBG/`KUL_FS_multiparc.sh` respectively). The reference is built by:

1. resampling `aseg.mgz` onto the anatomical grid (`mri_vol2vol --regheader`);
2. splitting cerebral WM and cortex by hemisphere (aseg labels 2/41 and 3/42);
3. counting lesion voxels per hemisphere to establish laterality;
4. eroding the contralesional WM by 4 passes, and excluding everything within 4
   dilations of the lesion.

**If either input is missing the parametric maps are still written** and only
the normalisation is skipped, with a note saying so. Re-running the script once
the segmentation exists adds the ratios without recomputing the fit.

Guards worth knowing about: a lesion that overlaps neither hemisphere, or a NAWM
mask under 100 voxels, skips normalisation rather than reporting an unstable
ratio; and a map whose NAWM median is not meaningfully positive gets `n/a` in
the summary instead of a normalised volume. `K2` routinely hits that last case,
which is expected — it is a leakage coefficient, not a perfusion quantity.

---

## Changes from the prototype

These are behavioural changes, not just refactoring. `-L` restores the
prototype's deconvolution if you need to reproduce old numbers.

**Bugs that stopped it running, or ran it wrong**

- `scipy.integrate.cumtrapz` / `simps` were removed in scipy ≥ 1.14. The fit
  could not import at all on a current environment.
- **`DenoiseImage -d 4` corrupted the first and last frames.** It treats time as
  a fourth spatial axis, so its patch and search neighbourhoods are truncated at
  the series boundaries and those frames came back systematically darkened —
  measured at ~13 % on frame 0 of a 50-frame series, as deep as the bolus
  itself, which is a spurious ΔR2\* of several 1/s in *every* voxel. The time
  axis is now mirror-padded before denoising and cropped afterwards.
- **The baseline window was fixed at frames 5–10.** If the bolus arrives before
  frame 10, S₀ is measured part-way down the bolus, ΔR2\* goes negative, and
  every map — MTT in particular — silently changes sign. The window is now
  detected from the mask-average time course, anchored on the bolus peak; the
  fit refuses to continue if the AIF area comes out non-positive.
- **TE was scraped from the `mrinfo` comments field** and had `"0.0"` glued in
  front of it. It now comes from the sidecar `EchoTime`, and a TE that looks
  like milliseconds is rejected outright.
- The AIF candidate filter required an absolute signal drop > 100 intensity
  units — scanner- and scaling-dependent. Where more than `num_candidates`
  voxels cleared it the ranking was identical anyway; where none did, it left an
  empty pool and crashed PCA. Removed.
- Hardcoded `/mnt/DATA1/aradwa0/...` lesion paths, and a `python ./Good_DSCLC_fit5.py`
  call that only worked from one directory.

**Method changes**

- **Deconvolution.** The prototype z-scored the AIF, dropped `dt` from the
  convolution operator and regularised with `S_inv = 1/(S + 1e-3)`. rCBF was
  therefore in arbitrary units and `MTT = rCBV/rCBF` was not in seconds. The
  default is now a `dt`-scaled AIF with truncated SVD (cutoff 0.2·S_max, the
  standard sSVD choice), and MTT is normalised by the AIF area so it comes out
  in seconds. On a phantom where only bolus amplitude varies, this recovers MTT
  as amplitude-invariant (ratio 1.01) while rCBF tracks the amplitude — the
  legacy path gave an rCBF lesion/normal ratio of 1.04, i.e. no contrast at all.
- **Distortion correction.** `antsIntermodalityIntrasubject.sh` required a
  template and a subject-to-template warp purely to emit template-space extras
  nothing downstream consumed, which made a full SyN registration to MNI a
  prerequisite of every run. Replaced by the equivalent rigid + affine + SyN
  directly to the anatomy, with the deformation restricted to the phase-encoding
  axis — which a susceptibility correction should have, and which the original
  did not.
- **Performance.** The per-voxel `linregress` leakage fit and the per-voxel
  deconvolution loop are vectorised. The leakage slope has a closed form because
  the regressor is shared across voxels, so the whole volume is one
  matrix-vector product.
- **Interpolation.** The 4D series is resampled with `LanczosWindowedSinc`,
  which preserves spatial resolution; its ringing is tolerable there because one
  transform is applied to every timepoint and the artifact largely cancels in
  the S(t)/S₀ ratio. Vendor maps use `Linear` instead — they are already-derived
  3D quantities with no ratio to cancel the ringing, which would otherwise
  produce negative rCBV at the lesion rim.
- `mri_synthstrip -g` is now conditional. The prototype always passed it, which
  aborts on CPU-only nodes *and* on GPU nodes whose FreeSurfer ships a CPU-only
  torch.

---

## Known deviations

**The leakage correction is not textbook Boxerman-Schmainda-Weisskoff.** BSW
regresses each tissue curve on a whole-brain *non-enhancing reference tissue*
curve and its integral. This implementation, following the prototype, regresses
on the integral of the *arterial input function*:

```
ΔR2*_tissue(t) ~ K2·∫AIF + K1
rCBV_corr      = ∫[ ΔR2*_tissue(t) − K2·∫AIF ]
```

This was left as-is deliberately rather than changed silently — it is the model
the prototype was validated against. It does mean the K1/K2 values are not
directly comparable to published BSW K2 maps, and that the correction's
behaviour differs where the AIF shape and the mean tissue curve shape diverge.
Worth revisiting before the numbers are used for anything beyond internal
comparison.

**rCBV is not normalised by the AIF area.** It is the raw integral of ΔR2\*, in
arbitrary units, which is standard for "relative" CBV and harmless because it is
read as a NAWM ratio. MTT does divide by the AIF area internally, which is what
puts it in seconds.

---

## Options

| Flag | Meaning |
|------|---------|
| `-p` | participant name (**required**) |
| `-d` | 4D DSC series (default: `BIDS/sub-X/perf/*_dsc.nii.gz`) |
| `-a` | anatomical reference and output space (default: `cT1w_reg2_T1w`, then `T1w`) |
| `-l` | lesion mask (default: `*_lesion_and_cavity.nii.gz`, then `lesion.nii.gz`) |
| `-F` | FreeSurfer subject dir with `mri/aseg.mgz` (default: VBG output, then `BIDS/derivatives/freesurfer/`) |
| `-I` | directory of vendor/console perfusion maps to co-register |
| `-e` | echo time in seconds (default: sidecar `EchoTime`) |
| `-r` | repetition time in seconds (default: sidecar `RepetitionTime`) |
| `-b` | baseline frames `START:END`, or `auto` (default) |
| `-P` | phase-encoding axis `i`/`j`/`k`/`none` (default: from the sidecar) |
| `-L` | legacy deconvolution (arbitrary-unit rCBF, MTT not in seconds) |
| `-y` | conda env instead of `$KUL_PYFMRI_ENV` |
| `-R` | redo — discard previous results and recompute |
| `-n` | threads (default 32) |
| `-v` | verbosity 0/1/2 (default 1) |

## Requirements

ANTs, MRtrix3, FreeSurfer (`mri_synthstrip`, `mri_vol2vol`), and the conda env
named by `$KUL_PYFMRI_ENV` (default `pyfMRI`, created by the installer's
`env-pyfmri` section). The fit needs numpy, scipy, nibabel, matplotlib and
scikit-learn — the last arrives transitively via nilearn.
