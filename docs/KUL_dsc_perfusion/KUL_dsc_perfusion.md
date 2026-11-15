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
DSC,^T2[*] DSC Perfusion,-,2,j
```

That is the Philips description `T2* DSC Perfusion<suffix>`: `[*]` matches the
literal `*` (a bare `*` is a wildcard), and `^` keeps out the console's
`Reg - T2* DSC Perfusion...` copy. Take the string from the SeriesDescription,
not from the dcm2niix file name in `tmp_dcm2bids/` (which reads
`T2_DSC_Perfusion...` and never matches). `study_config/_base/sequences.txt`
already carries the Philips and Siemens rows.

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
| 1 | Resolve inputs | DSC series, anatomical reference, lesion mask, TE/TR and PE axis from the BIDS sidecar |
| 2 | Prepare the anatomy | `mri_synthstrip`, then regrid to the DSC voxel size (registration target only) |
| 3 | Denoise | MP-PCA (`dwidenoise`); in-plane patches (e×e×3) when slices are thick; `-D` skips it |
| 4 | Motion correct | `antsMotionCorr`, rigid, BSpline — the only resampling the 4D series gets |
| 5 | Distortion correct | rigid + affine + SyN to the anatomy, deformation restricted to the PE axis (transform only) |
| 6 | Masks | brain mask and dilated lesion mask, both on the **native DSC grid** |
| 7 | Fit | `share/dsc/KUL_dsc_fit.py` on the native grid — see [The fit](#the-fit-method-v2) |
| 8 | Deliver | the 3D maps and the DSC coverage mask to the full-resolution T1w grid, **one** Linear interpolation (affine + SyN) |
| 9 | Vendor maps | optional (`-I`); maps not on the DSC grid are skipped |
| 10 | Reference ROIs | contralesional NAWM (normalisation) and contralesional cortex (QC) |
| 11 | ROI report | `share/dsc/KUL_dsc_report.py`: normalised maps, per-ROI measures, PDF report and the `REPORT/` montage — see [The ROI report](#the-roi-report). Always runs; without a NAWM reference it holds absolute values only |

### The fit (method v2)

| Stage | What it does | Evidence |
|---|---|---|
| Timing | Baseline = every steady-state frame (frame 0 dropped) up to 2 frames before bolus arrival; arrival = where the brain-average signal leaves its early pre-bolus level by more than max(3σ, 10 % of the bolus depth). End of first pass = recirculation trough of the reference curve. | More baseline frames, less CBV noise (Boxerman 1997 JMRI; ASFNR, Welker 2015 AJNR) |
| S0 | Mean over the whole baseline (noise), levelled to the signal just before arrival using the brain-average curve, so a slow pre-bolus drift does not become a ΔR2\* offset that the integral multiplies (`--no_drift_adjust` in the fit to disable). | phantom, below |
| Leakage | Boxerman–Schmainda–Weisskoff: per voxel ΔR2\*(t) = K1·ref(t) − K2·∫ref; corrected curve = ΔR2\* + K2·∫ref. `ref` = mean ΔR2\* of non-enhancing brain (outside the lesion, not T1/T2\*-leaky, not CSF, not vessels). The reference is shifted to each voxel's bolus arrival (arrival-time correction). K2 > 0 = T1-dominant leakage, K2 < 0 = T2\*-dominant. | BSW, Boxerman 2006 AJNR; 2020 consensus (Boxerman 2020 Neuro-Oncol); arrival-time correction, Leigh 2012 PLoS One |
| rCBV | Trapezoidal integral of the (corrected) ΔR2\* from 2 frames before arrival to the end (`-w full`, default) or to the end of the first pass + 3 frames (`-w firstpass`). Not AIF-normalised: rCBV does not depend on the AIF. | ASFNR (Welker 2015) |
| AIF | Cluster selection (Mouridsen 2006 MRM): the 10 % of voxels with the largest first-pass area, minus the 25 % roughest, k-means on area-normalised curves, the cluster with the earliest first moment (twice), then its 20 highest-peaking voxels. Averaged; no smoothing over the peak. Checked: the AIF must peak earlier and > 2× higher than tissue, otherwise a warning says rCBF/MTT/Tmax are unreliable. | Mouridsen 2006; Bjørnerud & Emblem 2010 JCBFM |
| rCBF, Tmax | Block-circulant SVD with per-voxel oscillation-index threshold (oSVD, OI 0.095), on the leakage-corrected curve — insensitive to bolus delay. Tmax = time of the residue maximum. | Wu 2003 MRM |
| MTT | CBV/CBF with CBV = ∫corrected/∫AIF over the same window (central volume theorem), seconds. | |
| TTP, TT0 | From the temporally smoothed (σ = 1 frame) corrected curve, searched only inside the first pass, sub-frame interpolated; TT0 walks back from the peak to max(3σ, 15 % of peak). 0 where no bolus is detectable (CNR < 3). | |
| CNR | Peak of the smoothed corrected curve / SD of its baseline: a per-voxel "is there a measurable bolus here" map. | |
| Smoothing | None on the primary maps. With `-s FWHM` (default 3 mm) extra `*_smooth` copies of rCBV and rCBF are written — normalised 3D Gaussian inside the mask, for reading/display only. | ASFNR: "a small kernel may be considered" |

QC written to `BIDS/derivatives/KUL_compute/sub-X/KUL_dsc_perfusion/fit/`:
`qc_timing.png` (baseline/arrival/window on the mean curve; AIF and reference
curves), `curves.tsv`, `aif_mask.nii.gz`, `reference_tissue_mask.nii.gz` and
`KUL_dsc_fit.json` (every detected frame index, drift, AIF checks, K2 of the
reference, median oSVD threshold, fraction of voxels with a detectable bolus).

### Outputs

`RESULTS/sub-{participant}/Perfusion/`:

```
sub-X_rCBV_corrected.nii.gz          sub-X_nrCBV_corrected.nii.gz          (NAWM-normalised)
sub-X_rCBV_uncorrected.nii.gz        sub-X_nrCBV_uncorrected.nii.gz
sub-X_rCBF.nii.gz                    sub-X_nrCBF.nii.gz
sub-X_rCBV_corrected_smooth.nii.gz   sub-X_nrCBV_corrected_smooth.nii.gz   (-s > 0; display only)
sub-X_rCBV_uncorrected_smooth.nii.gz sub-X_nrCBV_uncorrected_smooth.nii.gz
sub-X_rCBF_smooth.nii.gz             sub-X_nrCBF_smooth.nii.gz
sub-X_MTT.nii.gz    sub-X_Tmax.nii.gz    sub-X_TTP.nii.gz    sub-X_TT0.nii.gz
sub-X_K1.nii.gz     sub-X_K2.nii.gz      sub-X_CNR.nii.gz
sub-X_contralesional_NAWM_mask.nii.gz    sub-X_contralesional_GM_mask.nii.gz
sub-X_dsc_coverage_mask.nii.gz           where the maps are defined (the DSC slab)
sub-X_perfusion_report.tsv       ONE ROW PER ROI: the measures you read (see below)
sub-X_perfusion_roi_stats.tsv    every map x every ROI: n, mean, std, p10..p95, value and ratio to NAWM
sub-X_perfusion_stats.tsv        lesion vs NAWM mean/median/std/voxels (legacy format)
sub-X_perfusion_summary.tsv      lesion median, NAWM median, normalised ratio (legacy format)
sub-X_perfusion_reference.txt    how the reference was chosen, and the method used
ROIs/                            optional: drop extra masks / label maps here (see below)
```

`REPORT/`:

```
sub-X_06_DSC_report.pdf    p1 ROI table + acquisition/fit QC, p2 nrCBV/nrCBF/K2 montages
                           with the ROI contours, p3 the fit QC curves
sub-X_06_DSC_rCBV.png      nrCBV montage through the lesion (radiological), ROI contours
```

Units: `rCBF` in 1/s relative to the (partial-volume-scaled) AIF, `MTT`/`Tmax`
in seconds, `TTP`/`TT0` in seconds from the start of the series, `rCBV` in
arbitrary units (the integral of ΔR2\*, which is why it is normalised against
NAWM before being read), `K2` in 1/s, `K1` and `CNR` dimensionless. The
unsmoothed maps are the ones the stats, PACS and Karawun use.

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

### The ROI report

`share/dsc/KUL_dsc_report.py` produces every perfusion number, so the TSVs, the
PDF and the normalised maps always agree. It runs at the end of every
invocation; on a re-run only it (and the cheap mask steps) recompute, so adding
an ROI and re-running takes well under a minute.

**ROIs**, each used only when the file exists:

| ROI | Source | Note |
|---|---|---|
| `enhancing` | `Lesion/sub-X_hdglio_lesion_solid_tissue` | HD-GLIO contrast-enhancing label, morphologically closed — it **includes enclosed necrosis**, so read its percentiles, not just its median |
| `nonenhancing` | `Lesion/sub-X_hdglio_lesion_perilesional_tissue` | non-enhancing FLAIR abnormality |
| `cavity` | `Lesion/sub-X_resseg_cavity_only` | resection cavity |
| `lesion` | the lesion mask (`-l`, default `lesion_and_cavity`) | whole lesion; also feeds the legacy tables |
| *user* | `Perfusion/ROIs/*.nii.gz` | a binary mask is one ROI; an integer label map is one ROI per label, named from an ITK-SNAP label file next to it (`<stem>.txt`, or `<stem>_itksnap.txt` as written by `KUL_lesion_compartments.py`). Any grid — resampled nearest-neighbour |
| `NAWM` | contralesional NAWM | the normalisation reference |
| `contralesional_GM` | contralesional cortex (aseg 3/42 minus the dilated lesion) | QC only: cortex/NAWM rCBV should be clearly > 1 |

**Rules applied to every number**

* Only voxels inside the DSC slab count (`dsc_coverage_mask`): the maps are 0
  outside it, and a lesion reaching beyond the slab would otherwise be pulled
  towards 0. The coverage % is reported per ROI.
* MTT, Tmax, TTP and TT0 use only voxels with a detectable bolus (CNR ≥ 3);
  MTT/TTP/TT0 voxels the fit could not estimate (0) are excluded, not averaged in.
* Always the **unsmoothed** maps; `*_smooth` are display only.
* An ROI covering fewer than 10 native DSC voxels, under 90 % inside the slab, or
  with a bolus in under half its voxels is flagged.

**Measures in `perfusion_report.tsv`** (ratios = value / contralesional NAWM median)

| Column | Meaning | Evidence |
|---|---|---|
| `nrCBV_corr_median` | typical leakage-corrected rCBV ratio of the ROI | |
| `nrCBV_corr_p90`, `_p95` | reproducible "hotspot" — replaces a hand-placed max ROI | histogram measures more reproducible than max-ROI: Wetzel 2002 (PMID 12202717) |
| `pct_nrCBV_corr_gt_1`, `_gt_1.75` | % of the ROI above 1.0 and above 1.75 (`--thresholds` in the report script) | 1.75 high-grade cut-off: Law 2003 (PMID 14625221); threshold volume fraction cf. fractional tumour burden: Hu 2012 (PMID 22561797) |
| `nrCBV_uncorr_median` | without leakage correction, for comparison | |
| `nrCBF_median`, `_p90` | rCBF ratio | |
| `MTT_median_s`, `Tmax_median_s`, `K2_median` | absolute | |
| `volume_ml`, `dsc_coverage_pct`, `native_dsc_voxels`, `bolus_detected_pct`, `flags` | how much to trust the row | |

The NAWM row's `% > 1.75` is the background rate produced by voxel noise alone
(about 17 % in the validation subject with 60 frames and 22 baseline frames); a
lesion fraction is only meaningful well above it.

The PDF also lists the fit QC (baseline frames, AIF quality, drift, reference
K2, bolus detection) and warns when the acquisition is below the 2020 consensus
(Boxerman 2020, PMID 32516388) or the cortex/NAWM ratio is implausibly low.

---

## Method v2 (2026-11): why the fit changed

v1 maps looked noisier and more mottled than the Siemens console maps of the
same patient. A comparison on the native DSC grid, an independent audit and a
literature review traced it to:

1. **Leakage model.** v1 regressed each curve on the integral of a single-voxel
   AIF (`ΔR2* ~ K2·∫AIF + K1`) with no intravascular term, so K2 absorbed the
   normal recirculation (22–38 % of NAWM rCBV removed) and depended on which
   voxel was picked as AIF. Replaced by reference-tissue BSW.
2. **Baseline.** The 5 %-of-depth plateau rule against a noisy pre-peak maximum
   kept only 7 of ~21 pre-bolus frames; the S0 error was then multiplied by an
   81 s integral.
3. **Resampling.** The 4D series was Lanczos-resampled onto a grid with the
   anatomy's axes (tilted ~7° to the 4 mm slices + gap), then the maps were
   upsampled again — blurred, displaced, "blocky". The fit now runs on the
   native grid and the maps are interpolated once.
4. **AIF / CBF / timing.** A single noisy AIF voxel from a mostly venous
   candidate pool, a Savitzky–Golay filter as wide as the bolus, delay-sensitive
   truncated SVD, and TTP/TT0 taken over all frames (TT0 landed in pre-bolus
   noise in 96 % of voxels).

**Synthetic phantom** (known truth, 60 frames, TR 1.35 s, tSNR 40/80; ratios to
WM; v1 → v2):

| | truth | v1 | v2 |
|---|---|---|---|
| GM rCBV, bolus delayed 3 s (tSNR 80) | 2.00 | 1.52 | 1.92 |
| tumour rCBV, T1 leakage (tSNR 80) | 4.00 | 2.44 | 2.73 |
| GM rCBV with a 1 % pre-bolus drift (tSNR 40) | 2.00 | 1.56 | 1.94 |
| WM rCBV noise (robust CV, tSNR 40, drift) | — | 1.85 | 0.66 |
| TT0 (s; true arrival 29.7 s, 15 %-of-peak definition) | ~33 | 1.4 | 33–34 |

**Validation subject** (odd/even-frame split-half noise; Siemens not used as
truth): rCBV_corrected noise CV 0.44 → 0.19, GM–WM contrast-to-noise 1.24 →
2.07, lesion/NAWM stable at 1.46 across baseline, window and data variants.

## Changes from the prototype (method v1, kept for history)

These were the v1 changes relative to the original research prototype; several
were superseded by method v2 above (AIF ranking, truncated SVD, the 4D
resampling, N4). The v1 code is in git history before the method-v2 commit.

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
- **AIF candidates were ranked by absolute signal drop** (S₀ − S_min), which
  is S₀·(1 − S_min/S₀) and so scales with receive-coil sensitivity. With the
  64-channel coil this favoured bright voxels near the coil: on one series it
  picked a frontal-surface voxel at 2× the median S₀ (peak ΔR2\* 37/s, FWHM
  12 s) over a deep temporal artery (87/s, 9.5 s). Candidates are now ranked
  by relative drop (1 − S_min/S₀, i.e. what ΔR2\* sees), PCA/KMeans runs on
  S/S₀ curves, and voxels with S₀ below half the median are excluded.
  `-A absolute` (or `dsc_aif_ranking: absolute` in `study_config/run_dsc.txt`)
  restores the old ranking to reproduce earlier maps; the fit JSON records
  which one was used (`aif_ranking`).
- Hardcoded lesion paths, and a `python ./Good_DSCLC_fit5.py`
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

## Known limitations

- **BSW assumes the tissue's transit time matches the reference.** Where tumour
  MTT is much longer than normal tissue, BSW under-corrects: −14 % rCBV in the
  phantom for a non-leaky tumour with MTT 6 s vs 4 s, and T1-leaky tumour rCBV
  is still underestimated (2.7 vs true 4.0 with strong leakage). Bjørnerud 2011
  (JCBFM, PMID 21505483) addresses this via the residue function; not implemented.
- **MTT in seconds is biased high** (~1.5× in the phantom at tSNR 80): SVD
  regularisation underestimates CBF. Read rCBF/MTT as ratios to normal tissue.
- **rCBF depends on the AIF shape.** Check `qc_timing.png`; the fit warns when
  the AIF is not clearly arterial. rCBV does not use the AIF.
- **The protocol limits what any software can recover.** The 2020 consensus asks
  for ≥ 120 frames with 30–50 baseline frames, 3 mm slices without gap at 3T,
  and 30° flip angle without preload or 60° with a full-dose preload. The
  validation series had 60 frames, ~21 baseline frames and 4 mm + 1.2 mm gap.
- **Vendor maps are not ground truth.** Packages disagree on identical data
  (Kudo 2013 Radiology; Bell 2019) and need different rCBV cut-offs. The Siemens
  relCBVCorr map looks smoother than the v2 unsmoothed map because it is
  smoothed (and clips negatives); the `*_smooth` maps are the like-for-like
  view.
- **ROI numbers are only as good as the ROI.** The compartments come from the
  automatic segmentation, which the report does not second-guess: in the
  validation subject the 1.3 mL ResSeg "cavity" lies on the vascular posterior
  rim and has nrCBV 2.35, i.e. it is not fluid. Look at the contours on page 2
  of the PDF before quoting a compartment, and drop a corrected mask in
  `Perfusion/ROIs/` when needed. A DSC voxel (here 1.7 × 1.7 × 5.2 mm) is ~20×
  a T1w voxel, so small ROIs and thin rims are partial-volumed whatever their
  T1w-grid voxel count.
- **Validated on one patient plus a phantom.** A DSC digital reference object
  (Semmineh 2017 Tomography; Bell 2019) is the next check.

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
| `-w` | rCBV integration window: `full` (default), `firstpass`, or `START:END` (`dsc_window` in `run_dsc.txt`) |
| `-s` | FWHM (mm) of the extra `*_smooth` maps; default 3, `0` = none (`dsc_smooth_fwhm` in `run_dsc.txt`) |
| `-D` | skip MP-PCA denoising |
| `-P` | phase-encoding axis `i`/`j`/`k`/`none` (default: from the sidecar) |
| `-A`, `-L` | v1 options; accepted and ignored |
| `-y` | conda env instead of `$KUL_PYFMRI_ENV` |
| `-R` | redo — discard previous results and recompute |
| `-n` | threads (default 32) |
| `-v` | verbosity 0/1/2 (default 1) |

## Requirements

ANTs, MRtrix3, FreeSurfer (`mri_synthstrip`, `mri_vol2vol`), and the conda env
named by `$KUL_PYFMRI_ENV` (default `pyfMRI`, created by the installer's
`env-pyfmri` section). The fit needs numpy, scipy, nibabel, matplotlib and
scikit-learn — the last arrives transitively via nilearn.
