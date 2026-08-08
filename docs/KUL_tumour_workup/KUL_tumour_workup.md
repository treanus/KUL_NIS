# User guide — tumour work-up with perfusion, Brainlab and PACS

A practical, start-to-finish walkthrough of a glioma case: DICOM in, perfusion
maps, tract labels and PACS series out.

This guide covers the pieces added in August 2026 — DSC perfusion, the lesion as
a Brainlab label, the restored FAT1w QA volume, and the lesion/perfusion PACS
series. For the pipeline as a whole see
[KUL_clinical_fmridti](/docs/KUL_clinical_fmridti/KUL_clinical_fmridti.md); for
the perfusion internals see
[KUL_dsc_perfusion](/docs/KUL_dsc_perfusion/KUL_dsc_perfusion.md).

---

## 0. One-time setup

**Conda environments.** The perfusion fit runs in `$KUL_PYFMRI_ENV` (default
`pyfMRI`), created by the installer's `env-pyfmri` section. It needs numpy,
scipy, nibabel, matplotlib and scikit-learn; the last arrives via nilearn.
The script exits early with a clear message if the env is missing.

**Karawun.** Colour indices above 30 are clamped by stock karawun, so bundles
assigned 31–41 all render in the same colour in Brainlab — and so do all the
fMRI activation labels, which live at 51–63. Install the KU Leuven
extended-palette fork into `KarawunDev` and pin it. Indices 0–30 are unchanged
there, so existing scenes are unaffected. See `KUL_PALETTE_NOTES.md` in that
fork, and
[KUL_karawun_prepare](/docs/KUL_karawun_prepare/KUL_karawun_prepare.md) for the
full label-colour convention.

**Tell dcm2bids about your DSC sequence.** Add a line to
`study_config/sequences.txt` with enough of the series description to be
unambiguous:

```
# Identifier,search-string,task,mb,pe_dir,acq_label
DSC,T2_DSC_Perfusion,-,2,j
```

The `pe_dir` column is not cosmetic — it is written into the BIDS sidecar as
`PhaseEncodingDirection`, and the perfusion script restricts its distortion
correction to that axis. Leave it as `-` and it assumes `j`.

---

## 1. Convert

```bash
KUL_clinical_fmridti.sh -p JaneDoe -d DICOM/JaneDoe.zip -n 32
```

The DSC series lands in the `perf/` datatype directory, a sibling of
`anat/`, `func/` and `dwi/`:

```
BIDS/sub-JaneDoe/
├── anat/   ├── dwi/   ├── func/   └── perf/
                                       ├── sub-JaneDoe_dsc.nii.gz
                                       └── sub-JaneDoe_dsc.json
```

`perf/*dsc*` is added to `.bidsignore` — BIDS has no standardised DSC suffix
yet, so the validator would otherwise reject it.

**Check the sidecar before going further.** `EchoTime` and `RepetitionTime` must
be present and in **seconds**:

```bash
python3 -m json.tool BIDS/sub-JaneDoe/perf/sub-JaneDoe_dsc.json | grep -E "EchoTime|Repetition|PhaseEnc"
```

A TE in milliseconds is the classic sidecar error and would scale ΔR2\* by 1000.
The script rejects anything above 1.0 outright rather than producing quietly
wrong maps, but it is cheaper to notice here.

---

## 2. Run

The same command runs everything. Perfusion is picked up automatically whenever
a DSC series exists, the way fMRI and dMRI are; `-W` skips it.

```bash
KUL_clinical_fmridti.sh -p JaneDoe -d DICOM/JaneDoe.zip -n 32
```

Perfusion is scheduled **after** VBG/multiparc, because the contralesional NAWM
reference needs a FreeSurfer `aseg`. Expect roughly 20–40 minutes for the
perfusion step on clinical data.

Run it standalone if you only want to redo perfusion:

```bash
KUL_dsc_perfusion.sh -p JaneDoe -n 32          # reuses existing work
KUL_dsc_perfusion.sh -p JaneDoe -n 32 -R       # recompute from scratch
```

---

## 3. Review the figures, then export

DICOM export is deliberately a second, explicit pass. Look at the PNGs first,
then re-run with `-R <underlay>`:

```bash
KUL_clinical_fmridti.sh -p JaneDoe -R 1 -n 32     # 1 = cT1w underlay
```

`-R` produces the PACS DICOMs **and** prepares the Karawun folder. Put a donor
DICOM (a single slice from a high-resolution anatomical series) in
`Karawun/sub-JaneDoe/DICOM/` first — the same donor is used for both, so the
metadata stays consistent.

---

## 4. Where everything lands

```
RESULTS/sub-JaneDoe/
├── Anat/
│   ├── T1w.nii.gz                       reference space for everything below
│   ├── cT1w_reg2_T1w.nii.gz
│   └── sub-JaneDoe_lesion.nii.gz        binarised lesion, next to the underlays
├── Lesion/
│   └── sub-JaneDoe_lesion_and_cavity.nii.gz
├── Perfusion/
│   ├── sub-JaneDoe_nrCBV_corrected.nii.gz     NAWM-normalised — the map you read
│   ├── sub-JaneDoe_nrCBF.nii.gz
│   ├── sub-JaneDoe_rCBV_corrected.nii.gz      un-normalised
│   ├── sub-JaneDoe_MTT.nii.gz  _TTP  _TT0  _K1  _K2
│   ├── sub-JaneDoe_contralesional_NAWM_mask.nii.gz
│   ├── sub-JaneDoe_perfusion_summary.tsv      lesion vs NAWM, one row per map
│   ├── sub-JaneDoe_perfusion_stats.tsv        mean/median/std/voxels per ROI
│   └── sub-JaneDoe_perfusion_reference.txt    how the reference was chosen
├── PACS/
│   ├── Clinical_cT1w/    lesion + perfusion series
│   ├── fMRI_cT1w/        activation series
│   └── Tracto_cT1w/      tract series
└── Clinical_figures_cT1w/   PNGs behind the Clinical PACS series

Karawun/sub-JaneDoe/
├── T1w.nii.gz
├── FAT1w.nii.gz          FA-weighted T1w — the registration QA volume
├── labels/
│   ├── Lesion.nii.gz     the tumour, palette colour 50
│   ├── afMRI_*.nii.gz    one per task, colours 51-63
│   └── *_center.nii.gz   one per bundle
└── tck/
```

---

## 5. Reading the perfusion output

Start with `sub-JaneDoe_perfusion_summary.tsv`:

```
map               lesion_median   nawm_median   normalised_ratio
rCBV_corrected    34.51           18.37         1.88
rCBF              0.229           0.124         1.85
MTT               2.65            2.62          1.01
```

The **normalised ratio** is the number to report — lesion median over
contralesional NAWM median. `nrCBV_corrected` is the map form of the same thing.

Units, with the default deconvolution: `rCBF` in 1/s, `MTT`/`TTP`/`TT0` in
seconds, `rCBV` in arbitrary units (the integral of ΔR2\*, which is exactly why
it is read as a ratio).

`K2` normally shows `n/a` in the ratio column. That is correct, not a failure:
it is a leakage coefficient, not a perfusion quantity, and its NAWM median is
often near zero or negative. Normalising by it would turn noise into a large
number, so the script refuses.

**The PACS perfusion series use fixed thresholds**, not the auto `max/3` applied
to activation maps: **1.75** on nrCBV (the conventional high-grade glioma
cutoff) and **1.0** on nrCBF ("above contralesional normal white matter"). These
are normalised ratios, so a threshold carries a fixed clinical meaning that
`max/3` would discard. `-T` overrides.

**Always open `sub-JaneDoe_perfusion_reference.txt`.** It records which
hemisphere the reference came from and how the laterality was decided. A
bilateral or midline lesion is where this is most likely to pick something you
would not have.

---

## 6. Brainlab import

The `importTractography` command is printed at the end of the Karawun prep. Load
**both** anatomicals:

```bash
conda activate KarawunDev
importTractography -d Karawun/sub-JaneDoe/DICOM/*.dcm \
  -o Karawun/sub-JaneDoe/sub-JaneDoe_for_elements \
  -n Karawun/sub-JaneDoe/T1w.nii.gz Karawun/sub-JaneDoe/FAT1w.nii.gz \
  -t Karawun/sub-JaneDoe/tck/*.tck \
  -l Karawun/sub-JaneDoe/labels/*.gz
```

**Check FAT1w first.** It is √FA · T1w, so white matter tracts light up inside
the T1w anatomy. If the FA-to-T1w registration has slipped, every tract in the
scene is displaced the same way — and nothing else in the export would show it,
because the tracts are self-consistent with each other. This is the cheapest
high-value QA step in the whole pipeline, which is why it is worth the extra
volume.

If `FAT1w.nii.gz` is absent, `KUL_dwiprep_anat.sh` has not produced
`dwiprep/sub-JaneDoe/sub-JaneDoe/qa/fa_reg2T1w.nii.gz` yet.

---

## Troubleshooting

**"no DSC series found in BIDS/…"** — the `DSC` line in `sequences.txt` did not
match. Check the actual series description:

```bash
grep -ri "dsc\|perfusion" BIDS/sub-JaneDoe/perf/*.json
```
Matching is on `SeriesDescription` only, so make the search string specific
enough not to also catch derived console maps.

**"could not read EchoTime"** — pass it explicitly: `-e 0.030`.

**"TE = 30 looks like milliseconds"** — the sidecar is in ms. Pass `-e 0.030`,
and fix the conversion so it does not recur.

**"The AIF dR2\* curve has a non-positive area"** — the baseline window overlaps
the bolus, or the series is truncated. Look at `aif_raw.png` and
`aif_deltaR2star.png` in the fit directory, then pin the window: `-b 5:15`.
This is a guard, not a crash — it fires precisely so you do not get maps with
the sign silently inverted.

**"no lesion mask found" / "no FreeSurfer aseg.mgz found"** — the parametric
maps are still written and are usable; only the NAWM ratios are skipped. Re-run
`KUL_dsc_perfusion.sh -p JaneDoe` once segmentation exists and the ratios are
added without recomputing the fit.

**"the lesion mask does not overlap either cerebral hemisphere"** — the mask is
in the wrong space, or entirely in cerebellum/brainstem. Overlay it on
`RESULTS/sub-JaneDoe/Anat/T1w.nii.gz` and check.

**A dozen bundles are the same colour in Brainlab, or all the fMRI activations
are** — stock karawun clamps colour indices above 30. Use the extended-palette
fork; see
[KUL_karawun_prepare](/docs/KUL_karawun_prepare/KUL_karawun_prepare.md).

**`mri_synthstrip: no usable GPU, falling back to CPU`** — informational. It
happens on CPU-only nodes and on GPU nodes whose FreeSurfer ships a CPU-only
torch. Slower, same result.

---

## Known limitations

**The leakage correction is not textbook Boxerman-Schmainda-Weisskoff.** BSW
regresses each tissue curve on a whole-brain non-enhancing *reference tissue*
curve and its integral; this implementation regresses on the integral of the
*arterial input function*, following the original in-house prototype. K1/K2 are
therefore not directly comparable to published BSW maps. This was kept
deliberately rather than changed silently — see "Known deviations" in the
perfusion doc.

**rCBV is not normalised by the AIF area** — it is the raw ΔR2\* integral, which
is standard for "relative" CBV and harmless when read as a NAWM ratio. MTT does
divide by the AIF area internally, which is what puts it in seconds.

**Absolute perfusion values are not comparable across scanners or sessions.**
Only the NAWM-normalised ratios are.

**One DSC series per participant.** If several are present the first is used and
a warning is printed; select explicitly with `-d`.
