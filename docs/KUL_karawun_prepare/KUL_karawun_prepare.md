# KUL_karawun_prepare.sh

Prepares a participant's tracts, target VOIs and lesion for import into Brainlab
Elements, via [Karawun](https://github.com/DevelopmentalImagingMCRI/karawun).

```bash
KUL_karawun_prepare.sh -p JaneDoe -t 1 -r 3
```

It is normally invoked for you by `KUL_clinical_fmridti.sh -R`, which then adds
the fMRI activation labels on top. Run it directly when you want to rebuild the
Karawun folder without redoing the DICOM export.

| flag | |
|---|---|
| `-p` | participant name (**required**) |
| `-t` | 1 = tumour (default), 2 = ET DBS (DRT), 3 = Parkinson DBS (CSHDP) |
| `-r` | relative threshold, as a percentage of tract density |
| `-a` | use the ACT output |
| `-v` | show command output |

## What it writes

```
Karawun/sub-JaneDoe/
├── T1w.nii.gz          rescaled anatomical
├── FAT1w.nii.gz        sqrt(FA) * T1w — registration QA volume
├── tck/*.tck           one per bundle
├── labels/
│   ├── <Bundle>_center.nii.gz    one per bundle
│   ├── Lesion.nii.gz             type 1, if a lesion mask exists
│   ├── VIM_Left.nii.gz           type 2, from FreeSurfer thalamic nuclei
│   ├── VIM_Right.nii.gz
│   └── DISTAL_STN_MOTOR_*.nii.gz type 3
└── DICOM/              you put a donor DICOM here
```

**Load `FAT1w.nii.gz` in Brainlab alongside the T1w.** It is √FA · T1w, so white
matter lights up inside the anatomy. If the FA-to-T1w registration has slipped,
every tract in the scene is displaced the same way — and because the tracts stay
self-consistent with each other, nothing else in the export reveals it. This is
the cheapest high-value QA step in the whole pipeline.

---

## Acquisition requirement: one plane must have isotropic voxels

**This constrains your scan protocol, and it is checked nowhere upstream.**

`check_isotropy()` in karawun refuses any volume whose three voxel dimensions
all differ:

```python
spu = np.unique(np.around(spacing, 6))
if spu.shape[0] == 3:
    raise ValueError("No plane with isotropic voxels - stopping - ...")
```

Karawun picks a slice plane to write DICOMs into and needs the **in-plane**
voxels square. Two matching dimensions is enough; three distinct is fatal.

| voxel size | |
|---|---|
| 1 × 1 × 1 | fine |
| 0.9 × 0.833 × 0.833 | fine (in-plane isotropic) |
| 1.8 × 1.8 × 4 | fine (typical clinical anisotropic slice) |
| **0.9 × 0.86 × 4.2** | **rejected** |

Two things make this worth knowing in advance:

- **It is an acquisition constraint.** If the T1w was scanned with three
  different voxel dimensions, no downstream step fixes it — you would have to
  resample, which is lossy and changes the space every other result lives in.
- **Fixing the T1w fixes everything.** Every label and anatomical in the export
  is resampled onto the T1w grid, so they inherit its spacing. Satisfy it once,
  for the T1w, and the whole export is safe.

Stock karawun only raises this at `importTractography` — the very last step,
after tractography, VBG, fMRI and the export have all run.
`KUL_karawun_prepare.sh` therefore pre-checks the T1w and warns early, using the
same rounding karawun does so the two cannot disagree.

If you are stuck with an already-acquired anisotropic T1w:

```bash
mrgrid <T1w> regrid -voxel 1,1,1 <T1w_iso>
```

and re-run the pipeline against the resampled volume.

Note the comparison rounds to 6 decimals, so two nominally-equal dimensions that
differ only by floating-point drift still count as equal. Volumes produced by
`antsApplyTransforms -r <T1w>` copy the reference geometry exactly and are safe
by construction.

## Label colour convention

**This is the part to read before changing anything.**

A label's **voxel value is its colour**. Karawun's `lookup_cie()` takes that
value, indexes a palette with it, and writes the result into the DICOM as
`RecommendedDisplayCIELabValue`. There is no separate colour field. So two
labels that share a value are *indistinguishable* in the Brainlab scene, no
matter what they are called.

Values are therefore allocated in fixed, non-overlapping ranges:

| range | use | written by |
|---|---|---|
| **1–41** | known tracts, from `KUL_karawun_tract_meta` | this script |
| **16, 30** | thalamic VIM left / right (type 2) | this script |
| **23, 24** | DBS STN VOIs (type 3) | this script |
| **42–49** | tracts *not* in the table, auto-assigned | this script |
| **50** | lesion (type 1) | this script |
| **51–63** | fMRI activation labels | `KUL_clinical_fmridti.sh -R` |

`2, 6, 8, 10, 12, 14` are deliberately left free for future tract entries.

16/30 and 23/24 sit in gaps the tract table happens to leave unused, which is
why the target VOIs can live inside the low block without disturbing it.

### The 31-entry ceiling

**Stock karawun ships a 31-colour palette (indices 0–30, 0 = background).**
`lookup_cie()` clamps anything larger to the last entry and prints
`Error - too many labels` — into a log nobody reads during a clinical run. So on
stock karawun *everything from 31 upward renders in a single colour*, including
every fMRI activation label and 13 of the tract-table entries.

The KU Leuven fork extends the palette to 64 entries (indices 0–63), with
**1–30 byte-identical to upstream** so previously exported scenes are unaffected.
The appended colours were chosen by greedy farthest-point selection under
CIEDE2000 and are each ≥ ΔE 11.6 from every other entry. See
`KUL_PALETTE_NOTES.md` in that fork.

**Pin the fork for the `KarawunDev` environment, or the table above is fiction
above index 30.** Note `KarawunEnv` is a separate, non-editable conda-forge
install that will not pick it up.

### Left/right share a colour, on purpose

For tracts with good hemispheric separation (AF, FAT, IFOF, ILF, UF,
OR_occlobe, Cingulum, MdLF) left and right share one colour. Position already
shows laterality in the 3D view, so colour is spent on encoding tract *type*
instead — which is what actually needs disambiguating.

CST, ML and PyT_SMA keep separate L/R colours: they run close to the midline in
the brainstem, where the two sides sit near each other and position alone does
not tell them apart.

Splitting the merged pairs would need exactly the 8 free low indices, leaving
none for the VIM labels, and would change the colour of every right-sided bundle
relative to scenes already reviewed. Deferred deliberately.

### If you need to change a colour

- **Prefer a free index** (`2, 6, 8, 10, 12, 14`) over reassigning one in use.
  Changing a value in `KUL_karawun_tract_meta` changes what a surgeon sees for a
  bundle they may already know by colour.
- **Stay inside the range** for what you are adding. A tract taking an fMRI
  value, or vice versa, is a misread waiting to happen.
- **Check separation, not just uniqueness.** The stock palette is not uniformly
  spaced — its worst pair (indices 3 and 26) is ΔE 2.80, below the
  just-noticeable threshold, so those two are effectively the same colour despite
  being different indices. `tools/extend_palette.py` in the karawun fork
  computes these distances.
- The ranges are checked programmatically; if you change them, re-verify that
  tracts, VOIs, lesion and fMRI remain disjoint and nothing exceeds 63.

### Auto-assigned tracts wrap

Bundles not in `KUL_karawun_tract_meta` get colours from 42–49. The shipped
`tracks_list.txt` has 71 bundles, 34 of them outside the table, and total demand
(41 + 34 + 13 + 1 = 89) exceeds the 63 usable indices — so reuse is unavoidable.
It is confined *within* the auto range, and warns once when it wraps:

```
WARNING: more bundles outside the known-tract table than reserved auto colours
(42-49); colours will now repeat between unknown bundles.
```

Two unknown bundles sharing a colour is a mild annoyance. A bundle taking an fMRI
activation's colour is not. Add frequently used bundles to
`KUL_karawun_tract_meta` to give them a stable name and colour.

---

## Target VOIs: what they actually are

The two DBS target labels are **not** equivalent in provenance, and it matters
when reading them.

**STN (type 3), colours 23/24** — the DISTAL atlas motor STN
(`DISTAL_STN_motor_bilateral_in_FSL_6thgen_symm.nii.gz`, shipped in
`KUL_FWT_templates/`), warped into the subject through two concatenated
nonlinear transforms, and repurposed from the CSHDP tract's first inclusion VOI.
It is **symmetrised**: left and right are 230 voxels each, mirror images by
construction. Its accuracy is entirely the accuracy of the registration, and no
subject-specific STN segmentation is involved anywhere.

**VIM (type 2), colours 16/30** — labels 8129 `Left-VLp` and 8229 `Right-VLp`
from FreeSurfer's thalamic subnuclei segmentation (`segment_subregions
thalamus`, run by `KUL_FS_multiparc.sh`), resampled onto the Karawun grid.
Segmented from the subject's own T1w, so it carries real individual anatomy and
genuine left/right asymmetry. VLp is the standard FreeSurfer analogue of the VIM
target — that atlas has no nucleus literally named VIM.

A side with fewer than 10 voxels is dropped with a warning rather than exported
as a sliver.

The filename is globbed, not hardcoded: FS 8.x writes
`ThalamicNuclei.FSvoxelSpace.mgz`, FS 7.x wrote
`ThalamicNuclei.v12.T1.FSvoxelSpace.mgz`.

---

## Importing into Brainlab

The `importTractography` command is printed at the end of the run. Put a donor
DICOM (one slice from a high-resolution anatomical series) in
`Karawun/sub-{participant}/DICOM/` first — the same donor should be used for the
PACS export, so metadata stays consistent.

```bash
conda activate KarawunDev
importTractography -d Karawun/sub-JaneDoe/DICOM/*.dcm \
  -o Karawun/sub-JaneDoe/sub-JaneDoe_for_elements \
  -n Karawun/sub-JaneDoe/T1w.nii.gz Karawun/sub-JaneDoe/FAT1w.nii.gz \
  -t Karawun/sub-JaneDoe/tck/*.tck \
  -l Karawun/sub-JaneDoe/labels/*.gz
```

Things to check in the output:

- **`Error - too many labels`** — the palette fork is not active; everything
  above index 30 is rendering in one colour.
- **Anatomical intensities.** Volumes handed to karawun must stay within 16-bit:
  it writes `LargestImagePixelValue` (0028,0107) from the *original* intensities,
  not from the values it rescales the pixel data to, and that tag is US (max
  65535). Both `T1w.nii.gz` and `FAT1w.nii.gz` are rescaled here for that reason.
- **Whether Elements honours `RecommendedDisplayCIELabValue` at all.** Some
  viewers let the user recolour on import and ignore the hint — worth confirming
  with two or three labels before relying on the whole scheme.

## Requirements

MRtrix3, Karawun (the KU Leuven extended-palette fork for `KarawunDev`), and
FreeSurfer (`mri_vol2vol`) for the VIM labels. `KUL_FAT1w.py` needs
`dwiprep/sub-X/sub-X/qa/fa_reg2T1w.nii.gz`, i.e. `KUL_dwiprep_anat.sh` must have
run; without it the QA volume is skipped with a message.
