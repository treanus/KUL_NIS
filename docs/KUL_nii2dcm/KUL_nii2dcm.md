# KUL_nii2dcm.py

Converts a NIfTI volume, 3D-TIFF, or a directory of PNG slices (e.g. `mrview` screenshots) to a DICOM series, given a donor DICOM for patient/study metadata. Used by `KUL_clinical_fmridti.sh -R` to wrap review figures into DICOMs for PACS and Karawun/Brainlab.

## Usage

```
KUL_nii2dcm.py [options] <nifti|tiff|png_dir> <donor.dcm> <dicomdir>
```

Plain NIfTI/3D-TIFF (geometry taken from the input volume itself, no underlay needed):

```
KUL_nii2dcm.py -s "tractography_result" -n 900 result.nii.gz donor.dcm DICOM_out/
```

PNG directory, TRA or COR screenshots (position/orientation computed from the underlay NIfTI; `-p` gives pixel spacing when no underlay is supplied):

```
KUL_nii2dcm.py -s "FT_tract_TRA" -p 0.5 -u underlay_T1w.nii.gz -o TRA png_dir/ donor.dcm DICOM_out/
```

PNG directory, SAG screenshots matched to the donor's own frames (`-M`, no `-p` needed — spacing is copied from the donor):

```
KUL_nii2dcm.py -s "FT_tract_SAG" -u underlay_T1w.nii.gz -o SAG -M png_dir/ donor.dcm DICOM_out/
```

These are the exact patterns `KUL_clinical_fmridti.sh` uses when converting `mrview` screenshots to DICOM: SAG orientations always go through `-M` (donor-match mode), TRA/COR always go through `-u`/`-o` with an explicit `-p`.

### Using it outside the pipeline

There are two independent routes, depending on how much of the workflow you want.

**This script on its own.** `KUL_nii2dcm.py` has no dependency on the pipeline at all — give it a volume, a donor and an output directory. This is the whole tool for a one-off conversion:

```bash
KUL_nii2dcm.py -q -s "my_ADC_map" --units "10^-6 mm^2/s" adc.nii.gz donor.dcm DICOM_out/
```

**The pipeline's export workflow, on a results-only directory.** The rendering, drop folders, thresholding, orientations and series numbering live in `KUL_clinical_fmridti.sh -R`, not here. That block runs no preprocessing, and its preprocessing pre-flight checks are skipped, so it works on a directory that only holds `RESULTS/sub-{participant}/`:

```bash
cd /path/to/anything/with/RESULTS
mkdir -p RESULTS/sub-XYZ/DICOM                      # put one donor DICOM here
cp my_map.nii.gz RESULTS/sub-XYZ/PACS_input/series_quantitative/
KUL_clinical_fmridti.sh -p XYZ -t 1 -R 4 -O SAG
```

The drop folders are created on the first `-R`/`-F` run (and at the end of a normal pipeline run), each with a `README.txt`. Files must already be registered to the anatomical — they inherit the donor's Frame of Reference.

## Options

| Option | Description |
|---|---|
| `-v, --verbose` | Print all donor DICOM metadata tags before conversion. |
| `-s, --seriesdescription` | Series description written to the output DICOMs (default: `IKTsimple - KUL_NIS`). |
| `-n, --seriesnumber` | Series number written to the output DICOMs. |
| `-p, --pixelspacing MM` | In-plane pixel spacing in mm. Required for PNG input without `-u/--underlay`; corrects PACS distance measurements. |
| `-u, --underlay NIFTI` | NIfTI underlay used for `mrview` rendering; supplies correct spatial metadata (position/orientation) for PACS linking and multi-planar reconstruction. |
| `-o, --plane {TRA,SAG,COR}` | Orientation plane of the PNG screenshots. Required when `-u/--underlay` is given. |
| `-M, --match-donor` | SAG mode only: match each PNG to the spatially closest donor frame and copy its exact IPP/IOP/PixelSpacing, guaranteeing PACS alignment without geometry computation. |
| `-q, --quantitative` | NIfTI input: write a **measurable** series. See below. |
| `--label` | NIfTI input: integer label/segmentation map, written with no rescaling so each label keeps its value. Implies `-q`. |
| `--units TEXT` | Value units for `-q`, written to RescaleType (0028,1054), e.g. `ml/100g` or `ratio`. Default `US` (unspecified). |
| `--window-headroom F` | Multiply the top of the default display window by `F` (`-q` only; default 1.0 = plain p99.5). Raise it for maps whose bright tail is anatomy you want out of saturation — the pipeline passes 1.3 for CBV maps, whose tail is choroid plexus. |

## Quantitative series (`-q`)

Without `-q`, a NIfTI is rescaled to fill the int16 range and the scale factor is
thrown away, so the DICOM carries arbitrary counts: an ROI drawn on PACS reads a
number with no physical meaning.

`-q` instead records the mapping, so `stored × RescaleSlope + RescaleIntercept`
recovers the original value and PACS reports real units — the same thing the
scanner's own ADC maps give you. It also switches the SOP Class to **MR Image
Storage**, because Secondary Capture has no Modality LUT module and viewers
would ignore the rescale tags entirely.

```
KUL_nii2dcm.py -q --units ratio -s rCBV -n 10120 rCBV.nii.gz donor.dcm DICOM_out/
```

Details worth knowing:

- The intercept is fixed at 0 and the whole mapping carried by the slope. GDCM
  (which SimpleITK writes through) treats these tags as an instruction to
  inverse-rescale the pixel data and refuses any non-integer slope or intercept,
  so the tags are attached in a post-pass with **pydicom** — the one mode that
  needs it. Anchoring at zero also keeps zero meaning zero for a masked map.
- Precision is 1 part in 32767 of the largest magnitude present.
- Non-finite voxels are mapped to 0 (or the minimum, if 0 is outside the range)
  and reported. Values are rounded and clipped, never wrapped.
- A default WindowCenter/WindowWidth is written in **rescaled (real) units**, as
  DICOM requires — they are applied after the Modality LUT (PS3.3 C.11.2). Taken
  from the 2nd–99.5th percentile of the **foreground** (non-zero) voxels, so a
  brain-masked map — typically ~87% background — opens on its tissue rather than
  its mask, and the very high choroid-plexus/CSF voxels on a DSC map do not blow
  out the ventricles. It is only a starting point; window on PACS. `--label`
  spans the full label range instead, since labels are categorical.
- Do **not** pass `-u`/`-o` with `-q`; they are ignored with a note. Geometry
  comes from the volume itself, because the underlay path reproduces `mrview`'s
  screenshot display convention, which is correct for rendered PNGs and wrong
  for voxel data.

Positional arguments: `nifti` (nifti/3d-tiff file, or a directory of PNG slices), `donor` (a donor DICOM file supplying patient/study metadata and Frame of Reference UID), `dicomdir` (output directory; cleared and recreated on each run).

### Donor formats

Both layouts are handled: **classic** (one frame per file, a directory of them)
and **enhanced/multi-frame** (a single file holding all frames, which SimpleITK
returns as a 4-D image — it is collapsed to 3-D before any geometry is read).
Exercised so far on Philips classic (Achieva dStream, 200 files) and Philips
enhanced (SmartBrain, 100 frames in one file), and on Siemens. **GE has not been
tested**; it should take the same path, but that is an expectation rather than a
result.

Donate a slice from a **high-resolution anatomical** series, not a localiser.
Geometry is derived from the underlay rather than the donor, so a mismatched
donor FOV is no longer fatal, but the donor is still what ties the output to the
patient's study.

## Input handling

- **NIfTI + `-q`/`--label`**: written as a measurable series with the value mapping recorded — see [Quantitative series](#quantitative-series--q) above.
- **NIfTI/3D-TIFF, no `-q`**: rescaled to fill the int16 range and reoriented to LPS; geometry comes from the input volume itself. The scale factor is *not* recoverable from the output, so use `-q` for anything anyone will measure.
- **PNG directory + `-u`/`-o` (TRA or COR)**: per-slice position/orientation is computed geometrically from the underlay NIfTI's direction matrix, with the rendering margin around the `mrview` screenshot content cropped out based on FOV vs. PNG dimensions.
- **PNG directory + `-M` (SAG)**: each PNG is matched to the nearest donor frame by slice location and inherits that frame's exact spatial metadata instead of computed geometry.
- **PNG directory, no underlay**: falls back to unit spacing with no spatial metadata (positions only from array indices).

Donor text tags (patient name, institution, etc.) are decoded using the donor's DICOM Specific Character Set (0008,0005), covering the single-byte ISO 8859 families and UTF-8; unsupported multibyte code-extension sets fall back to Latin-1 rather than raising.

## Dependencies

| Package | Version | Used for |
|---|---|---|
| python3 | ≥ 3.8 | script requires `python3` explicitly (shebang is `#!/usr/bin/env python3`) |
| [SimpleITK](https://simpleitk.org/) | any recent | reading donor DICOM metadata/geometry, reading NIfTI/TIFF input, writing the output DICOM series (GDCM-based; no separate DICOM library needed). Handles both older SimpleITK builds (metadata returned as `bytes`) and newer 2.5.x builds (surrogate-escaped `str`) transparently. |
| [Pillow](https://python-pillow.org/) (PIL) | any recent | loading/cropping/resizing PNG input. Resolves the LANCZOS resampling filter under both the old (`Image.LANCZOS`, Pillow < 9.1) and new (`Image.Resampling.LANCZOS`, Pillow ≥ 9.1/10) APIs, so it works across distro-packaged Pillow versions without pinning. |
| [numpy](https://numpy.org/) | any recent | array/volume manipulation for PNG and TIFF paths |
| [pydicom](https://pydicom.github.io/) | ≥ 2.3 | **`-q`/`--label` only.** Attaches RescaleSlope/Intercept/Type after the pixel data is written. SimpleITK cannot: GDCM reads those tags as an instruction to inverse-rescale the pixel data on write and rejects any non-integer slope or intercept, which would quantise every fractional map to whole numbers. Every other mode runs without it. |

Everything except `-q`/`--label` goes through SimpleITK alone. All other imports (`argparse`, `sys`, `time`, `os`, `shutil`, `glob`, `hashlib`) are Python standard library.

The pipeline runs this script from a dedicated conda env so these versions are pinned rather than inherited from whatever `python3` is on PATH. [KUL_Linux_setup](https://github.com/Rad-dude/KUL_Linux_setup) creates it as part of a normal install:

```
./setup_environment.sh --only env-dicom
```

Or by hand, from the equivalent spec shipped here:

```
mamba env create -f share/envs/KUL_dicom.yml
```

`KUL_clinical_fmridti.sh` picks it up automatically (default name `KUL_dicom`, from `$KUL_DICOM_ENV`; override with `-m <env>`), and falls back to plain `python3` with a warning if it does not exist.

See also the **Requirements** section of the main [README](/README.md) and the [KUL_clinical_fmridti](/docs/KUL_clinical_fmridti/KUL_clinical_fmridti.md) doc, which invokes this script with `-u`/`-o` (or `-M` for SAG) after generating `mrview` screenshots.
