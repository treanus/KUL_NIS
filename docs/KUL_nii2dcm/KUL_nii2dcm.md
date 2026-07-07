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

Positional arguments: `nifti` (nifti/3d-tiff file, or a directory of PNG slices), `donor` (a donor DICOM file supplying patient/study metadata and Frame of Reference UID), `dicomdir` (output directory; cleared and recreated on each run).

## Input handling

- **NIfTI/3D-TIFF**: rescaled to int16 and reoriented to LPS; geometry comes from the input volume itself.
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

No `pydicom` or other DICOM library is required — all DICOM I/O goes through SimpleITK. All other imports (`argparse`, `sys`, `time`, `os`, `shutil`, `glob`, `hashlib`) are Python standard library.

See also the **Requirements** section of the main [README](/README.md) and the [KUL_clinical_fmridti](/docs/KUL_clinical_fmridti/KUL_clinical_fmridti.md) doc, which invokes this script with `-u`/`-o` (or `-M` for SAG) after generating `mrview` screenshots.
