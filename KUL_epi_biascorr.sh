#!/bin/bash
# Remove the receive-coil (B1-) intensity bias from a 3D/4D EPI series.
#
# Requires ANTs (N4BiasFieldCorrection) and FSL (fslmaths, fslstats)
#
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@uzleuven.be
#
# Dense head coils without Prescan Normalize give EPI a strong spatial shading
# (64ch HeadNeck: deep brain ~40% darker than cortex). That breaks everything
# that thresholds raw EPI intensity: SPM's implicit mask (mthresh) dropped the
# whole brain centre, and SUSAN's single brightness threshold smooths dark
# regions harder than bright ones. fMRIPrep N4-corrects its boldref but not
# desc-preproc_bold, so the series itself still carries the bias.
#
# The field is estimated once by N4 on <ref_img> (default: the temporal mean of
# <in>), rescaled to unit mean inside <mask> so intensities stay on the input's
# scale, and every volume is divided by it. Only a fixed spatial scaling: voxel
# time courses, t-statistics and ratio measures (dR2*, %BOLD) are unchanged.
#
# Kept standalone (no KUL_main_functions.sh) so `set -euo pipefail` scripts such
# as KUL_fmri_denoise.sh can call it too.
version="0.1 - dd 21/10/2026"

if [ "$#" -lt 3 ]; then
    cat <<USAGE >&2

`basename $0` removes the receive-coil bias field from a 3D/4D EPI image.

Usage:

  `basename $0` <in> <mask> <out> [out_ref] [ref_img]

     in       3D or 4D EPI (nii/nii.gz)
     mask     brain mask on the same grid (e.g. fMRIPrep desc-brain_mask)
     out      corrected image (float, masked)
     out_ref  optional: the corrected reference, on the same scale as <out>
              (use it as the SUSAN USAN image)
     ref_img  optional: 3D image to estimate the field on
              (default: temporal mean of <in>)

Also writes <out stem>_biasfield.nii.gz (unit mean inside the mask).

USAGE
    exit 1
fi

in="$1"; mask="$2"; out="$3"; out_ref="${4:-}"; ref="${5:-}"
field="${out%.nii*}_biasfield.nii.gz"

tmp=$(mktemp -d "${TMPDIR:-/tmp}/KUL_epi_n4_XXXXXX") || exit 1
trap 'rm -rf "$tmp"' EXIT

if [ -z "$ref" ]; then
    ref="$tmp/ref.nii.gz"
    fslmaths "$in" -Tmean "$ref" -odt float || exit 1
fi

N4BiasFieldCorrection -d 3 -i "$ref" -x "$mask" -s 2 -b [100] \
    -c [50x50x50x50,0.0] -o ["$tmp/ref_n4.nii.gz","$tmp/field.nii.gz"] \
    || { echo "ERROR: N4BiasFieldCorrection failed on $ref" >&2; exit 1; }

# fslstats pads its output with a space; quoted, fslmaths reads that as a filename
fmean=$(fslstats "$tmp/field.nii.gz" -k "$mask" -M | tr -d '[:space:]')
fslmaths "$tmp/field.nii.gz" -div "$fmean" "$field" -odt float || exit 1
fslmaths "$in" -div "$field" -mas "$mask" "$out" -odt float || exit 1
if [ -n "$out_ref" ]; then
    fslmaths "$ref" -div "$field" -mas "$mask" "$out_ref" -odt float || exit 1
fi
