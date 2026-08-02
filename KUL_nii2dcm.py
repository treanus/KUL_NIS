#!/usr/bin/env python3

# Convert nifti to dicom given a donor dicom image
# Stefan Sunaert - 27/02/2023
# Mainly based on SimpleITK - https://simpleitk.readthedocs.io/en/master/link_DicomSeriesFromArray_docs.html
# Extended by KUL_NIS: PNG directory input, correct PixelSpacing, spatial metadata for PACS linking and MPR

import SimpleITK as sitk
import argparse
import sys
import time
import os
import shutil
import glob
import hashlib
import numpy as np
from PIL import Image

# Pillow >=9.1 moved resampling filters under Image.Resampling and (>=10) drops
# some top-level aliases; resolve once so this works across distro Pillow versions.
try:
    LANCZOS = Image.Resampling.LANCZOS
except AttributeError:
    LANCZOS = Image.LANCZOS


# DICOM Specific Character Set (0008,0005) defined term -> Python codec.
# Covers the single-byte ISO 8859 families plus UTF-8, which is everything a
# clinical donor realistically uses. Multibyte code-extension sets that need a
# stateful ISO 2022 decoder (Japanese ISO 2022 IR 87/159, Korean IR 149, and
# GB18030 escape sequences) are deliberately NOT claimed here — Arabic does not
# need them (ISO 8859-6 is single-byte; modern Arabic is UTF-8).
_DICOM_CHARSET = {
    "": "ascii", "ISO_IR 6": "ascii",
    "ISO_IR 100": "latin-1",    # Latin-1  Western European (DICOM default extension)
    "ISO_IR 101": "iso8859-2",  # Latin-2  Central European
    "ISO_IR 109": "iso8859-3",  # Latin-3
    "ISO_IR 110": "iso8859-4",  # Latin-4
    "ISO_IR 144": "iso8859-5",  # Cyrillic
    "ISO_IR 127": "iso8859-6",  # Arabic
    "ISO_IR 126": "iso8859-7",  # Greek
    "ISO_IR 138": "iso8859-8",  # Hebrew
    "ISO_IR 148": "iso8859-9",  # Latin-5  Turkish
    "ISO_IR 192": "utf-8",      # Unicode UTF-8
    "GB18030": "gb18030", "GBK": "gbk",
}


def _resolve_codec(reader):
    """Pick a decode codec from Specific Character Set (0008,0005).

    Falls back to Latin-1 (the DICOM default) when the tag is absent or maps to a
    multibyte set we don't handle, so decoding always has a sane single-byte codec
    rather than raising.
    """
    for key in ("0008|0005", "0008|0005 "):
        if reader.HasMetaDataKey(key):
            term = reader.GetMetaData(key)
            if isinstance(term, bytes):
                term = term.decode("ascii", "ignore")
            # Multi-valued (code extensions): take the first component. Normalise
            # the "ISO 2022 IR nnn" spelling to the "ISO_IR nnn" table key; for the
            # single-byte sets the underlying tables are identical.
            term = str(term).split("\\")[0].strip().replace("ISO 2022 IR", "ISO_IR").strip()
            codec = _DICOM_CHARSET.get(term)
            if codec:
                return codec
    return "latin-1"


def _dcm_str(v, codec="latin-1"):
    """Coerce a donor tag value to a str SimpleITK.SetMetaData will accept.

    A non-UTF-8 tag value (accented Latin name, Cyrillic/Greek/Arabic, ...) surfaces
    differently across SimpleITK builds: older ones return raw ``bytes``; newer ones
    (2.5.x) return a surrogate-escaped ``str`` (e.g. 'H\\udcf4pital'). SetMetaData
    rejects BOTH with 'argument 3 of type std::string const &', so a version bump
    does not fix it — the trigger is the tag content, not the version.

    Three cases:
      * bytes                      -> decode with the donor `codec`.
      * str with surrogate escapes -> SimpleITK could not decode it; recover the
                                      original bytes and decode with `codec`.
      * clean str (no surrogates)  -> SimpleITK already decoded valid UTF-8
                                      correctly (e.g. an ISO_IR 192 Arabic name);
                                      trust it as-is. This passthrough is what keeps
                                      UTF-8 content from being re-decoded — and
                                      corrupted — by a single-byte codec.
    Errors fall back to replacement so it can never raise. NUL pad is stripped.
    """
    if v is None:
        return ""
    if isinstance(v, bytes):
        return v.decode(codec, "replace").rstrip("\x00")
    s = str(v)
    if any(0xD800 <= ord(c) <= 0xDFFF for c in s):
        return s.encode("utf-8", "surrogateescape").decode(codec, "replace").rstrip("\x00")
    return s.rstrip("\x00")


# Get and check commandline
parser = argparse.ArgumentParser(description="Convert a nifti, 3d-tiff, or PNG directory to dicom given a donor dicom image",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-v", "--verbose", action="store_true", help="increase verbosity")
parser.add_argument("-s", "--seriesdescription")
parser.add_argument("-n", "--seriesnumber")
parser.add_argument("-p", "--pixelspacing", type=float, default=None,
                    help="in-plane pixel spacing in mm (required for PNG input; corrects PACS distance measurements)")
parser.add_argument("-u", "--underlay", default=None,
                    help="NIfTI underlay used for mrview rendering; provides correct spatial metadata for PACS linking and MPR")
parser.add_argument("-o", "--plane", choices=["TRA", "SAG", "COR"], default=None,
                    help="orientation plane of the screenshots (TRA/SAG/COR); required when --underlay is given")
parser.add_argument("-M", "--match-donor", action="store_true",
                    help="SAG mode: match each PNG to the spatially closest donor frame and copy its exact "
                         "IPP/IOP/PixelSpacing — guarantees PACS alignment without any geometry computation")
parser.add_argument("nifti", help="nifti, 3d-tiff, or directory of PNG slices")
parser.add_argument("donor", help="dicom donor image")
parser.add_argument("dicomdir", help="dicom output directory")
args = parser.parse_args()


# Define functions
def _axis_offset_and_dir(vec, flip, spacing, size):
    """Return (corner_offset, display_dir) for one in-plane voxel axis.

    `flip` says whether pixel index 0 along this screen axis sits at voxel
    index N-1 rather than 0 (i.e. whether display order is reversed relative
    to voxel-index order for this axis).
    """
    if flip:
        offset = [(size - 1) * spacing * v for v in vec]
        direction = [-v for v in vec]
    else:
        offset = [0.0, 0.0, 0.0]
        direction = list(vec)
    return offset, direction


def compute_slice_geometry(underlay_geom, plane, slice_idx):
    """Return (image_position, row_dir, col_dir, slice_normal, slice_thickness) for a given slice.

    mrview renders with a fixed anatomical display convention (superior/anterior
    up, patient-right on the left) regardless of how a volume's own voxel axes
    are signed or how oblique its direction matrix is. Per axis, whether pixel
    index 0 sits at voxel index 0 or N-1 depends only on that axis's own
    dominant anatomical sign — verified empirically against mrview by rendering
    synthetic phantoms with anisotropic spacing, oblique (gantry-tilt-like)
    rotation, and each axis's sign independently flipped:
      i (L-R):  flips when i_dir[0] < 0
      j (A-P):  flips when j_dir[1] < 0
      k (S-I):  flips when k_dir[2] > 0
    This one rule holds unchanged across TRA/COR/SAG since each plane just
    picks 2 of the 3 axes for its row/col roles and the third as its normal.
    """
    origin, spacing, d, sz = underlay_geom
    i_dir = [d[0], d[3], d[6]]
    j_dir = [d[1], d[4], d[7]]
    k_dir = [d[2], d[5], d[8]]
    N_i, N_j, N_k = sz[0], sz[1], sz[2]

    i_off, i_disp = _axis_offset_and_dir(i_dir, i_dir[0] < 0, spacing[0], N_i)
    j_off, j_disp = _axis_offset_and_dir(j_dir, j_dir[1] < 0, spacing[1], N_j)
    k_off, k_disp = _axis_offset_and_dir(k_dir, k_dir[2] > 0, spacing[2], N_k)

    if plane == 'TRA':
        row_off, row_dir = i_off, i_disp
        col_off, col_dir = j_off, j_disp
        slice_off = [spacing[2] * slice_idx * v for v in k_dir]
        normal, thick = list(k_dir), spacing[2]
    elif plane == 'COR':
        row_off, row_dir = i_off, i_disp
        col_off, col_dir = k_off, k_disp
        slice_off = [spacing[1] * slice_idx * v for v in j_dir]
        normal, thick = list(j_dir), spacing[1]
    else:  # SAG
        row_off, row_dir = j_off, j_disp
        col_off, col_dir = k_off, k_disp
        slice_off = [spacing[0] * slice_idx * v for v in i_dir]
        normal, thick = list(i_dir), spacing[0]

    position = [origin[x] + row_off[x] + col_off[x] + slice_off[x] for x in range(3)]
    return position, row_dir, col_dir, normal, thick


def writeSlices(series_tag_values, new_img, out_dir, i, underlay_geom=None, plane=None):
    image_slice = new_img[:, :, i]

    # Tags shared by the series
    list(map(lambda tag_value: image_slice.SetMetaData(tag_value[0], tag_value[1]),
             series_tag_values))

    # Slice-specific date/time
    image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))
    image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))
    image_slice.SetMetaData("0020|0013", str(i))  # Instance Number

    if underlay_geom is not None and plane is not None:
        pos, row_dir, col_dir, normal, thick = compute_slice_geometry(underlay_geom, plane, i)
        image_slice.SetMetaData("0020|0032", "\\".join(f"{v:.6f}" for v in pos))
        image_slice.SetMetaData("0020|0037",
            "\\".join(f"{v:.6f}" for v in row_dir + col_dir))
        # Slice Location: signed distance along slice normal from origin
        slice_loc = sum(pos[x] * normal[x] for x in range(3))
        image_slice.SetMetaData("0020|1041", f"{slice_loc:.4f}")
        image_slice.SetMetaData("0018|0050", f"{thick:.4f}")  # Slice Thickness
    else:
        image_slice.SetMetaData("0020|0032",
            "\\".join(map(str, new_img.TransformIndexToPhysicalPoint((0, 0, i)))))

    writer.SetFileName(os.path.join(out_dir, str(i).rjust(6, '0') + ".dcm"))
    writer.Execute(image_slice)


# --- set inputs and check ---
donor_dcm = args.donor
if not os.path.exists(donor_dcm):
    print(donor_dcm + ' does not exist')
    sys.exit(1)
nifti_input = args.nifti
if not os.path.exists(nifti_input):
    print(nifti_input + ' does not exist')
    sys.exit(1)
img_input, img_ext = os.path.splitext(nifti_input)
img_ext = img_ext.lower()
if os.path.isdir(nifti_input):
    input_type = 'png_dir'
    print('Assuming input is a directory of PNG slices')
elif img_ext in ('.tiff', '.tif'):
    input_type = 'tiff'
    print('Assuming input is a 3d-tiff')
else:
    input_type = 'nifti'
    print('Assuming input is nifti')
dcm_output = args.dicomdir

# set defaults
seriesdesc   = args.seriesdescription if args.seriesdescription else 'IKTsimple - KUL_NIS'
seriesnumber = args.seriesnumber if args.seriesnumber else ''

# Read the donor DICOM for patient/study metadata
reader = sitk.ImageFileReader()
reader.SetFileName(donor_dcm)
reader.LoadPrivateTagsOn()
reader.ReadImageInformation()

# Codec for decoding donor text tags, from Specific Character Set (0008,0005).
_charset = _resolve_codec(reader)

if args.verbose:
    for k in reader.GetMetaDataKeys():
        try:
            print(f'({k}) = "{reader.GetMetaData(k)}"')
        except:
            pass

tags_to_copy = [
    "0010|0010",  # Patient Name
    "0010|0020",  # Patient ID
    "0010|0030",  # Patient Birth Date
    "0010|0040",  # Patient Sex
    # Study/Series/SOP UIDs are set explicitly in series_tag_values_b — do not copy from donor.
    # SOP Class UID (0002|0002 / 0008|0016) is likewise NOT copied from the donor: these outputs
    # are DERIVED/SECONDARY RGB captures, not real acquisitions, and copying the donor's own
    # SOP Class (e.g. MR Image Storage) onto RGB pixel data makes GDCM's DICOM writer silently
    # re-derive/reset Image Orientation Patient and Pixel Spacing from the image's own (default
    # identity) geometry, discarding whatever was explicitly set via SetMetaData below.
    "0020|0010",  # Study ID
    "0020|0052",  # Frame of Reference UID — links series for cursor alignment in PACS
    "0008|0020",  # Study Date
    "0008|0022",  # Acquisition Date
    "0008|0023",  # Content Date
    "0008|0030",  # Study Time
    "0008|0032",  # Acquisition Time
    "0008|0033",  # Content Time
    "0008|0050",  # Accession Number
    "0008|0060",  # Modality
    "0008|0080",  # Institution Name
]

# Secondary Capture Image Storage: what these RGB-rendered outputs actually are,
# regardless of the donor's own (real acquisition) SOP Class.
_SC_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.7"

# Load underlay geometry for correct spatial metadata (enables PACS linking and MPR)
underlay_geom = None
if args.underlay is not None and args.plane is not None:
    if os.path.exists(args.underlay):
        ref = sitk.ReadImage(args.underlay)
        underlay_geom = (ref.GetOrigin(), ref.GetSpacing(), ref.GetDirection(), ref.GetSize())
        print(f'Spatial reference: {args.underlay} (plane={args.plane})')
    else:
        print(f'Warning: underlay {args.underlay} not found — spatial metadata will be approximate')

# -------------------------------------------------------------------------
# DONOR-MATCH MODE: SAG PNGs matched to donor frames by SliceLocation.
# Copies exact IPP/IOP/PixelSpacing from the donor — no geometry computation.
# -------------------------------------------------------------------------
if args.match_donor and input_type == 'png_dir':
    print('Donor-match mode: reading donor 3D geometry...')

    # Read donor as a full 3D series using its parent directory so that per-frame
    # spatial metadata (IPP) is available for all slices.
    # sitk.ReadImage on a single DICOM file gives a 2D image even for a multi-slice
    # series, which breaks the per-frame matching below.
    donor_dir = os.path.dirname(os.path.abspath(donor_dcm))
    series_ids = sitk.ImageSeriesReader.GetGDCMSeriesIDs(donor_dir)
    if series_ids:
        series_files = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(donor_dir, series_ids[0])
        donor_3d = sitk.ReadImage(series_files)
    else:
        donor_3d = sitk.ReadImage(donor_dcm)
        if donor_3d.GetDimension() == 2:
            donor_3d = sitk.JoinSeries(donor_3d)
    d3 = donor_3d.GetDirection()
    sp3 = donor_3d.GetSpacing()
    or3 = donor_3d.GetOrigin()
    sz3 = donor_3d.GetSize()   # (n_cols, n_rows, n_slices)

    def _norm(v):
        n = (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5
        return [x / n for x in v]

    def _dot(a, b):
        return sum(a[i] * b[i] for i in range(3))

    def _cross(a, b):
        return [a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0]]

    donor_ps        = sp3[0]          # pixel spacing (assumed isotropic in-plane)
    donor_cols      = sz3[0]          # image width  per slice
    donor_rows      = sz3[1]          # image height per slice
    n_donor_frames  = sz3[2]

    # IOP for the output DICOMs comes from the underlay NIfTI (SAG orientation),
    # NOT from the donor — the donor is acquired axially, so its IOP is axial and
    # would produce completely wrong orientation for SAG output slices.
    # The donor provides only Frame of Reference UID, patient/study metadata,
    # and (below) the frame count/order + pixel matrix size to match.
    if underlay_geom is None:
        print('ERROR: donor-match mode requires --underlay for correct SAG orientation')
        sys.exit(1)
    nifti_or_tmp, nifti_sp_tmp, nifti_d_tmp, nifti_sz_tmp = underlay_geom
    _, out_row_dir, out_col_dir, out_slice_dir, _ = compute_slice_geometry(underlay_geom, 'SAG', 0)

    # donor_slice_dir is still needed for per-frame loc computation
    donor_row_raw   = [d3[0], d3[3], d3[6]]
    donor_col_raw   = [d3[1], d3[4], d3[7]]
    donor_row_dir   = _norm(donor_row_raw)
    col_orth_d2     = [donor_col_raw[i] - _dot(donor_row_dir, donor_col_raw) * donor_row_dir[i]
                       for i in range(3)]
    donor_col_dir   = _norm(col_orth_d2)
    donor_slice_dir = _cross(donor_row_dir, donor_col_dir)

    out_iop_str = "\\".join(f"{v:.6f}" for v in out_row_dir + out_col_dir)

    print(f'Donor 3D: {n_donor_frames} frames, {donor_cols}×{donor_rows} px, ps={donor_ps:.4f}mm')
    print(f'Output IOP (from NIfTI SAG): {out_iop_str}')

    # Build per-frame table: IPP and SliceLocation
    donor_frames = []
    for j in range(n_donor_frames):
        ipp = [or3[x] + donor_slice_dir[x] * sp3[2] * j for x in range(3)]
        loc = sum(ipp[x] * donor_slice_dir[x] for x in range(3))
        donor_frames.append({'ipp': ipp, 'loc': loc})
    donor_frames.sort(key=lambda f: f['loc'])

    # Build PNG table sorted by SliceLocation along the donor slice normal
    pngs = sorted(glob.glob(os.path.join(nifti_input, '*.png')))
    if not pngs:
        print(f'No PNG files found in {nifti_input}')
        sys.exit(1)

    if underlay_geom is not None:
        nifti_or, nifti_sp, nifti_d, _ = underlay_geom
        nifti_i_dir = [nifti_d[0], nifti_d[3], nifti_d[6]]
        base_loc  = sum(nifti_or[x] * donor_slice_dir[x] for x in range(3))
        loc_step  = sum(nifti_i_dir[x] * donor_slice_dir[x] for x in range(3)) * nifti_sp[0]
        png_table = sorted(
            [(i, base_loc + loc_step * i, pngs[i]) for i in range(len(pngs))],
            key=lambda t: t[1])
    else:
        png_table = [(i, i, pngs[i]) for i in range(len(pngs))]

    if len(png_table) != n_donor_frames:
        print(f'Warning: {len(png_table)} PNGs vs {n_donor_frames} donor frames — '
              f'using nearest-neighbour matching')

    # Compute rendering margin geometrically from the in-plane FOV and PNG dimensions.
    # Pixel-content thresholds are unreliable: both the rendering margin and the dark
    # brain exterior are near-zero with a black background (or both near-white with white).
    # mrview's mode-1 view uses ONE zoom for the whole volume (so anatomy doesn't visually
    # resize when scrolling between orthogonal planes), based on the single largest
    # physical extent across ALL 3 volume axes -- not just this plane's own 2 in-plane
    # dimensions. Verified empirically against mrview with synthetic phantoms.
    sample_img  = Image.open(pngs[len(pngs) // 2])
    PNG_W, PNG_H = sample_img.size   # width × height in pixels
    j_fov = donor_cols * donor_ps    # horizontal in-plane FOV (mm) for SAG = n_j * sp
    k_fov = donor_rows * donor_ps    # vertical   in-plane FOV (mm) for SAG = n_k * sp
    _global_max_fov = max(nifti_sz_tmp[0] * nifti_sp_tmp[0],
                          nifti_sz_tmp[1] * nifti_sp_tmp[1],
                          nifti_sz_tmp[2] * nifti_sp_tmp[2])
    scale = min(PNG_W, PNG_H) / _global_max_fov   # px/mm, matches mrview's fixed zoom
    content_W = round(j_fov * scale)             # expected rendered content width  (px)
    content_H = round(k_fov * scale)             # expected rendered content height (px)
    crop_left  = (PNG_W - content_W) // 2
    crop_right = crop_left + content_W           # exclusive end column
    crop_top   = (PNG_H - content_H) // 2
    crop_bot   = crop_top + content_H            # exclusive end row
    print(f'Geometric crop: PNG={PNG_W}×{PNG_H}  FOV={j_fov:.1f}×{k_fov:.1f}mm  '
          f'scale={scale:.3f}px/mm  content={content_W}×{content_H}px  '
          f'crop cols {crop_left}:{crop_right} rows {crop_top}:{crop_bot}  '
          f'→ resizing to {donor_cols}×{donor_rows}')

    # Common series tags (patient/study from donor reader)
    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")
    series_tag_values_a = [
        (k, _dcm_str(reader.GetMetaData(k), _charset))
        for k in tags_to_copy
        if reader.HasMetaDataKey(k)
    ]
    _series_hash = str(int(hashlib.md5(seriesdesc.encode()).hexdigest()[:8], 16))
    _study_uid = _dcm_str(reader.GetMetaData("0020|000d"), _charset) if reader.HasMetaDataKey("0020|000d") else \
                 _dcm_str(reader.GetMetaData("0020|000D"), _charset) if reader.HasMetaDataKey("0020|000D") else \
                 "1.2.826.0.1.3680043.2.1125." + modification_date + modification_time
    series_tag_values_b = [
        ("0002|0002", _SC_SOP_CLASS),
        ("0008|0016", _SC_SOP_CLASS),
        ("0008|0031", modification_time),
        ("0008|0021", modification_date),
        ("0008|0008", "DERIVED\\SECONDARY"),
        ("0020|000d", _study_uid),
        ("0020|000e", "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time + "." + _series_hash),
        ("0008|103e", seriesdesc),
        ("0020|0011", seriesnumber),
    ]
    series_tags = series_tag_values_a + series_tag_values_b

    # Clean output dir
    if os.path.exists(dcm_output):
        shutil.rmtree(dcm_output)
    os.makedirs(dcm_output, exist_ok=True)

    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()

    for out_idx, (png_idx, png_loc, png_path) in enumerate(png_table):
        # Load, crop rendering margin, resize to donor dimensions
        img_arr = np.array(Image.open(png_path).convert('RGB'))
        img_arr = img_arr[crop_top:crop_bot, crop_left:crop_right, :]
        img_resized = np.array(
            Image.fromarray(img_arr).resize((donor_cols, donor_rows), LANCZOS))
        # IPP = physical position of top-left pixel for this SAG slice, using
        # the same per-axis flip convention as out_row_dir/out_col_dir above
        # (verified against mrview's actual rendering, not assumed).
        position, _, _, normal, _ = compute_slice_geometry(underlay_geom, 'SAG', png_idx)

        slice_2d = sitk.GetImageFromArray(img_resized, isVector=True)
        slice_2d.SetSpacing([donor_ps, donor_ps])

        ipp_str = "\\".join(f"{v:.6f}" for v in position)
        slice_loc = sum(position[x] * normal[x] for x in range(3))

        for tag, val in series_tags:
            slice_2d.SetMetaData(tag, val)
        slice_2d.SetMetaData("0020|0037", out_iop_str)
        slice_2d.SetMetaData("0020|0032", ipp_str)
        slice_2d.SetMetaData("0020|1041", f"{slice_loc:.4f}")
        slice_2d.SetMetaData("0018|0050", f"{sp3[2]:.4f}")
        slice_2d.SetMetaData("0018|0088", f"{sp3[2]:.4f}")
        slice_2d.SetMetaData("0028|0030", f"{donor_ps:.6f}\\{donor_ps:.6f}")
        slice_2d.SetMetaData("0020|0013", str(out_idx))
        slice_2d.SetMetaData("0008|0012", modification_date)
        slice_2d.SetMetaData("0008|0013", modification_time)
        slice_2d.SetMetaData("0008|0018",
            "1.2.826.0.1.3680043.2.1125." + modification_date + modification_time + "." + str(int(hashlib.md5(seriesdesc.encode()).hexdigest()[:8], 16)) + "." + str(out_idx))

        out_path = os.path.join(dcm_output, str(out_idx).rjust(6, '0') + ".dcm")
        writer.SetFileName(out_path)
        writer.Execute(slice_2d)

    print(f'Wrote {len(png_table)} matched SAG DICOMs to {dcm_output}')
    sys.exit(0)

# -------------------------------------------------------------------------
# TRA/COR PNG MODE: corrected geometry from NIfTI + geometric crop/resize
# -------------------------------------------------------------------------
if input_type == 'png_dir' and underlay_geom is not None and args.plane in ('TRA', 'COR'):
    nifti_or, nifti_sp, nifti_d, nifti_sz = underlay_geom
    N_i, N_j, N_k = nifti_sz[0], nifti_sz[1], nifti_sz[2]

    if args.plane == 'TRA':
        h_fov, v_fov    = N_i * nifti_sp[0], N_j * nifti_sp[1]
        out_w,  out_h   = N_i, N_j
        out_sp_x, out_sp_y = nifti_sp[0], nifti_sp[1]
    else:  # COR
        h_fov, v_fov    = N_i * nifti_sp[0], N_k * nifti_sp[2]
        out_w,  out_h   = N_i, N_k
        out_sp_x, out_sp_y = nifti_sp[0], nifti_sp[2]

    pngs = sorted(glob.glob(os.path.join(nifti_input, '*.png')))
    if not pngs:
        print(f'No PNG files found in {nifti_input}')
        sys.exit(1)

    sample_img = Image.open(pngs[len(pngs) // 2])
    PNG_W, PNG_H = sample_img.size
    # mrview's mode-1 view uses ONE zoom for the whole volume (so anatomy doesn't
    # visually resize when scrolling between orthogonal planes), based on the single
    # largest physical extent across ALL 3 volume axes -- not just this plane's own
    # 2 in-plane dimensions. Verified empirically against mrview with synthetic phantoms.
    _global_max_fov = max(N_i * nifti_sp[0], N_j * nifti_sp[1], N_k * nifti_sp[2])
    scale     = min(PNG_W, PNG_H) / _global_max_fov
    content_W = round(h_fov * scale)
    content_H = round(v_fov * scale)
    crop_left  = (PNG_W - content_W) // 2
    crop_right = crop_left + content_W
    crop_top   = (PNG_H - content_H) // 2
    crop_bot   = crop_top  + content_H
    print(f'Geometric crop ({args.plane}): PNG={PNG_W}×{PNG_H}  FOV={h_fov:.1f}×{v_fov:.1f}mm  '
          f'scale={scale:.3f}px/mm  content={content_W}×{content_H}px  '
          f'crop cols {crop_left}:{crop_right} rows {crop_top}:{crop_bot}  '
          f'→ resizing to {out_w}×{out_h}')

    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")
    series_tag_values_a = [
        (k, _dcm_str(reader.GetMetaData(k), _charset)) for k in tags_to_copy if reader.HasMetaDataKey(k)
    ]
    _series_hash2 = str(int(hashlib.md5(seriesdesc.encode()).hexdigest()[:8], 16))
    _study_uid2 = _dcm_str(reader.GetMetaData("0020|000d"), _charset) if reader.HasMetaDataKey("0020|000d") else \
                  _dcm_str(reader.GetMetaData("0020|000D"), _charset) if reader.HasMetaDataKey("0020|000D") else \
                  "1.2.826.0.1.3680043.2.1125." + modification_date + modification_time
    series_tag_values_b = [
        ("0002|0002", _SC_SOP_CLASS),
        ("0008|0016", _SC_SOP_CLASS),
        ("0008|0031", modification_time),
        ("0008|0021", modification_date),
        ("0008|0008", "DERIVED\\SECONDARY"),
        ("0020|000d", _study_uid2),
        ("0020|000e", "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time + "." + _series_hash2),
        ("0008|103e", seriesdesc),
        ("0020|0011", seriesnumber),
    ]
    series_tags = series_tag_values_a + series_tag_values_b

    if os.path.exists(dcm_output):
        shutil.rmtree(dcm_output)
    os.makedirs(dcm_output, exist_ok=True)

    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()

    for i, png_path in enumerate(pngs):
        position, row_dir, col_dir, normal, thick = compute_slice_geometry(
            underlay_geom, args.plane, i)
        iop_str   = "\\".join(f"{v:.6f}" for v in row_dir + col_dir)
        ipp_str   = "\\".join(f"{v:.6f}" for v in position)
        slice_loc = sum(position[x] * normal[x] for x in range(3))

        img_arr    = np.array(Image.open(png_path).convert('RGB'))
        img_arr    = img_arr[crop_top:crop_bot, crop_left:crop_right, :]
        img_resized = np.array(
            Image.fromarray(img_arr).resize((out_w, out_h), LANCZOS))

        slice_img = sitk.GetImageFromArray(img_resized, isVector=True)
        slice_img.SetSpacing([out_sp_x, out_sp_y])

        for tag, val in series_tags:
            slice_img.SetMetaData(tag, val)
        slice_img.SetMetaData("0020|0037", iop_str)
        slice_img.SetMetaData("0020|0032", ipp_str)
        slice_img.SetMetaData("0020|1041", f"{slice_loc:.4f}")
        slice_img.SetMetaData("0018|0050", f"{thick:.4f}")
        slice_img.SetMetaData("0018|0088", f"{thick:.4f}")
        slice_img.SetMetaData("0028|0030", f"{out_sp_y:.6f}\\{out_sp_x:.6f}")
        slice_img.SetMetaData("0020|0013", str(i))
        slice_img.SetMetaData("0008|0012", modification_date)
        slice_img.SetMetaData("0008|0013", modification_time)
        # Explicit unique SOP Instance UID per slice — without this, nothing
        # sets one on these freshly-created images and KeepOriginalImageUIDOn()
        # has nothing to "keep", so the writer falls back to auto-generating
        # one per call. That fallback isn't guaranteed unique across a tight
        # loop of many slices (this is what the SAG donor-match path already
        # does explicitly below/above); a collision here means PACS silently
        # drops/overwrites same-UID slices, leaving the series too short for
        # MPR even though each individual slice looks fine on its own.
        slice_img.SetMetaData("0008|0018",
            "1.2.826.0.1.3680043.2.1125." + modification_date + modification_time + "." + _series_hash2 + "." + str(i))

        out_path = os.path.join(dcm_output, str(i).rjust(6, '0') + ".dcm")
        writer.SetFileName(out_path)
        writer.Execute(slice_img)

    print(f'Wrote {len(pngs)} {args.plane} DICOMs to {dcm_output}')
    sys.exit(0)

# -------------------------------------------------------------------------
# Standard mode (NIfTI / TIFF / PNG with computed geometry)
# -------------------------------------------------------------------------

# Read input image
if input_type == 'png_dir':
    pngs = sorted(glob.glob(os.path.join(nifti_input, '*.png')))
    if not pngs:
        print(f'No PNG files found in {nifti_input}')
        sys.exit(1)
    print(f'Reading {len(pngs)} PNG slices from {nifti_input}')
    frames = [np.array(Image.open(p).convert('RGB')) for p in pngs]
    volume = np.stack(frames, axis=0)           # (Z, H, W, 3)
    new_img = sitk.GetImageFromArray(volume, isVector=True)
    new_img.SetSpacing([1.0, 1.0, 1.0])
elif input_type == 'tiff':
    nii_img = sitk.ReadImage(nifti_input)
    new_img = nii_img
else:
    print('Converting the nifti to 16bit')
    nii_img = sitk.ReadImage(nifti_input)
    img_data = sitk.GetArrayFromImage(nii_img)
    max_val = np.amax(img_data)
    if max_val == 0:
        print('Error: input image is all zeros; cannot rescale to int16')
        sys.exit(1)
    img_int16 = (img_data * (np.iinfo(np.int16).max / max_val)).astype(np.int16)
    new_img = sitk.GetImageFromArray(img_int16)
    new_img.CopyInformation(nii_img)
    new_img = sitk.DICOMOrient(new_img, "LPS")

# Apply correct pixel spacing
if args.pixelspacing is not None:
    ps = args.pixelspacing
    new_img.SetSpacing([ps, ps, 1.0])
    print(f'Pixel spacing set to {ps:.4f} mm')

writer = sitk.ImageFileWriter()
writer.KeepOriginalImageUIDOn()

modification_time = time.strftime("%H%M%S")
modification_date = time.strftime("%Y%m%d")

# Image Orientation: use underlay geometry if available, else fall back to new_img direction
if underlay_geom is not None:
    _, row_dir, col_dir, _, _ = compute_slice_geometry(underlay_geom, args.plane, 0)
    orientation_str = "\\".join(f"{v:.6f}" for v in row_dir + col_dir)
else:
    d = new_img.GetDirection()
    orientation_str = "\\".join(map(str, (d[0], d[3], d[6], d[1], d[4], d[7])))

series_tag_values_a = [
    (k, _dcm_str(reader.GetMetaData(k), _charset))
    for k in tags_to_copy
    if reader.HasMetaDataKey(k)
]
series_tag_values_b = [
    ("0002|0002", _SC_SOP_CLASS),
    ("0008|0016", _SC_SOP_CLASS),
    ("0008|0031", modification_time),
    ("0008|0021", modification_date),
    ("0008|0008", "DERIVED\\SECONDARY"),
    ("0020|000e",
     "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time + "." + str(int(hashlib.md5(seriesdesc.encode()).hexdigest()[:8], 16))),
    ("0020|0037", orientation_str),
    ("0008|103e", seriesdesc),
    ("0020|0011", seriesnumber),
]
series_tag_values = series_tag_values_a + series_tag_values_b

print('Incorporating the following dicom tags:')
print(series_tag_values)

# Clean and make the output dir
if os.path.exists(dcm_output):
    shutil.rmtree(dcm_output)
os.makedirs(dcm_output, exist_ok=True)

# Write slices
list(map(
    lambda i: writeSlices(series_tag_values, new_img, dcm_output, i,
                          underlay_geom=underlay_geom, plane=args.plane),
    range(new_img.GetDepth())
))

sys.exit(0)
