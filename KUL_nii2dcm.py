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
parser.add_argument("-q", "--quantitative", action="store_true",
                    help="NIfTI input: write a measurable MR Image Storage series. Voxel values are mapped "
                         "to int16 and the mapping is recorded in RescaleSlope/RescaleIntercept, so an ROI "
                         "drawn on PACS reads the original units (rCBV, ALFF, ReHo, FA, ...). Without this "
                         "the values are rescaled to fill int16 and the scale factor is lost.")
parser.add_argument("--label", action="store_true",
                    help="NIfTI input: integer label/segmentation map. Values are written through unchanged "
                         "(slope 1, intercept 0) so each label keeps its identity. Implies -q.")
parser.add_argument("--units", default=None,
                    help="value units for -q, written to RescaleType (0028,1054), e.g. 'ml/100g' or 'ratio'. "
                         "Default 'US' (unspecified).")
parser.add_argument("--window-headroom", type=float, default=1.0,
                    help="multiply the top of the default display window by this factor (-q only). "
                         "1.0 = the plain p99.5 of foreground. Raise it (e.g. 1.3) for maps whose "
                         "bright tail is anatomy you want to keep out of saturation, such as the "
                         "choroid plexus on a CBV map.")
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


def _attach_modality_lut(out_dir, slope, intercept, units, window_center, window_width):
    """Add RescaleSlope/Intercept/Type (+ default window) to an written series.

    These cannot go through SimpleITK: GDCM interprets them as a request to
    inverse-rescale the pixel data while writing, and rejects non-integer slope
    or intercept. So the series is written first with the stored values it
    already has, and the Modality LUT is attached here as a pure metadata edit.

    Without these tags a PACS ROI reports raw stored integers instead of rCBV /
    ALFF / FA, which is the whole reason -q exists -- so a missing pydicom is a
    hard error rather than a silent downgrade.
    """
    try:
        import pydicom
    except ImportError:
        print('Error: -q/--label needs pydicom to attach RescaleSlope/RescaleIntercept.')
        print('       Without them PACS would report raw stored values, not real units.')
        print('       Install it into the DICOM env:')
        print('         mamba env create -f <KUL_NIS>/share/envs/KUL_dicom.yml')
        print('       (or: pip install pydicom), then re-run.')
        return False

    # Assigned as pre-formatted strings: left to itself pydicom writes the full
    # float repr, which overflows DS's 16-byte limit for a small slope.
    slope_s, intercept_s = _ds16(slope), _ds16(intercept)
    for f in sorted(glob.glob(os.path.join(out_dir, '*.dcm'))):
        ds = pydicom.dcmread(f)
        ds.RescaleSlope = slope_s
        ds.RescaleIntercept = intercept_s
        ds.RescaleType = units
        if window_center is not None:
            ds.WindowCenter = _ds16(window_center)
            ds.WindowWidth = _ds16(window_width)
        ds.save_as(f)
    return True


def _ds16(v):
    """Format a float as DS keeping as many significant digits as fit in 16 bytes.

    Used for RescaleSlope, where _ds()'s fixed 6 decimals would be ruinous: a
    slope of 0.00011174462045 would round to 0.000112, a 0.2% error on every
    voxel. Left to itself Python emits 0.00011174462045097997 (22 chars), which
    overflows DS.
    """
    v = float(v)
    for prec in range(15, 0, -1):
        s = f"{v:.{prec}g}"
        if len(s) <= 16:
            return s
    return f"{v:.6g}"[:16]


def _ds(v):
    """Format a float for a DICOM DS element.

    DS is limited to 16 bytes. Python's repr (what `str()` gives) can emit 18+
    characters for an ordinary oblique direction cosine -- '0.9999999999999998'
    -- which strict PACS and validators reject. Six decimals is well inside the
    limit and far finer than any geometry we have.
    """
    return f"{v:.6f}"


def writeSlices(series_tag_values, new_img, out_dir, i, underlay_geom=None, plane=None,
                series_uid=None):
    image_slice = new_img[:, :, i]

    # Tags shared by the series
    list(map(lambda tag_value: image_slice.SetMetaData(tag_value[0], tag_value[1]),
             series_tag_values))

    # Slice-specific date/time
    image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))
    image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))
    # Instance Number is 1-based by convention; 0 makes some viewers mislabel or
    # mis-sort the first slice.
    image_slice.SetMetaData("0020|0013", str(i + 1))

    # SOP Instance UID, derived from the series UID so it is unique and stable.
    # writer.KeepOriginalImageUIDOn() has nothing to keep on a freshly created
    # image, and the writer's own generation is not guaranteed unique across a
    # tight loop of many slices -- duplicates make PACS silently drop slices.
    # (The two PNG paths already do this; the standard path never did.)
    if series_uid is not None:
        image_slice.SetMetaData("0008|0018", f"{series_uid}.{i + 1}")

    if underlay_geom is not None and plane is not None:
        pos, row_dir, col_dir, normal, thick = compute_slice_geometry(underlay_geom, plane, i)
        image_slice.SetMetaData("0020|0032", "\\".join(_ds(v) for v in pos))
        image_slice.SetMetaData("0020|0037",
            "\\".join(_ds(v) for v in row_dir + col_dir))
        # Slice Location: signed distance along slice normal from origin
        slice_loc = sum(pos[x] * normal[x] for x in range(3))
        image_slice.SetMetaData("0020|1041", _ds(slice_loc))
        image_slice.SetMetaData("0018|0050", _ds(thick))  # Slice Thickness
    else:
        # Geometry straight from the volume's own direction matrix. This is the
        # correct path for actual voxel data: compute_slice_geometry above
        # reproduces mrview's *display* flip convention, which is right for
        # screenshots and wrong for a scalar map.
        sp = new_img.GetSpacing()
        d = new_img.GetDirection()
        row_dir = (d[0], d[3], d[6])   # direction of increasing column index
        col_dir = (d[1], d[4], d[7])   # direction of increasing row index
        normal = (d[2], d[5], d[8])
        pos = new_img.TransformIndexToPhysicalPoint((0, 0, i))
        image_slice.SetMetaData("0020|0032", "\\".join(_ds(v) for v in pos))
        image_slice.SetMetaData("0020|0037",
            "\\".join(_ds(v) for v in list(row_dir) + list(col_dir)))
        image_slice.SetMetaData("0020|1041", _ds(sum(pos[x] * normal[x] for x in range(3))))
        # PixelSpacing is [between rows, between columns] = [dy, dx], the
        # opposite order to SimpleITK's (x, y, z) spacing.
        image_slice.SetMetaData("0028|0030", f"{_ds(sp[1])}\\{_ds(sp[0])}")
        image_slice.SetMetaData("0018|0050", _ds(sp[2]))   # Slice Thickness
        image_slice.SetMetaData("0018|0088", _ds(sp[2]))   # Spacing Between Slices

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

# MR Image Storage, used only by -q/--label. Secondary Capture has no Modality
# LUT module, so RescaleSlope/RescaleIntercept are not part of its IOD and PACS
# ROI tools will not apply them — which is exactly what a measurable series
# needs. The concern documented above (donor SOP Class making GDCM re-derive
# geometry) does not apply here: this path writes single-channel int16 with
# geometry taken from the image's own direction matrix, not RGB captures.
_MR_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.4"

# --label is a special case of -q (no rescaling rather than a computed slope).
if args.label:
    args.quantitative = True

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

    # Force the donor to exactly 3 dimensions before any geometry is read.
    #
    # A single *multi-frame* (enhanced) DICOM -- what Philips exports, e.g. a
    # 100-frame SmartBrain localiser in one file -- reads back as 4-D
    # (cols, rows, frames, 1), so GetDirection() returns a 4x4 matrix. The
    # column extraction below indexes it as 3x3 (d3[0], d3[3], d3[6]), which on
    # a 4x4 picks up (0,0,0); normalising that divided by zero and killed the
    # conversion. SAG is the only orientation routed through donor-match mode,
    # so this presented as "sagittal DICOMs fail, the others are fine".
    while donor_3d.GetDimension() > 3:
        donor_3d = donor_3d[..., 0]
    if donor_3d.GetDimension() == 2:
        donor_3d = sitk.JoinSeries(donor_3d)
    d3 = donor_3d.GetDirection()
    sp3 = donor_3d.GetSpacing()
    or3 = donor_3d.GetOrigin()
    sz3 = donor_3d.GetSize()   # (n_cols, n_rows, n_slices)

    def _norm(v):
        n = (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5
        if n == 0:
            # Should be unreachable now the donor is forced to 3-D above, but a
            # degenerate direction must not surface as a bare ZeroDivisionError.
            print('ERROR: donor has a degenerate direction matrix '
                  f'({donor_3d.GetDimension()}-D, direction {d3}).')
            print('       Cannot derive geometry from this donor; use a different one.')
            sys.exit(1)
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
        # Informational only: every PNG below is written as its own slice
        # regardless of the donor's frame count -- geometry comes entirely from
        # the underlay (compute_slice_geometry), so a count mismatch here has
        # no effect on correctness.
        print(f'Note: {len(png_table)} PNGs vs {n_donor_frames} donor frames '
              f'(informational only -- all {len(png_table)} PNGs are written)')

    # Compute rendering margin geometrically from the in-plane FOV and PNG dimensions.
    # Pixel-content thresholds are unreliable: both the rendering margin and the dark
    # brain exterior are near-zero with a black background (or both near-white with white).
    # mrview's mode-1 view uses ONE zoom for the whole volume (so anatomy doesn't visually
    # resize when scrolling between orthogonal planes), based on the single largest
    # physical extent across ALL 3 volume axes -- not just this plane's own 2 in-plane
    # dimensions. Verified empirically against mrview with synthetic phantoms.
    sample_img  = Image.open(pngs[len(pngs) // 2])
    PNG_W, PNG_H = sample_img.size   # width × height in pixels

    # In-plane geometry comes from the UNDERLAY, not the donor.
    #
    # The PNGs are renders of the underlay, so their in-plane extent is the
    # underlay's: for SAG, j (A-P) horizontally and k (S-I) vertically. These
    # used to be taken from the donor (donor_cols * donor_ps), which is only
    # correct when the donor happens to share the underlay's field of view.
    # With, say, a 320x320 @ 1.09mm SmartBrain localiser (350mm FOV) donating
    # for a 64x64x64 @ 2mm underlay (128mm FOV), every SAG series was written
    # claiming a 350mm in-plane extent for 128mm of anatomy -- a 2.7x stretch,
    # clipped at the edge of the frame. TRA/COR never had this because they
    # derive all of it from the underlay; SAG now does the same.
    # The donor still supplies identity, Frame of Reference and frame matching.
    j_fov = nifti_sz_tmp[1] * nifti_sp_tmp[1]   # horizontal in-plane FOV (mm), A-P
    k_fov = nifti_sz_tmp[2] * nifti_sp_tmp[2]   # vertical   in-plane FOV (mm), S-I
    out_w, out_h = nifti_sz_tmp[1], nifti_sz_tmp[2]
    out_sp_x, out_sp_y = nifti_sp_tmp[1], nifti_sp_tmp[2]
    out_slice_sp = nifti_sp_tmp[0]              # sagittal step = underlay i spacing
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
          f'→ resizing to {out_w}×{out_h} @ {out_sp_x:.4f}×{out_sp_y:.4f}mm')

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
            Image.fromarray(img_arr).resize((out_w, out_h), LANCZOS))
        # IPP = physical position of top-left pixel for this SAG slice, using
        # the same per-axis flip convention as out_row_dir/out_col_dir above
        # (verified against mrview's actual rendering, not assumed).
        position, _, _, normal, _ = compute_slice_geometry(underlay_geom, 'SAG', png_idx)

        slice_2d = sitk.GetImageFromArray(img_resized, isVector=True)
        slice_2d.SetSpacing([out_sp_x, out_sp_y])

        ipp_str = "\\".join(f"{v:.6f}" for v in position)
        slice_loc = sum(position[x] * normal[x] for x in range(3))

        for tag, val in series_tags:
            slice_2d.SetMetaData(tag, val)
        slice_2d.SetMetaData("0020|0037", out_iop_str)
        slice_2d.SetMetaData("0020|0032", ipp_str)
        slice_2d.SetMetaData("0020|1041", f"{slice_loc:.4f}")
        # Slice thickness/spacing are the underlay's sagittal step, not the
        # donor's (sp3[2]) -- same reason as the in-plane geometry above.
        slice_2d.SetMetaData("0018|0050", f"{out_slice_sp:.4f}")
        slice_2d.SetMetaData("0018|0088", f"{out_slice_sp:.4f}")
        # PixelSpacing is [between rows, between columns] = [dy, dx].
        slice_2d.SetMetaData("0028|0030", f"{out_sp_y:.6f}\\{out_sp_x:.6f}")
        # Instance Number is 1-based by convention; 0 mis-sorts in some viewers
        # (same fix as writeSlices(), never propagated here until now).
        slice_2d.SetMetaData("0020|0013", str(out_idx + 1))
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
        # Instance Number is 1-based by convention; 0 mis-sorts in some viewers
        # (same fix as writeSlices(), never propagated here until now).
        slice_img.SetMetaData("0020|0013", str(i + 1))
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
    nii_img = sitk.ReadImage(nifti_input)

    # 4-D input is not supported. It used to die several lines later inside
    # CopyInformation with a SimpleITK dimension error; say so plainly instead.
    # Writing a timeseries properly means one series carrying every timepoint
    # with TemporalPositionIdentifier and per-timepoint acquisition times, which
    # this converter does not do -- see the hand-off notes.
    if nii_img.GetDimension() > 3:
        _sz = nii_img.GetSize()
        print(f'Error: {nifti_input} is {nii_img.GetDimension()}-D {_sz}; only 3-D volumes are supported.')
        print('       For a timeseries (e.g. a 4-D PWI), split it first:')
        print(f'         mrconvert "{nifti_input}" -coord 3 <index> vol.nii.gz')
        print('       and convert each volume as its own series.')
        sys.exit(1)

    img_data = sitk.GetArrayFromImage(nii_img).astype(np.float64)

    # Non-finite voxels would propagate through the arithmetic and cast to
    # arbitrary integers. Map them to 0 when 0 is inside the value range
    # (normal for a masked brain map), otherwise to the minimum.
    _nonfinite = ~np.isfinite(img_data)
    _n_nonfinite = int(_nonfinite.sum())
    _finite = img_data[~_nonfinite]
    if _finite.size == 0:
        print('Error: input image has no finite voxels')
        sys.exit(1)
    vmin, vmax = float(_finite.min()), float(_finite.max())
    if _n_nonfinite:
        _fill = 0.0 if vmin <= 0.0 <= vmax else vmin
        img_data[_nonfinite] = _fill
        print(f'Note: {_n_nonfinite} non-finite voxel(s) set to {_fill}')

    _i16 = np.iinfo(np.int16)
    if args.label:
        # Label maps must keep their exact values; any scaling destroys the
        # identity of each label.
        if vmin < _i16.min or vmax > _i16.max:
            print(f'Error: label values [{vmin}, {vmax}] do not fit in int16')
            sys.exit(1)
        if not np.allclose(img_data, np.round(img_data)):
            print('Warning: --label given but values are not integers; rounding')
        img_int16 = np.round(img_data).astype(np.int16)
        rescale_slope, rescale_intercept = 1.0, 0.0
        print(f'Label mode: {len(np.unique(img_int16))} distinct value(s), written unscaled')
    elif args.quantitative:
        # Scale about zero: stored = value / slope, intercept 0.
        #
        # The intercept is deliberately kept at exactly 0 rather than shifted to
        # pack the values into the full int16 range. GDCM's Rescaler asserts
        # `intercept == (int)intercept` when writing integer pixel data, so a
        # fractional intercept aborts the write outright:
        #   gdcmRescaler.cxx:66 An invalid logic behavior occurred
        # Anchoring at zero also keeps zero meaning zero, which matters for the
        # masked background of a brain map. The cost is at most one bit of
        # precision (1 part in 32767 of the largest magnitude), far below the
        # noise of any map this handles.
        _peak = max(abs(vmin), abs(vmax))
        if _peak > 0:
            rescale_slope = _peak / 32767.0
            rescale_intercept = 0.0
            img_int16 = np.clip(np.round(img_data / rescale_slope),
                                _i16.min, _i16.max).astype(np.int16)
        else:
            # All zeros: nothing to scale, and slope must stay non-zero.
            rescale_slope, rescale_intercept = 1.0, 0.0
            img_int16 = np.zeros(img_data.shape, dtype=np.int16)
        print(f'Quantitative mode: [{vmin:.6g}, {vmax:.6g}] -> int16, '
              f'slope={rescale_slope:.6g} intercept={rescale_intercept:.6g}')
    else:
        # Historical behaviour: fill the int16 range from 0..max. Kept for
        # backward compatibility, but now rounded and clipped rather than
        # truncated and wrapped (numpy wraps modulo 2**16 on overflow, so
        # strongly negative voxels used to come back as large positives), and
        # the scale factor is recorded instead of discarded.
        print('Converting the nifti to 16bit')
        if vmax == 0:
            print('Error: input image is all zeros; cannot rescale to int16')
            sys.exit(1)
        _scale = _i16.max / vmax
        img_int16 = np.clip(np.round(img_data * _scale),
                            _i16.min, _i16.max).astype(np.int16)
        rescale_slope, rescale_intercept = 1.0 / _scale, 0.0

    new_img = sitk.GetImageFromArray(img_int16)
    new_img.CopyInformation(nii_img)
    new_img = sitk.DICOMOrient(new_img, "LPS")

# Apply correct pixel spacing.
# Only the in-plane spacing is overridden: -p describes a rendered screenshot's
# pixel size and says nothing about slice separation, so overwriting the third
# component (as this used to) corrupted every slice position downstream.
if args.pixelspacing is not None:
    ps = args.pixelspacing
    _sp = new_img.GetSpacing()
    new_img.SetSpacing([ps, ps, _sp[2] if len(_sp) > 2 else 1.0])
    print(f'Pixel spacing set to {ps:.4f} mm (slice spacing kept at {_sp[2] if len(_sp) > 2 else 1.0})')

writer = sitk.ImageFileWriter()
writer.KeepOriginalImageUIDOn()

modification_time = time.strftime("%H%M%S")
modification_date = time.strftime("%Y%m%d")

# Image Orientation: use underlay geometry if available, else fall back to new_img direction
if underlay_geom is not None:
    _, row_dir, col_dir, _, _ = compute_slice_geometry(underlay_geom, args.plane, 0)
    orientation_str = "\\".join(_ds(v) for v in row_dir + col_dir)
else:
    d = new_img.GetDirection()
    orientation_str = "\\".join(_ds(v) for v in (d[0], d[3], d[6], d[1], d[4], d[7]))

series_uid = ("1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time
              + "." + str(int(hashlib.md5(seriesdesc.encode()).hexdigest()[:8], 16)))

# Study Instance UID: keep the output inside the donor's study. This was never
# set on this path -- neither copied nor generated -- so GDCM invented a fresh
# one and every series landed as an orphan study for the same patient, with
# cross-series cursor sync (which needs a shared frame of reference within a
# study) unable to engage. Both key spellings are probed, as elsewhere.
study_uid = None
for _k in ("0020|000d", "0020|000D"):
    if reader.HasMetaDataKey(_k):
        study_uid = _dcm_str(reader.GetMetaData(_k), _charset).strip()
        break
if not study_uid:
    study_uid = ("1.2.826.0.1.3680043.2.1125." + modification_date + ".2" + modification_time)
    print('Warning: donor has no StudyInstanceUID; generating one')

# A measurable series is an MR image, not a screenshot: Secondary Capture has no
# Modality LUT module, so PACS would ignore RescaleSlope/Intercept and report
# raw stored values in an ROI.
_is_quant = args.quantitative and input_type == 'nifti'
_sop_class = _MR_SOP_CLASS if _is_quant else _SC_SOP_CLASS

series_tag_values_a = [
    (k, _dcm_str(reader.GetMetaData(k), _charset))
    for k in tags_to_copy
    if reader.HasMetaDataKey(k)
]
series_tag_values_b = [
    ("0002|0002", _sop_class),
    ("0008|0016", _sop_class),
    ("0008|0031", modification_time),
    ("0008|0021", modification_date),
    ("0008|0008", "DERIVED\\SECONDARY"),
    ("0020|000d", study_uid),
    ("0020|000e", series_uid),
    ("0020|0037", orientation_str),
    ("0008|103e", seriesdesc),
    ("0020|0011", seriesnumber),
]
if _is_quant:
    # NOTE: the Modality LUT tags are deliberately NOT added here.
    #
    # SimpleITK writes through GDCM, and GDCM reads RescaleSlope/RescaleIntercept
    # from the dictionary as an instruction to *inverse-rescale* the pixel data
    # on the way out -- it assumes the array it was handed holds real-world
    # values and divides them down to stored values. Our array already holds
    # stored values, so that would scale them a second time. It also refuses any
    # non-integer slope or intercept outright:
    #   gdcmRescaler.cxx: An invalid logic behavior occurred slope == (int)slope
    # and an integer-only slope would quantise every fractional map to whole
    # numbers, defeating the entire point of a quantitative series.
    #
    # So the pixel data is written first, without these tags, and they are
    # attached afterwards by _attach_modality_lut(). See that function.
    pass
else:
    # Conversion Type is Type 1 (required) for Secondary Capture; it was missing.
    series_tag_values_b += [("0008|0064", "WSD")]       # Workstation

series_tag_values = series_tag_values_a + series_tag_values_b

print('Incorporating the following dicom tags:')
print(series_tag_values)

# Clean and make the output dir
if os.path.exists(dcm_output):
    shutil.rmtree(dcm_output)
os.makedirs(dcm_output, exist_ok=True)

# Write slices.
# For a quantitative series the geometry must come from the volume itself, never
# from compute_slice_geometry() -- that reproduces mrview's screenshot display
# convention, which is correct for rendered PNGs and wrong for voxel data.
_geom = None if _is_quant else underlay_geom
_plane = None if _is_quant else args.plane
if _is_quant and underlay_geom is not None:
    print('Note: -u/-o ignored in quantitative mode; geometry comes from the input volume')

list(map(
    lambda i: writeSlices(series_tag_values, new_img, dcm_output, i,
                          underlay_geom=_geom, plane=_plane,
                          series_uid=series_uid),
    range(new_img.GetDepth())
))

if _is_quant:
    # Default window from robust percentiles of the FOREGROUND stored values.
    #
    # Background is excluded on purpose. These maps are brain-masked, so ~87% of
    # every volume is exactly zero; taking percentiles over the whole array puts
    # p2 deep inside that zero mass and drags p98 down to just above it. On real
    # DSC data that produced a window of [0, 533] for rCBV whose tissue actually
    # runs to 1272, so a large part of the brain saturated white -- the series
    # opened far too dark with too much of it blown out.
    # In REAL units, not stored ones. Window Center/Width are applied after the
    # Modality LUT (DICOM PS3.3 C.11.2), so they live in the rescaled value
    # space. Computing them from the int16 pixel values made the error scale
    # with the slope: K2 (slope 2e-05) got a window 1128x its own data range and
    # rendered as a flat grey plane with the negative lobe invisible, while MTT
    # got one far too narrow and blew out. Only maps whose slope happens to be
    # near 1 looked approximately right.
    _real = img_int16.astype(np.float64) * rescale_slope + rescale_intercept
    _fg = _real[img_int16 != 0]
    if _fg.size < 100:          # not a masked map (or nearly empty): use everything
        _fg = _real
    if args.label:
        # Labels are categorical: span them all rather than clipping the tails.
        _lo, _hi = float(_fg.min()), float(_fg.max())
    else:
        # Asymmetric on purpose. The top is p99.5, not p98: on a DSC map the
        # choroid plexus (genuinely very vascular) and CSF (deconvolution
        # garbage) reach ~46x the cortical median, and a p98 top left the
        # ventricles as blown-out white blobs. p99.5 cuts the saturated voxels
        # from 2.0% to 0.5% of tissue. It costs some brightness -- cortex sits
        # at ~14% of the scale rather than ~21% -- which is the right trade for
        # a series meant to be re-windowed on PACS: brightening is easy,
        # recovering detail that was clipped away is not.
        _lo, _hi = np.percentile(_fg, 2), np.percentile(_fg, 99.5)
        # Optional extra headroom above p99.5. There is no statistic that
        # separates the maps that want it from the ones that don't -- on real
        # DSC data K2 and MTT have far heavier tails than rCBV, so a tail-based
        # rule would give exactly the wrong maps more room. It is a per-map
        # display preference, so it is passed in rather than inferred.
        if args.window_headroom and args.window_headroom != 1.0 and _hi > 0:
            _hi *= args.window_headroom
    if _hi <= _lo:
        _lo, _hi = float(_real.min()), float(_real.max())
    if _hi > _lo:
        _wc, _ww = (_lo + _hi) / 2.0, _hi - _lo
    else:
        _wc = _ww = None
    if not _attach_modality_lut(dcm_output, rescale_slope, rescale_intercept,
                                args.units if args.units else 'US', _wc, _ww):
        sys.exit(1)
    print(f'Modality LUT attached: slope={rescale_slope:.6g} '
          f'intercept={rescale_intercept:.6g} type={args.units if args.units else "US"}')

print(f'Wrote {new_img.GetDepth()} DICOMs to {dcm_output}')
sys.exit(0)
