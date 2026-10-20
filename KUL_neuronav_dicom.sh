#!/bin/bash
# Bash shell script to export fMRI maps as plain DICOM image series for
# neuronavigation systems that cannot read DICOM SEG objects (e.g. ANT Neuro
# Visor2 for TMS).
#
# Requires Mrtrix3, dcm2niix, dcmdump (dcmtk), zip, KUL_nii2dcm.py
#
# Why this exists: the Karawun export writes every label as a DICOM
# Segmentation object (SOP class SegmentationStorage / SurfaceSegmentationStorage,
# modality SEG). Brainlab reads those; most other navigation systems only import
# image series, so the T1w/FLAIR load and every label folder fails.
#
# This script writes, per stat map, two ordinary MR Image Storage series:
#   <name>_tmap_thr<T>   the map with everything below T set to 0 (real values
#                        kept via RescaleSlope, so an overlay can be windowed)
#   <name>_mask_thr<T>   binary 0/1 mask of map >= T
# Both borrow StudyInstanceUID and FrameOfReferenceUID from the donor (the
# Karawun T1w DICOM), so they land in the SAME study and frame as the
# anatomical the site already imported, and need no registration.
#
# Thresholding is one-sided (map >= T, T > 0), matching the Karawun fMRI
# labels, so negative values of a two-sided contrast are dropped. To export a
# negative contrast, negate the map first (mrcalc map.nii -neg neg_map.nii).
#
# @ KUL_NIS
version="0.2"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
cwd=$(pwd)
# Source the helpers WITHOUT this script's arguments: with a non-empty $1 they
# create KUL_LOG/<script>/ in the cwd -- the patient folder -- and print a
# welcome line, neither of which this standalone exporter wants.
_kul_args=("$@")
set --
source "$kul_main_dir/KUL_main_functions.sh"
set -- "${_kul_args[@]}"

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

$(basename "$0") exports fMRI maps as DICOM image series (thresholded t-map and
binary mask per map) for neuronavigation systems without DICOM SEG support.

Run from the clinical working directory (the one holding RESULTS/ and Karawun/).

Usage:

  $(basename "$0") -p <participant> [-T <threshold>] [OPT_ARGS]

Example:

  $(basename "$0") -p JohnDoe -T 5 -t afMRI_Sham_p001unc_k50=3

Threshold per map (map >= T, T must be a number > 0), first that applies:
     1. -t <mapname>=<value>
     2. <mapdir>/<mapname>.thresh, the same file that pins the threshold for
        the PACS export and Karawun labels (so all three agree)
     3. -T <value>
     A map with none of these is an error, as is a map with no voxel >= T.

Optional arguments:

     -T:  default threshold for maps without -t or a .thresh file
     -t:  per-map threshold as <mapname>=<value> (map name without .nii/.nii.gz);
          repeat for several maps
     -i:  directory with the stat maps (NIfTI, on the anatomical grid)
          default: RESULTS/sub-<p>/PACS_input/overlays
     -d:  donor DICOM series directory (the anatomical the site imports)
          default: Karawun/sub-<p>/sub-<p>_for_elements/T1w
     -a:  anatomical series directories to copy into the export, comma-separated
          default: the donor directory plus its siblings FLAIR_reg2_T1w and
          FAT1w, if present
     -o:  output directory
          default: Neuronav/sub-<p>
     -v:  show output from mrtrix commands

Output:

     <out>/DICOM/<series>/        one folder per DICOM series
     <out>/sub-<p>_neuronav.zip   the same, zipped for transfer (only if QA passes)
     <out>/QA.txt                 readback check of every generated series

Series numbers: t-maps 99001.., masks 99101.. -- a block no scanner or the PACS
export (donor SeriesNumber*100 + ...) uses. At most 99 maps per run.

Exit status is non-zero whenever anything was not exported as asked.

USAGE

	exit 1
}

# PixelData (7FE0,0010) present? Detected by its tag bytes, as in
# KUL_karawun_import_hint (KUL_clinical_fmridti.sh): per-slice Philips PACS
# exports put a Raw Data Storage object without pixels (often MR000001.dcm) in
# every series folder, which SimpleITK cannot read and an importer may reject.
_has_pixels() { LC_ALL=C grep -q -a -F $'\xe0\x7f\x10\x00' "$1"; }

# A threshold must be a plain positive number: map >= 0 (or a negative T)
# selects the whole zero background, and a non-number makes mrcalc fail.
_valid_thr() {
	[[ "$1" =~ ^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]] && awk -v t="$1" 'BEGIN {exit !(t > 0)}'
}

_die() { echo "Error: $*" >&2; exit 1; }


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1
default_thr=""
mapdir=""
donor_dir=""
anat_dirs=""
outdir=""
declare -A map_thr=()

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1
fi

while getopts "p:T:t:i:d:a:o:v" OPT; do

	case $OPT in
	p) #participant
		participant=$OPTARG
	;;
	T) #default threshold
		default_thr=$OPTARG
	;;
	t) #per-map threshold
		if [[ "$OPTARG" != *=* ]]; then
			echo "Option -t expects <mapname>=<value> (got '$OPTARG')" >&2
			exit 1
		fi
		map_thr["${OPTARG%%=*}"]="${OPTARG#*=}"
	;;
	i) #map dir
		mapdir=$OPTARG
	;;
	d) #donor dir
		donor_dir=$OPTARG
	;;
	a) #anatomical dirs
		anat_dirs=$OPTARG
	;;
	o) #output
		outdir=$OPTARG
	;;
	v) #verbose
		silent=0
	;;
	\?)
		echo "Invalid option: -$OPTARG" >&2
		echo
		Usage >&2
		exit 1
	;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		echo
		Usage >&2
		exit 1
	;;
	esac

done

if [ -z "$participant" ]; then
	echo "Option -p is required." >&2
	exit 2
fi

[ -z "$mapdir" ]    && mapdir="RESULTS/sub-${participant}/PACS_input/overlays"
[ -z "$donor_dir" ] && donor_dir="Karawun/sub-${participant}/sub-${participant}_for_elements/T1w"
[ -z "$outdir" ]    && outdir="Neuronav/sub-${participant}"
if [ -z "$anat_dirs" ]; then
	anat_dirs="$donor_dir"
	for _sib in FLAIR_reg2_T1w FAT1w; do
		_d="$(dirname "$donor_dir")/$_sib"
		[ -d "$_d" ] && anat_dirs="$anat_dirs,$_d"
	done
fi

# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
fi


# --- MAIN ---
printf "\n\n\n"

# --- Inputs ---------------------------------------------------------------
[ -d "$mapdir" ] || _die "map directory $mapdir not found"
maps=()
for f in "$mapdir"/*.nii "$mapdir"/*.nii.gz; do
	[ -f "$f" ] && maps+=("$f")
done
[ ${#maps[@]} -gt 0 ] || _die "no NIfTI maps in $mapdir"
[ ${#maps[@]} -le 99 ] || _die "${#maps[@]} maps; at most 99 per run (series numbers 99001-99099 / 99101-99199)"

[ -d "$donor_dir" ] || _die "donor directory $donor_dir not found"
donor_dcm=""
while IFS= read -r -d '' _f; do
	if _has_pixels "$_f"; then donor_dcm="$_f"; break; fi
done < <(find "$donor_dir" -maxdepth 1 -type f -not -name ".*" -not -iname "DICOMDIR" -print0 | sort -z)
[ -n "$donor_dcm" ] || _die "no DICOM with image data in $donor_dir (non-image objects such as Raw Data Storage cannot be a donor)"

for tool in dcm2niix dcmdump mrcalc mrstats mrinfo zip; do
	command -v $tool >/dev/null 2>&1 || _die "$tool not found"
done

# Resolve and validate every map's threshold, and check it leaves voxels,
# BEFORE converting anything: a run either exports every map as asked or
# exports nothing.
declare -A thr_of=() src_of=()
for f in "${maps[@]}"; do
	n=$(basename "$f"); n=${n%.gz}; n=${n%.nii}
	_tf="$mapdir/$n.thresh"
	if [ -n "${map_thr[$n]+x}" ]; then
		thr=${map_thr[$n]}; src="-t"
	elif [ -f "$_tf" ] && [ -n "$(head -1 "$_tf" | tr -d '[:space:]')" ]; then
		thr=$(head -1 "$_tf" | tr -d '[:space:]'); src="$(basename "$_tf")"
	elif [ -n "$default_thr" ]; then
		thr=$default_thr; src="-T"
	else
		_die "no threshold for map $n (give -t $n=<value>, a $n.thresh file, or -T)"
	fi
	_valid_thr "$thr" || _die "threshold '$thr' for map $n (from $src) is not a number > 0"
	# Voxels >= T: mean x count of the 0/1 mask (mrstats has no sum).
	nvox=$(mrcalc "$f" "$thr" -ge - -datatype uint8 2>/dev/null \
		| mrstats - -output mean -output count 2>/dev/null \
		| awk 'NF==2 {printf "%d\n", $1*$2+0.5}')
	[[ "$nvox" =~ ^[0-9]+$ ]] || _die "could not threshold map $n (mrcalc/mrstats failed on $f)"
	[ "$nvox" -gt 0 ] || _die "map $n has no voxels >= $thr (from $src); lower the threshold or remove the map"
	thr_of[$n]=$thr; src_of[$n]=$src
done
for n in "${!map_thr[@]}"; do
	[ -f "$mapdir/$n.nii" ] || [ -f "$mapdir/$n.nii.gz" ] \
		|| _die "-t given for $n, but $mapdir/$n.nii(.gz) does not exist"
done

# --- Python interpreter for KUL_nii2dcm.py (same lookup as KUL_clinical_fmridti.sh)
_nii2dcm_py=""
if command -v conda >/dev/null 2>&1; then
	_conda_base=$(conda info --base 2>/dev/null)
	if [ -n "$_conda_base" ] && [ -x "$_conda_base/envs/$KUL_DICOM_ENV/bin/python" ]; then
		_nii2dcm_py="$_conda_base/envs/$KUL_DICOM_ENV/bin/python"
	fi
fi
if [ -z "$_nii2dcm_py" ]; then
	for _sw in "${SOFTWARE_ROOT:-}" /opt/kul_software /usr/local/KUL_apps; do
		[ -z "$_sw" ] && continue
		if [ -x "$_sw/miniforge3/envs/$KUL_DICOM_ENV/bin/python" ]; then
			_nii2dcm_py="$_sw/miniforge3/envs/$KUL_DICOM_ENV/bin/python"
			break
		fi
	done
fi
[ -z "$_nii2dcm_py" ] && _nii2dcm_py="python3"
"$_nii2dcm_py" -c "import SimpleITK, PIL, numpy, pydicom" >/dev/null 2>&1 \
	|| _die "'$_nii2dcm_py' cannot import SimpleITK/Pillow/numpy/pydicom (see KUL_DICOM_ENV)"
_nii2dcm_script="$kul_main_dir/KUL_nii2dcm.py"

echo "KUL_neuronav_dicom v${version}"
echo "  maps:   $mapdir (${#maps[@]})"
echo "  donor:  $donor_dcm"
echo "  output: $outdir"

dcmdir="$outdir/DICOM"
logdir="$outdir/logs"
tmp="$outdir/tmp"
rm -rf "$dcmdir" "$tmp" "$logdir" "$outdir/QA.txt" "$outdir/sub-${participant}_neuronav.zip"
mkdir -p "$dcmdir" "$logdir" "$tmp"
trap 'rm -rf "$tmp"' EXIT


# --- Anatomicals: image objects copied as-is, they already import fine ------
# Only files that hold pixels: a scanner series folder also carries a Raw Data
# Storage object, which must not reach the navigation system's importer.
IFS=',' read -ra _anats <<< "$anat_dirs"
for a in "${_anats[@]}"; do
	[ -d "$a" ] || _die "anatomical series directory $a not found"
	_dst="$dcmdir/$(basename "$a")"
	mkdir -p "$_dst"
	_ncp=0; _nskip=0
	while IFS= read -r -d '' _f; do
		if _has_pixels "$_f"; then
			cp -p "$_f" "$_dst/" || _die "copying $_f failed"
			_ncp=$((_ncp+1))
		else
			_nskip=$((_nskip+1))
		fi
	done < <(find "$a" -maxdepth 1 -type f -not -name ".*" -not -iname "DICOMDIR" -print0)
	[ $_ncp -gt 0 ] || _die "no image DICOMs in anatomical directory $a"
	echo "  copied anatomical: $(basename "$a") ($_ncp files$([ $_nskip -gt 0 ] && echo ", skipped $_nskip without image data"))"
done


# --- Maps -> thresholded t-map + binary mask -------------------------------
series_list=()   # "<series dir>|<source nifti>|<kind>"
idx=0
for f in "${maps[@]}"; do
	idx=$((idx+1))
	n=$(basename "$f"); n=${n%.gz}; n=${n%.nii}
	thr=${thr_of[$n]}

	tmap_name="${n}_tmap_thr${thr}"
	mask_name="${n}_mask_thr${thr}"
	# Temporary files are named by index: map names may hold spaces or other
	# characters that later tools (dcm2niix -f) rewrite.
	t_mask="$tmp/map${idx}_mask.nii.gz"
	t_tmap="$tmp/map${idx}_tmap.nii.gz"

	mrcalc "$f" "$thr" -ge "$t_mask" -datatype uint8 -force \
		|| _die "mrcalc failed making the mask for $n"
	mrcalc "$f" "$t_mask" -mult "$t_tmap" -datatype float32 -force \
		|| _die "mrcalc failed making the t-map for $n"
	echo "  $n: threshold $thr (${src_of[$n]})"

	"$_nii2dcm_py" "$_nii2dcm_script" -q -s "$tmap_name" -n $((99000+idx)) \
		"$t_tmap" "$donor_dcm" "$dcmdir/$tmap_name" >"$logdir/map${idx}_tmap.log" 2>&1 \
		|| _die "KUL_nii2dcm.py failed for $tmap_name -- see $logdir/map${idx}_tmap.log"
	"$_nii2dcm_py" "$_nii2dcm_script" --label -s "$mask_name" -n $((99100+idx)) \
		"$t_mask" "$donor_dcm" "$dcmdir/$mask_name" >"$logdir/map${idx}_mask.log" 2>&1 \
		|| _die "KUL_nii2dcm.py failed for $mask_name -- see $logdir/map${idx}_mask.log"
	series_list+=("$dcmdir/$tmap_name|$t_tmap|tmap")
	series_list+=("$dcmdir/$mask_name|$t_mask|mask")
done

# Belt and braces: every map was checked above, so this cannot be short.
[ ${#series_list[@]} -eq $((2 * ${#maps[@]})) ] \
	|| _die "wrote ${#series_list[@]} map series for ${#maps[@]} maps (expected $((2 * ${#maps[@]})))"


# --- QA: read every generated series back and check it -------------------
# 1. same Study and Frame of Reference UID as the donor
# 2. same voxel grid (dims + transform) as the donor series, read back the same way
# 3. same values as the NIfTI it came from (mask exact; t-map within one
#    quantisation step, i.e. the RescaleSlope)
# Reading both back through dcm2niix means a geometry or left/right error in
# the conversion shows up as a transform mismatch against the anatomical the
# site sees, rather than at the navigation system.
qa="$outdir/QA.txt"
fail=0
# Top-level tags only (anchored at column 0): the Karawun T1w also carries the
# ORIGINAL scanner StudyInstanceUID nested in reference sequences. -Un prints
# UIDs as numbers rather than names (=MRImageStorage).
_tag() { dcmdump -Un "$2" 2>/dev/null | grep -m1 "^($1)" | sed 's/.*\[\(.*\)\].*/\1/'; }
_first() { find "$1" -maxdepth 1 -type f | sort | head -1; }
_readback() {  # $1 dicom dir, $2 output basename (safe chars only) -> echoes nifti path
	dcm2niix -z y -b n -f "$2" -o "$tmp" "$1" >/dev/null 2>&1
	ls "$tmp/$2"*.nii.gz 2>/dev/null | head -1
}
_same_grid() {  # $1 $2 niftis -> 0 if same size and transform within 1e-3
	[ "$(mrinfo -size "$1")" = "$(mrinfo -size "$2")" ] || return 1
	paste <(mrinfo -transform "$1" | tr -s ' ' '\n' | grep .) \
	      <(mrinfo -transform "$2" | tr -s ' ' '\n' | grep .) |
		awk '{d=$1-$2; if (d<0) d=-d; if (d>1e-3) bad=1} END {exit bad}'
}

d_study=$(_tag 0020,000d "$donor_dcm")
d_for=$(_tag 0020,0052 "$donor_dcm")
# Read back the donor from its (image-only) copy in the export.
donor_nii=$(_readback "$dcmdir/$(basename "$donor_dir")" "qa_donor")

{
echo "KUL_neuronav_dicom v${version} QA  $(date)"
echo "donor: $donor_dcm"
echo "  StudyInstanceUID    $d_study"
echo "  FrameOfReferenceUID $d_for"
echo
} > "$qa"

if [ -z "$d_study" ] || [ -z "$d_for" ]; then
	echo "FAIL  donor has no StudyInstanceUID/FrameOfReferenceUID" | tee -a "$qa"
	fail=1
fi
if [ -z "$donor_nii" ]; then
	echo "FAIL  could not read back the donor series with dcm2niix" | tee -a "$qa"
	fail=1
fi

qi=0
for entry in "${series_list[@]}"; do
	qi=$((qi+1))
	IFS='|' read -r sdir src kind <<< "$entry"
	sname=$(basename "$sdir")
	msg=""
	f1=$(_first "$sdir")

	[ "$(_tag 0008,0016 "$f1")" = "1.2.840.10008.5.1.4.1.1.4" ] || msg="$msg not-MRImageStorage"
	[ "$(_tag 0020,000d "$f1")" = "$d_study" ] || msg="$msg study-UID"
	[ "$(_tag 0020,0052 "$f1")" = "$d_for" ]   || msg="$msg FoR-UID"

	rb=$(_readback "$sdir" "qa_series${qi}")
	if [ -z "$rb" ]; then
		msg="$msg readback-failed"
	else
		[ -n "$donor_nii" ] && { _same_grid "$rb" "$donor_nii" || msg="$msg grid-vs-donor"; }
		if _same_grid "$rb" "$src"; then
			maxdiff=$(mrcalc "$rb" "$src" -sub -abs - 2>/dev/null | mrstats - -output max)
			if [ "$kind" = "mask" ]; then
				tol=0
			else
				slope=$(_tag 0028,1053 "$f1")
				tol=$(awk -v s="$slope" 'BEGIN {printf "%.10f", s*0.51}')
			fi
			if [ -z "$maxdiff" ]; then
				msg="$msg values(not-computed)"
			else
				awk -v d="$maxdiff" -v t="$tol" 'BEGIN {exit !(d<=t)}' || msg="$msg values(maxdiff=$maxdiff)"
			fi
		else
			msg="$msg grid-vs-source"
		fi
	fi

	if [ -z "$msg" ]; then
		echo "PASS  $sname" | tee -a "$qa"
	else
		echo "FAIL  $sname:$msg" | tee -a "$qa"
		fail=1
	fi
done

if [ $fail -ne 0 ]; then
	echo ""
	echo "QA FAILED -- do not send this export. See $qa"
	exit 1
fi


# --- Package --------------------------------------------------------------
zipf="$(cd "$outdir" && pwd)/sub-${participant}_neuronav.zip"
rm -f "$zipf"
(cd "$outdir" && zip -qr "$zipf" DICOM) || _die "zip failed"
echo ""
echo "Done. ${#series_list[@]} map series (${#maps[@]} maps) + ${#_anats[@]} anatomical series, QA passed."
echo "  $dcmdir"
echo "  $zipf"
