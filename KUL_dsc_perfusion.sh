#!/bin/bash
# Bash shell script to process DSC (dynamic susceptibility contrast) perfusion MRI
#  of tumour patients, and to normalise the resulting maps against contralesional NAWM
#
# Requires ANTs, MRtrix3, FreeSurfer (mri_synthstrip, mri_vol2vol) and the
#  conda env named by $KUL_PYFMRI_ENV (default 'pyfMRI')
#
# @ Ahmed Radwan - KU Leuven, Translational MRI - ahmed.radwan@kuleuven.be
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 08/08/2026
version="2.0"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` computes DSC perfusion maps (rCBV, rCBF, MTT, Tmax, TTP, TT0,
 K1, K2, CNR) and normalises them against contralesional normal-appearing white
 matter.

Method (v2): MP-PCA denoising and rigid motion correction on the native DSC
 grid; the fit runs on that grid (no 4D resampling); BSW leakage correction
 with a non-enhancing reference tissue curve and arrival-time correction;
 cluster-selected AIF; oSVD deconvolution. The 3D maps are then moved to T1w
 space with ONE interpolation. See docs/KUL_dsc_perfusion/KUL_dsc_perfusion.md.

It expects the 4D DSC series in BIDS (BIDS/sub-{participant}/perf/*_dsc.nii.gz,
 as written by KUL_dcm2bids.sh with a 'DSC' entry in study_config/sequences.txt),
 but any 4D NIfTI can be given with -d.

All maps are delivered in the participant's T1w space, i.e. the same space as
 everything else in RESULTS/sub-{participant}/, so they can be sent to PACS
 and Karawun alongside the fMRI and tractography results.

Normalisation (the lesion/NAWM ratios that are actually reported) additionally
 needs a lesion mask and a FreeSurfer aseg. Both are produced by
 KUL_clinical_fmridti.sh; without them the parametric maps are still written
 and the ROI report holds absolute values only, so it can be re-run later.

ROI report (always written): per-ROI measures for the HD-GLIO enhancing and
 non-enhancing compartments, the resection cavity and the whole lesion (when
 RESULTS/sub-{participant}/Lesion/ holds them), any mask or label map dropped
 in RESULTS/sub-{participant}/Perfusion/ROIs/, and the contralesional NAWM and
 cortex references -> Perfusion/*_perfusion_report.tsv, *_perfusion_roi_stats.tsv
 and REPORT/sub-{participant}_06_DSC_report.pdf. Re-running the script after
 adding ROIs recomputes only the report (the fit is reused).

Usage:

  `basename $0` <OPT_ARGS>

Examples:

  # inside a KUL_NIS study dir, after KUL_clinical_fmridti.sh has run
  `basename $0` -p JaneDoe -n 32

  # standalone, giving the DSC series and the anatomical reference explicitly
  `basename $0` -p JaneDoe -d /data/raw/dsc.nii.gz -a /data/raw/T1w.nii.gz \\
      -e 0.030 -r 1.5 -n 32

Required arguments:

     -p:  participant name

Optional arguments:

     -d:  4D DSC series (default: BIDS/sub-{participant}/perf/*_dsc.nii.gz)
     -a:  anatomical reference to register the DSC to, and the output space
            (default: RESULTS/sub-{participant}/Anat/cT1w_reg2_T1w.nii.gz,
             falling back to .../Anat/T1w.nii.gz)
     -l:  lesion mask in the -a space, for the NAWM normalisation
            (default: RESULTS/sub-{participant}/Lesion/sub-{participant}_lesion_and_cavity.nii.gz,
             falling back to .../Lesion/lesion.nii.gz)
     -F:  FreeSurfer subject dir holding mri/aseg.mgz, for the hemispheric WM
            (default: the KUL_VBG output, falling back to
             BIDS/derivatives/freesurfer/sub-{participant})
     -I:  directory with vendor/console-exported perfusion maps (nrCBVcorr,
            nrCBVuncorr, nrelCBFAIF, ...) to co-register for side-by-side QC.
            Optional; the in-house fit is what gets reported either way.
     -e:  echo time in seconds (default: EchoTime from the BIDS sidecar)
     -r:  repetition time in seconds (default: RepetitionTime from the sidecar,
            falling back to the 4th dimension spacing of the DSC series)
     -b:  pre-contrast baseline frames as START:END, or 'auto' (default)
     -w:  rCBV integration window: full (default; bolus arrival to the end of
            the series), firstpass (to the end of the first pass), or START:END
     -s:  FWHM in mm of the light 3D smoothing used for the EXTRA *_smooth
            rCBV/rCBF maps (default 3; 0 = none). The unsmoothed maps are always
            written and are the ones used for the stats.
     -D:  do not denoise (MP-PCA) the series before fitting
     -P:  phase-encoding axis of the DSC EPI, used to restrict the distortion
            correction: i, j, k, or none for unrestricted SyN
            (default: the axis of PhaseEncodingDirection in the BIDS sidecar,
             falling back to j)
     -A, -L:  accepted for backward compatibility and ignored (the v1 AIF
            ranking and legacy deconvolution are gone; the v1 code is in git
            history before the method-v2 commit)
     -y:  conda env to use instead of \$KUL_PYFMRI_ENV (default 'pyfMRI')
     -R:  redo, i.e. throw away previous results and recompute everything
     -n:  number of threads to use (default 32)
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ants_verbose=1
ncpu=32
verbose_level=1
redo=0
window="full"
smooth_fwhm=3
do_denoise=1
pe_axis="j"
pe_axis_given=0
baseline="auto"
in_dsc=""
anat_ref=""
lesion_mask=""
fs_dir=""
isp_dir=""
te_arg=""
tr_arg=""
pyfmri_env="$KUL_PYFMRI_ENV"

# Set required options
p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:d:a:l:F:I:e:r:b:w:s:P:y:n:v:A:DLR" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
			p_flag=1
		;;
		d) #DSC series
			in_dsc=$OPTARG
		;;
		a) #anatomical reference
			anat_ref=$OPTARG
		;;
		l) #lesion mask
			lesion_mask=$OPTARG
		;;
		F) #freesurfer subject dir
			fs_dir=$OPTARG
		;;
		I) #vendor/ISP maps
			isp_dir=$OPTARG
		;;
		e) #echo time
			te_arg=$OPTARG
		;;
		r) #repetition time
			tr_arg=$OPTARG
		;;
		b) #baseline frames
			baseline=$OPTARG
		;;
		P) #phase encoding axis
			pe_axis=$OPTARG
			pe_axis_given=1
		;;
		y) #conda env
			pyfmri_env=$OPTARG
		;;
		w) #integration window
			window=$OPTARG
		;;
		s) #smoothing FWHM for the extra *_smooth maps
			smooth_fwhm=$OPTARG
		;;
		D) #no denoising
			do_denoise=0
		;;
		A|L) #v1 options, kept so older callers do not break
			echo "NOTE: -$OPT is a v1 option and is ignored (method v2: cluster AIF, oSVD)." >&2
		;;
		R) #redo
			redo=1
		;;
		n) #ncpu
			ncpu=$OPTARG
		;;
		v) #verbose
			verbose_level=$OPTARG
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

fi

# check for required options
if [ $p_flag -eq 0 ] ; then
	echo
	echo "Option -p is required: give the BIDS name of the participant." >&2
	echo
	exit 2
fi
if [[ ! "$window" =~ ^(full|firstpass|[0-9]+:[0-9]+)$ ]]; then
	echo "Option -w must be full, firstpass or START:END (got '$window')." >&2
	exit 2
fi
if ! awk -v x="$smooth_fwhm" 'BEGIN { exit !(x ~ /^[0-9]*\.?[0-9]+$/) }'; then
	echo "Option -s must be a non-negative number of mm (got '$smooth_fwhm')." >&2
	exit 2
fi

KUL_LOG_DIR="KUL_LOG/${script}/sub-${participant}"
mkdir -p $KUL_LOG_DIR

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
	export MRTRIX_QUIET=1
	silent=1
	str_silent=" > /dev/null 2>&1"
	ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
	silent=0
	str_silent=""
	ants_verbose=1
fi

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$ncpu


# --- directories ---
globalresultsdir=$cwd/RESULTS/sub-$participant
derivativesdir=$cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_dsc_perfusion
resultsdir=$globalresultsdir/Perfusion

if [ $redo -eq 1 ]; then
	kul_echo "Option -R given: removing previous DSC results"
	rm -rf $derivativesdir $resultsdir
fi

workdir=$derivativesdir/work
fitdir=$derivativesdir/fit
mkdir -p $workdir $fitdir $resultsdir


# --- functions ---

# Float helpers. These deliberately do not use bc: mrstats reports small values
# in scientific notation (e.g. 4.59e-05) and bc cannot parse that at all -- it
# writes "syntax error" to stderr and prints nothing, so a test like
# [ "$(echo "$x > 0" | bc -l)" -eq 1 ] silently takes the false branch on
# exactly the values worth checking. awk handles the notation natively.
function KUL_float_gt {
	awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 > b + 0) }'
}

function KUL_float_div {
	awk -v a="$1" -v b="$2" 'BEGIN { if (b + 0 == 0) print "n/a"; else printf "%.4f\n", a / b }'
}

function KUL_float_pct {
	awk -v a="$1" -v b="$2" 'BEGIN { if (b + 0 == 0) print "0.00"; else printf "%.2f\n", 100 * a / b }'
}

# read a string field out of a BIDS json sidecar; prints nothing when absent
function KUL_json_string {
	local json="$1"
	local field="$2"
	[ -f "$json" ] || return 0
	python3 -c "
import json, sys
try:
    v = json.load(open('$json')).get('$field')
except Exception:
    sys.exit(0)
if isinstance(v, str):
    print(v)
" 2>/dev/null
}

# restrict the SyN deformation to the phase-encoding axis; an unrestricted SyN
# is free to warp anatomy that the EPI distortion never touched, which on a
# perfusion map shows up as displaced rCBV rather than as an obvious mis-registration
function KUL_dsc_set_pe_restriction {
	case $pe_axis in
		i) restrict_deformation="1x0x0" ;;
		j) restrict_deformation="0x1x0" ;;
		k) restrict_deformation="0x0x1" ;;
		none) restrict_deformation="1x1x1" ;;
		*)
			echo "Invalid phase-encoding axis '$pe_axis': expected i, j, k or none." >&2
			exit 2
		;;
	esac
	kul_echo "  distortion correction restricted to the ${pe_axis} axis ($restrict_deformation)"
}

# read a numeric field out of a BIDS json sidecar; prints nothing when absent
function KUL_json_field {
	local json="$1"
	local field="$2"
	[ -f "$json" ] || return 0
	python3 -c "
import json, sys
try:
    v = json.load(open('$json')).get('$field')
except Exception:
    sys.exit(0)
if isinstance(v, (int, float)):
    print(v)
" 2>/dev/null
}

# resolve the DSC series, the anatomical reference and the acquisition timing
function KUL_dsc_resolve_inputs {

	if [ -z "$in_dsc" ]; then
		local found=($(find $cwd/BIDS/sub-${participant} -name "*_dsc.nii.gz" -type f 2>/dev/null | sort))
		if [ ${#found[@]} -eq 0 ]; then
			echo "ERROR: no DSC series found in BIDS/sub-${participant} (looked for *_dsc.nii.gz)." >&2
			echo "  Add a 'DSC,<search-string>' line to study_config/sequences.txt and re-run" >&2
			echo "  KUL_dcm2bids.sh, or point at the file directly with -d." >&2
			exit 1
		fi
		in_dsc=${found[0]}
		if [ ${#found[@]} -gt 1 ]; then
			kul_echo "WARNING: ${#found[@]} DSC series found; using $(basename $in_dsc). Use -d to pick another."
		fi
	fi
	if [ ! -f "$in_dsc" ]; then
		echo "ERROR: DSC series '$in_dsc' does not exist." >&2
		exit 1
	fi

	# the DSC has to be 4D, otherwise everything below is meaningless
	local n_vols=$(mrinfo -size "$in_dsc" | awk '{print $4}')
	if [ -z "$n_vols" ] || [ "$n_vols" -lt 10 ]; then
		echo "ERROR: '$in_dsc' has ${n_vols:-no} volumes; a DSC time series is required." >&2
		exit 1
	fi
	kul_echo "  DSC series: $in_dsc (${n_vols} volumes)"

	if [ -z "$anat_ref" ]; then
		if [ -f "$globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz" ]; then
			anat_ref="$globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz"
		elif [ -f "$globalresultsdir/Anat/T1w.nii.gz" ]; then
			anat_ref="$globalresultsdir/Anat/T1w.nii.gz"
		else
			echo "ERROR: no anatomical reference found in $globalresultsdir/Anat." >&2
			echo "  Run KUL_clinical_fmridti.sh first, or give one with -a." >&2
			exit 1
		fi
	fi
	if [ ! -f "$anat_ref" ]; then
		echo "ERROR: anatomical reference '$anat_ref' does not exist." >&2
		exit 1
	fi
	kul_echo "  anatomical reference / output space: $anat_ref"

	# TE and TR: the sidecar is authoritative. The prototype scraped TE out of
	# the mrinfo comments field and glued "0.0" in front of it, which quietly
	# produced a wrong TE on any series whose comment format differed.
	local json="${in_dsc%.nii.gz}.json"
	if [ -z "$te_arg" ]; then
		te_arg=$(KUL_json_field "$json" EchoTime)
	fi
	if [ -z "$te_arg" ]; then
		echo "ERROR: could not read EchoTime from $json." >&2
		echo "  Give the echo time in seconds with -e (e.g. -e 0.030)." >&2
		exit 1
	fi
	if [ -z "$tr_arg" ]; then
		tr_arg=$(KUL_json_field "$json" RepetitionTime)
	fi
	if [ -z "$tr_arg" ]; then
		tr_arg=$(mrinfo -spacing "$in_dsc" | awk '{print $4}')
		kul_echo "  no RepetitionTime in the sidecar; using the 4th-dimension spacing"
	fi
	if [ -z "$tr_arg" ] || ! KUL_float_gt "$tr_arg" 0; then
		echo "ERROR: could not determine a valid TR; give it with -r." >&2
		exit 1
	fi

	# a TE in ms rather than s is the classic sidecar mistake, and it scales
	# dR2* by a factor of 1000 without changing anything visible on the maps
	if KUL_float_gt "$te_arg" 1; then
		echo "ERROR: TE = $te_arg looks like milliseconds; -e expects seconds (e.g. 0.030)." >&2
		exit 1
	fi
	kul_echo "  TE = ${te_arg} s, TR = ${tr_arg} s"

	# the lesion mask is needed before the fit now: it is kept out of the BSW
	# reference tissue and the AIF search (an enhancing tumour in the reference
	# would make the correction partly cancel itself)
	if [ -z "$lesion_mask" ]; then
		if [ -f "$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" ]; then
			lesion_mask="$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz"
		elif [ -f "$globalresultsdir/Lesion/lesion.nii.gz" ]; then
			lesion_mask="$globalresultsdir/Lesion/lesion.nii.gz"
		fi
	fi

	# PhaseEncodingDirection is written by KUL_dcm2bids.sh from the pe_dir
	# column of study_config/sequences.txt. Only the axis matters here -- the
	# sign says which way the distortion goes, not which axis it acts along.
	if [ $pe_axis_given -eq 0 ]; then
		local pe_dir=$(KUL_json_string "$json" PhaseEncodingDirection)
		if [ -n "$pe_dir" ]; then
			pe_axis="${pe_dir%-}"
			kul_echo "  phase-encoding direction from the sidecar: $pe_dir"
		else
			kul_echo "  no PhaseEncodingDirection in the sidecar; assuming '$pe_axis' (override with -P)"
		fi
	fi
}


# brain-extract the anatomical reference and regrid it to the DSC voxel size.
# This grid is only the FIXED image of the DSC-to-anatomy registration (the
# registration parameters below were validated on it). The fit itself runs on
# the native DSC grid since method v2: resampling the 4D series onto a grid with
# the anatomy's axes cut obliquely through the 4 mm slices (+ gap), blurred and
# displaced structure, and gave the "blocky" maps.
function KUL_dsc_prepare_anat {

	if [ ! -f "$anat_regrid" ]; then

		task_in="KUL_dsc_strip_anat"
		KUL_task_exec $verbose_level "DSC - brain extracting the anatomical reference" "01_anat_brain" || return 1

		local spacing=($(mrinfo -spacing "$in_dsc"))
		task_in="mrgrid -force -voxel ${spacing[0]},${spacing[1]},${spacing[2]} \
			$workdir/anat_brain.nii.gz regrid $anat_regrid; \
			mrcalc -force $anat_regrid 0 -gt $anat_regrid_mask"
		KUL_task_exec $verbose_level "DSC - regridding the anatomical reference" "02_anat_regrid" || return 1

	else
		kul_echo "  the anatomical reference is already prepared"
	fi
}


# mri_synthstrip wrapper. The prototype always passed -g, which aborts outright
# ("CUDA is not available") on a CPU-only node and, less obviously, also on GPU
# nodes whose FreeSurfer ships a CPU-only torch. So try the GPU once and
# remember the answer, rather than assuming either way.
synthstrip_gpu=-1  # -1 = not probed yet, 1 = use -g, 0 = CPU
function KUL_dsc_synthstrip {
	local in_img="$1"
	local out_img="$2"
	local out_mask="$3"

	if [ $synthstrip_gpu -ne 0 ]; then
		if mri_synthstrip -i "$in_img" -o "$out_img" -m "$out_mask" -g; then
			synthstrip_gpu=1
			return 0
		fi
		if [ $synthstrip_gpu -eq 1 ]; then
			return 1   # the GPU worked before, so this is a real failure
		fi
		synthstrip_gpu=0
		echo "mri_synthstrip: no usable GPU, falling back to CPU" >&2
	fi

	mri_synthstrip -i "$in_img" -o "$out_img" -m "$out_mask"
}

function KUL_dsc_strip_anat {
	KUL_dsc_synthstrip "$anat_ref" "$workdir/anat_brain.nii.gz" "$workdir/anat_brain_mask.nii.gz"
}

function KUL_dsc_strip_dscvol1 {
	mrconvert -force -coord 3 0 -axes 0,1,2 "$dsc_mc" "$dsc_vol1" && \
	KUL_dsc_synthstrip "$dsc_vol1" "$dsc_vol1_brain" "$dsc_vol1_mask"
}

# brain mask ON THE NATIVE DSC GRID: synthstrip of the temporal median,
# intersected with the anatomical brain mask pulled back through the inverse
# DSC-to-anat transform (removes extra-cerebral voxels that would otherwise
# enter the AIF pool). Written with the strides of the DSC volume so the fit
# can check that mask and series share one grid.
function KUL_dsc_make_brainmask {
	mrmath -force -axis 3 "$dsc_mc" median "$workdir/dsc_mc_median.nii.gz" && \
	KUL_dsc_synthstrip "$workdir/dsc_mc_median.nii.gz" \
		"$workdir/dsc_mc_median_brain.nii.gz" "$workdir/dsc_mc_median_mask.nii.gz" && \
	antsApplyTransforms -d 3 --float 1 -i "$workdir/anat_brain_mask.nii.gz" -r "$dsc_vol1" \
		-t "[${epic_prefix}0GenericAffine.mat,1]" -t "${epic_prefix}1InverseWarp.nii.gz" \
		-n Linear -o "$workdir/anat_mask_in_dsc.nii.gz" && \
	mrcalc -quiet "$workdir/anat_mask_in_dsc.nii.gz" 0.5 -gt "$workdir/dsc_mc_median_mask.nii.gz" -mult - | \
		mrconvert -quiet -force - -strides "$dsc_vol1" -datatype uint8 "$dsc_mask"
}

# the lesion (and cavity) on the native DSC grid, dilated by one voxel: kept out
# of the BSW reference tissue and of the AIF search
function KUL_dsc_make_exclusion {
	antsApplyTransforms -d 3 --float 1 -i "$lesion_mask" -r "$dsc_vol1" \
		-t "[${epic_prefix}0GenericAffine.mat,1]" -t "${epic_prefix}1InverseWarp.nii.gz" \
		-n Linear -o "$workdir/lesion_in_dsc.nii.gz" && \
	mrcalc -quiet "$workdir/lesion_in_dsc.nii.gz" 0.05 -gt - | \
		maskfilter -quiet - dilate -npass 1 - | \
		mrconvert -quiet -force - -strides "$dsc_vol1" -datatype uint8 "$dsc_excl"
}


# The ANTs calls below live in their own functions purely so that KUL_task_exec
# receives a bare function name. It expands "$task_in" unquoted, which would let
# the shell treat bracketed ANTs arguments -- MI[...], [out,warped], [0.1,3,0] --
# as glob character classes.
# Marchenko-Pastur PCA denoising (dwidenoise), not ANTs DenoiseImage -d 4.
#
# The command is named for diffusion data but nothing in the algorithm is
# diffusion-specific: it takes a plain 4D image (it has no -grad/-fslgrad
# options at all), reshapes each spatial patch into a voxel-by-volume matrix and
# fits the MP distribution to that matrix's eigenspectrum. An 80-frame DSC
# series satisfies its assumptions at least as well as a dMRI series does --
# baseline, bolus passage and recirculation are a handful of temporal
# components, and every frame shares TE/TR/coil so the noise is i.i.d. across
# volumes. Its patches are purely SPATIAL, with the volumes as matrix columns,
# so there is no temporal neighbourhood, no edge effect, and hence none of the
# mirror-padding the previous implementation needed.
#
# Measured on a validation subject (80 frames, 144x144x42), against DenoiseImage -d 4:
#
#                       ANTs NLM      MP-PCA      raw
#   runtime             23.1 min      21.5 s      --
#   proxy rCBV vs raw   -11.28%       -0.88%      --
#   rCBV corr with raw    0.9685       0.9970     1.0
#   bolus depth          32.00%       34.48%     34.48%
#   spatial sharpness     0.15344      0.17826    0.18003
#
# The NLM bias is spatial, not temporal: bolus FWHM (7 frames) and per-voxel TTP
# (median 23, IQR 1) are identical under both, so TTP/MTT were never at risk,
# but NLM blurs ~15% and that mixes deep-dipping voxels with shallow neighbours,
# costing ~11% of rCBV. MP-PCA removes less noise (baseline tSNR 42.6 -> 53.9,
# versus 140.2 for NLM), but most of NLM's advantage there is the blurring
# itself. If the maps ever need to be smoother, smooth them explicitly after
# the fit rather than buying it here at the cost of the amplitude.
#
# Patch shape: dwidenoise's default is the smallest cube with more voxels than
# volumes (5x5x5 for 60 frames). With thick slices (here 4 mm + 1.2 mm gap) that
# cube spans 26 mm through-plane, pooling tissue that has nothing in common. When
# the slice spacing is over 1.5x the in-plane spacing the patch is therefore
# e x e x 3, with e the smallest odd size giving >= 1.2x as many voxels as frames.
# On the validation subject MP-PCA lowered split-half rCBV noise by ~20% without
# changing lesion/NAWM; it remains optional (-D) since it has no DSC-specific
# validation in the literature.
function KUL_dsc_run_denoise {
	if [ $do_denoise -eq 0 ]; then
		mrconvert -force "$in_dsc" "$dsc_denoised"
		return
	fi
	local sp=($(mrinfo -spacing "$in_dsc"))
	local nv=$(mrinfo -size "$in_dsc" | awk '{print $4}')
	local extent=""
	if awk -v z="${sp[2]}" -v x="${sp[0]}" 'BEGIN { exit !(z > 1.5 * x) }'; then
		local e=$(awk -v n="$nv" 'BEGIN { e = int(sqrt(1.2 * n / 3)); if (e < 5) e = 5; if (e % 2 == 0) e++; while (e * e * 3 < 1.2 * n) e += 2; print e }')
		extent="-extent ${e},${e},3"
	fi
	dwidenoise -force -nthreads $ncpu $extent -noise "$workdir/dsc_noise.nii.gz" \
		"$in_dsc" "$dsc_denoised"
}

# -p BSpline, not the antsMotionCorr default of Linear. This is the first of two
# resamplings the series goes through (motion correction, then the distortion
# correction below), and linear interpolation is a triangular kernel whose
# frequency response rolls off badly -- it softens the data before the second,
# sinc-interpolated resample ever sees it, so the resolution that step is trying
# to preserve has already been given away. That blurring leaks arterial signal,
# which carries an order of magnitude more dR2* than tissue, into neighbouring
# voxels and biases peritumoral rCBV upward. Unlike a per-voxel scale factor,
# that leak does not cancel in the S(t)/S0 ratio.
#
# Rigid, not affine (method v2): head motion is rigid, and the 12-dof affine has
# shear/scale freedom that the bolus-driven intensity change can push around.
# Since v2 this is the ONLY resampling the 4D series gets.
function KUL_dsc_run_motioncorr {
	antsMotionCorr -d 3 -a "$dsc_denoised" -o "$workdir/dsc_avg.nii.gz" && \
	antsMotionCorr -d 3 -o [$workdir/dsc_mc_,$dsc_mc,$workdir/dsc_avg.nii.gz] \
		-m MI[$workdir/dsc_avg.nii.gz,$dsc_denoised,1,32,Random,0.50] \
		-i 20 -u 1 -e 1 -n 10 -t Rigid[0.1] -s 0 -f 1 \
		-p BSpline
}

# The moving image is the SKULL-STRIPPED first frame, not the raw one. The
# fixed image ($anat_regrid) is brain-only -- it is a regrid of the synthstripped
# anat -- so feeding the raw whole-head DSC frame as moving made
# --initial-moving-transform's centre-of-mass estimate compare a brain against a
# head: the DSC's skull, scalp and whatever neck falls in the slab all pull its
# COM inferiorly. The -x masks below do NOT rescue this; they restrict metric
# sampling only, and ANTs computes the COM initialisation from the image
# intensities themselves, before any mask is applied. On the validation subject the two brains
# started ~23 mm apart in both y and z, the bad initialisation put the affine
# stage in the wrong basin, and MI at shrink factor 8 never climbed out --
# final Dice against the anat brain mask was 0.39.
# NO --initial-moving-transform. The DSC and the anatomy come from the same
# session and already share a scanner frame of reference, so the identity is the
# correct starting point -- which is exactly what antsIntermodalityIntrasubject.sh
# did in tools/KUL_DSC_analysis/DSC_proc_script_WIP3.sh: grep it and you will
# find no initialisation of any kind.
#
# A centre-of-MASS init (the ",1" this used to pass) is actively harmful here.
# The DSC is a 42 x 3 mm slab and the anat is whole-head, so their brain centres
# of mass differ by ~23 mm in y and z even when the two are perfectly aligned.
# Matching those centres manufactures that offset as a displacement: the affine
# came out translating [-0.8, 46.5, 31.7] mm instead of [1.4, -5.8, 0.7].
#
# Measured on a validation subject, Dice of the DSC brain warped into anat space:
#
#   COM init + stripped moving                    0.2889   (69.6% outside)
#   COM init + raw moving (the original failure)  0.3865   (58.3% outside)
#   no init + stripped moving, params as before   0.8074   (17.5% outside)
#   no init + params below, raw moving            0.9673   ( 4.7% outside)
#   no init + params below, stripped moving       0.9668   ( 4.8% outside)  <- this
#
# The remaining 0.81 -> 0.97 comes from the parameters below, taken from the
# wrapper: histogram matching on, no metric masks, and smoothing sigmas in mm.
# That last one is worth keeping in mind before "tidying" them back to vox --
# on this 1.5x1.5x3 mm grid, 3vox is 4.5 mm in-plane but 9 mm through-plane.
#
# --restrict-deformation is retained: it costs nothing measurable (0.9673 with,
# 0.9676 without) and without it the SyN stage is an unconstrained warp rather
# than a susceptibility correction.
function KUL_dsc_run_epicorrect {
	antsRegistration --verbose $ants_verbose --dimensionality 3 --float 1 \
		--output [${epic_prefix},${epic_prefix}Warped.nii.gz] \
		--interpolation LanczosWindowedSinc \
		--winsorize-image-intensities [0.005,0.995] \
		--collapse-output-transforms 1 \
		--transform Affine[0.1] \
		--metric MI[$anat_regrid,$dsc_vol1_brain,1,32,Regular,0.25] \
		--convergence [1000x500x250x100,1e-6,10] \
		--shrink-factors 8x4x2x1 --smoothing-sigmas 4x2x1x0 \
		--use-histogram-matching 1 \
		--transform SyN[0.1,3,0] \
		--metric CC[$anat_regrid,$dsc_vol1_brain,1,4] \
		--convergence [50x50x20,1e-7,5] \
		--shrink-factors 4x2x1 --smoothing-sigmas 2x1x0mm \
		--use-histogram-matching 1 \
		--restrict-deformation $restrict_deformation
}

# Method v2 removed the step that resampled the 4D series onto the anatomy-axis
# grid (LanczosWindowedSinc) and divided out a static N4 field. The resample cut
# obliquely through 4 mm slices with gaps and displaced and blurred structure; the
# N4 field cancels in S(t)/S0 and only served a synthstrip call that now runs on
# the native median. The 3D maps are moved to the anatomy once, after the fit.


# denoise, motion correct, and estimate the DSC-to-anatomy transform
function KUL_dsc_preprocess {

	if [ ! -f "$dsc_denoised" ]; then
		task_in="KUL_dsc_run_denoise"
		KUL_task_exec $verbose_level "DSC - denoising the time series" "03_denoise" || return 1
	fi

	if [ ! -f "$dsc_mc" ]; then
		task_in="KUL_dsc_run_motioncorr"
		KUL_task_exec $verbose_level "DSC - motion correction" "04_motioncorr" || return 1
	fi

	if [ ! -f "${epic_prefix}1Warp.nii.gz" ]; then

		task_in="KUL_dsc_strip_dscvol1"
		KUL_task_exec $verbose_level "DSC - extracting the first volume" "05_dsc_vol1" || return 1

		# EPI distortion correction, as a direct rigid+affine+SyN to the
		# anatomy. The prototype used antsIntermodalityIntrasubject.sh, which
		# demands a template and a subject-to-template warp purely to emit
		# template-space extras nothing here consumes -- that made a full SyN
		# registration to MNI a prerequisite of every run. This does the same
		# EPI-to-anatomy registration directly, and adds the PE-axis
		# restriction that a susceptibility correction should have.
		task_in="KUL_dsc_run_epicorrect"
		KUL_task_exec $verbose_level "DSC - EPI distortion correction" "06_epi_correct" || return 1
	fi
}


# the masks the fit runs with, on the native DSC grid
function KUL_dsc_brainmask {

	if [ ! -f "$dsc_mask" ]; then
		task_in="KUL_dsc_make_brainmask"
		KUL_task_exec $verbose_level "DSC - computing the brain mask" "08_brainmask" || return 1
	else
		kul_echo "  the DSC brain mask already exists"
	fi

	if [ -n "$lesion_mask" ] && [ -f "$lesion_mask" ]; then
		if [ ! -f "$dsc_excl" ]; then
			task_in="KUL_dsc_make_exclusion"
			KUL_task_exec $verbose_level "DSC - lesion mask onto the DSC grid" "08b_lesion_exclusion" || return 1
		fi
	else
		kul_echo "  no lesion mask: the BSW reference and AIF search use the whole brain mask"
	fi
}


# the quantification itself
function KUL_dsc_fit {

	if [ -f "$fitdir/rCBV_corrected.nii.gz" ]; then
		kul_echo "  the DSC fit has already run"
		return 0
	fi

	local excl_opt=""
	[ -f "$dsc_excl" ] && excl_opt="--exclude_mask $dsc_excl"

	task_in="python3 ${kul_main_dir}/share/dsc/KUL_dsc_fit.py \
		$dsc_mc $dsc_mask $fitdir \
		--te $te_arg --tr $tr_arg --baseline $baseline --window $window \
		--smooth_fwhm_mm $smooth_fwhm $excl_opt"
	KUL_task_exec $verbose_level "DSC - leakage correction and perfusion fit" "09_fit" || return 1
}


# move the 3D maps from the native DSC grid to the full-resolution anatomy:
# one Linear interpolation through the composed affine + SyN transform. Linear,
# not sinc: these are derived maps, and sinc side lobes would ring at the lesion
# rim. A map is re-exported whenever the fit is newer than the delivered copy.
function KUL_dsc_maps_to_anat {

	local map
	for map in "${dsc_maps[@]}"; do
		local src="$fitdir/${map}.nii.gz"
		local dst="$resultsdir/sub-${participant}_${map}.nii.gz"
		if [ -f "$src" ] && { [ ! -f "$dst" ] || [ "$src" -nt "$dst" ]; }; then
			antsApplyTransforms -d 3 --float 1 --verbose $ants_verbose \
				-i "$src" -o "$dst" -r "$anat_ref" \
				-t "${epic_prefix}1Warp.nii.gz" -t "${epic_prefix}0GenericAffine.mat" \
				-n Linear
		fi
	done

	# where the maps are defined: the fit mask on the T1w grid. Outside it the
	# maps are 0, and the ROI stats must not count those voxels as values.
	if [ ! -f "$coverage_mask" ] || [ "$dsc_mask" -nt "$coverage_mask" ]; then
		antsApplyTransforms -d 3 --float 1 --verbose $ants_verbose \
			-i "$dsc_mask" -o "$workdir/coverage_linear.nii.gz" -r "$anat_ref" \
			-t "${epic_prefix}1Warp.nii.gz" -t "${epic_prefix}0GenericAffine.mat" \
			-n Linear
		mrcalc -quiet -force "$workdir/coverage_linear.nii.gz" 0.5 -gt "$coverage_mask" -datatype uint8
	fi
	kul_echo "  perfusion maps written to $resultsdir"
}


# co-register vendor/console maps, if the user pointed at any
function KUL_dsc_isp_maps {

	[ -z "$isp_dir" ] && return 0

	if [ ! -d "$isp_dir" ]; then
		echo "WARNING: -I '$isp_dir' is not a directory; skipping the vendor maps"
		return 0
	fi

	local isp_maps=($(find "$isp_dir" -maxdepth 1 -name "*.nii.gz" -type f | sort))
	if [ ${#isp_maps[@]} -eq 0 ]; then
		echo "WARNING: no NIfTI files in '$isp_dir'; skipping the vendor maps"
		return 0
	fi

	mkdir -p $resultsdir/vendor
	local dsc_size=$(mrinfo -size "$dsc_vol1" 2>/dev/null)
	local m n_done=0
	for m in "${isp_maps[@]}"; do
		local base=$(basename "$m" .nii.gz)
		local dst="$resultsdir/vendor/sub-${participant}_${base}_reg2_anat.nii.gz"
		[ -f "$dst" ] && [ "$dst" -nt "${epic_prefix}1Warp.nii.gz" ] && continue
		# only maps on the DSC acquisition grid share its transform; consoles also
		# export things like 512x512x1 bolus-curve screenshots
		local m_size=$(mrinfo -size "$m" 2>/dev/null | awk '{print $1, $2, $3}')
		if [ "$m_size" != "$dsc_size" ]; then
			kul_echo "  skipping vendor file $(basename $m): ${m_size// /x} is not the DSC grid (${dsc_size// /x})"
			continue
		fi
		n_done=$((n_done + 1))
		# vendor maps come off the console on the DSC grid, so the DSC-to-anat
		# transform applies to them unchanged.
		# Linear, not the LanczosWindowedSinc used for the 4D series: these are
		# already-derived parametric maps, so the sinc side lobes ring at the
		# lesion rim with no S(t)/S0 ratio to cancel them, and can hand back
		# negative rCBV on the very edge being read.
		antsApplyTransforms -d 3 --float 1 --verbose $ants_verbose \
			-i "$m" -o "$dst" -r "$anat_ref" \
			-t ${epic_prefix}1Warp.nii.gz -t ${epic_prefix}0GenericAffine.mat \
			-n Linear
	done
	kul_echo "  ${n_done} vendor map(s) registered into $resultsdir/vendor"
}


# resolve the lesion mask and the FreeSurfer aseg needed for normalisation
function KUL_dsc_resolve_normalisation_inputs {

	if [ -z "$lesion_mask" ]; then
		if [ -f "$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" ]; then
			lesion_mask="$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz"
		elif [ -f "$globalresultsdir/Lesion/lesion.nii.gz" ]; then
			lesion_mask="$globalresultsdir/Lesion/lesion.nii.gz"
		fi
	fi

	if [ -z "$fs_dir" ]; then
		local vbg_fs="$cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG/output_VBG/sub-${participant}/sub-${participant}_FS_output/sub-${participant}"
		local std_fs="$cwd/BIDS/derivatives/freesurfer/sub-${participant}"
		if [ -f "$vbg_fs/mri/aseg.mgz" ]; then
			fs_dir="$vbg_fs"
		elif [ -f "$std_fs/mri/aseg.mgz" ]; then
			fs_dir="$std_fs"
		fi
	fi

	if [ -z "$lesion_mask" ] || [ ! -f "$lesion_mask" ]; then
		echo ""
		echo "NOTE: no lesion mask found, so the NAWM normalisation is skipped."
		echo "  The perfusion maps in $resultsdir are complete and usable;"
		echo "  re-run this script once the lesion segmentation exists to add the ratios."
		return 1
	fi
	if [ -z "$fs_dir" ] || [ ! -f "$fs_dir/mri/aseg.mgz" ]; then
		echo ""
		echo "NOTE: no FreeSurfer aseg.mgz found, so the NAWM normalisation is skipped."
		echo "  It is produced by KUL_VBG (types 1-3) or KUL_FS_multiparc (types 4-6)."
		echo "  Re-run this script after those, or point at a subject dir with -F."
		return 1
	fi

	kul_echo "  lesion mask: $lesion_mask"
	kul_echo "  FreeSurfer aseg: $fs_dir/mri/aseg.mgz"
	return 0
}


# build the contralesional NAWM reference ROI (and the contralesional cortex,
# a QC reference). Every step only runs when its output is missing, but the
# laterality is always worked out, so a re-run still records it.
function KUL_dsc_make_nawm {

	# aseg.mgz lives on FreeSurfer's conformed grid; --regheader moves it onto
	# the anatomical grid using the shared scanner coordinates
	if [ ! -f $workdir/aseg_in_anat.nii.gz ]; then
		task_in="mri_vol2vol --mov $fs_dir/mri/aseg.mgz --targ $anat_ref \
			--regheader --o $workdir/aseg_in_anat.nii.gz --nearest"
		KUL_task_exec $verbose_level "DSC - resampling the aseg to the anatomical grid" "10_aseg" || return 1
	fi

	# FreeSurfer aseg labels: 2/41 cerebral WM, 3/42 cerebral cortex (L/R)
	if [ ! -f $workdir/Rhemi_GM.nii.gz ]; then
		task_in="mrcalc -force $workdir/aseg_in_anat.nii.gz 2 -eq $workdir/Lhemi_WM.nii.gz; \
			mrcalc -force $workdir/aseg_in_anat.nii.gz 41 -eq $workdir/Rhemi_WM.nii.gz; \
			mrcalc -force $workdir/aseg_in_anat.nii.gz 3 -eq $workdir/Lhemi_GM.nii.gz; \
			mrcalc -force $workdir/aseg_in_anat.nii.gz 42 -eq $workdir/Rhemi_GM.nii.gz"
		KUL_task_exec $verbose_level "DSC - splitting the hemispheres" "11_hemispheres" || return 1
	fi

	# the lesion masks are rebuilt every time: they are cheap, and the lesion
	# segmentation may have been edited since the last run
	task_in="mrcalc -force $workdir/Lhemi_GM.nii.gz $workdir/Lhemi_WM.nii.gz -add 0 -gt $workdir/Lhemi_mask.nii.gz; \
		mrcalc -force $workdir/Rhemi_GM.nii.gz $workdir/Rhemi_WM.nii.gz -add 0 -gt $workdir/Rhemi_mask.nii.gz; \
		mrcalc -force $lesion_mask 0 -gt $workdir/lesion_bin.nii.gz; \
		maskfilter -force -npass 4 $workdir/lesion_bin.nii.gz dilate $workdir/lesion_dil4.nii.gz"
	KUL_task_exec $verbose_level "DSC - preparing the hemisphere and lesion masks" "12_masks" || return 1

	local n_left=$(mrstats -output count -mask $workdir/Lhemi_mask.nii.gz -ignorezero $workdir/lesion_bin.nii.gz)
	local n_right=$(mrstats -output count -mask $workdir/Rhemi_mask.nii.gz -ignorezero $workdir/lesion_bin.nii.gz)
	n_left=${n_left:-0}
	n_right=${n_right:-0}

	if [ "$n_left" -eq 0 ] && [ "$n_right" -eq 0 ]; then
		echo "WARNING: the lesion mask does not overlap either cerebral hemisphere;"
		echo "  cannot determine laterality, so the NAWM normalisation is skipped."
		return 1
	fi

	local total=$((n_left + n_right))
	pct_left=$(KUL_float_pct "$n_left" "$total")
	pct_right=$(KUL_float_pct "$n_right" "$total")

	local ref_hemi
	if [ "$n_left" -gt "$n_right" ]; then
		lesion_side="left"
		contralesional="right"
		ref_hemi="R"
	else
		lesion_side="right"
		contralesional="left"
		ref_hemi="L"
	fi
	echo "  lesion is ${pct_left}% left / ${pct_right}% right -> sampling NAWM in the ${contralesional} hemisphere"

	# erode the reference WM away from the GM boundary and from any partial-volume
	# CSF, then exclude anything within 4 dilations of the lesion
	if [ ! -f "$nawm_mask" ]; then
		task_in="maskfilter -force -npass 4 $workdir/${ref_hemi}hemi_WM.nii.gz erode $workdir/ref_WM_ero4.nii.gz; \
			mrcalc -force $workdir/lesion_dil4.nii.gz 0 -eq $workdir/ref_WM_ero4.nii.gz -mult $nawm_mask"
		KUL_task_exec $verbose_level "DSC - building the contralesional NAWM mask" "13_nawm" || return 1
	fi

	# contralesional cortex, for the GM/WM rCBV sanity check in the report.
	# Not eroded: the cortical ribbon would vanish.
	if [ ! -f "$gm_mask" ]; then
		mrcalc -quiet -force $workdir/${ref_hemi}hemi_GM.nii.gz $workdir/lesion_dil4.nii.gz 0 -eq -mult \
			"$gm_mask" -datatype uint8
	fi

	local n_nawm=$(mrstats -output count -ignorezero $nawm_mask)
	if [ -z "$n_nawm" ] || [ "$n_nawm" -lt 100 ]; then
		echo "WARNING: the contralesional NAWM mask has only ${n_nawm:-0} voxels."
		echo "  That is too few for a stable reference; skipping the normalisation."
		rm -f $nawm_mask
		return 1
	fi
	echo "  NAWM reference: ${n_nawm} voxels"
	return 0
}


# record how the reference was chosen next to the numbers it produced
function KUL_dsc_write_reference {

	local nawm_line="$nawm_mask"
	[ $nawm_ok -eq 1 ] || nawm_line="none (no ratios; see the NOTE/WARNING above)"
	cat <<-EOT > $resultsdir/sub-${participant}_perfusion_reference.txt
	participant:            sub-${participant}
	DSC series:             $in_dsc
	anatomical reference:   $anat_ref
	lesion mask:            ${lesion_mask:-none}
	FreeSurfer aseg:        ${fs_dir:+$fs_dir/mri/aseg.mgz}
	lesion laterality:      ${pct_left:-n/a}% left / ${pct_right:-n/a}% right -> ${lesion_side:-n/a}
	NAWM reference side:    ${contralesional:-n/a} hemisphere
	NAWM mask:              $nawm_line
	DSC coverage mask:      $coverage_mask
	TE / TR:                ${te_arg} s / ${tr_arg} s
	method:                 v${dsc_method_version} (fit grid: native DSC; one interpolation to T1w)
	leakage correction:     BSW reference tissue + arrival-time correction
	rCBV integration:       ${window}
	AIF:                    cluster-selected (see fit/qc_timing.png)
	deconvolution:          oSVD (rCBF 1/s relative to the AIF, MTT s)
	smoothed maps:          $(KUL_float_gt "$smooth_fwhm" 0 && echo "*_smooth: ${smooth_fwhm} mm FWHM, display only" || echo "none")
	ROI report:             $resultsdir/sub-${participant}_perfusion_report.tsv
	fit provenance:         $fitdir/KUL_dsc_fit.json
	EOT
}


# ROI measures, the normalised maps and the report (share/dsc/KUL_dsc_report.py)
function KUL_dsc_report {

	local args="--perf_dir $resultsdir --participant $participant --anat $anat_ref \
		--coverage $coverage_mask --native_ref $fitdir/rCBV_corrected.nii.gz \
		--fit_json $fitdir/KUL_dsc_fit.json --qc_png $fitdir/qc_timing.png \
		--report_dir $cwd/REPORT"

	# the tumour compartments from KUL_anat_segment_tumor.sh, when they exist
	local L="$globalresultsdir/Lesion/sub-${participant}"
	[ -f "${L}_hdglio_lesion_solid_tissue.nii.gz" ] && \
		args="$args --roi enhancing=${L}_hdglio_lesion_solid_tissue.nii.gz"
	[ -f "${L}_hdglio_lesion_perilesional_tissue.nii.gz" ] && \
		args="$args --roi nonenhancing=${L}_hdglio_lesion_perilesional_tissue.nii.gz"
	[ -f "${L}_resseg_cavity_only.nii.gz" ] && \
		args="$args --roi cavity=${L}_resseg_cavity_only.nii.gz"
	[ -n "$lesion_mask" ] && [ -f "$lesion_mask" ] && args="$args --roi lesion=$lesion_mask"
	# operator-drawn ROIs / label maps
	[ -d "$resultsdir/ROIs" ] && args="$args --roi_dir $resultsdir/ROIs"
	if [ $nawm_ok -eq 1 ]; then
		args="$args --nawm $nawm_mask"
		[ -f "$gm_mask" ] && args="$args --ref_gm $gm_mask"
	fi

	mkdir -p $cwd/REPORT
	task_in="python3 ${kul_main_dir}/share/dsc/KUL_dsc_report.py $args"
	KUL_task_exec $verbose_level "DSC - ROI measures and report" "14_report" || return 1
	kul_echo "  ROI report: $resultsdir/sub-${participant}_perfusion_report.tsv"
	kul_echo "              $cwd/REPORT/sub-${participant}_06_DSC_report.pdf"
}


# --- MAIN ---

echo ""
echo "Starting $script v$version for sub-${participant}"

# conda env existence check (mirrors KUL_run_rsfMRI_networks.sh)
if ! conda env list | awk '{print $1}' | grep -qx "$pyfmri_env"; then
	echo "ERROR: conda env '$pyfmri_env' was not found (checked 'conda env list')." >&2
	echo "  Run the KUL_NIS installer's env-pyfmri section, or pass -y <env> to use a different one." >&2
	exit 1
fi

total_errorcount=0

# the maps KUL_dsc_fit.py produces, in the order they are reported; the
# *_smooth ones exist only with -s > 0
dsc_maps=(rCBV_corrected rCBV_uncorrected rCBF MTT Tmax TTP TT0 K1 K2 CNR \
	rCBV_corrected_smooth rCBV_uncorrected_smooth rCBF_smooth)

# intermediate files
anat_regrid="$workdir/anat_brain_regrid.nii.gz"
anat_regrid_mask="$workdir/anat_brain_regrid_mask.nii.gz"
dsc_denoised="$workdir/dsc_denoised.nii.gz"
dsc_mc="$workdir/dsc_denoised_mc.nii.gz"
dsc_vol1="$workdir/dsc_mc_vol1.nii.gz"
dsc_vol1_brain="$workdir/dsc_vol1_brain.nii.gz"
dsc_vol1_mask="$workdir/dsc_mc_vol1_brain_mask.nii.gz"
epic_prefix="$workdir/dsc_2_anat_"
dsc_mask="$workdir/dsc_brain_mask_native.nii.gz"
dsc_excl="$workdir/lesion_exclusion_native.nii.gz"
nawm_mask="$resultsdir/sub-${participant}_contralesional_NAWM_mask.nii.gz"
gm_mask="$resultsdir/sub-${participant}_contralesional_GM_mask.nii.gz"
coverage_mask="$resultsdir/sub-${participant}_dsc_coverage_mask.nii.gz"

# Work from method v1 (anatomy-grid fit, affine motion correction, cubic MP-PCA
# patches) must not be silently reused by v2, and maps exported from it must be
# replaced. A version stamp in the derivatives dir decides; the anatomy-space
# NAWM mask does not depend on the method and is kept.
dsc_method_version=2
if [ -d "$workdir" ] && [ "$(cat $derivativesdir/method_version 2>/dev/null)" != "$dsc_method_version" ]; then
	kul_echo "  previous DSC work is from an older method version: recomputing it"
	rm -rf "$workdir" "$fitdir"
	for _m in "${dsc_maps[@]}"; do
		rm -f "$resultsdir/sub-${participant}_${_m}.nii.gz" "$resultsdir/sub-${participant}_n${_m}.nii.gz"
	done
	rm -rf "$resultsdir/vendor"
	mkdir -p "$workdir" "$fitdir"
fi
echo "$dsc_method_version" > "$derivativesdir/method_version"

# STEP 1 - work out what we are processing
KUL_dsc_resolve_inputs
KUL_dsc_set_pe_restriction

# STEP 2 - prepare the anatomical reference
KUL_dsc_prepare_anat || exit 1

# STEP 3 - denoise, motion correct, distortion correct
KUL_dsc_preprocess || exit 1

# STEP 4 - the mask the fit runs in
KUL_dsc_brainmask || exit 1

# STEP 5 - the quantification
KUL_activate_conda_env "$pyfmri_env"
KUL_dsc_fit || exit 1
KUL_activate_conda_env base

# STEP 6 - deliver the maps in the anatomical (T1w) space
KUL_dsc_maps_to_anat

# STEP 7 - optional vendor/console maps
KUL_dsc_isp_maps

# STEP 8 - the contralesional NAWM reference, if we can build one
nawm_ok=0
if KUL_dsc_resolve_normalisation_inputs; then
	KUL_dsc_make_nawm && nawm_ok=1
fi
KUL_dsc_write_reference

# STEP 9 - ROI measures, normalised maps and the report (always; without a
# NAWM reference it holds absolute values only)
KUL_activate_conda_env "$pyfmri_env"
KUL_dsc_report || total_errorcount=$((total_errorcount + 1))
KUL_activate_conda_env base

echo ""
if [ $total_errorcount -eq 0 ]; then
	echo "Finished $script for sub-${participant} - results in $resultsdir"
else
	echo "Finished $script for sub-${participant} with $total_errorcount error(s) - check $KUL_LOG_DIR"
fi

exit 0
