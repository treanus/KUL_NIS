#!/bin/bash
# Bash shell script to prepare fMRI/DTI results for Brainlab Elements Server
#
# Requires Mrtrix3, Karawun
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 12/11/2021
version="0.3"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` is a script that prepares data for input into Brainlab Elements

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -t 3 -r 40

Required arguments:

     -p:  participant name

Optional arguments:

     -t:  processing type
        type 1: (DEFAULT) prepare for tumor patient
        type 2: prepare for a ET DBS patient (DRT)
        type 3: prepare for a Parkinson DBS patient (CSHDP)
     -r:  use a relative treshold (in percent of tract density)
     -a:  use the ACT output
     -v:  show output from commands

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ncpu=15
type=1
act_type=0
relative=0
treshold=0 

# Set required options
p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:r:va" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        t) #type
			type=$OPTARG
		;;
        r) #relative threshold
            relative=1
			threshold=$OPTARG
		;;
        v) #verbose
			silent=0
		;;
        a) #ACT
			act_type=1
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

# MRTRIX and others verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
    str_silent=" > /dev/null 2>&1" 
    ants_verbose=0
fi

if [ $act_type -eq 1 ] ; then
    ACT="_ACT"
else
    ACT=""
fi

#----- functions

function KUL_karawun_get_tract {
    local tract_dir="BIDS/derivatives/KUL_compute/sub-${participant}/FWT/sub-${participant}_TCKs_output/${tract_name_orig}_output"
    local fin_tck="${tract_dir}/${tract_name_orig}_fin_BT${ACT}_iFOD2.tck"
    local fin_map="${tract_dir}/${tract_name_orig}_fin_map_BT${ACT}_iFOD2.nii.gz"
    local use_tck=""
    local use_map=""

    if [ -f "$fin_tck" ]; then
        use_tck="$fin_tck"
    else
        use_tck=$(ls "${tract_dir}/${tract_name_orig}_filt"*"_BT${ACT}_iFOD2.tck" 2>/dev/null | grep -v "_inMNI" | sort -V | tail -1)
        if [ -n "$use_tck" ]; then
            echo "  ${tract_name_orig}: fin not found, using $(basename $use_tck)"
        fi
    fi

    if [ -f "$fin_map" ]; then
        use_map="$fin_map"
    else
        use_map=$(ls "${tract_dir}/${tract_name_orig}_filt"*"_map_BT${ACT}_iFOD2.nii.gz" "${tract_dir}/${tract_name_orig}_filt"*"_BT${ACT}_iFOD2_map.nii.gz" 2>/dev/null | grep -v "_inMNI" | sort -V | tail -1)
        if [ -n "$use_map" ]; then
            echo "  ${tract_name_orig}: fin map not found, using $(basename "$use_map")"
        fi
    fi

    if [ -n "$use_tck" ]; then
        cp "$use_tck" Karawun/sub-${participant}/tck/${tract_name_final}.tck

        if [ $type -eq 1 ]; then

            num_tck=$(tckstats -quiet -output count Karawun/sub-${participant}/tck/${tract_name_final}.tck)
            echo $num_tck
            tract_threshold=$(( num_tck*1/3/100 ))
            echo $tract_threshold

        fi

        if [ $relative -eq 1 ]; then

            num_tck=$(tckstats -quiet -output count Karawun/sub-${participant}/tck/${tract_name_final}.tck)
            echo "$tract_name_final has $num_tck streamlines"
            tract_threshold=$(( num_tck*threshold/tract_corr_threshold/100 ))
            echo "The compute threshold is: $tract_threshold (with a correction of $tract_corr_threshold)"

        fi

        if [ -n "$use_map" ]; then
            mrgrid "$use_map" \
                regrid -template Karawun/sub-${participant}/T1w.nii.gz \
                -interp linear \
                - | mrcalc - ${tract_threshold} -gt ${tract_color} -mul \
                Karawun/sub-${participant}/labels/${tract_name_final}_center.nii.gz -force
        else
            echo "  Warning: no map found for ${tract_name_orig}, label not generated"
        fi
    else
        echo "Does not exist: $fin_tck"
    fi
}

function KUL_karawun_get_voi {
    if [ -f BIDS/derivatives/KUL_compute/sub-${participant}/FWT/sub-${participant}_VOIs/${tract_name_orig}_VOIs/${tract_name_orig}_incs1/${tract_name_orig}_incs1_map.nii.gz ]; then
        mrgrid BIDS/derivatives/KUL_compute/sub-${participant}/FWT/sub-${participant}_VOIs/${tract_name_orig}_VOIs/${tract_name_orig}_incs1/${tract_name_orig}_incs1_map.nii.gz \
            regrid -template Karawun/sub-${participant}/T1w.nii.gz \
            -interp linear \
            - | mrcalc - ${voi_threshold} -gt ${voi_color} -mul \
            Karawun/sub-${participant}/labels/${voi_name_final}.nii.gz -force
    else
        echo "Does not exist: BIDS/derivatives/KUL_compute/sub-${participant}/FWT/sub-${participant}_VOIs/${tract_name_orig}_VOIs/${tract_name_orig}_incs1/${tract_name_orig}_incs1_map.nii.gz"
    fi
}

# Known tract display-name/color/threshold metadata, keyed by the FWT bundle name
# (tract_name_orig). Used to look up nice Brainlab labels/colors for bundles that
# are actually found (see the auto-discovery loop in MAIN below) instead of the
# previous approach of unconditionally attempting a hardcoded list of ~40 bundles
# regardless of which ones this patient's FWT config actually generated -- that
# meant any bundle not already in this list (e.g. a newly added FWT bundle, or a
# custom one) was silently never picked up, no matter what -t was set to (-t does
# not actually filter this list at all; it only changes the threshold formula for
# type 1, see KUL_karawun_get_tract above).
# Color convention: value = Brainlab's RecommendedDisplayCIELabValue index
# (karawun's lookup_cie(), 31-entry palette, indices 0-30 -- anything >30
# silently clamps to the same last color, so 30 is the real usable ceiling).
# For most tract families, left/right share the same color -- side is
# already obvious from spatial position in the 3D view, so color instead
# encodes tract *type*, which is what actually matters for interpretation,
# and halving the color budget leaves room for fMRI labels without
# colliding. Exception: tracts that pass through/near the brainstem or are
# otherwise high-stakes/eloquent (CST, ML, PyT_SMA) keep separate L/R
# colors -- the brainstem packs both sides close to midline with much less
# spatial separation than cortical/subcortical tracts, so color is doing
# real disambiguation work there that position alone doesn't replace.
declare -A KUL_karawun_tract_meta=(
    [AF_all_LT]="Arcuate_Fasc_Left|1|20|3"
    [AF_all_RT]="Arcuate_Fasc_Right|1|20|3"
    [CST_LT]="Corticospinal_Tract_Left|3|20|3"
    [CST_RT]="Corticospinal_Tract_Right|4|20|3"
    [CCing_LT]="Cingulum_cing_Left|5|20|3"
    [TCing_LT]="Cingulum_temporal_Left|5|20|3"
    [CCing_RT]="Cingulum_cing_Right|5|20|3"
    [TCing_RT]="Cingulum_temporal_Right|5|20|3"
    [FAT_LT]="FrontalAslant_Tract_Left|7|20|3"
    [FAT_RT]="FrontalAslant_Tract_Right|7|20|3"
    [IFOF_LT]="IFOF_Left|9|20|3"
    [IFOF_RT]="IFOF_Right|9|20|3"
    [ILF_LT]="InferiorLongitudinal_Fasc_Left|11|20|3"
    [ILF_RT]="InferiorLongitudinal_Fasc_Right|11|20|3"
    [UF_LT]="Uncinate_Fasc_Left|13|20|3"
    [UF_RT]="Uncinate_Fasc_Right|13|20|3"
    [OR_occlobe_LT]="Occiptal_Radition_Left|15|20|3"
    [OR_occlobe_RT]="Occiptal_Radition_Right|15|20|3"
    [ML_LT]="Medial_Lemniscus_Tract_Left|17|20|3"
    [ML_RT]="Medial_Lemniscus_Tract_Right|18|20|3"
    [MdLF_LT]="MiddleLongitudinal_Fasc_Left|29|20|3"
    [MdLF_RT]="MiddleLongitudinal_Fasc_Right|29|20|3"
    [Ant_Comm]="Anterior_Commissure|31|20|3"
    [Post_Comm]="Posterior_Commissure|32|20|3"
    [CC_Motor_Comm]="CC_Motor|33|20|3"
    [CC_Occipital_Comm]="CC_Occipital|34|20|3"
    [CC_Parietal_Comm]="CC_Parietal|35|20|3"
    [CC_PMandSM_Comm]="CC_PreMotor_SupplMotor|36|20|3"
    [CC_PreF_Comm]="CC_Prefrontal|37|20|3"
    [CC_Sensory_Comm]="CC_Sensory|38|20|3"
    [CC_Temporal_Comm]="CC_Temporal|39|20|3"
    [DRT_LT]="DRT_Left|19|20|1"
    [DRT_RT]="DRT_Right|20|20|1"
    [ThR_S1_LT]="S1VC_Left|21|10|4"
    [ThR_S1_RT]="S1VC_Right|22|10|4"
    [CSHDP_LT]="CSHDP_Left|25|40|1"
    [CSHDP_RT]="CSHDP_Right|26|40|1"
    [ThR_Ant_LT]="ThR_Ant_LT|27|20|1"
    [ThR_Ant_RT]="ThR_Ant_RT|28|20|1"
    # Placed at the end with their own unused colors (were previously
    # colliding with CST_LT/RT's 3/4, swapped no less) rather than reusing
    # any of the colors above, so existing/more heavily relied-on bundles
    # like CST keep their familiar colors undisturbed.
    [PyT_SMA_LT]="Pyramidal_Tract_Supplementary_Motor_Left|40|20|3"
    [PyT_SMA_RT]="Pyramidal_Tract_Supplementary_Motor_Right|41|20|3"
)

function KUL_karawun_auto_discover_tracts {
    # Mirrors the auto-discovery pattern KUL_clinical_fmridti.sh already uses for
    # screenshots/PACS export: glob over whatever *_output directories FWT actually
    # produced for this participant, rather than assuming a fixed list. Any bundle
    # found gets its known display name/color/thresholds from the table above if
    # present, or a sensible auto-assigned fallback (raw FWT name, next unused
    # color starting at 100) if it's a bundle this script doesn't recognise yet --
    # so a new or custom bundle in the FWT config still reaches Karawun instead of
    # being silently skipped.
    local tck_root="BIDS/derivatives/KUL_compute/sub-${participant}/FWT/sub-${participant}_TCKs_output"
    # 42 = one past the highest index the known-tract table uses (41), so
    # auto-assigned bundles start where the curated ones stop.
    #
    # This was 100, which karawun's lookup_cie() clamps to the last palette
    # entry -- so every auto-discovered bundle came out the same colour as every
    # other one, and as the 13 table entries that also sit above the palette
    # size. Against the stock 31-colour palette 42 still clamps, so this is no
    # worse today; against an extended palette it does the right thing.
    local next_auto_color=42

    if [ ! -d "$tck_root" ]; then
        echo "Does not exist: $tck_root"
        return
    fi

    for tck_outdir in "$tck_root"/*_output; do
        [ -d "$tck_outdir" ] || continue
        tract_name_orig=$(basename "$tck_outdir" _output)

        if [ -n "${KUL_karawun_tract_meta[$tract_name_orig]:-}" ]; then
            IFS='|' read -r tract_name_final tract_color tract_threshold tract_corr_threshold \
                <<< "${KUL_karawun_tract_meta[$tract_name_orig]}"
        else
            echo "  ${tract_name_orig}: not in the known tract table, using defaults (color ${next_auto_color})"
            tract_name_final="$tract_name_orig"
            tract_color=$next_auto_color
            tract_threshold=20
            tract_corr_threshold=3
            next_auto_color=$((next_auto_color + 1))
        fi

        KUL_karawun_get_tract
    done
}

#---- MAIN

mkdir -p Karawun/sub-${participant}/labels
mkdir -p Karawun/sub-${participant}/tck
mkdir -p Karawun/sub-${participant}/DICOM

T1w_in="RESULTS/sub-${participant}/Anat/T1w.nii.gz"
T1w_min=$(mrstats -output min $T1w_in)
T1w_max=$(mrstats -output max $T1w_in)
T1w_factor=$(scale=10; echo "($T1w_max-($T1w_min))/32767" | bc)
mrcalc $T1w_in $T1w_min -sub $T1w_factor -div Karawun/sub-${participant}/T1w.nii.gz -force

# FAT1w = sqrt(FA) * T1w (Goedemans et al., Imaging Neurosci 2024), loaded into
# Brainlab as a second anatomical alongside the T1w. It makes the FA-to-T1w
# registration directly inspectable, which is the QA step this whole folder
# depends on -- if that registration has slipped, every tract in the scene is
# wrong in the same direction and nothing else here would show it.
#
# KUL_FAT1w.py had become orphaned: nothing called it any more, so the file this
# block looks for never existed and FAT1w was silently always empty. Generate it
# here instead of assuming some earlier step did.
FAT1="BIDS/derivatives/KUL_compute/sub-${participant}/KUL_FAT1/FAT1w.nii.gz"
FA_reg2T1w="dwiprep/sub-${participant}/sub-${participant}/qa/fa_reg2T1w.nii.gz"

if [ ! -f "$FAT1" ]; then
    if [ -f "$FA_reg2T1w" ] && [ -f "$T1w_in" ]; then
        echo "Computing FAT1w (FA-weighted T1w) for Brainlab QA"
        # no -s: smoothing blurs exactly the tissue boundaries this image is
        # meant to let you check the registration against
        "${kul_main_dir}/KUL_FAT1w.py" -p "${participant}" || \
            echo "WARNING: KUL_FAT1w.py failed - continuing without the FAT1w QA volume"
    else
        echo "No $FA_reg2T1w (run KUL_dwiprep_anat.sh first) - skipping the FAT1w QA volume"
    fi
fi

if [ -f $FAT1 ]; then
    FAT1w="Karawun/sub-${participant}/FAT1w.nii.gz"

    # Rescale into 16-bit range instead of copying as-is. karawun writes
    # SmallestImagePixelValue/LargestImagePixelValue (0028,0106/0107) from the
    # *original* intensities, not from the values it rescales the pixel data to,
    # and both tags are US -- capped at 65535. FAT1w is sqrt(FA) * T1w, so it
    # comes out well above that (97575 on the first real case), and
    # importTractography aborts with
    #   "'H' format requires 0 <= number <= 65535 ... (0028,0107) US: 97575".
    # The T1w a few lines above is rescaled for the same reason.
    fat1_min=$(mrstats -output min "$FAT1")
    fat1_max=$(mrstats -output max "$FAT1")
    # awk, not bc: mrstats can report in scientific notation, which bc cannot parse
    fat1_factor=$(awk -v a="$fat1_max" -v b="$fat1_min" \
        'BEGIN { d = (a - b) / 32767; if (d <= 0) d = 1; printf "%.10f", d }')
    mrcalc "$FAT1" "$fat1_min" -sub "$fat1_factor" -div "$FAT1w" -force
else
    FAT1w=""
fi

# Lesion as a Brainlab label, so the tumour/cavity shows up in the same scene as
# the tracts. Colour 16 is free in the 0-30 palette (see the colour convention
# note above) and sits clear of the fMRI task labels, which count up from 1.
lesion_color=16
lesion_in=""
for _lesion_cand in \
    "RESULTS/sub-${participant}/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" \
    "RESULTS/sub-${participant}/Lesion/lesion.nii.gz"; do
    if [ -f "$_lesion_cand" ]; then
        lesion_in="$_lesion_cand"
        break
    fi
done

if [ -n "$lesion_in" ]; then
    echo "Karawun lesion label: $(basename "$lesion_in") -> colour ${lesion_color}"
    # GenericLabel, not mrgrid's nearest: it interpolates each label's indicator
    # and takes the argmax, so boundaries follow the anatomy instead of the grid
    # and no label is invented that wasn't in the input. Today this is a no-op --
    # Karawun's T1w.nii.gz is the RESULTS T1w with only its intensities rescaled,
    # so the grids are identical and any interpolator is exact -- but that is an
    # assumption about an upstream step, and nearest is the option that degrades
    # worst if it ever stops holding. No -t: identity transform, resample only.
    _lesion_tmp="$(mktemp -d)"
    mrcalc "$lesion_in" 0 -gt "${_lesion_tmp}/lesion_bin.nii.gz" -force -quiet && \
    antsApplyTransforms -d 3 --float 1 --verbose 0 \
        -i "${_lesion_tmp}/lesion_bin.nii.gz" \
        -r Karawun/sub-${participant}/T1w.nii.gz \
        -o "${_lesion_tmp}/lesion_rs.nii.gz" \
        -n GenericLabel && \
    mrcalc "${_lesion_tmp}/lesion_rs.nii.gz" 0 -gt ${lesion_color} -mult \
        Karawun/sub-${participant}/labels/Lesion.nii.gz -force -quiet
    rm -rf "${_lesion_tmp}"
else
    echo "No lesion mask found in RESULTS/sub-${participant}/Lesion/ - no Karawun lesion label"
fi


    KUL_karawun_auto_discover_tracts

    # CSHDP is the one FWT bundle used both as a VOI-derived label (the distal
    # STN-motor VOI itself) and, separately, as a full reconstructed tract -- the
    # VOI form isn't produced under TCKs_output/*_output so it can't be picked up
    # by the auto-discovery loop above and stays explicit here. Only relevant for
    # type 3 (DBS Parkinson/CSHDP) -- gated so tumor/DBS-ET runs don't print
    # spurious "Does not exist" lines for a VOI that was never computed for them.
    if [ $type -eq 3 ]; then
        tract_name_orig="CSHDP_LT"
        voi_name_final="DISTAL_STN_MOTOR_Left"
        voi_color=23
        voi_threshold=0.1
        KUL_karawun_get_voi

        tract_name_orig="CSHDP_RT"
        voi_name_final="DISTAL_STN_MOTOR_Right"
        voi_color=24
        voi_threshold=0.1
        KUL_karawun_get_voi
    fi

# give information
echo "See to it that Karawun/sub-${participant}/DICOM/ contains a donor DICOM"
echo "(a single file is enough) before importing - use one slice from a"
echo "high-resolution anatomical series (T1w, FLAIR, T2, ...). Use the same"
echo "donor for the PACS export too, so the metadata stays consistent."
echo "Then copy into terminal: "
echo "conda activate KarawunDev"
echo "importTractography -d Karawun/sub-${participant}/DICOM/*.dcm \
-o Karawun/sub-${participant}/sub-${participant}_for_elements \
-n Karawun/sub-${participant}/T1w.nii.gz $FAT1w \
-t Karawun/sub-${participant}/tck/*.tck \
-l Karawun/sub-${participant}/labels/*.gz"
