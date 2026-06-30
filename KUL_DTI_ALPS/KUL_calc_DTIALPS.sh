#!/bin/bash
# Bash shell script to calculate DTI_ALPS index using predefined 5mm radius spheres in MNI2009 space
#
# Requires Mrtrix3, FSL, ants
#
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@uzleuven.be
#

#set -x

v="v0.1 - dd 28/12/2025"

# Assumes that you've already run fmriprep, KUL_dwiprep.sh, KUL_dwiprep_anat.sh and KUL_dwiprep_MNI.sh
# Start by using KUL_clinical_fmridti.sh, then run KUL_dwiprep_MNI.sh and then this script
#   - Apply inv warps to bring MNI Spheres for DTI ALPS calculation to native dMRI space
#   - Calculate DTI ALPS and saves them to a text file per subject
# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl (for logging)
#
# this script uses "preprocessing control", i.e. if some steps are already processed it will skip these
KUL_ALPS_maindir=`dirname "$0"`
cwd=$(pwd)
# FUNCTIONS --------------
# function Usage
function Usage {
cat <<USAGE
`basename $0` Calculate DTI_ALPS index.
Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001 -n 6 

Required arguments:

     -p:  participant (anonymised name of the subject)

Optional arguments:

     -n:  number of cpu for parallelisation
     -v:  show output from mrtrix commands


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
ncpu=6
silent=1

# Set required options
p_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:n:vh" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            subj=$OPTARG
        ;;
        n) #parallel
            ncpu=$OPTARG
        ;;
        v) #verbose
            silent=0
        ;;
        h) #help
            Usage >&2
            exit 0
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
    echo "Option -p is required: give the anonymised name of a subject (this will create a directory subject_preproc with results)." >&2
    echo
    exit 2 
fi 

# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then 

    export MRTRIX_QUIET=1

fi

# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Some parallelisation
FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=log/logDTI_ALPS_${d}.txt

# --- MAIN ----------------

bids_subj=BIDS/sub-${subj}

ALPS_out=${cwd}/dwiprep/sub-${subj}/sub-${subj}/DTI_ALPS

dMRI_data=${cwd}/dwiprep/sub-${subj}/sub-${subj}/dwi_preproced.mif

native_fa=${cwd}/dwiprep/sub-${subj}/sub-${subj}/qa/fa.nii.gz

mkdir -p ${ALPS_out}

echo "id,Ass_L_Dx,Ass_L_Dz,Ass_R_Dx,Ass_R_Dz,Proj_L_Dx,Proj_L_Dy,Proj_R_Dx,Proj_R_Dy,ALPS_L,ALPS_R,ALPS_mean" > ${ALPS_out}/results.txt

## here we find the fmriprep transform from MNI to T1
function find_first_match {

    local glob="$1"
    local description="$2"
    local matches=()

    shopt -s nullglob
    matches=($glob)
    shopt -u nullglob

    if [ ${#matches[@]} -eq 0 ]; then
        echo "Error: could not locate ${description} (pattern: ${glob})" >&2
        exit 1
    fi

    printf '%s\n' "${matches[0]}"
}

MNI2T1transform=$(find_first_match "${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_*from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5" "fmriprep MNI-to-T1w transform")
dMRI2T1transform=${cwd}/dwiprep/sub-${subj}/sub-${subj}/dwi_reg/rigid_out0GenericAffine.mat

## first we bring the VOIs to native space
antsApplyTransforms -d 3 -i ${KUL_ALPS_maindir}/MNI_ROI_1_L.nii.gz -o ${ALPS_out}/NAT_ROI_1_L.nii.gz -r ${native_fa} -t [${dMRI2T1transform},1] -t ${MNI2T1transform} -n GenericLabel
antsApplyTransforms -d 3 -i ${KUL_ALPS_maindir}/MNI_ROI_2_L.nii.gz -o ${ALPS_out}/NAT_ROI_2_L.nii.gz -r ${native_fa} -t [${dMRI2T1transform},1] -t ${MNI2T1transform} -n GenericLabel
antsApplyTransforms -d 3 -i ${KUL_ALPS_maindir}/MNI_ROI_1_R.nii.gz -o ${ALPS_out}/NAT_ROI_1_R.nii.gz -r ${native_fa} -t [${dMRI2T1transform},1] -t ${MNI2T1transform} -n GenericLabel
antsApplyTransforms -d 3 -i ${KUL_ALPS_maindir}/MNI_ROI_2_R.nii.gz -o ${ALPS_out}/NAT_ROI_2_R.nii.gz -r ${native_fa} -t [${dMRI2T1transform},1] -t ${MNI2T1transform} -n GenericLabel

# probably better to collapse those to COMs and remake the spheres in native space, no?

## second we generate the DT based on available preprocessed dMRI data and calculate DTI_ALPS
# first we recalc the tensor
dwi2tensor -mask ${cwd}/dwiprep/sub-${subj}/sub-${subj}/dwi_mask.nii.gz ${dMRI_data} ${ALPS_out}/DT_4_ALPS.mif -force

# then we extract the elements of interest
mrconvert ${ALPS_out}/DT_4_ALPS.mif -coord 3 0 ${ALPS_out}/Dx.nii.gz -force
mrconvert ${ALPS_out}/DT_4_ALPS.mif -coord 3 1 ${ALPS_out}/Dy.nii.gz -force
mrconvert ${ALPS_out}/DT_4_ALPS.mif -coord 3 2 ${ALPS_out}/Dz.nii.gz -force
Dx=${ALPS_out}/Dx.nii.gz
Dy=${ALPS_out}/Dy.nii.gz
Dz=${ALPS_out}/Dz.nii.gz
 

id=${subj}

# # clean up the warped spheres (to make sure we only work further with the spheres)
# ImageMath 3 ${ALPS_out}/NAT_ROI_1_L_lcc.nii.gz GetLargestComponent ${ALPS_out}/NAT_ROI_1_L.nii.gz  
# ImageMath 3 ${ALPS_out}/NAT_ROI_2_L_lcc.nii.gz GetLargestComponent ${ALPS_out}/NAT_ROI_2_L.nii.gz
# ImageMath 3 ${ALPS_out}/NAT_ROI_1_R_lcc.nii.gz GetLargestComponent ${ALPS_out}/NAT_ROI_1_R.nii.gz
# ImageMath 3 ${ALPS_out}/NAT_ROI_2_R_lcc.nii.gz GetLargestComponent ${ALPS_out}/NAT_ROI_2_R.nii.gz

# get COMS
read cx1l cy1l cz1l <<< "$(fslstats "${ALPS_out}/NAT_ROI_1_L.nii.gz" -C)"
read cx2l cy2l cz2l <<< "$(fslstats "${ALPS_out}/NAT_ROI_2_L.nii.gz" -C)"
read cx1r cy1r cz1r <<< "$(fslstats "${ALPS_out}/NAT_ROI_1_R.nii.gz" -C)"
read cx2r cy2r cz2r <<< "$(fslstats "${ALPS_out}/NAT_ROI_2_R.nii.gz" -C)"

# remake the spheres
# Sphere 1 L
fslmaths "${native_fa}" -mul 0 -add 1 -roi "${cx1l}" 1 "${cy1l}" 1 "${cz1l}" 1 0 1 ${ALPS_out}/NAT_ROI_1_L_point.nii.gz -odt float
fslmaths ${ALPS_out}/NAT_ROI_1_L_point.nii.gz -kernel sphere 5 -fmean -thr 1e-5 -bin ${ALPS_out}/NAT_ROI_1_L_newsphere.nii.gz -odt char
fslmaths ${ALPS_out}/NAT_ROI_1_L_newsphere.nii.gz -mas ${cwd}/dwiprep/sub-${subj}/sub-${subj}/dwi_mask.nii.gz ${ALPS_out}/Proj_L_sphere.nii.gz
# Sphere 2 L
fslmaths "${native_fa}" -mul 0 -add 1 -roi "${cx2l}" 1 "${cy2l}" 1 "${cz2l}" 1 0 1 ${ALPS_out}/NAT_ROI_2_L_point.nii.gz -odt float
fslmaths ${ALPS_out}/NAT_ROI_2_L_point.nii.gz -kernel sphere 5 -fmean -thr 1e-5 -bin ${ALPS_out}/NAT_ROI_2_L_newsphere.nii.gz -odt char
fslmaths ${ALPS_out}/NAT_ROI_2_L_newsphere.nii.gz -mas ${cwd}/dwiprep/sub-${subj}/sub-${subj}/dwi_mask.nii.gz ${ALPS_out}/Ass_L_sphere.nii.gz
# Sphere 1 R
fslmaths "${native_fa}" -mul 0 -add 1 -roi "${cx1r}" 1 "${cy1r}" 1 "${cz1r}" 1 0 1 ${ALPS_out}/NAT_ROI_1_R_point.nii.gz -odt float
fslmaths ${ALPS_out}/NAT_ROI_1_R_point.nii.gz -kernel sphere 5 -fmean -thr 1e-5 -bin ${ALPS_out}/NAT_ROI_1_R_newsphere.nii.gz -odt char
fslmaths ${ALPS_out}/NAT_ROI_1_R_newsphere.nii.gz -mas ${cwd}/dwiprep/sub-${subj}/sub-${subj}/dwi_mask.nii.gz ${ALPS_out}/Proj_R_sphere.nii.gz
# Sphere 2 R
fslmaths "${native_fa}" -mul 0 -add 1 -roi "${cx2r}" 1 "${cy2r}" 1 "${cz2r}" 1 0 1 ${ALPS_out}/NAT_ROI_2_R_point.nii.gz -odt float
fslmaths ${ALPS_out}/NAT_ROI_2_R_point.nii.gz -kernel sphere 5 -fmean -thr 1e-5 -bin ${ALPS_out}/NAT_ROI_2_R_newsphere.nii.gz -odt char
fslmaths ${ALPS_out}/NAT_ROI_2_R_newsphere.nii.gz -mas ${cwd}/dwiprep/sub-${subj}/sub-${subj}/dwi_mask.nii.gz ${ALPS_out}/Ass_R_sphere.nii.gz


Ass_L=${ALPS_out}/Ass_L_sphere.nii.gz
Ass_R=${ALPS_out}/Ass_R_sphere.nii.gz
Proj_R=${ALPS_out}/Proj_R_sphere.nii.gz
Proj_L=${ALPS_out}/Proj_L_sphere.nii.gz

# defining command to trim trailing spaces
trim() { awk '{$1=$1}1'; }

Ass_L_Dx=$(mrstats $Dx -mask $Ass_L -output mean -quiet | trim)
Ass_L_Dz=$(mrstats $Dz -mask $Ass_L -output mean -quiet | trim)
Ass_R_Dx=$(mrstats $Dx -mask $Ass_R -output mean -quiet | trim)
Ass_R_Dz=$(mrstats $Dz -mask $Ass_R -output mean -quiet | trim)
Proj_L_Dx=$(mrstats $Dx -mask $Proj_L -output mean -quiet | trim)
Proj_L_Dy=$(mrstats $Dy -mask $Proj_L -output mean -quiet | trim)
Proj_R_Dx=$(mrstats $Dx -mask $Proj_R -output mean -quiet | trim)
Proj_R_Dy=$(mrstats $Dy -mask $Proj_R -output mean -quiet | trim)


ALPS_L=$(awk -v dxp="${Proj_L_Dx}" -v dxa="${Ass_L_Dx}" -v dyp="${Proj_L_Dy}" -v dza="${Ass_L_Dz}" 'BEGIN{num=(dxp+dxa); den=(dyp+dza); if(den==0){print "nan"} else {printf "%.6f", num/den}}')
ALPS_R=$(awk -v dxp="${Proj_R_Dx}" -v dxa="${Ass_R_Dx}" -v dyp="${Proj_R_Dy}" -v dza="${Ass_R_Dz}" 'BEGIN{num=(dxp+dxa); den=(dyp+dza); if(den==0){print "nan"} else {printf "%.6f", num/den}}')
ALPS_mean=$(awk -v a="${ALPS_L}" -v b="${ALPS_R}" 'BEGIN{ if(a=="nan" || b=="nan"){print "nan"} else {printf "%.6f", (a+b)/2}}')

read ALPS_L ALPS_R ALPS_mean <<< "$(python3 -c '
import sys, math
def f(x):
    try: return float(x.strip())
    except: return float("nan")

pLDx,aLDx,pLDy,aLDz,pRDx,aRDx,pRDy,aRDz = [f(x) for x in sys.argv[1:9]]

def alps(dxp,dxa,dyp,dza):
    den=dyp+dza
    num=dxp+dxa
    v=num/den if den!=0 and math.isfinite(num) and math.isfinite(den) else float("nan")
    return v

L=alps(pLDx,aLDx,pLDy,aLDz)
R=alps(pRDx,aRDx,pRDy,aRDz)
M=(L+R)/2 if math.isfinite(L) and math.isfinite(R) else float("nan")

def out(v): return "nan" if math.isnan(v) else f"{v:.6f}"
print(out(L), out(R), out(M))
' "$Proj_L_Dx" "$Ass_L_Dx" "$Proj_L_Dy" "$Ass_L_Dz" "$Proj_R_Dx" "$Ass_R_Dx" "$Proj_R_Dy" "$Ass_R_Dz")"

# DTI_ALPS = (mean (Dxx_proj, Dxx_ass))/(mean (Dyy_proj, Dzz_proj))

echo "$id,$Ass_L_Dx,$Ass_L_Dz,$Ass_R_Dx,$Ass_R_Dz,$Proj_L_Dx,$Proj_L_Dy,$Proj_R_Dx,$Proj_R_Dy,$ALPS_L,$ALPS_R,$ALPS_mean" >> ${ALPS_out}/results.txt

# done

