#!/bin/bash -e
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 09/11/2018 - alpha version
v="v0.1 - dd 09/11/2018"

# To Do
#  - register dwi to T1 with ants-syn
#  - fod calc msmt-5tt in stead of dhollander


# A few fixed (for now) parameters:

    # Specify additional options for FSL eddy
    # eddy_options="--data_is_shelled --slm=linear --repol "
    eddy_options="--slm=linear --repol "

    # Number of desired streamlines
    nods=2000


# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl (for logging)
#
# this script uses "preprocessing control", i.e. if some steps are already processed it will skip these

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs dMRI segmentation of the Dentato-rubro-thalamic tract in thalamus for DBS target selection.

Usage:

  `basename $0` -s subject <OPT_ARGS>

Example:

  `basename $0` -s pat001 -p 6 -d pat001.zip 

Required arguments:

     -s:  subject (anonymised name of the subject)

Optional arguments:

     -p:  number of cpu for parallelisation
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
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "s:p:d:vh" OPT; do

        case $OPT in
        s) #subject
            s_flag=1
            subj=$OPTARG
        ;;
        p) #parallel
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
if [ $s_flag -eq 0 ] ; then 
    echo 
    echo "Option -s is required: give the anonymised name of a subject (this will create a directory subject_preproc with results)." >&2
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

# Directory to write preprocessed data in
preproc=dwiprep/sub-${subj}

# Directory to put raw mif data in
raw=${preproc}/raw

# set up preprocessing & logdirectory
mkdir -p ${preproc}/raw
mkdir -p ${preproc}/log

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=log/log_${d}.txt



# SAY HELLO ---

kul_e2cl "Welcome to KUL_dwi_preproc $v - $d" ${preproc}/${log}

bids_subj=BIDS/"sub-$subj"/ses-tp1

# STEP 1 - CONVERSION of BIDS to MIF ---------------------------------------------

# test if conversion has been done
if [ ! -f ${preproc}/dwi_orig.mif ]; then

    kul_e2cl " Preparing datasets from BIDS directory..." ${preproc}/${log}

    # convert raw T1w data, using -strides 1:3 to get orientation correct for FSL
    bids_t1w="$bids_subj/anat/sub-${subj}_T1w.nii.gz"
    #mrconvert "$bids_t1w" ${raw}/T1w.nii.gz -strides 1:3 -force -nthreads $ncpu -clear_property comments -nthreads $ncpu

    # convert raw FLAIR data, using -strides 1:3 to get orientation correct for FSL
    bids_flair="$bids_subj/anat/sub-${subj}_FLAIR.nii.gz"
    #mrconvert "$bids_flair" ${raw}/FLAIR.nii.gz -strides 1:3 -force -nthreads $ncpu -clear_property comments -nthreads $ncpu

    # convert dwi
    bids_dwi_search="$bids_subj/dwi/sub-*_dwi.nii.gz"
    bids_dwi_found=$(ls $bids_dwi_search)
    
    number_of_bids_dwi_found=$(echo $bids_dwi_found | wc -w)

    if [ $number_of_bids_dwi_found -eq 1 ]; then

        kul_e2cl "   only 1 dwi dataset, scaling not necessary" ${preproc}/${log}
        dwi_base=${bids_dwi_found%%.*}
        mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
        -json_import ${dwi_base}.json ${preproc}/dwi_orig.mif -strides 1:3 -force -clear_property comments -nthreads $ncpu 

    else 
        
        kul_e2cl "   found $number_of_bids_dwi_found dwi datasets, scaling & catting" ${preproc}/${log}
        
        dwi_i=1
        for dwi_file in $bids_dwi_found; do
            dwi_base=${dwi_file%%.*}
        
            mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
            -json_import ${dwi_base}.json ${raw}/dwi_p${dwi_i}.mif -strides 1:3 -force -clear_property comments -nthreads $ncpu

            dwiextract -quiet -bzero ${raw}/dwi_p${dwi_i}.mif - | mrmath -axis 3 - mean ${raw}/b0s_p${dwi_i}.mif -force
        
            # read the median b0 values
            scale[dwi_i]=$(mrstats ${raw}/b0s_p1.mif -mask `dwi2mask ${raw}/dwi_p${dwi_i}.mif - -quiet` -output median)
            kul_e2cl "   dataset p${dwi_i} has ${scale[dwi_i]} as mean b0 intensity" ${preproc}/${log}

            #echo "scaling ${raw}/dwi_p${dwi_i}_scaled.mif"
            mrcalc -quiet ${scale[1]} ${scale[dwi_i]} -divide ${raw}/dwi_p${dwi_i}.mif -mult ${raw}/dwi_p${dwi_i}_scaled.mif -force

            ((dwi_i++))

        done 

        #echo "catting dwi_orig"
        mrcat ${raw}/dwi_p*_scaled.mif ${preproc}/dwi_orig.mif

    fi


else

    echo " Conversion has been done already... skipping to next step"

fi


# STEP 2 - DWI Preprocessing ---------------------------------------------
# dwidenoise
cd ${preproc}
mkdir -p dwi

# check if first 2 steps of dwi preprocessing are done 
if [ ! -f dwi/degibbs.mif ]; then

    kul_e2cl " Start part 1 of Preprocessing" ${log}

    # dwidenoise
    kul_e2cl "   dwidenoise..." ${log}
    dwidenoise dwi_orig.mif dwi/denoise.mif -noise dwi/noiselevel.mif -nthreads $ncpu -force

    # mrdegibbs
    kul_e2cl "   mrdegibbs..." ${log}
    mrdegibbs dwi/denoise.mif dwi/degibbs.mif -nthreads $ncpu -force
    rm dwi/denoise.mif

fi

# check if step 3 of dwi preprocessing is done (dwipreproc, i.e. motion and distortion correction takes very long)
if [ ! -f dwi/geomcorr.mif ]; then

    # motion and distortion correction using rpe_header
    kul_e2cl "   dwipreproc using rpe_header (this takes time!)..." ${log}

    dwiextract dwi_orig.mif -bzero - | dwiextract - -pe 0,1,0 raw/b0s_pe1.mif -force
    dwiextract dwi_orig.mif -bzero - | dwiextract - -pe 0,-1,0 raw/b0s_pe2.mif -force
    mrconvert raw/b0s_pe1.mif -coord 3 1:2 raw/b0s_pe1_first2.mif -force
    mrconvert raw/b0s_pe2.mif -coord 3 1:2 raw/b0s_pe2_first2.mif -force
    mrcat raw/b0s_pe1_first2.mif raw/b0s_pe2_first2.mif raw/se_epi_for_topup.mif -force

    dwipreproc dwi/degibbs.mif dwi/geomcorr.mif -rpe_header \
    -se_epi raw/se_epi_for_topup.mif -nthreads $ncpu -eddy_options "${eddy_options}" 
    
    #rm dwi/degibbs.mif

fi

# check if next 4 steps of dwi preprocessing are done
if [ ! -f dwi_preproced.mif ]; then

    kul_e2cl " Start part 2 of Preprocessing" ${log}

    # bias field correction
    kul_e2cl "    dwibiascorrect" ${log}
    dwibiascorrect -ants dwi/geomcorr.mif dwi/biascorr.mif -bias dwi/biasfield.mif -nthreads $ncpu -force 

    # upsample the images
    kul_e2cl "    upsampling resolution..." ${log}
    mrresize dwi/biascorr.mif -vox 1.3 dwi/upsampled.mif -nthreads $ncpu -force 
    rm dwi/biascorr.mif

    # copy to main directory for subsequent processing
    kul_e2cl "    saving..." ${log}
    mrconvert dwi/upsampled.mif dwi_preproced.mif -set_property comments "Preprocessed dMRI data." -nthreads $ncpu -force 
    rm dwi/upsampled.mif

    # create mask of the dwi data (note masking works best on low b-shells, if high b-shells are noisy)
    kul_e2cl "    creating mask of the dwi data (note masking works best on low b-shells, if high b-shells are noisy)..." ${log}
    #dwiextract -shells 0,200,500,1200 dwi_preproced.mif - | dwi2mask - dwi_mask.nii.gz -nthreads $ncpu -force 
    dwiextract -shells 0,2400 dwi_preproced.mif - | dwi2mask - dwi_mask.nii.gz -nthreads $ncpu -force
    
    # create mean b0 of the dwi data
    kul_e2cl "    creating mean b0 of the dwi data ..." ${log}
    dwiextract -quiet dwi_preproced.mif -bzero - | mrmath -axis 3 - mean dwi_b0.nii.gz -force 

    # create 2nd mask of the dwi data (using ants)
    #kul_e2cl "    create 2nd mask of the dwi data (using ants)..." ${log}
    #
    #antsBrainExtraction.sh -d 3 -a dwi_b0.nii.gz -e ../T2_template_and_tpms/mni_icbm152_t2_tal_nlin_asym_09a.nii \
    #    -m ../T2_template_and_tpms/mni_icbm152_t2_tal_nlin_asym_09a_mask.nii -o ./dwi_mask2_ -s nii.gz -u 1

else

    echo " Preprocessing already done, skipping"

fi



# STEP 3 - RESPONSE ---------------------------------------------
mkdir -p response
# response function estimation
if [ ! -f response/wm_response.txt ]; then
    kul_e2cl "   Calculating dwi2response..." ${log}
    dwi2response dhollander dwi_preproced.mif response/wm_response.txt response/gm_response.txt response/csf_response.txt -nthreads $ncpu -force 

else

    echo " dwi2response already done, skipping..."

fi

if [ ! -f response/wmfod.mif ]; then
    kul_e2cl "   Calculating dwi2fod..." ${log}
    dwi2fod msmt_csd dwi_preproced.mif response/wm_response.txt response/wmfod.mif response/gm_response.txt response/gm.mif \
        response/csf_response.txt response/csf.mif -mask dwi_mask.nii.gz -force -nthreads $ncpu 

else

    echo " dwi2fod already done, skipping..."

fi


# STEP 4 - DO QA ---------------------------------------------
# Make an FA/dec image

mkdir -p qa

if [ ! -f qa/dec.mif ]; then
    kul_e2cl "   Calculating FA/dec..." ${log}
    dwi2tensor dwi_preproced.mif dwi_dt.mif -force
    tensor2metric dwi_dt.mif -fa qa/fa.nii.gz -mask dwi_mask.nii.gz -force
    fod2dec response/wmfod.mif qa/dec.mif -force

    #mrconvert dwi/noiselevel.mif qa/noiselevel.nii.gz

fi

# STEP 5 - Anatomical Processing ---------------------------------------------
# Brain_extraction, Registration of dmri to T1, MNI Warping, freesurfer, 5tt
mkdir -p T1w
mkdir -p dwi_reg

bids_anat_search="${cwd}/${bids_subj}/anat/*_T1w.nii.gz"
echo $bids_anat_search
bids_anat=$(ls $bids_anat_search)
#echo $bids_anat

# bet the T1w using ants
if [ ! -f T1w/T1w_BrainExtractionBrain.nii.gz ]; then
    kul_e2cl " skull stripping the T1w using ants..." $log
    antsBrainExtraction.sh -d 3 -a $bids_anat -e $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz \
        -m $FSLDIR/data/standard/MNI152_T1_1mm_brain_mask.nii.gz -o T1w/T1w_ -u 1 

else

    echo " skull stripping of the T1w already done, skipping..."

fi

#mrresize -voxel 2 T1w/T1w_BrainExtractionBrain.nii.gz T1w/T1w_BrainExtractionBrain_LR.nii.gz -force
ants_b0=dwi_b0.nii.gz
ants_anat=T1w/T1w_BrainExtractionBrain.nii.gz
ants_type=dwi_reg/rigid

# register mean b0 to betted T1w (rigid)
if [ ! -f dwi_reg/rigid_outWarped.nii.gz ]; then

    kul_e2cl " registering the the dmri b0 to the betted T1w image (rigid)..." ${log}
    antsRegistration --verbose 1 --dimensionality 3 \
        --output [${ants_type}_out,${ants_type}_outWarped.nii.gz,${ants_type}_outInverseWarped.nii.gz] \
        --interpolation Linear \
        --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [$ants_anat,$ants_b0,1] \
        --transform Rigid[0.1] \
        --metric MI[$ants_anat,$ants_b0,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox

else

    echo " registering the T1w image to  (rigid) already done, skipping..."

fi

ants_type=dwi_reg/affine

# register mean b0 to betted T1w (affine)
if [ ! -f dwi_reg/affine_outWarped.nii.gz ]; then

    kul_e2cl " registering the the dmri b0 to the betted T1w image (affine)..." ${log}
    antsRegistration --verbose 1 --dimensionality 3 \
        --output [${ants_type}_out,${ants_type}_outWarped.nii.gz,${ants_type}_outInverseWarped.nii.gz] \
        --interpolation Linear \
        --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [$ants_anat,$ants_b0,1] \
        --transform Rigid[0.1] \
        --metric MI[$ants_anat,$ants_b0,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform Affine[0.1] \
        --metric MI[$ants_anat,$ants_b0,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox 

else

    echo " registering the T1w image to  (affine) already done, skipping..."

fi

# Apply the rigid transformation of the dMRI to T1 
#  to the wmfod and the preprocessed dMRI data

ConvertTransformFile 3 dwi_reg/rigid_out0GenericAffine.mat dwi_reg/rigid_out0GenericAffine.txt

transformconvert dwi_reg/rigid_out0GenericAffine.txt itk_import \
    dwi_reg/rigid_out0GenericAffine_mrtrix.txt -force

mrtransform dwi_preproced.mif -linear dwi_reg/rigid_out0GenericAffine_mrtrix.txt \
    dwi_preproced_reg2T1w.mif -nthreads $ncpu -force 
mrtransform response/wmfod.mif -linear dwi_reg/rigid_out0GenericAffine_mrtrix.txt \
    response/wmfod_reg2T1w.mif -nthreads $ncpu -force 



# DO QA ---------------------------------------------
# Make an FA/dec image


if [ ! -f qa/dec_reg2T1w.mif ]; then

    kul_e2cl "   Calculating FA/dec..." ${log}
    dwi2tensor dwi_preproced_reg2T1w.mif dwi_dt_reg2T1w.mif -force
    tensor2metric dwi_dt_reg2T1w.mif -fa qa/fa_reg2T1w.nii.gz -mask dwi_mask.nii.gz -force
    fod2dec response/wmfod_reg2T1w.mif qa/dec_reg2T1w.mif -force
    fod2dec response/wmfod_reg2T1w.mif qa/dec_reg2T1w_on_t1w.mif -contrast $ants_anat -force

fi

kul_e2cl "Finished " ${log}
exit 1

OLD-CODE:
ants_type=dwi_reg/syn
syn_convergence=1e-2

# register mean b0 to betted T1w (minimal application of SyN)
if [ ! -f dwi_reg/syn_outWarped.nii.gz ]; then

    kul_e2cl " registering the the dmri b0 to the betted T1w image (minimal application of SyN)..." ${log}
    antsRegistration --verbose 1 --dimensionality 3 \
        --output [${ants_type}_out,${ants_type}_outWarped.nii.gz,${ants_type}_outInverseWarped.nii.gz] \
        --interpolation Linear \
        --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [$ants_anat,$ants_b0,1] \
        --transform Rigid[0.1] \
        --metric MI[$ants_anat,$ants_b0,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform Affine[0.1] \
        --metric MI[$ants_anat,$ants_b0,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform SyN[0.1,3,0] \
        --metric CC[$ants_anat,$ants_b0,1,4] --convergence [100x70x50x20,1e-2,2] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox

else

    echo " registering the T1w image to  (rigid) already done, skipping..."

fi

warpinit $mrtransform_file dwi_reg/identity_warp[].nii -force

for i in {0..2}; do
  WarpImageMultiTransform 3 dwi_reg/identity_warp${i}.nii dwi_reg/mrtrix_warp${i}.nii \
    -R $ants_anat dwi_reg/syn_out1Warp.nii.gz dwi_reg/syn_out0GenericAffine.mat;
done;

warpcorrect dwi_reg/mrtrix_warp[].nii dwi_reg/mrtrix_warp_corrected.mif -nthreads $ncpu -force 

