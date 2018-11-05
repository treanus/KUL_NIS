#!/bin/bash -e
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#  Developed for Segmentation of the Dentato-rubro-thalamic tract in thalamus for DBS target selection
#   following the paper "Connectivity derived thalamic segmentation in deep brain stimulation for tremor"
#       of Akram et al. 2018 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5790021/pdf/main.pdf)
#  Project PI's: Stefan Sunaert & Bart Nuttin
#
# Requires Mrtrix3, FSL, ants, freesurfer
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 11/10/2018 - alpha version
v="v0.1 - dd 11/10/2018"

# To Do
#  - start from BIDS format?
#  - check requirements
#  - register dwi to T1 with ants-syn
#  - align b0s for calc of mean b0 and for topup
#  - fod calc msmt-5tt in stead of dhollander
#  - use HPC of KUL?
#  - how to import in neuronavigation?
#  - warp the resulted TH-* back into MNI space for group analysis 


# A few fixed (for now) parameters:

    # Specify additional options for FSL eddy
    eddy_options="--data_is_shelled --slm=linear --repol "

    # Number of desired streamlines
    nods=2000

    # tmp directory for temporary processing
    tmp=/tmp
# 


# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl (for logging)
#  - dcmtags (for reading specific parameters from dicom header)
#
# this script uses "preprocessing control", i.e. if some steps are already processed it will skip these

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh

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

 
Semi-Optional arguments (to be used on first run, but can be ommited if conversion of dicoms is already done):

     -d:  dicom_zip_file (the zip or tar.gz containing all your dicoms)

Optional arguments:

     -p:  number of cpu for parallelisation
     -b:  use BIDS format (needs to be there already)
     -v:  show output from mrtrix commands


USAGE

    exit 1
}

#dcmtags (for reading specific parameters from dicom header)
function dcmtags {

    # for calc of ess/trt are needed: FieldStrength (0018 0087), WaterFatShift, EPIFactor
    
    local dcm_file=$1 
    local out=${preproc}/log/${subj}_dcm_tags.txt

    local seriesdescr=$(dcminfo "$dcm_file" -tag 0008 103E | awk '{print $(NF)}')
    local waterfatshift=$(dcminfo "$dcm_file" -tag 2001 1022 | awk '{print $(NF)}')
    local fieldstrength=$(dcminfo "$dcm_file" -tag 0018 0087 | awk '{print $(NF)}')
    local epifactor=$(dcminfo "$dcm_file" -tag 2001 1013 | awk '{print $(NF)}')
    local water_fat_diff_ppm=3.3995
    local resonance_freq_mhz_tesla=42.576

    local water_fat_shift_hz=$(echo $fieldstrength $water_fat_diff_ppm $resonance_freq_mhz_tesla echo $fieldstrength $water_fat_diff_ppm | awk '{print $1 * $2 * $3}')
    
    #effective_echo_spacing_msec  = 1000 * WFS_PIXEL/(water_fat_shift_hz * (EPI_FACTOR + 1))
    ees_sec=$(echo $waterfatshift $water_fat_shift_hz $epifactor | awk '{print $1 / ($2 * ($3 + 1))}')
    
    #total_readout_time_fsl_msec      = EPI_FACTOR * effective_echo_spacing_msec;
    trt_sec=$(echo $epifactor $ees_sec | awk '{print $1 * $2 }')
    
    if [ ! -f $out ]; then
        echo -e "series \t field \t epi \t wfs \t ees_sec \t trt_sec " > $out
    fi
    echo -e "$seriesdescr \t $fieldstrength \t $epifactor \t $waterfatshift \t $ees_sec \t $trt_sec " >> $out

}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
ncpu=6
silent=1

# Set required options
s_flag=0
dcm_flag=0
bids_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "s:p:d:bvh" OPT; do

        case $OPT in
        s) #subject
            s_flag=1
            subj=$OPTARG
        ;;
        p) #parallel
            ncpu=$OPTARG
        ;;
        d) #dicom_zip_file
            dcm_flag=1
            dcm=$OPTARG
        ;;
        b) #Use BIDS format
            bids_flag=1
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

if [ $dcm_flag -eq 0 && $bids_flag -eq 0 ] ; then 

    if [ ! -f ${preproc}/dwi_orig.mif ]; then

        echo
        echo "For the dcm2nii conversion step provide the -d dicom_zip_file option please" >&2
        echo
        exit 2
    
    fi
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
preproc=${subj}_preproc

# Directory to put raw mif data in
raw=${preproc}/raw

# set up preprocessing & logdirectory
mkdir -p ${preproc}/raw
mkdir -p ${preproc}/log

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=log/log_${d}.txt



# SAY HELLO ---

kul_e2cl "Welcome to KUL_dwi_preproc $v - $d" ${preproc}/${log}


# STEP 1 - CONVERSION ---------------------------------------------
# test if conversion has been done
if [ ! -f ${preproc}/dwi_orig.mif ]; then

    if [ $bids_flag -eq 0 ]; then

    kul_e2cl " Converting datasets..." ${preproc}/${log}

    # take the content of the zip file with dicoms
    tar -ztvf ${dcm} > ${preproc}/log/${subj}_dcm_content.txt

    # find the part1 in the zip file with dicoms and extract & mrconvert it
    part1=$(grep part1 ${preproc}/log/${subj}_dcm_content.txt | head -n 1 | cut -c 49-)
    part1_file=$(grep part1 ${preproc}/log/${subj}_dcm_content.txt | awk 'NR==2' | cut -c 49-)
    echo "  part1 dicoms are in ${part1}"
    tar -C ${tmp} -xzf ${dcm} "${part1}"
    dcmtags "${tmp}/${part1_file}"
    trt_part1=$trt_sec
    echo 0 | mrconvert "${tmp}/${part1}" ${raw}/dwi_p1.mif -force -clear_property comments -nthreads $ncpu

    # check the number of directions in the dwi (part1)
    n_d_p1=$(mrinfo -size ${raw}/dwi_p1.mif | awk '{print $NF}')
    kul_e2cl "  sequence part 1 has ${n_d_p1} directions" ${preproc}/${log}

    # find the part2 in the zip file with dicoms and extract & mrconvert it
    part2=$(grep part2 ${preproc}/log/${subj}_dcm_content.txt | head -n 1 | cut -c 49-)
    part2_file=$(grep part2 ${preproc}/log/${subj}_dcm_content.txt | awk 'NR==2' | cut -c 49-)
    echo "  part2 dicoms are in ${part2}"
    tar -C ${tmp} -xzf ${dcm} "${part2}"
    dcmtags "${tmp}/${part2_file}"
    trt_part2=$trt_sec
    echo 1 | mrconvert "${tmp}/${part2}" ${raw}/dwi_p2.mif -force -clear_property comments -nthreads $ncpu

    # check the number of directions in the dwi (part2)
    n_d_p2=$(mrinfo -size ${raw}/dwi_p2.mif | awk '{print $NF}')
    kul_e2cl "  sequence part 2 has ${n_d_p2} directions" ${preproc}/${log}

    # find the revphase in the zip file with dicoms and extract & mrconvert it (if any)
    part3=$(grep "revphase" ${preproc}/log/${subj}_dcm_content.txt | head -n 1 | cut -c 49-)
    if [ "${part3}" == "" ]; then
        kul_e2cl "  there is no revphase" ${preproc}/${log}
    else
        echo "  revphase dicoms are in ${part3}"
        tar -C ${tmp} -xzf ${dcm} "${part3}"
        part3_file=$(grep "revphase" ${preproc}/log/${subj}_dcm_content.txt | awk 'NR==2' | cut -c 49-)
        dcmtags "${tmp}/${part3_file}"
        trt_part3=$trt_sec
        echo 2 | mrconvert "${tmp}/${part3}" ${raw}/dwi_p3.mif -force -clear_property comments -nthreads $ncpu

        n_d_p3=$(mrinfo -size ${raw}/dwi_p3.mif | awk '{print $NF}')
        kul_e2cl "  sequence rev_phase has ${n_d_p3} directions" ${preproc}/${log}
    fi

    # setup FSL/Mrtrix pe_info
    if [ "${part3}" == "" ]; then
        for i in `seq $n_d_p1`; do echo "0 1 0 $trt_part1"; done > ${raw}/pe_p1_1.txt
        for i in `seq $n_d_p1`; do printf '1 %.0s'; done > ${raw}/pe_p1_2.txt
        for i in `seq $n_d_p2`; do echo "0 -1 0 $trt_part2"; done > ${raw}/pe_p2_1.txt
        for i in `seq $n_d_p2`; do printf '1 %.0s'; done > ${raw}/pe_p2_2.txt
    else
        for i in `seq $n_d_p1`; do echo "0 1 0 $trt_part1"; done > ${raw}/pe_p1_1.txt
        for i in `seq $n_d_p1`; do printf '1 %.0s'; done > ${raw}/pe_p1_2.txt
        for i in `seq $n_d_p2`; do echo "0 1 0 $trt_part2"; done > ${raw}/pe_p2_1.txt
        for i in `seq $n_d_p2`; do printf '1 %.0s'; done > ${raw}/pe_p2_2.txt
        for i in `seq $n_d_p3`; do echo "0 -1 0 $trt_part3"; done > ${raw}/pe_p3_1.txt
        for i in `seq $n_d_p3`; do printf '1 %.0s'; done > ${raw}/pe_p3_2.txt
    fi

    # convert raw dwi data (second pass)
    echo 0 | mrconvert "${tmp}/${part1}" ${raw}/dwi_p1.mif -strides 1:3 -force -import_pe_eddy ${raw}/pe_p1_1.txt ${raw}/pe_p1_2.txt -nthreads $ncpu -clear_property comments -nthreads $ncpu
    rm -rf "${tmp}/${part1}"
    echo 1 | mrconvert "${tmp}/${part2}" ${raw}/dwi_p2.mif -strides 1:3 -force -import_pe_eddy ${raw}/pe_p2_1.txt ${raw}/pe_p2_2.txt -nthreads $ncpu -clear_property comments -nthreads $ncpu
    rm -rf "${tmp}/${part2}"
    if [ ! "${part3}" == "" ]; then
        echo 2 | mrconvert "${tmp}/${part3}" ${raw}/dwi_p3.mif -strides 1:3 -force -import_pe_eddy ${raw}/pe_p3_1.txt ${raw}/pe_p3_2.txt -nthreads $ncpu -clear_property comments -nthreads $ncpu
        rm -rf "${tmp}/${part3}"
    fi


    # convert raw T1w data, using -strides 1:3 to get orientation correct for FSL
    t1=$(grep "MPRAGE\|3DTFE" ${preproc}/log/${subj}_dcm_content.txt | head -n 1 | cut -c 49-) #'MPRAGE\|3DTFE' searches for either MPRAGE or 3DTFE
    t1_file=$(grep "MPRAGE\|3DTFE" ${preproc}/log/${subj}_dcm_content.txt | awk 'NR==2' | cut -c 49-)
    echo "  3D-T1w dicoms are in ${t1}"
    tar -C ${tmp} -xzf ${dcm} "${t1}"
    dcmtags "${tmp}/${t1_file}"
    echo 3 | mrconvert "${tmp}/${t1}" ${raw}/T1w.nii.gz -strides 1:3 -force -nthreads $ncpu -clear_property comments -nthreads $ncpu
    rm -rf "${tmp}/${t1}"

    # convert raw FLAIR data, using -strides 1:3 to get orientation correct for FSL
    fl=$(grep "FLAIR" ${preproc}/log/${subj}_dcm_content.txt | head -n 1 | cut -c 49-)
    fl_file=$(grep "FLAIR" ${preproc}/log/${subj}_dcm_content.txt | awk 'NR==2' | cut -c 49-)
    echo "  3D-FLAIR dicoms are in ${fl}"
    tar -C ${tmp} -xzf ${dcm} "${fl}"
    dcmtags "${tmp}/${fl_file}"
    echo 4 | mrconvert "${tmp}/${fl}" ${raw}/FLAIR.nii.gz -strides 1:3 -force -nthreads $ncpu -clear_property comments -nthreads $ncpu
    rm -rf "${tmp}/${fl}"

    # intensity matching
    kul_e2cl "  adjusting scaling factor of the native dwi datasets" ${preproc}/${log}

    if [ "${part3}" == "" ]; then
        mrcat -quiet ${raw}/dwi_p1.mif ${raw}/dwi_p2.mif - | mrstats - -mask `dwi2mask ${raw}/dwi_p1.mif - -quiet` -output median > ${raw}/median.txt  
    else
        mrcat -quiet ${raw}/dwi_p1.mif ${raw}/dwi_p2.mif ${raw}/dwi_p3.mif - | mrstats - -mask `dwi2mask ${raw}/dwi_p1.mif - -quiet` -output median > ${raw}/median.txt  
    fi

    dwiextract -quiet -bzero ${raw}/dwi_p1.mif - | mrmath -axis 3 - mean ${raw}/b0s_p1.mif -force
    dwiextract -quiet -bzero ${raw}/dwi_p2.mif - | mrmath -axis 3 - mean ${raw}/b0s_p2.mif -force
    if [ ! "${part3}" == "" ]; then
        dwiextract -quiet -bzero ${raw}/dwi_p3.mif - | mrmath -axis 3 - mean ${raw}/b0s_p3.mif -force
    fi

    # read the median b0 values
    scale0=$(mrstats ${raw}/b0s_p1.mif -mask `dwi2mask ${raw}/dwi_p1.mif - -quiet` -output median)
    kul_e2cl "   dataset p1 has $scale0 as mean b0 intensity$" ${preproc}/${log}
    scale1=$(mrstats ${raw}/b0s_p2.mif -mask `dwi2mask ${raw}/dwi_p2.mif - -quiet` -output median)
    kul_e2cl "   dataset p2 has $scale1 as mean b0 intensity" ${preproc}/${log}
    if [ ! "${part3}" == "" ]; then
        scale2=$(mrstats ${raw}/b0s_p3.mif -mask `dwi2mask ${raw}/dwi_p3.mif - -quiet` -output median)
        kul_e2cl "   dataset p3 has $scale2 as mean b0 intensity" ${preproc}/${log}
    fi

    # adjust the scale accroding to the median b0 values
    if [ "${part3}" == "" ]; then
        mrcalc -quiet $scale0 $scale1 -divide ${raw}/dwi_p2.mif -mult - | mrcat ${raw}/dwi_p1.mif - ${preproc}/dwi_orig.mif 
    else
        mrcalc -quiet $scale0 $scale1 -divide ${raw}/dwi_p2.mif -mult - | mrcat ${raw}/dwi_p1.mif - ${raw}/dwi_p1p2.mif 
        mrcalc -quiet $scale0 $scale2 -divide ${raw}/dwi_p3.mif -mult - | mrcat ${raw}/dwi_p1p2.mif - ${preproc}/dwi_orig.mif 
    fi

    else

    # BIDS
    kul_e2cl " Preparing datasets from BIDS directory..." ${preproc}/${log}
    
    bids_subj=BIDS/"sub-$subj"/ses-tp1/

    # convert raw T1w data, using -strides 1:3 to get orientation correct for FSL
    bids_t1w="$bids_subj/anat/sub-${subj}_T1w.nii.gz"
    #mrconvert "$bids_t1w" ${raw}/T1w.nii.gz -strides 1:3 -force -nthreads $ncpu -clear_property comments -nthreads $ncpu

    # convert raw FLAIR data, using -strides 1:3 to get orientation correct for FSL
    bids_flair="$bids_subj/anat/sub-${subj}_FLAIR.nii.gz"
    #mrconvert "$bids_flair" ${raw}/FLAIR.nii.gz -strides 1:3 -force -nthreads $ncpu -clear_property comments -nthreads $ncpu

    # convert dwi
    bids_dwi_search="$bids_subj/dwi/sub-*_dwi.nii.gz"

    bids_dwi_found=$(ls $bids_dwi_search)
    echo $bids_dwi_found
    
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
            echo $dwi_base
        
            mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
            -json_import ${dwi_base}.json ${raw}/dwi_p${dwi_i}.mif -strides 1:3 -force -clear_property comments -nthreads $ncpu

            dwiextract -quiet -bzero ${raw}/dwi_p${dwi_i}.mif - | mrmath -axis 3 - mean ${raw}/b0s_p${dwi_i}.mif -force
        
            # read the median b0 values
            scale[dwi_i]=$(mrstats ${raw}/b0s_p1.mif -mask `dwi2mask ${raw}/dwi_p${dwi_i}.mif - -quiet` -output median)
            kul_e2cl "   dataset p${dwi_i} has ${scale[dwi_i]} as mean b0 intensity$" ${preproc}/${log}

            echo "scaling ${raw}/dwi_p${dwi_i}_scaled.mif"
            mrcalc -quiet ${scale[1]} ${scale[dwi_i]} -divide ${raw}/dwi_p${dwi_i}.mif -mult ${raw}/dwi_p${dwi_i}_scaled.mif -force

            ((dwi_i++))

        done 

        echo "catting dwi_orig"
        mrcat ${raw}/dwi_p*_scaled.mif ${preproc}/dwi_orig.mif

    fi

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
    dwipreproc dwi/degibbs.mif dwi/geomcorr.mif -rpe_header -nthreads $ncpu -eddy_options "${eddy_options}" 
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
    mrresize dwi/biascorr.mif -vox 1 dwi/upsampled.mif -nthreads $ncpu -force 
    rm dwi/biascorr.mif

    # copy to main directory for subsequent processing
    kul_e2cl "    saving..." ${log}
    mrconvert dwi/upsampled.mif dwi_preproced.mif -set_property comments "Preprocessed dMRI data." -nthreads $ncpu -force 
    rm dwi/upsampled.mif

    # create mask of the dwi data (note masking works best on low b-shells, if high b-shells are noisy)
    kul_e2cl "    create mask of the dwi data (note masking works best on low b-shells, if high b-shells are noisy)..." ${log}
    #dwiextract -shells 0,200,500,1200 dwi_preproced.mif - | dwi2mask - dwi_mask.nii.gz -nthreads $ncpu -force 
    dwiextract -shells 0,2400 dwi_preproced.mif - | dwi2mask - dwi_mask.nii.gz -nthreads $ncpu -force

    # create 2nd mask of the dwi data (using ants)
    #kul_e2cl "    create 2nd mask of the dwi data (using ants)..." ${log}
    #dwiextract -quiet dwi_preproced.mif -bzero - | mrmath -axis 3 - mean dwi_b0.nii.gz -force 
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

# QA - make an FA/dec image
mkdir -p qa

if [ ! -f qa/dec.mif ]; then
    kul_e2cl "   Calculating FA/dec..." ${log}
    dwi2tensor dwi_preproced.mif dwi_dt.mif -force
    tensor2metric dwi_dt.mif -fa qa/fa.nii.gz -mask dwi_mask.nii.gz -force
    fod2dec response/wmfod.mif qa/dec.mif -force
    #fod2dec response/wmfod.mif qa/dec_t1w.mif -contrast T1w_brain_reg2_b0_HR.nii.gz -force

    mrconvert dwi/noiselevel.mif qa/noiselevel.nii.gz

fi


# STEP 4 - Anatomical Processing ---------------------------------------------
# Brain_extraction, Registration of dmri to T1, MNI Warping, freesurfer, 5tt
mkdir -p T1w
mkdir -p dwi_reg

# bet the T1w using ants
if [ ! -f T1w/T1w_BrainExtractionBrain.nii.gz ]; then
    kul_e2cl " Skull stripping the T1w using ants..." ${log}
    antsBrainExtraction.sh -d 3 -a ../${raw}/T1w.nii.gz -e $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz -m $FSLDIR/data/standard/MNI152_T1_1mm_brain_mask.nii.gz -o T1w/T1w_ -u 1 

else

    echo " Skull stripping of the T1w already done, skipping..."

fi





if 0; then 
# register betted T1w to mean b0 (rigid)
if [ ! -f T1w/T1w_brain_reg2_b0_deformed.nii.gz ]; then

    kul_e2cl " registering the T1w image to the dmri b0 (rigid)..." ${log}
    antsaffine.sh 3 dwi_b0.nii.gz T1w/T1w_BrainExtractionBrain.nii.gz T1w/T1w_brain_reg2_b0_ PURELY-RIGID
    
else

    echo " registering the T1w image to the dmri b0 (rigid) already done, skipping..."

fi

# register T1w to b0 (rigid), writing results
if [ ! -f T1w_brain_reg2_b0_HR.nii.gz ]; then
    kul_e2cl " registering the T1w image to the dmri b0 (rigid), writing data..." ${log}
    WarpImageMultiTransform 3 T1w/T1w_BrainExtractionBrain.nii.gz   T1w_brain_reg2_b0_HR.nii.gz -R T1w/T1w_BrainExtractionBrain.nii.gz  T1w/T1w_brain_reg2_b0_Affine.txt --use-BSpline
    WarpImageMultiTransform 3 ../${raw}/T1w.nii.gz                  T1w/T1w_full_reg2_b0_HR.nii -R ../${raw}/T1w.nii.gz                 T1w/T1w_brain_reg2_b0_Affine.txt --use-BSpline
fi


# registering the dmri fa (using ants syn) to the T1w
if [ ! -f dwi_reg/fa_reg2_T1w_brain_bj1_Warped.nii.gz ]; then

    kul_e2cl " registering the dmri fa (using ants syn) to the T1w..." ${log}
    antsRegistrationSyNQuick.sh -d 3 -m qa/fa.nii.gz -f T1w/T1w_brain_reg2_b0_deformed.nii.gz -x dwi_mask2.nii.gz T1w/-o dwi_reg/fa_reg2_T1w_brain_std_ -t s -n $ncpu 

    antsRegistrationSyNQuick.sh -d 3 -m qa/fa.nii.gz -f T1w/T1w_brain_reg2_b0_deformed.nii.gz -x dwi_mask2.nii.gz -o dwi_reg/fa_reg2_T1w_brain_stdj1_ -t s -j 1 -n $ncpu

    antsRegistrationSyNQuick.sh -d 3 -m qa/fa.nii.gz -f T1w/T1w_brain_reg2_b0_deformed.nii.gz -x dwi_mask2.nii.gz -o dwi_reg/fa_reg2_T1w_brain_rigid_ -t r -n $ncpu

    antsRegistrationSyNQuick.sh -d 3 -m qa/fa.nii.gz -f T1w/T1w_brain_reg2_b0_deformed.nii.gz -x dwi_mask2.nii.gz -o dwi_reg/fa_reg2_T1w_brain_b_ -t b -n $ncpu

    antsRegistrationSyNQuick.sh -d 3 -m qa/fa.nii.gz -f T1w/T1w_brain_reg2_b0_deformed.nii.gz -x dwi_mask2.nii.gz -o dwi_reg/fa_reg2_T1w_brain_bj1_ -t b -j 1-n $ncpu

else

    echo " registering of fa dMRI (using ants syn) to the T1w already done, skipping..."

fi


# warp the T1w to MNI space using ants
if [ ! -f T1w/T1w_brain_reg2_b0_MNI_Warped.nii.gz ]; then

    kul_e2cl " warping the T1w to MNI (1 mm) space using ants..." ${log}
    antsRegistrationSyN.sh -d 3 -f $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -m T1w_brain_reg2_b0_HR.nii.gz -o T1w/T1w_brain_reg2_b0_MNI_ -t s -j 1 -n $ncpu 

else

    echo " warping the T1w to MNI (1 mm) space already done, skipping..."

fi
fi

# freesurfer
mkdir -p freesurfer
if [ ! -f freesurfer/${subj}/scripts/recon-all.done ]; then

    kul_e2cl " performing freesurfer recon-all..." ${log}
    SUBJECTS_DIR=$(pwd)/freesurfer
    export SUBJECTS_DIR
    recon-all -subject $subj -i T1w/T1w_full_reg2_b0_HR.nii -all -openmp ${ncpu} -parallel 

else

    echo " freesurfer already done, skipping..."


fi

# 5tt segmentation & tracking
mkdir -p 5tt
if [ ! -f 5tt/5tt2gmwmi.nii.gz ]; then

    kul_e2cl " Performig 5tt..." ${log}
    5ttgen fsl T1w_brain_reg2_b0_HR.nii.gz 5tt/5ttseg.mif -premasked -nocrop -force -nthreads $ncpu #
    5ttcheck -masks 5tt/failed_5tt 5tt/5ttseg.mif -force -nthreads $ncpu #
    5tt2gmwmi 5tt/5ttseg.mif 5tt/5tt2gmwmi.nii.gz -force #

else

    echo " 5tt already done, skipping..."

fi


# STEP 5 - ROI Processing ---------------------------------------------
mkdir -p roi

# Warp the MNI ROIS into subject space (apply INVERSE warp using ants)
if [ ! -f atlas/TH-SMA_R.nii.gz ]; then
kul_e2cl " Warping the MNI ROIS into subjects space..." ${log}
WarpImageMultiTransform 3 ../ROIS/ROI_DENTATE_L.nii.gz roi/DENTATE_L.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/ROI_DENTATE_R.nii.gz roi/DENTATE_R.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/ROI_THALAMUS_L.nii.gz roi/THALAMUS_L.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/ROI_THALAMUS_R.nii.gz roi/THALAMUS_R.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/ROI_M1.nii.gz roi/M1_full.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz \
    -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/ROI_S1.nii.gz roi/S1_full.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz \
    -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/ROI_SMA_and_PMC.nii.gz roi/SMA_and_PMC_full.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz \
    -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 

# transect the S1, M1 and SMA_and_PMC ROIS with 5ttgen wm/gm interface
kul_e2cl " Intersecting ROIS with 5tt WM/GM..." ${log}
WarpImageMultiTransform 3 5tt/5tt2gmwmi.nii.gz 5tt/5tt2gmwmi_dwi.nii.gz -R roi/M1_full.nii.gz --reslice-by-header
fslmaths roi/M1_full.nii.gz -mas 5tt/5tt2gmwmi_dwi.nii.gz roi/M1.nii.gz
fslmaths roi/S1_full.nii.gz -mas 5tt/5tt2gmwmi_dwi.nii.gz roi/S1.nii.gz
fslmaths roi/SMA_and_PMC_full.nii.gz -mas 5tt/5tt2gmwmi_dwi.nii.gz roi/SMA_and_PMC.nii.gz

# Warp the Atlas ROIS into subjects space (apply INVERSE warp using ants)
mkdir -p atlas
kul_e2cl " Warping the Atlas ROIS into subjects space..." ${log}
WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/lh/Dentate.nii.gz atlas/TH-Dentate_L.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/lh/M1.nii.gz atlas/TH-M1_L.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/lh/S1.nii.gz atlas/TH-S1_L.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/lh/SMA.nii.gz atlas/TH-SMA_L.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 

WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/rh/Dentate.nii.gz atlas/TH-Dentate_R.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/rh/M1.nii.gz atlas/TH-M1_R.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/rh/S1.nii.gz atlas/TH-S1_R.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 
WarpImageMultiTransform 3 ../ROIS/Thalamic_DBS_Connectivity_atlas_Akram_2018/rh/SMA.nii.gz atlas/TH-SMA_R.nii.gz -R T1w/T1w_brain_reg2_b0_deformed.nii.gz -i T1w/T1w_brain_reg2_b0_MNI_0GenericAffine.mat T1w/T1w_brain_reg2_b0_MNI_1InverseWarp.nii.gz 

else

echo " Reverse warping of rois/atlas has been done already, skipping" 

fi

# freesurfer determined rois
if [ ! -f roi/WM_fs_R.nii.gz ]; then
kul_e2cl " Making the Freesurfer ROIS from subject space..." ${log}
#  project labels back into native T1w space
#mri_convert -rl freesurfer/$subj/mri/rawavg.mgz -rt nearest freesurfer/$subj/mri/aparc+aseg.mgz labels_from_FS.nii.gz
mri_convert -rl T1w/T1w_brain_reg2_b0_deformed.nii.gz -rt nearest freesurfer/$subj/mri/aparc+aseg.mgz labels_from_FS.nii.gz


# Extract relevant labels
# M1_R is 2024
fslmaths labels_from_FS -thr 2024 -uthr 2024 -bin roi/M1_fs_R
# M1_L is 1024
fslmaths labels_from_FS -thr 1024 -uthr 1024 -bin roi/M1_fs_L
# S1_R is 2022
fslmaths labels_from_FS -thr 2022 -uthr 2022 -bin roi/S1_fs_R
# S1_L is 1024
fslmaths labels_from_FS -thr 1022 -uthr 1022 -bin roi/S1_fs_L
# Thalamus_R is 49
fslmaths labels_from_FS -thr 49 -uthr 49 -bin roi/THALAMUS_fs_R
# Thalamus_L is 10
fslmaths labels_from_FS -thr 10 -uthr 10 -bin roi/THALAMUS_fs_L
# SMA_and_PMC_L are
# 1003    ctx-lh-caudalmiddlefrontal
# 1028    ctx-lh-superiorfrontal
fslmaths labels_from_FS -thr 1003 -uthr 1003 -bin roi/MFG_fs_R
fslmaths labels_from_FS -thr 1028 -uthr 1028 -bin roi/SFG_fs_R
fslmaths roi/MFG_fs_R -add roi/SFG_fs_R -bin roi/SMA_and_PMC_fs_R
# SMA_and_PMC_L are
# 2003    ctx-lh-caudalmiddlefrontal
# 2028    ctx-lh-superiorfrontal
fslmaths labels_from_FS -thr 2003 -uthr 2003 -bin roi/MFG_fs_L
fslmaths labels_from_FS -thr 2028 -uthr 2028 -bin roi/SFG_fs_L
fslmaths roi/MFG_fs_L -add roi/SFG_fs_L -bin roi/SMA_and_PMC_fs_L
# 41  Right-Cerebral-White-Matter
fslmaths labels_from_FS -thr 41 -uthr 41 -bin roi/WM_fs_R
# 2   Left-Cerebral-White-Matter
fslmaths labels_from_FS -thr 2 -uthr 2 -bin roi/WM_fs_L


else

echo " Making the Freesurfer ROIS has been done already, skipping" 

fi


# STEP 5 - Tractography Processing ---------------------------------------------

function kul_mrtrix_tracto_drt {

    for a in iFOD2 Tensor_Prob; do

        if [ ! -f ${tract}_${a}.nii.gz ]; then

            mkdir -p tracts_${a}

            # make the intersect string (this is the first of the seeds)
            intersect=${seeds%% *}

            kul_e2cl " Calculating $a ${tract} tract (all seeds with -select $nods, intersect with $intersect)" ${log} 

            # make the seed string
            local s=$(printf " -seed_image roi/%s.nii.gz"  "${seeds[@]}")
    
            # make the include string (which is same rois as seed)
            local i=$(printf " -include roi/%s.nii.gz"  "${seeds[@]}")

            # make the exclude string (which is same rois as seed)
            local e=$(printf " -exclude roi/%s.nii.gz"  "${exclude[@]}")

            # make the mask string 
            local m="-mask dwi_mask.nii.gz"

            if [ "${a}" == "iFOD2" ]; then

                # perform IFOD2 tckgen
                tckgen response/wmfod.mif tracts_${a}/${tract}.tck -algorithm $a -select $nods $s $i $e $m -nthreads $ncpu -force

            else

                # perform Tensor_Prob tckgen
                tckgen dwi_preproced.mif tracts_${a}/${tract}.tck -algorithm $a -cutoff 0.01 -select $nods $s $i $e $m -nthreads $ncpu -force

            fi

            # convert the tck in nii
            tckmap tracts_${a}/${tract}.tck tracts_${a}/${tract}.nii.gz -template T1w/T1w_brain_reg2_b0_deformed.nii.gz -force 

            # intersect the nii tract image with the thalamic roi
            fslmaths tracts_${a}/${tract}.nii -mas roi/${intersect}.nii.gz tracts_${a}/${tract}_masked
    
            # make a probabilistic image
            local m=$(mrstats -quiet tracts_${a}/${tract}_masked.nii.gz -output max)
            fslmaths tracts_${a}/${tract}_masked -div $m ${tract}_${a}

        fi
    
    done

}



# M1-Thalamic tracts
tract="TH-M1_R"
seeds=("THALAMUS_R" "M1")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-M1_L"
seeds=("THALAMUS_L" "M1")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

# M1_fs-Thalamic tracts
tract="TH-M1_fs_R"
seeds=("THALAMUS_fs_R" "M1_fs_R")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-M1_fs_L"
seeds=("THALAMUS_fs_L" "M1_fs_L")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

# S1-Thalamic tracts
tract="TH-S1_R"
seeds=("THALAMUS_R" "S1")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-S1_L"
seeds=("THALAMUS_L" "S1")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

# S1_fs-Thalamic tracts
tract="TH-S1_fs_R"
seeds=("THALAMUS_fs_R" "S1_fs_R")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-S1_fs_L"
seeds=("THALAMUS_fs_L" "S1_fs_L")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

# SMA_and_PMC-Thalamic tracts
tract="TH-SMA_and_PMC_R"
seeds=("THALAMUS_R" "SMA_and_PMC")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-SMA_and_PMC_L"
seeds=("THALAMUS_L" "SMA_and_PMC")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt  

# SMA_and_PMC-Thalamic tracts
tract="TH-SMA_and_PMC_R"
seeds=("THALAMUS_fs_R" "SMA_and_PMC_fs_L")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-SMA_and_PMC_L"
seeds=("THALAMUS_fs_L" "SMA_and_PMC_fs_L")
kul_mrtrix_tracto_drt  

# Dentato-Rubro_Thalamic tracts
tract="TH-DR_R"
seeds=("THALAMUS_R" "M1" "DENTATE_L")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-DR_L"
seeds=("THALAMUS_L" "M1" "DENTATE_R")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

kul_e2cl "KUL_dwi_preproc $v - finished processing" ${log}