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
v="v0.2 - dd 19/10/2018"

# To Do
#  - use 5ttgen with freesurfer
#  - register dwi to T1 with ants-syn
#  - fod calc msmt-5tt in stead of dhollander
#  - use HPC of KUL?
#  - how to import in neuronavigation?
#  - warp the resulted TH-* back into MNI space for group analysis 


# A few fixed (for now) parameters:

    # Number of desired streamlines
    # nods=2000
    # this has become a command line optional parameter

    # Maximum angle between successive steps for iFOD2
    theta=35

    # sift1 filtering
    # termination ratio - defined as the ratio between reduction in cost
    # function, and reduction in density of streamlines.
    # Smaller values result in more streamlines being filtered out.
    do_sift_th=10000 # when to do sift? (if more than 5000 streamlines in tract e.g.)
    term_ratio=0.5 # reduce by e.g. 50%

    # tmp directory for temporary processing
    tmp=/tmp

    # development for Donatienne
    Donatienne=0
# 


# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl (for logging)
#  - dcmtags (for reading specific parameters from dicom header)
#
# this script uses "preprocessing control", i.e. if some steps are already processed it will skip these

kul_main_dir=`dirname "$0"`
script=`basename "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs dMRI segmentation of the Dentato-rubro-thalamic tract in thalamus for DBS target selection.

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001 -n 6 

Required arguments:

     -p:  participant (anonymised name of the subject)

Optional arguments:

     -s:  session (of the participant)
     -o:  number of desired streamlines to select in tckgen (default nods=2000)
     -n:  number of cpu for parallelisation
     -v:  show output from mrtrix commands


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
nods=2000
ncpu=6
silent=1

# Set required options
p_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:o:n:vh" OPT; do

        case $OPT in
        p) #subject
            p_flag=1
            subj=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        o) #nods
            nods=$OPTARG
            #remove leading/trailing spaces
            #awk '{$nods=$nods;print}'
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
    echo "Option -s is required: give the anonymised name of a subject." >&2
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
log=log/log_${d}.txt


# --- MAIN ----------------

bids_subj=BIDS/sub-${subj}

# Either a session is given on the command line
# If not the session(s) need to be determined.
if [ $s_flag -eq 1 ]; then

    # session is given on the command line
    search_sessions=BIDS/sub-${subj}/ses-${ses}

else

    # search if any sessions exist
    search_sessions=($(find BIDS/sub-${subj} -type d | grep dwi))

fi    
 
num_sessions=${#search_sessions[@]}
    
echo "  Number of BIDS sessions: $num_sessions"
echo "    notably: ${search_sessions[@]}"


# ---- BIG LOOP for processing each session
for i in `seq 0 $(($num_sessions-1))`; do

# set up directories 
cd $cwd
long_bids_subj=${search_sessions[$i]}
#echo $long_bids_subj
bids_subj=${long_bids_subj%dwi}

# Create the Directory to write preprocessed data in
preproc=dwiprep/sub-${subj}/$(basename $bids_subj) 
#echo $preproc

# Directory to put raw mif data in
raw=${preproc}/raw


kul_e2cl " Start processing $bids_subj" ${preproc}/${log}

#---------- MAIN ---------------------------------------------------------------------------------------
echo $subj
echo $nods


# STEP 1 - PROCESSING  ---------------------------------------------
cd ${preproc}

# Where is the freesurfer parcellation? 
fs_aparc=${cwd}/freesurfer/sub-${subj}/${subj}/mri/aparc+aseg.mgz

# Where is the T1w anat?
ants_anat=T1w/T1w_BrainExtractionBrain.nii.gz

# Convert FS aparc back to original space
mkdir -p roi
fs_labels=roi/labels_from_FS.nii.gz
if [ ! -f $fs_labels ]; then
    mri_convert -rl $ants_anat -rt nearest $fs_aparc $fs_labels
fi

# 5tt segmentation & tracking
mkdir -p 5tt

if [ ! -f 5tt/5tt2gmwmi.nii.gz ]; then

    kul_e2cl " Performig 5tt..." ${log}
    #5ttgen fsl $ants_anat 5tt/5ttseg.mif -premasked -nocrop -force -nthreads $ncpu 
    #5ttgen freesurfer $fs_aparc 5tt/5ttseg.mif -nocrop -force -nthreads $ncpu
    5ttgen freesurfer $fs_labels 5tt/5ttseg.mif -nocrop -force -nthreads $ncpu
    
    5ttcheck -masks 5tt/failed_5tt 5tt/5ttseg.mif -force -nthreads $ncpu 
    5tt2gmwmi 5tt/5ttseg.mif 5tt/5tt2gmwmi.nii.gz -force 

else

    echo " 5tt already done, skipping..."

fi

# Extract relevant freesurfer determined rois
if [ ! -f roi/WM_fs_R.nii.gz ]; then

    kul_e2cl " Making the Freesurfer ROIS from subject space..." ${log}

    # M1_R is 2024
    fslmaths $fs_labels -thr 2024 -uthr 2024 -bin roi/M1_fs_R
    # M1_L is 1024
    fslmaths $fs_labels -thr 1024 -uthr 1024 -bin roi/M1_fs_L
    # S1_R is 2022
    fslmaths $fs_labels -thr 2022 -uthr 2022 -bin roi/S1_fs_R
    # S1_L is 1024
    fslmaths $fs_labels -thr 1022 -uthr 1022 -bin roi/S1_fs_L
    # Thalamus_R is 49
    fslmaths $fs_labels -thr 49 -uthr 49 -bin roi/THALAMUS_fs_R
    # Thalamus_L is 10
    fslmaths $fs_labels -thr 10 -uthr 10 -bin roi/THALAMUS_fs_L
    # SMA_and_PMC_L are
    # 1003    ctx-lh-caudalmiddlefrontal
    # 1028    ctx-lh-superiorfrontal
    fslmaths $fs_labels -thr 1003 -uthr 1003 -bin roi/MFG_fs_L
    fslmaths $fs_labels -thr 1028 -uthr 1028 -bin roi/SFG_fs_L
    fslmaths roi/MFG_fs_L -add roi/SFG_fs_L -bin roi/SMA_and_PMC_fs_L
    # SMA_and_PMC_L are
    # 2003    ctx-lh-caudalmiddlefrontal
    # 2028    ctx-lh-superiorfrontal
    fslmaths $fs_labels -thr 2003 -uthr 2003 -bin roi/MFG_fs_R
    fslmaths $fs_labels -thr 2028 -uthr 2028 -bin roi/SFG_fs_R
    fslmaths roi/MFG_fs_R -add roi/SFG_fs_R -bin roi/SMA_and_PMC_fs_R
    # 41  Right-Cerebral-White-Matter
    fslmaths $fs_labels -thr 41 -uthr 41 -bin roi/WM_fs_R
    # 2   Left-Cerebral-White-Matter
    fslmaths $fs_labels -thr 2 -uthr 2 -bin roi/WM_fs_L


else

    echo " Making the Freesurfer ROIS has been done already, skipping" 

fi

function KUL_antsApply_Transform {

    antsApplyTransforms -d 3 --float 1 \
    --verbose 1 \
    -i $input \
    -o $output \
    -r $reference \
    -t $transform \
    -n Linear
}

if [ ! -f roi/DENTATE_L.nii.gz ]; then

    kul_e2cl " Warping the SUIT3.3 atlas ROIS of the DENTATE to subject space..." ${log}
    # transform the T1w into MNI space using fmriprep data
    input=$ants_anat
    output=T1w/T1w_MNI152NLin2009cAsym.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
    reference=/KUL_apps/fsl/data/standard/MNI152_T1_1mm.nii.gz
    KUL_antsApply_Transform

    # inversly transform the T1w in MNI space to subject space (for double checking)
    input=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz
    output=T1w/T1w_test_inv_MNI_warp.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform

    # We get the Dentate rois out of MNI space, from the SUIT v3.3 atlas
    # http://www.diedrichsenlab.org/imaging/suit_download.htm
    # fslmaths Cerebellum-SUIT.nii -thr 30 -uthr 30 Dentate_R
    # fslmaths Cerebellum-SUIT.nii -thr 29 -uthr 29 Dentate_L
    input=${kul_main_dir}/atlasses/Local/Dentate_R.nii.gz
    output=roi/DENTATE_R.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform

    input=${kul_main_dir}/atlasses/Local/Dentate_L.nii.gz
    output=roi/DENTATE_L.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform

    if [ $Donatienne -eq 1 ]; then

    # We get the SN & PATUMEN rois out of MNI space, from Donatienne's PET data
    # We warp them back to individual subject space
    input=${cwd}/ROIS/rsn_l.nii
    output=roi/SUBNIG_L.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform

    input=${cwd}/ROIS/rputamen_l.nii
    output=roi/PUTAMEN_L.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform

    input=${cwd}/ROIS/rsn_r.nii
    output=roi/SUBNIG_R.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform

    input=${cwd}/ROIS/rputamen_r.nii
    output=roi/PUTAMEN_R.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform

    fi

else

    echo " Warping the SUIT3.3 atlas ROIS of the DENTATE to subject space has been done already, skipping" 

fi


# STEP 2 - Tractography  ---------------------------------------------

function kul_mrtrix_tracto_drt {

    #for a in iFOD2 Tensor_Prob; do
    for a in iFOD2; do
    
        # do the tracking
        echo tracts_${a}/${tract}.tck
        
        if [ ! -f tracts_${a}/${tract}.tck ]; then 

            mkdir -p tracts_${a}

            # make the intersect string (this is the first of the seeds)
            intersect=${seeds%% *}

            kul_e2cl " running tckgen of ${tract} tract with algorithm $a (all seeds with -select $nods, intersect with $intersect)" ${log}

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
                tckgen $wmfod tracts_${a}/${tract}.tck -algorithm $a -select $nods $s $i $e $m -angle $theta -nthreads $ncpu -force

            else

                # perform Tensor_Prob tckgen
                tckgen $dwi_preproced tracts_${a}/${tract}.tck -algorithm $a -cutoff 0.01 -select $nods $s $i $e $m -nthreads $ncpu -force

            fi
        
        else

            echo "  tckgen of ${tract} tract already done, skipping"

        fi

        # Check if any fibers have been found & log to the information file
        echo "   checking tracts_${a}/${tract}"
        local count=$(tckinfo tracts_${a}/${tract}.tck | grep count | head -n 1 | awk '{print $(NF)}')
        echo "$subj, $a, $tract, $count" >> tracts_info.csv

        # do further processing of tracts are found
        if [ ! -f tracts_${a}/MNI_Space_${tract}_${a}.nii.gz ]; then

            if [ $count -eq 0 ]; then

                # report that no tracts were found and stop further processing
                kul_e2cl "  no streamlines were found for the tracts_${a}/${tract}.tck" ${log}

            else

                # report how many tracts were found and continue processing
                echo "   $count streamlines were found for the tracts_${a}/${tract}.tck"
                
                echo "   generating subject/MNI space images"
                # convert the tck in nii
                tckmap tracts_${a}/${tract}.tck tracts_${a}/${tract}.nii.gz -template $ants_anat -force 

                # Warp the full tract image to MNI space
                input=tracts_${a}/${tract}.nii.gz
                output=tracts_${a}/MNI_Space_FULL_${tract}_${a}.nii.gz
                transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
                reference=/KUL_apps/fsl/data/standard/MNI152_T1_1mm.nii.gz
                KUL_antsApply_Transform

                # intersect the nii tract image with the thalamic roi
                fslmaths tracts_${a}/${tract}.nii -mas roi/${intersect}.nii.gz tracts_${a}/${tract}_masked

                # make a probabilistic image
                local m=$(mrstats -quiet tracts_${a}/${tract}_masked.nii.gz -output max)
                fslmaths tracts_${a}/${tract}_masked -div $m tracts_${a}/Subj_Space_${tract}_${a}

                # Warp the probabilistic image to MNI space
                input=tracts_${a}/Subj_Space_${tract}_${a}.nii.gz
                output=tracts_${a}/MNI_Space_${tract}_${a}.nii.gz
                transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
                reference=/KUL_apps/fsl/data/standard/MNI152_T1_1mm.nii.gz
                KUL_antsApply_Transform
                
                
                if [ $count -lt $do_sift_th ]; then
                
                    kul_e2cl "  NOT running tckshift since less than $do_sift_th streamlines" ${log}

                else
                    kul_e2cl "  running tckshift & generation subject/MNI space images" ${log}

                    # perform filtering on tracts with tcksift (version 1)
                    tcksift -term_ratio $term_ratio -act ${cwd}/dwiprep/sub-${subj}/5tt/5ttseg.mif tracts_${a}/${tract}.tck \
                        $wmfod tracts_${a}/sift1_${tract}.tck -nthreads $ncpu -force

                    # convert the tck in nii
                    tckmap tracts_${a}/sift1_${tract}.tck tracts_${a}/sift1_${tract}.nii.gz -template $ants_anat -force 

                    # intersect the nii tract image with the thalamic roi
                    fslmaths tracts_${a}/sift1_${tract}.nii -mas roi/${intersect}.nii.gz tracts_${a}/sift1_${tract}_masked

                    # make a probabilistic image
                    local m=$(mrstats -quiet tracts_${a}/sift1_${tract}_masked.nii.gz -output max)
                    fslmaths tracts_${a}/sift1_${tract}_masked -div $m tracts_${a}/Subj_Space_sift1_${tract}_${a}

                    # Warp the probabilistic image to MNI space
                    input=tracts_${a}/Subj_Space_sift1_${tract}_${a}.nii.gz
                    output=tracts_${a}/MNI_Space_sift1_${tract}_${a}.nii.gz
                    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
                    reference=/KUL_apps/fsl/data/standard/MNI152_T1_1mm.nii.gz
                    KUL_antsApply_Transform

                fi

            fi

        else
        
            echo "  tckshift & generation subject/MNI space images already done, skipping..."
        
        fi
    
    done

}

wmfod=response/wmfod_reg2T1w.mif
dwi_preproced=dwi_preproced_reg2T1w.mif

# Make an empty log file with information about the tracts
echo "subject, algorithm, tract, count" > tracts_info.csv

# M1_fs-Thalamic tracts
tract="TH-M1_fs_R_nods${nods}"
seeds=("THALAMUS_fs_R" "M1_fs_R")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-M1_fs_L_nods${nods}"
seeds=("THALAMUS_fs_L" "M1_fs_L")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

# S1_fs-Thalamic tracts
tract="TH-S1_fs_R_nods${nods}"
seeds=("THALAMUS_fs_R" "S1_fs_R")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-S1_fs_L_nods${nods}"
seeds=("THALAMUS_fs_L" "S1_fs_L")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

# SMA_and_PMC-Thalamic tracts
tract="TH-SMA_and_PMC_R_nods${nods}"
seeds=("THALAMUS_fs_R" "SMA_and_PMC_fs_R")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-SMA_and_PMC_L_nods${nods}"
seeds=("THALAMUS_fs_L" "SMA_and_PMC_fs_L")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt  

# Dentato-Rubro_Thalamic tracts
tract="TH-DR_R_nods${nods}"
seeds=("THALAMUS_fs_R" "M1_fs_R" "DENTATE_L")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-DR_L_nods${nods}"
seeds=("THALAMUS_fs_L" "M1_fs_L" "DENTATE_R")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

if [ $Donatienne -eq 1 ]; then
tract="NST_L_nods${nods}"
seeds=("SUBNIG_L" "PUTAMEN_L")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt 

tract="NST_R_nods${nods}"
seeds=("SUBNIG_R" "PUTAMEN_R")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 
fi


# Now prepare the data for iPlan
if [ ! -f for_iplan/TH_SMAPMC_R.hdr ]; then

    mkdir -p for_iplan

    # copy the tracts in analyze format
    fslmaths tracts_iFOD2/TH-DR_L_nods${nods} -s 0.5 -thr 4 -bin for_iplan/Tract_DRT_L
    fslchfiletype NIFTI_PAIR for_iplan/Tract_DRT_L for_iplan/Tract_DRT_L

    fslmaths tracts_iFOD2/TH-DR_R_nods${nods} -s 0.5 -thr 4 -bin for_iplan/Tract_DRT_R
    fslchfiletype NIFTI_PAIR for_iplan/Tract_DRT_R for_iplan/Tract_DRT_R

    # copy the T1w in analyze format
    cp T1w/T1w_BrainExtractionBrain.nii.gz for_iplan/anat.nii.gz
    fslchfiletype NIFTI_PAIR for_iplan/anat for_iplan/anat

    # copy the Thalamic probabilistic images as speudo fmri activation maps
    fslmaths tracts_iFOD2/Subj_Space_TH-DR_L_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_DRT_L
    fslchfiletype NIFTI_PAIR for_iplan/TH_DRT_L for_iplan/TH_DRT_L
    fslmaths tracts_iFOD2/Subj_Space_TH-DR_R_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_DRT_R
    fslchfiletype NIFTI_PAIR for_iplan/TH_DRT_R for_iplan/TH_DRT_R

    fslmaths tracts_iFOD2/Subj_Space_TH-M1_fs_L_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_M1_L
    fslchfiletype NIFTI_PAIR for_iplan/TH_M1_L for_iplan/TH_M1_L
    fslmaths tracts_iFOD2/Subj_Space_TH-M1_fs_R_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_M1_R
    fslchfiletype NIFTI_PAIR for_iplan/TH_M1_R for_iplan/TH_M1_R

    fslmaths tracts_iFOD2/Subj_Space_TH-S1_fs_L_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_S1_L
    fslchfiletype NIFTI_PAIR for_iplan/TH_S1_L for_iplan/TH_S1_L
    fslmaths tracts_iFOD2/Subj_Space_TH-S1_fs_R_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_S1_R
    fslchfiletype NIFTI_PAIR for_iplan/TH_S1_R for_iplan/TH_S1_R

    fslmaths tracts_iFOD2/Subj_Space_TH-SMA_and_PMC_L_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_SMAPMC_L
    fslchfiletype NIFTI_PAIR for_iplan/TH_SMAPMC_L for_iplan/TH_SMAPMC_L
    fslmaths tracts_iFOD2/Subj_Space_TH-SMA_and_PMC_R_nods${nods}_iFOD2.nii.gz -s 0.5 -thr 0.25 for_iplan/TH_SMAPMC_R
    fslchfiletype NIFTI_PAIR for_iplan/TH_SMAPMC_R for_iplan/TH_SMAPMC_R

    # clean up
    rm -rf for_iplan/*.nii.gz

fi

# ---- END BIG LOOP for processing each session
done

exit 0











---- EVERYTHING BELOW IS OLD - needs to be removed in beta version

# M1-Thalamic tracts
tract="TH-M1_R"
seeds=("THALAMUS_R" "M1")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-M1_L"
seeds=("THALAMUS_L" "M1")
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

# SMA_and_PMC-Thalamic tracts
tract="TH-SMA_and_PMC_R"
seeds=("THALAMUS_R" "SMA_and_PMC")
exclude="WM_fs_L"
kul_mrtrix_tracto_drt 

tract="TH-SMA_and_PMC_L"
seeds=("THALAMUS_L" "SMA_and_PMC")
exclude="WM_fs_R"
kul_mrtrix_tracto_drt  


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
