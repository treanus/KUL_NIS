#!/bin/bash -e
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#  Developed for generating major fiber bundles for presurgical mapping with Tensor_Prof abd iFOD2 msmt_CSD
#	for S61759.
#	 Project PI: Stefan Sunaert
#
# Requires Mrtrix3, FSL, ants, freesurfer
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@kuleuven.be
#
# v0.1 - dd 02/02/2019 - Dev - AR
v="v0.1 - dd 02/02/2019"

# To Do
#  - add iFOD2 and Tensor_Prob/Tensor_Det fiber tracking for: 
#  - CC (Genu, Body and Splenium ?) or all segments, VOF, AC, Fornix, TIF, AIF


# A few fixed (for now) parameters:

    # Number of desired streamlines
    # nods=2000
    # this has become a command line optional parameter

    # Maximum angle between successive steps for iFOD2
	# still harcoded but now defined separately for each bundle
	# added $stop for stopping or not, depending on the bundle tracked
	# need to also define another var for nods per bundle
	# another var also for thr for filtering (should be related to how many fibers I get from first tckgen e.g. 0.2% of 10k fibers would be 2 fibers)


    # sift1 filtering
    # termination ratio - defined as the ratio between reduction in cost
    # function, and reduction in density of streamlines.
    # Smaller values result in more streamlines being filtered out.
    do_sift_th=10000 # when to do sift? (if more than 5000 streamlines in tract e.g.)
    term_ratio=0.5 # reduce by e.g. 50%

    # tmp directory for temporary processing
    tmp=/tmp




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

`basename $0` performs fiber tractography.

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
# nods=2000
ncpu=6
silent=1

# Set required options
p_flag=0
s_flag=0

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
for current_session in `seq 0 $(($num_sessions-1))`; do

    # set up directories 
    cd $cwd
    long_bids_subj=${search_sessions[$current_session]}
    echo $long_bids_subj
    bids_subj=${long_bids_subj%dwi}

    # Create the Directory to write preprocessed data in
    preproc=dwiprep/sub-${subj}/$(basename $bids_subj) 
    echo $preproc

    kul_e2cl " Start processing $bids_subj" ${preproc}/${log}

    
    # STEP 1 - PROCESSING  ---------------------------------------------
    cd ${preproc}

    # Where is the freesurfer parcellation? 
    fs_aparc=${cwd}/freesurfer/sub-${subj}/${subj}/mri/aparc+aseg.mgz
	
	# there are other parcellation files we're interested in namely:
	# wm_parc 
	fs_wm_parc=${cwd}/freesurfer/sub-${subj}/${subj}/mri/wmparc.mgz
	fs_wm_lobes=${cwd}/freesurfer/sub-${subj}/${subj}/mri/wmparc.lobes.mgz
	
	# lobe specific wm segmentations from wm+parc
	if [ ! -f $fs_wm_lobes ]; then
	mri_annotation2label --subject ${subj} --sd ${cwd}/freesurfer/sub-${subj} --hemi lh --lobesStrict lobes
	mri_annotation2label --subject ${subj} --sd ${cwd}/freesurfer/sub-${subj} --hemi rh --lobesStrict lobes
	mri_aparc2aseg --s ${subj} --sd ${cwd}/freesurfer/sub-${subj}  --labelwm --hypo-as-wm --rip-unknown \
	  --volmask --o ${cwd}/freesurfer/sub-${subj}/${subj}/mri/wmparc.lobes.mgz --ctxseg aparc+aseg.mgz \
	  --annot lobes --base-offset 200
	else 
		
		echo " wm lobe specific labels already done, skipping ..."
		
	fi
	
    # Where is the T1w anat?
    ants_anat=T1w/T1w_BrainExtractionBrain.nii.gz
	ants_mask=T1w/T1w_BrainExtractionMask.nii.gz

	if [ ! -f $ants_mask ]; then
		fslmaths $ants_anat -bin $ants_mask
	else
		echo " T1 brain mask already done, skipping ..."
	fi


    # Convert FS aparc back to original space
    mkdir -p roi
    fs_labels=roi/labels_from_FS.nii.gz
    fs_wm_labels=roi/wm_labels_from_FS.nii.gz
	fs_wm_lobe_labels=roi/wm_lobes_labels_from_FS.nii.gz
    mri_convert -rl $ants_anat -rt nearest $fs_aparc $fs_labels
	mri_convert -rl $ants_anat -rt nearest $fs_wm_parc $fs_wm_labels
	mri_convert -rl $ants_anat -rt nearest $fs_wm_lobes $fs_wm_lobe_labels
	
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
    # Add the paracentral lobules for the rest of the SM cortex
    fslmaths $fs_labels -thr 2017 -uthr 2017 -bin roi/PCent_fs_R
    fslmaths $fs_labels -thr 1017 -uthr 1017 -bin roi/PCent_fs_L
	# superior parietal for STR
    fslmaths $fs_labels -thr 2029 -uthr 2029 -bin roi/SPL_fs_R
    fslmaths $fs_labels -thr 1029 -uthr 1029 -bin roi/SPL_fs_L
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
	# add brainstem as a VOI
	fslmaths $fs_labels -thr 16 -uthr 16 -bin roi/BStem
	# add IFG PTr and POp w/wo POr
	# 2020	ctx-rh-parstriangularis
	# 2019	ctx-rh-parsorbitalis
	# 2018	ctx-rh-parsopercularis
	# 1020	ctx-lh-parstriangularis
	# 1019	ctx-lh-parsorbitalis
	# 1018	ctx-lh-parsopercularis
	fslmaths $fs_labels -thr 2020 -uthr 2020 -bin roi/IFG_PTr_fs_R	
	fslmaths $fs_labels -thr 2018 -uthr 2018 -bin roi/IFG_POp_fs_R	
	fslmaths $fs_labels -thr 1020 -uthr 1020 -bin roi/IFG_PTr_fs_L	
	fslmaths $fs_labels -thr 1018 -uthr 1018 -bin roi/IFG_POp_fs_L
	# add FP and TP
	# 2032	ctx-rh-frontalpole
	# 2033	ctx-rh-temporalpole
	# 1032	ctx-lh-frontalpole
	# 1033	ctx-rh-temporalpole	
	fslmaths $fs_labels -thr 2032 -uthr 2032 -bin roi/FP_fs_L
	fslmaths $fs_labels -thr 2033 -uthr 2033 -bin roi/TP_fs_L	
	fslmaths $fs_labels -thr 1032 -uthr 1032 -bin roi/FP_fs_R	
	fslmaths $fs_labels -thr 1033 -uthr 1033 -bin roi/TP_fs_R	
	# add temporal gyri
	fslmaths $fs_labels -thr 1009 -uthr 1009 -bin roi/ITG_fs_L
	fslmaths $fs_labels -thr 2009 -uthr 2009 -bin roi/ITG_fs_R
	fslmaths $fs_labels -thr 1015 -uthr 1015 -bin roi/MTG_fs_L
	fslmaths $fs_labels -thr 2015 -uthr 2015 -bin roi/MTG_fs_R
	fslmaths $fs_labels -thr 1030 -uthr 1030 -bin roi/STG_fs_L
	fslmaths $fs_labels -thr 2030 -uthr 2030 -bin roi/STG_fs_R
	fslmaths $fs_labels -thr 1001 -uthr 1001 -bin roi/bSTG_fs_L	
	fslmaths $fs_labels -thr 2001 -uthr 2001 -bin roi/bSTG_fs_R
	fslmaths $fs_labels -thr 1031 -uthr 1031 -bin roi/SMG_fs_L
	fslmaths $fs_labels -thr 2031 -uthr 2031 -bin roi/SMG_fs_R
	# add Orbito_frontal ctx
	fslmaths $fs_labels -thr 1012 -uthr 1012 -bin roi/L_OF_fs_L
	fslmaths $fs_labels -thr 2012 -uthr 2012 -bin roi/L_OF_fs_R
	fslmaths $fs_labels -thr 1014 -uthr 1014 -bin roi/M_OF_fs_L
	fslmaths $fs_labels -thr 2014 -uthr 2014 -bin roi/M_OF_fs_R
	fslmaths $fs_labels -thr 1019 -uthr 1019 -bin roi/IFG_POr_fs_L
	fslmaths $fs_labels -thr 2019 -uthr 2019 -bin roi/IFG_POr_fs_R
	# add entorhinal cortex
	fslmaths $fs_labels -thr 1006 -uthr 1006 -bin roi/enR_fs_L
	fslmaths $fs_labels -thr 2006 -uthr 2006 -bin roi/enR_fs_R
	# add Insula
	# 1035	ctx-lh-insula
	# 2035 	ctx-rh-insula
	fslmaths $fs_labels -thr 1035 -uthr 1035 -bin roi/Ins_fs_L
	fslmaths $fs_labels -thr 2035 -uthr 2035 -bin roi/Ins_fs_R			
	# add corpus callosum
	# 251 - 255 CC Post,MidPost,Central,MidAnt,Ant
	fslmaths $fs_labels -thr 251 -uthr 251 -bin roi/CC_fs_post			
	fslmaths $fs_labels -thr 252 -uthr 252 -bin roi/CC_fs_midpost
	fslmaths $fs_labels -thr 253 -uthr 253 -bin roi/CC_fs_central
	fslmaths $fs_labels -thr 254 -uthr 254 -bin roi/CC_fs_midant
	fslmaths $fs_labels -thr 255 -uthr 255 -bin roi/CC_fs_ant
	fslmaths roi/CC_fs_ant.nii.gz -add roi/CC_fs_midant.nii.gz -add roi/CC_fs_central.nii.gz \
		-add roi/CC_fs_midpost.nii.gz -add roi/CC_fs_post.nii.gz -bin roi/CC_fs_all
	# add ACC
	# 1002	ctx-lh-caudalanteriorcingulate
	# 1026  ctx-lh-rostralanteriorcingulate
	# 2002	ctx-lh-caudalanteriorcingulate
	# 2026  ctx-lh-rostralanteriorcingulate
	fslmaths $fs_labels -thr 1002 -uthr 1002 -bin roi/cACC_fs_L
	fslmaths $fs_labels -thr 1026 -uthr 1026 -bin roi/rACC_fs_L
	fslmaths $fs_labels -thr 2002 -uthr 2002 -bin roi/cACC_fs_R	
	fslmaths $fs_labels -thr 2026 -uthr 2026 -bin roi/rACC_fs_R
	# add PCC
	# 1023	ctx-lh-posteriorcingulate
	# 2023	ctx-rh-posteriorcingulate
	# 1010	ctx-lh-isthmuscingulate
	# 2010	ctx-rh-isthmuscingulate
	fslmaths $fs_labels -thr 1023 -uthr 1023 -bin roi/PCC_fs_L
	fslmaths $fs_labels -thr 2023 -uthr 2023 -bin roi/PCC_fs_R
	fslmaths $fs_labels -thr 1010 -uthr 1010 -bin roi/iPCC_fs_L
	fslmaths $fs_labels -thr 2010 -uthr 2010 -bin roi/iPCC_fs_R
	# hippocampi
	# 17	Left-Hippocampus
	# 18	Left-Amygdala
	# 53	Right-Hippocampus
	# 54	Right-Amygdala
	fslmaths $fs_labels -thr 17 -uthr 17 -bin roi/Hippo_fs_L	
	fslmaths $fs_labels -thr 53 -uthr 53 -bin roi/Hippo_fs_R
	fslmaths $fs_labels -thr 18 -uthr 18 -bin roi/Amyg_fs_L
	fslmaths $fs_labels -thr 54 -uthr 54 -bin roi/Amyg_fs_R
	# Pericalcarine cortex
	# 2021	ctx-rh-pericalcarine
	# 1021	ctx-rh-pericalcarine
	fslmaths $fs_labels -thr 2021 -uthr 2021 -bin roi/periCalc_fs_R
	fslmaths $fs_labels -thr 1021 -uthr 1021 -bin roi/periCalc_fs_L
	# Ventral DC VOIs for exclude in OR tracking
	fslmaths $fs_labels -thr 28 -uthr 28 -bin roi/vDC_fs_L
	fslmaths $fs_labels -thr 60 -uthr 60 -bin roi/vDC_fs_R
	# Lateral ventricles may also be used as exclude ROIs
	fslmaths $fs_labels -thr 4 -uthr 4 -bin roi/Lat_V_fs_L
	fslmaths $fs_labels -thr 43 -uthr 43 -bin roi/Lat_V_fs_R
	# Caudates may also be used as exclude for ORs
	fslmaths $fs_labels -thr 11 -uthr 11 -bin roi/Caudate_fs_L
	fslmaths $fs_labels -thr 50 -uthr 50 -bin roi/Caudate_fs_R
	# Optic Chiasm
	fslmaths $fs_labels -thr 85 -uthr 85 -bin roi/OC_fs
	# WM tract of the limbic system
	fslmaths $fs_wm_labels -thr 3010 -uthr 3010 -bin roi/iPCC_wm_fs_L
	fslmaths $fs_wm_labels -thr 4010 -uthr 4010 -bin roi/iPCC_wm_fs_R
	fslmaths $fs_wm_labels -thr 3023 -uthr 3023 -bin roi/PCC_wm_fs_L
	fslmaths $fs_wm_labels -thr 4023 -uthr 4023 -bin roi/PCC_wm_fs_R
	# Putamina for exclusions
	fslmaths $fs_labels -thr 12 -uthr 12 -dilM -bin roi/putamen_dil_fs_L
	fslmaths $fs_labels -thr 51 -uthr 51 -dilM -bin roi/putamen_dil_fs_R
	# Putamina without dilation
	fslmaths $fs_labels -thr 12 -uthr 12 -dilM -bin roi/putamen_dil_fs_L
	fslmaths $fs_labels -thr 51 -uthr 51 -dilM -bin roi/putamen_dil_fs_R
	# accumbens
	fslmaths $fs_wm_labels -thr 26 -uthr 26 -dilM -bin roi/accumbens_dil_fs_L
	fslmaths $fs_wm_labels -thr 58 -uthr 58 -dilM -bin roi/accumbens_dil_fs_R
	# add the CSF tpm from ANTs for further cleanup (still trying this out)
	fslmaths $cwd/fmriprep/sub-${subj}/anat/sub-${subj}_label-CSF_probseg.nii.gz \
		-thr 0.2 $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/CSF_ants.nii.gz
	fslmaths $cwd/fmriprep/sub-${subj}/anat/sub-${subj}_label-CSF_probseg.nii.gz \
		-thr 0.2 -binv $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/CSF_ants_binv.nii.gz
	# add the GM tpm from ANTs for further cleanup (still trying this out)
	fslmaths $cwd/fmriprep/sub-${subj}/anat/sub-${subj}_label-GM_probseg.nii.gz \
		-thr 0.2 $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/GM_ants.nii.gz
	# Lobe specific wm labels
	fslmaths $fs_wm_lobe_labels -thr 3204 -uthr 3204 -bin roi/Occ_wm_fs_L
	fslmaths $fs_wm_lobe_labels -thr 4204 -uthr 4204 -bin roi/Occ_wm_fs_R
	fslmaths $fs_wm_lobe_labels -thr 3207 -uthr 3207 -bin roi/Ins_wm_fs_L
	fslmaths $fs_wm_lobe_labels -thr 4204 -uthr 4204 -bin roi/Ins_wm_fs_R
	fslmaths $fs_wm_lobe_labels -thr 3205 -uthr 3205 -bin roi/Temp_wm_fs_L
	fslmaths $fs_wm_lobe_labels -thr 4205 -uthr 4205 -bin roi/Temp_wm_fs_R
	fslmaths $fs_wm_lobe_labels -thr 3201 -uthr 3201 -bin roi/Front_wm_fs_L
	fslmaths $fs_wm_lobe_labels -thr 4201 -uthr 4201 -bin roi/Front_wm_fs_R
	fslmaths $fs_wm_lobe_labels -thr 3206 -uthr 3206 -bin roi/Par_wm_fs_L
	fslmaths $fs_wm_lobe_labels -thr 4206 -uthr 4206 -bin roi/Par_wm_fs_R
	# SM1 and paracentral lobule wm
	fslmaths $fs_wm_labels -thr 4022 -uthr 4022 -bin roi/S1_wm_fs_R	
	fslmaths $fs_wm_labels -thr 3022 -uthr 3022 -bin roi/S1_wm_fs_L	
	fslmaths $fs_wm_labels -thr 4024 -uthr 4024 -bin roi/M1_wm_fs_R	
	fslmaths $fs_wm_labels -thr 3024 -uthr 3024 -bin roi/M1_wm_fs_L	
	fslmaths $fs_wm_labels -thr 3017 -uthr 3017 -bin roi/PCent_wm_fs_R
	fslmaths $fs_wm_labels -thr 4017 -uthr 4017 -bin roi/PCent_wm_fs_L
	# SPL wm
	fslmaths $fs_wm_labels -thr 4029 -uthr 4029 -bin roi/SPL_wm_fs_R
	fslmaths $fs_wm_labels -thr 3029 -uthr 3029 -bin roi/SPL_wm_fs_L	
	# Temporal gyri wm
	fslmaths $fs_wm_labels -thr 4030 -uthr 4030 -bin roi/STG_wm_fs_R	
	fslmaths $fs_wm_labels -thr 3030 -uthr 3030 -bin roi/STG_wm_fs_L
	fslmaths $fs_wm_labels -thr 4015 -uthr 4015 -bin roi/MTG_wm_fs_R	
	fslmaths $fs_wm_labels -thr 3015 -uthr 3015 -bin roi/MTG_wm_fs_L
	fslmaths $fs_wm_labels -thr 3009 -uthr 3009 -bin roi/ITG_wm_fs_L
	fslmaths $fs_wm_labels -thr 4009 -uthr 4009 -bin roi/ITG_wm_fs_R
	fslmaths $fs_wm_labels -thr 3007 -uthr 3007 -bin roi/fusiform_wm_fs_L
	fslmaths $fs_wm_labels -thr 4007 -uthr 4007 -bin roi/fusiform_wm_fs_R
	# supramarginal gyri
	fslmaths $fs_wm_labels -thr 4031 -uthr 4031 -bin roi/SMG_wm_fs_R
	fslmaths $fs_wm_labels -thr 3031 -uthr 3031 -bin roi/SMG_wm_fs_L
	# Inferior frontal gyri
	fslmaths $fs_wm_labels -thr 4020 -uthr 4020 -bin roi/IFG_PTr_wm_fs_R	
	fslmaths $fs_wm_labels -thr 3020 -uthr 3020 -bin roi/IFG_PTr_wm_fs_L
	fslmaths $fs_wm_labels -thr 4018 -uthr 4018 -bin roi/IFG_POp_wm_fs_R
	fslmaths $fs_wm_labels -thr 3018 -uthr 3018 -bin roi/IFG_POp_wm_fs_L
	# Lateral orbitofrontal white matter
	fslmaths $fs_wm_labels -thr 3012 -uthr 3012 -bin roi/L_OF_wm_fs_L
	fslmaths $fs_wm_labels -thr 4012 -uthr 4012 -bin roi/L_OF_wm_fs_R
	# Temporal pole white matter
	fslmaths $fs_wm_labels -thr 3033 -uthr 3033 -bin roi/TP_wm_fs_L
	fslmaths $fs_wm_labels -thr 4033 -uthr 4033 -bin roi/TP_wm_fs_R
	# parahippocampal white matter 4016 rt and 3016 lt
	fslmaths $fs_wm_labels -thr 3016 -uthr 3016 -bin roi/Phip_wm_fs_L
	fslmaths $fs_wm_labels -thr 4016 -uthr 4016 -bin roi/Phip_wm_fs_R
	# fs unsegmented wm can also be useful
	fslmaths $fs_wm_labels -thr 5001 -uthr 5001 -bin roi/unseg_wm_fs_L
	fslmaths $fs_wm_labels -thr 5002 -uthr 5002 -bin roi/unseg_wm_fs_R
	# The CST, SMA/CST and ML need subsegmentation of the brainstem
	# This can be done using the JHU white matter labels ROIs

	# first we define the template and atlas for JHU_wm_labels
	JHU_temp=$FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz
	JHU_labels=$FSLDIR/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz

	# Extract brainstem ROIs
	# should make ones for including (with smooth and a low thr) and ones for excluding (with ero)
	fslmaths $JHU_labels -thr 8 -uthr 8 -s 1.5 -thr 0.3 -bin roi/LT_CST_JHU_roi_MNI.nii.gz
	fslmaths $JHU_labels -thr 7 -uthr 7 -s 1.5 -thr 0.3 -bin roi/RT_CST_JHU_roi_MNI.nii.gz
	fslmaths $JHU_labels -thr 10 -uthr 10 -s 1.5 -thr 0.3 -bin roi/LT_ML_JHU_roi_MNI.nii.gz
	fslmaths $JHU_labels -thr 9 -uthr 9 -s 1.5 -thr 0.3 -bin roi/RT_ML_JHU_roi_MNI.nii.gz

	# Warp MNI image to native and apply warps to individual ROIs
	antsRegistrationSyNQuick.sh -d 3 -m $JHU_temp -f $ants_anat -t s -n $ncpu -j 1 -o T1w/JHU_temp_in_native
	WarpImageMultiTransform 3 $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/LT_CST_JHU_roi_MNI.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/LT_CST_JHU_roi_native.nii.gz -R $ants_anat \
		$cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native1Warp.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native0GenericAffine.mat
	WarpImageMultiTransform 3 $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/RT_CST_JHU_roi_MNI.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/RT_CST_JHU_roi_native.nii.gz -R $ants_anat \
		$cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native1Warp.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native0GenericAffine.mat
	WarpImageMultiTransform 3 $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/LT_ML_JHU_roi_MNI.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/LT_ML_JHU_roi_native.nii.gz -R $ants_anat \
		$cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native1Warp.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native0GenericAffine.mat
	WarpImageMultiTransform 3 $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/RT_ML_JHU_roi_MNI.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/RT_ML_JHU_roi_native.nii.gz -R $ants_anat \
		$cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native1Warp.nii.gz $cwd/dwiprep/sub-${subj}/sub-${subj}/T1w/JHU_temp_in_native0GenericAffine.mat

	# delete CST rois in MNI space
	rm $cwd/dwiprep/sub-${subj}/sub-${subj}/roi/*JHU_roi_MNI.nii.gz
	
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

# if [ ! -f roi/DENTATE_L.nii.gz ]; then
#
#     kul_e2cl " Warping the SUIT3.3 atlas ROIS of the DENTATE to subject space..." ${log}
#     # transform the T1w into MNI space using fmriprep data
#     input=$ants_anat
#     output=T1w/T1w_MNI152NLin2009cAsym.nii.gz
#     transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
#     reference=/KUL_apps/fsl/data/standard/MNI152_T1_1mm.nii.gz
#     KUL_antsApply_Transform
#
    # inversly transform the T1w in MNI space to subject space (for double checking)
    input=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz
    output=T1w/T1w_test_inv_MNI_warp.nii.gz
    transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    reference=$ants_anat
    KUL_antsApply_Transform
#
#     # We get the Dentate rois out of MNI space, from the SUIT v3.3 atlas
#     # http://www.diedrichsenlab.org/imaging/suit_download.htm
#     # fslmaths Cerebellum-SUIT.nii -thr 30 -uthr 30 Dentate_R
#     # fslmaths Cerebellum-SUIT.nii -thr 29 -uthr 29 Dentate_L
#     input=${kul_main_dir}/atlasses/Local/Dentate_R.nii.gz
#     output=roi/DENTATE_R.nii.gz
#     transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
#     reference=$ants_anat
#     KUL_antsApply_Transform
#
#     input=${kul_main_dir}/atlasses/Local/Dentate_L.nii.gz
#     output=roi/DENTATE_L.nii.gz
#     transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
#     reference=$ants_anat
#     KUL_antsApply_Transform
#
# else
#
#     echo " Warping the SUIT3.3 atlas ROIS of the DENTATE to subject space has been done already, skipping"
#
# fi



# STEP 2 - Tractography  ---------------------------------------------

function kul_mrtrix_FT {

    #for a in iFOD2 Tensor_Prob; do
		# where the hell does iFOD2 come from here ? need to ask SS
    for a in iFOD2; do
    
        # do the tracking
        # echo tracts_${a}/${tract}.tck
        
        if [ ! -f tracts_${a}/${tract}.tck ]; then 

            mkdir -p tracts_${a}

            # make the intersect string (this is the first of the seeds)
            intersect=${seeds%% *}
            
            echo $log

            kul_e2cl " running tckgen of ${tract} tract with algorithm $a all seeds with -select $nods, intersect with $intersect " ${log}

            # make the seed string
            local s=$(printf " -seed_image roi/%s.nii.gz"  "${seeds[@]}")
    
            # make the include string (which is same rois as seed)
            local i=$(printf " -include roi/%s.nii.gz"  "${seeds[@]}")

            # make the exclude string (which is same rois as seed)
            local e=$(printf " -exclude roi/%s.nii.gz"  "${exclude[@]}")

            # make the mask string 
            local m="-mask dwi_preproced_reg2T1w_mask.nii.gz"
			local mask="dwi_preproced_reg2T1w_mask.nii.gz"
			# question, why do we continue to use the suboptimal mask, and not the T1 brain mask 
			
			# make the string that will be used to create the tractogram filter mask
			local filter=$(printf " -add roi/%s.nii.gz" "${seeds[@]}")

			# make the string to set minimum length of fibers
			local min_L=$(printf " -minlength 20")

            if [ "${a}" == "iFOD2" ]; then

                # perform IFOD2 tckgen
				# now using the stop option to terminate streamlines within include rois
				if echo "$tract" | grep "_R_" > /dev/null 2>&1 ; then 
					echo "it's right sided"
                	tckgen $wmfod tracts_${a}/${tract}.tck -algorithm $a -select $nods $s -include roi/WM_fs_R.nii.gz $i $e $m -angle $theta $min_L -nthreads $ncpu $act $stop -force;
				elif echo "$tract" | grep "_L_" > /dev/null 2>&1 ; then 
					 echo "it's left sided"
                	tckgen $wmfod tracts_${a}/${tract}.tck -algorithm $a -select $nods $s -include roi/WM_fs_L.nii.gz $i $e $m -angle $theta $min_L -nthreads $ncpu $act $stop -force;
				fi
				

            else

                # perform Tensor_Prob tckgen
				if echo "$tract" | grep "_R_" > /dev/null 2>&1 ; then 
					echo "it's right sided"
                	tckgen $wmfod tracts_${a}/${tract}.tck -algorithm $a -select $nods $s -include roi/WM_fs_R.nii.gz $i $e $m -angle $theta $min_L -nthreads $ncpu $act $stop -force;
				elif echo "$tract" | grep "_L_" > /dev/null 2>&1 ; then 
					 echo "it's left sided"
                	tckgen $wmfod tracts_${a}/${tract}.tck -algorithm $a -select $nods $s -include roi/WM_fs_L.nii.gz $i $e $m -angle $theta $min_L -nthreads $ncpu $act $stop -force;
				fi

            fi
        
        else

            echo "  tckgen of ${tract} tract already done, skipping"

        fi

        # Check if any fibers have been found & log to the information file
        echo "   checking tracts_${a}/${tract}"
        local count=$(tckinfo tracts_${a}/${tract}.tck -count | grep count | head -n 1 | awk '{print $(NF)}')
		# setting a relative thr for tract filtering using 0.3% as thr cutoff
		act_count="$(echo $count | cut -d':' -f2)"
		thr="$(bc <<< "scale = 2; (($act_count*0.3/100))")"
		
        echo "$subj, $a, $tract, $count" >> tracts_FT_info.csv

        # do further processing if tracts are found
		# if [ ! -f tracts_${a}/MNI_Space_${tract}_${a}.nii.gz ]; then
        if [ ! -f tracts_${a}/${tract}.nii.gz ]; then

            if [ $count -eq 0 ]; then

                # report that no tracts were found and stop further processing
                kul_e2cl "  no streamlines were found for the tracts_${a}/${tract}.tck" ${log}

            else

                # report how many tracts were found and continue processing
                echo "   $count streamlines were found for the tracts_${a}/${tract}.tck"
                
                echo "   generating subject/MNI space images"
                # convert the tck in nii
                tckmap tracts_${a}/${tract}.tck tracts_${a}/${tract}.nii.gz -template $ants_anat -force 
                
				# skipping SIFT for now
				
                    kul_e2cl "  NOT running tcksift since less than $do_sift_th streamlines" ${log}
					
                    tckmap tracts_${a}/${tract}.tck tracts_${a}/${tract}.nii.gz -template $ants_anat -force 
					
					fslmaths tracts_${a}/${tract}.nii.gz -thr $thr -s 1.5 -thr 0.2 -bin tracts_${a}/${tract}_bin_mask.nii.gz

					fslmaths tracts_${a}/${tract}_bin_mask.nii.gz $filter -mas $ants_mask -bin -mas $ants_CSF_binv tracts_${a}/${tract}_filter_mask.nii.gz
					
					# here we add tckedit with the resulting thresholded prob mask for filtered bundles
					# trying tract mask as include with brain mask rather than using the tract mask as -mask
					if echo "$tract" | grep "_R_" > /dev/null 2>&1 ; then 
						echo "it's right sided"
						tckedit -number $nods -include roi/WM_fs_R.nii.gz $i $e $min_L -mask tracts_${a}/${tract}_filter_mask.nii.gz -nthreads $ncpu -force \
							tracts_${a}/${tract}.tck tracts_${a}/${tract}_filt.tck
					elif echo "$tract" | grep "_L_" > /dev/null 2>&1 ; then 
						 echo "it's left sided"
 						tckedit -number $nods -include roi/WM_fs_L.nii.gz $i $e $min_L -mask tracts_${a}/${tract}_filter_mask.nii.gz -nthreads $ncpu -force \
 							tracts_${a}/${tract}.tck tracts_${a}/${tract}_filt.tck
					fi
					
					# To acquire CST/ML separately we need an if loop to run tckedit with specific ROIs for subsegmentation
					# if the CST/ML is found initially.
					nods1=20000
					nods2=10000
					
					if [ -f tracts_${a}/CSTML_all_L_nods${nods1}.tck ]; then
						
						if [ ! -f tracts_${a}/CST_L_nods${nods2}.tck ]; then
							
							# separate the L_CST
							tckedit -number $nods2 -include roi/WM_fs_L.nii.gz -include roi/BStem.nii.gz -include roi/S1_wm_fs_L.nii.gz \
								-include roi/S1_wm_fs_L.nii.gz -include roi/LT_CST_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_R.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/LT_ML_JHU_roi_native.nii.gz \
										-exclude roi/RT_CST_JHU_roi_native.nii.gz -exclude roi/RT_ML_JHU_roi_native.nii.gz \
											$min_L $m -nthreads $ncpu -force tracts_${a}/CSTML_all_L_nods${nods1}.tck tracts_${a}/CST_L_nods${nods2}.tck
		                    tckmap tracts_${a}/CST_L_nods${nods2}.tck tracts_${a}/CST_L_nods${nods2}.nii.gz -template $ants_anat -force 
							# count and thr for L_CST
							count_CST_L=$(tckinfo tracts_${a}/CST_L_nods${nods2}.tck -count | grep count | head -n 1 | awk '{print $(NF)}')
									# setting a relative thr for tract filtering using 0.2% as thr cutoff
									act_count_CST_L="$(echo $count_CST_L | cut -d':' -f2)"
									thr_CST_L="$(bc <<< "scale = 2; (($act_count_CST_L*0.2/100))")"
							fslmaths tracts_${a}/CST_L_nods${nods2}.nii.gz -thr $thr_CST_L -s 1.5 -thr 0.2 -bin tracts_${a}/CST_L_nods${nods2}_bin_mask.nii.gz
							fslmaths tracts_${a}/CST_L_nods${nods2}_bin_mask.nii.gz -add roi/BStem.nii.gz -add roi/M1_wm_fs_L.nii.gz \
								-add roi/S1_wm_fs_L.nii.gz -add roi/LT_CST_JHU_roi_native.nii.gz -mas $ants_mask -bin -mas $ants_CSF_binv \
									tracts_${a}/CST_L_nods${nods2}_filter_mask.nii.gz
							# filter L_CST
							tckedit -number ${nods2} -include roi/WM_fs_L.nii.gz -include roi/BStem.nii.gz -include roi/M1_wm_fs_L.nii.gz \
								-include roi/S1_wm_fs_L.nii.gz -include roi/LT_CST_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_R.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/LT_ML_JHU_roi_native.nii.gz \
										-exclude roi/RT_CST_JHU_roi_native.nii.gz -exclude roi/RT_ML_JHU_roi_native.nii.gz \
											-mask tracts_${a}/CST_L_nods${nods2}_filter_mask.nii.gz $min_L -nthreads $ncpu \
												-force tracts_${a}/CST_L_nods${nods2}.tck tracts_${a}/CST_L_nods${nods2}_filt.tck
						else
				            
							echo "  CST_L already done , skipping..."
							
						fi
							
						if [ ! -f tracts_${a}/ML_L_nods${nods2}.tck ]; then
							
							# separate the L_ML
							tckedit -number ${nods2} -include roi/WM_fs_L.nii.gz -include roi/BStem.nii.gz -include roi/M1_wm_fs_L.nii.gz \
								-include roi/S1_wm_fs_L.nii.gz -include roi/LT_ML_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_R.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/LT_CST_JHU_roi_native.nii.gz \
										-exclude roi/RT_ML_JHU_roi_native.nii.gz -exclude roi/RT_CST_JHU_roi_native.nii.gz \
											$min_L $m -nthreads $ncpu -force tracts_${a}/CSTML_all_L_nods${nods1}.tck tracts_${a}/ML_L_nods${nods2}.tck
							# count and thr for L_ML
							count_ML_L=$(tckinfo tracts_${a}/ML_L_nods${nods2}.tck -count | grep count | head -n 1 | awk '{print $(NF)}')
									# setting a relative thr for tract filtering using 0.2% as thr cutoff
									act_count_ML_L="$(echo $count_ML_L | cut -d':' -f2)"
									thr_ML_L="$(bc <<< "scale = 2; (($act_count_ML_L*0.2/100))")"
		                    tckmap tracts_${a}/ML_L_nods${nods2}.tck tracts_${a}/ML_L_nods${nods2}.nii.gz -template $ants_anat -force 
							fslmaths tracts_${a}/ML_L_nods${nods2}.nii.gz -thr $thr_ML_L -s 1.5 -thr 0.2 -bin tracts_${a}/ML_L_nods${nods2}_bin_mask.nii.gz
							fslmaths tracts_${a}/ML_L_nods${nods2}_bin_mask.nii.gz -add roi/BStem.nii.gz -add roi/M1_wm_fs_L.nii.gz \
								-add roi/S1_wm_fs_L.nii.gz -add roi/LT_ML_JHU_roi_native.nii.gz -mas $ants_mask -bin -mas $ants_CSF_binv \
									tracts_${a}/ML_L_nods${nods2}_filter_mask.nii.gz
							# filter L_ML
							tckedit -number $nods2 -include roi/WM_fs_L.nii.gz -include roi/BStem.nii.gz -include roi/M1_wm_fs_L.nii.gz \
								-include roi/S1_wm_fs_L.nii.gz -include roi/LT_ML_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_R.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/LT_CST_JHU_roi_native.nii.gz \
										-exclude roi/RT_ML_JHU_roi_native.nii.gz -exclude roi/RT_CST_JHU_roi_native.nii.gz \
											-mask tracts_${a}/ML_L_nods${nods2}_filter_mask.nii.gz $min_L -nthreads $ncpu \
												-force tracts_${a}/ML_L_nods${nods2}.tck tracts_${a}/ML_L_nods${nods2}_filt.tck
						else
								
							echo "  ML_L already done , skipping..."
								
						fi
							
					else
							
						echo "  CST/ML_L not yet generated , try later"
							
							
					fi
						
							
					if [ -f tracts_${a}/CSTML_all_R_nods${nods1}.tck ]; then
							
						if [ ! -f tracts_${a}/CST_R_nods${nods2}.tck ]; then
								
							# separate the R_CST
							tckedit -number ${nods2} -include roi/WM_fs_R.nii.gz -include roi/BStem.nii.gz -include roi/M1_wm_fs_R.nii.gz \
								-include roi/S1_wm_fs_R.nii.gz -include roi/RT_CST_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_L.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/RT_ML_JHU_roi_native.nii.gz \
										-exclude roi/LT_CST_JHU_roi_native.nii.gz -exclude roi/LT_ML_JHU_roi_native.nii.gz \
											$min_L $m -nthreads $ncpu -force tracts_${a}/CSTML_all_R_nods${nods1}.tck tracts_${a}/CST_R_nods${nods2}.tck
							# count and thr for R_CST
							count_CST_R=$(tckinfo tracts_${a}/CST_R_nods${nods2}.tck -count | grep count | head -n 1 | awk '{print $(NF)}')
									# setting a relative thr for tract filtering using 0.2% as thr cutoff
									act_count_CST_R="$(echo $count_CST_R | cut -d':' -f2)"
									thr_CST_R="$(bc <<< "scale = 2; (($act_count_CST_R*0.2/100))")"
		                    tckmap tracts_${a}/CST_R_nods${nods2}.tck tracts_${a}/CST_R_nods${nods2}.nii.gz -template $ants_anat -force 
							fslmaths tracts_${a}/CST_R_nods${nods2}.nii.gz -thr $thr_CST_R -s 1.5 -thr 0.2 -bin tracts_${a}/CST_R_nods${nods2}_bin_mask.nii.gz
							fslmaths tracts_${a}/CST_R_nods${nods2}_bin_mask.nii.gz -add roi/BStem.nii.gz -add roi/M1_wm_fs_R.nii.gz \
								-add roi/S1_wm_fs_R.nii.gz -add roi/RT_CST_JHU_roi_native.nii.gz -mas $ants_mask -bin -mas $ants_CSF_binv \
									tracts_${a}/CST_R_nods${nods2}_filter_mask.nii.gz
							# filter R_CST
							tckedit -number ${nods2} -include roi/WM_fs_R.nii.gz -include roi/BStem.nii.gz -include roi/M1_wm_fs_R.nii.gz \
								-include roi/S1_wm_fs_R.nii.gz -include roi/RT_CST_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_L.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/RT_ML_JHU_roi_native.nii.gz \
										-exclude roi/LT_CST_JHU_roi_native.nii.gz -exclude roi/LT_ML_JHU_roi_native.nii.gz \
											-mask tracts_${a}/CST_R_nods${nods2}_filter_mask.nii.gz $min_L -nthreads $ncpu \
												-force tracts_${a}/CST_R_nods${nods2}.tck tracts_${a}/CST_R_nods${nods2}_filt.tck

						else
				            
							echo "  CST_R already done , skipping..."
							
						fi
						

						# ML is only to the primary somatosensory cortex, so no need to use M1.
						if [ ! -f tracts_${a}/ML_R_nods${nods2}.tck ]; then
							
							# separate the R_ML
							tckedit -number ${nods2} -include roi/WM_fs_R.nii.gz -include roi/BStem.nii.gz -include roi/M1_wm_fs_R.nii.gz \
								-include roi/S1_wm_fs_R.nii.gz -include roi/RT_ML_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_L.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/RT_CST_JHU_roi_native.nii.gz \
										-exclude roi/LT_ML_JHU_roi_native.nii.gz -exclude roi/LT_CST_JHU_roi_native.nii.gz \
											$min_L $m -nthreads $ncpu -force tracts_${a}/CSTML_all_R_nods${nods1}.tck tracts_${a}/ML_R_nods${nods2}.tck
							# count and thr for R_ML
							count_ML_R=$(tckinfo tracts_${a}/ML_R_nods${nods2}.tck -count | grep count | head -n 1 | awk '{print $(NF)}')
									# setting a relative thr for tract filtering using 0.2% as thr cutoff
									act_count_ML_R="$(echo $count_ML_R | cut -d':' -f2)"
									thr_ML_R="$(bc <<< "scale = 2; (($act_count_ML_R*0.2/100))")"
		                    tckmap tracts_${a}/ML_R_nods${nods2}.tck tracts_${a}/ML_R_nods${nods2}.nii.gz -template $ants_anat -force 
							fslmaths tracts_${a}/ML_R_nods${nods2}.nii.gz -thr $thr_ML_R -s 1.5 -thr 0.2 -bin tracts_${a}/ML_R_nods${nods2}_bin_mask.nii.gz
							fslmaths tracts_${a}/ML_R_nods${nods2}_bin_mask.nii.gz -add roi/BStem.nii.gz -add roi/M1_wm_fs_R.nii.gz \
								-add roi/S1_wm_fs_R.nii.gz -add roi/RT_ML_JHU_roi_native.nii.gz -mas $ants_mask -bin -mas $ants_CSF_binv \
									tracts_${a}/ML_R_nods${nods2}_filter_mask.nii.gz
							# filter R_ML
							tckedit -number ${nods2} -include roi/WM_fs_R.nii.gz -include roi/BStem.nii.gz -include roi/M1_wm_fs_R.nii.gz \
								-include roi/S1_wm_fs_R.nii.gz -include roi/RT_ML_JHU_roi_native.nii.gz \
									-exclude roi/WM_fs_L.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/RT_CST_JHU_roi_native.nii.gz \
										-exclude roi/LT_ML_JHU_roi_native.nii.gz -exclude roi/LT_CST_JHU_roi_native.nii.gz \
											-mask tracts_${a}/ML_R_nods${nods2}_filter_mask.nii.gz $min_L -nthreads $ncpu \
												-force tracts_${a}/ML_R_nods${nods2}.tck tracts_${a}/ML_R_nods${nods2}_filt.tck
							
						else
			            
							echo "  ML_R already done , skipping..."
						
						fi
							
							# fslmaths tracts_${a}/${tract}.nii.gz -thr $thr -s 2 -thr 0.1 -bin tracts_${a}/${tract}_bin_mask.nii.gz
							# tckedit -number $nods -include roi/WM_fs_L.nii.gz -include roi/BStem.nii.gz -include roi/S1_wm_fs_L.nii.gz \
							# 	-include roi/S1_wm_fs_L.nii.gz -include roi/LT_CST_JHU_roi_native.nii.gz \
							# 		-exclude roi/WM_fs_R.nii.gz -exclude roi/CC_fs_all.nii.gz -exclude roi/LT_ML_JHU_roi_native.nii.gz \
							# 			-exclude roi/RT_CST_JHU_roi_native.nii.gz -exclude roi/RT_ML_JHU_roi_native.nii.gz \
							# 				-mask tracts_${a}/CSTML_all_L_nods$nods.tck -nthreads $ncpu -force tracts_${a}/CSTML_all_L_nods$nods.tck tracts_${a}/$CST_L_nods$nods.tck
							# tracts_${a}/CST_all_L_nods$nods.tck tracts_${a}/${tract}_filt.tck
						# Then we need the filtering step using the fiber density map
						
					else
						
						echo "  CST/ML_R not yet generated , try later"
						
						
					fi

            fi

        else
        
            echo "  tcksift & generation subject/MNI space images already done, skipping..."
        
        fi
    
    done

}

wmfod=response/wmfod_reg2T1w.mif
dwi_preproced=dwi_preproced_reg2T1w.mif
fs_5tt=5tt/5ttseg.mif
gmwmi=5tt/5tt2gmwmi.nii.gz
ants_CSF_binv=roi/CSF_ants_binv.nii.gz


# Make an empty log file with information about the tracts
echo "subject, algorithm, tract, count" > tracts_FT_info.csv

# Trackings for S61759 (radwan)
# will divide the trackings in two, the first will not use -stop, the second will.
# we keep tract the same but add theta and -stop as another var

# AF
# we'll use the IFG Pt + Pop and the SFG as seeds
# CC, BStem and contralateral WM as excludes
# for further refinement use lobe wm labels
# Using only WM labels for AF
#
nods=5000
tract="AF_R_nods${nods}"
seeds=("IFG_POp_wm_fs_R" "IFG_PTr_wm_fs_R" "STG_wm_fs_R" "MTG_wm_fs_R" "SMG_wm_fs_R")
exclude=("WM_fs_L" "BStem" "CC_fs_all" "Ins_fs_R" "putamen_dil_fs_R" "CSF_ants")
theta=60
stop=()
act=()
kul_mrtrix_FT

nods=5000
tract="AF_L_nods${nods}"
seeds=("IFG_POp_wm_fs_L" "IFG_PTr_wm_fs_L" "STG_wm_fs_L" "MTG_wm_fs_L" "SMG_wm_fs_L")
exclude=("WM_fs_R" "BStem" "CC_fs_all" "Ins_fs_L" "putamen_dil_fs_L" "CSF_ants")
theta=60
stop=()
act=()
kul_mrtrix_FT

# nods=8000
# tract="AF_L_GM_nods${nods}"
# seeds=("IFG_POp_wm_fs_L" "IFG_POp_fs_L" "IFG_PTr_wm_fs_L" "IFG_PTr_fs_L" "STG_wm_fs_L" "STG_fs_L" "MTG_wm_fs_L" "SMG_wm_fs_L")
# exclude=("WM_fs_R" "BStem" "CC_fs_all" "Ins_fs_L" "putamen_dil_fs_L")
# theta=60
# stop=()
# act=()
# kul_mrtrix_FT

nods=5000
tract="AF_R_STG_only_nods${nods}"
seeds=("IFG_POp_wm_fs_R" "IFG_PTr_wm_fs_R" "STG_wm_fs_R")
exclude=("WM_fs_L" "BStem" "CC_fs_all" "Ins_fs_R" "putamen_dil_fs_R" "THALAMUS_fs_R" "CSF_ants")
theta=60
stop=()
act=()
kul_mrtrix_FT

nods=5000
tract="AF_L_STG_only_nods${nods}"
seeds=("IFG_POp_wm_fs_L" "IFG_PTr_wm_fs_L" "STG_wm_fs_L")
exclude=("WM_fs_R" "BStem" "CC_fs_all" "Ins_fs_L" "putamen_dil_fs_L" "THALAMUS_fs_L" "CSF_ants")
theta=60
stop=()
act=()
kul_mrtrix_FT

# CST
# we'll use the S1 and M1 + BStem as seeds
# contrateral cerebral WM and CC as excludes
nods=20000
tract="CSTML_all_R_nods${nods}"
seeds=("BStem" "S1_wm_fs_R" "M1_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "Amyg_fs_R")
theta=60
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=20000
tract="CSTML_all_L_nods${nods}"
seeds=("BStem" "S1_wm_fs_L" "M1_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "Amyg_fs_L")
theta=60
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# nods=8000
# tract="pcCST_L_nods${nods}"
# seeds=("WM_fs_L" "BStem" "PCent_wm_fs_L")
# exclude=("WM_fs_R" "CC_fs_all")
# theta=60
# stop=()
# act=()
# kul_mrtrix_FT
#
# nods=8000
# tract="pcCST_R_nods${nods}"
# seeds=("WM_fs_R" "BStem" "PCent_wm_fs_R")
# exclude=("WM_fs_L" "CC_fs_all")
# theta=60
# stop=()
# act=()
# kul_mrtrix_FT

# SMA_PMC
nods=10000
tract="SMA_PMC_R_nods${nods}"
seeds=("SMA_and_PMC_fs_R" "BStem")
exclude=("WM_fs_L" "CC_fs_all")
theta=50
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=10000
tract="SMA_PMC_L_nods${nods}"
seeds=("SMA_and_PMC_fs_L" "BStem")
exclude=("WM_fs_R" "CC_fs_all")
theta=50
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# FAT
nods=6000
tract="FAT_L_nods${nods}"
seeds=("IFG_POp_wm_fs_L" "SFG_fs_L")
exclude=("IFG_PTr_wm_fs_L" "M1_fs_L" "WM_fs_R" "CC_fs_all" "Ins_fs_L" "putamen_dil_fs_L")
theta=60
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=6000
tract="FAT_R_nods${nods}"
seeds=("IFG_POp_wm_fs_R" "SFG_fs_R")
exclude=("IFG_PTr_wm_fs_R" "M1_fs_R" "WM_fs_L" "CC_fs_all" "Ins_fs_R" "putamen_dil_fs_R")
theta=60
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT
#
# ATR
# exclude the hippocampi, insulae, par and temp wm as well as vDC
# try using -stop and -act
nods=3000
tract="ATR_R_nods${nods}"
seeds=("THALAMUS_fs_R" "L_OF_wm_fs_R")
exclude=("BStem" "WM_fs_L" "CC_fs_all" "Temp_wm_fs_R" "putamen_dil_fs_R" "Ins_wm_fs_R" "vDC_fs_R")
theta=50
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=3000
tract="ATR_L_nods${nods}"
seeds=("THALAMUS_fs_L" "L_OF_wm_fs_L")
exclude=("BStem" "WM_fs_R" "CC_fs_all" "Temp_wm_fs_L" "putamen_dil_fs_L" "Ins_wm_fs_L" "vDC_fs_L")
theta=50
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# STR
# needs occ wm and vDC as excludes still
nods=3000
tract="STR_R_nods${nods}"
seeds=("THALAMUS_fs_R" "Par_wm_fs_R" "SPL_fs_R" "S1_fs_R")
exclude=("BStem" "WM_fs_L" "CC_fs_all" "Temp_wm_fs_R" "Front_wm_fs_R" "vDC_fs_L" "putamen_dil_fs_R" "accumbens_dil_fs_R")
theta=50
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=3000
tract="STR_L_nods${nods}"
seeds=("THALAMUS_fs_L" "Par_wm_fs_L" "SPL_fs_L" "S1_fs_L")
exclude=("BStem" "WM_fs_R" "CC_fs_all" "Temp_wm_fs_L" "Front_wm_fs_L" "vDC_fs_L" "putamen_dil_fs_L" "accumbens_dil_fs_L")
theta=50
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# OT
nods=1000
tract="OT_L_nods${nods}"
seeds=("THALAMUS_fs_L" "OC_fs")
exclude=("WM_fs_R" "CC_fs_all" "BStem" "Amyg_fs_L" "THALAMUS_fs_R" "Caudate_fs_L" "iPCC_wm_fs_L" "CSF_ants" "GM_ants")
theta=45
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=1000
tract="OT_R_nods${nods}"
seeds=("THALAMUS_fs_R" "OC_fs")
exclude=("WM_fs_L" "CC_fs_all" "BStem" "Amyg_fs_R" "THALAMUS_fs_L" "Caudate_fs_R" "iPCC_wm_fs_R" "CSF_ants" "GM_ants")
theta=45
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
# act=$(printf " -act %s" "${fs_5tt[@]}")
# another option here is to try with gm_ants as exclude
kul_mrtrix_FT

# OR
nods=5000
tract="OR_L_nods${nods}"
seeds=("THALAMUS_fs_L" "periCalc_fs_L" "Occ_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "vDC_fs_L" "BStem" "THALAMUS_fs_R" "Caudate_fs_L" "Lat_V_fs_L" "iPCC_wm_fs_L" "CSF_ants")
theta=45
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=5000
tract="OR_R_nods${nods}"
seeds=("THALAMUS_fs_R" "periCalc_fs_R" "Occ_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "vDC_fs_R" "BStem" "THALAMUS_fs_L" "Caudate_fs_R" "Lat_V_fs_R" "iPCC_wm_fs_R" "CSF_ants")
theta=60
stop=$(printf " -backtrack -crop_at_gmwmi")
act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# IFOF
# now trying without the temporal lobe as exclude
nods=6000
tract="IFOF_R_nods${nods}"
seeds=("Occ_wm_fs_R" "L_OF_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "BStem" "THALAMUS_fs_R" "Caudate_fs_R" "Phip_wm_fs_R")
theta=60
stop=()
act=()
kul_mrtrix_FT

nods=6000
tract="IFOF_L_nods${nods}"
seeds=("Occ_wm_fs_L" "L_OF_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "BStem" "THALAMUS_fs_L" "Caudate_fs_L" "Phip_wm_fs_L")
theta=60
stop=()
act=()
kul_mrtrix_FT

# ILF
# needs ITG as the include ROI in the temporal lobe
nods=6000
tract="ILF_R_nods${nods}"
seeds=("Occ_wm_fs_R" "ITG_wm_fs_R" "ITG_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "Front_wm_fs_R" "THALAMUS_fs_R" "Caudate_fs_R" "unseg_wm_fs_R" "Phip_wm_fs_R")
theta=45
stop=()
act=()
kul_mrtrix_FT

nods=6000
tract="ILF_L_nods${nods}"
seeds=("Occ_wm_fs_L" "ITG_wm_fs_L" "ITG_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "Front_wm_fs_L" "THALAMUS_fs_L" "Caudate_fs_L" "unseg_wm_fs_L" "Phip_wm_fs_L")
theta=45
stop=()
act=()
kul_mrtrix_FT

# MLF
# need to also add MLF with MTG as the include of temporal lobe.
nods=6000
tract="MLF_R_nods${nods}"
seeds=("Occ_wm_fs_R" "MTG_wm_fs_R" "MTG_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "Front_wm_fs_R" "THALAMUS_fs_R" "Caudate_fs_R" "unseg_wm_fs_R" "Phip_wm_fs_R")
theta=45
stop=()
act=()
kul_mrtrix_FT

nods=6000
tract="MLF_L_nods${nods}"
seeds=("Occ_wm_fs_L" "MTG_wm_fs_L" "MTG_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "Front_wm_fs_L" "THALAMUS_fs_L" "Caudate_fs_L" "unseg_wm_fs_L" "Phip_wm_fs_L")
theta=45
stop=()
act=()
kul_mrtrix_FT

# SLF
nods=6000
tract="SLF_R_nods${nods}"
seeds=("SMG_wm_fs_R" "IFG_POp_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "Putamen_dil_fs_R" "THALAMUS_fs_R" "Caudate_fs_R" "unseg_wm_fs_R" "Phip_wm_fs_R")
theta=45
stop=-stop
act=()
kul_mrtrix_FT

nods=6000
tract="SLF_L_nods${nods}"
seeds=("SMG_wm_fs_L" "IFG_POp_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "Putamen_dil_fs_L" "THALAMUS_fs_L" "Caudate_fs_L" "unseg_wm_fs_L" "Phip_wm_fs_L")
theta=45
stop=-stop
act=()
kul_mrtrix_FT

# UF
nods=5000
tract="UF_L_nods${nods}"
seeds=("L_OF_fs_L" "L_OF_wm_fs_L" "TP_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "BStem" "THALAMUS_fs_L" "Caudate_fs_L" "CSF_ants" "unseg_wm_fs_L" "accumbens_dil_fs_L")
theta=60
stop=()
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=5000
tract="UF_R_nods${nods}"
seeds=("L_OF_fs_R" "L_OF_wm_fs_R" "TP_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "BStem" "THALAMUS_fs_R" "Caudate_fs_R" "CSF_ants" "unseg_wm_fs_R" "accumbens_dil_fs_R")
theta=60
stop=()
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# TIF
nods=3000
tract="TIF_L_nods${nods}"
seeds=("Ins_fs_L" "Ins_wm_fs_L" "TP_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "BStem" "THALAMUS_fs_R" "Front_wm_fs_L" "CSF_ants" "unseg_wm_fs_L" "Amyg_fs_L")
theta=60
stop=()
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=3000
tract="TIF_R_nods${nods}"
seeds=("Ins_fs_R" "Ins_wm_fs_R" "TP_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "BStem" "THALAMUS_fs_R" "Front_wm_fs_R" "CSF_ants" "unseg_wm_fs_R" "Amyg_fs_R")
theta=60
stop=()
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# AIF
nods=2000
tract="AIF_L_nods${nods}"
seeds=("Ins_fs_L" "Amyg_fs_L" "TP_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "BStem" "THALAMUS_fs_R" "CSF_ants" "unseg_wm_fs_L")
theta=60
stop=()
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=2000
tract="AIF_R_nods${nods}"
seeds=("Ins_fs_R" "Amyg_fs_R" "TP_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "BStem" "THALAMUS_fs_R" "CSF_ants" "unseg_wm_fs_R")
theta=60
stop=()
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# Cingulum
nods=6000
tract="cCing_R_nods${nods}"
seeds=("cACC_fs_R" "rACC_fs_R" "PCC_fs_R" "iPCC_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "unseg_wm_fs_R" "STG_wm_fs_R")
theta=50
stop=-stop
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=6000
tract="cCing_L_nods${nods}"
seeds=("cACC_fs_L" "rACC_fs_L" "PCC_fs_L" "iPCC_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "unseg_wm_fs_L" "STG_wm_fs_L")
theta=50
stop=-stop
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=4000
# wm post. cing needs to be added here
tract="pCing_R_nods${nods}"
seeds=("Hippo_fs_R" "iPCC_fs_R" "Phip_wm_fs_R")
exclude=("WM_fs_L" "CC_fs_all" "unseg_wm_fs_R" "PCC_wm_fs_R" "fusiform_wm_fs_R")
theta=50
stop=-stop
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

nods=4000
tract="pCing_L_nods${nods}"
seeds=("Hippo_fs_L" "iPCC_fs_L" "Phip_wm_fs_L")
exclude=("WM_fs_R" "CC_fs_all" "unseg_wm_fs_L" "PCC_wm_fs_L" "fusiform_wm_fs_L")
theta=50
stop=-stop
act=()
# act=$(printf " -act %s" "${fs_5tt[@]}")
kul_mrtrix_FT

# Now prepare the data for iPlan
# if [ ! -f for_iplan/TH_SMAPMC_R.hdr ]; then
#
#     mkdir -p for_iplan
#
#     smooth_sigma=0.6
#
#     # copy the tracts in analyze format
#     fslmaths tracts_iFOD2/TH-DR_L_nods${nods} -s $smooth_sigma -thr 3 -bin for_iplan/Tract_DRT_L
#     fslchfiletype NIFTI_PAIR for_iplan/Tract_DRT_L for_iplan/Tract_DRT_L
#
#     fslmaths tracts_iFOD2/TH-DR_R_nods${nods} -s $smooth_sigma -thr 3 -bin for_iplan/Tract_DRT_R
#     fslchfiletype NIFTI_PAIR for_iplan/Tract_DRT_R for_iplan/Tract_DRT_R
#
#     # copy the T1w in analyze format
#     cp T1w/T1w_BrainExtractionBrain.nii.gz for_iplan/anat.nii.gz
#     fslchfiletype NIFTI_PAIR for_iplan/anat for_iplan/anat
#
#     # copy the Thalamic probabilistic images as speudo fmri activation maps
#     fslmaths tracts_iFOD2/Subj_Space_TH-DR_L_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_DRT_L
#     fslchfiletype NIFTI_PAIR for_iplan/TH_DRT_L for_iplan/TH_DRT_L
#     fslmaths tracts_iFOD2/Subj_Space_TH-DR_R_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_DRT_R
#     fslchfiletype NIFTI_PAIR for_iplan/TH_DRT_R for_iplan/TH_DRT_R
#
#     fslmaths tracts_iFOD2/Subj_Space_TH-M1_fs_L_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_M1_L
#     fslchfiletype NIFTI_PAIR for_iplan/TH_M1_L for_iplan/TH_M1_L
#     fslmaths tracts_iFOD2/Subj_Space_TH-M1_fs_R_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_M1_R
#     fslchfiletype NIFTI_PAIR for_iplan/TH_M1_R for_iplan/TH_M1_R
#
#     fslmaths tracts_iFOD2/Subj_Space_TH-S1_fs_L_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_S1_L
#     fslchfiletype NIFTI_PAIR for_iplan/TH_S1_L for_iplan/TH_S1_L
#     fslmaths tracts_iFOD2/Subj_Space_TH-S1_fs_R_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_S1_R
#     fslchfiletype NIFTI_PAIR for_iplan/TH_S1_R for_iplan/TH_S1_R
#
#     fslmaths tracts_iFOD2/Subj_Space_TH-SMA_and_PMC_L_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_SMAPMC_L
#     fslchfiletype NIFTI_PAIR for_iplan/TH_SMAPMC_L for_iplan/TH_SMAPMC_L
#     fslmaths tracts_iFOD2/Subj_Space_TH-SMA_and_PMC_R_nods${nods}_iFOD2.nii.gz -s $smooth_sigma -thr 0.3 for_iplan/TH_SMAPMC_R
#     fslchfiletype NIFTI_PAIR for_iplan/TH_SMAPMC_R for_iplan/TH_SMAPMC_R
#
#     # clean up
#     #rm -rf for_iplan/*.nii.gz
#
#fi

# ---- END BIG LOOP for processing each session
done

# write a file to indicate that dwiprep_drtdbs runned succesfully
#   this file will be checked by KUL_preproc_all

echo "done" > ../dwiprep_FT_is_done.log


kul_e2cl "   done KUL_dwiprep_drtdbs on participant $BIDS_participant" $log
