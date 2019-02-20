#!/bin/bash -e
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
#  Project PI's: Stefan Sunaert & Bart Nuttin
#
# Requires Mrtrix3, FSL, ants, freesurfer
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 14/02/2019 - alpha version
v="v0.1 - dd 14/02/2019"

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

    # development for Rose Bruffaerts
    Rose=0

# 


# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl (for logging)

kul_main_dir=`dirname "$0"`
script=`basename "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs dMRI fibertractography.

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001 -c study_config/tracto_tracts.csv -r study_config/tracto_rois.csv 

Required arguments:

     -p:  participant (anonymised name of the subject)
     -c:  tractography config file (what tracts to generate & settings)
     -r:  tractography file with ROIs

Optional arguments:

     -s:  session (of the participant)
     -o:  number of desired streamlines to select in tckgen (default nods=2000)
     -n:  number of cpu for parallelisation
     -v:  show output from mrtrix commands


USAGE

    exit 1
}

function kul_mrtrix_tracto {

    for a in iFOD2 Tensor_Prob; do
    #for a in iFOD2; do
    
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


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
nods=2000
ncpu=6
silent=1

# Set required options
p_flag=0
c_flag=0
r_flag=0
s_flag=0

if [ "$#" -lt 3 ]; then

    echo
    echo "Please specify all required options!"
    echo 

    Usage >&2
    exit 1

else

    while getopts "p:c:r:s:o:n:vh" OPT; do

        case $OPT in
        p) #subject
            p_flag=1
            subj=$OPTARG
        ;;
        c) #tracto-config
            c_flag=1
            tracts_config=$OPTARG
        ;;
        r) #rois-config
            r_flag=1
            rois_config=$OPTARG
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
if [ $c_flag -eq 0 ] ; then 
    echo 
    echo "Option -c is required: give the config file with tractography settings." >&2
    echo
    exit 2 
fi 
if [ $r_flag -eq 0 ] ; then 
    echo 
    echo "Option -r is required: give the config file with rois to create." >&2
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

    # Change the Directory to write preprocessed data in
    preproc=dwiprep/sub-${subj}/$(basename $bids_subj) 
    echo $preproc
    cd $preproc

    # STEP 1 - create the ROIS for fibertractography -------------------------------------------------------
    kul_e2cl " Creating ROIS for $bids_subj from" ${log}

    # Where is the freesurfer parcellation? 
    fs_aparc=${cwd}/freesurfer/sub-${subj}/${subj}/mri/aparc+aseg.mgz

    # Where is the T1w anat?
    ants_anat=T1w/T1w_BrainExtractionBrain.nii.gz

    # Where is fs_labels?
    fs_labels=roi/labels_from_FS.nii.gz


    # we read the config file (and it may be csv, tsv or ;-seperated)
    while IFS=$'\t,;' read -r roi_name from_atlas space label_id; do
 
        if [ "$roi_name" = "roi_name" ]; then
        
        echo "first line" > /dev/null 2>&1

        else

            if [ $space = "subject" ]; then

                echo " creating the $space space $roi_name ROI from $from_atlas using label_id $label_id..." 

                fslmaths $fs_labels -thr $label_id -uthr $label_id -bin roi/${roi_name}

            fi
        
        fi 

    done < ${cwd}/$rois_config


    # STEP 2 - perform fibertractography -------------------------------------------------------
    wmfod=response/wmfod_reg2T1w.mif
    dwi_preproced=dwi_preproced_reg2T1w.mif

    # Make an empty log file with information about the tracts
    echo "subject, algorithm, tract, count" > tracts_info.csv

    # we read the config file (and it may be csv, tsv or ;-seperated)
    while IFS=$'\t,;' read -r tract_name seed_rois include_rois exclude_rois algorithm paramaters; do

        echo "tract_name    = $tract_name"

        if [ "$tract_name" = "tract_name" ]; then
        
        echo "first line" > /dev/null 2>&1

        else

            tract=$tract_name
            seeds=($seed_rois)  
            exclude=$exclude_rois
            kul_mrtrix_tracto
        
        fi 

    done < ${cwd}/$tracts_config




done

exit


# Extract relevant freesurfer determined rois
if [ ! -f roi/WM_fs_R.nii.gz ]; then

    

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
    # 1018   ctx-lh-parsopercularis
    # 2018   ctx-rh-parsopercularis
    fslmaths $fs_labels -thr 1018 -uthr 1018 -bin roi/IFGparsopercularis_fs_L
    fslmaths $fs_labels -thr 2018 -uthr 2018 -bin roi/IFGparsopercularis_fs_R

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

    for a in iFOD2 Tensor_Prob; do
    #for a in iFOD2; do
    
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

if [ $Rose -eq 1 ]; then
    
    tract="Aslant_L_nods${nods}"
    seeds=("IFGparsopercularis_fs_L" "SFG_fs_L")  
    exclude="WM_fs_R"
    kul_mrtrix_tracto_drt 

    tract="Aslant_R_nods${nods}"
    seeds=("IFGparsopercularis_fs_R" "SFG_fs_R")  
    exclude="WM_fs_L"
    kul_mrtrix_tracto_drt 

fi



# ---- END BIG LOOP for processing each session
done

# write a file to indicate that dwiprep_drtdbs runned succesfully
#   this file will be checked by KUL_preproc_all

echo "done" > ../${script}_is_done.log


kul_e2cl "   done $script on participant $BIDS_participant" $log
