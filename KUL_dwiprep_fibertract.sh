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

    # sift1 filtering
    # termination ratio - defined as the ratio between reduction in cost
    # function, and reduction in density of streamlines.
    # Smaller values result in more streamlines being filtered out.
    do_sift_th=10000 # when to do sift? (if more than 5000 streamlines in tract e.g.)
    term_ratio=0.5 # reduce by e.g. 50%

    # tmp directory for temporary processing
    tmp=/tmp

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

     -f:  perform whole/full brain fibertractography first (and tckedit)
     -w:  which wmfod to use (default = dhollander_wmfod_norm_reg2T1w)
     -s:  session (of the participant)
     -n:  number of cpu for parallelisation
     -v:  show output from mrtrix commands


USAGE

    exit 1
}

function kul_mrtrix_tracto {

    if [ $f_flag -eq 1 ]; then
        
        local a=${algorithm}_WBFT
    
    else

        local a=$algorithm
    
    fi

    local d=tracts_${a}_from_${wmfod_select}
    
    mkdir -p ${d}
    
    kul_e2cl " running tckgen of ${tract} tract with algorithm $a all seeds with parameters $parameters" ${log}

    # make the seed string
    local s=$(printf " -seed_image roi/%s.nii.gz"  ${seeds[@]})
    
    # make the include string (which is same rois as seed)
    local i=$(printf " -include roi/%s.nii.gz"  ${include[@]})

    # make the exclude string (which is same rois as seed)
    local e=$(printf " -exclude roi/%s.nii.gz"  ${exclude[@]})

    # make the mask string 
    local m="-mask dwi_preproced_reg2T1w_mask.nii.gz"

  
    # do the tracking
            
    if [ ! -f ${d}/${tract}.tck ]; then 

        #echo ${a}            
        #pwd
        #echo $s
        #echo $i
        #echo $e
        #echo $m
        #echo $parameters

        if [ $f_flag -eq 1 ]; then

            #WBFT
            tckedit $i $e $m ${output_wbft}.tck ${d}/${tract}.tck -nthreads $ncpu -force

        else

            if [[ "${a}" =~ "iFOD" ]]; then

                # perform IFOD1 or 2 tckgen
                tckgen $wmfod ${d}/${tract}.tck -algorithm $a $parameters $s $i $e $m -nthreads $ncpu -force

            elif [[ "${a}" =~ "Tensor" ]]; then

                # perform Tensor_Prob or Tensor_Det tckgen
                tckgen $dwi_preproced ${d}/${tract}.tck -algorithm $a $parameters $s $i $e $m -nthreads $ncpu -force

            fi
        
        fi

    else

        echo "  tckgen of ${tract} tract already done, skipping"

    fi


    # Check if any fibers have been found & log to the information file
    echo "   checking ${d}/${tract}"
    local count=$(tckinfo ${d}/${tract}.tck | grep count | head -n 1 | awk '{print $(NF)}')
    echo "$subj, $a, $tract, $count" >> tracts_info.csv

        # do further processing of tracts are found
        if [ ! -f ${d}/MNI_Space_${tract}_${a}.nii.gz ]; then

            if [ $count -eq 0 ]; then

                # report that no tracts were found and stop further processing
                kul_e2cl "  no streamlines were found for the ${d}/${tract}.tck" ${log}

            else

                # report how many tracts were found and continue processing
                echo "   $count streamlines were found for the ${d}/${tract}.tck"
                
                echo "   generating subject/MNI space images"
                # convert the tck in nii
                tckmap ${d}/${tract}.tck ${d}/${tract}.nii.gz -template $ants_anat -force 

                # Warp the full tract image to MNI space
                input=${d}/${tract}.nii.gz
                output=${d}/MNI_Space_${tract}_${a}.nii.gz
                transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
                reference=$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz
                KUL_antsApply_Transform

                # make a probabilistic image in subject and MNI space
                local m=$(mrstats -quiet ${d}/${tract}.nii.gz -output max)
                #echo $m
                fslmaths ${d}/${tract} -div $m ${d}/Subj_Space_prob_${tract}_${a}
                fslmaths ${d}/MNI_Space_${tract}_${a}.nii.gz -div $m ${d}/MNI_Space_prob_${tract}_${a}

                # Warp the probabilistic image to MNI space
                #input=${d}/Subj_Space_${tract}_${a}.nii.gz
                #output=${d}/MNI_Space_${tract}_${a}.nii.gz
                #transform=${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
                #reference=$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz
                #KUL_antsApply_Transform
            
                # Make a smoothed version of the tracts for iPlan
                smooth_sigma=0.6
                prob_treshold=0.15
                fslmaths ${d}/Subj_Space_prob_${tract}_${a} -s $smooth_sigma -thr $prob_treshold ${d}/Subj_Space_prob_smooth_${tract}_${a}

            fi

        else
        
            echo "  tckshift & generation subject/MNI space images already done, skipping..."
        
        fi
    

}


function KUL_antsApply_Transform {

    # Fix bug in antsApplytransforms (add EOF at tranform file)

    cp $transform /tmp/transform_tmp.txt
    echo "" >> /tmp/transform_tmp.txt

    antsApplyTransforms -d 3 --float 1 \
    --verbose 1 \
    -i $input \
    -o $output \
    -r $reference \
    -t /tmp/transform_tmp.txt \
    -n Linear

    rm -rf /tmp/transform_tmp.txt
    
}

# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
ncpu=6
silent=1
wmfod_select=dhollander_wmfod_norm_reg2T1w
#wmfod_select=dhollander_wmfod_norm_noGM_reg2T1w

# Set required options
p_flag=0
c_flag=0
r_flag=0
s_flag=0
f_flag=0

if [ "$#" -lt 3 ]; then

    echo
    echo "Please specify all required options!"
    echo 

    Usage >&2
    exit 1

else

    while getopts "p:c:r:s:n:w:fvh" OPT; do

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
        f) #session
            f_flag=1
        ;;
        n) #parallel
            ncpu=$OPTARG
        ;;
        w) #wmfod
            wmfod_select=$OPTARG
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
fmriprep_subj=fmriprep/"sub-${subj}"

# Either a session is given on the command line
# If not the session(s) need to be determined.
if [ $s_flag -eq 1 ]; then

    # session is given on the command line
    search_sessions=BIDS/sub-${subj}/ses-${ses}
    xfm_search=($(find ${cwd}/${fmriprep_subj}/ses-${ses} -type f -name "*from-orig_to-T1w_mode-image_xfm*" ! -name "*gadolinium*"))
    num_xfm=${#xfm_search[@]}
    echo "  Xfm files: number : $num_xfm"
    echo "    notably: ${xfm_search[@]}"

else

    # search if any sessions exist
    search_sessions=($(find BIDS/sub-${subj} -type d | grep dwi))

    # Transforming the T1w to fmriprep space
    xfm_search=($(find ${cwd}/${fmriprep_subj} -type f -name "*from-orig_to-T1w_mode-image_xfm*" ! -name "*gadolinium*"))
    num_xfm=${#xfm_search[@]}
    echo "  Xfm files: number : $num_xfm"
    echo "    notably: ${xfm_search[@]}"

fi    



num_sessions=${#search_sessions[@]}
    
echo "  Number of BIDS sessions: $num_sessions"
echo "    notably: ${search_sessions[@]}"


# ---- BIG LOOP for processing each session
for current_session in `seq 0 $(($num_sessions-1))`; do

    # set up directories 
    cd $cwd
    long_bids_subj=${search_sessions[$current_session]}
    #echo $long_bids_subj
    bids_subj=${long_bids_subj%dwi}

    # Change the Directory to write preprocessed data in
    preproc=dwiprep/sub-${subj}/$(basename $bids_subj) 
    #echo $preproc
    cd $preproc

    # STEP 1 - create the ROIS for fibertractography -------------------------------------------------------
    kul_e2cl " Creating ROIS for $bids_subj from" ${log}

    # Where is the freesurfer parcellation? 
    fs_aparc=${cwd}/freesurfer/sub-${subj}/${subj}/mri/aparc+aseg.mgz

    # Where is the T1w anat?
    ants_anat=T1w/T1w_BrainExtractionBrain.nii.gz

    # Where is fs_labels?
    #fs_labels=roi/labels_from_FS.nii.gz
    #fs_wmlabels=roi/labels_wm_from_FS.nii.gz

    #fmriprep_anat="${cwd}/${fmriprep_subj}/anat/sub-${subj}_desc-preproc_T1w.nii.gz"

    # we read the config file (and it may be csv, tsv or ;-seperated)
    while IFS=$'\t,;' read -r roi_name from_atlas space label_id; do
 
        if [[ ! $roi_name == \#* ]]; then
        

            if [ ! -f roi/${roi_name}.nii.gz ]; then   

                if [ $space = "subject" ]; then

                    echo " creating the $space space $roi_name ROI from $from_atlas using label_id $label_id..." 
                    
                    labelArray=($label_id)
                    #echo ${labelArray[@]}
                    
                    for label_id_tmp in "${labelArray[@]}"
                    do
                        
                        #echo $label_id_tmp
                        atlas_file="$(echo -e "${from_atlas}" | sed -e 's/^[[:space:]]*//')"
                        echo "fslmaths roi/$atlas_file -thr $label_id_tmp -uthr $label_id_tmp -bin roi/${roi_name}_${label_id_tmp}"
                        fslmaths roi/$atlas_file -thr $label_id_tmp -uthr $label_id_tmp -bin roi/${roi_name}_${label_id_tmp}
                    
                    done 
                    fslmerge -t roi/${roi_name}_merged roi/${roi_name}_* 
                    fslmaths roi/${roi_name}_merged -Tmean -bin roi/${roi_name}

                elif [ $space = "mni" ]; then

                    echo " creating the $space space $roi_name ROI from $from_atlas..." 

                    input=${kul_main_dir}/atlasses/Local/${from_atlas}
                    input=${input//[[:blank:]]/}
                    output=roi/${roi_name}_tmp.nii.gz
                    transform="${cwd}/fmriprep/sub-${subj}/anat/sub-${subj}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
                    reference=$ants_anat
                    KUL_antsApply_Transform

                    input=roi/${roi_name}_tmp.nii.gz
                    output=roi/${roi_name}.nii.gz
                    transform=${xfm_search[$current_session]}
                    reference=$ants_anat
                    if [ -z "$transform" ];then
                        cp $input $output
                    else
                        KUL_antsApply_Transform
                    fi

                fi

            fi
        
        fi 

    done < ${cwd}/$rois_config


    # STEP 2 - perform fibertractography -------------------------------------------------------
    
    wmfod=response/${wmfod_select}.mif
    dwi_preproced=dwi_preproced_reg2T1w.mif
    dwi_mask=dwi_preproced_reg2T1w_mask.nii.gz
    #output_wbft=tracks_50_million
    #input_wbft=50000000
    
    output_wbft=WBFT_20_million_from_${wmfod_select}
    input_wbft=20000000
    
    
    # Check WBFT or direct tracking
    if [ $f_flag -eq 1 ]; then

        # WBFT
        wbft_options="-maxlen 250 -minlen 10 -select $input_wbft"

        echo " Performing Whole Brain Fiber Tractography first"

        if [ ! -f ${output_wbft}.tck ]; then
            tckgen $wmfod -seed_image $dwi_mask -mask $dwi_mask ${output_wbft}.tck $wbft_options -nthreads $ncpu
        fi

    fi

    # Make an empty log file with information about the tracts
    echo "subject, algorithm, tract, count" > tracts_info.csv

    # we read the config file (and it may be csv, tsv or ;-seperated)
    while IFS=$'\t,;' read -r tract_name seed_rois include_rois exclude_rois algorithm parameters; do

        #echo "tract_name    = $tract_name"

        if [[ ! $tract_name == \#* ]]; then
        
            tract=$tract_name
            seeds=($seed_rois)  
            include=($include_rois)
            exclude=($exclude_rois)
            kul_mrtrix_tracto
        
        fi 

    done < ${cwd}/$tracts_config


done

