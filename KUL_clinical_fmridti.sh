#!/bin/bash -e
# Bash shell script to analyse clinical fMRI/DTI
#
# Requires matlab fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 16/12/2020
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` is a batch analysis of clinical fMRI/DTI data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -d DICOM/JohnDoe.zip

Required arguments:

     -p:  participant name
     -t:  processing type
        type 1: do hd-glio and vbg (tumor with T1w, cT1w, T2w and FLAIR)
        type 2: do vbg with manual mask (tumor but missing one of T1w, cT1w, T2w and FLAIR)
        type 3: don't run hd-glio nor vbg (cavernoma, epilepsy, etc... cT1w)

Optional arguments:

     -d:  dicom zip file (or directory)
     -n:  number of cpu to use (default 15)
     -v:  show output from commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ncpu=15

# Set required options
p_flag=0
d_flag=0 
t_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:d:n:v" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        t) #type
			type=$OPTARG
            t_flag=1
		;;
        d) #dicomzip
			dicomzip=$OPTARG
            d_flag=1
		;;
        n) #ncpu
			ncpu=$OPTARG
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

fi

# check for required options
if [ $p_flag -eq 0 ] ; then
	echo
	echo "Option -p is required: give the BIDS name of the participant." >&2
	echo
	exit 2
fi

if [ $t_flag -eq 0 ] ; then
	echo
	echo "Option -t is required: give the analysis type." >&2
	echo
	exit 2
fi

# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
    str_silent=" > /dev/null 2>&1" 
fi

if [ $type -eq 1 ]; then
    hdglio=1
    vbg=1
elif [ $type -eq 2 ]; then
    hdglio=0
    vbg=1
elif [ $type -eq 3 ]; then
    hdglio=0
    vbg=0
fi

# --- functions ---

function KUL_antsApply_Transform {
    antsApplyTransforms -d 3 --float 1 \
        --verbose 1 \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n Linear
}

function KUL_convert2bids {
    # convert the DICOM to BIDS
    if [ ! -d "BIDS/sub-${participant}" ];then
        KUL_dcm2bids_new.sh -d $dicomzip -p ${participant} -c study_config/sequences.txt -e
    else
        echo "BIDS conversion already done"
    fi
}

function KUL_run_fmriprep {
    if [ ! -f fmriprep/sub-${participant}.html ]; then
        cp study_config/run_fmriprep.txt KUL_LOG/sub-${participant}_run_fmriprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_fmriprep.txt 
        rm -fr fmriprep_work_${participant}
    else
        echo "fmriprep already done"
    fi
}

function KUL_run_dwiprep {
    if [ ! -f dwiprep/sub-${participant}/dwiprep_is_done.log ]; then
        cp study_config/run_dwiprep.txt KUL_LOG/sub-${participant}_run_dwiprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_dwiprep.txt
        rm -f KUL_LOG/sub-${participant}_run_dwiprep.txt.bck
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_dwiprep.txt 
    else
        echo "dwiprep already done"
    fi
}

function KUL_run_freesurfer {
    if [ ! -f BIDS/derivatives/freesurfer/${participant}_freesurfer_is.done ]; then
        cp study_config/run_freesurfer.txt KUL_LOG/sub-${participant}_run_freesurfer.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_freesurfer.txt
        rm -f KUL_LOG/sub-${participant}_run_freesurfer.txt.bck
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_freesurfer.txt 
    else
        echo "freesurfer already done"
    fi
}

function KUL_compute_SPM_matlab {
    fmriresults="$computedir/RESULTS/stats_$shorttask"
    mkdir -p $fmriresults
    pcf="${scriptsdir}/stats_${shorttask}.m" #participant config file
    pjf="${scriptsdir}/stats_${shorttask}_job.m" #participant job file
    # get rid of - in filename, since this breaks -r in matlab
    pcf=${pcf/run-/run}
    pjf=${pjf/run-/run}
    #echo "$pcf -- $pjf"
    cp $tcf $pcf
    cp $tjf $pjf
    sed -i.bck "s|###JOBFILE###|$pjf|" $pcf
    sed -i.bck "s|###FMRIDIR###|$fmridatadir|" $pjf
    sed -i.bck "s|###FMRIFILE###|$fmrifile|" $pjf
    sed -i.bck "s|###FMRIRESULTS###|$fmriresults|" $pjf
    $matlab_exe -nodisplay -nosplash -nodesktop -r "run('$pcf');exit;"
            
    result=$computedir/RESULTS/MNI/${shorttask}_space-MNI152NLin6Asym.nii
    cp $fmriresults/spmT_0001.nii $result
            
    result_global=$cwd/RESULTS/sub-$participant/SPM_${shorttask}.nii
            
    # since SPM analysis was in MNI space, we transform back in native space
    input=$result
    output=$result_global
    transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5
    find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
    reference=${find_T1w[0]}
    #echo "input=$input"
    #echo "output=$output"
    #echo "transform=$transform"
    #echo "reference=$reference"
    KUL_antsApply_Transform

    gm_mask="$fmriprepdir/../anat/sub-${participant}_label-GM_probseg.nii.gz"
    gm_mask2=$computedir/RESULTS/gm_mask_${shorttask}.nii.gz
    gm_result_global=$cwd/RESULTS/sub-$participant/SPM_${shorttask}_gm.nii
    mrgrid $gm_mask regrid -template $result_global $gm_mask2
    gm_mask3=$computedir/RESULTS/gm_mask_smooth_${shorttask}.nii.gz
    mrfilter $gm_mask2 smooth $gm_mask3
    #mrcalc $result_global $gm_mask3 0.3 -gt -mul $gm_result_global

} 

function KUL_compute_SPM {
    #  setup variables
    computedir="$cwd/compute/SPM/sub-$participant"
    fmridatadir="$computedir/fmridata"
    scriptsdir="$computedir/scripts"
    fmriprepdir="fmriprep/sub-$participant/func"
    searchtask="_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii"
    matlab_exe=$(which matlab)
    #  the template files in KNT for SPM analysis
    tcf="$kul_main_dir/share/spm12/spm12_fmri_stats_1run.m" #template config file
    tjf="$kul_main_dir/share/spm12/spm12_fmri_stats_1run_job.m" #template job file

    mkdir -p $fmridatadir
    mkdir -p $scriptsdir
    mkdir -p $computedir/RESULTS/MNI

    # Provide the anatomy
    cp -f $fmriprepdir/../anat/sub-${participant}_desc-preproc_T1w.nii.gz $globalresultsdir/T1w.nii.gz
    gunzip -f $globalresultsdir/T1w.nii.gz

    if [ ! -f KUL_LOG/sub-${participant}_SPM.done ]; then
        echo "Preparing for SPM"
        tasks=( $(find $fmriprepdir -name "*${searchtask}.gz" -type f) )
        #echo ${tasks[@]}

        # we loop over the found tasks
        for task in ${tasks[@]}; do
            d1=${task#*_task-}
            shorttask=${d1%_space*}
            #echo "$task -- $shorttask"
            if [[ ! "$shorttask" = *"rest"* ]]; then
                echo " Analysing task $shorttask"
                fmrifile="${shorttask}${searchtask}"
                cp $fmriprepdir/*$fmrifile.gz $fmridatadir
                gunzip -f $fmridatadir/*$fmrifile.gz
                KUL_compute_SPM_matlab

                # do the combined analysis
                if [[ "$shorttask" == *"run-2" ]]; then
                    echo "this is run2, we run full analysis now" 
                    shorttask=${shorttask%_run-2}
                    #echo $shorttask
                    tcf="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs.m" #template config file
                    tjf="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs_job.m" #template job file
                    fmrifile="${shorttask}"
                    KUL_compute_SPM_matlab
                fi
            
            fi
        done
        touch KUL_LOG/sub-${participant}_SPM.done
    else
        echo "SPM analysis already done"
    fi
}

function KUL_segment_tumor {
    bidsdir="BIDS/sub-$participant"
    T1w=($(find $bidsdir -name "*T1w.nii.gz" ! -name "*gadolinium*" -type f ))
    nT1w=${#T1w[@]}
    echo "number of non-contrast T1w: $nT1w"
    cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
    ncT1w=${#cT1w[@]}
    echo "number of contrast enhanced T1w: $ncT1w"
    FLAIR=($(find $bidsdir -name "*FLAIR.nii.gz" -type f ))
    nFLAIR=${#FLAIR[@]}
    echo "number of FLAIR: $nFLAIR"
    T2w=($(find $bidsdir -name "*T2w.nii.gz" -type f ))
    nT2w=${#T2w[@]}
    echo "number of T2w: $nT2w"

    mkdir -p $globalresultsdir
    
    # this will segment the lesion automatically
    hdglioinputdir="compute/hdglio/sub-${participant}/input"
    hdgliooutputdir="compute/hdglio/sub-${participant}/output"
    if [ $hdglio -eq 1 ]; then

        if [ $nT1w -eq 1 ] && [ $ncT1w -eq 1 ] && [ $nFLAIR -eq 1 ] && [ $nT2w -eq 1 ];then
            mkdir -p $hdglioinputdir
            mkdir -p $hdgliooutputdir
            if [ ! -f "$hdgliooutputdir/volumes.txt" ]; then
                cp $T1w $hdglioinputdir/T1.nii.gz
                cp $cT1w $hdglioinputdir/CT1.nii.gz
                cp $FLAIR $hdglioinputdir/FLAIR.nii.gz
                cp $T2w $hdglioinputdir/T2.nii.gz
                if [[ $machine_type = "Darwin" ]];then
                    echo "Running HD-GLIO"
                    fslreorient2std $hdglioinputdir/T1.nii.gz $hdgliooutputdir/T1_reorient.nii.gz
                    fslreorient2std $hdglioinputdir/CT1.nii.gz $hdgliooutputdir/CT1_reorient.nii.gz
                    fslreorient2std $hdglioinputdir/T2.nii.gz $hdgliooutputdir/T2_reorient.nii.gz
                    fslreorient2std $hdglioinputdir/FLAIR.nii.gz $hdgliooutputdir/FLAIR_reorient.nii.gz
                    cd $hdgliooutputdir
                    # run hd bet
                    hd-bet -i T1_reorient.nii.gz -o t1_bet.nii.gz -s 1 -device cpu -mode fast -tta 0
                    hd-bet -i CT1_reorient.nii.gz -o ct1_bet.nii.gz -device cpu -mode fast -tta 0
                    hd-bet -i T2_reorient.nii.gz -o t2_bet.nii.gz -device cpu -mode fast -tta 0
                    hd-bet -i FLAIR_reorient.nii.gz -o flair_bet.nii.gz -device cpu -mode fast -tta 0
                    # register brain extracted images to t1, save matrix
                    flirt -in ct1_bet.nii.gz -out ct1_bet_reg.nii.gz -ref t1_bet.nii.gz -omat ct1_to_t1.mat -interp spline -dof 6 &
                    flirt -in t2_bet.nii.gz -out t2_bet_reg.nii.gz -ref t1_bet.nii.gz -omat t2_to_t1.mat -interp spline -dof 6 &
                    flirt -in flair_bet.nii.gz -out flair_bet_reg.nii.gz -ref t1_bet.nii.gz -omat flair_to_t1.mat -interp spline -dof 6 &
                    wait
                    # we are only interested in the matrices, delete the other output images
                    rm ct1_bet.nii.gz t2_bet.nii.gz flair_bet.nii.gz
                    rm ct1_bet_reg.nii.gz t2_bet_reg.nii.gz flair_bet_reg.nii.gz
                    # now apply the transformation matrices to the original images (pre hd-bet)
                    flirt -in CT1_reorient.nii.gz -out ct1_reg.nii.gz -ref t1_bet.nii.gz -applyxfm -init ct1_to_t1.mat -interp spline &
                    flirt -in T2_reorient.nii.gz -out t2_reg.nii.gz -ref t1_bet.nii.gz -applyxfm -init t2_to_t1.mat -interp spline &
                    flirt -in FLAIR_reorient.nii.gz -out flair_reg.nii.gz -ref t1_bet.nii.gz -applyxfm -init flair_to_t1.mat -interp spline &
                    wait
                    # now apply t1 brain mask to all registered images
                    fslmaths ct1_reg.nii.gz -mas t1_bet_mask.nii.gz CT1_reorient_reg_bet.nii.gz & # t1_bet_mask.nii.gz was generated by hd-bet (see above)
                    fslmaths t2_reg.nii.gz -mas t1_bet_mask.nii.gz T2_reorient_reg_bet.nii.gz & # t1_bet_mask.nii.gz was generated by hd-bet (see above)
                    fslmaths flair_reg.nii.gz -mas t1_bet_mask.nii.gz FLAIR_reorient_reg_bet.nii.gz & # t1_bet_mask.nii.gz was generated by hd-bet (see above)
                    wait
                    # run hd-glio
                    hd_glio_predict -t1 T1_reorient.nii.gz -t1c CT1_reorient_reg_bet.nii.gz -t2 T2_reorient_reg_bet.nii.gz \
                     -flair FLAIR_reorient_reg_bet.nii.gz -o OUTPUT_FILE.nii.gz
                    # change back
                    cd $cwd
                else
                    echo "Running HD-GLIO-AUTO using docker"
                    docker run --gpus all --mount type=bind,source=${cwd}/$hdglioinputdir,target=/input \
                     --mount type=bind,source=${cwd}/$hdgliooutputdir,target=/output \
                    jenspetersen/hd-glio-auto
                fi
                
                mrcalc $hdgliooutputdir/segmentation.nii.gz 1 -ge $globalresultsdir/lesion.nii

            else
                echo "HD-GLIO-AUTO already done"
            fi
        else
            echo "Not possible to run HD-GLIO-AUTO"
        fi 
    fi
}

function KUL_run_VBG {
    if [ $vbg -eq 1 ]; then
        vbg_test="lesion_wf/output_LWF/sub-${participant}/sub-${participant}_aparc+aseg.nii.gz"
        if [[ ! -f $vbg_test ]]; then
            echo "Starting KUL_VBG"
            KUL_VBG.sh -p ${participant} -l $globalresultsdir/lesion.nii -z T1 -b -B 1 -t -F -n $ncpu -v
            mkdir -p freesurfer
            ln -s ${cwd}/lesion_wf/output_LWF/sub-${participant}/sub-${participant}_FS_output/sub-${participant}/ freesurfer
            echo "done" > freesurfer/sub-${participant}_freesurfer_is.done
        else
            echo "KUL_VBG has already run"
        fi
    fi
}

function KUL_run_msbp {
    if [ ! -f KUL_LOG/sub-${participant}_MSBP.done ]; then

        echo " starting MSBP"

        # there seems tpo be a problem with docker if the fsaverage dir is a soft link; so we delete the link and hardcopy it
        rm -fr $cwd/BIDS/derivatives/freesurfer/fsaverage
        cp -r $FREESURFER_HOME/subjects/fsaverage $cwd/BIDS/derivatives/freesurfer/fsaverage

        my_cmd="docker run --rm -u $(id -u) -v $cwd/BIDS:/bids_dir \
         -v $cwd/BIDS/derivatives:/output_dir \
         -v $HOME/KUL_apps/freesurfer/license.txt:/opt/freesurfer/license.txt \
         sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
         --participant_label $participant --isotropic_resolution 1.0 --thalamic_nuclei \
         --brainstem_structures --skip_bids_validator --fs_number_of_cores $ncpu \
         --multiproc_number_of_cores $ncpu $str_silent"
        #echo $my_cmd
        eval $my_cmd
        
        touch KUL_LOG/sub-${participant}_MSBP.done
        
    else
        echo "MSBP already done"
    fi
}

function KUL_run_TCKSEG {

    echo " starting FWT VOI generation"
    my_cmd="KUL_FWT_make_VOIs.sh -p ${participant} \
     -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
     -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
     -c $cwd/study_config/trial_tracks_list_2.txt \
     -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
     -n $ncpu $str_silent"
    eval $my_cmd

    echo " starting FWT tracking"
    my_cmd="KUL_FWT_make_TCKs.sh -p ${participant} \
     -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
     -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
     -c $cwd/study_config/trial_tracks_list_2.txt \
     -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
     -T 1 -a iFOD2 \
     -Q -S \
     -n $ncpu $str_silent"
    eval $my_cmd
}

function KUL_compute_melodic {
# run FSL Melodic
computedir="$cwd/compute/FSL/melodic/sub-$participant"
fmridatadir="$computedir/fmridata"
scriptsdir="$computedir/scripts"
fmriprepdir="fmriprep/sub-$participant/func"
globalresultsdir=$cwd/RESULTS/sub-$participant
searchtask="_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii"

mkdir -p $fmridatadir
mkdir -p $computedir/RESULTS
mkdir -p $globalresultsdir

if [ ! -f KUL_LOG/sub-${participant}_melodic.done ]; then
    echo "Preparing for Melodic"
    tasks=( $(find $fmriprepdir -name "*${searchtask}.gz" -type f) )
    # we loop over the found tasks
    for task in ${tasks[@]}; do
        d1=${task#*_task-}
        shorttask=${d1%_space*}
        #echo "$task -- $shorttask"
        echo " Analysing task $shorttask"
        fmrifile="${shorttask}${searchtask}"
        cp $fmriprepdir/*$fmrifile.gz $fmridatadir
        gunzip $fmridatadir/*$fmrifile.gz
        fmriresults="$computedir/stats_$shorttask"
        mkdir -p $fmriresults
        melodic_in="$fmridatadir/sub-${participant}_task-$fmrifile"
        # find the TR
        tr=$(mrinfo $melodic_in -spacing | cut -d " " -f 4)
        # make model and contrast
        dyn=$(mrinfo $melodic_in -size | cut -d " " -f 4)
        t_glm_con="$kul_main_dir/share/FSL/fsl_glm.con"
        t_glm_mat="$kul_main_dir/share/FSL/fsl_glm_${dyn}dyn.mat"        
        # set dimensionality and model for rs-/a-fMRI
        if [[ $shorttask == *"rest"* ]]; then
            dim="--dim=15"
            model=""
        else
            dim=""
            model="--Tdes=$t_glm_mat --Tcon=$t_glm_con"
        fi
        
        #melodic -i Melodic/sub-Croes/fmridata/sub-Croes_task-LIP_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii -o test/ --report --Tdes=glm.mat --Tcon=glm.con
        melodic -i $melodic_in -o $fmriresults --report --tr=$tr --Oall $model $dim
        
        # now we compare to known networks
        mkdir -p $fmriresults/kul
        fslcc --noabs -p 3 -t .204 $kul_main_dir/atlasses/Yeo2011_rsfMRI_in_FSL_Space/yeo2011_7_liberal_combined.nii.gz \
         $fmriresults/melodic_IC.nii.gz > $fmriresults/kul/kul_networks.txt

        while IFS=$' ' read network ic stat; do
            echo $network
            network_name=$(sed "${network}q;d" $kul_main_dir/atlasses/Yeo2011_rsfMRI_in_FSL_Space/yeo2011_7_liberal_combined_networks.txt)
            echo $network_name
            icfile="$fmriresults/stats/thresh_zstat${ic}.nii.gz"
            network_file="$fmriresults/kul/melodic_${network_name}_ic${ic}.nii.gz"
            echo $icfile
            echo $network_file
            mrcalc $icfile 2 -gt $icfile -mul $network_file

            # since Melodic analysis was in MNI space, we transform back in native space
            input=$network_file
            output=$globalresultsdir/melodic_${shorttask}_${network_name}_ic${ic}.nii
            transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5
            find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
            reference=${find_T1w[0]}
            echo "input=$input"
            echo "output=$output"
            echo "transform=$transform"
            echo "reference=$reference"
            KUL_antsApply_Transform
        done < $fmriresults/kul/kul_networks.txt
    done
    touch KUL_LOG/sub-${participant}_melodic.done
else
    echo "Melodic analysis already done"
fi
}

# --- MAIN ---
globalresultsdir=$cwd/RESULTS/sub-$participant

# STEP 1 - BIDS conversion
KUL_convert2bids

# Run BIDS validation
if [ ! -f KUL_LOG/sub-${participant}_1_bidscheck.done ]; then 
    docker run -ti --rm -v ${cwd}/BIDS:/data:ro bids/validator /data

    read -p "Are you happy? (y/n) " answ
    if [[ ! "$answ" == "y" ]]; then
        exit 1
    else
        touch KUL_LOG/sub-${participant}_1_bidscheck.done
    fi
fi 

# STEP 2 - run fmriprep/dwiprep and continue
KUL_run_fmriprep &
KUL_run_dwiprep &

# STEP 3 - run HD-GLIO-AUTO
if [ $hdglio -eq 1 ];then
    KUL_segment_tumor
fi

# STEP 4 - run VBG+freesurfer or freesurfer only
if [ $vbg -eq 1 ];then
    KUL_run_VBG &
else
    KUL_run_freesurfer &
fi

# WAIT FOR ALL TO FINISH
wait

# STEP 5 - run SPM/melodic/msbp
KUL_dwiprep_anat.sh -p $participant -n $ncpu > /dev/null &
KUL_compute_SPM &
KUL_compute_melodic &
KUL_run_msbp &

wait 

KUL_run_TCKSEG

echo "Finished"

