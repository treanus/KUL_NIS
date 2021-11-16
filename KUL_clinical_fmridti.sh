#!/bin/bash -e
# Bash shell script to analyse clinical fMRI/DTI
#
# Requires matlab fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 11/11/2021
version="0.7"

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

Optional arguments:

     -t:  processing type
        type 1: (DEFAULT) do hd-glio and vbg (tumor with T1w, cT1w, T2w and FLAIR)
        type 2: do vbg with manual mask (tumor but missing one of T1w, cT1w, T2w and FLAIR; 
                    put lesion.nii in RESULTS/sub-{participant}/Anat)
        type 3: don't run hd-glio nor vbg (cavernoma, epilepsy, etc... cT1w)
     -d:  dicom zip file (or directory)
     -c:  make a backup and cleanup 
     -n:  number of cpu to use (default 15)
     -r:  redo certain steps (program will ask)
     -R:  make results ready
     -v:  show output from commands

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ants_verbose=1
ncpu=15
bc=0 
type=1
redo=0
results=0 

# Set required options
p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:d:n:Rrcv" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        t) #type
			type=$OPTARG
		;;
        d) #dicomzip
			dicomzip=$OPTARG
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
		c) #backup&clean
			bc=1
		;;
        r) #redo
			redo=1
		;;
        R) #make results
			results=1
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

# MRTRIX and others verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
    str_silent=" > /dev/null 2>&1" 
    ants_verbose=0
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


if [ $bc -eq 1 ]; then
    # clean some stuff
    clean_dwiprep="./dwiprep/sub-${participant}/sub-${participant}/kul_dwifsl* \
        ./dwiprep/sub-${participant}/sub-${participant}/raw \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi_orig* \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi_preproced.mif"

    rm -fr $clean_dwiprep

    #we backup everything
    bck_bids="./BIDS/sub-${participant}"
    bck_dicom="./DICOM/${participant}*"
    bck_derivatives_KUL_VBG="./BIDS/derivatives/KUL_VBG/sub-${participant}"
    bck_derivatives_KUL_compute="./BIDS/derivatives/KUL_compute/sub-${participant}"
    bck_derivatives_freesurfer="./BIDS/derivatives/freesurfer/sub-${participant}"
    bck_derivatives_cmp="./BIDS/derivatives/cmp/sub-${participant}"
    bck_derivatives_nipype="./BIDS/derivatives/nipype/sub-${participant}"
    bck_derivatives_ini="./BIDS/derivatives/sub-${participant}_anatomical_config.ini"
    bck_fmriprep="./fmriprep/sub-${participant}*"
    bck_dwiprep="./dwiprep/sub-${participant}"
    bck_karawun="./Karawun/sub-${participant}"
    bck_results="./RESULTS/sub-${participant}"
    bck_conf="./study_config"

    tar --ignore-failed-read -cvzf sub-${participant}.tar.gz $bck_bids \
        $bck_dicom $bck_derivatives_freesurfer $bck_derivatives_KUL_compute $bck_derivatives_KUL_VBG \
        $bck_derivatives_cmp $bck_derivatives_nipype $bck_derivatives_ini $bck_fmriprep $bck_dwiprep \
        $bck_results $bck_karawun \
        $bck_conf

    read -p "Are you sure the backup is complete and continue with delete? (y/n) " answ
    if [[ ! "$answ" == "y" ]]; then
        exit 1
    else  
        rm -fr $bck_bids \
            $bck_dicom $bck_derivatives_freesurfer $bck_derivatives_KUL_compute $bck_derivatives_KUL_VBG \
            $bck_derivatives_cmp $bck_derivatives_nipype $bck_derivatives_ini $bck_fmriprep $bck_dwiprep \
            $bck_results $bck_karawun
    fi

    exit 0
fi

if [ $results -eq 1 ];then

    ### under development 
    results_final_output="RESULTS/sub-${participant}/${participant}4silvia/for_PACS"
    mkdir -p $results_final_output
    #SPM
    read -p "Which SPM results (.e.g. TAAL_run-2) " answ_spm
    read -p "Which SPM threshold (.e.g. 7.5) " answ_thr
    SPM_orig="RESULTS/sub-${participant}/SPM/SPM_${answ_spm}.nii"
    SPM_output="${results_final_output}/tbfMRI_${answ_spm}_thr${answ_thr}.nii"
    mrcalc $SPM_orig $answ_thr -gt $SPM_orig -mul $SPM_output
    exit

    #Melodic
    melodic_network="visual"
    melodic_thr=3
    melodic_orig="RESULTS/sub-${participant}/Melodic/melodic*${Melodic_network}*.nii"
    melodic_output="${results_final_output}/rsfMRI_visual_thr${melodic_thr}.nii"
    mrmath $melodic_orig mean - | \
        mrcalc - $melodic_thr -gt - -mul $melodic_output
    exit

fi

# --- functions ---
function KUL_check_redo {
    if [ $redo -eq 1 ];then
        read -p "Redo: Melodic? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_melodic.done"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_melodic.done >/dev/null 2>&1
            rm -f ${cwd}/RESULTS/sub-${participant}/Melodic/*
            rm -fr ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/FSL_melodic/*
        fi
    fi
}

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
        KUL_dcm2bids.sh -d $dicomzip -p ${participant} -c study_config/sequences.txt -e -v
    else
        echo "BIDS conversion already done"
    fi
}

function KUL_check_data {
    
    mkdir -p $globalresultsdir

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
    SWI=($(find $bidsdir -name "*run-01_SWI.nii.gz" -type f ))
    nSWI=${#SWI[@]}
    SWIp=($(find $bidsdir -name "*run-02_SWI.nii.gz" -type f ))
    nSWIp=${#SWIp[@]}
    echo "number of SWI magnitude: $nSWI"
    echo "number of SWI phase: $nSWIp"

    # check the T1w
    if [ $nT1w -lt 1 ]; then
        echo "No T1w (without Gd) found. Fmriprep will not run."
        echo " Is the BIDS dataset correct?"
        read -p "Are you sure you want to continue? (y/n)? " answ
        if [[ ! "$answ" == "n" ]]; then
            exit 1
        fi
    if 

    # check hd-glio-auto requirements
    if [ $hdglio -eq 1 ]; then
        if [ $nT1w -lt 1 ] || [ $ncT1w -lt 1 ] || [ $nT2w -lt 1 ] || [ $nT1w -lt 1 ] && ; then
            echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
            echo " At least one is missing. Is the BIDS dataset correct?"
            read -p "Are you sure you want to continue? (y/n)? " answ
            if [[ ! "$answ" == "n" ]]; then
                exit 1
            fi
        fi
    if 

    # check the BIDS
    find_fmri=($(find ${cwd}/BIDS/sub-${participant} -name "*_bold.nii.gz"))
    n_fMRI=${#find_fmri[@]}
    if [ $n_fMRI -eq 0 ]; then
        echo "WARNING: no fMRI data"
    fi

    find_dwi=($(find ${cwd}/BIDS/sub-${participant} -name "*_dwi.nii.gz"))
    n_dwi=${#find_dwi[@]}
    if [ $n_dwi -eq 0 ]; then
        echo "WARNING: no dwi data"
    fi
    echo -e "\n\n"

}

function KUL_rigid_register {
    warp_field="${registeroutputdir}/${source_mri_label}_reg2_T1w"
    output_mri="${globalresultsdir}/Anat/${source_mri_label}_reg2_T1w.nii.gz"
    echo "Rigidly registering $source_mri to $target_mri"
    my_cmd="antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$warp_field,$output_mri] \
    --interpolation BSpline \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$target_mri,$source_mri,1] \
    --transform Rigid[0.1] \
    --metric MI[$target_mri,$source_mri,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox $str_silent"
    eval $my_cmd
    echo "Done rigidly registering $source_mri to $target_mri"
}

function KUL_run_fmriprep {
    if [ ! -f fmriprep/sub-${participant}.html ]; then
        
        # preparing for fmriprep
        cp study_config/run_fmriprep.txt KUL_LOG/sub-${participant}_run_fmriprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        if [ $n_fMRI -gt 0 ]; then
            fmriprep_options="--fs-no-reconall --use-aroma --use-syn-sdc "
        else
            fmriprep_options="--fs-no-reconall --anat-only "
        fi
        sed -i.bck "s/fmriprep_options: /fmriprep_options: ${fmriprep_options}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        
        # running fmriprep
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_fmriprep.txt 
        
        # cleaning the working directory
        rm -fr fmriprep_work_${participant}
        
        # copying the result to the global results dir
        cp -f fmriprep/sub-$participant/anat/sub-${participant}_desc-preproc_T1w.nii.gz $globalresultsdir/Anat/T1w.nii.gz
        gunzip -f $globalresultsdir/Anat/T1w.nii.gz
        
        # create a GM mask in the global results dir
        mrcalc fmriprep/sub-$participant/anat/sub-${participant}_dseg.nii.gz 1 -eq \
            fmriprep/sub-$participant/anat/sub-${participant}_dseg.nii.gz -mul - | \
            maskfilter - median - | \
            maskfilter - dilate $globalresultsdir/Anat/T1w_GM.nii.gz

    else
        echo "Fmriprep already done"
    fi
}

function KUL_run_dwiprep {
    if [ ! -f dwiprep/sub-${participant}/dwiprep_is_done.log ]; then
        cp study_config/run_dwiprep.txt KUL_LOG/sub-${participant}_run_dwiprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_dwiprep.txt
        rm -f KUL_LOG/sub-${participant}_run_dwiprep.txt.bck
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_dwiprep.txt 
    else
        echo "Dwiprep already done"
    fi
}

function KUL_run_freesurfer {
    if [ ! -f BIDS/derivatives/freesurfer/${participant}_freesurfer_is.done ]; then
        cp study_config/run_freesurfer.txt KUL_LOG/sub-${participant}_run_freesurfer.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_freesurfer.txt
        rm -f KUL_LOG/sub-${participant}_run_freesurfer.txt.bck
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_freesurfer.txt 
    else
        echo "Freesurfer already done"
    fi
}


function task_exec {

    echo "  " | tee -a ${prep_log} 
    
    echo ${task_in} | tee -a ${prep_log} 

    echo " Started @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${prep_log} 

    eval ${task_in} | tee -a ${prep_log} 2>&1 &

    echo " pid = $! basicPID = $$ " | tee -a ${prep_log}

    wait ${pid}

    ### STEFAN NEED TO DO: is the sleep needed, or can it be shorter?
    sleep 5

    if [ $? -eq 0 ]; then
        echo Success | tee -a ${prep_log}
    else
        echo Fail | tee -a ${prep_log}

        exit 1
    fi

    echo " Finished @  $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${prep_log} 

    echo "  " | tee -a ${prep_log} 

    unset task_in

}


function KUL_run_fastsurfer {

    echo
    echo "Hybrid parcellation flag is set, now starting FastSurfer/FreeSurfer hybrid recon-all based part of VBG"
    echo

    # make your log file
    prep_log="KUL_LOG/sub-${participant}_run_fastsurfer.txt" 
    if [[ ! -f ${prep_log} ]] ; then
        touch ${prep_log}
    else
        echo "${prep_log} already created"
    fi

    fs_output="${cwd}/BIDS/derivatives/freesurfer"
    output_d="${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/FastSurfer"
    #str_op="${output_d}/${participant}"
    #fasu_output="${str_op}fastsurfer"
    fasu_output=$output_d
    #T1_4_parc="${str_op}_T1_nat_4parc.mgz"
    T1_4_parc="${cwd}/$T1w"
    echo $T1_4_parc

    recall_scripts="${fs_output}/sub-${participant}/scripts"
    echo $recall_scripts
 

    #search_wf_mark4=($(find ${recall_scripts} -type f 2> /dev/null | grep recon-all.done))
    #echo $search_wf_mark4


    #FS_brain="${fs_output}/${participant}/mri/brainmask.mgz"

    #new_brain="${str_pp}_T1_Brain_4FS.mgz"

    task_in="mkdir -p ${fs_output} >/dev/null 2>&1"
    task_exec

    # Run recon-all and convert the T1 to .mgz for display
    # running with -noskulltrip and using brain only inputs
    # for recon-all
    # if we can run up to skull strip, break, fix with hd-bet result then continue it would be much better
    # if we can switch to fast-surf, would be great also
    # another possiblity is using recon-all -skullstrip -clean-bm -gcut -subjid <subject name>
    
    echo "starting recon-all stage 1"
    task_in="recon-all -i ${T1_4_parc} -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -autorecon1 -no-isrunning"
    task_exec
    echo "done recon-all stage 1"

    #task_in="mri_convert -rl ${fs_output}/${participant}/mri/brainmask.mgz ${T1_BM_4_FS} ${clean_BM_mgz}"
    #task_exec

    #task_in="mri_mask ${FS_brain} ${T1_BM_4_FS} ${new_brain} && mv ${new_brain} ${fs_output}/${participant}/mri/brainmask.mgz && cp \
    #${fs_output}/${participant}/mri/brainmask.mgz ${fs_output}/${participant}/mri/brainmask.auto.mgz"
    #task_exec

    #exit


    FaSu_loc=$(which run_fastsurfer.sh)
    #nvd_cu=$(nvcc --version)
    user_id_str=$(id -u $(whoami))
    T1_4_FaSu=$(basename ${T1_4_parc})
    nvram=$(echo $(nvidia-smi --query-gpu=memory.free --format=csv) | rev | cut -d " " -f2 | rev)
    if [[ ! -z ${nvram} ]]; then
        if [[ ${nvram} -lt 6000 ]]; then
            batch_fasu="4"
        elif [[ ${nvram} -gt 6500 ]] && [[ ${nvram} -lt 7000 ]]; then
            batch_fasu="6"
        elif [[ ${nvram} -gt 7000 ]]; then
            batch_fasu="8"
        fi
    else
        batch_fasu="4"
    fi


    if [[ ! -z ${FaSu_loc} ]]; then

        if [ -z ${nvram} -lt 4000 ]; then

            FaSu_cpu=" --no_cuda "
            echo " Running FastSurfer without CUDA " | tee -a ${prep_log}

        else

            FaSu_cpu=""
            echo " Running FastSurfer with CUDA " | tee -a ${prep_log}

        fi

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        task_in="run_fastsurfer.sh --t1 ${T1_4_parc} \
        --sid sub-${participant} --sd ${fasu_output} --parallel --threads ${ncpu} \
        --fs_license $FS_LICENSE --py python ${FaSu_cpu} --ignore_fs_version --batch ${batch_fasu}"

        task_exec

    else

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        echo "Local FastSurfer not found, switching to Docker version" | tee -a ${prep_log}
        T1_4_FaSu=$(basename ${T1_4_parc})

        if [[ ! -z ${nvd_cu} ]]; then

            FaSu_v="gpu"

        else

            FaSu_v="cpu"

        fi

        task_in="docker run -v ${output_d}:/data -v ${fasu_output}:/output \
        -v $FREESURFER_HOME:/fs60 --rm --user ${user_id_str} fastsurfer:${FaSu_v} \
        --fs_license /fs60/$(basename $FS_LICENSE) --sid sub-${participant} \
        --sd /output/ --t1 /data/${T1_4_FaSu} \
        --parallel --threads ${ncpu}"

        task_exec

    fi


    # time to copy the surfaces and labels from FaSu to FS dir
    # here we run FastSurfer first and 

    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/surf ${output_d}/${participant}_FS_output/${participant}/
    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/label ${output_d}/${participant}_FS_output/${participant}/
    
    cp -rf ${output_d}/sub-${participant}/surf/* $fs_output/sub-${participant}/surf/
    cp -rf ${output_d}/sub-${participant}/label/* $fs_output/sub-${participant}/surf/label/ 


    # task_in="recon-all -s ${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -all -noskullstrip"

    # task_exec

    task_in="recon-all -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -noskullstrip -no-isrunning -make all"
    task_exec

    exit

    task_in="mri_convert -rl ${fs_output}/${participant}/mri/brain.mgz ${T1_brain_clean} ${fs_output}/${participant}/mri/real_T1.mgz"

    task_exec

    task_in="mri_convert -rl ${fs_output}/${participant}/mri/brain.mgz -rt nearest ${Lmask_o} ${fs_output}/${participant}/mri/Lmask_T1_bin.mgz"

    task_exec

    fs_parc_mgz="${fs_output}/${participant}/mri/aparc+aseg.mgz"



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
            
    result_global=$cwd/RESULTS/sub-$participant/SPM/SPM_${shorttask}.nii
            
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
    gm_result_global=$cwd/RESULTS/sub-$participant/SPM/SPM_${shorttask}_gm.nii
    mrgrid $gm_mask regrid -template $result_global $gm_mask2
    gm_mask3=$computedir/RESULTS/gm_mask_smooth_${shorttask}.nii.gz
    mrfilter $gm_mask2 smooth $gm_mask3
    #mrcalc $result_global $gm_mask3 0.3 -gt -mul $gm_result_global

} 

function KUL_compute_SPM {
    #  setup variables
    computedir="$kulderivativesdir/sub-$participant/SPM"
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
    #cp -f $fmriprepdir/../anat/sub-${participant}_desc-preproc_T1w.nii.gz $globalresultsdir/Anat/T1w.nii.gz
    #gunzip -f $globalresultsdir/Anat/T1w.nii.gz

    if [ ! -f KUL_LOG/sub-${participant}_SPM.done ]; then
        echo "Computing SPM"
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
                my_cmd="KUL_compute_SPM_matlab $silent"
                eval $my_cmd

                # do the combined analysis
                if [[ "$shorttask" == *"run-2" ]]; then
                    echo "this is run2, we run full analysis now" 
                    shorttask=${shorttask%_run-2}
                    #echo $shorttask
                    tcf="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs.m" #template config file
                    tjf="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs_job.m" #template job file
                    fmrifile="${shorttask}"
                    my_cmd="KUL_compute_SPM_matlab $silent"
                    eval $my_cmd
                fi
            
            fi
        done
        touch KUL_LOG/sub-${participant}_SPM.done
        echo "Done computing SPM"
    else
        echo "SPM analysis already done"
    fi
}

function KUL_segment_tumor {
    
    # this will segment the lesion automatically
    hdglioinputdir="$kulderivativesdir/sub-${participant}/hdglio/input"
    hdgliooutputdir="$kulderivativesdir/sub-${participant}/hdglio/output"
    if [ $hdglio -eq 1 ]; then

        if [ $nT1w -eq 1 ] && [ $ncT1w -eq 1 ] && [ $nFLAIR -eq 1 ] && [ $nT2w -eq 1 ];then
            mkdir -p $hdglioinputdir
            mkdir -p $hdgliooutputdir
            if [ ! -f "$hdgliooutputdir/volumes.txt" ]; then
                cp $T1w $hdglioinputdir/T1.nii.gz
                cp $cT1w $hdglioinputdir/CT1.nii.gz
                cp $FLAIR $hdglioinputdir/FLAIR.nii.gz
                cp $T2w $hdglioinputdir/T2.nii.gz
                
                echo "Running HD-GLIO-AUTO using docker"
                docker run --gpus all --mount type=bind,source=$hdglioinputdir,target=/input \
                    --mount type=bind,source=$hdgliooutputdir,target=/output \
                jenspetersen/hd-glio-auto
                
                #mrcalc $hdgliooutputdir/segmentation.nii.gz 1 -ge $globalresultsdir/Anat/lesion.nii -force
                maskfilter $hdgliooutputdir/segmentation.nii.gz fill $globalresultsdir/Anat/lesion.nii -force
                mrcalc $hdgliooutputdir/segmentation.nii.gz 1 -eq $globalresultsdir/Anat/lesion_perilesional_oedema.nii -force
                mrcalc $hdgliooutputdir/segmentation.nii.gz 2 -eq $globalresultsdir/Anat/lesion_solid_tumour.nii -force
                mrcalc $globalresultsdir/Anat/lesion.nii $globalresultsdir/Anat/lesion_perilesional_oedema.nii -sub \
                    $globalresultsdir/Anat/lesion_solid_tumour.nii -sub $globalresultsdir/Anat/lesion_central_necrosis.nii

                echo "Done running HD-GLIO-AUTO using docker"

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
        vbg_test="${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG/output_VBG/sub-${participant}/sub-${participant}_T1_nat_filled.nii.gz"
        if [[ ! -f $vbg_test ]]; then
            echo "Computing KUL_VBG"
            mkdir -p ${cwd}/BIDS/derivatives/freesurfer/sub-${participant}
            mkdir -p ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG

            # Need to update to dev version
            my_cmd="KUL_VBG.sh -S ${participant} \
                -l $globalresultsdir/Anat/lesion.nii \
                -o ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG \
                -m ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG \
                -z T1 -b -B 1 -t -P 3 -n $ncpu $silent"       
            eval $my_cmd

            # copy the output of VBG to the derivatives freesurfer directory
            cp -r ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG/output_VBG/sub-${participant}/sub-${participant}_FS_output/sub-${participant} \
                BIDS/derivatives/freesurfer/
            #rm -fr ${cwd}/BIDS/derivatives/KUL_compute//sub-${participant}/KUL_VBG/sub-${participant}/sub-${participant}_FS_output/sub-${participant}/${participant}
            #ln -s ${cwd}/lesion_wf/output_LWF/sub-${participant}/sub-${participant}_FS_output/sub-${participant}/ freesurfer
            echo "done" > BIDS/derivatives/freesurfer/sub-${participant}_freesurfer_is.done
            echo "Done computing KUL_VBG"
        else
            echo "KUL_VBG has already run"
        fi
    fi
}

function KUL_run_msbp {
    if [ ! -f KUL_LOG/sub-${participant}_MSBP.done ]; then

        echo "Running MSBP"

        # there seems tpo be a problem with docker if the fsaverage dir is a soft link; so we delete the link and hardcopy it
        rm -fr $cwd/BIDS/derivatives/freesurfer/fsaverage
        cp -r $FREESURFER_HOME/subjects/fsaverage $cwd/BIDS/derivatives/freesurfer/fsaverage

        my_cmd="docker run --rm -u $(id -u) -v $cwd/BIDS:/bids_dir \
         -v $cwd/BIDS/derivatives:/output_dir \
         -v $FS_LICENSE:/opt/freesurfer/license.txt \
         sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
         --participant_label $participant --isotropic_resolution 1.0 --thalamic_nuclei \
         --brainstem_structures --skip_bids_validator --fs_number_of_cores $ncpu \
         --multiproc_number_of_cores $ncpu $str_silent"
        #echo $my_cmd
        eval $my_cmd
        
        echo "Done MSBP"
        touch KUL_LOG/sub-${participant}_MSBP.done
        
    else
        echo "MSBP already done"
    fi
}

function KUL_run_FWT {
    config="tracks_list.txt"
    if [ ! -f KUL_LOG/sub-${participant}_FWT.done ]; then
        echo "Running FWT VOI generation"
        my_cmd="KUL_FWT_make_VOIs.sh -p ${participant} \
        -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
        -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
        -c $cwd/study_config/${config} \
        -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
        -o $kulderivativesdir/sub-${participant}/FWT \
        -n $ncpu $str_silent"
        eval $my_cmd

        echo "Running FWT tracking"
        my_cmd="KUL_FWT_make_TCKs.sh -p ${participant} \
        -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
        -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
        -c $cwd/study_config/${config} \
        -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
        -o $kulderivativesdir/sub-${participant}/FWT \
        -T 1 -a iFOD2 \
        -Q -S \
        -n $ncpu $str_silent"
        eval $my_cmd

        ln -s $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*fin_map_BT_iFOD2.nii.gz $globalresultsdir/Tracto/
        ln -s $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*fin_BT_iFOD2.tck $globalresultsdir/Tracto/
        pdfunite $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*_output/Screenshots/*fin_BT_iFOD2_inMNI_screenshot2_niGB.pdf $globalresultsdir/Tracto/Tracts_Summary.pdf
        touch KUL_LOG/sub-${participant}_FWT.done
        
    else
        echo "FWT already done"
    fi
}

function KUL_compute_melodic {
# run FSL Melodic
computedir="$kulderivativesdir/sub-$participant/FSL_melodic"
fmridatadir="$computedir/fmridata"
scriptsdir="$computedir/scripts"
fmriprepdir="fmriprep/sub-$participant/func"
globalresultsdir=$cwd/RESULTS/sub-$participant
searchtask="_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii"

mkdir -p $fmridatadir
mkdir -p $computedir/RESULTS
mkdir -p $globalresultsdir

if [ ! -f KUL_LOG/sub-${participant}_melodic.done ]; then
    echo "Computing Melodic"
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
        fslcc --noabs -p 3 -t .204 $kul_main_dir/atlasses/Local/Sunaert2021/KUL_NIT_networks.nii.gz \
            $fmriresults/melodic_IC.nii.gz > $fmriresults/kul/kul_networks.txt

        while IFS=$' ' read network ic stat; do
            #echo $network
            network_name=$(sed "${network}q;d" $kul_main_dir/atlasses/Local/Sunaert2021/KUL_NIT_networks.txt)
            #echo $network_name
            icfile="$fmriresults/stats/thresh_zstat${ic}.nii.gz"
            network_file="$fmriresults/kul/melodic_${network_name}_ic${ic}.nii.gz"
            #echo $icfile
            #echo $network_file
            mrcalc $icfile 2 -gt $icfile -mul $network_file -force

            # since Melodic analysis was in MNI space, we transform back in native space
            input=$network_file
            output=$globalresultsdir/Melodic/melodic_${shorttask}_${network_name}_ic${ic}.nii
            transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5
            find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
            reference=${find_T1w[0]}
            #echo "input=$input"
            #echo "output=$output"
            #echo "transform=$transform"
            #echo "reference=$reference"
            KUL_antsApply_Transform
        done < $fmriresults/kul/kul_networks.txt
    done
    echo "Done computing Melodic"
    touch KUL_LOG/sub-${participant}_melodic.done
else
    echo "Melodic analysis already done"
fi
}


function KUL_register_anatomical_images {

    if [ ! -f KUL_LOG/sub-${participant}_anat_reg.done ]; then 
        target_mri=$T1w
        registeroutputdir="$kulderivativesdir/sub-${participant}/antsregister"
        mkdir -p $registeroutputdir

        if [ $ncT1w -gt 0 ];then
            source_mri_label="cT1w"
            source_mri=$cT1w
            KUL_rigid_register
        fi
        if [ $nT2w -gt 0 ];then
            source_mri_label="T2w"
            source_mri=$T2w
            KUL_rigid_register
        fi
        if [ $nFLAIR -gt 0 ];then
            source_mri_label="FLAIR"
            source_mri=$FLAIR
            KUL_rigid_register
        fi
        if [ $nSWI -gt 0 ];then
            source_mri_label="SWI"
            source_mri=$SWI
            KUL_rigid_register

            input=$SWIp
            transform="${registeroutputdir}/${source_mri_label}_reg2_T1w0GenericAffine.mat"
            output="${globalresultsdir}/Anat/${source_mri_label}_phase_reg2_T1w.nii.gz"
            reference=$target_mri
            KUL_antsApply_Transform
        fi
        touch KUL_LOG/sub-${participant}_anat_reg.done
    else 
        echo "Anatomical registration already done"
    fi
}


# --- MAIN ---
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
mkdir -p $kulderivativesdir
globalresultsdir=$cwd/RESULTS/sub-$participant
mkdir -p $globalresultsdir/Anat
mkdir -p $globalresultsdir/SPM
mkdir -p $globalresultsdir/Melodic
mkdir -p $globalresultsdir/Tracto

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


# Check if fMRI and/or dwi data are present and/or to redo some processing
echo "Starting KUL_clinical_fmridti"
KUL_check_data
KUL_check_redo

# STEP 2 - run fmriprep/dwiprep and continue
KUL_run_fmriprep &
if [ $n_dwi -gt 0 ];then
    KUL_run_dwiprep &
fi

# don't run too many AI tools (hd-bet on dwi and HD-GLIO-AUTO on structural) simultaneously on a 6GB GPU - wait a bit...
if [ ! -f dwiprep/sub-${participant}/dwiprep_is_done.log ]; then
    sleep 600
fi

# STEP 3 - run HD-GLIO-AUTO
if [ $hdglio -eq 1 ];then
    KUL_segment_tumor
fi


# STEP 3B - regsiter all anatomical other data to the T1w without contrast
KUL_register_anatomical_images &


# STEP 4 - run VBG+freesurfer or freesurfer only
if [ $vbg -eq 1 ];then
    KUL_run_VBG &
else
    #KUL_run_freesurfer &
    fast=0
    if [ $fast -eq 1 ];then
        KUL_run_fastsurfer
    else
        KUL_run_freesurfer
    fi
fi

# WAIT FOR ALL TO FINISH
wait


# STEP 5 - run SPM/melodic/msbp
KUL_run_msbp &
KUL_dwiprep_anat.sh -p $participant -n $ncpu > /dev/null &
if [ $n_fMRI -gt 0 ];then
    KUL_compute_SPM &  
    KUL_compute_melodic &
fi

wait 

# STEP 6 - run Fun With Tracts
if [ $n_dwi -gt 0 ];then
    KUL_run_FWT
fi

echo "Finished"
