#!/bin/bash
# Bash shell script to analyse clinical fMRI/DTI
#
# Requires matlab fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/12/2021
version="0.8"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

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
     -B:  make a backup and cleanup 
     -n:  number of cpu to use (default 15)
     -r:  redo certain steps (program will ask)
     -R:  make results ready
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)

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
verbose_level=1

# Set required options
p_flag=0
d_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:d:n:v:RBr" OPT; do

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
            d_flag=1
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
		B) #backup&clean
			bc=1
		;;
        r) #redo
			redo=1
		;;
        R) #make results
			results=1
		;;
        v) #verbose
            verbose_level=$OPTARG
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

KUL_LOG_DIR="KUL_LOG/${script}/sub-${participant}"
mkdir -p $KUL_LOG_DIR

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
	export MRTRIX_QUIET=1
    silent=1
    str_silent=" > /dev/null 2>&1" 
    ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
    silent=0
    str_silent="" 
    ants_verbose=1
fi

# Determine what to process depending on patient lesion type
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

# The BACKUP and clean option
if [ $bc -eq 1 ]; then
    # clean some stuff
    clean_dwiprep="./dwiprep/sub-${participant}/sub-${participant}/*dwifsl*tmp* \
        ./dwiprep/sub-${participant}/sub-${participant}/raw \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi \
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

# The make RESULTS option
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
        read -p "Redo: KUL_dwiprep? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_run_dwiprep.txt"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_run_dwiprep.txt >/dev/null 2>&1
            rm -rf ${cwd}/dwiprep/sub-${participant} >/dev/null 2>&1
        fi
        read -p "Redo: Melodic? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_melodic.done"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_melodic.done >/dev/null 2>&1
            rm -f ${cwd}/RESULTS/sub-${participant}/Melodic/*
            rm -fr ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/FSL_melodic/*
        fi
        read -p "Redo: KUL_VBG? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_VBG.log"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_VBG.log >/dev/null 2>&1
            rm -fr ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG/*
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
        if [ $d_flag -eq 1 ]; then
            KUL_dcm2bids.sh -d $dicomzip -p ${participant} -c study_config/sequences.txt -e -v
        else
            echo "Error: no dicom zip file given."
        fi
    else
        echo "BIDS conversion already done"
    fi
}

function KUL_check_data {
    
    mkdir -p $globalresultsdir
    echo -e "\n\nAn overview of the BIDS data:"
    bidsdir="BIDS/sub-$participant"
    T1w=($(find $bidsdir -name "*T1w.nii.gz" ! -name "*gadolinium*" -type f ))
    nT1w=${#T1w[@]}
    echo "  number of non-contrast T1w: $nT1w"
    cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
    ncT1w=${#cT1w[@]}
    echo "  number of contrast enhanced T1w: $ncT1w"
    FLAIR=($(find $bidsdir -name "*FLAIR.nii.gz" -type f ))
    nFLAIR=${#FLAIR[@]}
    echo "  number of FLAIR: $nFLAIR"
    FGATIR=($(find $bidsdir -name "*FGATIR.nii.gz" -type f ))
    nFGATIR=${#FGATIR[@]}
    echo "  number of FGATIR: $nFGATIR"
    T2w=($(find $bidsdir -name "*T2w.nii.gz" -type f ))
    nT2w=${#T2w[@]}
    echo "  number of T2w: $nT2w"
    SWI=($(find $bidsdir -name "*run-01_SWI.nii.gz" -type f ))
    nSWI=${#SWI[@]}
    SWIp=($(find $bidsdir -name "*run-02_SWI.nii.gz" -type f ))
    nSWIp=${#SWIp[@]}
    echo "  number of SWI magnitude: $nSWI"
    echo "  number of SWI phase: $nSWIp"

    # check the T1w
    if [ $nT1w -eq 0 ]; then
        echo "No T1w (without Gd) found. Fmriprep will not run."
        echo " Is the BIDS dataset correct?"
        read -p "Are you sure you want to continue? (y/n)? " answ
        if [[ ! "$answ" == "n" ]]; then
            exit 1
        fi
    fi 

    # check hd-glio-auto requirements
    if [ $hdglio -eq 1 ]; then
        if [ $nT1w -lt 1 ] || [ $ncT1w -lt 1 ] || [ $nT2w -lt 1 ] || [ $nT1w -lt 1 ]; then
            echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
            echo " At least one is missing. Is the BIDS dataset correct?"
            read -p "Are you sure you want to continue? (y/n)? " answ
            if [[ ! "$answ" == "n" ]]; then
                exit 1
            fi
        fi
    fi 

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
    #echo "Rigidly registering $source_mri to $target_mri"
    antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$warp_field,$output_mri] \
    --interpolation BSpline \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$target_mri,$source_mri,1] \
    --transform Rigid[0.1] \
    --metric MI[$target_mri,$source_mri,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    #echo "Done rigidly registering $source_mri to $target_mri"
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


function KUL_run_fastsurfer {
if [ ! -f KUL_LOG/sub-${participant}_FastSurfer.done ]; then
    #echo "Hybrid parcellation flag is set, now starting FastSurfer/FreeSurfer hybrid recon-all based part of VBG"

    # make your log file
    #prep_log="KUL_LOG/sub-${participant}_run_fastsurfer.txt" 
    #if [[ ! -f ${prep_log} ]] ; then
    #    touch ${prep_log}
    #else
    #    echo "${prep_log} already created"
    #fi
    kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"

    fs_output="${cwd}/BIDS/derivatives/freesurfer"
    #output_d="${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/FastSurfer"
    #str_op="${output_d}/${participant}"
    #fasu_output="${str_op}fastsurfer"
    fasu_output="${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/FastSurfer"
    #T1_4_parc="${str_op}_T1_nat_4parc.mgz"
    if [ $vbg -eq 1 ];then
        T1_4_parc=${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG/output_VBG/sub-${participant}/sub-Casier_T1_nat_4parc.mgz
    else
        T1_4_parc="${cwd}/$T1w"
    fi
    #echo $T1_4_parc
    #ls -l $T1_4_parc

    #recall_scripts="${fs_output}/sub-${participant}/scripts"
    #echo $recall_scripts
 
    mkdir -p ${fs_output} >/dev/null 2>&1
    mkdir -p ${fasu_output} >/dev/null 2>&1

    # Run recon-all and convert the T1 to .mgz for display
    # running with -noskulltrip and using brain only inputs
    # for recon-all
    # if we can run up to skull strip, break, fix with hd-bet result then continue it would be much better
    # if we can switch to fast-surf, would be great also
    # another possiblity is using recon-all -skullstrip -clean-bm -gcut -subjid <subject name>
    
    #echo "starting recon-all stage 1"
    task_in="recon-all -i ${T1_4_parc} -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -autorecon1 -no-isrunning"
    KUL_task_exec $verbose_level "FastSurfer part 1: recon-all stage 1" "$KUL_LOG_DIR/FastSurfer"
    #echo "done recon-all stage 1"



    #task_in="mri_convert -rl ${fs_output}/${participant}/mri/brainmask.mgz ${T1_BM_4_FS} ${clean_BM_mgz}"
    #task_exec

    #task_in="mri_mask ${FS_brain} ${T1_BM_4_FS} ${new_brain} && mv ${new_brain} ${fs_output}/${participant}/mri/brainmask.mgz && cp \
    #${fs_output}/${participant}/mri/brainmask.mgz ${fs_output}/${participant}/mri/brainmask.auto.mgz"
    #task_exec


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
        batch_fasu="2"
    fi


    if [[ ! -z ${FaSu_loc} ]]; then

        if [ ${nvram} -lt 5500 ]; then
            FaSu_cpu=" --no_cuda "
            FaSu_mode="cpu"
            #echo " Running FastSurfer without CUDA " | tee -a ${prep_log}
        else
            FaSu_cpu=""
            FaSu_mode="cuda-gpu"
            #echo " Running FastSurfer with CUDA " | tee -a ${prep_log}
        fi

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        task_in="run_fastsurfer.sh --t1 ${T1_4_parc} \
        --sid sub-${participant} --sd ${fasu_output} --fsaparc --parallel --threads ${ncpu} \
        --fs_license $FS_LICENSE --py python ${FaSu_cpu} --ignore_fs_version --batch ${batch_fasu}"
        kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"
        KUL_task_exec $verbose_level "FastSurfer part 2: Fastsurfer itself (script & $FaSu_mode mode)" "$KUL_LOG_DIR/FastSurfer"

    else

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        echo "Local FastSurfer not found, switching to Docker version" | tee -a ${prep_log}
        T1_4_FaSu=$(basename ${T1_4_parc})
        dir_4_FaSu=$(dirname ${T1_4_parc})

        if [ ${nvram} -lt 5500 ]; then
            FaSu_v="cpu"
        else
            FaSu_v="gpu"
        fi

        task_in="docker run -v ${dir_4_FaSu}:/data -v ${fasu_output}:/output \
        -v $FREESURFER_HOME:/fs60 --rm --user ${user_id_str} fastsurfer:${FaSu_v} \
        --fs_license /fs60/$(basename $FS_LICENSE) --sid sub-${participant} \
        --sd /output/ --t1 /data/${T1_4_FaSu} \
        --parallel --threads ${ncpu}"
        kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"
        KUL_task_exec $verbose_level "FastSurfer part 2: Fastsurfer itself (docker & $FaSu_mode mode)" "$KUL_LOG_DIR/FastSurfer"

    fi

    #fs_output="${cwd}/BIDS/derivatives/freesurfer"
    #fasu_output="${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/FastSurfer"
    

    # time to copy the surfaces and labels from FaSu to FS dirtask_exec
    # here we run FastSurfer first and 

    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/surf ${output_d}/${participant}_FS_output/${participant}/
    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/label ${output_d}/${participant}_FS_output/${participant}/
    
    #cp -rf ${output_d}/sub-${participant}/surf/* $fs_output/sub-${participant}/surf/
    #cp -rf ${output_d}/sub-${participant}/label/* $fs_output/sub-${participant}/surf/label/ 

    rsync -azv ${fasu_output}/sub-${participant}/ ${fs_output}/sub-${participant}/

    #task_in="recon-all -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -all -noskullstrip"
    #task_exec

    task_in="recon-all -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} \
        -parallel -no-isrunning -make all"
    KUL_task_exec $verbose_level "FastSurfer part 3: recon-all -make-all" "FastSurfer"

    #exit

    #task_in="mri_convert -rl ${fs_output}/${participant}/mri/brain.mgz ${T1_brain_clean} ${fs_output}/${participant}/mri/real_T1.mgz"

    #task_exec

    #task_in="mri_convert -rl ${fs_output}/${participant}/mri/brain.mgz -rt nearest ${Lmask_o} ${fs_output}/${participant}/mri/Lmask_T1_bin.mgz"

    #task_exec

    #fs_parc_mgz="${fs_output}/${participant}/mri/aparc+aseg.mgz"
    touch KUL_LOG/sub-${participant}_FastSurfer.done

else
    echo "Already done FastSurfer"
fi
}


function KUL_segment_tumor {
    
    # Segmentation of the tumor using HD-GLIO-AUTO
    
    # check if it needs to be performed
    if [ $hdglio -eq 1 ]; then

        hdglioinputdir="$kulderivativesdir/sub-${participant}/hdglio/input"
        hdgliooutputdir="$kulderivativesdir/sub-${participant}/hdglio/output"
        
        # only run if not yet done
        if [ ! -f "$globalresultsdir/Anat/lesion_central_necrosis_or_cyst.nii" ]; then

            # prepare the inputs
            mkdir -p $hdglioinputdir
            mkdir -p $hdgliooutputdir
            cp $T1w $hdglioinputdir/T1.nii.gz
            cp $cT1w $hdglioinputdir/CT1.nii.gz
            cp $FLAIR $hdglioinputdir/FLAIR.nii.gz
            cp $T2w $hdglioinputdir/T2.nii.gz
            
            # run HD-GLIO-AUTO using docker
            if [ ! -f /usr/local/KUL_apps/HD-GLIO-AUTO/scripts/run.py ]; then
                task_in="docker run --gpus all --mount type=bind,source=$hdglioinputdir,target=/input \
                    --mount type=bind,source=$hdgliooutputdir,target=/output \
                    jenspetersen/hd-glio-auto"
                hdglio_type="docker"
            else
                task_in="python /usr/local/KUL_apps/HD-GLIO-AUTO/scripts/run.py -i $hdglioinputdir -o $hdgliooutputdir"
                hdglio_type="local install"
            fi
            KUL_task_exec $verbose_level "HD-GLIO-AUTO using $hdglio_type" "2_hdglioauto"

            # compute some additional output
            task_in="maskfilter $hdgliooutputdir/segmentation.nii.gz dilate $hdgliooutputdir/lesion_dil5.nii.gz -npass 5 -force; \
                maskfilter $hdgliooutputdir/lesion_dil5.nii.gz fill $hdgliooutputdir/lesion_dil5_fill.nii.gz -force; \
                maskfilter $hdgliooutputdir/lesion_dil5_fill.nii.gz erode $globalresultsdir/Anat/lesion.nii -npass 5 -force; \
                mrcalc $hdgliooutputdir/segmentation.nii.gz 1 -eq $globalresultsdir/Anat/lesion_perilesional_oedema.nii -force; \
                mrcalc $hdgliooutputdir/segmentation.nii.gz 2 -eq $globalresultsdir/Anat/lesion_solid_tumour.nii -force" 
                #mrcalc $globalresultsdir/Anat/lesion.nii $globalresultsdir/Anat/lesion_perilesional_oedema.nii -sub \
                #    $globalresultsdir/Anat/lesion_solid_tumour.nii -sub $globalresultsdir/Anat/lesion_central_necrosis_or_cyst.nii -force"
            KUL_task_exec $verbose_level "compute lesion, oedema & solid parts" "2_hdglioauto"
            
        else
            echo "HD-GLIO-AUTO already done"
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

            task_in="KUL_VBG.sh -S ${participant} \
                -l $globalresultsdir/Anat/lesion.nii \
                -o ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG \
                -m ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG \
                -z T1 -b -B 1 -t -P 3 -n $ncpu"
            KUL_task_exec $verbose_level "KUL_VBG" "7_VBG"
            #wait

            # Need to update to dev version
            #my_cmd="KUL_VBG.sh -S ${participant} \
            #    -l $globalresultsdir/Anat/lesion.nii \
            #    -o ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG \
            #    -m ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG \
            #    -z T1 -b -B 1 -t -P 3 -n $ncpu $str_silent_VBG"       
            #eval $my_cmd

            # copy the output of VBG to the derivatives freesurfer directory
            #cp -r ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG/output_VBG/sub-${participant}/sub-${participant}_FS_output/sub-${participant} \
            #    BIDS/derivatives/freesurfer/
            #echo "Done computing KUL_VBG"
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

        task_in="docker run --rm -u $(id -u) -v $cwd/BIDS:/bids_dir \
         -v $cwd/BIDS/derivatives:/output_dir \
         -v $FS_LICENSE:/opt/freesurfer/license.txt \
         sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
         --participant_label $participant --isotropic_resolution 1.0 --thalamic_nuclei \
         --brainstem_structures --skip_bids_validator --fs_number_of_cores $ncpu \
         --multiproc_number_of_cores $ncpu"
        KUL_task_exec $verbose_level "MSBP" "7_msbp"

        echo "Done MSBP"
        touch KUL_LOG/sub-${participant}_MSBP.done
        
    else
        echo "MSBP already done"
    fi
}

function KUL_run_FWT {
    config="tracks_list.txt"
    if [ ! -f KUL_LOG/sub-${participant}_FWT.done ]; then

        task_in="KUL_FWT_make_VOIs.sh -p ${participant} \
        -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
        -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
        -c $cwd/study_config/${config} \
        -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
        -o $kulderivativesdir/sub-${participant}/FWT \
        -n $ncpu"
        KUL_task_exec $verbose_level "KUL_FWT voi generation" "FWTvoi"


        task_in="KUL_FWT_make_TCKs.sh -p ${participant} \
        -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
        -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
        -c $cwd/study_config/${config} \
        -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
        -o $kulderivativesdir/sub-${participant}/FWT \
        -T 1 -a iFOD2 \
        -f 1 \
        -Q -S \
        -n $ncpu"
        KUL_task_exec $verbose_level "KUL_FWT tract generation" "FWTtck"

        #ln -s $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*fin_map_BT_iFOD2.nii.gz $globalresultsdir/Tracto/
        #ln -s $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*fin_BT_iFOD2.tck $globalresultsdir/Tracto/
        mcp "$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*_fin_map_BT_iFOD2.nii.gz" \
            "$globalresultsdir/Tracto/Tract-csd_#2.nii.gz"
        mcp "$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*_fin_BT_iFOD2.tck" \
            "$globalresultsdir/Tracto/Tract-csd_#2.tck"
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

    if [ $silent -eq 1 ] ; then
        str_silent_melodic=" >> KUL_LOG/sub-${participant}_melodic.log"
    fi

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
        task_in="melodic -i $melodic_in -o $fmriresults --report --tr=$tr --Oall $model $dim"
        KUL_task_exec 0
        kul_log_file=KUL_LOG/sub-${participant}_melodic.log
        wait

        # now we compare to known networks
        mkdir -p $fmriresults/kul
        task_in="fslcc --noabs -p 3 -t .204 $kul_main_dir/atlasses/Local/Sunaert2021/KUL_NIT_networks.nii.gz \
            $fmriresults/melodic_IC.nii.gz > $fmriresults/kul/kul_networks.txt"
        KUL_task_exec 0
        kul_log_file=KUL_LOG/sub-${participant}_melodic.log
        wait

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
            output=$globalresultsdir/Melodic/rsfMRI_${shorttask}_${network_name}_ic${ic}.nii
            transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5
            find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
            reference=${find_T1w[0]}
            #echo "input=$input"
            #echo "output=$output"
            #echo "transform=$transform"
            #echo "reference=$reference"
            KUL_antsApply_Transform $str_silent_melodic
        done < $fmriresults/kul/kul_networks.txt
    done
    echo "Done computing Melodic"
    touch KUL_LOG/sub-${participant}_melodic.done
else
    echo "Melodic analysis already done"
fi
}


function KUL_register_anatomical_images {
    check="KUL_LOG/sub-${participant}_anat_reg.done"
    if [ ! -f $check ]; then

        target_mri=$T1w
        registeroutputdir="$kulderivativesdir/sub-${participant}/antsregister"
        mkdir -p $registeroutputdir

        if [ $ncT1w -gt 0 ];then
            source_mri_label="cT1w"
            source_mri=$cT1w
            task_in="KUL_rigid_register"
            KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "3_register_anat"
        fi
        if [ $nT2w -gt 0 ];then
            source_mri_label="T2w"
            source_mri=$T2w
            task_in="KUL_rigid_register"
            KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "3_register_anat"
        fi
        if [ $nFLAIR -gt 0 ];then
            source_mri_label="FLAIR"
            source_mri=$FLAIR
            task_in="KUL_rigid_register"
            KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "3_register_anat"
        fi
        if [ $nFGATIR -gt 0 ];then
            source_mri_label="FGATIR"
            source_mri=$FGATIR
            task_in="KUL_rigid_register"
            KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "3_register_anat"
        fi
        if [ $nSWI -gt 0 ];then
            source_mri_label="SWI"
            source_mri=$SWI
            task_in="KUL_rigid_register"
            KUL_task_exec $verbose_level "Rigidly registering the $source_mri_label to the T1w" "3_register_anat"

            input=$SWIp
            transform="${registeroutputdir}/${source_mri_label}_reg2_T1w0GenericAffine.mat"
            output="${globalresultsdir}/Anat/${source_mri_label}_phase_reg2_T1w.nii.gz"
            reference=$target_mri
            task_in="KUL_antsApply_Transform"
            KUL_task_exec $verbose_level "Applying the rigid registration of SWIm to SWIp too" "3_register_anat"
        fi
        touch $check
    else 
        echo "Anatomical registration already done"
    fi
}


# --- MAIN ---

# STEP 1 - BIDS conversion
KUL_convert2bids

KUL_check_participant

kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
mkdir -p $kulderivativesdir
globalresultsdir=$cwd/RESULTS/sub-$participant
mkdir -p $globalresultsdir/Anat
mkdir -p $globalresultsdir/SPM
mkdir -p $globalresultsdir/Melodic
mkdir -p $globalresultsdir/Tracto

if [ $KUL_DEBUG -gt 0 ]; then 
    echo "kulderivativesdir: $kulderivativesdir"
    echo "globalresultsdir: $globalresultsdir"
fi

# Run BIDS validation
check_in=${KUL_LOG_DIR}/1_bidscheck.done
if [ ! -f $check_in ]; then

    docker run -ti --rm -v ${cwd}/BIDS:/data:ro bids/validator /data

    read -p "Are you happy? (y/n) " answ
    if [[ ! "$answ" == "y" ]]; then
        exit 1
    else
        touch $check_in
    fi
fi


# Check if fMRI and/or dwi data are present and/or to redo some processing
echo "Starting KUL_clinical_fmridti"
KUL_check_data
KUL_check_redo

# STEP 1 - run HD-GLIO-AUTO
if [ $hdglio -eq 1 ];then
    KUL_segment_tumor
fi


# STEP 2 - run fmriprep and continue
KUL_run_fmriprep &


# STEP 3 - run dwiprep and continue
if [ $n_dwi -gt 0 ];then
    KUL_run_dwiprep &
fi


# STEP 4 - regsiter all anatomical other data to the T1w without contrast
KUL_register_anatomical_images &
wait


# STEP 5 & 6 - run SPM & melodic
if [ $n_fMRI -gt 0 ];then
    
    task_in="KUL_fmriproc_spm.sh -p $participant"
    KUL_task_exec $verbose_level "KUL_fmriproc_spm" "5_fmriproc_spm"

    KUL_compute_melodic &

fi
wait 


# STEP 7 - run VBG
if [ $vbg -eq 1 ];then
    KUL_run_VBG 
fi
wait


# STEP 8 - run SPM/melodic/msbp
#KUL_run_fastsurfer
KUL_run_freesurfer
wait 


# STEP 9 - run SPM/melodic/msbp
KUL_run_msbp
wait


# STEP 10 run dwiprep_anat
task_in="KUL_dwiprep_anat.sh -p $participant -n $ncpu"
KUL_task_exec $verbose_level "KUL_dwiprep_anat" "6_dwiprep_anat"



# STEP 11 - run Fun With Tracts
if [ $n_dwi -gt 0 ];then
    KUL_run_FWT
fi

echo "Finished"
