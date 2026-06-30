#!/bin/bash
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - UZ/KUL - ahmed.radwan@uzleuven.be
#
# v0.1 - dd 19/01/2019 - jurassic version
version="v0.2 - dd 05/12/2021"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions


# To Do
#  - make this thing work!
#  - This script should be able to do the following:
# Decide what kind of fMRI data you have in your BIDS dir
# carry out a single subject analysis of all this BOLD data automatically
# for rs-fMRI ICA, time series clustering (SLIC/NCUT), atlas based roi to roi analysis and connectome generation should be automated
# we will probably need to use some matlab and python functions and commands here and there
# For tb-fMRI we can do FEAT +/- all the previous
# Include an option for high resolution rendering in the three orthogonal planes and mosaic generation in TIFF/JPEG (slicer)
# For mTE-rsfMRI, after MEICA, run same as single echo rs-fMRI, check out tedana also

# assuming most of the preproc is taken care of by fmriprep


# -----------------------------------  MAIN  ---------------------------------------------


# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs an automated task based fMRI spm12 GLM analysis assuming:
    - a 30 seconds REST followed by 30 seconds TASK epochs
    - having run fmriprep with aroma

Note: requires matlab and spm12 installed

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001

Required arguments:

     -p:  participant


Optional arguments:
     
     -s:  session
     -v:  verbose (0=silent, 1=normal, 2=verbose; default=1)


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
verbose_level=1


# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:v:" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            participant=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
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

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
	export MRTRIX_QUIET=1
    ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
    ants_verbose=1
fi

# Functions --------------------------------------------------------
function KUL_antsApply_Transform {
    if [ $KUL_DEBUG -gt 0 ]; then
        echo "input=$input"
        echo "output=$output"
        echo "transform=$transform"
        echo "reference=$reference"
    fi
    antsApplyTransforms -d 3 --float 1 \
        --verbose $ants_verbose \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n Linear
}

function KUL_compute_SPM_matlab {
    
    # 
    fmriresults="$computedir/RESULTS/stats_$fmrifile"
    
    # clean a possible old result, 
    rm -rf $fmriresults
    mkdir -p $fmriresults

    # prepare the job and config files    
    spm_participant_config_file="${scriptsdir}/stats_${fmrifile}.m" #participant config file
    spm_participant_job_file="${scriptsdir}/stats_${fmrifile}_job.m" #participant job file
    # get rid of - in filename, since this breaks -r in matlab
    spm_participant_config_file=${spm_participant_config_file/run-/run}
    spm_participant_job_file=${spm_participant_job_file/run-/run}
    #echo "$spm_participant_config_file -- $spm_participant_job_file"
    cp $spm_template_config_file $spm_participant_config_file
    cp $spm_template_job_file $spm_participant_job_file
    sed -i.bck "s|###JOBFILE###|$spm_participant_job_file|" $spm_participant_config_file
    sed -i.bck "s|###FMRIDIR###|$fmridatadir|" $spm_participant_job_file
    sed -i.bck "s|###FMRIFILE###|$fmrifile|" $spm_participant_job_file
    sed -i.bck "s|###FMRIRESULTS###|$fmriresults|" $spm_participant_job_file
    sed -i.bck "s|###TR###|$TR|" $spm_participant_job_file
    rm -f "${spm_participant_config_file}.bck"
    rm -f "${spm_participant_job_file}.bck"


    # call matlab and execute
    cmd="$matlab_exe -nodisplay -nosplash -nodesktop -r \"run('$spm_participant_config_file');exit;\" $str_silent_SPM"
    #echo $cmd
    eval $cmd


    result=$computedir/RESULTS/MNI/${fmrifile}_space-MNI152NLin2009cAsym.nii
    cp $fmriresults/spmT_0001.nii $result
    
    global_result=${globalresultsdir}/afMRI_${fmrifile}.nii
            
    # since SPM analysis was in MNI space, we transform back in native space
    input=$result
    output=$global_result
    transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
    reference=${find_T1w[0]}
    KUL_antsApply_Transform

    # compute the gray matter mask
    #gm_mask="$fmriprepdir/anat/sub-${participant}_label-GM_probseg.nii.gz"
    #gm_mask2=$computedir/RESULTS/gm_mask_${fmrifile}.nii.gz
    #gm_mask3=$computedir/RESULTS/gm_mask_smooth_${fmrifile}.nii.gz
    #mrgrid $gm_mask regrid -template $global_result $gm_mask2 -force
    #mrfilter $gm_mask2 smooth $gm_mask3 -force

    # compute a gray matter masked SPM result
    #gm_result_global=${globalresultsdir}/SPM_${fmrifile}_gm.nii
    #mrcalc $global_result $gm_mask3 0.1 -gt -mul $gm_result_global -force

} 

# MAIN --------------------------------------------------------------
matlab_exe=$(which matlab)

if [ $KUL_DEBUG -gt 0 ]; then 
    echo "matlab lives at $matlab_exe"
fi

if [[ -z "$matlab_exe" ]]; then
    echo "Matlab is required but not found on path. Exitting"
    exit 1
fi


KUL_check_participant


#  setup variables
kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
computedir="$kulderivativesdir/sub-$participant/SPM"
fmridatadir="$computedir/fmridata"
scriptsdir="$computedir/scripts"
fmriprepdir="${cwd}/fmriprep/sub-$participant"
globalresultsdir=$cwd/RESULTS/sub-$participant/SPM


if [ $KUL_DEBUG -gt 0 ]; then 
    echo "kulderivativesdir: $kulderivativesdir"
    echo "fmridatadir: $fmridatadir"
    echo "fmriprepdir: $fmriprepdir"
    echo "globalresultsdir: $globalresultsdir"
fi

mkdir -p $fmridatadir
mkdir -p $scriptsdir
mkdir -p $computedir/RESULTS/MNI
mkdir -p $globalresultsdir

# fmriprep_output_type="_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii"

# edited by AR 04/11/2022
fmriprep_output_type="_space-MNI152NLin2009cAsym_desc-preproc_bold.nii"


if [ $verbose_level -lt 2 ] ; then
    str_silent_SPM=" >> KUL_LOG/$script/sub-${participant}_spm12.log"
fi
#echo $str_silent_SPM

# Provide the anatomy
#cp -f $fmriprepdir/../anat/sub-${participant}_desc-preproc_T1w.nii.gz $globalresultsdir/Anat/T1w.nii.gz
#gunzip -f $globalresultsdir/Anat/T1w.nii.gz

if [ ! -f KUL_LOG/sub-${participant}_SPM.done ]; then
    #echo "Computing SPM"
    
    # find the output of fmriprep
    fmriprep_match=($(find $fmriprepdir/func -name "*${fmriprep_output_type}.gz" -type f))

    # find the unique tasks
    tasks=()
    for match in ${fmriprep_match[@]}; do
        #echo ${match[@]}
        match_tmp1=${match[@]#*_task-}
        #echo ${match_tmp1[@]}
        match_tmp2=${match_tmp1[@]%_space*}
        #echo ${match_tmp2[@]}
        match_tmp3=${match_tmp2[@]%_run*}
        tasks=(${tasks[@]} $match_tmp3)
    done
    #echo ${tasks[@]}
    uniqe_tasks=($(for i in ${tasks[@]}; do echo $i; done | sort -u))
    #echo ${uniqe_tasks[@]}


    # we loop over the unique tasks
    for task in ${uniqe_tasks[@]}; do
  
        if [[ ! "$task" = *"rest"* ]]; then
            kul_echo " Analysing task $task"
            task_and_type_1="*${task}*${fmriprep_output_type}"
            #echo $task_and_type_1
            
            # find the number of runs
            runs_sharp=($(find $fmriprepdir/func -name "*${task_and_type_1}.gz" -type f))
            # edited by AR 04/11/2022
            # temporary smoothin solution - better to use susan
            for run_sharp in ${runs_sharp[@]}; do
                if [[ ! -f "$(dirname ${run_sharp})/$(basename ${run_sharp} .nii.gz)_smooth_3mm.nii.gz" ]]; then
                    fslmaths ${run_sharp} -s 3 $(dirname ${run_sharp})/$(basename ${run_sharp} .nii.gz)_smooth_6mm.nii.gz
                fi
            done

            # edited by AR 04/11/2022
            fmriprep_output_type_2=$(echo ${fmriprep_output_type} | cut -d "." -f1)
            task_and_type_2="*${task}*${fmriprep_output_type_2}_smooth_6mm.nii"
            #echo $task_and_type_2

            # edited by AR 04/11/2022
            runs=($(find $fmriprepdir/func -name "*${task_and_type_2}.gz" -type f))
            
            # unzip each run
            for run in ${runs[@]}; do
                #echo $run
                cp $run $fmridatadir
                shortrun=$(basename $run)
                kul_echo " gunzipping $shortrun"
                gunzip -f $fmridatadir/$shortrun
            done

            # determine the TR
            TR=$(mrinfo $fmridatadir/${shortrun%.gz} -spacing | awk '{print $(NF)}')
            kul_echo " the repetition time (TR) of $shortrun is: $TR"

            n_runs=${#runs[@]}
            #echo $n_runs     

            i_run=1
            for run in ${runs[@]}; do
            
                #echo $run
                spm_template_config_file="$kul_main_dir/share/spm12/spm12_fmri_stats_1run.m" #template config file
                spm_template_job_file="$kul_main_dir/share/spm12/spm12_fmri_stats_1run_job.m" #template job file
                if [ $n_runs -gt 1 ]; then
                    fmrifile="${task}_run-${i_run}"
                elif [ $n_runs -eq 1 ]; then
                    fmrifile="${task}"
                fi

                echo " computing $fmrifile"
                KUL_compute_SPM_matlab
                ((i_run++))
            
            done

            #  the template files in KNS for SPM analysis
            if [ $n_runs -gt 1 ]; then
                if [ $n_runs -eq 2 ]; then
                    spm_template_config_file="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs.m" #template config file
                    spm_template_job_file="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs_job.m" #template job file
                elif [ $n_runs -eq 3 ]; then
                    spm_template_config_file="$kul_main_dir/share/spm12/spm12_fmri_stats_3runs.m" #template config file
                    spm_template_job_file="$kul_main_dir/share/spm12/spm12_fmri_stats_3runs_job.m" #template job file
                else
                    "Error: Not yet defined more than 3 runs. Exitting"
                    exit 1
                fi

                fmrifile="${task}"
                echo " computing aggregate n=${n_runs} ${fmrifile}"
                KUL_compute_SPM_matlab
            fi

        fi
    done

    # cleanup
    #rm -rf $fmridatadir

    touch KUL_LOG/sub-${participant}_SPM.done
    echo "Done computing SPM"
else
    echo "SPM analysis already done"
fi
