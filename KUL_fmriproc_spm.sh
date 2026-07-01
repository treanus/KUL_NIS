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
     -S:  smoothing FWHM in mm for SUSAN (default: adaptive = mean voxel size)
            e.g. -S 6 for 6mm FWHM, -S 8 for 8mm FWHM
     -P:  FWE-corrected p-value threshold for Bizzi thresholding (default: 0.01)
            e.g. -P 0.01, -P 0.005, -P 0.001
     -v:  verbose (0=silent, 1=normal, 2=verbose; default=1)


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
verbose_level=1
smooth_fwhm=0  # 0 = adaptive (mean voxel size)
pfwe=0.01


# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:s:S:P:v:" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            participant=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        S) #smoothing FWHM
            smooth_fwhm=$OPTARG
        ;;
        P) #FWE p-value threshold
            pfwe=$OPTARG
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

}

function KUL_threshold_SPM_Bizzi {

    if [ ! -f "$fmriresults/SPM.mat" ]; then
        kul_echo " Bizzi thresholding skipped: SPM.mat not found in $fmriresults"
        return
    fi

    spm_bizzi_script="${scriptsdir}/thresh_Bizzi_${fmrifile}.m"
    spm_bizzi_script=${spm_bizzi_script/run-/run}

    cp "$kul_main_dir/share/spm12/spm12_threshold_Bizzi.m" "$spm_bizzi_script"
    sed -i.bck "s|###FMRIRESULTS###|$fmriresults|" "$spm_bizzi_script"
    sed -i.bck "s|###PFWE###|$pfwe|g" "$spm_bizzi_script"
    rm -f "${spm_bizzi_script}.bck"

    cmd="$matlab_exe -nodisplay -nosplash -nodesktop -r \"run('$spm_bizzi_script');exit;\" $str_silent_SPM"
    eval $cmd

    # build FWE tag matching the MATLAB output filename logic (e.g. 0.01->FWE01, 0.005->FWE005)
    pfwe_tag=$(echo "$pfwe" | sed 's/0\.//' | sed 's/0*$//')
    fwe_tag="FWE${pfwe_tag}_k50"

    transform="${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
    find_T1w=($(find ${cwd}/BIDS/sub-${participant}/anat/ -name "*_T1w.nii.gz" ! -name "*gadolinium*"))
    reference=${find_T1w[0]}

    for thresh_tag in "p001unc_k50" "${fwe_tag}"; do
        thresh_nii="$fmriresults/spmT_0001_${thresh_tag}.nii"
        if [ -f "$thresh_nii" ]; then
            mni_result="$computedir/RESULTS/MNI/${fmrifile}_${thresh_tag}_space-MNI152NLin2009cAsym.nii"
            cp "$thresh_nii" "$mni_result"
            input="$mni_result"
            output="${globalresultsdir}/afMRI_${fmrifile}_${thresh_tag}.nii"
            KUL_antsApply_Transform
        fi
    done
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

if [ ! -f KUL_LOG/sub-${participant}_SPM.done ]; then
    
    # find the output of fmriprep
    fmriprep_match=($(find $fmriprepdir/func -name "*${fmriprep_output_type}.gz" -type f))

    # find the unique tasks
    tasks=()
    for match in ${fmriprep_match[@]}; do
        match_tmp1=${match[@]#*_task-}
        match_tmp2=${match_tmp1[@]%_space*}
        match_tmp3=${match_tmp2[@]%_run*}
        tasks=(${tasks[@]} $match_tmp3)
    done
    uniqe_tasks=($(for i in ${tasks[@]}; do echo $i; done | sort -u))


    # we loop over the unique tasks
    for task in ${uniqe_tasks[@]}; do
  
        if [[ ! "$task" = *"rest"* ]]; then
            kul_echo " Analysing task $task"
            task_and_type_1="*${task}*${fmriprep_output_type}"

            # find the number of runs
            runs_sharp=($(find $fmriprepdir/func -name "*${task_and_type_1}.gz" -type f))
            # SUSAN edge-preserving smoothing: adaptive sigma (mean voxel size), bt = 2/3 masked median
            for run_sharp in ${runs_sharp[@]}; do
                smooth_out="$(dirname ${run_sharp})/$(basename ${run_sharp} .nii.gz)_smooth_6mm.nii.gz"
                if [[ ! -f "${smooth_out}" ]]; then
                    spacing=($(mrinfo ${run_sharp} -spacing))
                    if (( $(echo "$smooth_fwhm > 0" | bc -l) )); then
                        sigma=$(echo "$smooth_fwhm / 2.3548" | bc -l)
                    else
                        sigma=$(echo "(${spacing[0]} + ${spacing[1]} + ${spacing[2]}) / 3" | bc -l)
                    fi
                    # brightness threshold from masked median (2/3 rule)
                    mask_file="$(dirname ${run_sharp})/$(basename ${run_sharp} _desc-preproc_bold.nii.gz)_desc-brain_mask.nii.gz"
                    if [[ -f "${mask_file}" ]]; then
                        p50=$(fslstats ${run_sharp} -k ${mask_file} -p 50)
                    else
                        p50=$(fslstats ${run_sharp} -p 50)
                    fi
                    bt=$(echo "$p50 * 0.66666" | bc -l)
                    # boldref as USAN edge reference
                    boldref="$(dirname ${run_sharp})/$(basename ${run_sharp} _desc-preproc_bold.nii.gz)_desc-coreg_boldref.nii.gz"
                    if [[ -f "${boldref}" ]]; then
                        susan ${run_sharp} ${bt} ${sigma} 3 1 1 ${boldref} ${bt} ${smooth_out}
                    else
                        susan ${run_sharp} ${bt} ${sigma} 3 1 0 ${smooth_out}
                    fi
                fi
            done

            # edited by AR 04/11/2022
            fmriprep_output_type_2=$(echo ${fmriprep_output_type} | cut -d "." -f1)
            task_and_type_2="*${task}*${fmriprep_output_type_2}_smooth_6mm.nii"

            # edited by AR 04/11/2022
            runs=($(find $fmriprepdir/func -name "*${task_and_type_2}.gz" -type f))
            
            # unzip each run
            for run in ${runs[@]}; do
                cp $run $fmridatadir
                shortrun=$(basename $run)
                kul_echo " gunzipping $shortrun"
                gunzip -f $fmridatadir/$shortrun
            done

            # determine the TR
            TR=$(mrinfo $fmridatadir/${shortrun%.gz} -spacing | awk '{print $(NF)}')
            kul_echo " the repetition time (TR) of $shortrun is: $TR"

            n_runs=${#runs[@]}

            i_run=1
            for run in ${runs[@]}; do

                spm_template_config_file="$kul_main_dir/share/spm12/spm12_fmri_stats_1run.m" #template config file
                spm_template_job_file="$kul_main_dir/share/spm12/spm12_fmri_stats_1run_job.m" #template job file
                if [ $n_runs -gt 1 ]; then
                    run_id=$(basename "${run}" | grep -oP '(?<=_run-)\d+' | head -1)
                    fmrifile="${task}_run-${run_id}"
                elif [ $n_runs -eq 1 ]; then
                    fmrifile="${task}"
                fi

                echo " computing $fmrifile"
                KUL_compute_SPM_matlab
                KUL_threshold_SPM_Bizzi
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
                KUL_threshold_SPM_Bizzi
            fi

        fi
    done

    touch KUL_LOG/sub-${participant}_SPM.done
    echo "Done computing SPM"
else
    echo "SPM analysis already done"
fi
