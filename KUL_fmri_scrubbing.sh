#!/bin/bash 
# 
# Bash shell script to scrub resting-state functional MRI from fmriprep
#
# Requires Mrtrix3, csvkit
#
# @ Rob Colaes - KUL - rob.colaes@kuleuven.be
#
version="v0.1 - dd 20/05/2025"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd, mrtrix3new & $log_dir is made in main_functions


# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs fMRI scrubbing.

Usage:

  `basename $0` -p subject <OPT_ARGS>

Example:

  `basename $0` -p pat001 -n 6 

Required arguments:

     -p:  participant (anonymised name of the subject)
    OR 
     -a:  automatic mode (do all participants in BIDS folder)

Optional arguments:

     -s:  session (of the participant)
     -d:  n SD for dvars for the threshold: mean_dvars + n * SD_dvars
     -f:  framewise displacement threshold
     -l:  minimal number of volumes in remaining timeseries (5 minutes recommended, calculate this using your TR)
     -c:  space of the fmri (eg. MNI152NLin2009cAsym)
     -r:  denoised with https://github.com/arielletambini/denoiser (0/1)? (highly recommended to denoise before scrubbing, assumes output is in fmriprep: *_NR.nii.gz)
     -n:  number of cpu for parallelisation
     -o:  output directory
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
auto=0
ncpu=6
silent=1
verbose_level=1

dvars_threshold=3
fd_threshold=0.7
min_timepoints=334
denoised=0
space_flag=0

# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:n:s:d:f:l:r:c:o:v" OPT; do

        case $OPT in
        p) #participant
            p_flag=1
            participant=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        n) #parallel
            ncpu=$OPTARG
        ;;
        f) #dvars_threshold
            dvars_threshold=$OPTARG
        ;;
        f) #fd_threshold
            fd_threshold=$OPTARG
        ;;
        l) #min_timepoints
            min_timepoints=$OPTARG
        ;;
        o) #output_dir
            output_dir=$OPTARG
        ;;
        c) #space
            space_flag=1
            space=$OPTARG
        ;;
        r) #denoised
            denoised=$OPTARG
        ;;
        v) #verbose
            verbose_level=$OPTARG
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

# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Some parallelisation
FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

d=$(date "+%Y-%m-%d_%H-%M-%S")
#log=log/log_${d}.txt
log=${log_dir}/${script}_${d}.log

# --- MAIN ----------------

echo "subject fd_mean fd_vol_kept dvars_mean dvars_threshold dvars_vol_kept total_vol_kept" > $output_dir/scrubbing_information.txt

if [ $auto -eq 0 ]; then
    SUBJECTS=$participant
else
    exit
fi

for sub in ${SUBJECTS[@]}; do

    echo $sub

    if [ $s_flag -eq 1 ]; then
        timeseries=fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-rest_desc-confounds_timeseries.tsv
    else
        timeseries=fmriprep/sub-${sub}/func/sub-${sub}_task-rest_desc-confounds_timeseries.tsv
    fi

    if [ ! -f "$timeseries" ]; then
        continue  # Skip the rest of the loop if the file doesn't exist
    fi

    #### Framewise displacement: get indices
    fd=($(csvcut -t -c "framewise_displacement" $timeseries | tail -n +2))

    # Empty array that we add with volumes to keep:
    fd_keep=()

    for i in "${!fd[@]}"; do

        fd_decimal=$(echo "${fd[$i]}" | awk '{ printf "%f\n", $1 }')

        if (( $(echo "$fd_decimal < $fd_threshold" | bc -l) )); then

            # CAUTION: the way we coded here makes that volume 1 has index 0 -> which is the way it should be for mrconvert as it counts from 0.
            # Make sure to check and adapt this if necessary
            fd_keep+=($i)
            
        fi
    done

    num_fd_keep=${#fd_keep[@]}
    echo "Volumes kept by fd criterium = $num_fd_keep"

    # Let's also calculate some mean FD for plotting outliers
    n=0
    sum=0
    for num in "${fd[@]:1}"; do
        sum=$(echo "$sum + $num" | bc -l)
        ((n++))
    done
    mean_fd=$(echo "$sum / $n" | bc -l)
    echo $mean_fd

    #### Dvars: get indices not to scrub (threshold = mean(dvars)+3SD(dvars)) :
    dvars=($(csvcut -t -c "dvars" $timeseries | tail -n +2))

    # Mean DVARS
    n=0
    sum=0
    for num in "${dvars[@]:1}"; do
        sum=$(echo "$sum + $num" | bc -l)
        ((n++))
    done

    mean_dvars=$(echo "$sum / $n" | bc -l)
    echo $mean_dvars

    # SD DVARS
    sum_sq_diff=0
    for num in "${dvars[@]:1}"; do
        diff=$(echo "$num - $mean_dvars" | bc -l)
        sq_diff=$(echo "$diff * $diff" | bc -l)
        sum_sq_diff=$(echo "$sum_sq_diff + $sq_diff" | bc -l)
    done

    variance=$(echo "$sum_sq_diff / $n" | bc -l)
    std_dev=$(echo "scale=4; sqrt($variance)" | bc -l)

    # Threshold:
    threshold=$(echo "$mean_dvars + $dvars_threshold * $std_dev" | bc -l)

    # Empty array that we add with volumes to keep:
    dvars_keep=()

    for i in "${!dvars[@]}"; do

        dvars_decimal=$(echo "${dvars[$i]}" | awk '{ printf "%f\n", $1 }')

        if (( $(echo "$dvars_decimal < $threshold" | bc -l) )); then

            # CAUTION: the way we coded here makes that volume 1 has index 0 -> which is the way it should be for mrconvert as it counts from 0
            # Make sure to check and adapt this if necessary
            dvars_keep+=($i)
            
        fi
    done

    num_dvars_keep=${#dvars_keep[@]}
    echo "Volumes kept by standardized dvars criterium = $num_dvars_keep"

    #### CHECK TOTAL NUMBER OF VOLUMES -> needs to more than 5 minutes
    # To have at least 5 minutes we have to have minimal 334 volumes!
    mutual_array=($(comm -12 <(echo "${fd_keep[@]}" | tr ' ' '\n' | sort) <(echo "${dvars_keep[@]}" | tr ' ' '\n' | sort)))

    num_volumes_keep=${#mutual_array[@]}

    echo "$sub $mean_fd $num_fd_keep $mean_dvars $threshold $num_dvars_keep $num_volumes_keep" >> $output_dir/scrubbing_information.txt

    #### FINAL STEP: Keep only the volumes of our mutual_array with mrconvert

    # Convert mutual_array into comma seperated array, necessary for mrconvert
    final_array=$(IFS=,; echo "${mutual_array[*]}")

    # Define input based on options
    if [ $denoised -eq 1 ]; then
        end_dir=desc-preproc_bold_NR
    else
        end_dir=desc-preproc_bold
    fi

    if [ $space_flag -eq 1 ]; then
        space_dir=space-${space}_
    else
        space_dir=""
    fi

    if [ $s_flag -eq 1 ]; then
        input=fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-rest_${space_dir}${end_dir}.nii.gz
        echo $input
    else
        input=fmriprep/sub-${sub}/func/sub-${sub}_task-rest_${space_dir}${end_dir}.nii.gz
        echo $input
    fi

    # Define output
    if [ $s_flag -eq 1 ]; then
        output=$output_dir/sub-${sub}_${ses}_BOLD_final.nii.gz
    else
        output=$output_dir/sub-${sub}_BOLD_final.nii.gz
    fi

    # Only keep 'good' volumes
    if (( $(echo "$num_volumes_keep > $min_timepoints" | bc -l) )); then

       task_in="mrconvert $input $output -coord 3 ${final_array[@]} -force "
       KUL_task_exec $verbose_level "Scrubbing part 2: extracting volumes" "2_scrubbing"

    fi

done
