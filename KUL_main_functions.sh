#!/bin/bash
# Bash shell script to:
#  - define global functions 
#  - define defaults
#  - execute startup
#
# @ stefan.sunaert@kuleuven.be - v0.2 - 18/11/2021
#   update 26/08/2024


# MAIN FUNCTION - Function task_exec ###################################################################################
#  - obligatory variable to set
#       task_in (a command string that needs to be evaluated)
#  - facultative command line options
#       kul_verbose_level (0=silent, 1=normal; 2=verbose; 1=default)
#       kul_short_name (what to display as process)
#       kul_log_file (the path of the log file)
#
#  it will return $total_errorcount in the calling script (a sum of all errors that happened)
# 
# Example 1
#  task_in="a_big_long_process"
#  and also given are: 0 
#  and also given is: "KUL VBG" 
#  and also given is: "KUL_LOG/my_big_long_process"
#   It will run "process a_big_long_process"
#   Do it silently (the 0)
#   Display "KUL VBG" as the running process
#   Log to the files KUL_LOG/my_big_long_process.log and KUL_LOG/my_big_long_process.error.log
#
# Example 2
#  task_in="another big process"
#   Run process "another big process"
#   It will give output to the terminal and to log files (because default is kul_verbose_level 1)
#   Take the first 20 characters of "task_in" and display that (kul_short_name is generated automatically)
#   Log to the files another_big_process.log and another_big_process.error.log (/ -> _, and spaces too)
function KUL_task_exec {

    # task_in needs to be defined
    if [[ -z "$task_in" ]]; then
        echo "task_in is not defined, exitting"
        exit 1
    fi

    # task_participant, if undefined, and only 1 task_in, set to current participant, otherwise exit
    if [[ -z "$task_participant" ]]; then
        if [ ${#task_in[@]} -eq 1 ]; then
            task_participant[0]=$participant
        else
            echo "task_participant is not defined, exitting"
            exit 1
        fi
    fi

    local pidsArray=() # pids to wait for, separated by semi-colon
    local procsArray=() # name of procs to wait for, separated by semi-colon     
    local log_ttime=0 # local time instance for comparison
    local seconds_begin=$SECONDS # Seconds since the beginning of the script
    local exec_time=0 # Seconds since the beginning of this function
    local pidCount # number of given pids
    local c # counter for pids/procsArray
    local errorcount=0 # Number of pids that finished with errors

    local local_task_in
    local local_n_tasks
    local local_task_in_name

    local local_main_logdir=${cwd}/KUL_LOG/${script}

    
    # get the input variables; a check below wil put default values if $1/$2/$3 are empty
    local kul_verbose_level="$1"
    local kul_process_name="$2"
    local kul_log_files="$3"


    if [ $KUL_DEBUG -eq 1 ]; then
        echo "kul_verbose_level: $kul_verbose_level"
        echo "kul_process_name: $kul_process_name"
        echo "kul_log_files: $kul_log_files"
        echo "script: $script"
        echo "task_in: ${task_in[@]}"
        echo "task_participant: ${task_participant[@]}"
    fi

    
    # test the input to the function and of unset, set a default
    if [[ -z "$kul_verbose_level" ]]; then
        kul_verbose_level=1
    fi

    local_n_tasks=0
    for local_task_in in "${task_in[@]}"; do

        if [ $KUL_DEBUG -eq 1 ]; then
            echo "local_task_in: $local_task_in"
            echo "local_n_tasks: $local_n_tasks"
        fi
        #remove the double (or more) spaces from task_in
        local_task_in=$(echo $local_task_in | tr -s ' ')
        

        if [[ -z "$kul_process_name" ]]; then
            task_in_name[$local_n_tasks]="$(echo ${local_task_in:0:20} [sub-${task_participant[local_n_tasks]}])"
        else
            task_in_name[$local_n_tasks]="$kul_process_name [sub-${task_participant[local_n_tasks]}]"
        fi

        local local_main_logdir_participant="$local_main_logdir/sub-${task_participant[local_n_tasks]}"

        if [ ! -d "$local_main_logdir_participant" ]; then
            mkdir "$local_main_logdir_participant"
        fi

        if [[ -z "$kul_log_files" ]]; then 
            task_in_name_nospaces_tmp="${task_in_name[$local_n_tasks]// /_}"
            task_in_name_nospaces="${task_in_name_nospaces_tmp////_}"
            kul_log_file[$local_n_tasks]="$local_main_logdir_participant/"${task_in_name_nospaces}".log"
            kul_errorlog_file[$local_n_tasks]="$local_main_logdir_participant/"${task_in_name_nospaces}".error.log"
            kul_log_command[$local_n_tasks]="$local_main_logdir_participant/"${task_in_name_nospaces}".command"
        else
            kul_log_file[$local_n_tasks]="$local_main_logdir_participant/"${kul_log_files}.log""
            kul_errorlog_file[$local_n_tasks]="$local_main_logdir_participant/"${kul_log_files}.error.log""
            kul_log_command[$local_n_tasks]="$local_main_logdir_participant/"${kul_log_files}.command""
        fi

        ### STEP 2 - execute the task_in
        if [ $kul_verbose_level -lt 2 ]; then 

            #echo "using >"
            tput dim
            local task_in_final="($local_task_in)  >>${kul_log_file[$local_n_tasks]} 2>>${kul_errorlog_file[$local_n_tasks]}"
            #echo $task_in_final
            eval ${task_in_final} &
            tput sgr0
        
        else
            
            #echo "using tee"
            tput dim
            local task_in_final="($local_task_in)  > >(tee -a ${kul_log_file[$local_n_tasks]}) 2> >(tee -a ${kul_errorlog_file[$local_n_tasks]})"
            eval ${task_in_final} &
            tput sgr0
        fi
        
        # set the pids, first we get the pid by "$!", then feed it into an array
        task_in_pid="$!"
        pidsArray+=($task_in_pid)
        procsArray+=("${task_in_name[$local_n_tasks]}")
        #echo "procsArray: ${procsArray[$local_n_tasks]}"

        ### STEP 3 - give some information
        if [ $kul_verbose_level -gt 0 ]; then
            tput bold
            echo "KUL_task_exec: ${task_in_name[$local_n_tasks]}... started @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${kul_log_file[$local_n_tasks]}
            tput sgr0
        fi
        if [ $kul_verbose_level -eq 2 ]; then 
            tput dim
            echo -e "   The task_in command: ${local_task_in}" | tee -a ${kul_log_file[$local_n_tasks]}
            tput sgr0
        fi
        # log the command itself
        echo ${local_task_in} > ${kul_log_command[$local_n_tasks]}

        ((local_n_tasks++))

    done

    ### STEP 4 - keep checking if the process is still running    
    pidCount=${#pidsArray[@]}
    #echo "  pidCount: $pidCount"
    #echo "  pidsArray: ${pidsArray[@]}"
    #echo "  procsArray: ${procsArray[@]}"

    while [ ${#pidsArray[@]} -gt 0 ]; do

        newPidsArray=()
        newProcsArray=()
        c=0

        total_exec_time=$(($SECONDS - $script_start_time))
        #echo $total_exec_time
        total_exec_time_min=$(echo "scale=2; $total_exec_time/60" | bc)

        ## Log a standby message every hour
        every_time=1201
        exec_time=$(($SECONDS - $seconds_begin))
        if [ $((($exec_time + 1) % $every_time)) -eq 0 ]; then
            if [ $log_ttime -ne $exec_time ]; then
                log_ttime=$exec_time
                log_min=$(echo "scale=1; $log_ttime/60" | bc)
                if [ $kul_verbose_level -gt 0 ]; then
                    tput dim
                    echo "  Current tasks [${procsArray[@]}] still running after $log_min minutes with pids [${pidsArray[@]}]."
                    echo "    Total script time: $total_exec_time_min minutes"
                    tput sgr0
                fi
            fi
        fi


        for pid in "${pidsArray[@]}"; do
            #echo "pid: $pid"
            #echo "proc: ${procsArray[c]}"
            if kill -0 $pid > /dev/null 2>&1; then
                newPidsArray+=($pid)
                #echo "newPidsArray: ${newPidsArray[@]}"
                newProcsArray+=("${procsArray[c]}")
                #echo "newProcsArray: ${newProcsArray[@]}"
            else
                wait $pid
                result=$?
                #echo "result: $result"
                if [ $result -ne 0 ]; then

                    errorcount=$((errorcount+1))
                    fail_exec_time_seconds=$(($SECONDS - $seconds_begin))
                    fail_exec_time_minutes=$(echo "scale=1; $final_exec_time_seconds/60" | bc)
                    tput bold; tput setaf 1
                    echo "  *** WARNING! **** Process ${procsArray[c]} with pid $pid might have failed after $fail_exec_time_minutes minutes. (with exitcode [$result]). Check the ${kul_errorlog_file[$c]} log-file" | tee -a ${kul_errorlog_file[$c]}
                    tput sgr0
                
                else
                    
                    final_exec_time_seconds=$(($SECONDS - $seconds_begin))
                    final_exec_time_minutes=$(echo "scale=2; $final_exec_time_seconds/60" | bc)
                    if [ $kul_verbose_level -gt 0 ]; then
                        tput setaf 2
                        #echo "c: $c"
                        #echo "procsArray: ${procsArray[c]}"
                        echo " ${procsArray[c]} finished successfully after $final_exec_time_minutes minutes" | tee -a ${kul_log_file[$c]}
                        echo "    Total script time: $total_exec_time_min minutes"
                        tput sgr0
                    fi
                fi
            fi
            c=$((c+1))
        done


        pidsArray=("${newPidsArray[@]}")
        procsArray=("${newProcsArray[@]}")
        sleep 1

    done

    ### STEP 5 - return the status of execution
    if [ $errorcount -eq 0 ]; then
        if [ $kul_verbose_level -eq 2 ]; then
            echo -e "Success" | tee -a ${kul_log_file}
        fi
    else
        echo -e "Fail" | tee -a ${kul_log_file}
    fi

    unset task_in
    unset task_participant

    if [[ ! -z $total_errorcount ]]; then
        total_errorcount=$(($total_errorcount + $errorcount))
        if [ $kul_verbose_level -eq 2 ]; then
	        echo -e "total_errorcount: $total_errorcount\n\n\n"
        fi
    fi

    return $(( errorcount > 0 ? 1 : 0 ))

}


# Resolve the first file matching a glob pattern; exit with an error if none found
function find_first_match {
    local glob="$1"
    local description="$2"
    local matches=()

    shopt -s nullglob
    matches=($glob)
    shopt -u nullglob

    if [ ${#matches[@]} -eq 0 ]; then
        echo "Error: could not locate ${description} (pattern: ${glob})" >&2
        echo "Transform file does not exist: 1" >&2
        exit 1
    fi

    printf '%s\n' "${matches[0]}"
}


# Function to check if a Conda environment is activated
function is_conda_env_activated() {
    # Check if the current environment name matches the desired one
    [[ "$CONDA_DEFAULT_ENV" == "$ENV_TO_ACTIVATE" ]]
}

function KUL_activate_conda_env {
    ENV_TO_ACTIVATE=$1
    # function to activate a certain conda env
    # if the base is active activate the requeste one
    # if another is active, go back to base
    echo "Current conda env: $CONDA_DEFAULT_ENV"
    #echo "Available environments:"
    #conda info --envs
    echo "Conda asked to activate: $ENV_TO_ACTIVATE"

    # Check if the desired Conda environment is activated
    if ! is_conda_env_activated; then
        echo "Activating the Conda environment: $ENV_TO_ACTIVATE"
        if [[ $ENV_TO_ACTIVATE == 'base' ]]; then 
            echo "Switching back to the base Conda environment: $ENV_TO_ACTIVATE"
            conda deactivate
        else
        # Activate the environment
            source activate "$ENV_TO_ACTIVATE"
            ACTIVATED=true
        fi
    else
        echo "Conda environment $ENV_TO_ACTIVATE is already activated."
        ACTIVATED=false
    fi
    
}

# MAIN FUNCTION - kul_echo ######################################################################################
# echo loud or silent
function kul_echo {
    local info_to_log=$1
    #echo "previous verbose_level: $verbose_level"
    if [[ -z "$verbose_level" ]]; then
        #echo "setting verbose_level = 2"
        verbose_level=2
    fi
    #if [ $silent -eq 0 ];then
    #    echo $1
    #fi
    if [ $verbose_level -eq 1 ]; then
        #echo "log: $log"
        echo ${info_to_log[@]} >> ${log}
    elif [ $verbose_level -eq 2 ]; then
        echo ${info_to_log[@]} | tee -a ${log}
    fi
}



# MAIN FUNCTION - function kul_e2cl to echo to console & log file with matlab tic/toc behavior ##################
# parameters for logging
log_every_seconds=120
function kul_e2cl {

    x=0
    if [ -z ${elapsed_s+x} ]; then
        #echo "var is unset" 
        elapsed_s=0
        old_elapsed_s=0
    fi
    
    local old_elapsed_s=$elapsed_s
    local b=''
    local n=''

    
    # calculate time in minutes since start of the script
    now=$(date +%s)
    elapsed_s=$(($now-$start))
    elapsed_m=$(echo "scale=2; $elapsed_s/60" | bc) # use bc to get it in float minutes
    diff_s=$(($elapsed_s-$old_elapsed_s))
    diff_m=$(echo "scale=2; $diff_s/60" | bc) # use bc to get it in float minutes
    if [ "$diff_s" -gt $log_every_seconds ]; then
        echo "          computation took ${diff_m} minutes ${b}(since start: ${elapsed_m} minutes)${n}"
        echo "          computation took ${diff_m} minutes (since start: ${elapsed_m})" >> $2
    fi

    # now do the logging
    echo "${b}$1${n}"
    echo $1 >> $2

}


function KUL_check_participant {
    #echo ${cwd}/BIDS/sub-${participant}
    if [ ! -d ${cwd}/BIDS/sub-${participant} ]; then
        echo "Error: participant ${participant} does not exists."
        exit 1
    fi
}


### MAIN ############################################################################################

# Upon sourcing this script these variables are set
kul_main_dir=$(dirname "$0")
cwd=$(pwd)
script_start_time=$SECONDS

# The KUL_DEBUG variable
# export KUL_DEBUG=1 in terminal to debug
if [[ -z "$KUL_DEBUG" ]]; then
    KUL_DEBUG=0
elif [ $KUL_DEBUG -eq 1 ]; then 
    echo "KUL_DEBUG is set to $KUL_DEBUG"
elif [ $KUL_DEBUG -gt 1 ]; then
    echo "KUL_DEBUG is now set to $KUL_DEBUG"
    set -x
fi



machine_type=$(uname)

# -- Execute global startup --

# timestamp
start=$(date +%s)

# Define directory/files to log in 
log_dir=${cwd}/KUL_LOG/$script
d=$(date "+%Y-%m-%d_%H-%M-%S")
log=$log_dir/main_log_${d}.txt

# create preprocessing & log directory/files
if [ ! -z "$1" ];then
    mkdir -p $log_dir
    # -- Say Welcome --
    command_line_options=$@
    #if []
    echo "Welcome to $script, version $version, invoked with parameters $command_line_options"

fi

# Check the mrtrix3 version
mrtrix_version_major=$(mrconvert | head -1 | cut -d'.' -f1 | cut -d' ' -f2)
mrtrix_version_minor=$(mrconvert | head -1 | cut -d'.' -f2)
mrtrix_version_revision_major=$(mrconvert | head -1 | cut -d'.' -f3 | cut -d'-' -f 1)
mrtrix_version_revision_minor=$(mrconvert | head -1 | cut -d'-' -f2)
if [ $mrtrix_version_revision_major -eq 2 ]; then
	mrtrix3new=0
	kul_echo "you are using an older version of MRTrix3 $mrtrix_version_revision_major"
	kul_echo "this is not supported. Exiting"
	exit 1
elif [ $mrtrix_version_revision_major -eq 3 ] && [ $mrtrix_version_revision_minor -lt 200 ]; then
	mrtrix3new=1
	kul_echo "you are using a new version of MRTrix3 $mrtrix_version_revision_major $mrtrix_version_revision_minor but not the latest"
elif [ $mrtrix_version_revision_major -ge 3 ] && [ $mrtrix_version_revision_minor -gt 200 ]; then
	mrtrix3new=2
	kul_echo "you are using the newest version of MRTrix3 $mrtrix_version_revision_major $mrtrix_version_revision_minor"
elif [ $mrtrix_version_revision_major -ge 3 ]; then
    Has_legacy=$(dwi2mask | head -18 | grep legacy | rev | cut -d "," -f3 | rev)
    if [ ! -z ${Has_legacy} ]; then
        mrtrix3new=2
	    kul_echo "you are using the newest version of MRTrix3 $mrtrix_version_revision_major $mrtrix_version_revision_minor"
    fi
else 
	kul_echo "cannot find correct mrtrix versions - exitting"
	exit 1
fi