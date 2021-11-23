#!/bin/bash
# Bash shell script to:
#  - define global functions used by all /DATA/fmri_pats/BIDS/derivatives/KUL_compute/sub-Casier/FastSurfer/sub-Casier/scripts/lh.processing.cmdfsub-scripts
#  - define defaults
#  - execute startup
#
# @ prof.sunaert@gmail.com - v0.2 - 18/11/2021


# MAIN FUNCTION - Function task_exec ###################################################################################
#  - obligatory variable to set
#       task_in (a command string that needs to be evaluated)
#  - facultative command line options
#       kul_verbose_level (0=silent, 1=normal; 2=verbose; 1=default)
#       kul_short_name (what to display as process)
#       kul_log_file (the path of the log file)
# Example 1
#  task_in="a_big_long_process 0 1 "KUL VBG" KUL_LOG/my_big_long_process"
#   Run "process a_big_long_process"
#   Do it silently
#   Display "KUL VBG" as the running process
#   Log to the files KUL_LOG/my_big_long_process.log and KUL_LOG/my_big_long_process.error.log
# Example 2
#  task_in="another big process"
#   Run process "another_big_process"
#   Have output to the terminal and to log files (default kul_verbose_level 1)
#   Take the first 20 characters of "task_in" and display that (kul_short_name is generated automatically)
#   Log to the files KUL_LOG_kul_short_name.log and KUL_LOG_kul_short_name.error.log (/ -> _, and spaces too)
function KUL_task_exec {

    local kul_verbose_level="$1"
    local kul_process_name="$2"
    local kul_log_files="$3"

    # export KUL_DEBUG=1 in terminal to debug
    if [[ -z "$KUL_DEBUG" ]]; then
        KUL_DEBUG=0
    fi

    if [ $KUL_DEBUG -eq 1 ]; then
        echo $kul_verbose_level
        echo $kul_process_name
        echo $kul_log_files
        echo $script
    fi

    #local pidsArray=${task_in_pid[@]} # pids to wait for, separated by semi-colon
    #local procsArray=${task_in_name[@]} # name of procs to wait for, separated by semi-colon 
    local pidsArray=() # pids to wait for, separated by semi-colon
    local procsArray=() # name of procs to wait for, separated by semi-colon     
    local log_ttime=0 # local time instance for comparaison
    local seconds_begin=$SECONDS # Seconds since the beginning of the script
    local exec_time=0 # Seconds since the beginning of this function
    local errorcount=0 # Number of pids that finished with errors
    local pidCount # number of given pids
    local c # counter for pids/procsArray
    
    ### STEP 1 - test the input to the function and of unset, set a default
    if [[ -z "$kul_verbose_level" ]]; then
        kul_verbose_level=1
    fi

    if [[ -z "$kul_process_name" ]]; then
        task_in_name="$(echo ${task_in:0:20})"
    else
        task_in_name="$kul_process_name"
    fi

    if [[ -z "$kul_log_files" ]]; then 
        task_in_name_nospaces="${task_in_name// /_}"
        task_in_name_nospaces="${task_in_name_nospaces////_}"
        kul_log_file="KUL_LOG/${script}/"${task_in_name_nospaces}".log"
        kul_errorlog_file="KUL_LOG/${script}/"${task_in_name_nospaces}".error.log"
    else
        kul_log_file="${kul_log_files}.log"
        kul_errorlog_file="${kul_log_files}.error.log"
    fi

    ### STEP 2 - execute the task_in
    # to
    # implement multiple task_in (see preproc_all)

    if [ $kul_verbose_level -eq 0 ]; then 
    
        local task_in_final="$task_in  1>>${kul_log_file} 2>>${kul_errorlog_file}"
        eval ${task_in_final} &
    
    else
        
        local task_in_final="$task_in  > >(tee -a ${kul_log_file}) 2> >(tee -a ${kul_errorlog_file})"
        eval ${task_in_final} &
         
    fi
    
    # set the pids, first we get the pid by "$!", then feed it into an array
    task_in_pid="$!"
    pidsArray+=($task_in_pid)
    procsArray+=("$task_in_name")


    ### STEP 3 - give some information
    if [ $kul_verbose_level -gt 0 ]; then
        tput bold
        echo "${task_in_name}... started @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${kul_log_file}
        tput sgr0
    fi
    if [ $kul_verbose_level -eq 2 ]; then 
        tput dim
        echo -e "   The task_in command: ${task_in}" | tee -a ${kul_log_file}
        tput sgr0
    fi


    ### STEP 4 - keep checking if the process is still running    
    pidCount=${#pidsArray[@]}
    #echo "  pidCount: $pidCount"
    #echo "  pidsArray: ${pidsArray[@]}"

    while [ ${#pidsArray[@]} -gt 0 ]; do

        newPidsArray=()
        newProcsArray=()
        c=0

        for pid in "${pidsArray[@]}"; do
            #echo "pid: $pid"
            #echo "proc: ${procsArray[c]}"
            if kill -0 $pid > /dev/null 2>&1; then
                newPidsArray+=($pid)
                #echo "newPidsArray: ${newPidsArray[@]}"
                newProcsArray+=(${procsArray[c]})
                #echo "newProcsArray: ${newProcsArray[@]}"
            else
                wait $pid
                result=$?
                #echo "result: $result"
                if [ $result -ne 0 ]; then

                    errorcount=$((errorcount+1))
                    fail_exec_time_seconds=$(($SECONDS - $seconds_begin))
                    fail_exec_time_minutes=$(($final_exec_time_seconds / 60))
                    tput bold; tput setaf 1
                    echo "  *** WARNING! **** Process ${procsArray[c]} with pid $pid might have failed after $fail_exec_time_minutes minutes. (with exitcode [$result]). Check the ${kul_errorlog_file} log-file" | tee -a ${kul_errorlog_file}
                    tput sgr0
                
                else
                    
                    final_exec_time_seconds=$(($SECONDS - $seconds_begin))
                    final_exec_time_minutes=$(($final_exec_time_seconds / 60))
                    if [ $kul_verbose_level -gt 0 ]; then
                        echo "$task_in_name finished successfully after $final_exec_time_minutes minutes" | tee -a ${kul_log_file}
                    fi
                fi
            fi
            c=$((c+1))
        done

        ## Log a standby message every hour
        every_time=1201
        exec_time=$(($SECONDS - $seconds_begin))
        if [ $((($exec_time + 1) % $every_time)) -eq 0 ]; then
            if [ $log_ttime -ne $exec_time ]; then
                log_ttime=$exec_time
                log_min=$((log_ttime / 60))
                if [ $kul_verbose_level -gt 0 ]; then
                    echo "  Current tasks [${procsArray[@]}] still running after $log_min minutes with pids [${pidsArray[@]}]."
                fi
            fi
        fi

        pidsArray=("${newPidsArray[@]}")
        procsArray=("${newProcsArray[@]}")
        sleep 1

    done

    ### STEP 5 - return the status of execution 
    if [ $errorcount -eq 0 ]; then
        if [ $kul_verbose_level -eq 2 ]; then 
            echo -e "Success\n\n" | tee -a ${kul_log_file}
        fi
    else
        echo "Fail" | tee -a ${kul_log_file}
        #exit 1
    fi

    unset task_in

    # return errorcount

}



# MAIN FUNCTION - kul_echo ######################################################################################
# echo loud or silent
function kul_echo {
    if [ $silent -eq 0 ];then
        echo $1
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
    #local b=$(tput bold)
    #local n=$(tput sgr0)
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



### OTHER STUFF ############################################################################################
machine_type=$(uname)
#echo $machine_type

# Check the mrtrix3 version
mrtrix_version_major=$(mrconvert | head -1 | cut -d'.' -f1 | cut -d' ' -f2)
mrtrix_version_minor=$(mrconvert | head -1 | cut -d'.' -f2)
mrtrix_version_revision_major=$(mrconvert | head -1 | cut -d'.' -f3 | cut -d'-' -f 1)
mrtrix_version_revision_minor=$(mrconvert | head -1 | cut -d'-' -f2)
 
# -- Set global defaults --
silent=1
tmp=/tmp


# -- Execute global startup --

# timestamp
start=$(date +%s)

# Directory to write preprocessed data in, i.e $preproc
preproc=KUL_LOG

# Define directory/files to log in 
log_dir=${preproc}/$script
d=$(date "+%Y-%m-%d_%H-%M-%S")
log=$log_dir/main_log_${d}.txt

# create preprocessing & log directory/files
if [ ! -z "$1" ];then
    mkdir -p $log_dir
    # -- Say Welcome --
    command_line_options=$@
    #kul_e2cl "Welcome to $script, version $version, invoked with parameters $command_line_options" $log
    kul_echo "Welcome to $script, version $version, invoked with parameters $command_line_options"
    #echo "   starting at $d"
fi
