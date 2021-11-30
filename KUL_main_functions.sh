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

    #local pidsArray=${task_in_pid[@]} # pids to wait for, separated by semi-colon
    #local procsArray=${task_in_name[@]} # name of procs to wait for, separated by semi-colon 
    local pidsArray=() # pids to wait for, separated by semi-colon
    local procsArray=() # name of procs to wait for, separated by semi-colon     
    local log_ttime=0 # local time instance for comparaison
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
            task_in_name[$local_n_tasks]="$(echo ${local_task_in:0:20} [instance $local_n_tasks])"
        else
            task_in_name[$local_n_tasks]="$kul_process_name [instance $local_n_tasks]"
        fi

        if [ ! -d "$local_main_logdir" ]; then
            mkdir "$local_main_logdir"
        fi

        if [[ -z "$kul_log_files" ]]; then 
            task_in_name_nospaces_tmp="${task_in_name[$local_n_tasks]// /_}"
            task_in_name_nospaces="${task_in_name_nospaces_tmp////_}"
            kul_log_file[$local_n_tasks]="$local_main_logdir/"${task_in_name_nospaces}".log"
            kul_errorlog_file[$local_n_tasks]="$local_main_logdir/"${task_in_name_nospaces}".error.log"
        else
            kul_log_file[$local_n_tasks]="$local_main_logdir/"${kul_log_files}_[$local_n_tasks].log""
            kul_errorlog_file[$local_n_tasks]="$local_main_logdir/"${kul_log_files}_[$local_n_tasks].error.log""
        fi

        ### STEP 2 - execute the task_in
        # to
        # implement multiple task_in (see preproc_all)

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
                log_min=$(echo "scale=2; $log_ttime/60" | bc)
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
                    fail_exec_time_minutes=$(echo "scale=2; $final_exec_time_seconds/60" | bc)
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
        #exit 1
    fi

    unset task_in

    # return errorcount

    if [[ ! -z $total_errorcount ]]; then
        total_errorcount=$(($total_errorcount + $errorcount))
        if [ $kul_verbose_level -eq 2 ]; then 
	        echo -e "total_errorcount: $total_errorcount\n\n\n"
        fi
    fi

}



# MAIN FUNCTION - kul_echo ######################################################################################
# echo loud or silent
function kul_echo {
    #echo "previous verbose_level: $verbose_level"
    if [[ -z "$verbose_level" ]]; then
        echo "setting verbose_level = 2"
        verbose_level=2
    fi
    #if [ $silent -eq 0 ];then
    #    echo $1
    #fi
    if [ $verbose_level -eq 1 ]; then
        #echo "log: $log"
        echo "$1" >> ${log}
    elif [ $verbose_level -eq 2 ]; then
        echo $1 | tee -a ${log}
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



### MAIN ############################################################################################

# Upon sourcing this script these variables are set
kul_main_dir=$(dirname "$0")
cwd=$(pwd)
script_start_time=$SECONDS

# The KUL_DEBUG variable
# export KUL_DEBUG=1 in terminal to debug
if [[ -z "$KUL_DEBUG" ]]; then
    KUL_DEBUG=0
elif [ $KUL_DEBUG -eq 1]; then
    set -x
fi


machine_type=$(uname)
#echo $machine_type

# Check the mrtrix3 version
mrtrix_version_major=$(mrconvert | head -1 | cut -d'.' -f1 | cut -d' ' -f2)
mrtrix_version_minor=$(mrconvert | head -1 | cut -d'.' -f2)
mrtrix_version_revision_major=$(mrconvert | head -1 | cut -d'.' -f3 | cut -d'-' -f 1)
mrtrix_version_revision_minor=$(mrconvert | head -1 | cut -d'-' -f2)
 
# -- Set global defaults --
#silent=1
#tmp=/tmp


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
