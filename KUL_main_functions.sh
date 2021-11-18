#!/bin/bash
# Bash shell script to:
#  - define global functions used by all sub-scripts
#  - define defaults
#  - execute startup
#


# Function task_exec
# Inputs
#  - obligatory
#       task_in (a command string that needs to be evaluated)
#  - facultative
#       kul_verbous_level (1= default)
#       kul_log_file (the path of the log file)
function KUL_task_exec {

    local kul_verbous_level=$1

    if [[ -z "$kul_verbous_level" ]]; then
        kul_verbous_level=1
    fi
    
    if [[ -z "$kul_log_file" ]]; then 
        kul_log_file="/dev/null"
    fi

    ### TODO
    # implement multiple task_in (see preproc_all)
    task_in_short=$(echo ${task_in:0:20})
    eval ${task_in} | tee -a ${kul_log_file} &
    task_in_pid="$!"

    if [ $kul_verbous_level -gt 0 ]; then 
        echo -e "\n${task_in}" | tee -a ${kul_log_file} 
        echo "Starting [\"${task_in_short}...\"] @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${kul_log_file} 
    fi

    local pidsArray=${task_in_pid[@]} # pids to wait for, separated by semi-colon
    local procsArray=${task_in_short[@]} # name of procs to wait for, separated by semi-colon     
    #local exit_on_error="false"
    #local soft_alert=0 # Does a soft alert need to be triggered, if yes, send an alert once 
    local log_ttime=0 # local time instance for comparaison
    local seconds_begin=$SECONDS # Seconds since the beginning of the script
    local exec_time=0 # Seconds since the beginning of this function
    #local retval=0 # return value of monitored pid process
    local errorcount=0 # Number of pids that finished with errors
    local pidCount # number of given pids
    local c # counter for pids/procsArray

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
                    echo "  *** WARNING! **** Process ${procsArray[c]} with pid $pid FAILED (with exitcode [$result]). Check the log-file"
                else
                    if [ $kul_verbous_level -gt 0 ]; then 
                        echo "Process ${procsArray[c]} with pid $pid finished successfully @ $(date "+%Y-%m-%d_%H-%M-%S") (with exitcode [$result])."
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
                echo "  Current tasks [${procsArray[@]}] still running after $log_min minutes with pids [${pidsArray[@]}]."
            fi
        fi

        pidsArray=("${newPidsArray[@]}")
        procsArray=("${newProcsArray[@]}")
        sleep 1

    done

    if [ $errorcount -eq 0 ]; then
        if [ $kul_verbous_level -gt 0 ]; then 
            echo Success | tee -a ${kul_log_file}
        fi
    else
        echo Fail | tee -a ${kul_log_file}
        exit 1
    fi

    unset task_in

}


# parameters for logging
log_every_seconds=120

# echo loud or silent
function kul_echo {
    if [ $silent -eq 0 ];then
        echo $1
    fi
}

# -- function kul_e2cl to echo to console & log file with matlab tic/toc behavior ---
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
preproc=KUL_LOG/${subj}

# Define directory/files to log in 
log_dir=${preproc}/log/$script
d=$(date "+%Y-%m-%d_%H-%M-%S")
log=$log_dir/main_log_${d}.txt

# create preprocessing & log directory/files
if [ ! -z "$1" ];then
    mkdir -p $log_dir
    # -- Say Welcome --
    command_line_options=$@
    kul_e2cl "Welcome to $script, version $version, invoked with parameters $command_line_options" $log
    echo "   starting at $d"
fi
