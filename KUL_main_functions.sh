#!/bin/bash -e
# Bash shell script to:
#  - define global functions used by all sub-scripts
#  - define defaults
#  - execute startup
#

# -- function kul_e2cl to echo to console & log file with matlab tic/toc behavior ---
function kul_e2cl {

    x=0
    if [ -z ${elapsed_s+x} ]; then
        #echo "var is unset" 
        elapsed_s=0
        old_elapsed_s=0
    fi
    
    local old_elapsed_s=$elapsed_s
    local b=$(tput bold)
    local n=$(tput sgr0)
    
    # calculate time in minutes since start of the script
    now=$(date +%s)
    elapsed_s=$(($now-$start))
    elapsed_m=$(echo "scale=2; $elapsed_s/60" | bc) # use bc to get it in float minutes
    diff_s=$(($elapsed_s-$old_elapsed_s))
    diff_m=$(echo "scale=2; $diff_s/60" | bc) # use bc to get it in float minutes
    if [ "$diff_s" -gt 1 ]; then
        echo "          computation took ${diff_m} minutes ${b}(since start: ${elapsed_m} minutes)${n}"
        echo "          computation took ${diff_m} minutes (since start: ${elapsed_m})" >> $2
    fi

    # now do the logging
    echo "${b}$1${n}"
    echo $1 >> $2

}





# -- Set global defaults --

silent=1
tmp=/tmp




# -- Execute global startup --

# timestamp
start=$(date +%s)

# Directory to write preprocessed data in, i.e $preproc
preproc=KUL_preproc/${subj}

# Define directory/files to log in 
log_dir=${preproc}/log/$script

# create preprocessing & log directory/files
mkdir -p $log_dir

# main log file naming
d=$(date "+%Y-%m-%d_%H-%M-%S")
log=$log_dir/main_log_${d}.txt


# -- Say Welcome --
kul_e2cl "Welcome to $script $v - $d" $log
