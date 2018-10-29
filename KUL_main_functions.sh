#!/bin/bash -e
# Bash shell script to define main used functions
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
        echo "${b}          computation took ${diff_m} minutes${n} (since start: ${elapsed_m} minutes)"
        echo "          computation took ${diff_m} minutes (since start: ${elapsed_m})" >> $2
    fi

    # now do the logging
    echo "${b}$1${n}"
    echo $1 >> $2

}

