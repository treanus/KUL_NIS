#!/bin/bash -e
# Bash shell script to:
#  - define global functions used by all sub-scripts
#  - define defaults
#  - execute startup
#

# parameters for version_checking
mrtrix_version_needed=133
dcm2niix_version_needed=20180622
dcm2bids_version_needed=4
fsl_version_needed=6

# parameters for logging
log_every_seconds=120


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
    if [ "$diff_s" -gt $log_every_seconds ]; then
        echo "          computation took ${diff_m} minutes ${b}(since start: ${elapsed_m} minutes)${n}"
        echo "          computation took ${diff_m} minutes (since start: ${elapsed_m})" >> $2
    fi

    # now do the logging
    echo "${b}$1${n}"
    echo $1 >> $2

}

# check version of mrtrix3
mrtrix_version=$(mrconvert -version | head -n 1 | cut -d'-' -f 2)
if [ $mrtrix_version -lt $mrtrix_version_needed ]; then

    echo "Your mrtrix3 RC3 subversion is $mrtrix_version"
    echo "You need mrtrix3 RC3 subversion => $mrtrix_version_needed"
    exit 2

fi

# check version of dcm2niix
dcm2niix_version=$(dcm2niix | grep version | cut -d'.' -f 3 | cut -c -8)
if [ $dcm2niix_version -lt $dcm2niix_version_needed ]; then

    echo "Your version of dcm2nixx is $dcm2niix_version"
    echo "You need dcm2nixx version more recent than $dcm2niix_version_needed"
    exit 2

fi


# check version of dcm2bids
dcm2bids_version=$(dcm2bids -h | grep version | head -n 1 | cut -d'.' -f 2)
if [ $dcm2bids_version -lt $dcm2bids_version_needed ]; then

    echo "Your version of dcm2bids is $dcm2bids_version"
    echo "You need dcm2bids version equal or more than $dcm2bids_version_needed"
    exit 2

fi

# check version of fsl
fsl_version=$(flirt -version | cut -d' ' -f 3 | cut -d'.' -f 1)
if [ $fsl_version -lt $fsl_version_needed ]; then

    echo "Your version of FSL is $fsl_version"
    echo "You need FSL version equal or more than $fsl_version_needed"
    exit 2

fi

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
command_line_options=$@
kul_e2cl "Welcome to $script, version $v, invoked with parameters $command_line_options" $log
echo "   starting at $d"
