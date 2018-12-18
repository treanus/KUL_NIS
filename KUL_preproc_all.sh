#!/bin/bash -e
# @ Stefan Sunaert & Ahmed Radwan- UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 06/11/2018 - alpha version
# v0.2a - dd 12/12/2018 - preparing for beta release 0.2
v="v0.2a - dd 12/12/2018"

# This is the main script of the KUL_NeuroImaging_Toools
#
# Description:
#    This script preprocces an entire study (multiple subjects) with structural, functional and diffusion data at Stefan's lab
#      It will:
#       - convert dicom files to BIDS format
#       - perform mriqc on structural and functional data
#       - perform fmriprep on structural and functional data
#       - perform freesurfer on the structural data (only T1w for now) 
#       - perform mrtix3 and related processing on dMRI data
#       - optionally:
#           - perform combined structural and dMRI data analysis (depends on fmriprep, freesurfer and mrtrix3 above)
#           - perfrom dbsdrt (automated tractography of the dentato-rubro-thalamic tract) on dMRI + structural data (depends on all above)
#
#   Requirements:
#       A correct installation of your mac (for now, maybe later also a hpc) at the lab
#           - including:
#               - dcm2niix (in KUL_apps)
#               - dcm2bids (jooh fork, using pip)
#               - docker
#               - freesurfer (in KUL_apps)
#               - mrtrix (in KUL_apps)  
#               - last but not least, a correct installation of up-to-date KUL_NeuroImaging_Tools (in KUL_apps)
#               - correct setup of your .bashrc and .bash_profile
#
#  It depends on a major config file, e.g. "study_config/subjects_and_options.csv" in which one informs the script:
#       What and how (options) to perform:
#               - mriqc (yes/no) 
#                   (no options implemented yet)
#               - fmriprep (yes/no), and specifies options:
#                   all fmriprep options may be given,
#                   e.g.:
#                       --anat-only (to only process structural)
#               - freesurfer (yes/no) 
#                   (no options implemented yet)
#               - KUL_dwiprep processing, i.e. a full mrtrix processing pipeline (yes/no)
#                   options may be e.g.:
#                       --slm=linear --repol (to provide to eddy)
#               - KUL_dwiprep_anat processing (yes/no) 
#                   (no options implemented yet)
#               - KUL_dwiprep_dbsdrt processing (yes/no)
#                       option nods e.g. 4000
#






# To do:
# - update the description above (section DESCRIPTION) of what this script does exactly!
#
#       - other ideas:
#               - add KUL_dcm2bids in the loop of processing (was implemented, but temporarily out again)
#               - add processing for fmri stats
#               - add processing for automated tracking of major tracts (similar to tractseg e.g.)





# Source KUL_main_functions
# KUL_main_functions will:
#  - Say Welcome
#  - check wether all necessary software is installed (and exit if needed)
#  - provide some general functions like logging
kul_main_dir=$(dirname "$0")
source $kul_main_dir/KUL_main_functions.sh
script=$(basename "$0")
cwd=$(pwd)




# Start with defining local functions


# A Function to provide Usage information
#   - gives information about the script
function Usage {

cat <<USAGE

`basename $0` preproccesses an entire study

Usage:

  `basename $0` -c config_file -b bids_dir

Required arguments:

     -c:  description of the subjects and settings for processing
     -b:  bids directory

Optional arguments:
    
     -n:  number of cores to use (distrubuted over mriqc/fmriprep/freesurfer/etc...)
     -m:  max memory (in gigabytes) available in docker
     -t:  temporary directory (default = /tmp)
     -r:  reset docker (clean the images and download new ones)
     -v:  verbose

Example:

  `basename $0` -c study_config/subjects_and_options.csv -b BIDS -n 6 -m 12 -t /scratch -v 
    
    uses "study_config/subjects_and_options.csv" to 
        - reads the subjects (participants) on which to do processing
        - reads what processing (mriqc/fmriprep/freesurfer/etc...) to do on those
        - reads the options to give to mriqc/fmriprep/etc...
    uses the (already converted) BIDS data in directory "BIDS"
    uses 6 cores in total for all processes (distrubuted over mriqc/fmriprep/freesurfer/etc...)
    uses & memory of 12 GB
        - for fmriprep & mriqc 
        - (set this option equal to, or slightly less than what you specify in your docker prefences)
    specifies that temporary data are written to /scratch
    spits out more verbose logging to the terminal

USAGE

    exit 1
}


# A Function to start mriqc processing (in parallel)
function task_mriqc_participant {

# check if already performed mriqc
mriqc_dir_to_check=mriqc/sub-${BIDS_participant}

if [ ! -d $mriqc_dir_to_check ]; then

    mriqc_log=${preproc}/log/mriqc/${BIDS_participant}.txt

    kul_e2cl " started (in parallel) mriqc on participant $BIDS_participant (using $ncpu_mriqc cores, logging to $mriqc_log)" $log

    local task_mriqc_cmd=$(echo "docker run --read-only --tmpfs /run --tmpfs /tmp --rm \
 -v ${cwd}/${bids_dir}:/data:ro -v ${cwd}/mriqc:/out \
 poldracklab/mriqc:latest \
 --participant_label $BIDS_participant \
 --n_procs $ncpu_mriqc --ants-nthreads $ncpu_mriqc_ants --mem_gb $mem_gb --no-sub \
 /data /out participant \
 > $mriqc_log 2>&1 ") 

    echo "   using cmd: $task_mriqc_cmd"

    # now we start the parallel job
    eval $task_mriqc_cmd &
    mriqc_pid="$!"
    echo " mriqc pid is $mriqc_pid"

    sleep 2

    #kul_e2cl "   done mriqc on participant $BIDS_participant" $log

else
        
    echo " mriqc of participant $BIDS_participant already done, skipping..."

fi

}


# A function to start fmriprep processing (in parallel)
function task_fmriprep {

# check if already performed fmriprep
fmriprep_file_to_check=fmriprep/sub-${BIDS_participant}.html

if [ ! -f $fmriprep_file_to_check ]; then

    fmriprep_log=${preproc}/log/fmriprep/${BIDS_participant}.txt

    kul_e2cl " started (in parallel) fmriprep on participant ${BIDS_participant}... (with options $fmriprep_options, using $ncpu_fmriprep cores, logging to $fmriprep_log)" ${log}

    local task_fmriprep_cmd=$(echo "docker run --rm \
 -v ${cwd}/${bids_dir}:/data \
 -v ${cwd}:/out \
 -v ${cwd}/fmriprep_work:/scratch \
 -v /KUL_apps/freesurfer/license.txt:/opt/freesurfer/license.txt \
 poldracklab/fmriprep:latest \
 --participant_label ${BIDS_participant} \
 -w /scratch \
 --nthreads $ncpu_fmriprep --omp-nthreads $ncpu_fmriprep_ants \
 --mem_mb $mem_mb \
 --fs-no-reconall \
 $fmriprep_options \
 --notrack \
 /data /out \
 participant \
 > $fmriprep_log  2>&1") 

    echo "   using cmd: $task_fmriprep_cmd"

    # Now start the parallel job
    eval $task_fmriprep_cmd &
    fmriprep_pid="$!"
    echo " fmriprep pid is $fmriprep_pid"

    sleep 2
    
    #kul_e2cl "   done fmriprep on participant $BIDS_participant" $log

else
        
    echo " fmriprep of participant $BIDS_participant already done, skipping..."

fi

}

# A Function to start freesurfer processing (in parallel)
function task_freesurfer {

# check if already performed freesurfer
freesurfer_file_to_check=freesurfer/sub-${BIDS_participant}/${BIDS_participant}/scripts/recon-all.done
        
if [ ! -f  $freesurfer_file_to_check ]; then

    freesurfer_log=${preproc}/log/freesurfer/${BIDS_participant}.txt

    kul_e2cl " started (in parallel) freesurfer recon-all on participant ${BIDS_participant}... (using $ncpu_freesurfer cores, logging to $freesurfer_log)" ${log}

    mkdir -p freesurfer

    bids_subj=${bids_dir}/sub-${BIDS_participant}/ses-tp1
    bids_anat=$(ls ${bids_subj}/anat/*_T1w.nii.gz)
    #echo $bids_anat

    SUBJECTS_DIR=${cwd}/freesurfer/sub-${BIDS_participant}

    #start clean
    rm -rf $SUBJECTS_DIR
    mkdir -p $SUBJECTS_DIR
    export SUBJECTS_DIR

    local task_freesurfer_cmd=$(echo "recon-all -subject $BIDS_participant -i $bids_anat -all -openmp $ncpu_freesurfer \
 -parallel > $freesurfer_log 2>&1 ")

    echo "   using cmd: $task_freesurfer_cmd"

    eval $task_freesurfer_cmd &
    freesurfer_pid="$!"
    echo "   freesurfer pid is $freesurfer_pid"
    
    sleep 2

    #kul_e2cl "   done freesufer on participant $BIDS_participant" $log

else

    echo " freesurfer of subjet $BIDS_participant already done, skipping..."
        
fi

}



# A function to start KUL_dwiprep processing (in parallel)
function task_KUL_dwiprep {

# check if already performed KUL_dwiprep
dwiprep_file_to_check=dwiprep/sub-${BIDS_participant}/qa/dec.mif

#FLAG we still need to implement topup_options

if [ ! -f  $dwiprep_file_to_check ]; then

    dwiprep_log=${preproc}/log/dwiprep/dwiprep_${BIDS_participant}.txt

    kul_e2cl " started (in parallel) KUL_dwiprep on participant ${BIDS_participant}... (using $ncpu_dwiprep cores, logging to $dwiprep_log)" ${log}

    local task_dwiprep_cmd=$(echo "KUL_dwiprep.sh -p ${BIDS_participant} -n $ncpu_dwiprep -d $dwipreproc_options -e \"${eddy_options} \" -v \
 > $dwiprep_log 2>&1 ")

    echo "   using cmd: $task_dwiprep_cmd"

    # Now we start the parallel job
    eval $task_dwiprep_cmd &
    dwiprep_pid="$!"
    echo " KUL_dwiprep pid is $dwiprep_pid"

    sleep 2

    #kul_e2cl "   done KUL_dwiprep on participant $BIDS_participant" $log

else

    echo " KUL_dwiprep of participant $BIDS_participant already done, skipping..."
        
fi

}




# A Function to start KUL_dwiprep_anat processing
function task_KUL_dwiprep_anat {

# check if already performed KUL_dwiprep_anat
dwiprep_anat_file_to_check=dwiprep/sub-${BIDS_participant}/qa/dec_reg2T1w_on_t1w.mif

if [ ! -f  $dwiprep_anat_file_to_check ]; then

    dwiprep_anat_log=${preproc}/log/dwiprep/dwiprep_anat_${BIDS_participant}.txt

    kul_e2cl " performing KUL_dwiprep_anat on subject ${BIDS_participant}... (using $ncpu cores, logging to $dwiprep_anat_log)" ${log}

    KUL_dwiprep_anat.sh -p ${BIDS_participant} -n $ncpu -v \
        > $dwiprep_anat_log 2>&1 

    kul_e2cl "   done KUL_dwiprep_anat on participant $BIDS_participant" $log

else

    echo " KUL_dwiprep_anat of subjet $BIDS_participant already done, skipping..."
        
fi

}




# A Function to start KUL_dwiprep_drtdbs processing
function task_KUL_dwiprep_drtdbs {

# check if already performed KUL_dwiprep_drtdbs
dwiprep_drtdbs_file_to_check=dwiprep/sub-${BIDS_participant}/qa/void

if [ ! -f  $dwiprep_drtdbs_file_to_check ]; then

    dwiprep_drtdbs_log=${preproc}/log/dwiprep/dwiprep_drtdbs_${BIDS_participant}.txt

    kul_e2cl " performing KUL_dwiprep_drtdbs on subject ${BIDS_participant}... (using $ncpu cores, logging to $dwiprep_drtdbs_log)" ${log}

    local task_dwiprep_drtdbs_cmd=$(echo "KUL_dwiprep_drtdbs.sh -p ${BIDS_participant} -n $ncpu -v -o $drtdbs_options -v")

    echo "   using cmd: $task_dwiprep_drtdbs_cmd"
    
    eval $task_dwiprep_drtdbs_cmd

    kul_e2cl "   done KUL_dwiprep_drtdbs on participant $BIDS_participant" $log

else

    echo " KUL_dwiprep_drtdbs of subjet $BIDS_participant already done, skipping..."
        
fi

}

# end of local function --------------





# MAIN STARTS HERE

# Set some defaults
silent=1
mem_gb=16
tmp=/tmp

# Set flags
conf_flag=0
bids_flag=0
tmp_flag=0
cpu_flag=0
mem_flag=0
docker_reset_flag=0

# Check command line options, and return function Usage if required options are not given
if [ "$#" -lt 2 ]; then
    Usage >&2
    exit 1

else

    while getopts "c:b:n:m:t:rvh" OPT; do

        case $OPT in
        c) #config_file
            conf_flag=1
            conf=$OPTARG
        ;;
        b) #bids_dir
            bids_flag=1
            bids_dir=$OPTARG
        ;;
        n) #ncpu
            cpu_flag=1
            ncpu=$OPTARG
        ;;
        m) #bids_dir
            mem_flag=1
            mem_gb=$OPTARG
        ;;
        t) #temporary directory
            tmp_flag=1
            tmp=$OPTARG
        ;;
        r) #reset docker
            docker_reset_flag=1
        ;;
        v) #verbose
            silent=0
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
if [ $conf_flag -eq 0 ] ; then 
    echo 
    echo "Option -c is required: give the path to the file that describes the subjects" >&2
    echo
    exit 2 

elif [ ! -f $conf ] ; then
    echo 
    echo "The config file $conf does not exist"
    echo
    exit 2
fi 

if [ $bids_flag -eq 0 ] ; then 
    echo 
    echo "Option -o is required: give the path to the BIDS directory" >&2
    echo
    exit 2 
fi 

if [ $cpu_flag -eq 0 ] ; then 
    echo 
    echo "Option -p is required: give the number of cpu's to use" >&2
    echo
    exit 2 
fi 


# INITIATE ---



# ----------- MAIN ----------------------------------------------------------------------------------

if [ $silent -eq 0 ]; then
    echo "  The script you are running has basename `basename "$0"`, located in dirname $kul_main_dir"
    echo "  The present working directory is `pwd`"
fi

# ---------- SET MAIN DEFAULTS ---
# set mem_mb for mriqc
gb=1024
mem_mb=$(echo $mem_gb $gb | awk '{print $1 * $2 }')

# ---------- PROCESS CONTROL & LOAD BALANCING --------
# We will be running 4 preprocessings in parallel: mriqc, fmriprep, freesurfer & KUL_dwiprep
# We need to do some load balancing #FLAG, needs optimisation, a.o. if some processes finished already!

# set number of cores for task mriqc
load_mriqc=33 # higher number means less cpu need (mriqc does not need much)
ncpu_mriqc=$(((($ncpu/$load_mriqc))+1))
ncpu_mriqc_ants=$(((($ncpu/$load_mriqc))+1))

# set number of cores for task fmriprep
load_fmriprep=5
ncpu_fmriprep=$(((($ncpu/$load_fmriprep))+1))
ncpu_fmriprep_ants=$(((($ncpu/$load_fmriprep))+1))

# set number of cores for task freesurfer
load_freesurfer=3
ncpu_freesurfer=$(((($ncpu/$load_freesurfer))+1))

# set number of cores for task KUL_dwiprep
load_dwiprep=2
ncpu_dwiprep=$(((($ncpu/$load_dwiprep))+1))


# Ask if docker needs to be reset
if [ $docker_reset_flag -eq 1 ];then
    docker system prune -a
fi


# ----------- STEP 1 - CONVERT TO BIDS ---
#kul_e2cl "Performing KUL_multisubjects_dcm2bids... " $log
#KUL_multisubjects_dcm2bids.sh -d DICOM -c $conf -o $bids_dir -e
# TODO:
#  - this will be changed: KUL_multisubjects_dcm2bids will become obsolete
#  - instead we will call KUL_dcm2bids for each subject in the loop below.


# ----------- STEP 2 - Preprocess each subject with mriqc, fmriprep, freesurfer and KUL_dwiprep ---
# set up logging directories and clean left over fmriprep_work directory
# TODO:
#  - this should best go into the task_*
mkdir -p ${preproc}/log/mriqc
mkdir -p ${preproc}/log/freesurfer
mkdir -p ${preproc}/log/dwiprep
mkdir -p ${preproc}/log/fmriprep
rm -fr ${cwd}/fmriprep_work



# we read the config file (and it may be csv, tsv or ;-seperated)
while IFS=$'\t,;' read -r BIDS_participant EAD dicom_zip config_file do_mriqc mriqc_options do_fmriprep fmriprep_options do_freesurfer freesurfer_options do_dwiprep dwipreproc_options topup_options eddy_options do_dwiprep_anat anat_options do_dwiprep_drtdbs drtdbs_options; do
    
    
    if [ "$dicom_zip" = "dicom_zip" ]; then
        
        echo "first line" > /dev/null 2>&1

    else

        kul_e2cl "Performing preprocessing of subject $BIDS_participant... " $log

        # check if already performed dcm2bids
        dcm2bids_dir_to_check=BIDS/sub-${BIDS_participant}
        
        if [ ! -d  $dcm2bids_dir_to_check ]; then

            kul_e2cl " Converting dicom to BIDS for subject $BIDS_participant... " $log
            KUL_dcm2bids.sh -p ${BIDS_participant} -d ${dicom_zip} -c ${config_file} -o BIDS -v
        
        else

            echo " BIDS conversion for subject $BIDS_participant already done, skipping... "

        fi
        
        kul_e2cl "  Now starting (depending on your config-file) mriqc, fmriprep, freesurfer and KUL_dwiprep... " $log
        echo "   note: further processing with KUL_dwiprep_anat, KUL_dwiprep_drtdbs depend on fmriprep, freesurfer and KUL_dwiprep (which need to run fully)"

        if [ $silent -eq 0 ]; then
            
            echo "  if this script fails, please check your configuration file (given to -c); for now this was what was defined:"
            echo "    BIDS_participant: $BIDS_participant"
            echo "    EAD: $EAD"
            echo "    dicom_zip: $dicom_zip"
            echo "    config_file: $config_file"
            echo "    do_mriqc: $do_mriqc"
            echo "    mriqc_options: $mriqc_options"
            echo "    do_fmriprep: $do_fmriprep"
            echo "    fmriprep_options: $fmriprep_options"
            echo "    do_freesurfer: $do_freesurfer"
            echo "    freesurfer_options: $freesurfer_options"
            echo "    do_dwiprep: $do_dwiprep"
            echo "    dwipreproc_options: $dwipreproc_options"
            echo "    topup_options: $topup_options"
            echo "    eddy_options: $eddy_options"
            echo "    do_dwiprep_anat: $do_dwiprep_anat"
            echo "    anat_options: $anat_options"
            echo "    do_dwiprep_drtdbs: $do_dwiprep_drtdbs"
            echo "    drtdbs_options: $drtdbs_options"
        
        fi


        if [ $do_mriqc -eq 1 ]; then
            
            task_mriqc_participant 

        fi

        if [ $do_fmriprep -eq 1 ]; then
            
            task_fmriprep 

        fi

        if [ $do_dwiprep -eq 1 ]; then
            
            task_KUL_dwiprep

        fi

        if [ $do_freesurfer -eq 1 ]; then
            
            task_freesurfer 

        fi

        # wait for mriqc, fmriprep, freesurfer and KUL_dwiprep to finish
        kul_e2cl " waiting for processes mriqc, fmriprep, freesurfer and KUL_dwiprep for subject $BIDS_participant to finish before continuing with further processing... (this can take hours!)... " $log
        wait $mriqc_pid $fmriprep_pid $dwiprep_pid $freesurfer_pid

        kul_e2cl " processes mriqc, fmriprep, freesurfer and KUL_dwiprep for subject $BIDS_participant have finished" $log

        # clean up after jobs finished
        rm -fr ${cwd}/fmriprep_work


        # Here we could also have fMRI statistical analysis e.g.
        # task_KUL_fmri_model # needs to be made

        # Here we could also have rsfMRI processing e.g.
        # task_KUL_fmri_melodic # needs to be made

        # continue with KUL_dwiprep_anat, which depends on finished data from freesurfer, fmriprep & KUL_dwiprep
        if [ $do_dwiprep_anat -eq 1 ]; then
            task_KUL_dwiprep_anat
        fi 

        # Here we could also have some whole brain tractography processing e.g.
        # task_KUL_mrtix_wb_tckgen # needs to be made
        # task_KUL_mrtrix_tractsegment # needs to be made
        
        # continue with KUL_dwiprep_drtdbs
        if [ $do_dwiprep_drtdbs -eq 1 ]; then
            echo $drtdbs_options
            task_KUL_dwiprep_drtdbs
        fi

    fi

done < $conf


# ----------- STEP 3 - Compute mriqc group summary ---

kul_e2cl "Performing mriqc group summary" $log

# check if already performed group mriqc
mriqc_file_to_check=mriqc/group_bold.html

if [ ! -f $mriqc_file_to_check ]; then
    
    docker run --read-only --tmpfs /run --tmpfs /tmp --rm \
        -v ${cwd}/${bids_dir}:/data:ro -v ${cwd}/mriqc:/out \
        poldracklab/mriqc:latest \
        /data /out group

else

    kul_e2cl " group mriqc already done, skipping..." $log

fi

# ----------- STEP 4 - Compute dwiprep group summary ---

cat dwiprep/sub-*/tracts_info.csv > dwiprep/group_tracts_info.csv




kul_e2cl "Finished all... " $log
