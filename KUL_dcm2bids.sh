#!/bin/bash -e
# Bash shell script to convert dicoms to bids format
#
# Requires dcm2bids (jooh fork), dcm2niix, Mrtrix3
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 26/10/2018 - alpha version
v="v0.1 - dd 26/10/2018"

# TODO
#  - make it work for multiple vendors
#  - wrap around for multiple subjects



# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl from KUL_main_functions (for logging)
#  - kul_dcmtags (for reading specific parameters from dicom header)

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh


# BEGIN LOCAL FUNCTIONS --------------

# --- function Usage ---
function Usage {

cat <<USAGE

`basename $0` converts dicoms to bids format

Usage:

  `basename $0` -d dicom_dir -s subject (participant) -c config_file -o bids_dir <OPT_ARGS>

  Depends on a config file that defines parameters with sequence information, e.g.

  Identifier,search-string,task,mb,pe_dir
  T1w,MPRAGE
  FLAIR,FLAIR
  func,rsfMRI,rest,8,j
  func,tb_fMRI,nback,2,j
  dwi,part1,-,2,j
  dwi,part2,-,2,j-

  explains that the T1w scan should be found by the search string "MPRAGE"
  func by rsfMRI, and has multiband_factor 8, and pe_dir = j

Example:

  `basename $0` -s pat001 -d pat001.zip -c definitions_of_sequences.txt -o BIDS

Required arguments:

     -d:  dicom_zip_file (the zip or tar.gz containing all your dicoms)
     -s:  subject (anonymised name of the participant)
     -c:  definitions of sequences (T1w=MPRAGE,dwi=seq, etc...)
     -o:  bids directory

Optional arguments:

     -t:  temporary directory (default = /tmp)
     -e:  copy task-*_events.tsv from config to BIDS dir
     -v:  verbose 

USAGE

    exit 1
}



# --- function kul_dcmtags (for reading specific parameters from dicom header & calculating missing BIDS parameters) ---
function kul_dcmtags {

    local dcm_file=$1 

    local out=${preproc}/log/${subj}_dcm_tags.txt
    local seriesdescr=$(dcminfo "$dcm_file" -tag 0008 103E | awk '{print $(NF)}')


    # 1/ Calculate ess/trt; needed are : FieldStrength, WaterFatShift, EPIFactor
    
        local waterfatshift=$(dcminfo "$dcm_file" -tag 2001 1022 | awk '{print $(NF)}')
        local fieldstrength=$(dcminfo "$dcm_file" -tag 0018 0087 | awk '{print $(NF)}')
        local epifactor=$(dcminfo "$dcm_file" -tag 2001 1013 | awk '{print $(NF)}')
        local water_fat_diff_ppm=3.3995
        local resonance_freq_mhz_tesla=42.576

        local water_fat_shift_hz=$(echo $fieldstrength $water_fat_diff_ppm $resonance_freq_mhz_tesla echo $fieldstrength $water_fat_diff_ppm | awk '{print $1 * $2 * $3}')
    
        #effective_echo_spacing_msec  = 1000 * WFS_PIXEL/(water_fat_shift_hz * (EPI_FACTOR + 1))
        ees_sec=$(echo $waterfatshift $water_fat_shift_hz $epifactor | awk '{print $1 / ($2 * ($3 + 1))}')

        #total_readout_time_fsl_msec      = EPI_FACTOR * effective_echo_spacing_msec;
        trt_sec=$(echo $epifactor $ees_sec | awk '{print $1 * $2 }')
    
    # 2/ Calculate slice SliceTiming

        #function SliceTime=KUL_slicetiming(MB, NS, TR)
        #%MB = 1; % Multiband factor
        #%NS = 30; % Number of Slices
        #%TR = 1.7; % TR in seconds
        #st = repmat(0:TR/(NS/MB):TR-.0000001,1,MB);
        #SliceTime = sprintf('%.8f,' , st);
        #SliceTime = ['"SliceTiming": [' SliceTime(1:end-1) ']'];
        #end     

        multiband_factor=$mb

        #NS
        local number_of_slices=$(dcminfo "$dcm_file" -tag 2001 1018 | awk '{print $(NF)}')

        #TR (in milliseconds)
        local repetion_time_msec=$(dcminfo "$dcm_file" -tag 0018 0080 | awk '{print $(NF)}')

        #single_slice_time (in seconds)
        local single_slice_time=$(echo $repetion_time_msec $number_of_slices $multiband_factor | awk '{print $1 / ($2 / $3) / 1000}')

    
        # number of excitations given multiband
        local e=$(echo $number_of_slices $multiband_factor | awk '{print ($1 / $2) -1 }')
        

        slit=0

        for (( c=1; c<=$e; c++ )); do 
            sl=$(echo $c $single_slice_time | awk '{print $1 * $2}')
            slit="$slit, $sl"
        done

        slit2=$slit
    
        rep=$(echo $multiband_factor | awk '{print $1 - 1}')

        for (( c=1; c<=$rep; c++ )); do 
            slit2="$slit2, $slit"
        done
    
        
        slice_time=[$slit2]
        

    if [ $silent -eq 0 ]; then
        echo "   the dicom file we are reading = $dcm_file"
        echo "     series = $seriesdescr"
        echo "      fieldstrength = $fieldstrength"
        echo "      waterfatshift =  $waterfatshift"
        echo "      epifactor = $epifactor"
        echo "        calulated ees  = $ees_sec"
        echo "        calculated trt  = $trt_sec"
        echo "      number of slices = $number_of_slices"
        echo "      repetion_time_msec = $repetion_time_msec"
        echo "      multiband_factor = $mb"
        echo "        calculated single_slice_time = $single_slice_time"
        echo "        number of excitations - 1 = $e"
        echo "        for 1 multiband = $slit"
        echo "        complete slice_time = $slice_time"
    fi

    if [ ! -f $out ]; then
        echo -e "series \t field \t epi \t wfs \t ees_sec \t trt_sec \t #slices \t repetion_time_msec \t multiband_factor" > $out
    fi
    echo -e "$seriesdescr \t $fieldstrength \t $epifactor \t $waterfatshift \t $ees_sec \t $trt_sec \t $number_of_slices \t $repetion_time_msec \t $multiband_factor" >> $out

    
}

function kul_untar {

    # find the search_string in the zip file with dicoms and extract

    echo "Searching for $identifier using $search_string"
            
    seq=$(grep $search_string ${preproc}/log/${subj}_dcm_content.txt | head -n 1 | cut -c 49-)
    seq_file=$(grep $search_string ${preproc}/log/${subj}_dcm_content.txt | awk 'NR==2' | cut -c 49-)
    if [ "$seq_file" = "" ]; then

        kul_e2cl "  $identifier dicoms are NOT FOUND" $log
        seq_found=0

    else

        seq_found=1
        echo "  $identifier dicoms are in ${seq}"

        echo "  Untarring $dcm to $tmp"
        mkdir -p ${tmp}/$subj
        tar -C ${tmp}/$subj -xzf ${dcm} "${seq}"

    fi

}


# END LOCAL FUNCTIONS --------------




# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
silent=1
tmp=/tmp

# Set required options
subj_flag=0
dcm_flag=0
conf_flag=0
bids_flag=0
tmp_flag=0
events_flag=0

if [ "$#" -lt 3 ]; then
    Usage >&2
    exit 1

else

    while getopts "d:s:c:o:vet" OPT; do

        case $OPT in
        s) #subject
            subj_flag=1
            subj=$OPTARG
        ;;
        c) #config_file
            conf_flag=1
            conf=$OPTARG
        ;;
        d) #dicom_zip_file
            dcm_flag=1
            dcm=$OPTARG
        ;;
        o) #bids output directory
            bids_flag=1
            bids_output=$OPTARG
        ;;
        v) #verbose
            silent=0
        ;;
        e) #events for task based fmri
            events_flag=1
        ;;
        t) #temporary directory
            tmp_flag=1
            tmp=$OPTARG
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
if [ $dcm_flag -eq 0 ] ; then 
    echo 
    echo "Option -d is required: give the file with your raw dicoms either .zip of tar.gz" >&2
    echo
    exit 2 
fi 

if [ $subj_flag -eq 0 ] ; then 
    echo 
    echo "Option -s is required: give the anonymised name of a subject this will create a directory subject_preproc with results." >&2
    echo
    exit 2 
fi 

if [ $conf_flag -eq 0 ] ; then 
    echo 
    echo "Option -c is required: give the path to the file that describes the sequences" >&2
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
    echo "Option -o is required: give the path to the BODS directory" >&2
    echo
    exit 2 
fi 

# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Directory to write preprocessed data in
preproc=${subj}_preproc

# set up preprocessing & logdirectory
mkdir -p ${preproc}/log

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=${preproc}/log/log_${d}.txt



# ----------- SAY HELLO ----------------------------------------------------------------------------------

kul_e2cl "Welcome to KUL_dcm2bids $v - $d" $log

if [ $silent -eq 0 ]; then
    echo "  The script you are running has basename `basename "$0"`, located in dirname $kul_main_dir"
    echo "  The present working directory is `pwd`"
fi

# dump the content of the zip file with dicoms
kul_e2cl "  reading the compressed dicom file $dcm" $log
tar -ztvf ${dcm} > ${preproc}/log/${subj}_dcm_content.txt


# we will untar the relevant bits of dicoms using the config file
# create empty bids description
bids=""

# we read the config file
while IFS=, read identifier search_string task mb pe_dir; do
    
    if [ $identifier = "T1w" ]; then 
        
        kul_untar

        sub_bids=$(cat <<EOF
        {
        "dataType": "anat",
        "suffix": "T1w",
        "criteria": {
            "in": {
            "SeriesDescription": "${search_string}",
            "ImageType": "ORIGINAL"
                }
            }
        })
        #echo $sub_bids

        bids="$bids,$sub_bids"

    fi

    if [ $identifier = "FLAIR" ]; then 
        
        kul_untar

        sub_bids=$(cat <<EOF
        {
        "dataType": "anat",
        "suffix": "FLAIR",
        "criteria": {
            "in": {
            "SeriesDescription": "${search_string}",
            "ImageType": "ORIGINAL"
                }
            }
        })
        #echo $sub_bids

        bids="$bids,$sub_bids"

    fi

    if [ $identifier = "func" ]; then 

        kul_untar

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${tmp}/$subj/${seq_file}"

            sub_bids=$(cat <<EOF
            {
            "dataType": "func",
            "suffix": "bold",
            "criteria": {
                "in": {
                "SeriesDescription": "${search_string}",
                "ImageType": "ORIGINAL"
                }
            },
            "customHeader": {
                "TaskName": "${task}",
                "EffectiveEchoSpacing": ${ees_sec},
                "TotalReadoutTime": ${trt_sec},
                "MultibandAccelerationFactor": ${mb},
                "PhaseEncodingDirection": "${pe_dir}",
                "SliceTiming": $slice_time
            }
            })
            #echo $sub_bids

            bids="$bids,$sub_bids"

        fi

    fi

    if [ $identifier = "dwi" ]; then 

        kul_untar

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${tmp}/$subj/${seq_file}"

            sub_bids=$(cat <<EOF
            {
            "dataType": "dwi",
            "suffix": "dwi",
            "criteria": {
                "in": {
                "SeriesDescription": "${search_string}",
                "ImageType": "ORIGINAL"
                }
            },
            "customHeader": {
                "EffectiveEchoSpacing": ${ees_sec},
                "TotalReadoutTime": ${trt_sec},
                "MultibandAccelerationFactor": ${mb},
                "PhaseEncodingDirection": "${pe_dir}",
                "SliceTiming": $slice_time
            }
            })
            
            bids="$bids,$sub_bids"

        fi

    fi


done < $conf

# make the full bids_conf string and write it to file
bids_conf=$(cat <<EOF
{
  "descriptions": [
      ${bids:1}
        ]
}
EOF)

echo $bids_conf  > ${preproc}/log/${subj}_bids_config.json 

# invoke dcm2bids
kul_e2cl "Calling dcm2bids..." $log
dcm2bids  -d "${tmp}/$subj" -p $subj -c ${preproc}/log/${subj}_bids_config.json -o $bids_output
rm -rf $bids_output/tmp_dcm2bids

# copying task based events.tsv to BIDS directory
if [ $events_flag -eq 1 ]; then
    kul_e2cl "Copying task based events.tsv to BIDS directory" $log
    cp config/task-*_events.tsv $bids_output
fi

# clean up
rm -rf "${tmp}/$subj"

kul_e2cl "Finished" $log