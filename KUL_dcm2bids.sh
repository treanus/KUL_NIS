#!/bin/bash -e

# set -x

# Bash shell script to convert dicoms to bids format
#
# Requires dcm2bids, dcm2niix, Mrtrix3
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - KUL - ahmed.radwan@kuleuven.be
#
v="v0.8 - dd 08/09/2020"

# Notes
#  - NOW USES https://github.com/UNFmontreal/Dcm2Bids
#  - works for GE/Siemens/Philips
#  - wrap around for multiple subjects: use KUL_multisubjects_dcm2bids


# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl from KUL_main_functions (for logging)
#  - kul_dcmtags (for reading specific parameters from dicom header)

# source general functions
kul_main_dir=`dirname "$0"`
script=`basename "$0"`
source "${kul_main_dir}/KUL_main_functions.sh"
cwd=$(pwd)

# BEGIN LOCAL FUNCTIONS --------------

# --- function Usage ---
function Usage {

cat <<USAGE

`basename $0` converts dicoms to bids format

Usage:

  `basename $0` -d dicom_dir -p subject (participant) -c config_file -o bids_dir <OPT_ARGS>

  Depends on a config file that defines parameters with sequence information, e.g.

  For Philips dicom we need to manually specify the mb and pe_dir:
    # Identifier,search-string,fmritask/inteded_for,mb,pe_dir,acq_label
    # Structural scans
    T1w,T1_PRE
    cT1w,T1_Post
    FLAIR,3D_FLAIR
    T2w,3D_T2
    SWI,SWI
    MTI,mtc_2dyn
    # functional scans
    func,rsfMRI_MB6,rest,6,j,singleTE
    sbref,rsfMRI_SBREF,rest,1,j,singleTE
    func,MB_mTE,rest,4,j,multiTE
    sbref,mTE_SBREF,rest,4,j,multiTE
    func,MB2_hand,HAND,2,j
    func,MB2_lip,LIP,2,j
    func,MB2_nback,nback,2,j
    sbref,MB2_SBREF,HAND nback,1,j
    # fmap: 'task' is now 'IntendedFor' (the func images/tasks that the B0_map is used for SDC)
    fmap,B0_map,[HAND LIP nback]
    # dMRI
    dwi,p1_b1200,-,3,j-,b1200
    dwi,p2_b0,-,3,j-,b0
    dwi,p3_b2500,-,3,j-,b2500
    dwi,p4_b2500,-,3,j,rev
    # ASL support is very limited (not BIDS compliant for now)
    # Indentifier, search-string 
    ASL,pCASL

  explains that the T1w scan should be found by the search string "T1_PRE"
  func by rsfMRI, and has multiband_factor 6, and pe_dir = j
  the sbref will be used for the tb_fMRI for both tasks (hands and nback)
  
  For Siemens and GE dicom it can be as simple as:
    # Identifier,search-string,fmritask/inteded_for,mb,pe_dir,acq_label
    # Structural scans
    T1w,T1_PRE
    cT1w,T1_Post
    FLAIR,3D_FLAIR
    T2w,3D_T2
    # functional scans
    func,rsfMRI_MB6,rest,-,-,singleTE
    sbref,rsfMRI_SBREF,rest,-,-,singleTE
    func,MB_mTE,rest,-,-,multiTE
    sbref,mTE_SBREF,rest,-,-,multiTE
    func,MB2_hand,HAND,-,-
    func,MB2_lip,LIP,-,-
    func,MB2_nback,nback,-,-
    sbref,MB2_SBREF,HAND nback,-,-
    # dMRI
    dwi,p1_b1200,-,-,-,b1200
    dwi,p2_b0,-,-,-,b0
    dwi,p3_b2500,-,-,-,b2500
    dwi,p4_b2500,-,-,-,rev
    # fmap, SWI, MTC and ASL have not been tested on Siemens or GE
    
Example:

  `basename $0` -p pat001 -d pat001.zip -c definitions_of_sequences.txt -o BIDS

Required arguments:

     -d:  dicom_zip_file (the zip or tar.gz containing all your dicoms, or directory containing dicoms)
     -p:  participant (anonymised name of the subject in bids convention)
     -c:  definitions of sequences (T1w=MPRAGE,dwi=seq, etc..., see above)

Optional arguments:

     -o:  bids directory
     -s:  session (for longitudinal study with multiple timepoints)
     -t:  temporary directory (default = /tmp)
     -e:  copy task-*_events.tsv from config to BIDS dir
     -a:  further anonymise the subject by using pydeface (takes much longer)
     -v:  verbose 

USAGE

    exit 1
}

# check if jsontool is installed and install it if not
if [[ $(which jsontool) ]]; then

    echo "  jsontool already installed, good" $log

else

    echo "  jsontool not installed, installing it with pip using pip install jsontool" $log
    pip install jsontool

fi
# check if pydeface is installed and install it if not
if [[ $(which pydeface) ]]; then

    echo "  pydeface already installed, good" $log

else

    echo "  pydeface not installed, installing it with pip using pip install jsontool" $log
    pip install pydeface

fi
# check if the correct dcm2bids is installed and install it if not
if [[ $(which dcm2bids_scaffold) ]]; then

    echo "  dcm2bids already installed, good" $log

else

    echo "  dcm2bids not installed, installing it with pip using pip install dcm2bids" $log
    pip uninstall Dcm2Bids
    pip install dcm2bids

fi

# --- function kul_dcmtags (for reading specific parameters from dicom header & calculating missing BIDS parameters) ---
function kul_dcmtags {

    local dcm_file=$1 

    local out=$final_dcm_tags_file
    
    # 0/ Read out standard tags for logging
    local seriesdescr=$(dcminfo "$dcm_file" -tag 0008 103E | cut -c 13-)
    local manufacturer=$(dcminfo "$dcm_file" -tag 0008 0070 | cut -c 13-)
    local software=$(dcminfo "$dcm_file" -tag 0018 1020 | cut -c 13-)
    local imagetype=$(dcminfo "$dcm_file" -tag 0008 0008 2>/dev/null | cut -c 13- | head -n 1)
    local patid=$(dcminfo "$dcm_file" -tag 0010 0020 | cut -c 13-)
    local pixelspacing=$(dcminfo "$dcm_file" -tag 0028 0030 2>/dev/null | cut -c 13- | head -n 1)
    local slicethickness=$(dcminfo "$dcm_file" -tag 0018 0088 2>/dev/null | cut -c 13- | head -n 1)
    local acquisitionMatrix=$(dcminfo "$dcm_file" -tag 0018 1310 2>/dev/null | cut -c 13- | head -n 1)
    local FovAP=$(dcminfo "$dcm_file" -tag 2005 1074 | cut -c 13-)
    local FovFH=$(dcminfo "$dcm_file" -tag 2005 1075 | cut -c 13-)
    local FovRL=$(dcminfo "$dcm_file" -tag 2005 1076 | cut -c 13-)
    # local echonumber=$(dcminfo "$dcm_file" -tag 0018 0086 2>/dev/null | cut -c 13- | head -n 1)
    # need to add local echonumber or something similar for mTE (0018,0086)     

    # Now we need to determine what vendor it is.
    #   Philips needs all the following calculations
    #   Siemens works out of the box
    #   GE most recent version also seem to work fine

    #echo $manufacturer
    if [ "$manufacturer" = "SIEMENS" ]; then

        slicetime_provided_by_vendor=1
        ees_trt_provided_by_vendor=1

    #elif [ "$manufacturer" = 'GE ]

        # need to be tested

    else

        slicetime_provided_by_vendor=0
        ees_trt_provided_by_vendor=0

    fi


    # 1/ Calculate ess/trt; needed are : FieldStrength, WaterFatShift, EPIFactor
         
    # Note: only calculate it when it is provided (sometimes this has been thrown away by anonymising the dicom-data)
    # We check whether needed tags exist
        
    tags_are_present=1

    test_waterfatshift=$(dcminfo "$dcm_file" -tag 2001 1022)
    if [ -z "$test_waterfatshift" ]; then
        tags_are_present=0
        local waterfatshift="empty"
    else
        local waterfatshift=$(dcminfo "$dcm_file" -tag 2001 1022 | awk '{print $(NF)}')
    fi

    test_fieldstrength=$(dcminfo "$dcm_file" -tag 0018 0087)
    if [ -z "$test_fieldstrength" ]; then
        tags_are_present=0
        local fieldstrength="empty"
    else
        local fieldstrength=$(dcminfo "$dcm_file" -tag 0018 0087 | awk '{print $(NF)}')
    fi

    test_epifactor=$(dcminfo "$dcm_file" -tag 2001 1013)
    if [ -z "$test_epifactor" ]; then
        tags_are_present=0
        local epifactor="empty"
    else
        local epifactor=$(dcminfo "$dcm_file" -tag 2001 1013 | awk '{print $(NF)}')
    fi


    if [ $tags_are_present -eq 0 ]; then

        ees_sec="empty"
        trt_sec="empty"

    else

        local water_fat_diff_ppm=3.3995
        local resonance_freq_mhz_tesla=42.576

        local water_fat_shift_hz=$(echo $fieldstrength $water_fat_diff_ppm $resonance_freq_mhz_tesla echo $fieldstrength $water_fat_diff_ppm | awk '{print $1 * $2 * $3}')
    
        #effective_echo_spacing_msec  = 1000 * WFS_PIXEL/(water_fat_shift_hz * (EPI_FACTOR + 1))
        ees_sec=$(echo $waterfatshift $water_fat_shift_hz $epifactor | awk '{print $1 / ($2 * ($3 + 1))}')

        #total_readout_time_fsl_msec      = EPI_FACTOR * effective_echo_spacing_msec;
        trt_sec=$(echo $epifactor $ees_sec | awk '{print $1 * $2 }')
    
    fi

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

    # if mb is not found in the config file or is set to 0
    # simply set it to 1 and run
    if [[ -z ${multiband_factor} ]] || [[ ${multiband_factor} == 0 ]]; then
        multiband_factor=1
    fi
    

    tags_are_present=1

    test_number_of_slices=$(dcminfo "$dcm_file" -tag 2001 1018)
    if [ -z "$test_number_of_slices" ]; then
        tags_are_present=0
        local number_of_slices="empty"
    else
        local number_of_slices=$(dcminfo "$dcm_file" -tag 2001 1018 | awk '{print $(NF)}')
    fi

    test_repetion_time_msec=$(dcminfo "$dcm_file" -tag 0018 0080 )
    if [ -z "$test_repetion_time_msec" ]; then
        tags_are_present=0
        local repetion_time_msec="empty"
    else
        local repetion_time_msec=$(dcminfo "$dcm_file" -tag 0018 0080 | awk '{print $(NF)}')
    fi
    
    test_slice_scan_order=$(dcminfo "$dcm_file" -tag 2005 1081 )
    if [ -z "$test_slice_scan_order" ]; then
        #tags_are_present=0
        local slice_scan_order="empty"
    else
        local slice_scan_order=$(dcminfo "$dcm_file" -tag 2005 1081 | head -n 1 | awk '{print $(NF)}')
    fi
    

    if [ $tags_are_present -eq 1 ]; then

        if [ ! $multiband_factor = "" ];then

            #single_slice_time (in seconds)
            local single_slice_time=$(echo $repetion_time_msec $number_of_slices $multiband_factor | awk '{print $1 / ($2 / $3) / 1000}')
    
            # number of excitations given multiband
            # e = n. of excitations/slices per band
            local e=$(echo $number_of_slices $multiband_factor | awk '{print ($1 / $2) -1 }')
            local spb=$((${e}+1));
            echo $e
            echo $spb

            echo "${slice_scan_order}"
        
            # here we need to adapt to account for different slice orders
            # e.g. 
            if [[ "${slice_scan_order}" == "rev. central" ]]; then
                # this is a bit different from interleaved... namely we split it into 2 gps
                # lower group is regular ascending and second group is regular descending
                half_e=$(echo "scale=2;(${spb}/2)" | bc | awk '{print int($1+0.5)}')
                
                for (( zc=0; zc<${half_e}; zc++ )); do 
                    sl1=$(echo $zc $single_slice_time | awk '{print $1 * $2}')
                    if [[ -z ${slit1} ]]; then 
                        slit1="${sl1}"
                    else
                        slit1="${slit1}, ${sl1}"
                    fi
                done

                for (( zx=${e}; zx>=${half_e}; zx-- )); do 
                    sl2=$(echo $zx $single_slice_time | awk '{print $1 * $2}')
                    if [[ -z ${slit2} ]]; then 
                        slit2="${sl2}"
                    else
                        slit2="${slit2}, ${sl2}"
                    fi
                done

                slit="${slit1}, ${slit2}"
                echo "${slit}"

                for iz in ${!tmp_order[@]}; do
                    sl=$(echo $((${tmp_order[$iz]})) ${single_slice_time} | awk '{print $1 * $2}');
                    echo ${sl}
                    if [[ -z ${slit} ]]; then 
                        slit="${sl}"
                        echo ${slit}
                    else
                        slit="${slit}, ${sl}"
                        echo ${slit}
                    fi
                done
                
                # interleaved is ready!
            elif [[ "${slice_scan_order}" == "interleaved" ]]; then
                
                declare -a tmp_order
                step=$(echo "sqrt(${spb})" | bc)
                unset tmp_order curr bh ik slgp; 
                slgp=0; 
                declare -a tmp_order; 
                tmp_order[0]=0;
                for ik in $(seq 1 ${e}); do 
                    bh=$((${ik}-1)); 
                    curr=$((${tmp_order[$bh]}+${step})); 
                    if [[ ${curr} -gt ${e} ]]; then 
                        ((slgp++)); 
                        curr=${slgp}; 
                    fi; 
                    tmp_order[$ik]=${curr}; 
                    
                done; 

                for iz in ${!tmp_order[@]}; do
                    sl=$(echo $((${tmp_order[$iz]})) ${single_slice_time} | awk '{print $1 * $2}');
                    echo ${sl}
                    if [[ -z ${slit} ]]; then 
                        slit="${sl}"
                        echo ${slit}
                    else
                        slit="${slit}, ${sl}"
                        echo ${slit}
                    fi
                done
                
            elif [[ "${slice_scan_order}" == "FH" ]]; then

                for (( c=0; c<=$e; c++ )); do 
                    sl=$(echo $c $single_slice_time | awk '{print $1 * $2}')
                    if [[ -z ${slit} ]]; then 
                        slit="${sl}"
                    else
                        slit="${slit}, ${sl}"
                    fi
                done

            elif [[ "${slice_scan_order}" == "HF" ]]; then

                for (( c=${e}; c>=0; c-- )); do 
                    sl=$(echo $c $single_slice_time | awk '{print $1 * $2}')
                    if [[ -z ${slit} ]]; then 
                        slit="${sl}"
                    else
                        slit="${slit}, ${sl}"
                    fi
                done

            elif [[ "${slice_scan_order}" == "default" ]]; then

                echo ${slice_scan_order}
                a=$((${spb} %2));
                echo ${a}
                # this is still untested
                if [[ "${spb}" -le 6 ]]; then
                    step=2;
                    hlpp=$(($((${e}+1))/${step}));
                    lpsl=0;
                    lpal=${lpsl};
                    hpsl=${e};
                    hpal=${hpsl};
                    order=0;

                    for ii in $(seq 0 1 ${spb}); do
                        if [[ ${lpal} -lt $((${hlpp}-1)) ]]; then
                            tmp_order[${order}]=${lpal};
                            lpal=$((${lpal}+${step}));
                            ((order++))
                        elif [[ ${hpal} -ge $((${hlpp}-1)) ]]; then
                            tmp_order[${order}]=${hpal};
                            hpal=$((${hpal}-${step}));
                            ((order++))
                        else
                            lpal=$((${lpsl}+1))
                            hpal=$((${hpsl}-1))
                        fi
                    done

                    # We will not add a 1 as done in the matlab version but we iterate over spb not e also
                    for iz in ${!tmp_order[@]}; do
                        sl=$(echo $((${tmp_order[$iz]})) ${single_slice_time} | awk '{print $1 * $2}');
                        if [[ -z ${slit} ]]; then 
                            slit="${sl}"
                        else
                            slit="${slit}, ${sl}"
                        fi
                    done

                # this is still untested
                elif [[ "${spb}" == 8 ]]; then
                    declare -a tmp_order
                    step=$(echo "sqrt(${spb})" | bc)
                    echo "step is ${step}"
                    unset tmp_order curr bh ik slgp; 
                    slgp=0; 
                    declare -a tmp_order; 
                    tmp_order[0]=0;
                    echo ${tmp_order[@]}; 
                    for ik in $(seq 1 ${e}); do 
                        echo " ik is ${ik}"; 
                        bh=$((${ik}-1)); 
                        echo "bh is ${bh}"; 
                        curr=$((${tmp_order[$bh]}+${step})); 
                        echo "tmp_order of bh is ${tmp_order[$bh]}; 
                        echo "initially tmp_order of ik is ${tmp_order[$ik]}; 
                        echo "step is ${step}"; 
                        if [[ ${curr} -gt ${e} ]]; then 
                            echo "slgp is ${slgp}"; 
                            ((slgp++)); 
                            echo "inceremented slgp is ${slgp}"; 
                            curr=${slgp}; 
                            echo "now curr = ${slgp}"; 
                        fi; 
                        tmp_order[$ik]=${curr}; 
                        echo "tmp_order of ik is ${tmp_order[$ik]}"; 
                    done; 
                    echo ${tmp_order[@]}

                    for iz in ${!tmp_order[@]}; do
                        sl=$(echo $((${tmp_order[$iz]})) ${single_slice_time} | awk '{print $1 * $2}');
                        if [[ -z ${slit} ]]; then 
                            slit="${sl}"
                        else
                            slit="${slit}, ${sl}"
                        fi
                    done

                else

                    # if we have an odd no. of slices per band
                    if [[ ${a} == 1 ]]; then
                        
                        for zc in $(seq 0 2 ${e} ); do 
                            sl1=$(echo $zc $single_slice_time | awk '{print $1 * $2}')
                            if [[ -z ${slit1} ]]; then 
                                slit1="${sl1}"
                            else
                                slit1="${slit1}, ${sl1}"
                            fi
                        done

                        for zx in $(seq 1 2 ${e} ); do 
                            sl2=$(echo $zx $single_slice_time | awk '{print $1 * $2}')
                            if [[ -z ${slit2} ]]; then 
                                slit2="${sl2}"
                            else
                                slit2="${slit2}, ${sl2}"
                            fi
                        done

                        slit="${slit1}, ${slit2}"
                        echo "${slit}"
                    
                    # if we have an odd no. of slices per band
                    elif [[ ${a} == 0 ]]; then
                        step=2;
                        hlpp=$(($((${e}+1))/${step}));
                        lpsl=0;
                        lpal=${lpsl};
                        hpsl=${e};
                        hpal=${hpsl};
                        order=0;

                        for ii in $(seq 0 1 ${spb}); do
                            if [[ ${lpal} -lt $((${hlpp}-1)) ]]; then
                                tmp_order[${order}]=${lpal};
                                lpal=$((${lpal}+${step}));
                                ((order++))
                            elif [[ ${hpal} -ge $((${hlpp}-1)) ]]; then
                                tmp_order[${order}]=${hpal};
                                hpal=$((${hpal}-${step}));
                                ((order++))
                            else
                                lpal=$((${lpsl}+1))
                                hpal=$((${hpsl}-1))
                            fi
                        done

                        # We will not add a 1 as done in the matlab version but we iterate over spb not e also
                        for iz in ${!tmp_order[@]}; do
                            sl=$(echo $((${tmp_order[$iz]})) ${single_slice_time} | awk '{print $1 * $2}');
                            if [[ -z ${slit} ]]; then 
                                slit="${sl}"
                            else
                                slit="${slit}, ${sl}"
                            fi
                        done

                    fi
                fi

            fi

            slit2=$slit
    
            rep=$(echo $multiband_factor | awk '{print $1 - 1}')

            for (( c=1; c<=$rep; c++ )); do 
                slit2="$slit2, $slit"
            done
    
        
            slice_time=[$slit2]
        
        fi
    
    else

        slice_time="empty"

    fi


    if [ $silent -eq 0 ]; then
        echo "   patid = $patid"
        echo "   the dicom file we are reading = $dcm_file"
        echo "   manufacturer = $manufacturer"
        echo "   software version = $software"
        echo "   imagetype = $imagetype"
        echo "   acquisitionMatrix = $acquisitionMatrix"
        echo "   FovAP = $FovAP"
        echo "   FovFH = $FovFH"
        echo "   FovRL = $FovRL"
        echo "   pixelspacing = $pixelspacing"
        echo "   slicethickness = $slicethickness"
        echo "     series = $seriesdescr"
        echo "      fieldstrength = $fieldstrength"
        echo "      waterfatshift =  $waterfatshift"
        echo "      epifactor = $epifactor"
        echo "        calulated ees  = $ees_sec"
        echo "        calculated trt  = $trt_sec"
        echo "      number of slices = $number_of_slices"
        echo "      repetion_time_msec = $repetion_time_msec"
        echo "      slice_scan_order = $slice_scan_order"
        
        if [ ! $multiband_factor = "" ];then
            echo "      multiband_factor = $mb"
            echo "        calculated single_slice_time = $single_slice_time"
            echo "        number of excitations - 1 = $e"
            echo "        for 1 multiband = $slit"
            echo "        complete slice_time = $slice_time"
        fi

    fi

    if [ ! -f $out ]; then
        echo -e "participant,session,dcm_file,manufacturer,software_version,series_descr,imagetype,fieldstrength,acquisitionMatrix,FovAP,FovFH,FovRL,pixelspacing,slicethickness,epifactor,wfs,ees_sec,trt_sec,nr_slices,slice_scan_order,repetion_time_msec,multiband_factor" > $out
    fi
    echo -e "$subj,${sess},$dcm_file,$manufacturer,$software,$seriesdescr,$imagetype,$fieldstrength,$acquisitionMatrix,$FovAP,$FovFH,$FovRL,$pixelspacing,$slicethickness,$epifactor,$waterfatshift,$ees_sec,$trt_sec,$number_of_slices,$slice_scan_order,$repetion_time_msec,$multiband_factor" >> $out

    
}

function kul_find_relevant_dicom_file {

    kul_e2cl "  Searching for ${identifier} using search_string $search_string" $log

    # find the search_string in the dicom dump_file            
    # search for search_string in dump_file, find ORIGINAL, remove dicom tags, sort, take first line, remove trailing space 
    seq_file=$(grep $search_string $dump_file | grep ORIGINAL - | cut -f1 -d"[" | sort | head -n 1 | sed -e 's/[[:space:]]*$//')

    if [ "$seq_file" = "" ]; then

        kul_e2cl "    ${identifier} dicoms are NOT FOUND" $log
        seq_found=0

    else

        seq_found=1
        kul_e2cl "    a relevant ${identifier} dicom is $(basename "${seq_file}") " $log

    fi

}


# END LOCAL FUNCTIONS --------------




# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
sess=""
bids_output=BIDS

# Set flags
subj_flag=0
sess_flag=0
dcm_flag=0
conf_flag=0
bids_flag=0
tmp_flag=0
events_flag=0
n_sbref_tasks=0
silent=1
anon=0

if [ "$#" -lt 4 ]; then
    Usage >&2
    exit 1

else

    while getopts "c:d:p:o:s:t:aveh" OPT; do

        case $OPT in
        d) #dicom_zip_file
            dcm_flag=1
            dcm=$OPTARG
        ;;
        p) #participant
            subj_flag=1
            subj=$OPTARG
        ;;
        c) #config_file
            conf_flag=1
            conf=$OPTARG
        ;;
        o) #bids output directory
            bids_flag=1
            bids_output=$OPTARG
        ;;
        s) #session
            sess_flag=1
            sess=$OPTARG
        ;;
        a) #pydeface
            anon=1
        ;;
        v) #verbose
            silent=0
        ;;
        e) #events for task based fmri
            events_flag=1
        ;;
        t) #temporary directory
            tmp_flag=1
            tempo=$OPTARG
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
    echo "Option -p is required: give the anonymised name of a subject this will create a directory subject_preproc with results." >&2
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


# INITIATE ---

# main log file naming
d=$(date "+%Y-%m-%d_%H-%M-%S")
log=$log_dir/${subj}_main_log_${d}.txt

# file for initial dicom tags
dump_file=${log_dir}/${subj}_${sess}_initial_dicom_info.txt

# file with final dicom tags
final_dcm_tags_file=${log_dir}/${subj}_${sess}_final_dicom_info.csv

# location of bids_config_json_file
bids_config_json_file=${log_dir}/${subj}_${sess}_bids_config.json

# location of dcm2niix_log_file
dcm2niix_log_file=$log_dir/${subj}_${sess}_dcm2niix_log_file.txt


 if [[ $tmp_flag -eq 1 ]] ; then

    tmp=${cwd}/${tempo}

 else

    tmp="/tmp/${subj}"

 fi

rm -fr ${tmp}

# exit

# remove previous existances to start fresh
rm -f $dump_file
rm -f $final_dcm_tags_file
rm -f $bids_config_json_file
# rm -fr ${cwd}/${tmp}/$subj

# ----------- SAY HELLO ----------------------------------------------------------------------------------


if [[ $silent -eq 0 ]]; then
    echo "  The script you are running has basename `basename "$0"`, located in dirname $kul_main_dir"
    echo "  The present working directory is `pwd`"
fi

# uncompress the zip file with dicoms or link the directory to tmp

# clear the /tmp directory

mkdir -p ${tmp}

if [[ -d "$dcm" ]]; then

    # it is a directory
    kul_e2cl "  you gave the directory $dcm as input; linking to to $tmp/$subj" $log
    ln -s "${cwd}/${dcm}" $tmp/$subj

else

    kul_e2cl "  uncompressing the zip file $dcm to $tmp/$subj" $log
    # Check the extention of the archive
    arch_ext="${dcm##*.}"
    #echo $arch_ext


    if [[ $arch_ext = "zip" ]]; then 

        unzip -q -o ${dcm} -d ${tmp}

    elif [[] $arch_ext = "tar" ]]; then

        tar --strip-components=5 -C ${tmp} -xzf ${dcm}

    fi

fi

# dump the dicom tags of all dicoms in a file
kul_e2cl "  brute force extraction of some relevant dicom tags of all dicom files of subject $subj into file $dump_file" $log

# check if a DICOMIR file exists and archive it (we do a brute force extraction and don't need the DICOMDIR file)
test_DICOMDIR=($(find -L ${tmp} -type f -name "DICOMDIR" ))
#echo "DICOMDIR = $test_DICOMDIR"
if [[ $test_DICOMDIR == "" ]]; then
    echo "  OK there is no DICOMDIR"
else
    gzip $test_DICOMDIR
fi

# Do the bruce force extract
echo hello > $dump_file

task(){
    dcm1=$(dcminfo "$dcm_file" -tag 0008 103E -tag 0008 0008 -tag 0008 0070 -nthreads 4 2>/dev/null | tr -s '\n' ' ')
    echo "$dcm_file" $dcm1 >> $dump_file
}

(
find -L ${tmp} -type f | 
while IFS= read -r dcm_file; do
    task 
done
)

kul_e2cl "    done reading dicom tags of $dcm" $log


# create empty bids description
# bids=""

# we read the config file

declare -a sub_bids

while IFS=, read identifier search_string task mb pe_dir acq_label; do

    bs=$(( $bs + 1))


 if [[ ! ${identifier} == \#* ]]; then

    if [[ ${identifier} == "T1w" ]]; then 
        
        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_T1='{"dataType": "anat", "modalityLabel": "T1w", "criteria": {  
             "SeriesDescription": "*'${search_string}'*"}}'

            sub_bids_[$bs]=$(echo ${sub_bids_T1} | python -m json.tool )

        fi

    fi

    if [[ ${identifier} == "cT1w" ]]; then 
        
        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_T1='{"dataType": "anat", "modalityLabel": "T1w", "criteria": {  
             "SeriesDescription": "*'${search_string}'*"},
            "customLabels": "ce-gadolinium",
            "sidecarChanges": {"KUL_dcm2bids": "yes","ContrastBolusIngredient": "gadolinium"}
            }'

            sub_bids_[$bs]=$(echo ${sub_bids_T1} | python -m json.tool )

        fi

    fi

    if [[ ${identifier} == "T2w" ]]; then 
        
        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_T2='{"dataType": "anat", "modalityLabel": "T2w", "criteria": { 
             "SeriesDescription": "*'${search_string}'*"}}'

            sub_bids_[$bs]=$(echo ${sub_bids_T2} | python -m json.tool)

        fi

    fi

    if [[ ${identifier} == "PD" ]]; then 
        
        kul_find_relevant_dicom_file
The complete ASL

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_PD='{"dataType": "anat", "modalityLabel": "PD", "criteria": {  
             "SeriesDescription": "*'${search_string}'*"}}'

            sub_bids_[$bs]=$(echo ${sub_bids_PD} | python -m json.tool)

        fi

    fi

    if [[ ${identifier} == "SWI" ]]; then 
        
        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_SWI='{"dataType": "anat", "modalityLabel": "SWI", "criteria": {  
             "SeriesDescription": "*'${search_string}'*"}}'

            sub_bids_[$bs]=$(echo ${sub_bids_SWI} | python -m json.tool)

        fi

    fi

    if [[ ${identifier} == "MTI" ]]; then 
        
        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_SWI='{"dataType": "anat", "modalityLabel": "MTI", "criteria": {  
             "SeriesDescription": "*'${search_string}'*","ImageType": [
                "ORIGINAL","PRIMARY","M","FFE","M","FFE"]}}'

            sub_bids_[$bs]=$(echo ${sub_bids_SWI} | python -m json.tool)

        fi

    fi

    if [[ ${identifier} == "ASL" ]]; then 
        
        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_SWI='{"dataType": "perf", "modalityLabel": "asl", 
                "criteria": {  
                    "SeriesDescription": "*'${search_string}'*",
                    "ImageType": [
                        "ORIGINAL","PRIMARY","PERFUSION","NONE","REAL"
                    ]}}'

            sub_bids_[$bs]=$(echo ${sub_bids_SWI} | python -m json.tool)

        fi

    fi

    if [[ ${identifier} == "FLAIR" ]]; then 
        
        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_FL='{"dataType": "anat", "modalityLabel": "FLAIR", "criteria": { 
             "SeriesDescription": "*'${search_string}'*"}}'

            sub_bids_[$bs]=$(echo ${sub_bids_FL} | python -m json.tool)

        fi     

    fi

    if [[ ${identifier} == "fmap" ]]; then 
        
        kul_find_relevant_dicom_file
        
        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"
            
            sub_bids_fm1='
            {
                "dataType": "fmap",
                "modalityLabel": "magnitude",
                "criteria": 
                    {
                    "SeriesDescription": "*'${search_string}'*",
                    "EchoNumber": 1
                    }
            }'
            
            #echo $sub_bids_fm1
            sub_bids_fm1a=$(echo ${sub_bids_fm1} | python -m json.tool)
            
            sub_bids_fm2='
            {
                "dataType": "fmap",
                "modalityLabel": "fieldmap",
                "criteria": 
                    {
                    "SeriesDescription": "'*${search_string}'*",
                    "EchoNumber": 2
                    },
                "sidecarChanges":
                {"Units": "Hz","IntendedFor": "##REPLACE_ME_INTENDED_FOR##"}
            }'
            
            #echo $sub_bids_fm2
            sub_bids_fm2a=$(echo ${sub_bids_fm2} | python -m json.tool)

            sub_bids_[$bs]="${sub_bids_fm1a},${sub_bids_fm2a}"

            echo $sub_bids_[$bs]

            fmap_task=$task
            intended_tasks_array=($task)

        fi     

    fi

    if [[ ${identifier} == "func" ]]; then 

        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            # remove any whitespaces
            task_nospace="$(echo -e "${task}" | tr -d '[:space:]')"

            sub_bids_fu1='{"dataType": "func","modalityLabel": 
            "bold","criteria": {"SeriesDescription": "*'${search_string}'*"},
            "customLabels": "task-'${task_nospace}''

            # add an acq_label if any
            if [ "$acq_label" = "" ];then
                sub_bids_fu1b='"',
            else
                sub_bids_fu1b='_acq-'${acq_label}'",'
            fi

            sub_bids_fu1c='"sidecarChanges": {"KUL_dcm2bids": "yes","TaskName": "'${task}'"'

            # for siemens (& ge?) ess/trt is not necessary as sidecarChanges
            # also not for philips, when it cannot be calculated
            #echo "ess_trt_provided_by_vendor: $ees_trt_provided_by_vendor"
            #echo "ees_sec: $ees_sec"
            if [ $ees_trt_provided_by_vendor -eq 1 ]  || [ $ees_sec = "empty" ]; then

                if [ $ees_trt_provided_by_vendor -eq 1 ]; then
                    kul_e2cl "   It's a SIEMENS, ees/trt are in the dicom-header" $log
                else
                    kul_e2cl "   It's NOT original dicom data (anonymised?): ees/trt could not be calculated" $log
                fi
            
                sub_bids_fu2=""

            else

                kul_e2cl "   It's a PHILIPS, ees/trt are were calculated" $log

                sub_bids_fu2=',"EffectiveEchoSpacing": '${ees_sec}',"TotalReadoutTime":
                '${trt_sec}',"MultibandAccelerationFactor": '${mb}',"PhaseEncodingDirection": "'${pe_dir}'"'
                
                    
            fi        
                    
            # for siemens (& ge?) slicetiming is not necessary as sidecarChanges
            # also not for philips, when it cannot be calculated
            #echo "slicetime_provided_by_vendor: $slicetime_provided_by_vendor"
            #echo "slice_time: $slice_time"
            if [ $slicetime_provided_by_vendor -eq 1 ]  || [ "$slice_time" = "empty" ]; then
                
                if [ $slicetime_provided_by_vendor -eq 1 ]; then
                    kul_e2cl "   It's a SIEMENS, slicetiming is in the dicom-header" $log
                else
                    kul_e2cl "   It's NOT original dicom data (anonymised?): slicetiming could not be calculated" $log
                fi

                sub_bids_fu3='}}'

            else

                kul_e2cl "   It's a PHILIPS, slicetiming was calculated" $log

                sub_bids_fu3=',"SliceTiming": '${slice_time}'}}'
            

            fi

            sub_bids_[$bs]=$(echo ${sub_bids_fu1}${sub_bids_fu1b}${sub_bids_fu1c}${sub_bids_fu2}${sub_bids_fu3} | python -m json.tool)

        fi

    fi

    if [[ ${identifier} == "sbref" ]]; then 

        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            # Take care of fact that 1 sbref could be used for multiple funcs
            # we convert the first here, and copy later (see below)
            sbref_tasks=($task)
            n_sbref_tasks=${#sbref_tasks[@]}
            sbref_task1=${sbref_tasks[0]}
            #echo $sbref_task1

            sub_bids_sb1='{"dataType": "func","modalityLabel": 
            "sbref","criteria": {"SeriesDescription": "*'${search_string}'*"}, 
            "customLabels": "task-'${sbref_task1}''
            
            if [ "$acq_label" = "" ];then
                sub_bids_sb1b='"',
            else
                sub_bids_sb1b='_acq-'${acq_label}'",'
            fi
            
            sub_bids_sb1c='"sidecarChanges": {"KUL_dcm2bids": "yes","TaskName": "'${sbref_task1}'"'

            # for siemens (& ge?) ess/trt is not necessary as sidecarChanges
            # also not for philips, when it cannot be calculated
            #echo "ess_trt_provided_by_vendor: $ees_trt_provided_by_vendor"
            #echo "ees_sec: $ees_sec"
            if [ $ees_trt_provided_by_vendor -eq 1 ]  || [ $ees_sec = "empty" ]; then

                if [ $ees_trt_provided_by_vendor -eq 1 ]; then
                    kul_e2cl "   It's a SIEMENS, ees/trt are in the dicom-header" $log
                else
                    kul_e2cl "   It's NOT original dicom data (anonymised?): ees/trt could not be calculated" $log
                fi
            
                sub_bids_sb2=""

            else

                kul_e2cl "   It's a PHILIPS, ees/trt are were calculated" $log

                sub_bids_sb2=',"EffectiveEchoSpacing": '${ees_sec}',"TotalReadoutTime":
                '${trt_sec}',"MultibandAccelerationFactor": '${mb}',"PhaseEncodingDirection": "'${pe_dir}'"'
                
                    
            fi        
                    
            # for siemens (& ge?) slicetiming is not necessary as sidecarChanges
            # also not for philips, when it cannot be calculated
            #echo "slicetime_provided_by_vendor: $slicetime_provided_by_vendor"
            #echo "slice_time: $slice_time"
            if [ $slicetime_provided_by_vendor -eq 1 ]  || [ "$slice_time" = "empty" ]; then
                
                if [ $slicetime_provided_by_vendor -eq 1 ]; then
                    kul_e2cl "   It's a SIEMENS, slicetiming is in the dicom-header" $log
                else
                    kul_e2cl "   It's NOT original dicom data (anonymised?): slicetiming could not be calculated" $log
                fi

                sub_bids_sb3='}}'

            else

                kul_e2cl "   It's a PHILIPS, slicetiming was calculated" $log

                sub_bids_sb3=',"SliceTiming": '${slice_time}'}}'
            

            fi

            sub_bids_[$bs]=$(echo ${sub_bids_sb1}${sub_bids_sb1b}${sub_bids_sb1c}${sub_bids_sb2}${sub_bids_sb3} | python -m json.tool)

        fi

    fi

    if [[ ${identifier} == "dwi" ]]; then 

        kul_find_relevant_dicom_file

        if [ $seq_found -eq 1 ]; then

            # read the relevant dicom tags
            kul_dcmtags "${seq_file}"

            sub_bids_dw1='{"dataType": "dwi","modalityLabel": "dwi",
            "criteria": {"SeriesDescription": "*'${search_string}'*"},'


            if [ "$acq_label" = "" ];then
                sub_bids_dw1b=""
            else
                sub_bids_dw1b='"customLabels": "acq-'${acq_label}'",'
            fi


            sub_bids_dw1c='"sidecarChanges": {"KUL_dcm2bids": "yes"'
            

            # for siemens (& ge?) ess/trt is not necessary as sidecarChanges
            # also not for philips, when it cannot be calculated
            #echo "ees_trt_provided_by_vendor: $ees_trt_provided_by_vendor"
            #echo "ees_sec: $ees_sec"
            if [ $ees_trt_provided_by_vendor -eq 1 ]  || [ $ees_sec = "empty" ]; then

                if [ $ees_trt_provided_by_vendor -eq 1 ]; then
                    kul_e2cl "   It's a SIEMENS, ees/trt are in the dicom-header" $log
                else
                    kul_e2cl "   It's NOT original dicom data (anonymised?): ees/trt could not be calculated" $log
                fi
            
                sub_bids_dw2=""

            else

                kul_e2cl "   It's a PHILIPS, ees/trt were calculated" $log

                sub_bids_dw2=',"EffectiveEchoSpacing": '${ees_sec}', "TotalReadoutTime": '${trt_sec}', 
                "MultibandAccelerationFactor": '${mb}', "PhaseEncodingDirection": "'${pe_dir}'"'
                    
            fi
                    
            # for siemens (& ge?) slicetiming is not necessary as sidecarChanges
            # also not for philips, when it cannot be calculated
            #echo "slicetime_provided_by_vendor: $slicetime_provided_by_vendor"
            #echo "slice_time: $slice_time"
            if [ $slicetime_provided_by_vendor -eq 1 ]  || [ "$slice_time" = "empty" ]; then

                if [ $slicetime_provided_by_vendor -eq 1 ]; then
                    kul_e2cl "   It's a SIEMENS, slicetiming is in the dicom-header" $log
                else
                    kul_e2cl "   It's NOT original dicom data (anonymised?): slicetiming could not be calculated" $log
                fi
            
                sub_bids_dw3='}}'    

            else

                kul_e2cl "   It's a PHILIPS, slicetiming was calculated" $log

                sub_bids_dw3=', "SliceTiming": '${slice_time}'}}'
            

            fi

            sub_bids_[$bs]=$(echo ${sub_bids_dw1}${sub_bids_dw1b}${sub_bids_dw1c}${sub_bids_dw2}${sub_bids_dw3} | python -m json.tool)

        fi

    fi

 fi 

done < $conf

# make the full bids_conf string and write it to file

bids_conf=""

for bf in ${!sub_bids_[@]}; do

    #echo "now generating dcm2bids config files"

    if [[ ! ${bids_conf} ]]; then

        #echo "first field"

        bids_conf="${sub_bids_[$bf]}"

    else

        #echo "Series ${bf}"

        bids_conf="${bids_conf}, ${sub_bids_[$bf]}"

    fi

done

# WRITE THE .JSON FILE USED BY dcm2bids
# Note: we also set dcm2niix options here and set pydeface
json_anon=""
if [ $anon -eq 1 ]; then
    $json_anon='"defaceTpl": "pydeface --outfile {dstFile} {srcFile}",'
fi
bids_conf_str="{${json_anon}
 \"dcm2niixOptions\": \"-b y -ba y -z y -i y -f '%3s_%f_%p_%t'\",
 \"descriptions\":[ ${bids_conf} ] }"
#echo ${bids_conf_str}
echo ${bids_conf_str} | python -m json.tool > ${bids_config_json_file}


# MAIN HERE - WE RUN dcm2bids - HERE
# invoke dcm2bids
kul_e2cl "  Calling dcm2bids... (for the output see $dcm2niix_log_file)" $log
if [ ! -d BIDS ];then
    mkdir -p BIDS
    cd BIDS
    dcm2bids_scaffold
    echo "tmp_dcm2bids/*" > .bidsignore
    echo "*/anat/*SWI*" >> .bidsignore
    echo "*/anat/*MTI*" >> .bidsignore
    echo "*/perf/*asl*" >> .bidsignore
    cd ..
fi

if [ $sess_flag -eq 1 ]; then
    dcm2bids_session=" -s ${sess} "
else
    dcm2bids_session=""
fi

dcm2bids  -d "${tmp}" -p $subj $dcm2bids_session -c $bids_config_json_file \
    -o $bids_output -l DEBUG --clobber > $dcm2niix_log_file


# Multi Echo func needs extra work. dcm2bids does not convert these correctly. "run" needs to be "echo"

if [[ ${sess} = "" ]] ; then 
    ses_long=""
else
    ses_long="/ses-${sess}"  
fi

me_file=($(grep EchoNumber ${bids_output}/sub-${subj}${ses_long}/func/*.json 2> /dev/null | awk -F ':' '{print $1}'))
me_echo=($(grep EchoNumber ${bids_output}/sub-${subj}${ses_long}/func/*.json 2> /dev/null | awk -F ':' '{print $3}' | cut -c 2 )) 

if [[ ${me_file} = "" ]] ; then 

    echo " No Multiecho fMRI data found "

else

    echo " Multiecho fMRI data found "

    n_multi_echo=${#me_echo[@]}

    for echo_number in $(seq 0 $(($n_multi_echo-1)) ) ; do 

        me_file_before_run=$(echo ${me_file[$echo_number]} | awk -F '_run-' '{print $1}')
        me_file_after_run=$(echo ${me_file[$echo_number]} | awk -F '_run-' '{print $2}')
        me_file_after_run=${me_file_after_run:2}
        cmd_json="mv ${me_file[$echo_number]} ${me_file_before_run}_echo-${me_echo[$echo_number]}${me_file_after_run}"
        cmd_nii=$(echo $cmd_json | perl -p -e 's/json/nii.gz/g')

        eval $cmd_json
        eval $cmd_nii

    done


fi


# Update the Intended For of the fmaps
# Here we define the Intended For according to the BIDS specs
#  It tells fmriprep how to use the fmap, notably for which func(s)
#  fmap_task variable (1 strings or space-separated string) is used 

if [[ ${fmap_task} ]] ; then 

    for intended_task in "${intended_tasks_array[@]}"; do
        #echo $intended_task
        #echo $cwd
        search_runs_of_task=($(find ${cwd}/${bids_output}/sub-${subj}/func -type f | grep task-${intended_task} | grep nii.gz))
        #echo ${search_runs_of_task[@]}

        n_runs=${#search_runs_of_task[@]}
        echo "  we found $n_runs of task $intended_task"

        for run_func in ${search_runs_of_task[@]}; do
                    
            if [[ $sess = "" ]]; then

                intended_for_string="func\/$(basename $run_func)"
                    
            else

                intended_for_string="ses-${sess}\/func\/$(basename $run_func)"
                    
            fi
                    
            full_intended_for_string="$full_intended_for_string, \"${intended_for_string}\""
                
        done
            
    done

    full_intended_for_string="[ ${full_intended_for_string:1} ]"      

    # Now we replace it in the json file
    echo "  NOTE: we set the following string as Intended_For in the fieldmap: ${full_intended_for_string}"

    if [[ $sess = "" ]]; then

        perl  -pi -e "s/\"##REPLACE_ME_INTENDED_FOR##\"/${full_intended_for_string}/g" ${bids_output}/sub-${subj}/fmap/sub-${subj}_fieldmap.json
                    
    else

        perl  -pi -e "s/\"##REPLACE_ME_INTENDED_FOR##\"/${full_intended_for_string}/g" ${bids_output}/sub-${subj}/ses-${sess}/fmap/sub-${subj}_fieldmap.json
                    
    fi

else 

    echo " No fmap tasks given "

fi

# Take care of fact that 1 sbref could be used for multiple funcs
# we convert the first above, now copy for each task
i=$((n_bref_tasks-1))
if [ $n_sbref_tasks -gt 1 ];then
    for t in ${sbref_tasks[@]:0-$i}; do
        cp  ${bids_output}/sub-${subj}/func/sub-${subj}_task-${sbref_tasks[0]}_sbref.json \
            ${bids_output}/sub-${subj}/func/sub-${subj}_task-${t}_sbref.json
        cp  ${bids_output}/sub-${subj}/func/sub-${subj}_task-${sbref_tasks[0]}_sbref.nii.gz \
            ${bids_output}/sub-${subj}/func/sub-${subj}_task-${t}_sbref.nii.gz
    done
fi


# copying task based events.tsv to BIDS directory
if [ $events_flag -eq 1 ]; then
    test_events_exist=$(ls -l *conf*/task-*_events.tsv | grep "No such file")
    #echo $test_events_exist
    if [ "$test_events_exist" = "" ]; then
        kul_e2cl "Copying task based events.tsv to BIDS directory" $log
        cp *conf*/task-*_events.tsv $bids_output
    fi  
fi

# clean up
cleanup="rm -fr ${tmp}"
echo ${cleanup}
eval ${cleanup}

# Fix BIDS validation
echo "This BIDS was made using KUL_NeuroImagingTools" >> ${bids_output}/README
sed -i.bck 's/"Funding": ""/"Funding": [""]/' ${bids_output}/dataset_description.json
rm ${bids_output}/dataset_description.json.bck

# Run BIDS validation
docker run -ti --rm -v ${cwd}/${bids_output}:/data:ro bids/validator /data


kul_e2cl "Finished $script" $log