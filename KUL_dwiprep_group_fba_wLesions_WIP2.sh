#!/bin/bash
set -x
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3 
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# @ Ahmed Radwan - KUL - ahmed.radwan@kuleuven.be
#
# v1.0 - dd 19/11/2021 - beta version
v="v1.1 - dd 12/07/2023"

# Changes made by AR:
# removed dependence on FWT, FS and MSBP for subjects!


# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)
ncpu_foreach=4
# suffix="_reg2T1w"
#suffix=""

#select_shells="0 700 1000 2000"
select_shells=""

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs a group fixel based analysis

Usage:

  `basename $0` -g group_name <OPT_ARGS>

Example:

  `basename $0` -g group_first_32 -n 6 -t "pat01 pat02 pat03 pat04 pat05 con01 con02 con03 con04 con05" 

Required arguments:

     -g: group_name

Optional arguments:
     
     -a:  algorithm: single shell single tissue (ssst), single shell multi tissue (ssmt), multi-shell multi-tissue (default=msmt)
     -t:  subjects used for population_template (useful if you have more than 30 to 40 subjects, otherwise the template building takes very long)
     -n:  number of cpu for parallelisation
     -s:  suffix of input dwis. Options are: 1- to leave it blank (to use native dMRI space processed dMRIs), 2- _reg2T1w (to use the processed dMRIs in native T1 space)
     -w:  specify the workflow. Options are: 1- Whole brain fixel-based (WB-FBA), 2- Bundle specific fixel-based (BS_FBA), 3- Whole brain peaks based (PBA), 4- Bundle specific peaks based (BS_PBA)
     *** NB: If opting for a bundle specific workflow (2 or 4) please make sure you have already generated FreeSurfer and MultiScale Brain Parcellator (v.1.1.1) outputs and included them in the BIDS/derivatives folder
     -L: lesion masks flag. Options are: 1- if KUL_dwiprep.sh was run with -L flag and dwi_mask_minLesion.nii.gz files were generated, 
                                         2- if lesion masks are included in a separate folder (naming should follow BIDS convention), lesion masks must be in processed dMRI space
                                         The lesion masks folder (named Lesion_mask) must be placed in the same directory as BIDS and dwiprep and contain files named as "sub-*_dwi_lesion_mask.nii.gz"
                                         and should be binary i.e. encoding the lesion as 1 and everything else as 0
     -v:  show output from mrtrix commands
     -S:  Create a symmetrical population template by using also flipped versions of each input fod (for template creation only)


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
ncpu=6 # default if option -n is not given
silent=1 # default if option -v is not given
algo=msmt

# Set required options
g_flag=0
t_flag=0
w_flag=0
a_flag=0
Symm_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "n:g:t:a:L:w:vS" OPT; do

        case $OPT in
        n) #ncpu
            ncpu=$OPTARG
        ;;
        g) #group_name
            group_name="$OPTARG"
            g_flag=1
        ;;
        t) #templatesubjects
            templatesubjects="$OPTARG"
            t_flag=1
        ;;
        a) #algorithm
            a_flag=1
            algo="$OPTARG"
        ;;
        s) #suffix 
            suffix="$OPTARG"
        ;;
        L) #Lesions 
            L_flag="$OPTARG"
        ;;
        w) # workflow
            w_flag="$OPTARG"
        ;;
        v) #verbose
            silent=0
        ;;
        S) #symmetrical along X
            Symm_flag=1
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
if [ $g_flag -eq 0 ] ; then 
    echo 
    echo "Option -g is required: give the group name for the analysis." >&2
    echo
    exit 2 
fi

if [ $a_flag -eq 0 ] ; then 
    echo 
    echo "Option -a is not specified using the multi-shell multitissue (msmt) CSD by defailt." >&2
    echo
else
    echo 
    if [[ $algo == "ssmt" ]]; then
        echo "Option -a is specified using the single-shell two-tissue ($algo) CSD." >&2
    elif [[ $algo == "msmt" ]]; then
        echo "Option -a is specified using the multi-shell multi-tissue ($algo) CSD." >&2
    elif [[ $algo == "ssst" ]]; then
        echo "Option -a is specified using the single-shell single-tissue ($algo) CSD." >&2       
    fi
    echo
fi 


# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then 

    export MRTRIX_QUIET=1

fi

# check whether this analysis will be symmetrical along x axis or not
if [ $Symm_flag = 0 ]; then
    echo "Option -S was not selected, we will not impose symmtery across x-axis in population template." >&2       
else 
    echo "Option -S was selected, we will do impose symmtery across x-axis in population template." >&2       
fi

# check suffix and use blank if not set
# Add a third condition here if needed
if [ -z $suffix ]; then
    suffix=""
    echo "No suffix is specified for dMRI data, a blank suffix will be used"
    echo "We will use the processed dMRIs in native diffusion space"
else
    echo "dMRI suffix is specified as ${suffix}"
    echo "We will use the processed dMRIs in native T1 space"
fi


# check lesion masks optional flag
if [[ ${L_flag} == 0 ]]; then
    echo "No lesion masks specified"
    echo "We will do a typical FBA not accounting for lesions"
elif [[ ${L_flag} == 1 ]]; then
    echo "Lesion masks will be used"
    echo "You specified -L 1, meaning that we expect KUL_dwiprep.sh to have been run with lesion masks, and the dwi_mask_minLesion.nii.gz files have been generated"
elif [[ ${L_flag} == 2 ]]; then
    echo "Lesion masks will be used"
    echo "You specified -L 2, meaning that we expect a folder called Lesion_masks to be present in ${cwd} with BIDS naming for each file"
fi

# check workflow optional flag
if [[ ${w_flag} == 0 ]]; then
    echo "Workflow flag not set, running default whole brain fba"
elif [[ ${w_flag} == 1 ]]; then
    echo "Workflow flag set to 1, running default whole brain fba"
elif [[ ${w_flag} == 2 ]]; then
    echo "Workflow flag set to 2, running bundle specific fba"
elif [[ ${w_flag} == 3 ]]; then
    echo "Workflow flag set to 3, running whole brain pba"
elif [[ ${w_flag} == 4 ]]; then
    echo "Workflow flag set to 4, running bundle specific pba"
fi

# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Some parallelisation
FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=log/log_${d}.txt

# --- MAIN ----------------

# make dirs
mkdir -p dwiprep/${group_name}/fba/subjects
if [ "$algo" = "ssst" ]; then 
    mkdir -p dwiprep/${group_name}/fba/dwiintensitynorm/dwi_input
    mkdir -p dwiprep/${group_name}/fba/dwiintensitynorm/mask_input
fi

cd dwiprep/${group_name}/fba

# find the preproced mifs
if [ ! -f data_prep.done ]; then 
    echo "   Preparing data in dwiprep/${group_name}/fba/"
    # finding all subjects in the dwiprep folder
    search_subjects=($(find ${cwd}/dwiprep/sub-* -type f | grep dwi_preproced${suffix}.mif | sort ))
    # num_sessions=${#search_subjects[@]}
    for i in ${search_subjects[@]}
    do
        sub=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
        ses=$(echo $i | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
        if [[ ! -z ${ses} ]]; then
            s=${sub}_${ses}
        else
            s=${sub}
        fi
        mkdir -p ${cwd}/dwiprep/${group_name}/fba/subjects/${s}
        ln -sfn $i ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}.mif
        if [ "$algo" = "ssst" ]; then 
            ln -sfn $i ${cwd}/dwiprep/${group_name}/fba/dwiintensitynorm/dwi_input/${s}_dwi_preproced${suffix}.mif
        fi

    done 
    
    # # find the preproced masks
    # # need to make sure this is correct in case KUL_dwiprep_anat.sh is used in advance
    # search_subjects=($(find ${cwd}/dwiprep -type f | grep dwi_mask${suffix}.nii.gz | sort ))
    # num_subjects=${#search_subjects[@]}

    # for i in ${search_subjects[@]}
    # do

    #     sub=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
    #     ses=$(echo $i | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
    #     # s=${sub}_${ses}
    #     if [[ ! -z ${ses} ]]; then
    #         s=${sub}_${ses}
    #     else
    #         s=${sub}
    #     fi

    #     if [ "$algo" = "ssst" ]; then 
        
    #         mrconvert $i ${cwd}/dwiprep/${group_name}/fba/dwiintensitynorm/mask_input/${s}_dwi_preproced${suffix}.mif -force
        
    #     fi

    # done

    # find the preproced masks min lesion if -L 1 is given
    # ATM this is only tested for native DWI space not for anat or MNI space
    # We will add support for KUL_dwiprep_anat and KUL_dwiprep_MNI later on

    if [[ ! -f mask.done ]]; then
        search_subjects=($(find ${cwd}/dwiprep -type f | grep dwi_mask${suffix}.nii.gz | sort ))
        if [[ "${L_flag}" == 1 ]]; then
            search_subjects_Lmask=($(find ${cwd}/dwiprep -type f | grep dwi_mask${suffix}_minLesion.nii.gz | sort ))
        elif [[ "${L_flag}" == 2 ]]; then
            search_subjects_Lmask=($(ls -f ${cwd}/Lesion_masks/sub-*_dwi_lesion_mask.nii.gz | sort ))
        fi
        num_subjects=${#search_subjects[@]}

        # exit 2

        for i in ${!search_subjects[@]}
        do
            sub=$(echo ${search_subjects[$i]} | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
            ses=$(echo ${search_subjects[$i]} | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
            
            if [[ ! -z ${ses} ]]; then
                s=${sub}_${ses}
            else
                s=${sub}
            fi

            # regular brain masks
            mkor=$(dirname ${search_subjects[$i]} | rev | cut -d "/" -f1 | rev)
            # brain masks min lesions
            for m in ${!search_subjects_Lmask[@]}; do 
                if [[ "${L_flag}" == 1 ]]; then
                    mkLe=$(dirname ${search_subjects_Lmask[$m]} | rev | cut -d "/" -f1 | rev)
                    bml=${search_subjects_Lmask[$m]}
                elif [[ "${L_flag}" == 2 ]]; then
                    mkLe=$(basename ${search_subjects_Lmask[$m]} | cut -d "_" -f1)
                    bml=${search_subjects_Lmask[$m]}
                fi
                # find out whether we'll use whole brain mask or brain mask min lesion
                if [[ ! "${mkor}" == "${mkLe}" ]]; then
                    echo " This subject ${search_subjects[$i]} has no brain_mask_minLesion.nii.gz file"
                    echo " We will use the regular brain mask for ${search_subjects[$i]}"
                    mrconvert ${search_subjects[$i]} ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif -force
                elif [[ "${mkor}" == "${mkLe}" ]]; then
                    echo " This subject ${search_subjects[$i]} has a lesion mask"
                    echo " We will use the lesion excluded brain mask for ${search_subjects[$i]} for registration"
                    mrconvert ${search_subjects[$i]} ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif -force
                    # also handle the brain mask minL work - replacing all zeros with nans
                    if [[ "${L_flag}" == 1 ]]; then
                        mrconvert ${bml} ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL.mif -force
                        mrconvert ${bml} - | mrthreshold - -abs 0.0 -comparison gt -nan ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif -force
                    elif [[ "${L_flag}" == 2 ]]; then
                        mrcalc ${bml} 0 -gt -1 -mult ${search_subjects[$i]} -add 0 -gt ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL.mif -force
                        mrcalc ${bml} 0 -gt -1 -mult ${search_subjects[$i]} -add 0 -gt - | mrthreshold - -abs 0.0 -comparison gt -nan ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif -force
                    fi
                fi
            done
            # All masking stuff is handled here - replacing all zeros with nans
            if [ "$algo" = "ssst" ]; then 
                mrconvert ${search_subjects[$i]} ${cwd}/dwiprep/${group_name}/fba/dwiintensitynorm/mask_input/${s}_dwi_preproced${suffix}.mif -force
            fi
        done
        # if [ $? -eq 0 ]; then
        #     echo "done" > mask.done
        # fi
    else
        echo "   Computing new brain mask images already done"
    fi

    # exit 2

    if [ "$algo" = "ssmt" ] || [ "$algo" = "msmt" ]; then 
        # find the response functions - CSF
        search_subjects=($(find ${cwd}/dwiprep -type f | grep dhollander_csf_response.txt | sort ))
        num_subjects=${#search_subjects[@]}
        for i in ${search_subjects[@]}
        do
            #s=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
            sub=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
            ses=$(echo $i | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
            # s=${sub}_${ses}
            if [[ ! -z ${ses} ]]; then
                s=${sub}_${ses}
            else
                s=${sub}
            fi
            ln -sfn $i ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dhollander_csf_response.txt
        done

        # find the response functions - GM
        search_subjects=($(find ${cwd}/dwiprep -type f | grep dhollander_gm_response.txt | sort ))
        num_subjects=${#search_subjects[@]}
        for i in ${search_subjects[@]}
        do
            #s=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
            sub=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
            ses=$(echo $i | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
            # s=${sub}_${ses}
            if [[ ! -z ${ses} ]]; then
                s=${sub}_${ses}
            else
                s=${sub}
            fi
            ln -sfn $i ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dhollander_gm_response.txt
        done

        # find the response functions - WM
        search_subjects=($(find ${cwd}/dwiprep -type f | grep dhollander_wm_response.txt | sort ))
        num_subjects=${#search_subjects[@]}
        for i in ${search_subjects[@]}
        do
            #s=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
            sub=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
            ses=$(echo $i | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
            # s=${sub}_${ses}
            if [[ ! -z ${ses} ]]; then
                s=${sub}_${ses}
            else
                s=${sub}
            fi
            ln -sfn $i ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dhollander_wm_response.txt
        done
    fi
    echo "done" > data_prep.done
else
    echo "   Preparing data in dwiprep/${group_name}/fba/ already done"
fi

# Option to select certain shells from the data
if [ "$select_shells" = "" ]; then 
    echo "No shell selection, just continue"
else
    echo "Shells $select_shells will now be used in further analysis"
    echo "NOT YET IMPLEMENTED!!!!"
fi

# STEP 1 - Intensity Normalisation (only for ST data)
if [ "$algo" = "ssst" ]; then 
    if [ ! -f dwiintensitynorm/fa_template_wm_mask.mif ]; then
        echo "   Doing Intensity Normalisation"
        dwinromalise group dwiintensitynorm/dwi_input/ dwiintensitynorm/mask_input/ \
        dwiintensitynorm/dwi_output/ dwiintensitynorm/fa_template.mif \
        dwiintensitynorm/fa_template_wm_mask.mif -nthreads $ncpu -force
        mrinfo dwiintensitynorm/dwi_output/* -property dwi_norm_scale_factor > CHECK_dwi_norm_scale_factor.txt
    else
        echo "   Intensity Normalisation already done"
    fi

    # Adding a subject
    # dwi2tensor new_subject/dwi_denoised_unringed_preproc_unbiased.mif -mask new_subject/dwi_temp_mask.mif - | tensor2metric - -fa - | mrregister -force \
    # ../dwiintensitynorm/fa_template.mif - -mask2 new_subject/dwi_temp_mask.mif -nl_scale 0.5,0.75,1.0 -nl_niter 5,5,15 -nl_warp - /tmp/dummy_file.mif | \
    # mrtransform ../dwiintensitynorm/fa_template_wm_mask.mif -template new_subject/dwi_denoised_unringed_preproc_unbiased.mif -warp - - | dwinormalise \ 
    # new_subject/dwi_denoised_unringed_preproc_unbiased.mif - ../dwiintensitynorm/dwi_output/new_subject.mif

fi

# Link back the normalised data
cd ${cwd}/dwiprep/${group_name}/fba/subjects

if [ "$algo" = "ssst" ]; then 
    for_each -info -force * : ln -sfn ${cwd}/dwiprep/${group_name}/fba/dwiintensitynorm/dwi_output/PRE_dwi_preproced${suffix}.mif \
    ${cwd}/dwiprep/${group_name}/fba/subjects/IN/dwi_preproced${suffix}_normalised.mif
fi

# STEP 2 - Computing an (average) white matter response function
if [ ! -f ../average_response.done ]; then 
    echo "   Computing an (average) white matter response function"
    if [ "$algo" = "ssst" ]; then 
        for_each -info -force * : dwi2response tournier IN/dwi_preproced${suffix}_normalised.mif \
        IN/response.txt -nthreads $ncpu -force
        responsemean */response.txt ../group_average_response.txt
    else
        responsemean */dhollander_wm_response.txt ../group_average_response_wm.txt
        responsemean */dhollander_gm_response.txt ../group_average_response_gm.txt
        responsemean */dhollander_csf_response.txt ../group_average_response_csf.txt
    fi

    if [ $? -eq 0 ]; then
        echo "done" > ../average_response.done
    fi
else
    echo "   Computing of an (average) white matter response function already done"
fi

# Use same masks generated by KUL_dwiprep - already handled above
# foreach * : dwi2mask IN/dwi_denoised_unringed_preproc_unbiased_normalised_upsampled.mif IN/dwi_mask_upsampled.mif

# STEP 3 - Fibre Orientation Distribution estimation (spherical deconvolution)
# see https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/st_fibre_density_cross-section.html
# Note that dwi2fod csd can be used, however here we use dwi2fod msmt_csd (even with single shell data) to benefit from the hard non-negativity constraint, which has been observed to lead to more robust outcomes

if [ ! -f ../fod_estimation.done ]; then
    echo "   Performing FOD estimation"
    if [ "$algo" = "ssst" ]; then 
        for_each -force -info -nthreads ${ncpu_foreach} * : dwiextract IN/dwi_preproced${suffix}_normalised.mif - \
        \| dwi2fod msmt_csd - ../group_average_response.txt IN/wmfod.mif \
        -mask IN/dwi_preproced${suffix}_mask.mif -info -force -nthreads 2
    elif [ "$algo" = "ssmt" ]; then 
        for_each -force -info -nthreads ${ncpu_foreach} * : dwi2fod msmt_csd IN/dwi_preproced${suffix}.mif \
        ../group_average_response_wm.txt IN/wmfod_nogm.mif \
        ../group_average_response_csf.txt IN/csf_nogm.mif \
        -mask IN/dwi_preproced${suffix}_mask.mif -info -force -nthreads 2
    elif [ "$algo" = "msmt" ]; then 
        for_each -force -info -nthreads ${ncpu_foreach} * : dwi2fod msmt_csd IN/dwi_preproced${suffix}.mif \
        ../group_average_response_wm.txt IN/wmfod.mif \
        ../group_average_response_gm.txt IN/gm.mif \
        ../group_average_response_csf.txt IN/csf.mif \
        -mask IN/dwi_preproced${suffix}_mask.mif -info -force -nthreads 2
    fi
    if [ $? -eq 0 ]; then
        echo "done" > ../fod_estimation.done
    fi
fi

# STEP 3b - for multi-tissue only - Joint bias field correction and intensity normalisation
if [ "$algo" = "msmt" ]; then 
    if [ ! -f ../mtnormalise.done ]; then
        for_each -info -force -nthreads ${ncpu_foreach} * : mtnormalise IN/wmfod.mif IN/wmfod_norm.mif \
        IN/gm.mif IN/gm_norm.mif IN/csf.mif IN/csf_norm.mif \
        -mask IN/dwi_preproced${suffix}_mask.mif
        if [ $? -eq 0 ]; then
            echo "done" > ../mtnormalise.done
        fi
    fi
elif [ "$algo" = "ssmt" ]; then 
    if [ ! -f ../mtnormalise.done ]; then
        for_each -info -force -nthreads ${ncpu_foreach} * : mtnormalise IN/wmfod_nogm.mif IN/wmfod_norm.mif \
        IN/csf_nogm.mif IN/csf_norm.mif \
        -mask IN/dwi_preproced${suffix}_mask.mif
        if [ $? -eq 0 ]; then
            echo "done" > ../mtnormalise.done
        fi
    fi
fi

# exit 2

# adding volume fraction calculation
search_subjects=($(ls -d ${cwd}/dwiprep/sub-*/sub-*))
if [[ ! -f "../vf_calc.done" ]]; then
    for i in ${!search_subjects[@]}; do
        s=$(basename ${search_subjects[$i]} | cut -d '-' -f2)
        # echo ${s}
        if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/vf_norm.mif" ]]; then
            mrconvert -coord 3 0 ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/wmfod_norm.mif - | \
            mrcat ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/csf_norm.mif ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/gm_norm.mif - \
            ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/vf_norm.mif
        fi
    done
    if [ $? -eq 0 ]; then
        echo "done" > ../vf_calc.done
    fi
fi
# exit 2

# STEP3C - for the bundle specific workflow, we need to bring the wmfod_norm back to dwiprep
# we also need to rename the response folder and create a new one

if [ ! "$algo" = "ssst" ]; then 
    search_fods=($(ls -f ${cwd}/dwiprep/${group_name}/fba/subjects/*/wmfod_norm.mif))
else
    search_fods=($(ls -f ${cwd}/dwiprep/${group_name}/fba/subjects/*/wmfod.mif))
fi
if [[ ! ${#search_subjects[@]} == ${#search_fods[@]} ]]; then
    echo " N subjects ${#search_subjects[@]} not equal to N of wmfod_norm.mif files ${#search_fods[@]}, exiting."
    exit 2
fi

# exit 2

# STEP3D
# generating T1 images in dMRI space
mkdir -p ${cwd}/BIDS/derivatives/cmp
search_subjects=($(ls -d ${cwd}/dwiprep/sub-*/sub-*))
ncpu_split=$(echo "scale=0; $ncpu/$ncpu_foreach" | bc)
ny=0
# we parallelize this over subjects as well - use ncpu_foreach for intersubject parallel limit, and ncpu_split for intrasubject parallel limit

for r in ${!search_subjects[@]}; do
    si=$(basename ${search_subjects[$r]} | cut -d "_" -f1 | cut -d "-" -f2)
    mkdir -p ${cwd}/BIDS/derivatives/cmp/sub-${si}
    if [[ ! -f "${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain.nii.gz" ]]; then 
        hd-bet -i ${cwd}/BIDS/sub-${si}/anat/sub-${si}_T1w.nii.gz \
        -o ${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain.nii.gz -mode accurate -s 1
    fi
done

unset r si

for r in ${!search_subjects[@]}; do
    ((ny++))
    ((ny=${ny}%${ncpu_foreach}))
    si=$(basename ${search_subjects[$r]} | cut -d "_" -f1 | cut -d "-" -f2)
    # to remove dependence on KUL_FWT_make_VOIs per subject
    # first convert the FS brain to nifti
    
    # first generate warp from FS brain to UKBB template
    if [[ ! -f "${cwd}/BIDS/derivatives/cmp/sub-${si}/fa_2_UKBB_vT1w_template.nii.gz" ]]; then

        ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ncpu_split}
        antsRegistrationSyN.sh -d 3 -f /usr/local/KUL_apps/KUL_FWT/KUL_FWT_templates/T1_preunbiased.nii.gz \
        -m ${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain.nii.gz \
        -x /usr/local/KUL_apps/KUL_FWT/KUL_FWT_templates/T1_preunbiased.nii.gz,${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain_mask.nii.gz \
        -o ${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain_2_UKBB_ -t s -n ${ncpu_split} \
        && antsIntermodalityIntrasubject.sh -d 3 -i ${cwd}/dwiprep/sub-${si}/sub-${si}/qa/fa.nii.gz -r ${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain.nii.gz \
        -x ${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain_mask.nii.gz -t 3 -w ${cwd}/BIDS/derivatives/cmp/sub-${si}/sub-${si}_T1w_brain_2_UKBB_ \
        -T /usr/local/KUL_apps/KUL_FWT/KUL_FWT_templates/T1_preunbiased.nii.gz -o ${cwd}/BIDS/derivatives/cmp/sub-${si}/fa_2_UKBB_vT1w_ &
    fi

    # second generate warp from T1 to FOD


    # KUL_FWT_make_VOIs.sh -p ${si} -M ${cwd}/BIDS/derivatives/cmp/sub-${si}/anat/sub-${si}_ses-01_label-L2018_desc-scale3_atlas.nii.gz \
    # -F ${cwd}/BIDS/derivatives/freesurfer/sub-${si}/mri/aparc+aseg.mgz \
    # -c ${cwd}/FWT_config.txt -d ${cwd}/dwiprep/sub-${si}/sub-${si} -o ${cwd}/BIDS/derivatives/FWT/sub-${si} -n ${ncpu_split} &

    if [[ $ny == 0 ]]; then
        wait
    fi
done
wait

# exit 2

echo " Finished with FWT_make_VOIs runs back to FBA script now "

# STEP 4 - Generate a study-specific unbiased FOD template
# This will be different if -L was used or not
# we could also include a flipping in RL step and creation of a perfectly symmetrical template
mkdir -p ../template
mkdir -p ../template/fod_input
mkdir -p ../template/gm_input
mkdir -p ../template/csf_input
mkdir -p ../template/mask_input
if [[ ! "${L_flag}" == 0 ]]; then
    mkdir -p ../template/nan_mask_input
fi

# declare -a links
templatesubjects_a=(${templatesubjects})
echo ${templatesubjects_a[@]}

if [ ! -f ../pop_template.done ]; then
    echo "   Generating FOD template"
    # search_sessions=($(find ${cwd}/dwiprep/${group_name}/fba/subjects | grep wmfod_norm.mif | sort ))
    search_sessions=($(ls -f ${cwd}/dwiprep/${group_name}/fba/subjects/*/wmfod_norm.mif))
    for t in ${!search_sessions[@]}; do
        s=$(echo ${search_sessions[$t]} | rev | cut -d '/' -f2 | rev)
        # If the template flag was used
        if [ $t_flag -eq 1 ]; then
            # Don't link subjects not given in -t
            for hb in ${!templatesubjects_a[@]}; do
                # echo "${templatesubjects_a[$hb]}"
                if [[ "${s}" == "${templatesubjects_a[$hb]}" ]]; then
                    ln -sfn ${search_sessions[$t]} ${cwd}/dwiprep/${group_name}/fba/template/fod_input/${s}_wmfod_norm.mif
                    # always need to check for gm and csf fods in case one or both of them don't exist (ssst and ssmt)
                    if [[ -f "$(dirname ${search_sessions[$t]})/gm_norm.mif" ]]; then
                        ln -sfn $(dirname ${search_sessions[$t]})/gm_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/gm_input/${s}_gm_norm.mif
                    fi
                    if [[ -f "$(dirname ${search_sessions[$t]})/csf_norm.mif" ]]; then
                        ln -sfn $(dirname ${search_sessions[$t]})/csf_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/csf_input/${s}_csf_norm.mif
                    fi
                    ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
                    ${cwd}/dwiprep/${group_name}/fba/template/mask_input/${s}_mask.mif
                    # here there's no diff between -L 1 and -L 2, we use the existence of _mask_minL.mif as it's already present and easier to use to iterate over subjects
                    if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}_mask.mif" ]]; then
                        if [[ -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif" ]]; then 
                            ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif \
                            ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}_mask.mif
                        else
                            # if the mask_minL isn't present simply use the whole brain mask
                            ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
                            ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}_mask.mif
                        fi
                    fi
                fi
            done
        else
            # If the template flag was not used
            ln -sfn ${search_sessions[$t]} ${cwd}/dwiprep/${group_name}/fba/template/fod_input/${s}_wmfod_norm.mif
            if [[ -f "$(dirname ${search_sessions[$t]})/gm_norm.mif" ]]; then
                ln -sfn $(dirname ${search_sessions[$t]})/gm_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/gm_input/${s}_gm_norm.mif
            fi
            if [[ -f "$(dirname ${search_sessions[$t]})/csf_norm.mif" ]]; then
                ln -sfn $(dirname ${search_sessions[$t]})/csf_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/csf_input/${s}_csf_norm.mif
            fi
            ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
            ${cwd}/dwiprep/${group_name}/fba/template/mask_input/${s}_mask.mif
            # here there's no diff between -L 1 and -L 2, we use the existence of _mask_minL.mif as it's already present and easier to use to iterate over subjects
            # now the Lmasks also have nans in them
            if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}_mask.mif" ]]; then
                if [[ -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif" ]]; then 
                    ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif \
                    ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}_mask.mif
                else
                    # if the mask_minL isn't present simply use the whole brain mask - now this is nan padded as well
                    mrthreshold ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif -abs 0.0 -comparison gt -nan ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}_mask.mif
                    # the ln -sfn step below is no longer needed right?
                    # ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
                    # ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}_mask.mif
                fi
            fi
        fi

        # add flipping in x-axis if required
        if [[ $Symm_flag == 1 ]]; then
            if [ $t_flag -eq 1 ]; then
                # Don't flip subjects not given in -t
                for fb in ${!templatesubjects_a[@]}; do
                    if [[ "${s}" == "${templatesubjects_a[$fb]}" ]]; then
                        mrtransform -flip 0 ${search_sessions[$t]} ${cwd}/dwiprep/${group_name}/fba/template/fod_input/${s}f_wmfod_norm.mif -force -nthreads ${ncpu} -reorient_fod 0
                        if [[ -f "$(dirname ${search_sessions[$t]})/gm_norm.mif" ]]; then
                            mrtransform -flip 0 $(dirname ${search_sessions[$t]})/gm_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/gm_input/${s}f_gm_norm.mif -force -nthreads ${ncpu} -reorient_fod 0
                        fi
                        if [[ -f "$(dirname ${search_sessions[$t]})/csf_norm.mif" ]]; then
                            mrtransform -flip 0 $(dirname ${search_sessions[$t]})/csf_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/csf_input/${s}f_csf_norm.mif -force -nthreads ${ncpu} -reorient_fod 0
                        fi
                        mrtransform -interp nearest -flip 0 ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
                        ${cwd}/dwiprep/${group_name}/fba/template/mask_input/${s}f_mask.mif -force -nthreads ${ncpu}
                        # here there's no diff between -L 1 and -L 2, we use the existence of _mask_minL.mif as it's already present and easier to use to iterate over subjects
                        if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}f_mask.mif" ]]; then
                            if [[ -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif" ]]; then 
                                mrtransform -interp nearest -flip 0 ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif \
                                ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}f_mask.mif -force -nthreads ${ncpu}
                            else
                                # if the mask_minL isn't present simply use the whole brain mask
                                mrtransform -interp nearest -flip 0 ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
                                ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}f_mask.mif -force -nthreads ${ncpu}
                            fi
                        fi
                    fi
                done
            else
                mrtransform -flip 0 ${search_sessions[$t]} ${cwd}/dwiprep/${group_name}/fba/template/fod_input/${s}f_wmfod_norm.mif -force -nthreads ${ncpu} -reorient_fod
                if [[ -f "$(dirname ${search_sessions[$t]})/gm_norm.mif" ]]; then
                    mrtransform -flip 0 $(dirname ${search_sessions[$t]})/gm_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/gm_input/${s}f_gm_norm.mif -force -nthreads ${ncpu} -reorient_fod
                fi
                if [[ -f "$(dirname ${search_sessions[$t]})/csf_norm.mif" ]]; then
                    mrtransform -flip 0 $(dirname ${search_sessions[$t]})/csf_norm.mif ${cwd}/dwiprep/${group_name}/fba/template/csf_input/${s}f_csf_norm.mif -force -nthreads ${ncpu} -reorient_fod
                fi
                mrtransform -interp nearest -flip 0 ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
                ${cwd}/dwiprep/${group_name}/fba/template/mask_input/${s}f_mask.mif -force -nthreads ${ncpu}
                # here there's no diff between -L 1 and -L 2, we use the existence of _mask_minL.mif as it's already present and easier to use to iterate over subjects
                # now the Lmasks also have nans in them
                if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}f_mask.mif" ]]; then 
                    if [[ -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif" ]]; then
                        mrtransform -interp nearest -flip 0 ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif \
                        ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}f_mask.mif -force -nthreads ${ncpu}
                    else
                        # if the mask_minL isn't present simply use the whole brain mask - now this is nan padded as well
                        mrthreshold ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif \
                        -abs 0.0 -comparison gt -nan ${cwd}/dwiprep/${group_name}/fba/template/nan_mask_input/${s}f_mask.mif
                    fi
                fi
            fi
        fi
    done

    unset hb
    srch_subs=($(ls -d ${cwd}/BIDS/derivatives/cmp/sub-*))
    mkdir ${cwd}/dwiprep/${group_name}/fba/template/T1s
    for xi in ${!srch_subs[@]}; do
        if [ $t_flag -eq 1 ]; then
            for hb in ${!templatesubjects_a[@]}; do
                so=$(basename ${srch_subs[$xi]})
                sd=$(echo ${so} | cut -d "-" -f2)
                # bring all subject T1s to template FOD space
                # we use the T1 to keep the skulls
                # use fa as a target for now
                if [[ "${sd}" == "${templatesubjects_a[$hb]}" ]]; then
                    if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}_T1w_in_FOD_native.nii.gz" ]]; then
                        antsApplyTransforms -d 3 -i ${cwd}/BIDS/${so}/anat/${so}_T1w.nii.gz -o ${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}_T1w_in_FOD_native.nii.gz \
                        -r ${cwd}/dwiprep/${so}/${so}/qa/fa.nii.gz -t \[${srch_subs[$xi]}/fa_2_UKBB_vT1w_0GenericAffine.mat,1\] \
                        -t ${srch_subs[$xi]}/fa_2_UKBB_vT1w_1InverseWarp.nii.gz
                    fi

                    if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}f_T1w_in_FOD_native.nii.gz" ]] && [[ ${Symm_flag} == 1 ]]; then
                        mrtransform -force -flip 0 ${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}_T1w_in_FOD_native.nii.gz ${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}f_T1w_in_FOD_native.nii.gz
                    fi
                fi
            done
        else
            so=$(basename ${srch_subs[$xi]})
            sd=$(echo ${so} | cut -d "-" -f2)
            # bring all subject T1s to template FOD space
            # we use the T1 to keep the skulls
            # use fa as a target for now
            if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}_T1w_in_FOD_native.nii.gz" ]]; then
                antsApplyTransforms -d 3 -i ${cwd}/BIDS/${so}/anat/${so}_T1w.nii.gz -o ${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}_T1w_in_FOD_native.nii.gz \
                -r ${cwd}/dwiprep/${so}/${so}/qa/fa.nii.gz -t \[${srch_subs[$xi]}/fa_2_UKBB_vT1w_0GenericAffine.mat,1\] \
                -t ${srch_subs[$xi]}/fa_2_UKBB_vT1w_1InverseWarp.nii.gz
            fi

            if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}f_T1w_in_FOD_native.nii.gz" ]] && [[ ${Symm_flag} == 1 ]]; then
                mrtransform -force -flip 0 ${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}_T1w_in_FOD_native.nii.gz ${cwd}/dwiprep/${group_name}/fba/template/T1s/${so}f_T1w_in_FOD_native.nii.gz
            fi
        
        fi
    done

    # could also add the work here for making flipped copies
    # flip all fod (wm, gm, csf) and T1s along x-axis and add to template input folders
    # if [[ $Symm_flag = 1 ]]; then
    #     mrtransform -flip 0
    # fi
    
    # switched to a multi-contrast population template and added generation of template mask
    # hopefully this still works even with the nans to replace the zeros!
    ### Need to test this - AR 27/06
    # could easily also create template of FA, MD, RD, AD, etc. both with and without flipping!
    # we could even feed the dwis in here too!!
    if [[ ${L_flag} == 0 ]]; then
        mask_opt="-mask_dir ../template/mask_input"
    else
        mask_opt="-mask_dir ../template/nan_mask_input"
    fi
    if [[ ${algo} == "ssst" ]]; then
        pop_temp_ins="../template/fod_input ../template/wmfod_template.mif ../template/T1s ../template/template_T1w.mif"
    elif [[ ${algo} == "ssmt" ]]; then
        pop_temp_ins="../template/fod_input ../template/wmfod_template.mif ../template/csf_input ../template/csf_template.mif ../template/T1s ../template/template_T1w.mif"
    elif [[ ${algo} == "msmt" ]]; then
        pop_temp_ins="../template/fod_input ../template/wmfod_template.mif ../template/gm_input ../template/gm_template.mif ../template/csf_input ../template/csf_template.mif ../template/T1s ../template/template_T1w.mif"
    fi
    population_template -linear_no_pause -aggregate median -scratch ../template ${mask_opt} -nocleanup -voxel_size 1.25 -nthreads $ncpu -template_mask ../template/template_mask.mif \
        ${pop_temp_ins}
    mrconvert ../template/wmfod_template.mif -coord 3 0 -axes 0,1,2 ../template/wmfod_template_vol1.nii.gz
    mrconvert ../template/template_mask.mif ../template/template_mask.nii.gz
    if [ $? -eq 0 ]; then
        echo "done" > ../pop_template.done
    fi
else
    echo "   FOD template already generated"
fi

# exit 2
# Register all subject FOD images to the FOD template
# foreach -${ncpu_foreach} * : mrregister IN/wmfod.mif -mask1 IN/dwi_mask_upsampled.mif ../template/wmfod_template.mif -nl_warp IN/subject2template_warp.mif IN/template2subject_warp.mif
if [ ! -f ../fod_reg2template.done ]; then
    # If -L flag is used we need to do the warping differently (to use cfm, which doesn't seem to work with mrregister)
    # then we need to convert the results to an mrtrix3 compatible format
    # as described here: https://community.mrtrix.org/t/registration-using-transformations-generated-from-other-packages/2259
    if [[ ! ${L_flag} == 0 ]]; then
        search_sessions=($(ls -f ${cwd}/dwiprep/${group_name}/fba/subjects/*/wmfod_norm.mif))
        echo " Registering all subject FOD images to the FOD template"
        for l in ${!search_sessions[@]}; do
            s=$(echo ${search_sessions[$l]} | rev | cut -d '/' -f2 | rev)
            if [[ ! -f "$(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_deformation.mif" ]]; then
                mrconvert ${search_sessions[$l]} -coord 3 0 -axes 0,1,2 $(dirname ${search_sessions[$l]})/wmfod_norm_vol1.nii.gz -force
                # remove nans again to avoid interfering with registration work
                if [[ -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL.mif" ]]; then
                    mask_4_temp=${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL.mif
                else
                    mask_4_temp=${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif
                fi
                # registration   
                if [[ ! -f "$(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_warped.nii.gz" ]]; then
                    mrregister -type rigid_affine_nonlinear -nl_warp $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_deformation.mif \
                    $(dirname ${search_sessions[$l]})/template_2_wmfod_norm_vol1_deformation.mif \
                    -mask1 ${mask_4_temp} \
                    -mask2 ${cwd}/dwiprep/${group_name}/fba/template/template_mask.nii.gz \
                    -transformed $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_warped.nii.gz \
                    $(dirname ${search_sessions[$l]})/wmfod_norm_vol1.nii.gz ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz -nthreads ${ncpu} -force
                    # -nl_warp_full $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_warp_full.mif
                fi

            fi
            # forwards
            # warpinit $(dirname ${search_sessions[$l]})/wmfod_norm_vol1.nii.gz $(dirname ${search_sessions[$l]})/identity_warp[].nii.gz -force
            # for j in {0..2}; do
            #     # forward
            #     antsApplyTransforms -d 3 -i $(dirname ${search_sessions[$l]})/identity_warp[$j].nii.gz -o $(dirname ${search_sessions[$l]})/subj_2_template_mrtrix_warp${j}.nii.gz -r ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz \
            #     -t $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_1Warp.nii.gz -t $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_0GenericAffine.mat --default-value 21474836472
            # done
            # # forward (will use same name as original script for resulting warps)
            # warpcorrect $(dirname ${search_sessions[$l]})/subj_2_template_mrtrix_warp[].nii.gz $(dirname ${search_sessions[$l]})/subject2template_warp.mif -marker 21474836472 -force
            # # backward
            # warpinit ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz $(dirname ${search_sessions[$l]})/inv_identity_warp[].nii.gz -force
            # for x in {0..2}; do
            #     # backward
            #     antsApplyTransforms -d 3 -i $(dirname ${search_sessions[$l]})/inv_identity_warp[$x].nii.gz -o $(dirname ${search_sessions[$l]})/template_2_subject_mrtrix_warp${x}.nii.gz -r $(dirname ${search_sessions[$l]})/wmfod_norm_vol1.nii.gz \
            #     -t \[$(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_0GenericAffine.mat,1\] -t $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_1InverseWarp.nii.gz --default-value 21474836475
            # done
            # # backward
            # warpcorrect $(dirname ${search_sessions[$l]})/template_2_subject_mrtrix_warp[].nii.gz $(dirname ${search_sessions[$l]})/template2subject_warp.mif -marker 21474836475 -force
            
            if [[ ! -f "$(dirname ${search_sessions[$l]})/dwi_preproced${suffix}_mask_in_Template.mif" ]]; then
                # warpcorrect $(dirname ${search_sessions[$l]})/subj_2_template_mrtrix_warp[].nii.gz $(dirname ${search_sessions[$l]})/subj_2_template_mrtrix_warp_corrected.mif -marker 2147483647
                # should follow this with mrtransform - VVVI DO NOT USE -DATATYPE BIT when transforming the FODs!
                mrtransform -template ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz -warp $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_deformation.mif \
                ${search_sessions[$l]} $(dirname ${search_sessions[$l]})/wmfod_in_template_space.mif -nthreads $ncpu -force -reorient_fod 0
                # also transform the nan brain masks and min lesion ones
                if [[ -f "${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif" ]]; then
                    mrtransform -template ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz -warp $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_deformation.mif \
                    -datatype bit ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask_minL_nan.mif $(dirname ${search_sessions[$l]})/dwi_preproced${suffix}_mask_minL_nan_in_Template.mif -nthreads $ncpu -force -reorient_fod 0
                    if [[ ! -f "$(dirname ${search_sessions[$l]})/wmfod_in_template_space_naned.mif" ]]; then
                        mrthreshold $(dirname ${search_sessions[$l]})/dwi_preproced${suffix}_mask_minL_nan_in_Template.mif -abs 0.0 -comparison gt -nan $(dirname ${search_sessions[$l]})/dwi_preproced${suffix}_mask_minL_renan_in_Template.mif -force
                        mrcalc $(dirname ${search_sessions[$l]})/wmfod_in_template_space.mif \
                        $(dirname ${search_sessions[$l]})/dwi_preproced${suffix}_mask_minL_renan_in_Template.mif -mult $(dirname ${search_sessions[$l]})/wmfod_in_template_space_naned.mif -force
                    fi
                fi
                mrtransform -template ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz -warp $(dirname ${search_sessions[$l]})/wmfod_norm_vol1_2_template_deformation.mif \
                -datatype bit ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/dwi_preproced${suffix}_mask.mif $(dirname ${search_sessions[$l]})/dwi_preproced${suffix}_mask_in_Template.mif -nthreads $ncpu -force -reorient_fod 0
            fi
        done
    else
        # at this point ssst, ssmt, and msmt are the same no?
        for_each -info -force -nthreads ${ncpu_foreach} * : mrregister -mask1 IN/dwi_preproced${suffix}_mask.mif -mask2 ${cwd}/dwiprep/${group_name}/fba/template/template_mask.nii.gz \
        -nl_warp IN/subject2template_warp.mif IN/template2subject_warp.mif -nthreads $ncpu -force IN/wmfod_norm.mif ../template/wmfod_template.mif
    fi

    if [ $? -eq 0 ]; then
        echo "done" > ../fod_reg2template.done
    fi
else 
    echo " Registration of all subject FOD images to the FOD template already done "
fi

# Insert NANs into the FODs in template space

# now need to create a fake BIDS dir
mkdir -p ${cwd}/dwiprep/${group_name}/fba/template/BIDS
mkdir -p ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives
mkdir -p ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/freesurfer
mkdir -p ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/Temp_4FWT
mkdir -p ${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp
mkdir -p ${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat

if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat/sub-Temp_T1w.nii.gz" ]];then
    mrconvert ${cwd}/dwiprep/${group_name}/fba/template/template_T1w.mif ${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat/sub-Temp_T1w.nii.gz -force
fi

cp ${cwd}/BIDS/*.json ${cwd}/dwiprep/${group_name}/fba/template/BIDS/
cp ${cwd}/BIDS/*.tsv ${cwd}/dwiprep/${group_name}/fba/template/BIDS/
cp ${cwd}/BIDS/sub-$(basename $(dirname ${search_sessions[0]}))/anat/sub-$(basename $(dirname ${search_sessions[0]}))_T1w.json \
${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat/sub-Temp_T1w.json
cp ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/Temp_4FWT/fod_template_firstvol.nii.gz
cp ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template.mif ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/Temp_4FWT/fod_template.mif
cp ${cwd}/dwiprep/${group_name}/fba/template/template_mask.nii.gz ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/Temp_4FWT/brain_fod_temp_mask.nii.gz
# if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat/sub-Temp_T1w.json" ]]; then
#     mrinfo ${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat/sub-Temp_T1w.nii.gz -json_all ${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat/sub-Temp_T1w.json -force
# fi
# need to copy T1 json - done
# need to copy rest of BIDS folder content - done (should be)

# then run FS and MSBP for the template
# double check if this will work with output dir including subject folder
# recon-all -i ${cwd}/dwiprep/${group_name}/fba/template/BIDS/sub-Temp/anat/sub-Temp_T1w.nii.gz \
# -s sub-Temp -sd ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/freesurfer/sub-Temp -parallel -openmp ${ncpu} -all

# run MSBP docker
if [[ ! -f "${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/cmp/sub-Temp/anat/sub-Temp_label-L2018_desc-scale3_atlas.nii.gz" ]]; then
    docker run -u $(id -u) -it --rm -v ${cwd}/dwiprep/${group_name}/fba/template/BIDS:/bids_dir \
    -v ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives:/output_dir \
    -v $FREESURFER_HOME/license.txt:/opt/freesurfer/license.txt \
    sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
    --participant_label Temp --isotropic_resolution 1.0 --thalamic_nuclei \
    --brainstem_structures --skip_bids_validator --fs_number_of_cores ${ncpu} \
    --multiproc_number_of_cores ${ncpu} 2>&1 >> ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/MSBP_Template_log.txt
fi

# then run KUL_FWT for the template
# need to also make a dir that contains wmfod.mif, gm.mif, and csf.mif, brain_mask.nii.gz
KUL_FWT_make_VOIs_4Temp.sh -p Temp -M ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/cmp/sub-Temp/anat/sub-Temp_label-L2018_desc-scale3_atlas.nii.gz \
-F ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/freesurfer/sub-Temp/mri/aparc+aseg.mgz \
-d ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/Temp_4FWT -c ${cwd}/FWT_config.txt -n ${ncpu} \
-o ${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output
KUL_FWT_make_TCKs_4Temp.sh -p Temp -M ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/cmp/sub-Temp/anat/sub-Temp_label-L2018_desc-scale3_atlas.nii.gz \
-F ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/freesurfer/sub-Temp/mri/aparc+aseg.mgz \
-d ${cwd}/dwiprep/${group_name}/fba/template/BIDS/derivatives/Temp_4FWT -c ${cwd}/FWT_config.txt -n ${ncpu} \
-o ${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output -T 1 -S -f 1 -Q

# exit 2

# Compute a white matter template analysis fixel mask
if [ ! -d ../template/fixel_mask ]; then
    echo "   Compute a white matter template analysis fixel mask"
    fod2fixel -mask ../template/template_mask.mif -fmls_peak_value 0.10 ../template/wmfod_template.mif ../template/fixel_mask -nthreads $ncpu -force
else
    echo " Computation of a whole brain white matter template analysis fixel mask already done "
fi

# Template bundle specific fixels - just for show
# create template bundle fixels only if the user asked for them
Tcks_4fba=($(ls -f ${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output/sub-Temp_TCKs_output/*/*fin_BT_iFOD2.tck))
Tckmaps_4fba=($(ls -f ${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output/sub-Temp_TCKs_output/*/*fin_map_BT_iFOD2.nii.gz))
search_sessions=($(ls -f ${cwd}/dwiprep/${group_name}/fba/subjects/*/wmfod_norm.mif))
if [[ $w_flag == 2 ]] || [[ $w_flag == 4 ]]; then
    for T in ${!Tcks_4fba[@]}; do
            tck_name=$(basename ${Tcks_4fba[$T]} _fin_BT_iFOD2.tck)
            mkdir -p ../template/${tck_name}
        if [[ ! -f "../template/${tck_name}/${tck_name}_temp_FODs/${tck_name}_temp_fixels/afd.mif" ]]; then
            # tck2fixel ${Tcks_4fba[$T]} ../template/fixel_mask ../template/${tck_name}_temp_fixels ${tck_name}_temp_fixels.mif -nthreads $ncpu -force
            # first smooth the map, then binarize and add nans
            maskfilter ${Tckmaps_4fba[$T]} dilate - | mrthreshold - -abs 0.0 -comparison gt $(dirname ${Tckmaps_4fba[$T]})/$(basename ${Tckmaps_4fba[$T]} .nii.gz)_dil_nan.mif -nthreads $ncpu -force
            # now we use this to mask the FOD of the template (we get an FOD segment per bundle then)
            mrcalc ../template/wmfod_template.mif $(dirname ${Tckmaps_4fba[$T]})/$(basename ${Tckmaps_4fba[$T]} .nii.gz)_dil_nan.mif -mult ../template/${tck_name}/${tck_name}_temp_FODs.mif -force
            # now we convert said FOD segments to fixels (Adding all measures we wanted)
            fod2fixel -mask $(dirname ${Tckmaps_4fba[$T]})/$(basename ${Tckmaps_4fba[$T]} .nii.gz)_dil_nan.mif ../template/${tck_name}/${tck_name}_temp_FODs.mif \
            ../template/${tck_name}/${tck_name}_temp_fixels/ -fmls_peak_value 0.10 -afd afd.mif -peak_amp peak_amp.mif -disp disp.mif -force -nthreads ${ncpu}
            fixelconnectivity ../template/${tck_name}/${tck_name}_temp_fixels ${Tcks_4fba[$T]} ../template/${tck_name}/${tck_name}_temp_conn_matrix
            # create dirs for afd, disp, peak_amp here per bundle?
            mkdir -p ../template/${tck_name}/afd_work
            mkdir -p ../template/${tck_name}/afd_smooth
            mkdir -p ../template/${tck_name}/disp_work
            mkdir -p ../template/${tck_name}/disp_smooth
            mkdir -p ../template/${tck_name}/peak_amp_work
            mkdir -p ../template/${tck_name}/peak_amp_smooth
            mkdir -p ../template/${tck_name}_fin
            mkdir -p ../template/${tck_name}_fin/afd_fin
            mkdir -p ../template/${tck_name}_fin/disp_fin
            mkdir -p ../template/${tck_name}_fin/peak_amp_fin
            
            # copy dirs and index to fin folder
            cp ../template/${tck_name}/${tck_name}_temp_fixels/index.mif ../template/${tck_name}_fin/afd_fin/index.mif
            cp ../template/${tck_name}/${tck_name}_temp_fixels/index.mif ../template/${tck_name}_fin/disp_fin/index.mif
            cp ../template/${tck_name}/${tck_name}_temp_fixels/index.mif ../template/${tck_name}_fin/peak_amp_fin/index.mif
            cp ../template/${tck_name}/${tck_name}_temp_fixels/directions.mif ../template/${tck_name}_fin/afd_fin/directions.mif
            cp ../template/${tck_name}/${tck_name}_temp_fixels/directions.mif ../template/${tck_name}_fin/disp_fin/directions.mif
            cp ../template/${tck_name}/${tck_name}_temp_fixels/directions.mif ../template/${tck_name}_fin/peak_amp_fin/directions.mif
        else
            echo " Computation of a ${tck_name} white matter template analysis fixel mask already done "
        fi
    done
fi

# exit 2

# subjects fod2fixel, gen afd, disp and peak_amp (non-reoriented fixels))
for ww in ${!search_subjects[@]}; do
    si=$(basename ${search_subjects[$ww]} | cut -d "_" -f1 | cut -d "-" -f2)
    # The lesion mask is taken into account when creating fixels from fods
    if [[ -f "$(dirname ${search_sessions[$ww]})/dwi_preproced${suffix}_mask_minL_renan_in_Template.mif" ]]; then
        mask_4fod2fix=$(dirname ${search_sessions[$ww]})/dwi_preproced${suffix}_mask_minL_renan_in_Template.mif
    else
        mask_4fod2fix=$(dirname ${search_sessions[$ww]})/dwi_preproced${suffix}_mask_in_Template.mif
    fi
    if [[ -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space_naned.mif" ]] && [[ ! -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_noreori" ]]; then
        # Whole brain naned (if needed) fixels FODs and fixels are created here for each subject
        # then convert to fixels with all the fixins!
        fod2fixel -mask ${mask_4fod2fix} -fmls_peak_value 0.10 $(dirname ${search_sessions[$ww]})/wmfod_in_template_space_naned.mif \
        ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_noreori -nthreads ${ncpu} -force -afd afd.mif -disp disp.mif -peak_amp peak_amp.mif
    elif [[ ! -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space_naned.mif" ]] && [[ -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space.mif" ]] \
    && [[ ! -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_noreori" ]]; then
        fod2fixel -mask ${mask_4fod2fix} -fmls_peak_value 0.10 $(dirname ${search_sessions[$ww]})/wmfod_in_template_space.mif \
        ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels -nthreads ${ncpu} -force -afd afd.mif -disp disp.mif -peak_amp peak_amp.mif
    elif [[ -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_noreori" ]] || [[ -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_noreori" ]]; then
        echo " Fixels already generated for ${si}, skipping"
    else
        echo " we have a problem, no appropriate wmfod maps in template space found for ${si}, exiting "
        exit 2
    fi

    # reorient the fixels!
    # at this point we don't really care whether the output is naned or not so we make no distinction in naming
    if [[ -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_noreori" ]] && [[ ! -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_reori" ]]; then
        fixelreorient ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_noreori $(dirname ${search_sessions[$ww]})/wmfod_norm_vol1_2_template_deformation.mif \
        ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_reori -nthreads $ncpu -force
        # get fixelcorrespondence done also
        if [[ ! -f "../template/peak_amp/${si}_peak_amp.mif" ]]; then
            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_reori/afd.mif ../template/fixel_mask ../template/afd ${si}_afd.mif
            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_reori/disp.mif ../template/fixel_mask ../template/disp ${si}_disp.mif
            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_reori/peak_amp.mif ../template/fixel_mask ../template/peak_amp ${si}_peak_amp.mif
        fi
    elif [[ -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_noreori" ]] && [[ ! -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_reori" ]]; then
        fixelreorient ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_noreori $(dirname ${search_sessions[$ww]})/wmfod_norm_vol1_2_template_deformation.mif \
        ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_reori -nthreads $ncpu -force
        if [[ ! -f "../template/peak_amp/${si}_peak_amp.mif" ]]; then
            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_reori/afd.mif ../template/fixel_mask ../template/afd ${si}_afd.mif
            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_reori/disp.mif ../template/fixel_mask ../template/disp ${si}_disp.mif
            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_reori/peak_amp.mif ../template/fixel_mask ../template/peak_amp ${si}_peak_amp.mif
        fi        
    elif [[ -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_reori" ]] || [[ -d "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_fixels_naned_reori" ]]; then
        echo " Fixels already reoriented and corresponded to template for ${si}, skipping"
    fi

    # this should only run if tract specific work was chosen by the user
    if [[ $w_flag == 2 ]] || [[ $w_flag == 4 ]]; then
        for T in ${!Tcks_4fba[@]}; do
            tck_name=$(basename ${Tcks_4fba[$T]} _fin_BT_iFOD2.tck)
            if [[ ! -f "../template/${tck_name}_fin/peak_amp_fin/${si}_peak_amp.mif" ]]; then
                if [[ -f "$(dirname ${search_sessions[$ww]})/dwi_preproced${suffix}_mask_minL_renan_in_Template.mif" ]]; then
                    rmasked_bundle=$(dirname ${search_sessions[$ww]})/${tck_name}_mask_minL_renan_in_Template.mif
                else
                    rmasked_bundle=$(dirname ${search_sessions[$ww]})/${tck_name}_mask_in_Template.mif
                fi
                if [[ ! -f "../template/${tck_name}/${tck_name}_FODs.mif" ]]; then
                    # the naned suffix given to the output of this step is to denote the nans added outside of the bundle mask
                    if [[ -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space_naned.mif" ]] && [[ ! -f "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_reori/${si}_afd.mif" ]]; then
                        input_FOD="$(dirname ${search_sessions[$ww]})/wmfod_in_template_space_naned.mif"
                    elif [[ ! -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space_naned.mif" ]] && [[ -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space.mif" ]] \
                    && [[ ! -f "${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_reori/${si}_afd.mif" ]]; then
                        input_FOD="$(dirname ${search_sessions[$ww]})/wmfod_in_template_space.mif"
                    elif [[ ! -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space_naned.mif" ]] && [[ ! -f "$(dirname ${search_sessions[$ww]})/wmfod_in_template_space.mif" ]]; then
                        echo " we have a problem, neither the wmfod_in_template_space.mif nor the wmfod_in_template_space_naned.mif were found for this subject ${si}, exiting"
                        exit 2
                    fi
                    # apply the dilated and naned template bundle mask to subject FODs directly
                    if [[ -f "${input_FOD}" ]]; then
                        mrcalc $(dirname ${Tckmaps_4fba[$T]})/$(basename ${Tckmaps_4fba[$T]} .nii.gz)_dil_nan.mif ${mask_4fod2fix} -mult ${rmasked_bundle} -force
                        mrcalc ${input_FOD} ${rmasked_bundle} -mult $(dirname ${search_sessions[$ww]})/${si}_${tck_name}_FODs_naned.mif -force
                        # generate naned fixels for each bundle of interest
                        fod2fixel -mask $(dirname ${Tckmaps_4fba[$T]})/$(basename ${Tckmaps_4fba[$T]} .nii.gz)_dil_nan.mif -fmls_peak_value 0.10 $(dirname ${search_sessions[$ww]})/${si}_${tck_name}_FODs_naned.mif \
                        ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_noreori -nthreads ${ncpu} -force -afd afd.mif -disp disp.mif -peak_amp peak_amp.mif
                        # reorient the fixels to match those of the template
                        fixelreorient ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_noreori $(dirname ${search_sessions[$ww]})/wmfod_norm_vol1_2_template_deformation.mif \
                        ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_reori -nthreads $ncpu -force
                        # fixelcorrespondence per parametric map
                        if [[ ! -f "../template/${tck_name}/${si}_${tck_name}_temp_fixels_peak_amp ${si}_peak_amp.mif" ]]; then 
                            # 08/07
                            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_reori/afd.mif ../template/${tck_name}/${tck_name}_temp_fixels ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_afd ${si}_afd.mif -nthreads $ncpu -force
                            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_reori/disp.mif ../template/${tck_name}/${tck_name}_temp_fixels ../template/${tck_name}/disp_work/${si}_${tck_name}_temp_fixels_disp ${si}_disp.mif -nthreads $ncpu -force
                            fixelcorrespondence ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/${si}_${tck_name}_fixels_naned_reori/peak_amp.mif ../template/${tck_name}/${tck_name}_temp_fixels ../template/${tck_name}/peak_amp_work/${si}_${tck_name}_temp_fixels_peak_amp ${si}_peak_amp.mif -nthreads $ncpu -force
                        fi

                        if [[ ! -d "../template/${tck_name}/peak_amp_work/${si}_${tck_name}_temp_fixels_peak_amp_renaned" ]]; then
                            # this needs to be done only once (using afd) 
                            # then the resulting mask applied to disp, peak_amp and afd using fixelcrop
                            voxel2fixel ${rmasked_bundle} ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_afd ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_mask ${si}_mask.mif -nthreads $ncpu -force
                            # smooth the fixels first, then crop... do we need to smooth the fixel mask too?
                            fixelfilter -matrix ../template/${tck_name}/${tck_name}_temp_conn_matrix ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_afd smooth ../template/${tck_name}/afd_smooth/${si}_${tck_name}_temp_fixels_afd_smooth
                            fixelfilter -matrix ../template/${tck_name}/${tck_name}_temp_conn_matrix ../template/${tck_name}/disp_work/${si}_${tck_name}_temp_fixels_disp smooth ../template/${tck_name}/disp_smooth/${si}_${tck_name}_temp_fixels_disp_smooth
                            fixelfilter -matrix ../template/${tck_name}/${tck_name}_temp_conn_matrix ../template/${tck_name}/peak_amp_work/${si}_${tck_name}_temp_fixels_peak_amp smooth ../template/${tck_name}/peak_amp_smooth/${si}_${tck_name}_temp_fixels_peak_amp_smooth

                            mrcalc ../template/${tck_name}/afd_smooth/${si}_${tck_name}_temp_fixels_afd_smooth/${si}_afd.mif ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_mask/${si}_mask.mif -mult ../template/${tck_name}/afd_smooth/${si}_${tck_name}_temp_fixels_afd_smooth/${si}_afd_naned.mif -force -nthreads $ncpu
                            mrcalc ../template/${tck_name}/disp_smooth/${si}_${tck_name}_temp_fixels_disp_smooth/${si}_disp.mif ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_mask/${si}_mask.mif -mult ../template/${tck_name}/disp_smooth/${si}_${tck_name}_temp_fixels_disp_smooth/${si}_disp_naned.mif -force -nthreads $ncpu
                            mrcalc ../template/${tck_name}/peak_amp_smooth/${si}_${tck_name}_temp_fixels_peak_amp_smooth/${si}_peak_amp.mif ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_mask/${si}_mask.mif -mult ../template/${tck_name}/peak_amp_smooth/${si}_${tck_name}_temp_fixels_peak_amp_smooth/${si}_peak_amp_naned.mif -force -nthreads $ncpu
                            
                            # copy fin outputs to fin folder
                            cp ../template/${tck_name}/afd_smooth/${si}_${tck_name}_temp_fixels_afd_smooth/${si}_afd_naned.mif ../template/${tck_name}_fin/afd_fin/${si}_afd.mif
                            cp ../template/${tck_name}/disp_smooth/${si}_${tck_name}_temp_fixels_disp_smooth/${si}_disp_naned.mif ../template/${tck_name}_fin/disp_fin/${si}_disp.mif
                            cp ../template/${tck_name}/peak_amp_smooth/${si}_${tck_name}_temp_fixels_peak_amp_smooth/${si}_peak_amp_naned.mif ../template/${tck_name}_fin/peak_amp_fin/${si}_peak_amp.mif
                            # should fixelfilter come before or after fixelcrop?
                        fi
                        # can add if statement for fixel2peaks here
                    fi
                fi
            fi         
        done
    fi
done

# okay this works, now we need fixel correspondence?

# exit 2
# bundle specific fc, logfd and fdc
# should restructure output folders for afd, disp and peak_amp to be ../template/afd_tracts, ../template/disp_tracts, ../template/peak_amp_tracts
# and put everything where it belongs
# could probably do this outside of any for loops, just an if statement for w_flag == 2 or 4
# 08/07/2023 (above lines need to be done still)
### 08/07 19:46 -- This sucks! we need to generate one directory per subject per bundle to make sure no mistakes happen when applying masks to maps!

if [[ $w_flag == 2 ]] || [[ $w_flag == 4 ]]; then
    unset T ww
    # calc fc, logfc and fdc per subject per bundle
    for T in ${!Tcks_4fba[@]}; do
        tck_name=$(basename ${Tcks_4fba[$T]} _fin_BT_iFOD2.tck)
        mkdir -p ../template/${tck_name}/fc_work
        mkdir -p ../template/${tck_name}/fc_smooth
        mkdir -p ../template/${tck_name}_fin
        mkdir -p ../template/${tck_name}_fin/fc_fin
        mkdir -p ../template/${tck_name}_fin/fdc_fin
        mkdir -p ../template/${tck_name}_fin/log_fc_fin

        cp ../template/${tck_name}/${tck_name}_temp_fixels/index.mif ../template/${tck_name}_fin/fc_fin/index.mif
        cp ../template/${tck_name}/${tck_name}_temp_fixels/index.mif ../template/${tck_name}_fin/fdc_fin/index.mif
        cp ../template/${tck_name}/${tck_name}_temp_fixels/index.mif ../template/${tck_name}_fin/log_fc_fin/index.mif
        cp ../template/${tck_name}/${tck_name}_temp_fixels/directions.mif ../template/${tck_name}_fin/fc_fin/directions.mif
        cp ../template/${tck_name}/${tck_name}_temp_fixels/directions.mif ../template/${tck_name}_fin/fdc_fin/directions.mif
        cp ../template/${tck_name}/${tck_name}_temp_fixels/directions.mif ../template/${tck_name}_fin/log_fc_fin/directions.mif
        
        # mkdir -p ../template/${tck_name}/fc_${tck_name}/fc_renaned
        # mkdir -p ../template/${tck_name}/fc_${tck_name}/fdc_renaned
        # mkdir -p ../template/${tck_name}/fc_${tck_name}/fdc
        # mkdir -p ../template/${tck_name}/log_fc_${tck_name}
        # mkdir -p ../template/${tck_name}/log_fc_${tck_name}/log_fc_renaned
        for ww in ${!search_subjects[@]}; do
            si=$(basename ${search_subjects[$ww]} | cut -d "_" -f1 | cut -d "-" -f2)
            rmasked_bundle=$(dirname ${search_sessions[$ww]})/${tck_name}_mask_minL_renan_in_Template.mif
            # mkdir -p ../template/${tck_name}/fc_work/${si}_fc_work
            if [[ ! -f "../template/${tck_name}_fin/log_fc_fin/${si}_log_fc.mif" ]]; then
                # calc fc
                warp2metric $(dirname ${search_sessions[$ww]})/wmfod_norm_vol1_2_template_deformation.mif -fc ../template/${tck_name}/${tck_name}_temp_fixels ../template/${tck_name}/fc_work/${si}_fc_work ${si}_fc.mif -force
                # calc log of fc
                mrcalc ../template/${tck_name}/fc_work/${si}_fc_work/${si}_fc.mif -log ../template/${tck_name}/fc_work/${si}_fc_work/${si}_log_fc.mif -force
                # calc fdc
                mrcalc ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_afd/${si}_afd.mif ../template/${tck_name}/fc_work/${si}_fc_work/${si}_fc.mif -mult ../template/${tck_name}/fc_work/${si}_fc_work/${si}_fdc.mif -force
                # this should be replaced by fixelcrop using same mask from afd... they're the same dimensions no?
                # voxel2fixel ${rmasked_bundle} ../template/${tck_name}/fc_work/${si}_fc_work ../template/${tck_name}/log_fc_${tck_name}/${si}_fc_work_renaned ${si}.mif -nthreads $ncpu -force
                fixelfilter -matrix ../template/${tck_name}/${tck_name}_temp_conn_matrix ../template/${tck_name}/fc_work/${si}_fc_work/ smooth ../template/${tck_name}/fc_smooth/${si}_fc_smooth
                
                # renan the results
                mrcalc ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_fc.mif ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_mask/${si}_mask.mif -mult ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_fc_naned.mif
                mrcalc ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_fdc.mif ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_mask/${si}_mask.mif -mult ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_fdc_naned.mif
                mrcalc ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_log_fc.mif ../template/${tck_name}/afd_work/${si}_${tck_name}_temp_fixels_mask/${si}_mask.mif -mult ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_log_fc_naned.mif

                # copy fin results to fin folder
                cp ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_fc_naned.mif ../template/${tck_name}_fin/fc_fin/${si}_fc.mif
                cp ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_fdc_naned.mif ../template/${tck_name}_fin/fdc_fin/${si}_fdc.mif
                cp ../template/${tck_name}/fc_smooth/${si}_fc_smooth/${si}_log_fc_naned.mif ../template/${tck_name}_fin/log_fc_fin/${si}_log_fc.mif
            fi
        done
    done
fi
# exit 2

unset i x s

for o in ${!Tcks_4fba[@]}; do
    # mkdir -p ../template/${tck_name}_segs
    tck_name=$(basename ${Tcks_4fba[$o]} _fin_BT_iFOD2.tck)
    voxel2fixel -nthreads ${ncpu} -force ${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output/sub-Temp_TCKs_output/${tck_name}_output/QQ/${tck_name}_fin_BT_iFOD2_rs1c_segments_inMNI/labels_map_inTemp.nii.gz \
    ${cwd}/dwiprep/${group_name}/fba/template/${tck_name}/${tck_name}_temp_fixels ${cwd}/dwiprep/${group_name}/fba/template/${tck_name}/${tck_name}_temp_bundle_segs_fixels_map ${tck_name}_temp_segs.mif
    cp ${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output/sub-Temp_TCKs_output/${tck_name}_output/QQ/${tck_name}_fin_BT_iFOD2_rs1c_segments_inMNI/labels_map_inTemp.nii.gz \
    ${cwd}/dwiprep/${group_name}/fba/template/${tck_name}/${tck_name}_temp_bundle_segs_voxels_map.nii.gz

    echo "sub-ID, segment_ID, AFD, log_fc, fdc, disp" >> ../template/mean_scores_fba_${tck_name}.txt
done

${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output/sub-Temp_TCKs_output/*/*fin_map_BT_iFOD2.nii.gz

for ww in ${!search_subjects[@]}; do
    si=$(basename ${search_subjects[$ww]} | cut -d "_" -f1 | cut -d "-" -f2)
    for x in ${!Tcks_4fba[@]}; do
        tck_name=$(basename ${Tcks_4fba[$x]} _fin_BT_iFOD2.tck)
        bundle_seg_map="${cwd}/dwiprep/${group_name}/fba/template/FWT_Temp_output/sub-Temp_TCKs_output/${tck_name}_output/QQ/${tck_name}_fin_BT_iFOD2_rs1c_segments_inMNI/labels_map_inTemp.nii.gz"
        for y in {1..101}; do
            # mrcalc ${bundle_seg_map} ${y} -eq ../template/${tck_name}_segs/segment_${y}.mif
            afd=$(mrstats -mask $(mrcalc ${cwd}/dwiprep/${group_name}/fba/template/${tck_name}/${tck_name}_temp_bundle_segs_fixels_map/${tck_name}_temp_segs.mif ${y} -eq - ) \
            -ignorezero -output mean ../template/${tck_name}_fin/afd_fin/${si}_afd.mif)
            log_fc=$(mrstats -mask $(mrcalc ${cwd}/dwiprep/${group_name}/fba/template/${tck_name}/${tck_name}_temp_bundle_segs_fixels_map/${tck_name}_temp_segs.mif ${y} -eq - ) \
            -ignorezero -output mean ../template/${tck_name}_fin/log_fc_fin/${si}_log_fc.mif)
            fdc=$(mrstats -mask $(mrcalc ${cwd}/dwiprep/${group_name}/fba/template/${tck_name}/${tck_name}_temp_bundle_segs_fixels_map/${tck_name}_temp_segs.mif ${y} -eq - ) \
            -ignorezero -output mean ../template/${tck_name}_fin/fdc_fin/${si}_fdc.mif)
            disp=$(mrstats -mask $(mrcalc ${cwd}/dwiprep/${group_name}/fba/template/${tck_name}/${tck_name}_temp_bundle_segs_fixels_map/${tck_name}_temp_segs.mif ${y} -eq - ) \
            -ignorezero -output mean ../template/${tck_name}_fin/disp_fin/${si}_disp.mif)
            echo "${si}, ${y}, ${afd}, ${log_fc}, ${fdc}, ${disp}" >> ../template/mean_scores_fba_${tck_name}.txt

        done
    done
done

exit 2

# change the transforms below to match mrregister output - or forget about it!
# Warp FOD images to template space
# if [ ! -f ../fod_warp.done ]; then
#     echo "   Warping FOD images to template space"
#     # without reorienting
#     for_each -info -force -nthreads ${ncpu_foreach} * : mrtransform IN/wmfod_norm.mif -template ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz -warp_full $(dirname ${search_sessions[$l]})/subject2template_warp.mif \
#     IN/fod_in_template_space_NOT_REORIENTED.mif -reorient_fod 0 -nthreads $ncpu -force
#     # with reorienting
#     for_each -info -force -nthreads ${ncpu_foreach} * : mrtransform IN/wmfod_norm.mif -template ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz -warp_full $(dirname ${search_sessions[$l]})/subject2template_warp.mif \
#     IN/fod_in_template_space_REORIENTED.mif -reorient_fod 1 -nthreads $ncpu -force
#     # if [ $? -eq 0 ]; then
#     #     echo "done" > ../fod_warp.done
#     # fi
# else
#     echo "   Warping FOD images to template space already done"
# fi  

# # this is redundant - we will use the (thresholded?) template bundle fixels to sample subjects AFD directly
# # sampling the subject specific fods using the template bundle fixels (minimize variance across subject)
# if [[ ! -f ../sub_bundle_fixels.done ]]; then
#     for o in ${!search_subjects[@]}; do
#         se=$(basename ${search_subjects[$o]} | cut -d "_" -f1 | cut -d "-" -f2)
#         for z in ${!Tcks_4fba[@]}; do
#             tck_name=$(basename ${Tcks_4fba[$z]} _fin_BT_iFOD2.tck)
#             ind_fix="${cwd}/dwiprep/${group_name}/fba/subjects/${se}/${se}_${tck_name}_fixels_naned"           
#             if [[ -d "${cwd}/dwiprep/${group_name}/fba/subjects/${se}/${se}_fixels_naned" ]]; then
#                 # new steps are as follow 07/07
                
#                 tck2fixel ${Tcks_4fba[$z]} ${cwd}/dwiprep/${group_name}/fba/subjects/${se}/${se}_fixels_naned ${ind_fix} ${se}_${tck_name}_temp_fixels.mif -nthreads $ncpu -force
#             else
#                 tck2fixel ${Tcks_4fba[$z]} ${cwd}/dwiprep/${group_name}/fba/subjects/${se}/${se}_fixels ${ind_fix} ${se}_${tck_name}_temp_fixels.mif -nthreads $ncpu -force
#             fi
#         done
#     done
# fi

# think about what needs to be added to list of voxelwise parametric maps
# Make FA/ADC images in template space
if [ ! -f ../fa_adc_warp.done ]; then
    # find the FA in subject space
    search_sessions=($(find ${cwd}/dwiprep/sub-* -type f | grep qa/fa${suffix}.nii.gz | sort ))
    # num_sessions=${#search_sessions[@]}
    for i in ${search_sessions[@]}; do
        sub=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
        ses=$(echo $i | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
        if [[ ! -z ${ses} ]]; then
            s=${sub}_${ses}
        else
            s=${sub}
        fi
        echo $i
        echo $s
        mrconvert $i ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/FA_subj_space.mif -force
    done
    # find the ADC in subject space
    search_sessions=($(find ${cwd}/dwiprep/sub-* -type f | grep qa/adc${suffix}.nii.gz | sort ))
    # num_sessions=${#search_sessions[@]}
    for i in ${search_sessions[@]}; do
        # s=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
        # s=$(echo $i | awk -F 'subjects/' '{print $2}' | awk -F '/' '{print $1}')
        sub=$(echo $i | awk -F 'sub-' '{print $2}' | awk -F '/' '{print $1}')
        ses=$(echo $i | awk -F 'ses-' '{print $2}' | awk -F '/' '{print $1}')
        # s=${sub}_${ses}
        if [[ ! -z ${ses} ]]; then
            s=${sub}_${ses}
        else
            s=${sub}
        fi
        mrconvert $i ${cwd}/dwiprep/${group_name}/fba/subjects/${s}/ADC_subj_space.mif -force
    done
    echo "   Warping FA/ADC images to template space"
    for_each -info -force -nthreads ${ncpu_foreach} * : mrtransform IN/FA_subj_space.mif -template ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz \
    -warp_full IN/wmfod_norm_vol1_2_template_warp.mif IN/FA_in_template_space.nii.gz -nthreads $ncpu -force
    for_each -info -force -nthreads ${ncpu_foreach} * : mrtransform IN/ADC_subj_space.mif -template ${cwd}/dwiprep/${group_name}/fba/template/wmfod_template_vol1.nii.gz \
    -warp_full IN/wmfod_norm_vol1_2_template_warp.mif IN/ADC_in_template_space.nii.gz -nthreads $ncpu -force
    if [ $? -eq 0 ]; then
        echo "done" > ../fa_adc_warp.done
    fi
    mkdir -p ../template/fa
    #ln -sfn ${cwd}/dwiprep/${group_name}/fba/dwiintensitynorm/dwi_output
    for_each -info -force -nthreads ${ncpu_foreach} * : ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/IN/FA_in_template_space.nii.gz ${cwd}/dwiprep/${group_name}/fba/template/fa/sub_IN_FA.nii.gz
    mkdir -p ../template/adc
    for_each -info -force -nthreads ${ncpu_foreach} * : ln -sfn ${cwd}/dwiprep/${group_name}/fba/subjects/IN/ADC_in_template_space.nii.gz ${cwd}/dwiprep/${group_name}/fba/template/adc/sub_IN_ADC.nii.gz
else
    echo "   Warping FA/ADC images to template space already done"
fi  

exit 2

# # Segment FOD images to estimate fixels and their apparent fibre density (FD)
# if [ ! -f ../fod_segment.done ]; then
#     echo "   Segment FOD images to estimate fixels and their apparent fibre density (FD)"
#     # whole brain WM
#     for_each -info -force -nthreads ${ncpu_foreach} * : fod2fixel -mask ../template/template_mask.mif IN/wmfod_in_template_space.mif \
#     IN/fixel_in_template_space_NOT_REORIENTED -afd fd.mif -nthreads $ncpu -force
#     # bundle specific - need input masks in template space per bundle with lesion voxels set to zero
#     for_each -info -force -nthreads ${ncpu_foreach} * : fod2fixel -mask ${cwd}/dwiprep/${group_name}/fba/subjects/${si}/dwi_preproced${suffix}_mask_minL.mif \
#     IN/fixel_in_template_space_NOT_REORIENTED -afd fd.mif -nthreads $ncpu -force
#     if [ $? -eq 0 ]; then
#         echo "done" > ../fod_segment.done
#     fi
# else
#     echo "   Segmenting of FOD images to estimate fixels and their apparent fibre density (FD) already done"
# fi

# # Reorient fixels
# #foreach * : fixelreorient IN/fixel_in_template_space_NOT_REORIENTED IN/subject2template_warp.mif IN/fixel_in_template_space
# if [ ! -f ../fod_reor_fixels.done ]; then
#     # this is for the whole brain WM fixels
#     echo "   Reorient fixels"
#     for_each -info -force -nthreads ${ncpu_foreach} * : fixelreorient IN/fixel_in_template_space_NOT_REORIENTED IN/subject2template_warp.mif \
#     IN/fixel_in_template_space -nthreads $ncpu -force
#     # this is for the bundle specific WM fixels
#     if [ $? -eq 0 ]; then
#         echo "done" > ../fod_reor_fixels.done
#     fi
# else
#     echo "   Reorient fixels already done"
# fi

# # Assign subject fixels to template fixels
# # foreach * : fixelcorrespondence IN/fixel_in_template_space/fd.mif ../template/fixel_mask ../template/fd PRE.mif
# # Note: do NOT run in PARALLEL
# if [ ! -f ../assign_fixels.done ]; then
#     echo "   Assign subject fixels to template fixels"
#     # for_each -info -force -nthreads ${ncpu_foreach} * : fixelcorrespondence IN/fixel_in_template_space/fd.mif \
#     # ../template/fixel_mask ../template/fd PRE.mif -force
#     # for_each -info -force -nthreads 1 * : Mau_Leuv_Cai - this line was mistakenly erased!!!!!
#     if [ $? -eq 0 ]; then
#         echo "done" > ../assign_fixels.done
#     fi
# else
#     echo " Assign subject fixels to template fixels already done"
# fi


# Compute the fibre cross-section (FC) metric
# IT IS possible that the script will exit at this stage with an error
# This seems to relate to the presence of an index file in the template/fc folder
# IT runs just fine if simply restarted and doesn't quit with fd or fdc
# THESE steps are in accordance with the guide on https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/st_fibre_density_cross-section.html
if [ ! -f ../compute_fc.done ]; then
    echo "   Compute the fibre cross-section (FC) metric"
    # for_each -info -force -nthreads ${ncpu_foreach} * : warp2metric IN/subject2template_warp.mif -fc ../template/fixel_mask ../template/fc IN.mif -force
    for_each -info -force -nthreads 1 * : warp2metric IN/subject2template_warp.mif -fc ../template/fixel_mask ../template/fc IN.mif -force
    if [ $? -eq 0 ]; then
        echo "done" > ../compute_fc.done
    fi
else
    echo " Compute the fibre cross-section (FC) metric already done"
fi


if [ ! -f ../compute_log_fc.done ]; then
    echo "   Compute the fibre cross-section LOG-(FC) metric"
    mkdir -p ../template/log_fc
    cp ../template/fc/index.mif ../template/fc/directions.mif ../template/log_fc
    for_each -info -force * : mrcalc ../template/fc/IN.mif -log ../template/log_fc/IN.mif -force
   if [ $? -eq 0 ]; then
        echo "done" > ../compute_log_fc.done
    fi
else
    echo " Compute the fibre cross-section LOG-(FC) metric already done"
fi

# Compute a combined measure of fibre density and cross-section (FDC)
if [ ! -f ../compute_fdc.done ]; then

    echo " Compute a combined measure of fibre density and cross-section (FDC)"
    mkdir -p ../template/fdc
    cp ../template/fc/index.mif ../template/fdc
    cp ../template/fc/directions.mif ../template/fdc
    for_each -info -nthreads ${ncpu_foreach} * : mrcalc ../template/fd/IN.mif ../template/fc/IN.mif -mult ../template/fdc/IN.mif -force
   if [ $? -eq 0 ]; then
        echo "done" > ../compute_fdc.done
    fi
else
    echo " Compute the fibre cross-section LOG-(FC) metric already done"
fi

# Perform whole-brain fibre tractography on the FOD template
cd ../template
if [ ! -f ../tckgen.done ]; then

    n=20000000

    echo " Perform whole-brain fibre tractography on the FOD template"
    tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_image template_mask.mif \
     -mask template_mask.mif -select $n -cutoff 0.10 tracks_20_million.tck
    
    if [ $? -eq 0 ]; then
        echo "done" > ../tckgen.done
    fi
    
else
    
    echo " Whole brain fibre tractography already done"

fi



# Reduce biases in tractogram densities
# tcksift tracks_20_million.tck wmfod_template.mif tracks_2_million_sift.tck -term_number 200000
if [ ! -f ../tcksift.done ]; then

    echo " Reduce biases in tractogram densities"

    n=2000000

    tcksift tracks_20_million.tck wmfod_template.mif tracks_2_million_sift.tck -term_number $n
    
    if [ $? -eq 0 ]; then
        echo "done" > ../tcksift.done
    fi
    
else
    
    echo " Reduce biases in tractogram densities already done"

fi


