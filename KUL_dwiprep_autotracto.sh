#!/bin/bash

# set -x

# This script is meant for automated segmentation of whole brain tractograms in .tck

# Ahmed Radwan, Louise Emsell, Stefan Sunaert
# @ Ahmed Radwan: ahemd.radwan@kuleuven.be
# @ Louise Emsell: louise.emsell@kuleuven.be
# @ Stefan Sunaert: stefan.sunaert@uzleuven.be

# v: 1.0 
# 19/12/2019

# we have .tck whole brain files for tckedit.

## date and ncpu

d=$(date "+%Y-%m-%d_%H-%M-%S");

ncpu="18"

subjs=("s3");



FS_Ls=("M1_GM_LT"  "M1_WM_LT"  "M1_GM_RT"  "M1_WM_RT"  \
"S1_GM_LT"  "S1_WM_LT"  "S1_GM_RT"  "S1_WM_RT"  \
"IFG_PTr_GM_LT"  "IFG_PTr_WM_LT"  "IFG_PTr_GM_RT"  "IFG_PTr_WM_RT"  \
"IFG_POp_GM_LT"  "IFG_POp_WM_LT"  "IFG_POp_GM_RT"  "IFG_POp_WM_RT"  \
"IFG_POr_GM_LT"  "IFG_POr_WM_LT"  "IFG_POr_GM_RT"  "IFG_POr_WM_RT"  \
"STG_GM_LT"  "STG_WM_LT"  "STG_GM_RT"  "STG_WM_RT"  \
"bSTS_GM_LT"  "bSTS_WM_LT"  "bSTS_GM_RT"  "bSTS_WM_RT"  \
"MTG_GM_LT"  "MTG_WM_LT"  "MTG_GM_RT"  "MTG_WM_RT"  \
"ITG_GM_LT"  "ITG_WM_LT"  "ITG_GM_RT"  "ITG_WM_RT"  \
"SMG_GM_LT"  "SMG_WM_LT"  "SMG_GM_RT"  "SMG_WM_RT"  \
"Insula_GM_LT"  "Insula_WM_LT"  "Insula_GM_RT"  "Insula_WM_RT"  \
"Amyg_LT"  "Amyg_RT"  "Put_LT"  "Put_RT"  "Pall_LT"  "Pall_RT"  \
"Thal_LT"  "Thal_RT"  "Caud_LT"  "Caud_RT"  \
"Phippo_GM_LT"  "Phippo_WM_LT"  "Phippo_GM_RT"  "Phippo_WM_RT"  \
"cACC_GM_LT"  "cACC_WM_LT"  "cACC_GM_RT"  "cACC_WM_RT"  \
"rACC_GM_LT"  "rACC_WM_LT"  "rACC_GM_RT"  "rACC_WM_RT"  \
"PCC_GM_LT"  "PCC_WM_LT"  "PCC_GM_RT"  "PCC_WM_RT"  \
"iPCC_GM_LT"  "iPCC_WM_LT"  "iPCC_GM_RT"  "iPCC_WM_RT"  \
"Hippo_LT"  "Hippo_RT"  "vDC_LT"  "vDC_RT"  \
"Fusi_GM_LT"  "Fusi_WM_LT"  "Fusi_GM_RT"  "Fusi_WM_RT"  \
"periCalc_GM_LT"  "periCalc_WM_LT" "periCalc_GM_RT"  "periCalc_WM_RT"  \
"Front_lobeGM_LT"  "Front_lobeWM_LT"  "Front_lobeGM_RT"  "Front_lobeWM_RT"  \
"Occ_lobeGM_LT"  "Occ_lobeWM_LT"  "Occ_lobeGM_RT"  "Occ_lobeWM_RT"
"Temp_lobeGM_LT"  "Temp_lobeWM_LT"  "Temp_lobeGM_RT"  "Temp_lobeWM_RT"
"Cing_lobeGM_LT"  "Cing_lobeWM_LT"  "Cing_lobeGM_RT"  "Cing_lobeWM_RT"  \
"Pari_lobeGM_LT"  "Pari_lobeWM_LT"  "Pari_lobeGM_RT"  "Pari_lobeWM_RT"
"Acc_LT"  "Acc_RT" "LatOF_GM_LT"  "LatOF_WM_LT"  "LatOF_GM_RT"  "LatOF_WM_RT"  \
"MedOF_GM_LT"  "MedOF_WM_LT"  "MedOF_GM_RT"  "MedOF_WM_RT"  \
"FrontP_GM_LT"  "FrontP_WM_LT"  "FrontP_GM_RT"  "FrontP_WM_RT"  \
"TempP_GM_LT"  "TempP_WM_LT"  "TempP_GM_RT"  "TempP_WM_RT"  \
"CC_post"  "CC_midpost"  "CC_central"  "CC_midant"  "CC_ant"  \
"Unseg_WM_LT"  "Unseg_WM_RT"  "GMhemi_LT"  "GMhemi_RT"  "BStem"  "WMhemi_LT"  "WMhemi_RT"
"Cerebellum_WM_LT"  "Cerebellum_WM_RT"  "Cerebellum_GM_LT"  "Cerebellum_GM_RT"
"SPL_GM_LT"  "SPL_WM_LT"  "SPL_GM_RT"  "SPL_WM_RT"  \
"IPL_GM_LT"  "IPL_WM_LT"  "IPL_GM_RT"  "IPL_WM_RT"  \
"SFG_GM_LT"  "SFG_WM_LT"  "SFG_GM_RT"  "SFG_WM_RT"  \
"cMFG_GM_LT"  "cMFG_WM_LT"  "cMFG_GM_RT"  "cMFG_WM_RT"  \
"rMFG_GM_LT"  "rMFG_WM_LT"  "rMFG_GM_RT"  "rMFG_WM_RT"  \
"Pc_GM_LT"  "Pc_WM_LT"  "Pc_GM_RT"  "Pc_WM_RT"  \
"Fornix"  "TTG_GM_LT"  "TTG_WM_LT"  "TTG_GM_RT"  "TTG_WM_RT"  \
"Cuneus_GM_LT"  "Cuneus_WM_LT"  "Cuneus_GM_RT"  "Cuneus_WM_RT"  \
"LatOcc_GM_LT"  "LatOcc_WM_LT"  "LatOcc_GM_RT"  "LatOcc_WM_RT");

FS_vals=("1024"  "3024"  "2024"  "4024"  \
"1022"  "3022"  "2022"  "4022"  \
"1020"  "3020"  "2020"  "4020"  \
"1018"  "3018"  "2018"  "4018"  \
"1019"  "3019"  "2019"  "4019"  \
"1030"  "3030"  "2030"  "4030"  \
"1001"  "3001"  "2001"  "4001"
"1015"  "3015"  "2015"  "4015"  \
"1009"  "3009"  "2009"  "4009"  \
"1031"  "3031"  "2031"  "4031"  \
"1035"  "3035"  "2035"  "4035"  \
"18"  "54"  "12"  "51"  "13"  "52" \
"10"  "49"  "11"  "50"  \
"1016"  "3016"  "2016"  "4016"  \
"1002"  "3002"  "2002"  "4002"  \
"1026"  "3026"  "2026"  "4026"  \
"1023"  "3023"  "2023"  "4023"  \
"1010"  "3010"  "2010"  "4010"  \
"17"  "53"  "28"  "60"  \
"1007"  "3007"  "2007"  "4007"  \
"1021"  "3021"  "2021"  "4021"  \
"1001"  "3001"  "2001"  "4001"  \
"1004"  "3004"  "2004"  "4004"  \
"1005"  "3005"  "2005"  "4005"  \
"1003"  "3003"  "2003"  "4003"  \
"1006"  "3006"  "2006"  "4006"  \
"26"  "58"  "1012"  "3012"  "2012"  "4012"  \
"1014"  "3014"  "2014"  "4014"  \
"1032"  "3032"  "2032"  "4032"  \
"1033"  "3033"  "2033"  "4033"  \
"251"  "252"  "253"  "254"  "255"  \
"5001"  "5002"  "3"  "42"  "16"  "2"  "41"
"7"  "46"  "8"  "47"  \
"1029"  "3029"  "2029"  "4029"  \
"1008"  "3008"  "2008"  "4008"  \
"1028"  "3028"  "2028"  "4028"  \
"1003"  "3003"  "2003"  "4003"  \
"1027"  "3027"  "2027"  "4027"  \
"1025"  "3025"  "2025"  "4025"  \
"250"  "1034"  "3034"  "2034"  "4034"  \
"1005"  "3005"  "2005"  "4005"  \
"1011"  "3011"  "2011"  "4011");

# Misc vars

work_dir=$(pwd)

ROIs_dir="$(pwd)/FS_ROIs"

submit_dir="$(pwd)/Submission"

Seg_tck_dir="${submit_dir}/TCK_seg_out"

tck_out_dir="${submit_dir}/TCK_output"

tckmap_out_dir="${submit_dir}/TCKmap_output"

scr_grab_out_dir="${submit_dir}/Screen_grab_output"

tck_fun_logs="${sub_log}"

mkdir -p ${ROIs_dir}

mkdir -p ${submit_dir}

mkdir -p ${Seg_tck_dir}

mkdir -p ${tck_out_dir}

mkdir -p ${tckmap_out_dir}

mkdir -p ${scr_grab_out_dir}

# Template vars

temps_dir="${work_dir}/UKBB_BStem_template"

UKBB_temp="${temps_dir}/T1_preunbiased.nii.gz"

UKBB_temp_mask="${temps_dir}/T1_UKBB_brain_mask.nii.gz"

UKBB_BStem_labels_orig=("PTX_cst_l.nii.gz"  "PTX_cst_r.nii.gz"  "PTX_ml_l.nii.gz"  "PTX_ml_r.nii.gz"  \
"AC_midline.nii.gz"  "PC_midline.nii.gz"  "IFOF_exclude_midline.nii.gz"  "FAT_exclude_midline.nii.gz")

declare -a UKBB_BStem_inSubj

declare -a UKBB_BSL_thr

declare -a UKBB_BS_labels_r

#####

# Functions

# Exec_all function

function task_exec {

    echo "-------------------------------------------------------------" >> ${sub_log}

    echo ${task_in} >> ${sub_log}

    echo " Started @ $(date "+%Y-%m-%d_%H-%M-%S")" >> ${sub_log}

    eval ${task_in} >> ${sub_log} 2>&1 &

    echo " pid = $! basicPID = $$ " >> ${sub_log}

    wait ${pid}

    sleep 5

    if [ $? -eq 0 ]; then

        echo Success >> ${sub_log}

    else

        echo Fail >> ${sub_log}

        exit 1

    fi

    echo " Finished @ $(date "+%Y-%m-%d_%H-%M-%S")" >> ${sub_log}

    echo "-------------------------------------------------------------" >> ${sub_log}
    
    echo "" >> ${sub_log}

    unset task_in

}


# 1- General purpose thresholding function using mri_binarize

function KUL_mri_thr {

    echo "thresholded roi output in mgz, minimal interpolation" >> ${sub_log}

    task_in="mri_binarize --i ${input} --min ${min} --max ${max} --o ${output}"

    task_exec

}

# 2- mgz2nii conversion with mri_convert

function KUL_mri_mgz2nii {

    task_in="mri_convert -rl ${ref} -rt ${interp} -i ${mgz} -o ${nii}"

    task_exec

}

# 3- Function to make includes for tckedit

# ROIs are just a name, they're all VOIs

function KUL_mk_ROIs {

    ROI_o_str=${ROI_out##*/}

    srch_roi=($(find ${s_ROIs_dir} -type f | grep "${ROI_o_str}"))

    if [[ ! ${srch_roi} ]]; then
    
        task_in="mrcalc -force -nthreads ${ncpu} ${input1} ${input2} -add 0 -gt ${ROI_out} -force -nthreads ${ncpu}"

        task_exec

    else

        echo "ROI ${roi_o_str} already created, skipping " >> ${sub_log}

    fi


}

# tckedit function

# might need to add the filtering step also.

function KUL_tck_edit {

    task_in="tckedit -nthreads ${ncpu} -force -mask ${mask} -maxlength ${maxL} -minlength ${minL} ${wb_tck} ${includes_str} ${excludes_str} ${tck_out}"

    task_exec

}

# this function will take care of getting the UKBB template tracts to subject space and creating BStem ROIs

function KUL_get_UKBB_ROIs {

    # UKBB_temp2subj="${s_ROIs_dir}/UKBB_2_${subjs[$sub]}"

    UKBB_temp2subj="${s_ROIs_dir}/UKBB_2_FA_${subjs[$sub]}"

    UKBB_temp2subj_str=${UKBB_temp2subj##*/}

    srch_UKBB2S=($(find ${s_ROIs_dir} -type f | grep "${UKBB_temp2subj_str}_Warped.nii.gz"));

    if [[ ! ${srch_UKBB2S} ]]; then
    
        # need to bring the UKBB template to Diffusion space and not the subject's T1
        # so use this for the UKBB2Subect FA
        task_in="antsRegistration -d 3 -o [${UKBB_temp2subj}_,${UKBB_temp2subj}_Warped.nii.gz,${UKBB_temp2subj}_InverseWarped.nii.gz] \
        -x [${T1_brain_mask},${FA_nii_brain_mask},NULL] \
        -m MI[${UKBB_temp},${FA_nii},1,32,Regular,0.50] \
        -c [1000x500x250x0,1e-7,5] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 -u 1 -v 1\
        -m Mattes[${UKBB_temp},${FA_nii},1,64,Regular,0.50] -c [200x100x50,1e-7,5] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0mm \
        -u 1 -z 1 --winsorize-image-intensities [0.005, 0.995]"
        
        # task_in="antsRegistrationSyN.sh -d 3 -f ${UKBB_temp} -m ${T1_brain_nii} -x ${UKBB_temp_mask},${T1_brain_mask} -j 1 -t s -n 8 -o ${UKBB_temp2subj}_"

        task_exec

    else 

        echo "UKBB 2 Subj reg already done for ${subjs[$sub]}" >> ${sub_log}

    fi

    # CSTs need a thr of 75000

    # MLs need a thr of 100000

    srch_UKBBLs=($(find ${s_ROIs_dir} -type f | grep "FAT_exclude_midline_thr_in_${subjs[$sub]}_r.nii.gz"));

    if [[ ! ${srch_UKBBLs} ]]; then

        for lab in ${!UKBB_BStem_labels_orig[@]}; do

            thrs=("75000"  "75000"  "100000"  "100000"  "0.5"  "0.5"  "0.5"  "0.5");

            orig_L="${UKBB_BStem_labels_orig[$lab]}"

            ### edit this to work for each side

            if [[ ${orig_L} == *"_l."* ]]; then

                echo "working on left side UKBB rois" >> ${sub_log}

                UKBB_BSL_str="${orig_L::${#orig_L}-9}_LT"

            elif [[ ${orig_L} == *"_r."* ]]; then

                echo "working on right side UKBB rois" >> ${sub_log}

                UKBB_BSL_str="${orig_L::${#orig_L}-9}_RT"

            else

                echo "working on midline ROIs for Fornix, AC and PC" >> ${sub_log}

                UKBB_BSL_str="${orig_L::${#orig_L}-7}"

            fi

            UKBB_BStem_inSubj[$lab]="${UKBB_BSL_str}_thr_in_${subjs[$sub]}.nii.gz"

            UKBB_BSL_thr="${UKBB_BSL_str}_thr.nii.gz"

            UKBB_BS_labels_r[$lab]="${UKBB_BSL_str}_thr_in_${subjs[$sub]}_r.nii.gz"

            # UKBB_BStem_inSubj[$lab]="${UKBB_BSL_str}_thr${thrs[$lab]}k_in_${subj[$i]}.nii.gz"
            # UKBB_BS_labels_r[$lab]="${UKBB_BSL_str}_thr${thrs[$lab]}k_in_${subj[$i]}.nii.gz"
            # thresholds are 75k, 75k, 100k, 100k
            # zeroes are added in command

            # task_in="mrcalc -force -nthreads ${ncpu} -abs ${thrs[$lab]}000 ${temps_dir}/${orig_L} ${work_dir}/${subjs[$sub]}_${UKBB_BSL_thr[$lab]}"

            task_in="mrthreshold -ignorezero -force -nthreads ${ncpu} -abs ${thrs[$lab]} ${temps_dir}/${orig_L} - | maskfilter - -force -nthreads ${ncpu} -npass 2 dilate ${s_ROIs_dir}/${UKBB_BSL_thr}"

            task_exec

            # Apply the inverse Warp and affine transform to bring UKBB labels to diffusion space

            # task_in="WarpImageMultiTransform 3 ${s_ROIs_dir}/${UKBB_BSL_thr} ${s_ROIs_dir}/${UKBB_BStem_inSubj[$lab]} -R ${T1_brain_nii} -i \
            # ${UKBB_temp2subj}_0GenericAffine.mat ${UKBB_temp2subj}_1InverseWarp.nii.gz --use-NN"

            task_in="WarpImageMultiTransform 3 ${s_ROIs_dir}/${UKBB_BSL_thr} ${s_ROIs_dir}/${UKBB_BStem_inSubj[$lab]} -R ${FA_nii} -i \
            ${UKBB_temp2subj}_0GenericAffine.mat ${UKBB_temp2subj}_1InverseWarp.nii.gz"

            task_exec

            if [[ ${orig_L} == *"_cst_"* ]]; then

                echo "It's a CST atlas ROI, masking with BStem" >> ${sub_log}
                
                task_in="ImageMath 3 ${s_ROIs_dir}/${UKBB_BS_labels_r[$lab]} m ${s_ROIs_dir}/${UKBB_BStem_inSubj[$lab]} \
                ${s_ROIs_dir}/BStemr.nii.gz"

                task_exec

            elif [[ ${orig_L} == *"_ml_"* ]]; then

                echo "It's an ML atlas ROI, masking with BStem" >> ${sub_log}

                task_in="ImageMath 3 ${s_ROIs_dir}/${UKBB_BS_labels_r[$lab]} m ${s_ROIs_dir}/${UKBB_BStem_inSubj[$lab]} \
                ${s_ROIs_dir}/BStemr.nii.gz"

                task_exec

            elif [[ ${orig_L} == *"_exclude_"* ]]; then

                echo "It's an exclude ROI, no masking needed " >> ${sub_log}

                task_in="mv ${s_ROIs_dir}/${UKBB_BStem_inSubj[$lab]} ${s_ROIs_dir}/${UKBB_BS_labels_r[$lab]}"

                task_exec

            else

                echo "It's not a BStem ROI, no masking needed" >> ${sub_log}

                task_in="maskfilter -force -nthreads 6 -npass 2 ${s_ROIs_dir}/${UKBB_BStem_inSubj[$lab]} dilate ${s_ROIs_dir}/${UKBB_BS_labels_r[$lab]}"

                task_exec

            fi

        done

    else

        echo "UKBB labels already warped to subject space, skipping " >> ${sub_log}

    fi


}

function TCK_gen_n_filt {

    tck_name="${subjs[$sub]}_${tck_in}_${sides[$s]}"

    srch_tck_name=($(find ${s_tck_dir} -type f | grep "${tck_name}"))

    if [[ ! ${srch_tck_name} ]]; then

        mkdir -p "${s_tck_dir}/${tck_name}_filtering"

        tck_init="${s_tck_output_dir}/${tck_name}/${tck_name}_initial.tck"

        # added for reporting results and uploading

        mkdir -p "${s_tck_output_dir}/${tck_name}"

        mkdir -p "${s_tck_map_out_dir}/${tck_name}"

        mkdir -p "${s_tck_grabs_dir}/${tck_name}"

        mkdir -p "${s_tck_output_dir}/${tck_name}_ROIs"

        if [[ ${sides[$s]} == *"LT"* ]]; then

            echo "It's a left side bundle, excluding right hemi WM" >> ${sub_log}

            hemi_exclude="WMhemi_RTr.nii.gz"

        elif [[ ${sides[$s]} == *"RT"* ]]; then

            echo "It's a right side bundle, excluding left hemi WM" >> ${sub_log}

            hemi_exclude="WMhemi_LTr.nii.gz"

        else

            echo "It's a commissural bundle, no extra hemispheric excludes" >> ${sub_log}

            hemi_exclude=""

        fi

        includes_str=$(printf " -include ${s_ROIs_dir}/%s"  "${includes[@]}")

        
        # if the hemi_exclude var is empty, don't use it!

        if [[ -z ${hemi_exclude} ]]; then

            hemi_exclude_str=""

        else

            hemi_exclude_str=$(printf " -exclude ${s_ROIs_dir}/%s"  "${hemi_exclude}")

            hemi_exc_srch=($(find ${s_tck_output_dir}/${tck_name}_ROIs -type f | grep "${hemi_exclude}"))

            if [[ ! ${hemi_exc_srch} ]]; then
            
                task_in="cp ${s_ROIs_dir}/${hemi_exclude} ${s_tck_output_dir}/${tck_name}_ROIs/${tck_name}_exclude_${hemi_exclude}"

                task_exec

            else
            
                echo "${hemi_exclude} already copied" >> ${sub_log}

            fi

        fi

        if [[ ${#excludes[@]} -lt 1 ]]; then

            excludes_str=""

        else

            excludes_str=$(printf " -exclude ${s_ROIs_dir}/%s"  "${excludes[@]}")

        fi

        task_in="tckedit -nthreads ${ncpu} -force -mask ${mask1} -maxlength ${maxL} -minlength ${minL} ${wb_tck} ${includes_str} ${excludes_str} ${hemi_exclude_str} ${tck_init}"

        task_exec

        task_in="mrview -load ${FA_nii} -interpolation 1 -tractography.load ${tck_init} -focus false \
        -intensity 0,1 -fullscreen -noannotations -capture.folder ${s_tck_grabs_dir}/${tck_name} -capture.prefix ${tck_name}_init_ \
        -fov 280  -tractography.thickness 0.2 -tractography.opacity 0.6 -imagevisible 0 -mode 3 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab \
        -imagevisible 1 -fov 200 -tractography.thickness 0.2 -tractography.opacity 0.6 -mode 4 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab -exit"

        task_exec

        tck_map_1="${tck_name}_map_FD.nii.gz"

        tck_map_1_bin="${tck_name}_map_init_bin.nii.gz"

        tck_ends="${tck_name}_map_EO.nii.gz"

        task_in="tckmap -precise -force -nthreads ${ncpu} -template ${T1_brain_nii} ${tck_init} ${s_tck_dir}/${tck_name}_filtering/${tck_map_1}"

        task_exec

        task_in="mrcalc -nthreads ${ncpu} -force ${s_tck_dir}/${tck_name}_filtering/${tck_map_1} 0 -gt ${s_tck_map_out_dir}/${tck_name}/${tck_map_1_bin}"

        task_exec

        task_in="tckmap -ends_only -force -nthreads ${ncpu} -template ${T1_brain_nii} ${tck_init} ${s_tck_dir}/${tck_name}_filtering/${tck_ends}"

        task_exec

        declare -a incs_ends_only

        for inc in ${!includes[@]}; do

            inc_4FILT="${includes[$inc]}"
            
            inc_str_4FILT="${inc_4FILT::${#inc_4FILT}-10}_${tck_name}_end_only.nii.gz"

            task_in="mrcalc -force -nthreads 8 ${s_tck_dir}/${tck_name}_filtering/${tck_ends} ${s_ROIs_dir}/${includes[$inc]} \
            -mult 0 -gt - | maskfilter - dilate -force -nthreads 6 -npass 4 ${s_tck_dir}/${tck_name}_filtering/${inc_str_4FILT}"

            task_exec

            incs_ends_only[$inc]="${inc_str_4FILT}"

            cp_inc_srch=($(find ${s_tck_output_dir}/${tck_name}_ROIs -type f | grep "${tck_name}_include_${includes[$inc]}"))

            if [[ ! ${cp_inc_srch} ]]; then
            
                task_in="cp ${s_ROIs_dir}/${includes[$inc]} ${s_tck_output_dir}/${tck_name}_ROIs/${tck_name}_include_${includes[$inc]}"

                task_exec

            else
            
                echo "${includes[$inc]} already copied" >> ${sub_log}

            fi

        done

        incs_ends_only_str=$(printf " -include ${s_tck_dir}/${tck_name}_filtering/%s"  "${incs_ends_only[@]}")


        for exc in ${!excludes[@]}; do

            cp_exc_srch=($(find ${s_tck_output_dir}/${tck_name}_ROIs -type f | grep "${tck_name}_exclude_${excludes[$exc]}"))

            if [[ ! ${cp_exc_srch} ]]; then
            
                task_in="cp ${s_ROIs_dir}/${excludes[$exc]} ${s_tck_output_dir}/${tck_name}_ROIs/${tck_name}_exclude_${excludes[$exc]}"

                task_exec

            else
            
                echo "${excludes[$exc]} already copied" >> ${sub_log}

            fi
        
        done

        # count=$(tckinfo ${tck_init} -count | grep count | head -n 1 | awk '{print $(NF)}')

        b=($(fslstats ${s_tck_dir}/${tck_name}_filtering/${tck_map_1} -R))

        echo ${b[@]} >> ${sub_log}

        a=($( echo ${b[1]} | cut -c1-6))

        echo ${a[@]} >> ${sub_log}

        # setting a relative thr for tract filtering using 3% as thr cutoff for prob tracking
        # setting a relative thr for tract filtering using 1% as thr cutoff for det tracking

        # should add an if loop here, if the wb_tck is probabilistic then we do the thr based on FD, if it's deterministic then do it based on the count
        
        # might be better to define the thr within the bc command as a variable

        count=($(tckstats -force -nthreads ${ncpu} -output count ${tck_init} -quiet ))

        echo ${count} >> ${sub_log}

        if [[ ${count} -ge 100 ]]; then 

            thr_prc="5"

        else

            thr_prc="2"

        fi

        echo "for this bundle ${tck_name}, we used a thr of ${thr_prc} for tractogram filtering, this was set based on streamlines count of the initial bundle ${count}" >> ${sub_log}

        thr="$(bc <<< "scale = 5; ((${a}*("${thr_prc}"/100)))")"; echo "thr for ${tck_map_1} is ${thr}" >> ${sub_log}

        tck_out_map_filt="${s_tck_dir}/${tck_name}_filter_mask.nii.gz"

        # filtering dir, specific for each subject, each tract

        local temp1="${s_tck_dir}/${tck_name}_filtering/${tck_name}_map_thr.nii.gz"

        local temp2="${s_tck_dir}/${tck_name}_filtering/${tck_name}_map_thr_s3.nii.gz"

        local temp3="${s_tck_dir}/${tck_name}_filtering/${tck_name}_map_thr_s3_thr.nii.gz"

        local temp4="${s_tck_dir}/${tck_name}_filtering/${tck_name}_map_thr_s3_thr_GLC.nii.gz"

        local temp4h="${s_tck_dir}/${tck_name}_filtering/${tck_name}_filter_mask_rh.nii.gz"

        local temp4_inv="${s_tck_dir}/${tck_name}_filtering/${tck_name}_filter_mask_r_inv.nii.gz"

        local sum_all_incs="${s_tck_dir}/${tck_name}_filtering/${tck_name}_all_includes.nii.gz"

        local sum_all_incs_smoothed="${s_tck_dir}/${tck_name}_filtering/${tck_name}_all_includes_s3.nii.gz"

        # thr -> smooth -> thr again -> GLC

        task_in="mrthreshold -abs ${thr} -nthreads ${ncpu} -force ${s_tck_dir}/${tck_name}_filtering/${tck_map_1} ${temp1} ; mrfilter ${temp1} smooth -fwhm 3 -force -nthreads ${ncpu} ${temp2} \
        ; mrcalc ${temp2} 0.15 -ge ${temp3} -force -nthreads ${ncpu} ; ImageMath 3 ${temp4} GetLargestComponent ${temp3}"

        task_exec

        echo "Include ROIs for tckedit are : "  >> ${sub_log}

        echo ${#includes[@]} >> ${sub_log}

        echo "FT filtering ROIs for tckedit are : "  >> ${sub_log}

        echo ${#includes_filt[@]} >> ${sub_log}

        # adding all includes together to make one mask

        if [[ ${#includes_filt[@]} -gt 2 ]]; then

            echo "it's more than 2 includes " >> ${sub_log}

            rois_2_add1=$(printf " ${s_ROIs_dir}/%s"  "${includes_filt[0]}"  "${includes_filt[1]} -add")

            rois_2_add2=$(printf " ${s_ROIs_dir}/%s -add" "${includes_filt[@]:2}")

            all_rois_2_add="${rois_2_add1}  ${rois_2_add2}"

            echo ${all_rois_2_add}  >> ${sub_log}

        elif [[ ${#includes_filt[@]} -le 2 ]]; then 

            echo "it's not more than 2 includes " >> ${sub_log}

            all_rois_2_add=$(printf " ${s_ROIs_dir}/%s"  "${includes_filt[0]}"  "${includes_filt[1]} -add")

            echo ${all_rois_2_add}  >> ${sub_log}

        fi

        task_in="mrcalc -force -nthreads ${ncpu} ${all_rois_2_add} ${temp4} -add ${temp4h} ; mrcalc ${temp4h} 0 -gt ${sum_all_incs} -force -nthreads ${ncpu} \
        ; mrfilter ${sum_all_incs} smooth -fwhm 3 -force -nthreads ${ncpu} ${sum_all_incs_smoothed}"

        task_exec

        # we will try making an inverted version of the tract filter mask

        task_in="mrcalc -force -nthreads ${ncpu} ${sum_all_incs} -neg 1 -add ${temp4_inv}"

        task_exec

        tck_filt_EO1="${s_tck_dir}/${tck_name}_filt_EO.tck"

        tck_filt_EO2="${s_tck_dir}/${tck_name}_filt_EO_wm_masked.tck"

        ## edited for final output and reporting

        tck_filt1="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD.tck"

        tck_filt2="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_WM_masked.tck"

        tck_filt3="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_GMWMincs.tck"

        tck_map_out1="${s_tck_map_out_dir}/${tck_name}/${tck_name}_FD_filt_map_bin.nii.gz"

        tck_map_out2="${s_tck_map_out_dir}/${tck_name}/${tck_name}_FD_filt_WMmasked_map_bin.nii.gz"

        tck_map_out3="${s_tck_map_out_dir}/${tck_name}/${tck_name}_FD_filt_GMWMincs_map_bin.nii.gz"

        map_1_3_diff="${s_tck_map_out_dir}/${tck_name}/${tck_name}_tckmap_GMo_GMWM_filt_diff.nii.gz"

        map_1_init_diff="${s_tck_map_out_dir}/${tck_name}/${tck_name}_tckmap_init_GMfilt_diff.nii.gz"

        TCK_dump_1="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_Lengths_dump.txt"

        TCK_dump_2="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_WMmasked_Lengths_dump.txt"

        TCK_report_1="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_Lengths_report.txt"

        TCK_report_2="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_WMmasked_Lengths_report.txt"

        TCK_histo_1="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_Lengths_histo.txt"

        TCK_histo_2="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_WMmasked_Lengths_histo.txt"

        ####

        includes_filt_str=$(printf " -include ${s_ROIs_dir}/%s"  "${includes_filt[@]}")

        # using GM labels only for initial filtering

        task_in="tckedit -nthreads ${ncpu} -force -mask ${sum_all_incs_smoothed} -maxlength ${maxL} -minlength ${minL} \
        ${wb_tck} ${includes_str} -exclude ${temp4_inv} ${excludes_str} ${hemi_exclude_str} ${tck_filt1}"

        task_exec

        task_in="tckinfo -force -nthreads ${ncpu} ${tck_filt1} -quiet -force 2>&1 >> ${TCK_report_1} ; \
        echo "" >> ${TCK_report_1} ; echo ${tck_name}_FD_filt_stats: >>  ${TCK_report_1} ; \
        tckstats -force -nthreads ${ncpu} -ignorezero ${tck_filt1} -histogram ${TCK_histo_1} -dump ${TCK_dump_1} \
        -explicit -quiet -force 2>&1 >> ${TCK_report_1}"

        task_exec

        task_in="tckedit -nthreads ${ncpu} -force -mask ${sum_all_incs_smoothed} -maxlength ${maxL} -minlength ${minL} \
        ${wb_tck} ${includes_filt_str} -exclude ${temp4_inv} ${excludes_str} ${hemi_exclude_str} ${tck_filt3}"

        task_exec

        task_in="tckmap -force -nthreads ${ncpu} -template ${FA_nii} ${tck_filt1} - | mrcalc - -force -nthreads ${ncpu} 0 -gt ${tck_map_out1}"

        task_exec

        task_in="tckmap -force -nthreads ${ncpu} -template ${FA_nii} ${tck_filt3} - | mrcalc - -force -nthreads ${ncpu} 0 -gt ${tck_map_out3}"

        task_exec

        # calculate difference between the two filtered masks (GM and GMWM)

        task_in="mrcalc -force -nthreads ${ncpu}  ${tck_map_out3}  ${tck_map_out1} -subtract 0 -gt  ${map_1_3_diff}"

        task_exec

        # calculate difference between GM filtered and initial

        task_in="mrcalc -force -nthreads ${ncpu}  ${s_tck_map_out_dir}/${tck_name}/${tck_map_1_bin}  ${tck_map_out1} -subtract 0 -gt  ${map_1_init_diff}"

        task_exec

        # Generating screen grabs with MRView

        task_in="mrview -load ${FA_nii} -interpolation 1 -tractography.load ${tck_filt1} -focus false \
        -intensity 0,1 -fullscreen -noannotations -capture.folder ${s_tck_grabs_dir}/${tck_name} -capture.prefix ${tck_name}_FD_filt_p1_ \
        -fov 280  -tractography.thickness 0.2 -tractography.opacity 0.6 -imagevisible 0 -mode 3 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab \
        -imagevisible 1 -fov 200 -tractography.thickness 0.2 -tractography.opacity 0.6 -mode 4 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab -exit"

        task_exec

        task_in="tckedit -nthreads ${ncpu} -force -mask ${T1_wm_mask_inFA} \
        -maxlength ${maxL} -minlength ${minL} ${tck_filt1} -exclude ${temp4_inv} ${excludes_str} ${hemi_exclude_str} ${tck_filt2}"

        task_exec

        task_in="tckinfo -force -nthreads ${ncpu} ${tck_filt2} -quiet -force 2>&1 >> ${TCK_report_2} ; \
        echo "" >> ${TCK_report_2} ; echo ${tck_name}_FD_filt_WMmasked_stats: >>  ${TCK_report_2} ; \
        tckstats -force -nthreads ${ncpu} -ignorezero ${tck_filt2} -histogram ${TCK_histo_2} -dump ${TCK_dump_2} \
        -explicit -quiet -force 2>&1 >> ${TCK_report_2}"

        task_exec

        task_in="tckmap -force -nthreads ${ncpu} -template ${FA_nii} ${tck_filt2} - | mrcalc - -force -nthreads ${ncpu} 0 -gt ${tck_map_out2}"

        task_exec

        task_in="mrview -load ${FA_nii} -interpolation 1 -tractography.load ${tck_filt2} -focus false \
        -intensity 0,1 -fullscreen -noannotations -capture.folder ${s_tck_grabs_dir}/${tck_name} -capture.prefix ${tck_name}_FD_filt_wm_masked_p1_ \
        -fov 280  -tractography.thickness 0.2 -tractography.opacity 0.6 -imagevisible 0 -mode 3 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab \
        -imagevisible 1 -fov 200 -tractography.thickness 0.2 -tractography.opacity 0.6 -mode 4 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab -exit"

        task_exec

        task_in="tckedit -nthreads ${ncpu} -force -mask ${sum_all_incs_smoothed} -maxlength ${maxL} -minlength ${minL} \
        ${wb_tck} ${incs_ends_only_str} -exclude ${temp4_inv} ${excludes_str} ${hemi_exclude_str} ${tck_filt_EO1}"

        task_exec

        task_in="tckedit -nthreads ${ncpu} -force -mask ${T1_wm_mask_inFA} -maxlength ${maxL} -minlength ${minL} \
        ${tck_filt_EO1} -exclude ${temp4_inv} ${excludes_str} ${hemi_exclude_str} ${tck_filt_EO2}"

        task_exec


    else

        echo "${tck_name} already done, skipping to next" >> ${sub_log}

        tck_init="${s_tck_output_dir}/${tck_name}/${tck_name}_initial.tck"

        tck_filt1="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD.tck"

        tck_filt2="${s_tck_output_dir}/${tck_name}/${tck_name}_filt_FD_WM_masked.tck"

        # for the initial version

        srch_SG_init=($(find ${s_tck_grabs_dir}/${tck_name} -type f | grep "${tck_name}_init_"))

        # if loops below generate screenshots if not done already, in case the bundles were already segmented

        if [[ ! ${srch_SG_init} ]]; then

            task_in="mrview -load ${FA_nii} -interpolation 1 -tractography.load ${tck_init} -focus false \
            -intensity 0,1 -fullscreen -noannotations -capture.folder ${s_tck_grabs_dir}/${tck_name} -capture.prefix ${tck_name}_init_ \
            -fov 280  -tractography.thickness 0.2 -tractography.opacity 0.6 -imagevisible 0 -mode 3 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab \
            -imagevisible 1 -fov 200 -tractography.thickness 0.2 -tractography.opacity 0.6 -mode 4 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab -exit"

            task_exec

        else

            echo "${srch_SG_init} already generated" >> ${sub_log}

        fi

        # for the filt_1 version

        srch_SG_filt_1=($(find ${s_tck_grabs_dir}/${tck_name} -type f | grep "${tck_name}_FD_filt_p1_"))

        if [[ ! ${srch_SG_filt_1} ]]; then

            task_in="mrview -load ${FA_nii} -interpolation 1 -tractography.load ${tck_filt1} -focus false \
            -intensity 0,1 -fullscreen -noannotations -capture.folder ${s_tck_grabs_dir}/${tck_name} -capture.prefix ${tck_name}_FD_filt_ \
            -fov 280  -tractography.thickness 0.2 -tractography.opacity 0.6 -imagevisible 0 -mode 3 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab \
            -imagevisible 1 -fov 200 -tractography.thickness 0.2 -tractography.opacity 0.6 -mode 4 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab -exit"

            task_exec

        else

            echo "${srch_SG_filt_1} already generated" >> ${sub_log}

        fi

        # for the filt_2 version

        srch_SG_filt_2=($(find ${s_tck_grabs_dir}/${tck_name} -type f | grep "${tck_name}_FD_filt_wm_masked_p1_"))

        if [[ ! ${srch_SG_filt_1} ]]; then

            task_in="mrview -load ${FA_nii} -interpolation 1 -tractography.load ${tck_filt2} -focus false \
            -intensity 0,1 -fullscreen -noannotations -capture.folder ${s_tck_grabs_dir}/${tck_name} -capture.prefix ${tck_name}_FD_filt_wm_masked_ \
            -fov 280  -tractography.thickness 0.2 -tractography.opacity 0.6 -imagevisible 0 -mode 3 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab \
            -imagevisible 1 -fov 200 -tractography.thickness 0.2 -tractography.opacity 0.6 -mode 4 -plane 0 -capture.grab -plane 1 -capture.grab -plane 2 -capture.grab -exit"

            task_exec

        else

            echo "${srch_SG_filt_1} already generated" >> ${sub_log}

        fi

    fi

}


# Let's do this, for loop for patients follows



for sub in ${!subjs[@]}; do

    # for sub in 0; do

    ## make the log file

    sub_log="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_${d}_tck_fun_log.txt"

    s_ROIs_dir="${ROIs_dir}/${subjs[$sub]}"

    s_tck_dir="${Seg_tck_dir}/${subjs[$sub]}"

    s_tck_output_dir="${tck_out_dir}/${subjs[$sub]}"

    s_tck_map_out_dir="${tckmap_out_dir}/${subjs[$sub]}"

    s_tck_grabs_dir="${scr_grab_out_dir}/${subjs[$sub]}"

    if [[ ! -f ${sub_log} ]]; then

        touch ${sub_log}

    else

        echo "${sub_log} already created" >> ${sub_log}

    fi

    mkdir -p ${s_ROIs_dir}

    mkdir -p ${s_tck_dir}

    # added for reporting results and uploading

    mkdir -p ${s_tck_output_dir}

    mkdir -p ${s_tck_map_out_dir}

    mkdir -p ${s_tck_grabs_dir}

    ####

    echo "" >> ${sub_log}; echo "Now working on ${subjs[$sub]}" >> ${sub_log}

    echo "" >> ${sub_log}; echo ${work_dir} >> ${sub_log}; echo "" >> ${sub_log}

    #####

    wb_tck="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_sift_tracking-probabilistic.tck"

    # wb_tck="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_tracking-deterministic.tck"

    UKBB_2_subj="${work_dir}/${subjs[$sub]}/UKBB_2_brain"

    FS_dir="${work_dir}/FS_parcs"

    aparc_mgz="${FS_dir}/${subjs[$sub]}/mri/aparc+aseg.mgz" # needs changing

    wmparc_mgz="${FS_dir}/${subjs[$sub]}/mri/wmparc.mgz"

    aseg_mgz="${FS_dir}/${subjs[$sub]}/mri/aseg.mgz"
    
    lobes_mgz="${FS_dir}/${subjs[$sub]}/mri/FS_lobes.mgz"

    Fornix_aseg="${FS_dir}/${subjs[$sub]}/mri/FS_Fx_aseg.mgz"

    wm_seg_mgz="${FS_dir}/${subjs[$sub]}/mri/wm.seg.mgz"

    T1_nii="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_T1.nii.gz"

    T1_brain_mask="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_mask-brain.nii.gz"

    T1_brain_nii="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_T1_brain.nii.gz"

    FA_nii="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_fa.nii.gz"

    FA_nii_brain_mask="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_fa_brain_mask.nii.gz"

    sub_FA2T1="${s_ROIs_dir}/${subjs[$sub]}_fa2T1_"

    sub_T1inFA="${s_ROIs_dir}/${subjs[$sub]}_T1inFA_Warped.nii.gz"

    sub_FAinT1="${s_ROIs_dir}/${subjs[$sub]}_FAinT1_Warped.nii.gz"

    # T1_csf_mask="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_mask-csf.nii.gz"

    # T1_gm_mask="${work_dir}/${subjs[$sub]}/${subjs[$sub]}_mask-gm.nii.gz"

    # need to use WM from FS recon-all instead of the one provided with the data

    # T1_wm_mask="${work_dir}/${subjs[$sub]}/mask-wm.nii.gz"

    T1_wm_mask_mgz="${s_ROIs_dir}/${subjs[$sub]}_FS_wm_mask.mgz"

    T1_wm_mask_nii="${s_ROIs_dir}/${subjs[$sub]}_FS_wm_mask.nii.gz"

    T1_wm_mask_inFA_int="${s_ROIs_dir}/${subjs[$sub]}_T1_wm_mask_inFA_int.nii.gz"

    T1_wm_mask_inFA="${s_ROIs_dir}/${subjs[$sub]}_T1_wm_mask_inFA.nii.gz"

    # make T1 brain image

    task_in="ImageMath 3 ${T1_brain_nii} m ${T1_nii} ${T1_brain_mask}"

    task_exec

    # declare -a CST_FS_LsP

    task_in="fslmaths ${FA_nii} -bin ${FA_nii_brain_mask}"

    task_exec

    # need to register the T1s to the FA map

    srch_lobes=($(find ${FS_dir}/${subjs[$sub]}/mri -type f | grep "FS_lobes.mgz"))

    if [[ ! ${srch_lobes}  ]]; then

        task_in="mri_annotation2label --subject ${subjs[$sub]} --sd ${FS_dir} --hemi lh --lobesStrict lobes"

        task_exec

    	task_in="mri_annotation2label --subject ${subjs[$sub]} --sd ${FS_dir} --hemi rh --lobesStrict lobes"

        task_exec   

        task_in="mri_aparc2aseg --s ${subjs[$sub]} --sd ${FS_dir}  --rip-unknown  --volmask --o ${lobes_mgz}  --annot lobes --labelwm  --hypo-as-wm"

        task_exec

    else

        echo "FS lobes parcellation already created, skipping "  >> ${sub_log}

    fi

    srch_Fx_aseg=($(find ${FS_dir}/${subjs[$sub]}/mri -type f | grep "FS_Fx_aseg.mgz"))

    if [[ ! ${srch_Fx_aseg}  ]]; then

        # this command generates a new CC and Fornix segmentation using the aseg.auto.noCCseg.mgz label file

        task_in="mri_cc -aseg aseg.auto_noCCseg.mgz -o FS_Fx_aseg.mgz -sdir ${FS_dir} -f -force ${subjs[$sub]}"

        task_exec

    else

        echo "FS lobes parcellation already created, skipping "  >> ${sub_log}

    fi

    srch_FA2T1=($(find ${s_ROIs_dir} -type f | grep "${sub_FA2T1}Warped.nii.gz"));

    if [[ ! ${srch_FA2T1} ]]; then

        # warping using MI metric and Affine + mattes metric and SyN 
        # task_in="antsRegistration -d 3 -o [${sub_FA2T1},${sub_FA2T1}Warped.nii.gz,${sub_FA2T1}InverseWarped.nii.gz] \
        # -x [${T1_brain_mask},${FA_nii_brain_mask},NULL] \
        # -m MI[${T1_brain_nii},${FA_nii},1,64,Regular,0.50] \
        # -c [1000x500x250x0,1e-7,5] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 -u 1 -v 1 \
        # -m mattes[${T1_brain_nii},${FA_nii},1,64,Regular,0.75] -c [200x100x50,1e-7,5] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0mm -u 1 -z 1 --winsorize-image-intensities [0.005, 0.995]"

        # warping with Affine only
        task_in="antsRegistrationSyN.sh -d 3 -n 8 -f ${T1_brain_nii} -m ${FA_nii} -t a -x ${T1_brain_mask},${FA_nii_brain_mask} -o ${sub_FA2T1}"

        task_exec

    else

        echo "ANTs intermodality registration of FA to T1 already done, skipping " >> ${sub_log}

    fi

    # make the WM mask from FS recon-all output
    # wm.aseg is already containing all WM voxels of interest, no need to sum up labels, simply convert that and use it.

    task_in="mri_binarize --i ${wm_seg_mgz}  --min 0.1 --o ${T1_wm_mask_mgz} ; mri_convert -rl ${T1_brain_nii} -rt interpolate ${T1_wm_mask_mgz} -o ${T1_wm_mask_nii}"

    task_exec

    # Affine registration of T1 brain and WM mask to dMRI space

    task_in="WarpImageMultiTransform 3 ${T1_brain_nii} ${sub_T1inFA} -R ${FA_nii} -i ${sub_FA2T1}0GenericAffine.mat"

    task_exec

    task_in="WarpImageMultiTransform 3 ${T1_wm_mask_nii} ${T1_wm_mask_inFA_int} -R ${FA_nii} -i ${sub_FA2T1}0GenericAffine.mat"

    task_exec

    task_in="mrcalc -force -nthreads ${ncpu} ${T1_wm_mask_inFA_int} 0.15 -gt ${T1_wm_mask_inFA}"

    task_exec

    # isolate the labels with a for loop

    for L in ${!FS_Ls[@]}; do

        srch_ROIs=($(find ${s_ROIs_dir} -type f | grep "${FS_Ls[$L]}r.nii.gz"));

        if [[ ! ${srch_ROIs} ]] ; then

            # decide which volumes label to use as source

            if  [[ ${FS_Ls[$L]} == *"_GM_"* ]]  ; then

                echo "it's a grey matter label, using ${aparc_mgz}" >> ${sub_log}

                input=${aparc_mgz}

            elif [[ ${FS_Ls[$L]} == *"_WM_"* ]]; then

                echo "it's a white matter label, using ${wmparc_mgz}" >> ${sub_log}

                input=${wmparc_mgz}

            elif [[ ${FS_Ls[$L]} == *"hemi"* ]]; then

                echo "it's a hemispheric label, using ${aseg_mgz}" >> ${sub_log}

                input=${aseg_mgz}

            elif [[ ${FS_Ls[$L]} == *"lobe"* ]]; then

                echo "it's a lobe label, using ${lobes_mgz}" >> ${sub_log}

                input=${lobes_mgz}

            elif [[ ${FS_Ls[$L]} == *"Fornix"* ]]; then

                echo "it's the Fornix label, using ${Fornix_aseg}" >> ${sub_log}

                input=${Fornix_aseg}

            else

                echo "it's ${FS_Ls[$L]} using ${aparc_mgz}" >> ${sub_log}

                input=${aparc_mgz}

            fi

            min=${FS_vals[$L]}

            max=${min}

            output="${s_ROIs_dir}/${FS_Ls[$L]}.mgz"

            # isolate rois

            KUL_mri_thr

            ref=${T1_brain_nii}

            interp="interpolate"

            mgz=${output}

            nii="${output::${#output}-4}.nii.gz"

            # convert ROIs from mgz to nii.gz

            KUL_mri_mgz2nii

            nii_r="${output::${#output}-4}r.nii.gz"

            # non-linear warps applied to labels from T1 to FA space
            # task_in="WarpImageMultiTransform 3 ${nii} ${nii_r} -R ${FA_nii} -i ${sub_FA2T1}0GenericAffine.mat ${sub_FA2T1}1InverseWarp.nii.gz"

            # trying Affine only
            task_in="WarpImageMultiTransform 3 ${nii} ${nii_r} -R ${FA_nii} -i ${sub_FA2T1}0GenericAffine.mat"

            task_exec

            # ROIs_nii_preped[$L]=${nii}


        else 

            echo "FS label ${FS_Ls[$L]} already prepped" >> ${sub_log}

        fi

    done

    # get UKBB and other labels to subject space

    KUL_get_UKBB_ROIs

    ## here we need to make the includes and excludes and track each bundle

    # first we create cerebellar masks

    input1="${s_ROIs_dir}/Cerebellum_GM_LTr.nii.gz"; input2="${s_ROIs_dir}/Cerebellum_WM_LTr.nii.gz"; ROI_out="${s_ROIs_dir}/Cerebellum_GMWM_LTr.nii.gz"

    KUL_mk_ROIs

    input1="${s_ROIs_dir}/Cerebellum_GM_RTr.nii.gz"; input2="${s_ROIs_dir}/Cerebellum_WM_RTr.nii.gz"; ROI_out="${s_ROIs_dir}/Cerebellum_GMWM_RTr.nii.gz"

    KUL_mk_ROIs

    # defining sides

    sides=("LT"  "RT")

    # for loop to track for each side separately

    for s in ${!sides[@]}; do

        # global vars for Tracking

        mask1=${T1_brain_mask}

        maxL="200"

        minL="50"

        wb_tck=${wb_tck}

        # For CST

        # sum the M1 and S1 GM masks to use as one include for tckgen

        input1="${s_ROIs_dir}/M1_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/S1_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SM1_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # create M1 and S1 GM+WM masks and sum them up for filtering

        input1="${s_ROIs_dir}/M1_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/M1_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/M1_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/S1_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/S1_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/S1_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # Sum both S1 and M1

        input1="${s_ROIs_dir}/M1_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/S1_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SM1_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # includes will be used for tracking, includes_filt will be used for tract filterting

        includes=("SM1_GM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz")

        includes_filt=("SM1_GMWM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz")

        declare -a excludes 

        excludes=("Cerebellum_GMWM_LTr.nii.gz"  "Cerebellum_GMWM_RTr.nii.gz")

        tck_in="CSTML"

        TCK_gen_n_filt

        unset includes includes_filt excludes tck_in

        # separate CST and MLs

        # might need to also generate CST fibers to M1 and CST fibers to S1 separately

        # For CST only

        if [[ ${sides[$s]} == *"LT"* ]]; then

            echo "It's a left side bundle, excluding right hemi WM" >> ${sub_log}

            excludes=("PTX_cst_RT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_LT_thr_in_${subjs[$sub]}_r.nii.gz" "PTX_ml_RT_thr_in_${subjs[$sub]}_r.nii.gz")

            includes=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

            includes_filt=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

        elif [[ ${sides[$s]} == *"RT"* ]]; then

            echo "It's a right side bundle, excluding left hemi WM" >> ${sub_log}

            excludes=("PTX_cst_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_RT_thr_in_${subjs[$sub]}_r.nii.gz")

            includes=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

            includes_filt=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

        fi

        includes+=("SM1_GM_${sides[$s]}r.nii.gz")

        includes_filt+=("SM1_GMWM_${sides[$s]}r.nii.gz")

        tck_in="CST"

        TCK_gen_n_filt

        unset includes includes_filt excludes tck_in

        # for ML only

        if [[ ${sides[$s]} == *"LT"* ]]; then

            echo "It's a left side bundle, excluding right hemi WM" >> ${sub_log}

            excludes=("PTX_cst_RT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_RT_thr_in_${subjs[$sub]}_r.nii.gz" "PTX_cst_LT_thr_in_${subjs[$sub]}_r.nii.gz" \
            "Cerebellum_GMWM_LTr.nii.gz"  "Cerebellum_GMWM_RTr.nii.gz")

            includes=("PTX_ml_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

            includes_filt=("PTX_ml_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

        elif [[ ${sides[$s]} == *"RT"* ]]; then

            echo "It's a right side bundle, excluding left hemi WM" >> ${sub_log}

            excludes=("PTX_cst_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_cst_RT_thr_in_${subjs[$sub]}_r.nii.gz" \
            "Cerebellum_WM_LTr.nii.gz"  "Cerebellum_WM_RTr.nii.gz")

            includes=("PTX_ml_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

            includes_filt=("PTX_ml_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

        fi
        
        includes+=("S1_GM_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz")
        
        includes_filt+=("S1_GM_${sides[$s]}r.nii.gz")

        tck_in="ML"

        TCK_gen_n_filt

        unset includes includes_filt excludes tck_in

        # For AF

        # for tracking

        # sum up the Inferior frontal gyral GM ROIs

        input1="${s_ROIs_dir}/IFG_PTr_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IFG_POp_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IFG_POpTr_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # sum up the temporal GM ROIs

        task_in="mrcalc ${s_ROIs_dir}/STG_GM_${sides[$s]}r.nii.gz  ${s_ROIs_dir}/bSTS_GM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SMG_GM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/STG_SMG_GM_${sides[$s]}r.nii.gz -force -nthreads ${ncpu}"

        task_exec        

        # for filtering

        # isolate Inferior frontal gyral GM and WM ROIs

        input1="${s_ROIs_dir}/IFG_PTr_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IFG_PTr_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IFG_PTr_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/IFG_POp_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IFG_POp_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IFG_POp_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # sum up GMWM IFG ROIs

        input1="${s_ROIs_dir}/IFG_PTr_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IFG_POp_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IFG_POpTr_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # isolate temporal GM and WM ROIs

        input1="${s_ROIs_dir}/STG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/STG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/STG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/bSTS_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/bSTS_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/bSTS_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/SMG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/SMG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SMG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Phippo_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Phippo_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Phippo_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Insula_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Insula_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Insula_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # eroding the Insula GMWM ROIs to keep fibers passing above it

        task_in="maskfilter ${s_ROIs_dir}/Insula_GMWM_${sides[$s]}r.nii.gz erode -npass 2 ${s_ROIs_dir}/Insula_GMWM_${sides[$s]}r_ero2.nii.gz -force"

        task_exec

        # sum up all temporal GM and WM ROIs

        task_in="mrcalc ${s_ROIs_dir}/STG_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/bSTS_GMWM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SMG_GMWM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/STG_SMG_GMWM_${sides[$s]}r.nii.gz -force -nthreads ${ncpu}"

        task_exec        

        # make the dilated putamen mask

        task_in="mrfilter ${s_ROIs_dir}/Put_${sides[$s]}r.nii.gz smooth -fwhm 3 -force -nthreads ${ncpu} - | mrcalc - 0.15 -ge ${s_ROIs_dir}/s3Put_${sides[$s]}r.nii.gz -force -nthreads ${ncpu}"

        task_exec

        # make CC all ROI

        task_in="mrcalc ${s_ROIs_dir}/CC_postr.nii.gz ${s_ROIs_dir}/CC_midpostr.nii.gz -add ${s_ROIs_dir}/CC_centralr.nii.gz -add \
        ${s_ROIs_dir}/CC_midantr.nii.gz -add ${s_ROIs_dir}/CC_antr.nii.gz -add 0 -gt ${s_ROIs_dir}/CC_allr.nii.gz -force -nthreads ${ncpu}"

        task_exec

        # define tracking and filtering variables

        includes=("IFG_POpTr_GM_${sides[$s]}r.nii.gz"  "STG_SMG_GM_${sides[$s]}r.nii.gz")

        includes_filt=("IFG_POpTr_GMWM_${sides[$s]}r.nii.gz"  "STG_SMG_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "s3Put_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz"  \
        "Phippo_GMWM_${sides[$s]}r.nii.gz"  "Insula_GMWM_${sides[$s]}r_ero2.nii.gz")

        tck_in="AF"

        TCK_gen_n_filt

        unset includes includes_filt excludes tck_in

        # For SLF1 - connecting the SPL to the SFG

        # need to add posterior cingulate cortex as exclude

        input1="${s_ROIs_dir}/SFG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/SFG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SFG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/SPL_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/SPL_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SPL_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/PCC_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/PCC_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/PCC_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        includes=("SFG_GM_${sides[$s]}r.nii.gz"  "SPL_GM_${sides[$s]}r.nii.gz")

        includes_filt=("SFG_GMWM_${sides[$s]}r.nii.gz"  "SPL_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "s3Put_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz"  \
        "Phippo_GMWM_${sides[$s]}r.nii.gz" "Insula_GMWM_${sides[$s]}r_ero2.nii.gz" "Temp_lobeWM_${sides[$s]}r.nii.gz" "PCC_GMWM_${sides[$s]}r.nii.gz")

        tck_in="SLF_1"

        TCK_gen_n_filt

        unset includes includes_filt excludes tck_in

        # For SLF2 - connecting the IPL/STG to the MFG/IFG

        # for tracking

        # sum the caudal and rostral MFG GM ROIs

        input1="${s_ROIs_dir}/cMFG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/rMFG_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/crMFG_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # sum the IPL and SMG GM ROIs

        input1="${s_ROIs_dir}/IPL_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/SMG_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IPL_SMG_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # for filtering

        # isolate and sum up the Inferior parietal lobe labels

        input1="${s_ROIs_dir}/IPL_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IPL_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IPL_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/IPL_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/SMG_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IPL_SMG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # isolate and sum up the middle frontal lobe labels

        input1="${s_ROIs_dir}/cMFG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/cMFG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/cMFG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/rMFG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/rMFG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/rMFG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/cMFG_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/rMFG_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/crMFG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        includes=("crMFG_GM_${sides[$s]}r.nii.gz"  "IPL_SMG_GM_${sides[$s]}r.nii.gz")

        includes_filt=("crMFG_GMWM_${sides[$s]}r.nii.gz"  "IPL_SMG_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "s3Put_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz"  \
        "Phippo_GMWM_${sides[$s]}r.nii.gz" "Insula_GMWM_${sides[$s]}r_ero2.nii.gz"  "IFG_POpTr_GMWM_${sides[$s]}r.nii.gz")

        tck_in="SLF_2"

        TCK_gen_n_filt

        unset tck_in includes includes_filt excludes

        # For SLF3 - connecting the angular gyrus to IFG

        includes=("IFG_POpTr_GM_${sides[$s]}r.nii.gz"  "IPL_SMG_GM_${sides[$s]}r.nii.gz")

        includes_filt=("IFG_POpTr_GMWM_${sides[$s]}r.nii.gz"  "IPL_SMG_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "s3Put_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz"  \
        "Phippo_GMWM_${sides[$s]}r.nii.gz" "Insula_GMWM_${sides[$s]}r_ero2.nii.gz")

        tck_in="SLF_3"

        TCK_gen_n_filt

        # for Temporal Cingulum 

        input1="${s_ROIs_dir}/iPCC_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/iPCC_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/iPCC_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/PCC_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/PCC_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/PCC_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Fusi_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Fusi_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Fusi_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # includes=("Hippo_${sides[$s]}r.nii.gz"  "Phippo_GM_${sides[$s]}r.nii.gz"  "iPCC_GM_${sides[$s]}r.nii.gz")

        includes=("Hippo_${sides[$s]}r.nii.gz"  "Phippo_GMWM_${sides[$s]}r.nii.gz"  "iPCC_GMWM_${sides[$s]}r.nii.gz")

        includes_filt=("Hippo_${sides[$s]}r.nii.gz"  "Phippo_GMWM_${sides[$s]}r.nii.gz"  "iPCC_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "PCC_GMWM_${sides[$s]}r.nii.gz"  "Unseg_WM_${sides[$s]}r.nii.gz"  "Fusi_GMWM_${sides[$s]}r.nii.gz")

        tck_in="Cing_temp"

        TCK_gen_n_filt

        # for Dorsal Cingulum

        input1="${s_ROIs_dir}/rACC_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/rACC_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/rACC_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/cACC_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/cACC_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/cACC_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # includes=("cACC_GM_${sides[$s]}r.nii.gz"  "PCC_GM_${sides[$s]}r.nii.gz"  "iPCC_GM_${sides[$s]}r.nii.gz")

        includes=("cACC_GMWM_${sides[$s]}r.nii.gz"  "PCC_GMWM_${sides[$s]}r.nii.gz"  "iPCC_GMWM_${sides[$s]}r.nii.gz")

        includes_filt=("cACC_GMWM_${sides[$s]}r.nii.gz"  "PCC_GMWM_${sides[$s]}r.nii.gz"  "iPCC_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "STG_GMWM_${sides[$s]}r.nii.gz"  "Unseg_WM_${sides[$s]}r.nii.gz"  "rACC_GMWM_${sides[$s]}r.nii.gz")

        tck_in="Cing_cing"

        TCK_gen_n_filt

        # for subgenual Cingulum

        input1="${s_ROIs_dir}/MedOF_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/MedOF_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/MedOF_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        # includes=("cACC_GM_${sides[$s]}r.nii.gz"  "rACC_GM_${sides[$s]}r.nii.gz"  "MedOF_GM_${sides[$s]}r.nii.gz")

        includes=("cACC_GMWM_${sides[$s]}r.nii.gz"  "rACC_GMWM_${sides[$s]}r.nii.gz"  "MedOF_GMWM_${sides[$s]}r.nii.gz")

        includes_filt=("cACC_GMWM_${sides[$s]}r.nii.gz"  "rACC_GMWM_${sides[$s]}r.nii.gz"  "MedOF_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "PCC_GMWM_${sides[$s]}r.nii.gz"  "Unseg_WM_${sides[$s]}r.nii.gz"  "iPCC_GMWM_${sides[$s]}r.nii.gz"  "STG_GMWM_${sides[$s]}r.nii.gz")

        tck_in="Cing_subgenu"

        TCK_gen_n_filt

        # for ILF

        input1="${s_ROIs_dir}/TempP_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/TempP_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/TempP_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs        

        input1="${s_ROIs_dir}/ITG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/ITG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/ITG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/MTG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/MTG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/MTG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/MTG_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/ITG_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/MTG_ITG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/MTG_GM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/ITG_GM_${sides[$s]}r.nii.gz -add \
        ${s_ROIs_dir}/TempP_GM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/MITG_TP_GM_${sides[$s]}r.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/MTG_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/ITG_GMWM_${sides[$s]}r.nii.gz -add \
        ${s_ROIs_dir}/TempP_GMWM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/MITG_TP_GMWM_${sides[$s]}r.nii.gz"

        task_exec

        input1="${s_ROIs_dir}/Occ_lobeGM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Occ_lobeWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Occ_lobeGMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        includes=("Occ_lobeGM_${sides[$s]}r.nii.gz"  "MITG_TP_GM_${sides[$s]}r.nii.gz")

        includes_filt=("Occ_lobeGMWM_${sides[$s]}r.nii.gz"  "MITG_TP_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "Unseg_WM_${sides[$s]}r.nii.gz"  "Front_lobeGM_${sides[$s]}r.nii.gz"  "STG_GM_${sides[$s]}r.nii.gz")

        tck_in="ILF_trial"

        TCK_gen_n_filt

        # for MLF

        # try with precuneus, cuneus, LOCC, STG, Tempolar pole, angular (IPL), SPL (GM only!) and don't exclude MTG or ITG.

        input1="${s_ROIs_dir}/Cuneus_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Cuneus_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Cuneus_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/LatOcc_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/LatOcc_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/LatOcc_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Pc_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Pc_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Pc_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/STG_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/TempP_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/STG_TP_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/STG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/TempP_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/STG_TP_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/SPL_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IPL_GMwM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SPL_IPL_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs
        
        input1="${s_ROIs_dir}/SPL_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IPL_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SPL_IPL_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Cing_lobeGM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Cing_lobeWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Cing_lobeGMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        task_in="mrcalc ${s_ROIs_dir}/Cuneus_GM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/Pc_GM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SPL_IPL_GM_${sides[$s]}r.nii.gz -add \
        ${s_ROIs_dir}/LatOcc_GM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/PO_MLF_incs_GM_${sides[$s]}r.nii.gz -force -nthreads ${ncpu}"

        task_exec        

        task_in="mrcalc ${s_ROIs_dir}/Cuneus_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/Pc_GMWM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SPL_IPL_GMWM_${sides[$s]}r.nii.gz -add \
        ${s_ROIs_dir}/LatOcc_GMWM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/PO_MLF_incs_GMWM_${sides[$s]}r.nii.gz -force -nthreads ${ncpu}"

        task_exec        

        includes=("PO_MLF_incs_GM_${sides[$s]}r.nii.gz"  "STG_TP_GM_${sides[$s]}r.nii.gz")

        includes_filt=("SPL_IPL_GMWM_${sides[$s]}r.nii.gz"  "STG_MTG_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "Unseg_WM_${sides[$s]}r.nii.gz"  "Front_lobeWM_${sides[$s]}r.nii.gz"  "Cing_lobeGMWM_${sides[$s]}r.nii.gz"  \
        "Phippo_GMWM_${sides[$s]}r.nii.gz"  "Insula_GMWM_${sides[$s]}r_ero2.nii.gz"  "s3Put_${sides[$s]}r.nii.gz")

        tck_in="MLF_trial"

        TCK_gen_n_filt

        # for IFOF

        input1="${s_ROIs_dir}/LatOF_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/LatOF_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/LatOF_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs        

        input1="${s_ROIs_dir}/MedOF_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/LatOF_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/MedLatOF_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/MedOF_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/LatOF_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/MedLatOF_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/FrontP_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/FrontP_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/FrontP_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/Pc_GM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/Occ_lobeGM_${sides[$s]}r.nii.gz \
        -add ${s_ROIs_dir}/Fusi_GM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SPL_IPL_GM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/Par_Occ_incs_IFOF_GM_${sides[$s]}r.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/Pc_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/Occ_lobeGMWM_${sides[$s]}r.nii.gz \
        -add ${s_ROIs_dir}/Fusi_GMWM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SPL_IPL_GMWM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/Par_Occ_incs_IFOF_GMWM_${sides[$s]}r.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/MedLatOF_GM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/IFG_POpTr_GM_${sides[$s]}r.nii.gz \
        -add ${s_ROIs_dir}/FrontP_GM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/IFG_POpTr_FrontP_MLOF_GM_${sides[$s]}r.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/MedLatOF_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/IFG_POpTr_GMWM_${sides[$s]}r.nii.gz \
        -add ${s_ROIs_dir}/FrontP_GMWM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/IFG_POpTr_FrontP_MLOF_GMWM_${sides[$s]}r.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/MedLatOF_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/FrontP_GMWM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/FrontP_MLOF_GMWM_${sides[$s]}r.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/MedLatOF_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/FrontP_GMWM_${sides[$s]}r.nii.gz -add 0 -gt ${s_ROIs_dir}/FrontP_MLOF_GMWM_${sides[$s]}r.nii.gz"

        task_exec

        includes=("Par_Occ_incs_IFOF_GM_${sides[$s]}r.nii.gz"  "MedLatOF_GM_${sides[$s]}r.nii.gz")

        includes_filt=("Par_Occ_incs_IFOF_GMWM_${sides[$s]}r.nii.gz"  "MedLatOF_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "Cing_lobeGMWM_${sides[$s]}r.nii.gz"  "M1_GMWM_${sides[$s]}r.nii.gz"  "TempP_GMWM_${sides[$s]}r.nii.gz"  "ITG_GMWM_${sides[$s]}r.nii.gz"  "Acc_${sides[$s]}r.nii.gz" \
        "IFOF_exclude_midline_thr_in_${subjs[$sub]}_r.nii.gz"  "Pall_${sides[$s]}r.nii.gz"  "Amyg_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz"  "Hippo_${sides[$s]}r.nii.gz"  "vDC_${sides[$s]}r.nii.gz")

        tck_in="IFOF"

        TCK_gen_n_filt

        # for UF

        # test it out with SFG and MFG

        input1="${s_ROIs_dir}/IFG_POpTr_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/MedLatOF_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IFG_POpTr_MedLatOF_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/IFG_POpTr_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/MedLatOF_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IFG_POpTr_MedLatOF_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/IFG_POpTr_MedLatOF_GM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/crMFG_GM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SFG_GM_${sides[$s]}r.nii.gz  -add \
        0 -gt ${s_ROIs_dir}/MLOFIMSFG_GM_${sides[$s]}r.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/IFG_POpTr_MedLatOF_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/crMFG_GMWM_${sides[$s]}r.nii.gz -add ${s_ROIs_dir}/SFG_GMWM_${sides[$s]}r.nii.gz  -add \
        0 -gt ${s_ROIs_dir}/MLOFIMSFG_GMWM_${sides[$s]}r.nii.gz"

        task_exec

        includes=("TempP_GM_${sides[$s]}r.nii.gz"  "IFG_POpTr_MedLatOF_GM_${sides[$s]}r.nii.gz")

        includes_filt=("TempP_GMWM_${sides[$s]}r.nii.gz"  "IFG_POpTr_MedLatOF_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "Unseg_WM_${sides[$s]}r.nii.gz"  "Acc_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz"  "Caud_${sides[$s]}r.nii.gz"  "rACC_GMWM_${sides[$s]}r.nii.gz")

        tck_in="UF"

        TCK_gen_n_filt

        unset includes includes_filt excludes tck_in

        includes=("TempP_GM_${sides[$s]}r.nii.gz"  "MLOFIMSFG_GM_${sides[$s]}r.nii.gz")

        includes_filt=("TempP_GMWM_${sides[$s]}r.nii.gz"  "MLOFIMSFG_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "Unseg_WM_${sides[$s]}r.nii.gz"  "Acc_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz"  "Caud_${sides[$s]}r.nii.gz"  "rACC_GMWM_${sides[$s]}r.nii.gz")

        tck_in="UF_trial"

        TCK_gen_n_filt

        # for OR

        input1="${s_ROIs_dir}/PeriCalc_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/PeriCalc_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/PeriCalc_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Front_lobeGM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Front_lobeWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Front_lobeGMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        task_in="mrcalc -force -nthreads ${ncpu} ${s_ROIs_dir}/iPCC_GMWM_${sides[$s]}r.nii.gz ${s_ROIs_dir}/Thal_${sides[$s]}r.nii.gz -subtract 0 -gt ${s_ROIs_dir}/iPCC_GMWM_min_Thal_${sides[$s]}r.nii.gz"

        task_exec

        includes=("Occ_lobeGM_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz") 

        includes_filt=("Occ_lobeGMWM_${sides[$s]}r.nii.gz"  "Thal_${sides[$s]}r.nii.gz") 

        excludes=("CC_allr.nii.gz"  "vDC_${sides[$s]}r.nii.gz"  "BStemr.nii.gz"  "Phippo_GMWM_${sides[$s]}r.nii.gz"  \
        "Hippo_${sides[$s]}r.nii.gz"  "Front_lobeGMWM_${sides[$s]}r.nii.gz"  "iPCC_GMWM_min_Thal_${sides[$s]}r.nii.gz")

        if [[ ${sides[$s]} == *"LT"* ]]; then

            echo "It's the left OR, excluding right Thalamus" >> ${sub_log}

            excludes+=("Thal_RTr.nii.gz")

        elif [[ ${sides[$s]} == *"RT"* ]]; then

            echo "It's the right OR, excluding left Thalamus" >> ${sub_log}

            excludes+=("Thal_LTr.nii.gz")

        fi

        tck_in="OR"

        TCK_gen_n_filt

        unset includes_filt includes excludes tck_in

        # for Frontal Aslant Tract

        input1="${s_ROIs_dir}/IFG_POr_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/IFG_POr_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/IFG_POr_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        includes=("SFG_GM_${sides[$s]}r.nii.gz"  "IFG_POp_GM_${sides[$s]}r.nii.gz")

        includes_filt=("SFG_GMWM_${sides[$s]}r.nii.gz"  "IFG_POp_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "BStemr.nii.gz"  "Insula_GMWM_${sides[$s]}r_ero2.nii.gz"  \
        "MedLatOF_GMWM_${sides[$s]}r.nii.gz"  "IFG_POr_GMWM_${sides[$s]}r.nii.gz"  "Caud_${sides[$s]}r.nii.gz"  "FAT_exclude_midline_thr_in_${subjs[$sub]}_r.nii.gz")

        tck_in="FAT"

        TCK_gen_n_filt

        unset tck_in includes includes_filt excludes

        # for PoPT

        input1="${s_ROIs_dir}/Pari_lobeGM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Pari_lobeWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Pari_lobeGMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Pari_lobeGM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Occ_lobeGM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Pari_Occ_lobeGM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Pari_lobeGMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Occ_lobeGMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Pari_Occ_lobeGMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        input1="${s_ROIs_dir}/Temp_lobeGM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Temp_lobeWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Temp_lobeGMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs        

        if [[ ${sides[$s]} == *"LT"* ]]; then

            echo "It's a left side bundle, excluding right hemi WM" >> ${sub_log}

            excludes=("PTX_cst_RT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_LT_thr_in_${subjs[$sub]}_r.nii.gz" "PTX_ml_RT_thr_in_${subjs[$sub]}_r.nii.gz")

            includes=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

            includes_filt=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

        elif [[ ${sides[$s]} == *"RT"* ]]; then

            echo "It's a right side bundle, excluding left hemi WM" >> ${sub_log}

            excludes=("PTX_cst_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_RT_thr_in_${subjs[$sub]}_r.nii.gz")

            includes=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

            includes_filt=("PTX_cst_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz")

        fi

        includes+=("Pari_Occ_lobeGM_${sides[$s]}r.nii.gz")

        includes_filt+=("Pari_Occ_lobeGMWM_${sides[$s]}r.nii.gz")

        excludes+=("CC_allr.nii.gz"  "Front_lobeGMWM_${sides[$s]}r.nii.gz"  "Temp_lobeGMWM_${sides[$s]}r.nii.gz"  \
        "Cerebellum_GMWM_LTr.nii.gz"  "Cerebellum_GMWM_RTr.nii.gz")

        tck_in="POPT"

        TCK_gen_n_filt

        unset tck_in includes includes_filt excludes

        # # for TIF

        # includes=("Insula_GM_${sides[$s]}r.nii.gz"  "TempP_GM_${sides[$s]}r.nii.gz")

        # includes_filt=("Insula_GMWM_${sides[$s]}r.nii.gz"  "TempP_GMWM_${sides[$s]}r.nii.gz")

        # excludes=("CC_allr.nii.gz"  "Front_lobeGMWM_${sides[$s]}r.nii.gz" )

        # tck_in="TIF"

        # TCK_gen_n_filt

        # unset tck_in includes includes_filt excludes

        # for Fornix - still needs some fine tuning

        task_in="maskfilter ${s_ROIs_dir}/Thal_${sides[$s]}r.nii.gz erode -npass 4 ${s_ROIs_dir}/Thal_${sides[$s]}r_ero4.nii.gz -force"

        task_exec

        input1="${s_ROIs_dir}/Phippo_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Phippo_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Phippo_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs

        includes=("Fornixr.nii.gz"  "Hippo_${sides[$s]}r.nii.gz")

        includes_filt=("Fornixr.nii.gz"  "Hippo_${sides[$s]}r.nii.gz")

        excludes=("Amyg_${sides[$s]}r.nii.gz"  "CC_midantr.nii.gz"  "Phippo_GM_${sides[$s]}r.nii.gz" \
        "s3Put_${sides[$s]}r.nii.gz"  "BStemr.nii.gz"  "Thal_${sides[$s]}r_ero4.nii.gz"  "Front_lobeGM_${sides[$s]}r.nii.gz")

        tck_in="Fx"

        TCK_gen_n_filt

        unset tck_in includes includes_filt excludes

        # includes=("Fx_subAC_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz"  "Fx_midline_thr_in_${subjs[$sub]}_r.nii.gz")

        # includes_filt=("Fx_subAC_${sides[$s]}_thr_in_${subjs[$sub]}_r.nii.gz"  "Fx_midline_thr_in_${subjs[$sub]}_r.nii.gz")

        # excludes=("s3Put_${sides[$s]}r.nii.gz"  "Acc_${sides[$s]}r.nii.gz"  "Amyg_${sides[$s]}r.nii.gz"  "Phippo_GMWM_${sides[$s]}r.nii.gz"  \
        # "Thal_${sides[$s]}r_ero4.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  "CC_centralr.nii.gz"  "CC_midantr.nii.gz"  "BStemr.nii.gz"  "Temp_lobeGMWM_${sides[$s]}r.nii.gz")

        # tck_in="Fx_2"

        # TCK_gen_n_filt

        # will use manual ROIs for this one , brought over to subject space from the UKBB template

        # for the Thalamic radiations

        # Frontal thalamic radiation

        # includes=("Thal_${sides[$s]}r.nii.gz"  "Front_lobeGM_${sides[$s]}r.nii.gz")

        # includes_filt=("Thal_${sides[$s]}r.nii.gz"  "Front_lobeGMWM_${sides[$s]}r.nii.gz")

        # excludes=("Temp_lobeGMWM_${sides[$s]}r.nii.gz"  "Pari_lobeGMWM_${sides[$s]}r.nii.gz"  \
        # "Occ_lobeGMWM_${sides[$s]}r.nii.gz"   "BStemr.nii.gz"   "vDC_${sides[$s]}r.nii.gz")

        # tck_in="Front_Thal_rad"

        # TCK_gen_n_filt

        # unset includes includes_filt excludes

        # Anterior thalamic radiation

        # input1="${s_ROIs_dir}/MedLatOF_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/SFG_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SFG_MLOF_GM_${sides[$s]}r.nii.gz"

        # KUL_mk_ROIs

        # input1="${s_ROIs_dir}/MedLatOF_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/SFG_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/SFG_MLOF_GMWM_${sides[$s]}r.nii.gz"

        # KUL_mk_ROIs

        # includes=("Thal_${sides[$s]}r.nii.gz"  "SFG_MLOF_GM_${sides[$s]}r.nii.gz")

        # includes_filt=("Thal_${sides[$s]}r.nii.gz"  "SFG_MLOF_GMWM_${sides[$s]}r.nii.gz")

        # excludes=("Temp_lobeGMWM_${sides[$s]}r.nii.gz"  "Pari_lobeGMWM_${sides[$s]}r.nii.gz"  \
        # "Occ_lobeGMWM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz"  "vDC_${sides[$s]}r.nii.gz")

        # tck_in="Ant_Thal_rad"

        # TCK_gen_n_filt

        # unset includes includes_filt excludes

        # sensory thalamic radiation

        # includes=("Thal_${sides[$s]}r.nii.gz"  "S1_GM_${sides[$s]}r.nii.gz")

        # includes_filt=("Thal_${sides[$s]}r.nii.gz"  "S1_GMWM_${sides[$s]}r.nii.gz")

        # excludes=("CC_allr.nii.gz"  "Occ_lobeGMWM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz"  "vDC_${sides[$s]}r.nii.gz")

        # tck_in="S1_Thal_rad"

        # TCK_gen_n_filt

        # unset includes includes_filt excludes

        # # motor thalamic radiation

        # includes=("Thal_${sides[$s]}r.nii.gz"  "M1_GM_${sides[$s]}r.nii.gz")

        # includes_filt=("Thal_${sides[$s]}r.nii.gz"  "M1_GMWM_${sides[$s]}r.nii.gz")

        # excludes=("Temp_lobeGMWM_${sides[$s]}r.nii.gz"  "CC_allr.nii.gz"  "Pari_lobeGMWM_${sides[$s]}r.nii.gz" \
        # "Occ_lobeGMWM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz"  "vDC_${sides[$s]}r.nii.gz"  )

        # tck_in="M1_Thal_rad"

        # TCK_gen_n_filt

        # unset includes includes_filt excludes

        # temporal thalamic radiation

        # includes=("Thal_${sides[$s]}r.nii.gz"  "Temp_lobeGM_${sides[$s]}r.nii.gz")

        # includes_filt=("Thal_${sides[$s]}r.nii.gz"  "Temp_lobeGMWM_${sides[$s]}r.nii.gz")

        # excludes=("CC_allr.nii.gz"  "Pari_lobeGMWM_${sides[$s]}r.nii.gz" \
        # "Occ_lobeGMWM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz"  "vDC_${sides[$s]}r.nii.gz"  )

        # tck_in="Temporal_Thal_rad"

        # TCK_gen_n_filt

        # unset includes includes_filt excludes

        # auditory thalamic radiation

        input1="${s_ROIs_dir}/TTG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/TTG_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/TTG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs        

        input1="${s_ROIs_dir}/TTG_GMWM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/STG_GMWM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/STG_TTG_GMWM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs     

        input1="${s_ROIs_dir}/TTG_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/STG_GM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/STG_TTG_GM_${sides[$s]}r.nii.gz"

        KUL_mk_ROIs     

        includes=("Thal_${sides[$s]}r.nii.gz"  "STG_TTG_GM_${sides[$s]}r.nii.gz")

        includes_filt=("Thal_${sides[$s]}r.nii.gz"  "STG_TTG_GMWM_${sides[$s]}r.nii.gz")

        excludes=("CC_allr.nii.gz"  "Occ_lobeGM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz")

        tck_in="Audio_Thal_rad"

        TCK_gen_n_filt

        # unset includes includes_filt excludes

        # Parietal thalamic radiation

        # includes=("Thal_${sides[$s]}r.nii.gz"  "Pari_lobeGM_${sides[$s]}r.nii.gz")

        # includes_filt=("Thal_${sides[$s]}r.nii.gz"  "Pari_lobeGMWM_${sides[$s]}r.nii.gz")

        # excludes=("CC_allr.nii.gz"  "BStemr.nii.gz")

        # tck_in="Parietal_Thal_rad"

        # TCK_gen_n_filt

        # unset includes includes_filt excludes

        # for the cerebellar pathways

        # Superior cerebellar peduncle

        # input1="${s_ROIs_dir}/Cerebellum_GM_${sides[$s]}r.nii.gz"; input2="${s_ROIs_dir}/Cerebellum_WM_${sides[$s]}r.nii.gz"; ROI_out="${s_ROIs_dir}/Cerebellum_GMWM_${sides[$s]}r.nii.gz"

        # KUL_mk_ROIs

        # includes=("Thal_${sides[$s]}r.nii.gz"  "Pari_lobeGM_${sides[$s]}r.nii.gz")

        # includes_filt=("Thal_${sides[$s]}r.nii.gz"  "Pari_lobeGMWM_${sides[$s]}r.nii.gz")

        # excludes=("CC_allr.nii.gz"  "Front_lobeGMWM_${sides[$s]}r.nii.gz"  "Temp_lobeGMWM_${sides[$s]}r.nii.gz" \
        # "Occ_lobeGMWM_${sides[$s]}r.nii.gz"  "BStemr.nii.gz"  "vDC_${sides[$s]}r.nii.gz"  )

        # tck_in="SCP"

        # TCK_gen_n_filt

        # unset includes includes_filt excludes

        # Inferior cerebellar peduncle



    done

    # commissural bundles must be done separately.

    # for CC

    # Track fibers to each lobe separately :)

    unset sides

    mask1=${T1_brain_mask}

    maxL="200"

    minL="30"

    wb_tck=${wb_tck}

    # CC frontal fibers

    includes=("CC_allr.nii.gz"  "Front_lobeGM_LTr.nii.gz"  "Front_lobeGM_RTr.nii.gz")

    includes_filt=("CC_allr.nii.gz"  "Front_lobeGMWM_LTr.nii.gz"  "Front_lobeGMWM_RTr.nii.gz")

    excludes=("Temp_lobeGM_RTr.nii.gz"  "Temp_lobeGM_LTr.nii.gz"  "Pari_lobeGM_LTr.nii.gz"  "Pari_lobeGM_RTr.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  \
    "Occ_lobeGM_LTr.nii.gz"  "Occ_lobeGM_RTr.nii.gz"  "BStemr.nii.gz"  "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz"  "Fornixr.nii.gz" )

    tck_in="CC_Frontal"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    # CC fibers to the temporal lobes without PHippo

    includes=("CC_allr.nii.gz"  "Temp_lobeGM_LTr.nii.gz"  "Temp_lobeGM_RTr.nii.gz")

    includes_filt=("CC_allr.nii.gz"  "Temp_lobeGMWM_LTr.nii.gz"  "Temp_lobeGMWM_RTr.nii.gz")

    excludes=("Front_lobeGM_RTr.nii.gz"  "Front_lobeGM_LTr.nii.gz"  "BStemr.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  \
    "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz"  "Phippo_GM_LTr.nii.gz"  "Phippo_GM_RTr.nii.gz")

    tck_in="CC_Temporal"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    # CC temporal fibers with PHippo

    includes=("CC_allr.nii.gz"  "Temp_lobeGM_LTr.nii.gz"  "Temp_lobeGM_RTr.nii.gz")

    includes_filt=("CC_allr.nii.gz"  "Temp_lobeGMWM_LTr.nii.gz"  "Temp_lobeGMWM_RTr.nii.gz")

    excludes=("Front_lobeGM_RTr.nii.gz"  "Front_lobeGM_LTr.nii.gz"  "BStemr.nii.gz"  \
    "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz")

    tck_in="CC_Temporal_wPhippo"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    ## CC parietal fibers

    includes=("CC_allr.nii.gz"  "Pari_lobeGM_LTr.nii.gz"  "Pari_lobeGM_RTr.nii.gz")

    includes_filt=("CC_allr.nii.gz"  "Pari_lobeGMWM_LTr.nii.gz"  "Pari_lobeGMWM_RTr.nii.gz")

    excludes=("AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  "Temp_lobeGM_RTr.nii.gz"  "Temp_lobeGM_LTr.nii.gz"  "BStemr.nii.gz"  "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz")

    tck_in="CC_Parietal"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    # CC fibers to the occipital lobes

    includes=("CC_allr.nii.gz"  "Occ_lobeGM_LTr.nii.gz"  "Occ_lobeGM_RTr.nii.gz")

    includes_filt=("CC_allr.nii.gz"  "Occ_lobeGMWM_LTr.nii.gz"  "Occ_lobeGMWM_RTr.nii.gz")

    excludes=("Front_lobeGM_RTr.nii.gz"  "Front_lobeGM_LTr.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  \
    "Temp_lobeGM_LTr.nii.gz"  "Temp_lobeGM_RTr.nii.gz"  "BStemr.nii.gz"  "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz")

    tck_in="CC_Occipital"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    # For MCPs and Transverse Pontine fibers

    # for the cerebellar pathways

    # includes=("Cerebellum_GM_LTr.nii.gz"  "Cerebellum_GM_LTr.nii.gz"  "BStemr.nii.gz")

    # includes_filt=("Cerebellum_GMWM_LTr.nii.gz"  "Cerebellum_GMWM_LTr.nii.gz"  "BStemr.nii.gz")

    # excludes=("CC_allr.nii.gz"  "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz")

    # tck_in="MCP"

    # TCK_gen_n_filt

    # unset includes includes_filt excludes

    # CC fibers to the cingulate lobes (not running this)

    # includes=("CC_allr.nii.gz"  "Cing_lobeGM_LTr.nii.gz"  "Cing_lobeGM_RTr.nii.gz")

    # includes_filt=("CC_allr.nii.gz"  "Cing_lobeGMWM_LTr.nii.gz"  "Cing_lobeGMWM_RTr.nii.gz")

    # excludes=("BStemr.nii.gz"  "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz")

    # tck_in="CC_Cingulate_inc"

    # TCK_gen_n_filt

    # unset includes includes_filt excludes

    # for AC

    includes=("Temp_lobeGM_LTr.nii.gz"  "Temp_lobeGM_RTr.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz")

    includes_filt=("Temp_lobeGMWM_LTr.nii.gz"  "Temp_lobeGMWM_RTr.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz")

    excludes=("CC_allr.nii.gz"  "BStemr.nii.gz")

    tck_in="AC_temp"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    includes=("Temp_lobeGM_LTr.nii.gz"  "Temp_lobeGM_RTr.nii.gz"  "Put_LTr.nii.gz"  "Put_RTr.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz")

    includes_filt=("Put_LTr.nii.gz"  "Put_RTr.nii.gz"   "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz")

    excludes=("CC_allr.nii.gz"  "BStemr.nii.gz")

    tck_in="AC_put_temp"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    # for PC
    # use the PC ROI defined manually + the parietal lobes
    # add occipital lobes to the includes of the PC and remove th vDC as an exclude, too restrictive
    # try UKBB labels as excludes

    input1="${s_ROIs_dir}/Temp_lobeGM_RTr.nii.gz"; input2="${s_ROIs_dir}/Occ_lobeGM_RTr.nii.gz"; ROI_out="${s_ROIs_dir}/Temp_Occ_GM_RTr.nii.gz"

    KUL_mk_ROIs        

    input1="${s_ROIs_dir}/Temp_lobeGM_LTr.nii.gz"; input2="${s_ROIs_dir}/Occ_lobeGM_LTr.nii.gz"; ROI_out="${s_ROIs_dir}/Temp_Occ_GM_LTr.nii.gz"

    KUL_mk_ROIs        

    input1="${s_ROIs_dir}/Temp_lobeGMWM_RTr.nii.gz"; input2="${s_ROIs_dir}/Occ_lobeGM_RTr.nii.gz"; ROI_out="${s_ROIs_dir}/Temp_Occ_GMWM_RTr.nii.gz"

    KUL_mk_ROIs        

    input1="${s_ROIs_dir}/Temp_lobeGMWM_LTr.nii.gz"; input2="${s_ROIs_dir}/Occ_lobeGMWM_LTr.nii.gz"; ROI_out="${s_ROIs_dir}/Temp_Occ_GMWM_LTr.nii.gz"

    KUL_mk_ROIs        

    includes=("Pari_lobeGM_LTr.nii.gz"  "Pari_lobeGM_RTr.nii.gz"  "PC_midline_thr_in_${subjs[$sub]}_r.nii.gz"   "BStemr.nii.gz")

    includes_filt=("Pari_lobeGMWM_LTr.nii.gz"  "Pari_lobeGMWM_RTr.nii.gz"  "PC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  "BStemr.nii.gz")

    # excludes=("CC_allr.nii.gz"  "vDC_LTr.nii.gz"  "vDC_RTr.nii.gz")

    excludes=("CC_allr.nii.gz"  "PTX_cst_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_cst_RT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_LT_thr_in_${subjs[$sub]}_r.nii.gz"  "PTX_ml_RT_thr_in_${subjs[$sub]}_r.nii.gz")

    tck_in="PC"

    TCK_gen_n_filt

    unset includes includes_filt excludes

    input1="${s_ROIs_dir}/Hippo_LTr.nii.gz"; input2="${s_ROIs_dir}/Hippo_RTr.nii.gz"; ROI_out="${s_ROIs_dir}/Bil_Hippos_r.nii.gz"

    KUL_mk_ROIs

    # for the Fornix bilateral version

    # includes=("Fornixr.nii.gz"  "Bil_Hippos_r.nii.gz")

    # includes_filt=("Fornixr.nii.gz"  "Bil_Hippos_r.nii.gz")

    # excludes=("Acc_LTr.nii.gz"  "Acc_RTr.nii.gz"  "Amyg_LTr.nii.gz"  "Amyg_RTr.nii.gz"  "CC_centralr.nii.gz"  "CC_midantr.nii.gz"  "Phippo_GMWM_LTr.nii.gz"  "Phippo_GMWM_RTr.nii.gz" \
    # "s3Put_LTr.nii.gz"  "s3Put_RTr.nii.gz"  "BStemr.nii.gz"  "Thal_LTr_ero4.nii.gz"   "Thal_RTr_ero4.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  "Front_lobeGMWM_LTr.nii.gz"  "Front_lobeGMWM_RTr.nii.gz")

    # tck_in="Fx_Bil_hippos"

    # TCK_gen_n_filt

    # unset tck_in includes includes_filt excludes

    # for the Fornix Right to left

    # includes=("Fornixr.nii.gz"  "Hippo_LTr.nii.gz"  "Hippo_LTr.nii.gz")

    # includes_filt=("Fornixr.nii.gz"  "Hippo_LTr.nii.gz"  "Hippo_LTr.nii.gz")

    # excludes=("Acc_LTr.nii.gz"  "Acc_RTr.nii.gz"  "Amyg_LTr.nii.gz"  "Amyg_RTr.nii.gz"  "CC_centralr.nii.gz"  "CC_midantr.nii.gz"  "Phippo_GMWM_LTr.nii.gz"  "Phippo_GMWM_RTr.nii.gz" \
    # "s3Put_LTr.nii.gz"  "s3Put_RTr.nii.gz"  "BStemr.nii.gz"  "Thal_LTr_ero4.nii.gz"   "Thal_RTr_ero4.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz"  "Front_lobeGMWM_LTr.nii.gz"  "Front_lobeGMWM_RTr.nii.gz")

    # tck_in="Fx_comm"

    # TCK_gen_n_filt

    # for the Audio commissure

    includes=("CC_allr.nii.gz"  "STG_TTG_GM_LTr.nii.gz"  "STG_TTG_GM_RTr.nii.gz")

    includes_filt=("CC_allr.nii.gz"  "STG_TTG_GMWM_LTr.nii.gz"  "STG_TTG_GMWM_RTr.nii.gz")

    excludes=("Occ_lobeGMWM_LTr.nii.gz"  "Occ_lobeGMWM_RTr.nii.gz"  "BStemr.nii.gz"  "AC_midline_thr_in_${subjs[$sub]}_r.nii.gz")

    tck_in="Audio_commissure"

    TCK_gen_n_filt

    unset includes includes_filt excludes


    # should add warping step to MNI

done