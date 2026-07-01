#!/bin/bash

participant=$1
conn_dir=CONN/conn_project01/results/firstlevel/SBC_01
sources_list=${conn_dir}/_list_sources.txt
fmriprep_dir=fmriprep/sub-${participant}
transform=${fmriprep_dir}/anat/sub-${participant}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
reference=${fmriprep_dir}/anat/sub-${participant}_desc-preproc_T1w.nii.gz
ref_fMRI=CONN/conn_project01/results/firstlevel/SBC_01/BETA_Subject001_Condition001_Source001.nii
T1w_GM=${fmriprep_dir}/anat/sub-${participant}_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz
function KUL_antsApply_Transform {
        echo "input=$input"
        echo "output=$output"
        echo "transform=$transform"
        echo "reference=$reference"
    antsApplyTransforms -d 3 --float 1 \
        --verbose $ants_verbose \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n Linear
}

function average_network {
    si=''
    for s in ${sources[@]}; do
        # extract the source name and image
        search_string="Source"$s
        echo $search_string
        line=$(grep "$search_string" "$sources_list")
        if [[ $line =~ networks\.(.*)\ \((.*)\) ]]; then
            extracted_text="${BASH_REMATCH[1]}"
            text_before_dot=$(echo "$extracted_text" | cut -d'.' -f1)
            formatted_text=$(echo "$extracted_text" | tr ' ' '_' | tr -d '()' | tr '.' '_')
        fi

        si="$si ${conn_dir}/BETA_Subject001_Condition001_Source${s}.nii "
        
        # extract the individual connectivity maps
        mrcalc ${conn_dir}/BETA_Subject001_Condition001_Source${s}.nii 0.3 -gt CONN/T1w_GM_resampled.nii.gz -mult \
            CONN/$formatted_text.nii.gz -force
        input=CONN/$formatted_text.nii.gz
        output=CONN/${formatted_text}_space-subject.nii.gz
        KUL_antsApply_Transform
        KUL_mrview_figure.sh -u $reference \
            -o CONN/${formatted_text}_space-subject.nii.gz \
            -t 2 -p $participant


    done
    echo $si
    net=$text_before_dot
    mrmath $si mean - | mrcalc - 0.2 -gt CONN/T1w_GM_resampled.nii.gz -mult \
        CONN/${net}_average.nii.gz -force
    input=CONN/${net}_average.nii.gz
    output=CONN/${net}_average_space-subject.nii.gz
    KUL_antsApply_Transform

    KUL_mrview_figure.sh -u $reference \
            -o CONN/${net}_average_space-subject.nii.gz \
            -t 2 -p $participant

}

# Make the T1w_GM same dimensions as fMRI (GM-only masking — default)
mrgrid $T1w_GM -template $ref_fMRI regrid CONN/T1w_GM_resampled.nii.gz -force
# DMN
sources=(001 002 003 004)
average_network

# SM
sources=(005 006 007)
average_network

sources=(008 009 010 011)
net='VIS'
average_network

sources=(012 013 014 015 016 017 018)
net='SAL'
average_network

sources=(023 024 025 026)
net='FP'
average_network

sources=(027 028 029 030)
net='LANGUAGE'
average_network
