cwd=$(pwd)

mkdir -p dwiprep/group/fba/dwiintensitynorm/dwi_input
mkdir -p dwiprep/group/fba/dwiintensitynorm/mask_input



#cd dwiprep/group/fba/dwiintensitynorm/dwi_input/

# find the preproced mifs
search_sessions=($(find dwiprep -type f | grep dwi_preproced_reg2T1w.mif))
num_sessions=${#search_sessions[@]}
    
echo "  Fixel based input dwi: $num_sessions"
echo "    notably: ${search_sessions[@]}"

for i in ${search_sessions[@]}
do

    s=$(echo $i | cut -d'/' -f 3)
    cp $i ${cwd}/dwiprep/group/fba/dwiintensitynorm/dwi_input/${s}_dwi_preproced_reg2T1w.mif

done

# find the preproced masks
search_sessions=($(find dwiprep -type f | grep dwi_preproced_reg2T1w_mask.nii.gz))
num_sessions=${#search_sessions[@]}
    
echo "  Fixel based input mask: $num_sessions"
echo "    notably: ${search_sessions[@]}"

for i in ${search_sessions[@]}
do

    s=$(echo $i | cut -d'/' -f 3)
    mrconvert $i dwiprep/group/fba/dwiintensitynorm/mask_input/${s}_dwi_preproced_reg2T1w_mask.mif

done

exit 0 

# Intensity Normalisation
dwiintensitynorm dwiprep/group/fba/dwiintensitynorm/dwi_input/ dwiprep/group/fba/dwiintensitynorm/mask_input/ dwiprep/group/fba/dwiintensitynorm/dwi_output/ dwiprep/group/fba/dwiintensitynorm/fa_template.mif dwiprep/group/fba/dwiintensitynorm/fa_template_wm_mask.mif -nthreads 32


# Adding a subject
#dwi2tensor new_subject/dwi_denoised_unringed_preproc_unbiased.mif -mask new_subject/dwi_temp_mask.mif - | tensor2metric - -fa - | mrregister -force \
# ../dwiintensitynorm/fa_template.mif - -mask2 new_subject/dwi_temp_mask.mif -nl_scale 0.5,0.75,1.0 -nl_niter 5,5,15 -nl_warp - /tmp/dummy_file.mif | \
# mrtransform ../dwiintensitynorm/fa_template_wm_mask.mif -template new_subject/dwi_denoised_unringed_preproc_unbiased.mif -warp - - | dwinormalise \ 
# new_subject/dwi_denoised_unringed_preproc_unbiased.mif - ../dwiintensitynorm/dwi_output/new_subject.mif

