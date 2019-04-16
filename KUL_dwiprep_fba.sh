cwd=$(pwd)

mkdir -p dwiprep/group/fba/dwiintensitynorm/dwi_input
mkdir -p dwiprep/group/fba/dwiintensitynorm/mask_input



cd dwiprep/group/fba

# find the preproced mifs
search_sessions=($(find ${cwd}/dwiprep -type f | grep dwi_preproced_reg2T1w.mif))
num_sessions=${#search_sessions[@]}
    
echo "  Fixel based input dwi: $num_sessions"
echo "    notably: ${search_sessions[@]}"

for i in ${search_sessions[@]}
do

    s=$(echo $i | cut -d'/' -f 3)
    cp $i ${cwd}/dwiprep/group/fba/dwiintensitynorm/dwi_input/${s}_dwi_preproced_reg2T1w.mif

done

# find the preproced masks
search_sessions=($(find ${cwd}/dwiprep -type f | grep dwi_preproced_reg2T1w_mask.nii.gz))
num_sessions=${#search_sessions[@]}
    
echo "  Fixel based input mask: $num_sessions"
echo "    notably: ${search_sessions[@]}"

for i in ${search_sessions[@]}
do

    s=$(echo $i | cut -d'/' -f 3)
    mrconvert $i ${cwd}/dwiprep/group/fba/dwiintensitynorm/mask_input/${s}_dwi_preproced_reg2T1w_mask.mif

done


# Intensity Normalisation
dwiintensitynorm dwiintensitynorm/dwi_input/ dwiintensitynorm/mask_input/ dwiintensitynorm/dwi_output/ dwiintensitynorm/fa_template.mif dwiintensitynorm/fa_template_wm_mask.mif -nthreads 32


# Adding a subject
#dwi2tensor new_subject/dwi_denoised_unringed_preproc_unbiased.mif -mask new_subject/dwi_temp_mask.mif - | tensor2metric - -fa - | mrregister -force \
# ../dwiintensitynorm/fa_template.mif - -mask2 new_subject/dwi_temp_mask.mif -nl_scale 0.5,0.75,1.0 -nl_niter 5,5,15 -nl_warp - /tmp/dummy_file.mif | \
# mrtransform ../dwiintensitynorm/fa_template_wm_mask.mif -template new_subject/dwi_denoised_unringed_preproc_unbiased.mif -warp - - | dwinormalise \ 
# new_subject/dwi_denoised_unringed_preproc_unbiased.mif - ../dwiintensitynorm/dwi_output/new_subject.mif

exit 0

# Computing an (average) white matter response function
foreach * : dwi2response tournier IN/dwi_denoised_unringed_preproc_unbiased_normalised.mif IN/response.txt

average_response */response.txt ../group_average_response.txt

# Fibre Orientation Distribution estimation (spherical deconvolution)
foreach * : dwiextract IN/dwi_denoised_unringed_preproc_unbiased_normalised_upsampled.mif - \| dwi2fod msmt_csd - ../group_average_response.txt IN/wmfod.mif -mask IN/dwi_mask_upsampled.mif

# Generate a study-specific unbiased FOD template
mkdir -p ../template/fod_input
mkdir ../template/mask_input

