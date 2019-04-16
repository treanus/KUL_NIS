cwd=$(pwd)
ncpu=32

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
dwiintensitynorm dwiintensitynorm/dwi_input/ dwiintensitynorm/mask_input/ dwiintensitynorm/dwi_output/ dwiintensitynorm/fa_template.mif dwiintensitynorm/fa_template_wm_mask.mif -nthreads $ncpu


# Adding a subject
#dwi2tensor new_subject/dwi_denoised_unringed_preproc_unbiased.mif -mask new_subject/dwi_temp_mask.mif - | tensor2metric - -fa - | mrregister -force \
# ../dwiintensitynorm/fa_template.mif - -mask2 new_subject/dwi_temp_mask.mif -nl_scale 0.5,0.75,1.0 -nl_niter 5,5,15 -nl_warp - /tmp/dummy_file.mif | \
# mrtransform ../dwiintensitynorm/fa_template_wm_mask.mif -template new_subject/dwi_denoised_unringed_preproc_unbiased.mif -warp - - | dwinormalise \ 
# new_subject/dwi_denoised_unringed_preproc_unbiased.mif - ../dwiintensitynorm/dwi_output/new_subject.mif


# Computing an (average) white matter response function
mkdir -p ${cwd}/dwiprep/group/fba/response
foreach -${ncpu} ${cwd}/dwiprep/group/fba/dwiintensitynorm/dwi_output/*.mif : dwi2response tournier IN \
${cwd}/dwiprep/group/fba/response/PRE_response.txt

average_response ${cwd}/dwiprep/group/fba/dwiintensitynorm/dwi_output/*response.txt ${cwd}/dwiprep/group/fba/group_average_response.txt

# Compute new brain mask images
mkdir -p ${cwd}/dwiprep/group/fba/mask
foreach -${ncpu} ${cwd}/dwiprep/group/fba/dwiintensitynorm/dwi_output/*.mif : dwi2mask IN \
${cwd}/dwiprep/group/fba/mask/PRE_mask.mif

# Fibre Orientation Distribution estimation (spherical deconvolution)
# see https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/st_fibre_density_cross-section.html
#  Note that dwi2fod csd can be used, however here we use dwi2fod msmt_csd (even with single shell data) to benefit from the hard non-negativity constraint, 
#  which has been observed to lead to more robust outcomes:
mkdir -p ${cwd}/dwiprep/group/fba/fod
foreach -${ncpu} ${cwd}/dwiprep/group/fba/dwiintensitynorm/dwi_output/*.mif : dwiextract IN - \
\| dwi2fod msmt_csd - ${cwd}/dwiprep/group/fba/group_average_response.txt ${cwd}/dwiprep/group/fba/fod/PRE_wmfod.mif \
-mask ${cwd}/dwiprep/group/fba/mask/PRE_mask.mif

# Generate a study-specific unbiased FOD template
mkdir -p ${cwd}/dwiprep/group/fba/template
population_template {cwd}/dwiprep/group/fba/fod -mask_dir ${cwd}/dwiprep/group/fba/mask ${cwd}/dwiprep/group/fba/template/wmfod_template.mif -voxel_size 1.3





