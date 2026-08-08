#!/bin/bash
set -x

# to do 23/03
# figure out which outputs need to go to outdir

# Function to display usage
usage() {
    echo "Usage: $0 -s <subject_name> -i <input_dir> -w <work_dir> -o <output_dir> -d <in_DSC> -t <DSCT1c> -l <template_brain>"
    echo "  -s: Subject name or ID"    
    echo "  -i: Input directory"
    echo "  -w: Working directory"
    echo "  -o: Output directory"
    echo "  -d: Input DSC file"
    echo "  -t: DSCT1c file"
    echo "  -l: (Optional) high resolution T1w post-contrast file. If not provided, a default will be used."
    exit 1
}

# Parse command-line arguments
while getopts ":s:i:w:o:d:t:l:" opt; do
    case $opt in
        s) subj="$OPTARG" ;;
        i) indir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        d) in_DSC="$OPTARG" ;;
        t) DSCT1c="$OPTARG" ;;
        l) template_brain="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check for missing arguments
if [[ -z "$subj" || -z "$indir" || -z "$workdir" || -z "$outdir" || -z "$in_DSC" || -z "$DSCT1c" ]]; then
    echo "Error: Missing arguments!"
    usage
fi

temp_flag=0

# set lesion source flag
if [[ ! -z "$template_brain" ]]; then
    temp_flag=1
fi

# Create directories
#mkdir -p "$indir"
mkdir -p "$workdir"
mkdir -p "$outdir"

# Logging
logfile="$workdir/processing_log.txt"
echo "Starting processing at $(date)" > "$logfile"

# vars
DSC_vol1="$workdir/DSC_vol1.nii.gz"
DSC_vol1_brain="$workdir/DSC_vol1_brain.nii.gz"
DSC_vol1_brain_mask="$workdir/DSC_vol1_brain_mask.nii.gz"
DSC_vol1_brain_mask_d2="$workdir/DSC_vol1_brain_mask_d2.nii.gz"

in_DSC_denoised="$workdir/in_DSC_denoised.nii.gz"
in_DSC_noise="$workdir/in_DSC_noise.nii.gz"

in_DSCT1c_denoised="$workdir/in_DSCT1c_denoised.nii.gz"
in_DSCT1c_noise="$workdir/in_DSCT1c_noise.nii.gz"
in_DSCT1c_dnbc="$workdir/in_DSCT1c_dnbc.nii.gz"
in_DSCT1c_dnbc_brain="$workdir/in_DSCT1c_dnbc_brain.nii.gz"
DSCT1c_brain_regrid="$workdir/DSCT1c_brain_regrid.nii.gz"


# if [[ ! -f $in_DSCT1c_dnbc ]]; then
#     in_DSCT1c_dnbc="$workdir/in_DSCT1c_denoised_dnbc.nii.gz"
#     in_DSCT1c_dnbc_brain="$workdir/in_DSCT1c_denoised_dnbc_brain.nii.gz"

# fi

in_DSCT1c_synthparc="$workdir/synthparc_in_DSCT1c.nii.gz"
in_DSCT1c_synthseg="$workdir/synthseg_in_DSCT1c.nii.gz"
in_DSCT1c_synthseg_rg="$workdir/synthseg_in_DSCT1c_rg.nii.gz"
in_DSCT1c_synthparc_rg="$workdir/synthparc_in_DSCT1c_rg.nii.gz"
in_DSCT1c_5tt="$workdir/in_DSCT1c_5tt.nii.gz"
Rhemi_WM="$workdir/Rhemi_WM.nii.gz"
Rhemi_GM="$workdir/Rhemi_GM.nii.gz"
Lhemi_WM="$workdir/Lhemi_WM.nii.gz"
Lhemi_GM="$workdir/Lhemi_GM.nii.gz"
Rhemi_mask="$workdir/Rhemi_mask.nii.gz"
Lhemi_mask="$workdir/Lhemi_mask.nii.gz"
Rhemi_lesion_part="$workdir/Lhemi_by_Lmask.nii.gz"
Lhemi_lesion_part="$workdir/Lhemi_by_Lmask.nii.gz"
Rhemi_WM_ero4="$workdir/Lhemi_WM_ero4.nii.gz"
Lhemi_WM_ero4="$workdir/Lhemi_WM_ero4.nii.gz"
NAWM_contralesion="$outdir/Contralesional_NAWM_mask.nii.gz"

in_DSC_denoised_avg="$workdir/in_DSC_denoised_avg.nii.gz"
in_DSC_denoised_mc="$workdir/in_DSC_denoised_mc.nii.gz"
in_DSC_denoised_mc_vol1="$workdir/in_DSC_denoised_mc_vol1.nii.gz"
DSCT1c_brain_regrid_mask="$workdir/DSCT1c_brain_regrid_mask.nii.gz"
DSC_antsII_EPICorrected="$workdir/DSC_antsII_EPICorrected.nii.gz"   
DSC_pp="$outdir/DSC_pp.nii.gz"
DSC_pp_median="$outdir/DSC_pp_median.nii.gz"
DSC_outbrainmask="$outdir/DSC_pp_brain_mask.nii.gz"

template_brain=${template_brain:-"/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"}


if [[ ! -f "$DSC_pp" ]]; then

    # Step 1: Convert and extract DSC vol1, extract brains
    echo "Step 1: Extracting DSC vol1 and brains..." | tee -a "$logfile"

    mrconvert -coord 3 0 -axes 0,1,2 "$in_DSC" "$DSC_vol1" >> "$logfile" 2>&1
    mri_synthstrip -i "$DSC_vol1" -o "$DSC_vol1_brain" -m "$DSC_vol1_brain_mask" -g >> "$logfile" 2>&1
    mri_synthstrip -i "$DSCT1c" -o "$workdir/DSCT1c_brain.nii.gz" -m "$workdir/DSCT1c_brain_mask.nii.gz" -g >> "$logfile" 2>&1
    # maskfilter -npass 2 "$DSC_vol1_brain_mask" dilate "$DSC_vol1_brain_mask_d2" >> "$logfile" 2>&1
    # ImageMath 3 "$DSC_vol1_brain_mask_d2" MD "$DSC_vol1_brain_mask" 2 >> "$logfile" 2>&1

    # Step 2: Denoise inputs
    echo "Step 2: Denoising inputs..." | tee -a "$logfile"
    
    if [[ ! -f "$in_DSCT1c_denoised" ]]; then 

        DenoiseImage -d 4 -i "$in_DSC" -o "[$in_DSC_denoised,$in_DSC_noise]" -v 1 >> "$logfile" 2>&1

        wait

        DenoiseImage -d 3 -i "$DSCT1c" -o "[$in_DSCT1c_denoised,$in_DSCT1c_noise]" -v 1 >> "$logfile" 2>&1

        wait

    fi
    # Step 3: Preprocess T1c further
    echo "Step 3: Preprocessing T1c..." | tee -a "$logfile"
    
    N4BiasFieldCorrection -d 3 -i "$in_DSCT1c_denoised" -o "$in_DSCT1c_dnbc" -x "$workdir/DSCT1c_brain_mask.nii.gz" -t "[0.15,0.01,200]" -r 1 >> "$logfile" 2>&1
    mrcalc "$in_DSCT1c_dnbc" "$workdir/DSCT1c_brain_mask.nii.gz" -mult "$in_DSCT1c_dnbc_brain" >> "$logfile" 2>&1
    aa=($(mrinfo -spacing "$DSC_vol1"))

    mrgrid -voxel "${aa[0]},${aa[1]},${aa[2]}" "$in_DSCT1c_dnbc_brain" regrid "$DSCT1c_brain_regrid" >> "$logfile" 2>&1

    # Step 4: Start registrations
    echo "Step 4: Starting registrations..." | tee -a "$logfile"
    if [[ ! -f "$workdir/DSC_LRT1c_2_HRT1c_Warped.nii.gz" ]]; then
        antsRegistrationSyN.sh -d 3 -m "$DSCT1c_brain_regrid" -f "$template_brain" -t s -o "$workdir/DSCT1c_2_temp_brain_" -n 80 >> "$logfile" 2>&1
        antsRegistrationSyN.sh -d 3 -m "$DSCT1c_brain_regrid" -f "$in_DSCT1c_dnbc_brain" -t r -o "$workdir/DSC_LRT1c_2_HRT1c_" -n 80 >> "$logfile" 2>&1
    fi
    # Step 5: Motion correct the DSC data
    echo "Step 5: Motion correcting DSC data..." | tee -a "$logfile"
    
    if [[ ! -f "$in_DSC_denoised_mc" ]]; then
        antsMotionCorr -d 3 -a "$in_DSC_denoised" -o "$in_DSC_denoised_avg" >> "$logfile" 2>&1
        antsMotionCorr -d 3 -o "[$in_DSC_denoised_mc,$in_DSC_denoised_mc,$in_DSC_denoised_avg]" \
            -m MI["$in_DSC_denoised_avg", "$in_DSC_denoised", 1, 32, Random, 0.50] -i 20 -u 1 -e 1 -n 10 -t Affine[0.005] -s 0 -f 1 >> "$logfile" 2>&1
    fi

    # Step 6: EPI distortion correction
    echo "Step 6: EPI distortion correction..." | tee -a "$logfile"
    
    mrconvert -coord 3 0 -axes 0,1,2 "$in_DSC_denoised_mc" "$in_DSC_denoised_mc_vol1" >> "$logfile" 2>&1
    mrcalc "$DSCT1c_brain_regrid" 0 -gt "$DSCT1c_brain_regrid_mask" >> "$logfile" 2>&1
    antsIntermodalityIntrasubject.sh -d 3 -i "$in_DSC_denoised_mc_vol1" -r "$DSCT1c_brain_regrid" -T "$template_brain" -x "$DSCT1c_brain_regrid_mask" \
    -w "$workdir/DSCT1c_2_temp_brain_" -o "$workdir/DSC_antsII_EPIC_" -t 3 >> "$logfile" 2>&1

    antsApplyTransforms -d 3 -e 3 -i "$in_DSC_denoised_mc" -o "$DSC_antsII_EPICorrected" -r "$DSCT1c_brain_regrid" -t "$workdir/DSC_antsII_EPIC_1Warp.nii.gz" \
        -t "$workdir/DSC_antsII_EPIC_0GenericAffine.mat" -n LanczosWindowedSinc >> "$logfile" 2>&1

    # Step 7: Bias field correction
    echo "Step 7: Bias field correction..." | tee -a "$logfile"
    # DSCT1c_brain_regrid_mask_d2="$workdir/DSCT1c_brain_regrid_mask_d2.nii.gz"

    # maskfilter -npass 2 "$DSCT1c_brain_regrid_mask" dilate "$DSCT1c_brain_regrid_mask_d2" >> "$logfile" 2>&1
    N4BiasFieldCorrection -d 4 -i "$DSC_antsII_EPICorrected" -o "$DSC_pp" -r 1 >> "$logfile" 2>&1

fi

pyDSC_out="$outdir/DSLC_fit_output_2"

if [[ ! -f "${pyDSC_out}/rCBV_corrected.nii.gz" ]]; then 

    mrmath -axis 3 "$DSC_pp" median "$DSC_pp_median" -force 
    mri_synthstrip -i "$DSC_pp_median" -o "$workdir/DSCpp_median_brain.nii.gz" -m "$workdir/DSCpp_median_brain_mask.nii.gz" -g >> "$logfile" 2>&1
    maskfilter -npass 3 "$workdir/DSCpp_median_brain_mask.nii.gz" dilate "$DSC_outbrainmask" -force

    # Step 8: Run leakage correction fitting in python
    echo "Step 8: Running leakage correction fitting in Python..." | tee -a "$logfile"
    # Call your R script here with appropriate inputs
    # we need to find N dim of D3 (time) and TE
    tsize=($(mrinfo -spacing "$in_DSC"))
    TE="0.0$(mrinfo -property comments "$in_DSC" | cut -d ";" -f1 | cut -d "=" -f2)"
    # Run the R script for fitting
    echo "Running Python script for fitting..." | tee -a "$logfile"
    # fix this to correct version and update inputs

    python ./Good_DSCLC_fit5.py --te=${TE} --tr=${tsize[3]} ${DSC_pp} $workdir/DSCpp_median_brain_mask.nii.gz ${pyDSC_out}/

fi

cp "$workdir/DSCT1c_brain_regrid.nii.gz" "$outdir/DSCT1c_brain_regrid.nii.gz"

## breaking off here for lesion fixes

# exit

if [[ ! -f "$NAWM_contralesion" ]]; then

    ## uncomment the following if statement when all lesion masks are present
    # if [[ -f "$(dirname $template_brain)/../../Masks_4_MR/sub-${subj}_AR_seg.nii.gz" ]]; then
        if [[ $temp_flag == 1 ]]; then 
            # Apply transforms to lesion_seg
            # could define lesion mask as input argument also
            if [[ -f "/mnt/DATA1/aradwa0/rad_dude/Data_w_DSC_portal_maps/niis/Extra_Lmasks/${subj}_Lmask.nii.gz" ]]; then 
                cp /mnt/DATA1/aradwa0/rad_dude/Data_w_DSC_portal_maps/niis/Extra_Lmasks/${subj}_Lmask.nii.gz \
                    $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz
            else
                antsApplyTransforms -d 3 -i $(dirname $template_brain)/../../Masks_4_MR/sub-${subj}_AR_seg.nii.gz \
                    -o $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz \
                    -r "$in_DSCT1c_dnbc_brain" -t "[$workdir/DSCT1c_2_temp_brain_0GenericAffine.mat,1]" \
                    -t "$workdir/DSCT1c_2_temp_brain_1InverseWarp.nii.gz" \
                    -t "$workdir/DSC_LRT1c_2_HRT1c_0GenericAffine.mat" -n MultiLabel
            fi    
        else
            # lesion mask is already in HR T1 space so we just sample
            echo "Lesion mask already in T1 space so we continue"
            cp /mnt/DATA1/aradwa0/rad_dude/Data_w_DSC_portal_maps/niis/Extra_Lmasks/${subj}_Lmask.nii.gz \
                $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz
        fi

        # Step 9: Register resulting maps to HR T1c
        echo "Step 9: Registering resulting maps to HR T1c..." | tee -a "$logfile"
        SC_LRT1c_2_HRT1c="$workdir/DSC_LRT1c_2_HRT1c_"
        pyrCBVcorr="$pyDSC_out/rCBV_corrected.nii.gz"
        pyrCBVcorr_inHRT1c="$outdir/pyrCBVcorr_inHRT1c.nii.gz"
        pyrCBVuncorr="$pyDSC_out/rCBV_uncorrected.nii.gz"
        pyrCBVuncorr_inHRT1c="$outdir/pyrCBVuncorr_inHRT1c.nii.gz"
        pyrCBF="$pyDSC_out/rCBF.nii.gz"
        pyrCBF_inHRT1c="$outdir/pyrCBF_inHRT1c.nii.gz"

        antsApplyTransforms -d 3 -i "$pyrCBVcorr" -o "$pyrCBVcorr_inHRT1c" -r "$in_DSCT1c_dnbc_brain" \
            -t "${SC_LRT1c_2_HRT1c}0GenericAffine.mat" -n LanczosWindowedSinc >> "$logfile" 2>&1
        antsApplyTransforms -d 3 -i "$pyrCBVuncorr" -o "$pyrCBVuncorr_inHRT1c" -r "$in_DSCT1c_dnbc_brain" \
            -t "${SC_LRT1c_2_HRT1c}0GenericAffine.mat" -n LanczosWindowedSinc >> "$logfile" 2>&1
        antsApplyTransforms -d 3 -i "$pyrCBF" -o "$pyrCBF_inHRT1c" -r "$in_DSCT1c_dnbc_brain" \
            -t "${SC_LRT1c_2_HRT1c}0GenericAffine.mat" -n LanczosWindowedSinc >> "$logfile" 2>&1

        # Step 10: Register DSC portal maps to T1c
        echo "Step 10: Registering DSC portal maps to T1c..." | tee -a "$logfile"
        in_DSCLC_portal="${indir}/sub-${subj}_T2_DSC_Perfusion_MB2_S2.2.nii.gz"
        in_DSCLC_portal_vol1="${workdir}/sub-${subj}_T2_DSC_portal_vol1.nii.gz"
        in_DSCLC_portal_vol1_brain="${workdir}/sub-${subj}_T2_DSC_portal_vol1_brain.nii.gz"
        in_DSCLC_portal_vol1_brain_mask="${workdir}/sub-${subj}_T2_DSC_portal_vol1_brain_mask.nii.gz"
        mrconvert -coord 3 0 -axes 0,1,2 ${in_DSCLC_portal} ${in_DSCLC_portal_vol1} -force
        mri_synthstrip -i "$in_DSCLC_portal_vol1" -o "$in_DSCLC_portal_vol1_brain" \
            -m "$in_DSCLC_portal_vol1_brain_mask" -g >> "$logfile" 2>&1
        antsIntermodalityIntrasubject.sh -d 3 -i $in_DSCLC_portal_vol1_brain -r "$DSCT1c_brain_regrid" -T "$template_brain" \
            -x $DSCT1c_brain_regrid_mask -w "$workdir/DSCT1c_2_temp_brain_" -o "$workdir/PortalDSC_2_LRT1c_EPIC_" -t 3

        # Repeat for each portal DSC map as needed
        ISPrCBVcorr=$(ls -f $indir/sub-${subj}_nrCBVcorr*.nii.gz)
        ISPrCBVcorr_inHRT1c="$outdir/ISPrCBVcorr_inHRT1c.nii.gz"
        ISPrCBVuncorr=$(ls $indir/sub-${subj}_nrCBVuncorr*.nii.gz)
        ISPrCBVuncorr_inHRT1c="$outdir/ISPrCBVuncorr_inHRT1c.nii.gz"
        ISPrCBF=$(ls $indir/sub-${subj}_nrelCBFAIF*.nii.gz)
        ISPrCBF_inHRT1c="$outdir/ISPrCBF_inHRT1c.nii.gz"

        antsApplyTransforms -d 3 -i "$ISPrCBVcorr" -o "$ISPrCBVcorr_inHRT1c" -r "$in_DSCT1c_dnbc_brain" -t "$workdir/PortalDSC_2_LRT1c_EPIC_1Warp.nii.gz" \
            -t "$workdir/PortalDSC_2_LRT1c_EPIC_0GenericAffine.mat" -n LanczosWindowedSinc >> "$logfile" 2>&1
        antsApplyTransforms -d 3 -i "$ISPrCBVuncorr" -o "$ISPrCBVuncorr_inHRT1c" -r "$in_DSCT1c_dnbc_brain" -t "$workdir/PortalDSC_2_LRT1c_EPIC_1Warp.nii.gz" \
            -t "$workdir/PortalDSC_2_LRT1c_EPIC_0GenericAffine.mat" -n LanczosWindowedSinc >> "$logfile" 2>&1
        antsApplyTransforms -d 3 -i "$ISPrCBF" -o "$ISPrCBF_inHRT1c" -r "$in_DSCT1c_dnbc_brain" -t "$workdir/PortalDSC_2_LRT1c_EPIC_1Warp.nii.gz" \
            -t "$workdir/PortalDSC_2_LRT1c_EPIC_0GenericAffine.mat" -n LanczosWindowedSinc >> "$logfile" 2>&1

        # Step 11: Sample the parametric maps using lesion and synthseg derived VOIs
        if [[ ! -f "$in_DSCT1c_synthparc" ]]; then 
            mri_synthseg --i "$in_DSCT1c_dnbc" --o "$in_DSCT1c_synthparc" --parc --threads 20
            mri_synthseg --i "$in_DSCT1c_dnbc" --o "$in_DSCT1c_synthseg" --threads 20
        fi
        antsApplyTransforms -d 3 -i ${in_DSCT1c_synthseg} -r "$in_DSCT1c_dnbc" \
        -o ${in_DSCT1c_synthseg_rg} -n MultiLabel

        antsApplyTransforms -d 3 -i ${in_DSCT1c_synthparc} -r "$in_DSCT1c_dnbc" \
        -o ${in_DSCT1c_synthparc_rg} -n MultiLabel

        # maybe not needed actually - similarly for the parcellation eh
        5ttgen freesurfer ${in_DSCT1c_synthparc_rg} ${in_DSCT1c_5tt} -force

        # based on the LUT of the segmentation
        # we determine relative laterality, and sample the NAWM of the contralesional side
        # or the side with less lesion overlap, with a distance determined by e.g. 4x dilated mask

        mrcalc ${in_DSCT1c_synthseg_rg} 41 -eq $Rhemi_WM -force
        mrcalc ${in_DSCT1c_synthseg_rg} 2 -eq $Lhemi_WM -force
        mrcalc ${in_DSCT1c_synthseg_rg} 42 -eq $Rhemi_GM -force
        mrcalc ${in_DSCT1c_synthseg_rg} 3 -eq $Lhemi_GM -force

        mrcalc $Lhemi_GM $Lhemi_WM -add 0 -gt $Lhemi_mask -force
        mrcalc $Rhemi_GM $Rhemi_WM -add 0 -gt $Rhemi_mask -force

        mrcalc $Lhemi_mask $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz -mult $Lhemi_lesion_part -force
        mrcalc $Rhemi_mask $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz -mult $Rhemi_lesion_part -force

        Lesion_Lvxcount=$(mrstats -output count -mask $Lhemi_mask -ignorezero $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz)
        Lesion_Rvxcount=$(mrstats -output count -mask $Rhemi_mask -ignorezero $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz)
        Lesion_totalvxcount=$(mrstats -output count -ignorezero $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz)

        # Calculate percentage of lesion in each hemisphere
        PrcntL=$(echo "scale=4; 100 * $Lesion_Lvxcount / $Lesion_totalvxcount" | bc)
        PrcntR=$(echo "scale=4; 100 * $Lesion_Rvxcount / $Lesion_totalvxcount" | bc)

        # Compare lesion percentages and define laterality
        # dilate lesion mask 4x to make sure it is not included in NAWM mask
        maskfilter -npass 4 $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz dilate \
            $workdir/sub-${subj}_Lesion_seg_inDSC_dil4.nii.gz -force

        if (( $(echo "$PrcntL > $PrcntR" | bc -l) )); then
            Lesion_laterality="Left"
            echo "Left or predominantly left sided lesion"
            maskfilter -npass 4 $Rhemi_WM erode $Rhemi_WM_ero4 -force
            mrcalc $workdir/sub-${subj}_Lesion_seg_inDSC_dil4.nii.gz 0 -eq $Rhemi_WM_ero4 -mult \
                $NAWM_contralesion -force
        else    
            Lesion_laterality="Right"
            echo "Right or predominantly right sided lesion"
            maskfilter -npass 4 $Lhemi_WM erode $Lhemi_WM_ero4 -force
            mrcalc $workdir/sub-${subj}_Lesion_seg_inDSC_dil4.nii.gz 0 -eq $Lhemi_WM_ero4 -mult \
                $NAWM_contralesion -force
        fi

fi
    # # sample lesion rCBV only in lesion mask 
    # mrcalc $outdir/sub-${subj}_Lesion_seg_inDSC.nii.gz 0 -gt - | maskfilter - dilate \
    #     $workdir/sub-${subj}_Lesion_seg_inDSCT1c_dil2.nii.gz -force

    # # is it almost truly bilateral? 
    # # Calculate absolute difference between Percent_L and Percent_R
    # diff=$(echo "$PrcntL - $PrcntR" | bc -l)
    # abs_diff=$(echo "if ($diff < 0) -1 * $diff else $diff" | bc -l)
    # threshold=10.0

    # ##
    # echo "Lesion is more dominant in the: $Lesion_laterality hemisphere"

    # echo "Left hemisphere lesion volume: ${PrcntL}%"
    # echo "Right hemisphere lesion volume: ${PrcntR}%"

    # if (( $(echo "$abs_diff < $threshold" | bc -l) )); then
    #     echo "fyi the lesion is bilateral or almost bilateral"
    # fi

    # # sample the parametric maps using the masks we made
    # # rCBV
    # L_out1=$(mrstats -mask $workdir/sub-${subj}_Lesion_seg_inDSCT1c_dil2.nii.gz  "$ISPrCBVcorr_inHRT1c")
    # L_out2=$(mrstats -mask $NAWM_contralesion  "$ISPrCBVcorr_inHRT1c" -quiet)

    # L_out3=$(mrstats -mask $workdir/sub-${subj}_Lesion_seg_inDSCT1c_dil2.nii.gz  "$ISPrCBVuncorr_inHRT1c" -quiet)
    # L_out4=$(mrstats -mask $NAWM_contralesion  "$ISPrCBVcorr_inHRT1c" -quiet)

    # L_out5=$(mrstats -mask $workdir/sub-${subj}_Lesion_seg_inDSCT1c_dil2.nii.gz  "$pyrCBVcorr_inHRT1c" -quiet)
    # L_out6=$(mrstats -mask $NAWM_contralesion  "$pyrCBVcorr_inHRT1c" -quiet)

    # L_out7=$(mrstats -mask $workdir/sub-${subj}_Lesion_seg_inDSCT1c_dil2.nii.gz  "$pyrCBVuncorr_inHRT1c" -quiet)
    # L_out8=$(mrstats -mask $NAWM_contralesion  "$pyrCBVuncorr_inHRT1c" -quiet)

    # #rCBF
    # L_out9=$(mrstats -mask $workdir/sub-${subj}_Lesion_seg_inDSCT1c_dil2.nii.gz  "$ISPrCBF_inHRT1c" -quiet)
    # L_out10=$(mrstats -mask $NAWM_contralesion  "$ISPrCBF_inHRT1c" -quiet)

    # L_out11=$(mrstats -mask $workdir/sub-${subj}_Lesion_seg_inDSCT1c_dil2.nii.gz  "$pyrCBF_inHRT1c" -quiet)
    # L_out12=$(mrstats -mask $NAWM_contralesion  "$pyrCBF_inHRT1c" -quiet)

    # echo "ISPrCBVcorr_lesion ${L_out1}" | tee -a $outdir/DSC_values.txt
    # echo "ISPrCBVcorr_NAWM ${L_out2}" | tee -a $outdir/DSC_values.txt
    # echo "ISPrCBVuncorr_lesion ${L_out3}" | tee -a $outdir/DSC_values.txt
    # echo "ISPrCBVuncorr_NAWM ${L_out4}" | tee -a $outdir/DSC_values.txt
    # echo "pyrCBVcorr_lesion ${L_out5}" | tee -a $outdir/DSC_values.txt
    # echo "pyrCBVcorr_NAWM ${L_out6}" | tee -a $outdir/DSC_values.txt
    # echo "pyrCBVuncorr_lesion ${L_out7}" | tee -a $outdir/DSC_values.txt
    # echo "pyrCBVuncorr_NAWM ${L_out8}" | tee -a $outdir/DSC_values.txt
    # echo "ISPrCBF_lesion ${L_out9}" | tee -a $outdir/DSC_values.txt
    # echo "ISPrCBF_NAWM ${L_out10}" | tee -a $outdir/DSC_values.txt
    # echo "pyrCBF_lesion ${L_out11}" | tee -a $outdir/DSC_values.txt
    # echo "pyrCBF_NAWM ${L_out12}" | tee -a $outdir/DSC_values.txt

    # echo "Processing complete! Output saved to $outdir" | tee -a "$logfile"

# fi