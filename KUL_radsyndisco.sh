#!/bin/bash

# @ AR 04/02/2021
# this script will apply synb0disco without docker or singularity

# to do:
# inset exec_function
# insert loggin
# insert input arg parse
# insert function path, mrtrix, FS, ANTs finders

TOPUP=1


# This script needs to know
# where ANTs, FS, FSL and Synb0DISCO live

# define some vars

radwd="$(pwd)"

# timestamp
start=$(date +%s)
d=$(date "+%Y-%m-%d_%H-%M-%S")

# logging
rsbzd_log="${radwd}/KUL_RadSynDisco_${d}.txt";

if [[ ! -f ${rsbzd_log} ]] ; then

    touch ${rsbzd_log}

else

    echo "${rsbzd_log} already created"

fi

## exec all func

function task_exec {

    echo "-------------------------------------------------------------" | tee -a ${rsbzd_log}

    echo ${task_in} | tee -a ${rsbzd_log}

    echo " Started @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${rsbzd_log}

    eval ${task_in} 2>&1 | tee -a ${rsbzd_log} &

    echo " pid = $! " | tee -a ${rsbzd_log}

    wait ${pid}

    sleep 5

    echo "exit status $?" | tee -a ${rsbzd_log}

    echo " Finished @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${rsbzd_log}

    echo "-------------------------------------------------------------" | tee -a ${rsbzd_log}

    echo "" | tee -a ${rsbzd_log}

    unset task_in

}


######

# Setup pytorch thru conda
# source /extra/pytorch/bin/activate
# only the pytorch sourcing is needed
# this can be done with conda also

a=($(which conda))
b=($(dirname ${a}))
c="${b}/../etc/profile.d/conda.sh"

sbzd_p1=($(which prepare_input_AR.sh))
sbzd_pD=($(dirname ${sbzd_p1}))
sbzd_p2="${sbzd_pD}/.."

# source conda and activate virtual env
# source ${c} # Or path to where your conda is
# conda activate radsyndisco # or conda virtual env of your choosing
# conda create --name radsyndisco; conda activate radsyndisco
# pip install torch
# pip install torchvision
# pip install numpy scipy matplotlib ipython jupyter pandas sympy nose
# pip install 'nibabel==2.5.2'


# make the output dir

######
# Start of actual script

# Prepare input

if [[ ! -f "${radwd}/OUTPUTS/b0_d_nonlin_atlas_2_5.nii.gz" ]]; then

    task_in="prepare_input_AR.sh ${radwd}/INPUTS/b0.nii.gz ${radwd}/INPUTS/T1.nii.gz ${radwd}/INPUTS/T1_mask.nii.gz ${sbzd_p2}/atlases/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz \
    ${sbzd_p2}/atlases/mni_icbm152_t1_tal_nlin_asym_09c_2_5.nii.gz ${radwd}/OUTPUTS 64"

    task_exec

fi

# exit 2

# Run inference
NUM_FOLDS=5

if [[ ! -f "${radwd}/OUTPUTS/b0_u_lin_atlas_2_5_FOLD_5.nii.gz" ]]; then

    for i in $(seq 1 ${NUM_FOLDS}); do 

        echo Performing inference on FOLD: "${i}" | tee -a ${rsbzd_log}

        # need to try this with py3.8
        # it was using python3.6 originally
        # task_in="python ${sbzd_p2}/src/inference.py ${radwd}/OUTPUTS/T1_norm_lin_atlas_2_5.nii.gz ${radwd}/OUTPUTS/b0_d_lin_atlas_2_5.nii.gz ${radwd}/OUTPUTS/b0_u_lin_atlas_2_5_FOLD_${i}.nii.gz ${sbzd_p2}/src/train_lin/num_fold_?_total_folds_${NUM_FOLDS}_seed_1_num_epochs_100_lr_0.0001_betas_*_1e-05_num_epoch_*.pth"

        python ${sbzd_p2}/src/inference.py ${radwd}/OUTPUTS/T1_norm_lin_atlas_2_5.nii.gz ${radwd}/OUTPUTS/b0_d_lin_atlas_2_5.nii.gz ${radwd}/OUTPUTS/b0_u_lin_atlas_2_5_FOLD_${i}.nii.gz ${sbzd_p2}/src/train_lin/num_fold_"${i}"_total_folds_"${NUM_FOLDS}"_seed_1_num_epochs_100_lr_0.0001_betas_\(0.9\,\ 0.999\)_weight_decay_1e-05_num_epoch_*.pth | tee -a ${rsbzd_log}


    done

    # Take mean
    echo "Taking ensemble average" | tee -a ${rsbzd_log}

    task_in="fslmerge -t ${radwd}/OUTPUTS/b0_u_lin_atlas_2_5_merged.nii.gz ${radwd}/OUTPUTS/b0_u_lin_atlas_2_5_FOLD_*.nii.gz"
    task_exec

    task_in="fslmaths ${radwd}/OUTPUTS/b0_u_lin_atlas_2_5_merged.nii.gz -Tmean ${radwd}/OUTPUTS/b0_u_lin_atlas_2_5.nii.gz"

    task_exec

fi

# Apply inverse xform to undistorted b0
echo "Applying inverse xform to undistorted b0" | tee -a ${rsbzd_log}

temp_dir="${sbzd_pD}/../atlases"

task_in="antsApplyTransforms -d 3 -i ${radwd}/OUTPUTS/b0_u_lin_atlas_2_5.nii.gz -r ${radwd}/INPUTS/b0.nii.gz -n BSpline -t [${temp_dir}/mni_2_5_2_1_0GenericAffine.mat,0] -t [${radwd}/OUTPUTS/T1_2_template_0GenericAffine.mat,1] -t [${radwd}/OUTPUTS/ANTs_II_d_0GenericAffine.mat,1] -o ${radwd}/OUTPUTS/b0_u.nii.gz"

task_exec

# Smooth image
echo "Applying slight smoothing to distorted b0" | tee -a ${rsbzd_log}

task_in="fslmaths ${radwd}/INPUTS/b0.nii.gz -s 1.15 ${radwd}/OUTPUTS/b0_d_smooth.nii.gz"

task_exec

# conda deactivate
# conda virtual env with py3.6 + numpy 1.16.4
# this is to ensure eddyquad will work
# conda activate py3.6_np1.16.4 

if [[ $TOPUP -eq 1 ]]; then
    # Merge results and run through topup

    if [[ ! -f "${radwd}/OUTPUTS/b0_all_topup.nii.gz" ]]; then
        echo Running topup | tee -a ${rsbzd_log}
        
        task_in="fslmerge -t ${radwd}/OUTPUTS/b0_all.nii.gz ${radwd}/OUTPUTS/b0_d_smooth.nii.gz ${radwd}/OUTPUTS/b0_u.nii.gz"

        task_exec

        task_in="topup -v --imain=${radwd}/OUTPUTS/b0_all.nii.gz \
            --datain=${radwd}/INPUTS/acqparams.txt --config=b02b0.cnf \
            --iout=${radwd}/OUTPUTS/b0_all_topup.nii.gz --out=${radwd}/OUTPUTS/topup \
            --fout=${radwd}/OUTPUTS/topup_fieldmap.nii.gz
            --subsamp=1,1,1,1,1,1,1,1,1 \
            --miter=10,10,10,10,10,20,20,30,30 \
            --lambda=0.00033,0.000067,0.0000067,0.000001,0.00000033,0.000000033,0.0000000033,0.000000000033,0.00000000000067 \
            --scale=0"

        task_exec

    fi

fi


# Done
echo "FINISHED!!!" | tee -a ${rsbzd_log}