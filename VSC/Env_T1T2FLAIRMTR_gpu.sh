echo "Activating conda environment mri_brain"
conda activate mri_brain

# activate Pytorch
#module load PyTorch/1.0.1-intel-2018a

# setup freesurfer development version of dd 10/03/2021
echo "Activating Freesurfer7 development versions centos7 of 10/03/2021"
export FREESURFER_HOME=$VSC_DATA/apps/freesurfer
export SUBJECTS_DIR=$FREESURFER_HOME/subjects
export FS_LICENSE=$VSC_DATA/apps/freesurfer_license/license.txt
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# necessary for samseg
module load cuDNN/7.4.1-CUDA-10.0.130

# Setup FastSurfer
export FASTSURFER_HOME=$VSC_DATA/apps/FastSurfer

# load KNT
echo "Loading KUL_NeuroImaging_Tools"
export PATH=${VSC_DATA}/apps/KUL_NeuroImaging_Tools:$PATH

# load MRtrix3
echo "Loading MRtrix3"
module load MRtrix/3.0.2-foss-2018a-Python-3.6.4

# load Ants
echo "Loading Ants"
module load ANTs/2.3.1-foss-2018a-Python-2.7.14
export ANTSPATH=/apps/leuven/skylake/2018a/software/ANTs/2.3.1-foss-2018a-Python-2.7.14/bin

# Loading Pytorch overwrites my conda env path 8-(
#export PATH="/data/leuven/322/vsc32269/miniconda3/envs/mri_brain/bin":$PATH

