#!/bin/bash -l

#PBS -l nodes=1:ppn=##CPU##:gpus=##GPU##
#PBS -A ##LP##
#PBS -l walltime=##WALLTIME##
##PARTITION##

export SINGULARITY_CACHEDIR=$VSC_SCRATCH/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR
export TMPDIR=$VSC_SCRATCH/tmp
mkdir -p $TMPDIR

#-------------------------------------------------


#load modules
module purge
module load FreeSurfer/6.0.0-centos6_x86_64
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export FS_LICENSE=${FREESURFER_HOME}/license.txt
module load FSL/6.0.1-foss-2018a
. ${FSLDIR}/etc/fslconf/fsl.sh
module load MRtrix/3.0.2-foss-2018a-Python-3.6.4 
module load Python/3.6.4-foss-2018a
module load CUDA/9.1.85
module load ANTs/2.3.1-foss-2018a-Python-2.7.14 

export OMP_NUM_THREADS=##CPU##

#add the path to neuroimaging tools
PATH=${VSC_DATA}/apps/KUL_NeuroImaging_Tools:${PATH}

##COMMAND##
