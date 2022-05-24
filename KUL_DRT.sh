#!/bin/bash
# Bash shell script to analyse DTI
#
# Requires fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 18/05/2022
version="0.1"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` is a batch analysis of DTI data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -d DICOM/JohnDoe.zip

Required arguments:

     -p:  participant name

Optional arguments:

     -n:  number of cpu to use (default 15)
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ants_verbose=1
ncpu=15
bc=0 
type=1
redo=0
results=0 
verbose_level=1

# Set required options
p_flag=0
d_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:n:v:" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
        v) #verbose
            verbose_level=$OPTARG
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
if [ $p_flag -eq 0 ] ; then
	echo
	echo "Option -p is required: give the BIDS name of the participant." >&2
	echo
	exit 2
fi

KUL_LOG_DIR="KUL_LOG/${script}/sub-${participant}"
mkdir -p $KUL_LOG_DIR

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
	export MRTRIX_QUIET=1
    silent=1
    str_silent=" > /dev/null 2>&1" 
    ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
    silent=0
    str_silent="" 
    ants_verbose=1
fi

function KUL_run_fmriprep {
    if [ ! -f fmriprep/sub-${participant}.html ]; then
        
        # preparing for fmriprep
        cp study_config/run_fmriprep.txt KUL_LOG/sub-${participant}_run_fmriprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        fmriprep_options="--fs-no-reconall --anat-only "
        sed -i.bck "s/fmriprep_options: /fmriprep_options: ${fmriprep_options}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        
        # running fmriprep
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_fmriprep.txt 
        
        # cleaning the working directory
        #rm -fr fmriprep_work_${participant}

        # copy to the RESULTS dir


    else
        echo "Fmriprep already done"
    fi
}

function KUL_run_msbp {

    if [ ! -f KUL_LOG/sub-${participant}_MSBP.done ]; then

        echo "Running MSBP"

        # there seems tpo be a problem with docker if the fsaverage dir is a soft link; so we delete the link and hardcopy it
        rm -fr $cwd/BIDS/derivatives/freesurfer/fsaverage
        cp -r $FREESURFER_HOME/subjects/fsaverage $cwd/BIDS/derivatives/freesurfer/fsaverage

        task_in="docker run --rm -u $(id -u) -v $cwd/BIDS:/bids_dir \
         -v $cwd/BIDS/derivatives:/output_dir \
         -v $FS_LICENSE:/opt/freesurfer/license.txt \
         sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
         --participant_label $participant --isotropic_resolution 1.0 --thalamic_nuclei \
         --brainstem_structures --skip_bids_validator --fs_number_of_cores $ncpu \
         --multiproc_number_of_cores $ncpu"
        KUL_task_exec $verbose_level "MSBP" "7_msbp"

        echo "Done MSBP"
        touch KUL_LOG/sub-${participant}_MSBP.done
        
    else
        echo "MSBP already done"
    fi
}

function KUL_run_FWT {
    config="tracks_list.txt"
    if [ ! -f KUL_LOG/sub-${participant}_ses-${ses}_FWT.done ]; then

        task_in="KUL_FWT_make_VOIs.sh -p ${participant} \
        -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
        -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
        -c $cwd/study_config/${config} \
        -d $cwd/dwiprep/sub-${participant}/ses-${ses} \
        -o $kulderivativesdir/sub-${participant}/ses-${ses}/FWT \
        -n $ncpu"
        KUL_task_exec $verbose_level "KUL_FWT voi generation" "FWTvoi"


        task_in="KUL_FWT_make_TCKs.sh -p ${participant} \
        -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
        -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
        -c $cwd/study_config/${config} \
        -d $cwd/dwiprep/sub-${participant}/ses-${ses} \
        -o $kulderivativesdir/sub-${participant}/ses-${ses}/FWT \
        -T 1 -a iFOD2 \
        -f 1 \
        -Q -S \
        -n $ncpu"
        KUL_task_exec $verbose_level "KUL_FWT tract generation" "FWTtck"

        #ln -s $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*fin_map_BT_iFOD2.nii.gz $globalresultsdir/Tracto/
        #ln -s $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*fin_BT_iFOD2.tck $globalresultsdir/Tracto/
        #mcp "$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*_fin_map_BT_iFOD2.nii.gz" \
        #    "$globalresultsdir/Tracto/Tract-csd_#2.nii.gz"
        #mcp "$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*_fin_BT_iFOD2.tck" \
        #    "$globalresultsdir/Tracto/Tract-csd_#2.tck"
        #pdfunite $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*_output/Screenshots/*fin_BT_iFOD2_inMNI_screenshot2_niGB.pdf $globalresultsdir/Tracto/Tracts_Summary.pdf
        touch KUL_LOG/sub-${participant}_ses-${ses}_FWT.done
        
    else
        echo "FWT already done"
    fi
}

################## MAIN ##############################
# GLOBAL defs
globalresultsdir=$cwd/RESULTS/sub-$participant

KUL_check_participant

kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
mkdir -p $kulderivativesdir
mkdir -p $globalresultsdir/Anat
mkdir -p $globalresultsdir/SPM
mkdir -p $globalresultsdir/Melodic
mkdir -p $globalresultsdir/Tracto

if [ $KUL_DEBUG -gt 0 ]; then 
    echo "kulderivativesdir: $kulderivativesdir"
    echo "globalresultsdir: $globalresultsdir"
fi

# Run BIDS validation
check_in=${KUL_LOG_DIR}/1_bidscheck.done
if [ ! -f $check_in ]; then

    docker run -ti --rm -v ${cwd}/BIDS:/data:ro bids/validator /data

    read -p "Are you happy? (y/n) " answ
    if [[ ! "$answ" == "y" ]]; then
        exit 1
    else
        touch $check_in
    fi
fi

# Check if fMRI and/or dwi data are present and/or to redo some processing
echo "Starting KUL_clinical_fmridti"

# STEP 1 - define major files
T1w=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_desc-preproc_T1w.nii.gz

# STEP 2 - run fmriprep and continue
KUL_run_fmriprep &


# STEP 3 - run dwiprep and continue
if [ ! -f dwiprep/sub-${participant}/dwiprep_is_done.log ]; then
    cp study_config/run_dwiprep.txt KUL_LOG/sub-${participant}_run_dwiprep.txt
    sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_dwiprep.txt
    rm -f KUL_LOG/sub-${participant}_run_dwiprep.txt.bck
    KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_dwiprep.txt 
else
    echo "Dwiprep already done"
fi

wait


# STEP 4 - run msbp (including freesurfer)
mkdir -p ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single/sub-${participant}/anat
cp ${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_desc-preproc_T1w.json \
    ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single/sub-${participant}/anat/sub-${participant}_T1w.json
cp ${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_desc-preproc_T1w.nii.gz \
    ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single/sub-${participant}/anat/sub-${participant}_T1w.nii.gz
cp -rf $cwd/BIDS/.bidsignore ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single
cp -rf $cwd/BIDS/CHANGES ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single
cp -rf $cwd/BIDS/participants.json ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single
cp -rf $cwd/BIDS/participants.tsv ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single
cp -rf $cwd/BIDS/README ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single
cp -rf $cwd/BIDS/dataset_description.json ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single
mkdir -p ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single/code
mkdir -p ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single/derivatives
mkdir -p ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single/sourcedata

if [ ! -f KUL_LOG/sub-${participant}_MSBP.done ]; then

        echo "Running MSBP"

        # there seems tpo be a problem with docker if the fsaverage dir is a soft link; so we delete the link and hardcopy it
        rm -fr $cwd/BIDS/derivatives/freesurfer/fsaverage
        cp -r $FREESURFER_HOME/subjects/fsaverage $cwd/BIDS/derivatives/freesurfer/fsaverage


        docker run --rm -u $(id -u) -v ${cwd}/BIDS/derivatives/KUL_compute/BIDS_single:/bids_dir \
         -v $cwd/BIDS/derivatives:/output_dir \
         -v $FS_LICENSE:/opt/freesurfer/license.txt \
         sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
         --participant_label $participant --isotropic_resolution 1.0 --thalamic_nuclei \
         --brainstem_structures --skip_bids_validator --fs_number_of_cores $ncpu \
         --multiproc_number_of_cores $ncpu
        #KUL_task_exec $verbose_level "MSBP" "7_msbp"

        echo "Done MSBP"
        touch KUL_LOG/sub-${participant}_MSBP.done
        
    else
        echo "MSBP already done"
    fi


# STEP 5 - run KUL_dwiprep_anat
task_in="KUL_dwiprep_anat.sh -p $participant -n $ncpu -m"
KUL_task_exec $verbose_level "KUL_dwiprep_anat" "6_dwiprep_anat"

# STEP 6 - run FWT
sessions=("T0 T1 T2")
for s in $sessions; do
    echo "Running FWT on session $s"
    ses=$s
    KUL_run_FWT
done
