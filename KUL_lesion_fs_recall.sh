#!/bin/bash

# set -x

# @ Ahmed Radwan ahmed.radwan@kuleuven.be # # v 1.0 - dd 20/03/2019 - dev (Alpha)
# @ Stefan Sunaert stefan.sunaert@kuleuven.be



v="1.0 - dd 20/03/2019"

# This script is meant for allowing a decent recon-all output in the presence of a large brain lesion 
# It is not a final end-all solution but a rather crude and simplistic work around 
# The main idea is to replace the lesion with a hole and fill the hole with information from normal hemisphere # maintains subject specificity and diseased hemisphere information but replaces lesion tissue with sham brain 
# should be followed by a loop calculating overlap between fake labels resulting from sham brain with actual lesion 
#  To do: # 
# - A similar approach for bilateral lesions could still work, using a population template perhaps ?
# - make sure there is only one of each modality in the BIDS folder, since we're not using a BIDS validator
# - add longitudinal option
# - add MNI images as input to recon-all also
# - test -a -t -f options without -b
# - test -s option in presence of multiple sessions
# - add N4 bias correction
# - improve lesion fill patch match to native brain (elastic reg ?)



# Description: 
# Info and Instructions: 
# This functiona in its current state is a WIP for S61759 AR,SK,EC,TT,PD,SS 
# - KUL_Lesion_FSrecall takes as input T1, T2, FLAIR and a lesion mask 
# - The end result is running Freesurfer recon-all on lesioned brains 
# - We start with unprocessed data 
# - At least 2 modalities (T1 + * ) are required plus a lesion mask 
# - The idea is to use warped healthy subject specific brain tissue to fill the lesioned area. 
# - Takes as input: subject_label, BIDS_dir, lesion mask file path & name 
# Steps: 
# - The input images are rigidly aligned to the T1, T1 affine 2 MNI, then all to MNI 
# - The MNI images are flipped in LR # % - The lesion mask in binarized and inverted in MNI + either smoothed or dilated 
# - This is followed by a rough brain extraction in MNI 
# - The lesion & brain masks are combined with apriori left right split masks in MNI 
# - Left & right hemispheres are split from the rough BET in MNI images 
# - Determine if lesion is R or L, and warp as follows 
# - Use lesioned_side-lesion_mask as target mask (fixed image MNI unflip brain) 
# - Use lesioned_side_w_lesion mask as input mask (source image MNI flip brain) 
# - Lesion fill is created out of the warped healthy to unhealthy brain tissue 
# - Lesion fill is inserted in place of lesion in original MNI images 
# - Inverse warp to native space T1 
# - Use both native T1&T2 and MNIT1&T2 to run recon-all as if we have two sessions with two modalities each. 
# - Calculate percentage volume overlap with each lobe. 
# - Move results to KUL_NITs location # % @ Ahmed Radwan (ahmed.radwan@kuleuven.be) 07/03/2019 
# Requires: FSL, Freesurfer, T1 and T2 WIs and lesion mask out of itk-snap 
# to add a lesion argument later 
# - lesion mask name should include space, modality, and adhere to BIDS naming scheme 


# ----------------------------------- MAIN --------------------------------------------- 
# this script uses "preprocessing control", i.e. if some steps are already processed it will skip these

kul_lesion_dir=`dirname "$0"`
script=`basename "$0"`
# source $kul_main_dir/KUL_main_functions.sh
cwd=($(pwd))

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` preps structural images with lesions and runs recon-all.

Usage:

  `basename $0` -p subject <OPT_ARGS> -l <OPT_ARGS> -z <OPT_ARGS> -b 
  
  or
  
  `basename $0` -p subject <OPT_ARGS> -a <OPT_ARGS> -b <OPT_ARGS> -c <OPT_ARGS> -l <OPT_ARGS> -z <OPT_ARGS>  
  
Examples:

  `basename $0` -p pat001 -b -n 6 -l /fullpath/lesion_T1w.nii.gz -z T1 -o /fullpath/output
  `basename $0` -p pat001 -n 6 -a /fullpath/T1w.nii.gz -b /fullpath/T2w.nii.gz -l /fullpath/lesion_T2w.nii.gz -z T2
  `basename $0` -p pat001 -n 6 -a /fullpath/T1w.nii.gz -c /fullpath/flair.nii.gz -l /fullpath/lesion_flair.nii.gz -z FLAIR
	

Required arguments:

	-p:  BIDS participant name (anonymised name of the subject without the "sub-" prefix)
	-b:  if data is in BIDS
	-l:  full path and file name to lesion mask file per session
	-z:  space of the lesion mask used (T1, T2, or FLAIR)
	-a:  full path and file name to T1
	-t:  full path and file name to T2
	-f:  full path and file name to T2 FLAIR
	-m:  full path to intermediate output dir
	-o:  full path to output dir (if not set reverts to default output ./lesion_wf_output)
	(This workflow require at least 2 modalities T1 + T2 or FLAIR)

Optional arguments:

	-s:  session (of the participant)
	-n:  number of cpu for parallelisation
	-v:  show output from mrtrix commands
	-h:  prints help menu


USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
ncpu=6
silent=1

# Set required options
p_flag=0
bids_flag=0
s_flag=0
l_flag=0
l_spaceflag=0
t1_flag=0
t2_flag=0
flair_flag=0
o_flag=0
m_flag=0
n_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "p:a:t:f:l:z:s:o:m:n:bvh" OPT; do

        case $OPT in
        p) #subject
            p_flag=1
            subj=$OPTARG
        ;;
        b) #BIDS or not ?
            bids_flag=1
        ;;
        a) #T1 WIs
            t1_flag=1
			t1_orig=$OPTARG
        ;;
        t) #T2 WIs
            t2_flag=1
			t2_orig=$OPTARG
        ;;
        f) #Flair WIs
            flair_flag=1
			flair_orig=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        l) #lesion_mask
            l_flag=1
            L_mask=$OPTARG
		;;
	    z) #lesion_mask
	        l_spaceflag=1
	        L_mask_space=$OPTARG		
	    ;;
	    m) #intermediate output
			m_flag=1
			wf_dir=$OPTARG		
        ;;
	    o) #output
			o_flag=1
			out_dir=$OPTARG		
        ;;
        n) #parallel
			n_flag=1
            ncpu=$OPTARG
        ;;
        v) #verbose
            silent=0
        ;;
        h) #help
            Usage >&2
            exit 0
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

# check for required inputs and define your workflow accordingly

if [[ $p_flag -eq 0 || $l_flag -eq 0 || $l_spaceflag -eq 0 ]]; then
	
    echo
    echo "Inputs -p -lesion -lesion_space must be set." >&2
    echo
    exit 2
	
else
	
	echo "Inputs are -p " $subj " -lesion " $L_mask " -lesion_space " $L_mask_space
	
fi
	
if [[ $bids_flag -eq 1 && $s_flag -eq 0 ]]; then
		
		# bids flag defined but not session flag
	    search_sessions=($(find ${cwd}/BIDS/sub-${subj} -type d | grep anat));
		num_sessions=${#search_sessions[@]};
		ses_long="";
		
		
		if [[ $num_sessions -eq 1 ]]; then 
			
			echo " we have one session in the BIDS dir, this is good."
			
			# now we need to search for the images
			# here we also need to search for the images
			# then also find which modalities are available and set wf accordingly
			
			search_T1=($(find $search_sessions -type f | grep T1w.nii.gz));
			search_T2=($(find $search_sessions -type f | grep T2w.nii.gz));
			search_FLAIR=($(find $search_sessions -type f | grep FLAIR.nii.gz));
			
			if [[ $search_T1 ]]; then
				
				T1_orig=$search_T1
				echo " We found T1 WIs " $T1_orig
				
			else
				
				echo " no T1 WIs found in BIDS dir, exiting"
				exit 2
				
			fi
			
			if [[ $search_T2 && ! $search_FLAIR ]]; then 
				
				wf=1;
				T2_orig=$search_T2;
				echo " We found also T2 WIs, but no FLAIR WIs "
				
			elif [[ ! $search_T2 && $search_FLAIR ]]; then 
				
				wf=2;
				FLAIR_orig=$search_FLAIR;
				echo " We found also FLAIR WIs, but no T2 WIs "
				
			elif [[ $search_T2 && $search_FLAIR ]]; then 
				
				wf=3;
				T2_orig=$search_T2;
				FLAIR_orig=$search_FLAIR;
				echo " We found also T2 WIs, and FLAIR WIs "
				
			else 
				
				echo " This script requires at least T1 WIs + either T2 or FLAIR WIs, exiting."
				exit 2
				
			fi
			
		else 
			
			echo " There's a problem with sessions in BIDS dir. "
			echo " Please double check your data structure &/or specify one session with -s if you have multiple ones. "
			exit 2
			
		fi
		
elif [[ $bids_flag -eq 1 && $s_flag -eq 1 ]]; then
		
		# this is fine
		search_sessions=($(find ${cwd}/BIDS/sub-${subj}_ses-${ses} -type d | grep anat));
		num_sessions=1;
		ses_long=ses-0${num_sessions}_;
		# here we also need to search for the images
		# then also find which modalities are available and set wf accordingly
		
		if [[ $num_sessions -eq 1 ]]; then 
			
			echo " One session " $ses " specified in BIDS dir, good."
			# now we need to search for the images
			# here we also need to search for the images
			# then also find which modalities are available and set wf accordingly
			
			search_T1=($(find $search_sessions -type f | grep T1w.nii.gz));
			search_T2=($(find $search_sessions -type f | grep T2w.nii.gz));
			search_FLAIR=($(find $search_sessions -type f | grep flair.nii.gz));
			
			if [[ $search_T1 ]]; then
				
				T1_orig=$search_T1;
				echo " We found T1 WIs " $T1_orig
				
			else
				
				echo " no T1 WIs found in BIDS dir, exiting "
				exit 2
				
			fi
			
			if [[ $search_T2 && ! $search_FLAIR ]]; then 
				
				wf=1;
				T2_orig=$search_T2;
				echo " We found also T2 WIs, but no FLAIR WIs "
				
			elif [[ ! $search_T2 && $search_FLAIR ]]; then 
				
				wf=2;
				FLAIR_orig=$search_FLAIR;
				echo " We found also FLAIR WIs, but no T2 WIs "
				
			elif [[ $search_T2 && $search_FLAIR ]]; then 
				
				wf=3;
				T2_orig=$search_T2;
				FLAIR_orig=$search_FLAIR;
				echo " We found also T2 WIs, and FLAIR WIs "
				
			else 
				
				echo " This script requires at least T1 WIs + either T2 or FLAIR WIs, exiting."
				exit 2
				
			fi
		
elif [[ $bids_flag -eq 0 && $s_flag -eq 0 ]]; then
		
		# this is fine if T1 and T2 and/or flair are set
		# find which ones are set and define wf accordingly
		ses_long="";
		
		if [[ $t1_flag ]]; then
			
			T1_orig=$t1_orig
			
		else
			
			echo " No T1 WIs specified, exiting. "
			exit 2
			
		fi
		
		if [[ $seach_T1 && $search_T2 && ! $search_FLAIR ]]; then 
			
			wf=1;
			T2_orig=$search_T2;
			echo " We also have T2 WIs, but no FLAIR WIs "
			
		elif [[ $seach_T1 && ! $search_T2 && $search_FLAIR ]]; then 
			
			wf=2;
			FLAIR_orig=$search_FLAIR;
			echo " We also have FLAIR WIs, but no T2 WIs "
			
		elif [[ $seach_T1 && $search_T2 && $search_FLAIR ]]; then 
			
			wf=3;
			T2_orig=$search_T2;
			FLAIR_orig=$search_FLAIR;
			echo " We also have T2 WIs, and FLAIR WIs "
			
		else 
			
			echo " This script requires at least T1 WIs + either T2 or FLAIR WIs, exiting."
			exit 2
			
		fi
		
elif [[ $bids_flag -eq 0 && $s_flag -eq 1 ]]; then
			
		echo " Wrong optional arguments, we can't have sessions without BIDS, exiting."
		exit 2
		
fi

fi

function_path=($(which KUL_lesion_fs_recall.sh | rev | cut -d"/" -f2- | rev))


# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Some parallelisation

if [[ $n_flag -eq 0 ]]; then

ncpu=6

echo " -n flag not set, using default 6 threads. "

else

echo " -n flag set, using " ${ncpu} " threads."

fi

# FSLPARALLEL=$ncpu; export FSLPARALLEL
#
# OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

d=$(date "+%Y-%m-%d_%H-%M-%S");
log=log/log_${d}.txt;


# --- MAIN ----------------
# here we give session invariate variables (e.g. template images)
# only in for WIP - R 09-03-2019
# subj=BAC

# The necessary priors
# need to use multiple templates here...
# will need to include a folder with priors

MNI_T1=${function_path}/Templates/MNI_T1.nii.gz

MNI_T1_brain=${function_path}/Templates/MNI_T1_brain.nii.gz

MNI_T2=${function_path}/Templates/MNI_T2.nii.gz

MNI_T2_brain=${function_path}/Templates/MNI_T2_brain.nii.gz

MNI_FLAIR=${function_path}/Templates/MNI_FLAIR.nii.gz

MNI_FLAIR_brain=${function_path}/Templates/MNI_FLAIR_brain.nii.gz

MNI_brain_mask=${function_path}/Templates/MNI_brain_mask.nii.gz

MNI_rl=${function_path}/Templates/MNI_RL1c_labels.nii.gz

# Either a session is given on the command line
# If not the session(s) need to be determined.


# ---- BIG LOOP for processing each session
#
for current_session in `seq 0 $(($num_sessions-1))`; do

	cd $cwd

    long_bids_subj=${search_sessions[$current_session]}

    echo $long_bids_subj

	bids_subj=${long_bids_subj%anat}

	echo $bids_subj

	lesion_wf=${cwd}/lesion_wf

	# output

	if [[ $o_flag -eq 1 ]]; then
		output=$out_dir

	else

		output=${lesion_wf}/lesion_wf_output/sub-PT01${ses_long}

	fi

	# intermediate folder

	if [[ $m_flag -eq 1 ]]; then

		preproc=$wf_dir

	else

		preproc=${lesion_wf}/lesion_wf_preproc/sub-PT01${ses_long}

	fi

	echo $lesion_wf

ROIs=${output}/ROIs
overlap=${output}/overlap

# make your dirs
mkdir -p ${preproc}
mkdir -p ${output}
mkdir -p ${ROIs}
mkdir -p ${overlap}


# The lesion as provided (this depends on the session)

L_mask_bin=${L_mask::${#L_mask}-7}.nii.gz

fslmaths $L_mask -bin $L_mask_bin

L_mask_orig=${L_mask_bin}

L_mask_binv_in_T1=${preproc}/sub-${subj}_${ses_long}L_mask_binv_T1.nii.gz

L_mask_in_MNI=${preproc}/sub-${subj}_${ses_long}L_mask_in_MNI.nii.gz

L_mask_in_MNI_dil=${preproc}/sub-${subj}_${ses_long}L_mask_dil_in_MNI.nii.gz

L_mask_in_MNI_dil_binv=${preproc}/sub-${subj}_${ses_long}L_mask_dil_binv_in_MNI.nii.gz

MNI_mask_min_lesion=${preproc}/sub-${subj}_${ses_long}MNI_mask_min_lesion.nii.gz

lesion_left_hemi_overlap=${preproc}/sub-${subj}_${ses_long}lesion_left_hemi_overlap.nii.gz

lesion_right_hemi_overlap=${preproc}/sub-${subj}_${ses_long}lesion_right_hemi_overlap.nii.gz

# vars for images
T2_in_T1=${preproc}/sub-${subj}_${ses_long}T2_in_T1_Warped.nii.gz

FLAIR_in_T1=${preproc}/sub-${subj}_${ses_long}FLAIR_in_T1_Warped.nii.gz

T2_to_T1_affine=${preproc}/sub-${subj}_${ses_long}T2_in_T1_0GenericAffine.mat

FLAIR_to_T1_affine=${preproc}/sub-${subj}_${ses_long}FLAIR_in_T1_0GenericAffine.mat

T1_in_MNI=${preproc}/sub-${subj}_${ses_long}T1_in_MNI_Warped.nii.gz

MNI_brain_mask_in_nat=${preproc}/sub-${subj}_${ses_long}MNI_brain_mask_in_T1nat.nii.gz

MNI_brain_mask_in_nat_min_lesion=${preproc}/sub-${subj}_${ses_long}MNI_brain_mask_in_T1nat_minlesion.nii.gz

T1_in_MNI_brain=${preproc}/sub-${subj}_${ses_long}T1_in_MNI_BrainExtractionBrain.nii.gz

subj_brain_mask_in_MNI=${preproc}/sub-${subj}_${ses_long}T1_in_MNI_BrainExtractionMask.nii.gz

T1_in_MNI_brain_flip=${preproc}/sub-${subj}_${ses_long}flip_T1_in_MNI_brain.nii.gz

T1_in_MNI_SyN_brain_flip=${preproc}/sub-${subj}_${ses_long}flip_T1_brain_in_MNI_SyN_Warped.nii.gz

T2_in_MNI=${preproc}/sub-${subj}_${ses_long}T2_in_MNI_Warped.nii.gz

T2_in_MNI_brain=${preproc}/sub-${subj}_${ses_long}T2_in_MNI_brain.nii.gz

T2_in_MNI_brain_flip=${preproc}/sub-${subj}_${ses_long}flip_T2_in_MNI_brain.nii.gz

T2_in_MNI_SyN_brain_flip=${preproc}/sub-${subj}_${ses_long}flip_T2_brain_in_MNI_SyN_Warped.nii.gz

FLAIR_in_MNI=${preproc}/sub-${subj}_${ses_long}FLAIR_in_MNI_Warped.nii.gz

FLAIR_in_MNI_brain=${preproc}/sub-${subj}_${ses_long}FLAIR_in_MNI_brain.nii.gz

FLAIR_in_MNI_brain_flip=${preproc}/sub-${subj}_${ses_long}flip_FLAIR_in_MNI_brain.nii.gz

FLAIR_in_MNI_SyN_brain_flip=${preproc}/sub-${subj}_${ses_long}flip_FLAIR_brain_in_MNI_SyN_Warped.nii.gz

T1_to_MNI_affine=${preproc}/sub-${subj}_${ses_long}T1_in_MNI_0GenericAffine.mat

T1_with_lesion_fill_MNI=${output}/sub-${subj}_${ses_long}MNI_T1_sim_lesionless.nii.gz

T2_with_lesion_fill_MNI=${output}/sub-${subj}_${ses_long}MNI_T2_sim_lesionless.nii.gz

FLAIR_with_lesion_fill_MNI=${output}/sub-${subj}_${ses_long}MNI_FLAIR_sim_lesionless.nii.gz

T1_with_lesion_fill_T1nat=${output}/sub-${subj}_${ses_long}T1nat_T1_sim_lesionless_Warped.nii.gz

T2_with_lesion_fill_T1nat=${output}/sub-${subj}_${ses_long}T1nat_T2_sim_lesionless_Warped.nii.gz

FLAIR_with_lesion_fill_T1nat=${output}/sub-${subj}_${ses_long}T1nat_FLAIR_sim_lesionless_Warped.nii.gz

MNI_left=${preproc}/MNI_left_bin.nii.gz

MNI_right=${preproc}/MNI_right_bin.nii.gz


# kul_e2cl " Start processing $bids_subj" ${preproc}/${log}

cd $preproc
		
# determine which workflow we need to apply and run it.
# this is distributed in nested if loops
# processing control flags are called wf_mark, listed below.


wf_mark1=${preproc}"/first_part_done.done"

wf_mark2=${preproc}"/second_part_done.done"

wf_mark3=${preproc}"/third_part_done.done"

search_wf_mark1=($(find ${preproc} -type f | grep first_part_done.done));

search_wf_mark2=($(find ${preproc} -type f | grep second_part_done.done));

search_wf_mark3=($(find ${preproc} -type f | grep third_part_done.done));

if [[ ! $search_wf_mark1 ]] ; then

antsRegistrationSyNQuick.sh -d 3 -f $MNI_T1 -m $T1_orig -t a -o ${preproc}/sub-${subj}_${ses_long}T1_in_MNI_ -n ${ncpu}

# This screws me over big time! let's do a respectable BET

# fslmaths $T1_in_MNI -mas $subj_brain_mask_in_MNI $T1_in_MNI_brain

antsBrainExtraction.sh -d 3 -a ${T1_in_MNI} -e ${MNI_T1} -m ${MNI_brain_mask} -o ${preproc}/sub-${subj}_${ses_long}T1_in_MNI_ -u 1

sleep 2

WarpImageMultiTransform 3 $subj_brain_mask_in_MNI $MNI_brain_mask_in_nat -R $T1_orig -i $T1_to_MNI_affine

fslmaths $MNI_brain_mask_in_nat -bin $MNI_brain_mask_in_nat

fslswapdim $T1_in_MNI_brain -x y z $T1_in_MNI_brain_flip

fslorient -swaporient $T1_in_MNI_brain_flip

echo " first part done. " >> $wf_mark1

else
	
	echo " fist part already done, skipping. "
	
fi


if  [[ ! $search_wf_mark2 ]] ; then
		
	if [[ $wf -eq 1 ]]; then

		antsRegistrationSyNQuick.sh -d 3 -f $T1_orig -m $T2_orig -t a -o $preproc/sub-${subj}_${ses_long}T2_in_T1_ -n ${ncpu}

		WarpImageMultiTransform 3 $T2_in_T1 $T2_in_MNI -R $MNI_T1 $T1_to_MNI_affine

		fslmaths $T2_in_MNI -mas $subj_brain_mask_in_MNI $T2_in_MNI_brain

		if [[ $L_mask_space == "T1" ]] ; then

			L_mask_in_T1=$L_mask_orig

		elif [[ $L_mask_space == "T2" ]] ; then

			L_mask_in_T2=$L_mask_orig

			L_mask_in_T1=$preproc/sub-${subj}_${ses_long}L_mask_in_T2_to_T1.nii.gz

			WarpImageMultiTransform 3 $L_mask_in_T2 $L_mask_in_T1 -R $T1_orig $T2_to_T1_affine

			fslmaths $L_mask_in_T1 -bin $L_mask_in_T1

		fi

		# flip T2 in MNI

		fslswapdim $T2_in_MNI_brain -x y z $T2_in_MNI_brain_flip

		fslorient -swaporient $T2_in_MNI_brain_flip

		WarpImageMultiTransform 3 $L_mask_in_T1 $L_mask_in_MNI -R $MNI_T1 $T1_to_MNI_affine

		fslmaths $L_mask_in_MNI -dilM -bin $L_mask_in_MNI_dil

		fslmaths $L_mask_in_MNI_dil -binv $L_mask_in_MNI_dil_binv

		fslmaths $subj_brain_mask_in_MNI -mas $L_mask_in_MNI_dil_binv $MNI_mask_min_lesion

		# SyN warp the flipped brains to the respective MNI template

		antsRegistrationSyNQuick.sh -d 3 -f $MNI_T1_brain -m $T1_in_MNI_brain_flip -t s \
			-x $MNI_brain_mask,$MNI_mask_min_lesion -o $preproc/sub-${subj}_${ses_long}flip_T1_brain_in_MNI_SyN_ -n ${ncpu}

		antsRegistrationSyNQuick.sh -d 3 -f $MNI_T2_brain -m $T2_in_MNI_brain_flip -t s \
			-x $MNI_brain_mask,$MNI_mask_min_lesion -o $preproc/sub-${subj}_${ses_long}flip_T2_brain_in_MNI_SyN_ -n ${ncpu}

	elif [[ $wf -eq 2 ]] ; then

		antsRegistrationSyNQuick.sh -d 3 -f $T1_orig -m $FLAIR_orig -t a -o $preproc/sub-${subj}_${ses_long}FLAIR_in_T1_ -n ${ncpu}

		WarpImageMultiTransform 3 $FLAIR_in_T1 $FLAIR_in_MNI -R $MNI_T1 $T1_to_MNI_affine

		fslmaths $FLAIR_in_MNI -mas $subj_brain_mask_in_MNI $FLAIR_in_MNI_brain

		if [[ $L_mask_space == "T1" ]] ; then

			L_mask_in_T1=$L_mask_orig

		elif [[ $L_mask_space == "FLAIR" ]] ; then

			L_mask_in_FLAIR=$L_mask_orig

			L_mask_in_T1=$preproc/sub-${subj}_${ses_long}L_mask_in_FLAIR_to_T1.nii.gz

			WarpImageMultiTransform 3 $L_mask_in_FLAIR $L_mask_in_T1 -R $T1_orig $FLAIR_to_T1_affine

			fslmaths $L_mask_in_T1 -bin $L_mask_in_T1

		fi

		# flip FLAIR in MNI

		fslswapdim $FLAIR_in_MNI_brain -x y z $FLAIR_in_MNI_brain_flip

		fslorient -swaporient $FLAIR_in_MNI_brain_flip

		WarpImageMultiTransform 3 $L_mask_in_T1 $L_mask_in_MNI -R $MNI_T1 $T1_to_MNI_affine

		fslmaths $L_mask_in_MNI -dilM -bin $L_mask_in_MNI_dil

		fslmaths $L_mask_in_MNI_dil -binv $L_mask_in_MNI_dil_binv

		fslmaths $subj_brain_mask_in_MNI -mas $L_mask_in_MNI_dil_binv $MNI_mask_min_lesion

		# SyN warp the flipped brains to the respective MNI template

		antsRegistrationSyNQuick.sh -d 3 -f $MNI_T1_brain -m $T1_in_MNI_brain_flip -t s \
			-x $MNI_brain_mask,$MNI_mask_min_lesion -o $preproc/sub-${subj}_${ses_long}flip_T1_brain_in_MNI_SyN_ -n ${ncpu}

		antsRegistrationSyNQuick.sh -d 3 -f $MNI_FLAIR_brain -m $FLAIR_in_MNI_brain_flip -t s \
			-x $MNI_brain_mask,$MNI_mask_min_lesion -o $preproc/sub-${subj}_${ses_long}flip_FLAIR_brain_in_MNI_SyN_ -n ${ncpu}


	elif [[ $wf -eq 3 ]] ; then

		antsRegistrationSyNQuick.sh -d 3 -f $T1_orig -m $FLAIR_orig -t a -o $preproc/sub-${subj}_${ses_long}FLAIR_in_T1_ -n ${ncpu}

		antsRegistrationSyNQuick.sh -d 3 -f $T1_orig -m $T2_orig -t a -o $preproc/sub-${subj}_${ses_long}T2_in_T1_ -n ${ncpu}

		WarpImageMultiTransform 3 $T2_in_T1 $T2_in_MNI -R $MNI_T1 $T1_to_MNI_affine

		WarpImageMultiTransform 3 $FLAIR_in_T1 $FLAIR_in_MNI -R $MNI_T1 $T1_to_MNI_affine

		fslmaths $T1_in_MNI -mas $subj_brain_mask_in_MNI $T1_in_MNI_brain

		fslmaths $T2_in_MNI -mas $subj_brain_mask_in_MNI $T2_in_MNI_brain

		fslmaths $FLAIR_in_MNI -mas $subj_brain_mask_in_MNI $FLAIR_in_MNI_brain

		if [[ $L_mask_space == "T1" ]] ; then

			L_mask_in_T1=$L_mask_orig

		elif [[ $L_mask_space == "T2" ]] ; then

			L_mask_in_T2=$L_mask_orig

			L_mask_in_T1=$preproc/sub-${subj}_${ses_long}L_mask_in_T2_to_T1.nii.gz

			WarpImageMultiTransform 3 $L_mask_in_T2 $L_mask_in_T1 -R $T1_orig $T2_to_T1_affine

			fslmaths $L_mask_in_T1 -bin $L_mask_in_T1

		elif [[ $L_mask_space == "FLAIR" ]] ; then

			L_mask_in_FLAIR=$L_mask_orig

			L_mask_in_T1=$preproc/sub-${subj}_${ses_long}L_mask_in_FLAIR_to_T1.nii.gz

			WarpImageMultiTransform 3 $L_mask_in_FLAIR $L_mask_in_T1 -R $T1_orig $FLAIR_to_T1_affine

			fslmaths $L_mask_in_T1 -bin $L_mask_in_T1

		fi

		# flip T2 and FLAIR in MNI

		fslswapdim $T2_in_MNI_brain -x y z $T2_in_MNI_brain_flip

		fslorient -swaporient $T2_in_MNI_brain_flip

		fslswapdim $FLAIR_in_MNI_brain -x y z $FLAIR_in_MNI_brain_flip

		fslorient -swaporient $FLAIR_in_MNI_brain_flip

		WarpImageMultiTransform 3 $L_mask_in_T1 $L_mask_in_MNI -R $MNI_T1 $T1_to_MNI_affine

		fslmaths $L_mask_in_MNI -dilM -bin $L_mask_in_MNI_dil

		fslmaths $L_mask_in_MNI_dil -binv $L_mask_in_MNI_dil_binv

		fslmaths $subj_brain_mask_in_MNI -mas $L_mask_in_MNI_dil_binv $MNI_mask_min_lesion

		# SyN warp the flipped brains to the respective MNI template

		antsRegistrationSyNQuick.sh -d 3 -f $MNI_T1_brain -m $T1_in_MNI_brain_flip -t s \
			-x $MNI_brain_mask,$MNI_mask_min_lesion -o $preproc/sub-${subj}_${ses_long}flip_T1_brain_in_MNI_SyN_ -n ${ncpu}

		antsRegistrationSyNQuick.sh -d 3 -f $MNI_T2_brain -m $T2_in_MNI_brain_flip -t s \
			-x $MNI_brain_mask,$MNI_mask_min_lesion -o $preproc/sub-${subj}_${ses_long}flip_T2_brain_in_MNI_SyN_ -n ${ncpu}

		antsRegistrationSyNQuick.sh -d 3 -f $MNI_FLAIR_brain -m $FLAIR_in_MNI_brain_flip -t s \
			-x $MNI_brain_mask,$MNI_mask_min_lesion -o $preproc/sub-${subj}_${ses_long}flip_FLAIR_brain_in_MNI_SyN_ -n ${ncpu}

	fi
	
	echo " second part done " >> $wf_mark2
	
	
else 
	
	echo " second part already done, skipping. "
	
	L_mask_in_T1=($(find ${preproc} -type f | grep _to_T1.nii.gz));
	
	
fi


	# all workflows collapse to this same step at this point


	T1_fake_left_hemi=${preproc}/T1_fake_left_hemi.nii.gz
	T2_fake_left_hemi=${preproc}/T2_fake_left_hemi.nii.gz
	FLAIR_fake_left_hemi=${preproc}/FLAIR_fake_left_hemi.nii.gz
	
	T1_fake_right_hemi=${preproc}/T1_fake_right_hemi.nii.gz
	T2_fake_right_hemi=${preproc}/T2_fake_right_hemi.nii.gz
	FLAIR_fake_right_hemi=${preproc}/FLAIR_fake_right_hemi.nii.gz
	
	T1_fake_left_hemi_mask=${preproc}/T1_fake_left_hemi_mask.nii.gz
	T2_fake_left_hemi_mask=${preproc}/T2_fake_left_hemi_mask.nii.gz
	FLAIR_fake_left_hemi_mask=${preproc}/FLAIR_fake_left_hemi_mask.nii.gz

	T1_fake_right_hemi_mask=${preproc}/T1_fake_right_hemi_mask.nii.gz
	T2_fake_right_hemi_mask=${preproc}/T2_fake_right_hemi_mask.nii.gz
	FLAIR_fake_right_hemi_mask=${preproc}/FLAIR_fake_right_hemi_mask.nii.gz
	
	T1_real_left_hemi=${preproc}/T1_real_left_hemi.nii.gz
	T2_real_left_hemi=${preproc}/T2_real_left_hemi.nii.gz
	FLAIR_real_left_hemi=${preproc}/FLAIR_real_left_hemi.nii.gz
	
	T1_real_right_hemi=${preproc}/T1_real_right_hemi.nii.gz
	T2_real_right_hemi=${preproc}/T2_real_right_hemi.nii.gz
	FLAIR_real_right_hemi=${preproc}/FLAIR_real_right_hemi.nii.gz
	
	T1_real_left_hemi_mask_min_lesion=${preproc}/T1_real_left_hemi_mask_min_lesion.nii.gz
	T2_real_left_hemi_mask_min_lesion=${preproc}/T2_real_left_hemi_mask_min_lesion.nii.gz
	FLAIR_real_left_hemi_mask_min_lesion=${preproc}/FLAIR_real_left_hemi_mask_min_lesion.nii.gz
	
	T1_real_right_hemi_mask_min_lesion=${preproc}/T1_real_right_hemi_mask_min_lesion.nii.gz
	T2_real_right_hemi_mask_min_lesion=${preproc}/T2_real_right_hemi_mask_min_lesion.nii.gz
	FLAIR_real_right_hemi_mask_min_lesion=${preproc}/FLAIR_real_right_hemi_mask_min_lesion.nii.gz
	
	T1_left_fake_2_real_hemi_SyN=${preproc}/sub-${subj}_${ses_long}T1_left_fake_2_real_hemiWarped.nii.gz
	T2_left_fake_2_real_hemi_SyN=${preproc}/sub-${subj}_${ses_long}T2_left_fake_2_real_hemiWarped.nii.gz
	FLAIR_left_fake_2_real_hemi_SyN=${preproc}/sub-${subj}_${ses_long}FLAIR_left_fake_2_real_hemiWarped.nii.gz

	T1_right_fake_2_real_hemi_SyN=${preproc}/sub-${subj}_${ses_long}T1_right_fake_2_real_hemiWarped.nii.gz
	T2_right_fake_2_real_hemi_SyN=${preproc}/sub-${subj}_${ses_long}T2_right_fake_2_real_hemiWarped.nii.gz
	FLAIR_right_fake_2_real_hemi_SyN=${preproc}/sub-${subj}_${ses_long}FLAIR_right_fake_2_real_hemiWarped.nii.gz
	
	T1_lesion_fill=${preproc}/sub-${subj}_${ses_long}lesion_fill_T1.nii.gz
	T2_lesion_fill=${preproc}/sub-${subj}_${ses_long}lesion_fill_T2.nii.gz
	FLAIR_lesion_fill=${preproc}/sub-${subj}_${ses_long}lesion_fill_FLAIR.nii.gz


	# determine lesion side
	# we have a problem here!!!!


if [[ ! $search_wf_mark3 ]]; then 
		
	fslmaths $MNI_rl -thr 100 -uthr 100 -bin $MNI_left
	fslmaths $MNI_left -mas $L_mask_in_MNI_dil $lesion_left_hemi_overlap
	overlap_left="$(fslstats $lesion_left_hemi_overlap -V | head -c 1)"

	fslmaths $MNI_rl -thr 1 -uthr 1 -bin $MNI_right
	fslmaths $MNI_right -mas $L_mask_in_MNI_dil $lesion_right_hemi_overlap
	overlap_right="$(fslstats $lesion_right_hemi_overlap -V | head -c 1)"

	fslmaths $L_mask_in_T1 -binv $L_mask_binv_in_T1
	fslmaths $MNI_brain_mask_in_nat -mas $L_mask_binv_in_T1 $MNI_brain_mask_in_nat_min_lesion

	# generate lesion fill and inject it into images in MNI, the hemisphere we manipulate depends on lesion laterality


	if [[ ! $overlap_left -eq 0 ]]; then
		echo "it's a left sided lesion"

		fslmaths $T1_in_MNI_SyN_brain_flip -mas $MNI_left $T1_fake_left_hemi

		fslmaths $T1_fake_left_hemi -bin $T1_fake_left_hemi_mask

		fslmaths $T1_in_MNI_brain -mas $MNI_left $T1_real_left_hemi

		fslmaths $T1_real_left_hemi -bin -mas $L_mask_in_MNI_dil_binv $T1_real_left_hemi_mask_min_lesion

		antsRegistrationSyN.sh -d 3 -f $T1_real_left_hemi -m $T1_fake_left_hemi -t s \
			-x $T1_real_left_hemi_mask_min_lesion,$T1_fake_left_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}T1_left_fake_2_real_hemi -n ${ncpu}

		fslmaths $T1_left_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $T1_lesion_fill

		fslmaths $T1_in_MNI -mas $L_mask_in_MNI_dil_binv -add $T1_lesion_fill $T1_with_lesion_fill_MNI

		antsRegistrationSyN.sh -d 3 -f $T1_orig -m $T1_with_lesion_fill_MNI -t s \
			-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_T1_sim_lesionless_ -n ${ncpu}

		if [[ $wf -eq 1 ]] ; then

			fslmaths $T2_in_MNI_SyN_brain_flip -mas $MNI_left $T2_fake_left_hemi

			fslmaths $T2_fake_left_hemi -bin $T2_fake_left_hemi_mask
			
			fslmaths $T2_in_MNI_brain -mas $MNI_left $T2_real_left_hemi

			fslmaths $T2_real_left_hemi -bin -mas $L_mask_in_MNI_dil_binv $T2_real_left_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $T2_real_left_hemi -m $T2_fake_left_hemi -t s \
				-x $T2_real_left_hemi_mask_min_lesion,$T2_fake_left_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}T2_left_fake_2_real_hemi -n ${ncpu}

			fslmaths $T2_left_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $T2_lesion_fill

			fslmaths $T2_in_MNI -mas $L_mask_in_MNI_dil_binv -add $T2_lesion_fill $T2_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $T2_in_T1 -m $T2_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_T2_sim_lesionless_ -n ${ncpu}

		elif [[ $wf -eq 2 ]] ; then

			fslmaths $FLAIR_in_MNI_SyN_brain_flip -mas $MNI_left $FLAIR_fake_left_hemi

			fslmaths $FLAIR_fake_left_hemi -bin $FLAIR_fake_left_hemi_mask

			fslmaths $FLAIR_in_MNI_brain -mas $MNI_left $FLAIR_real_left_hemi

			fslmaths $FLAIR_real_left_hemi -bin -mas $L_mask_in_MNI_dil_binv $FLAIR_real_left_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $FLAIR_real_left_hemi -m $FLAIR_fake_left_hemi -t s \
				-x $FLAIR_real_left_hemi_mask_min_lesion,$FLAIR_fake_left_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}FLAIR_left_fake_2_real_hemi -n ${ncpu}

			fslmaths $FLAIR_left_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $FLAIR_lesion_fill

			fslmaths $FLAIR_in_MNI -mas $L_mask_in_MNI_dil_binv -add $FLAIR_lesion_fill $FLAIR_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $FLAIR_in_T1 -m $FLAIR_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_FLAIR_sim_lesionless_ -n ${ncpu}

		elif [[ $wf -eq 3 ]] ; then

			fslmaths $T2_in_MNI_SyN_brain_flip -mas $MNI_left $T2_fake_left_hemi

			fslmaths $T2_fake_left_hemi -bin $T2_fake_left_hemi_mask

			fslmaths $T2_in_MNI_brain -mas $MNI_left $T2_real_left_hemi

			fslmaths $T2_real_left_hemi -bin -mas $L_mask_in_MNI_dil_binv $T2_real_left_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $T2_real_left_hemi -m $T2_fake_left_hemi -t s \
				-x $T2_real_left_hemi_mask_min_lesion,$T2_fake_left_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}T2_left_fake_2_real_hemi -n ${ncpu}

			fslmaths $T2_left_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $T2_lesion_fill

			fslmaths $T2_in_MNI -mas $L_mask_in_MNI_dil_binv -add $T2_lesion_fill $T2_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $T2_in_T1 -m $T2_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_T2_sim_lesionless_ -n ${ncpu}

			fslmaths $FLAIR_in_MNI_SyN_brain_flip -mas $MNI_left $FLAIR_fake_left_hemi

			fslmaths $FLAIR_fake_left_hemi -bin $FLAIR_fake_left_hemi_mask

			fslmaths $FLAIR_in_MNI_brain -mas $MNI_left $FLAIR_real_left_hemi

			fslmaths $FLAIR_real_left_hemi -bin -mas $L_mask_in_MNI_dil_binv $FLAIR_real_left_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $FLAIR_real_left_hemi -m $FLAIR_fake_left_hemi -t s \
				-x $FLAIR_real_left_hemi_mask_min_lesion,$FLAIR_fake_left_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}FLAIR_left_fake_2_real_hemi -n ${ncpu}

			fslmaths $FLAIR_left_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $FLAIR_lesion_fill

			fslmaths $FLAIR_in_MNI -mas $L_mask_in_MNI_dil_binv -add $FLAIR_lesion_fill $FLAIR_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $FLAIR_in_T1 -m $FLAIR_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_FLAIR_sim_lesionless_ -n ${ncpu}

		fi

	elif [[ ! $overlap_right -eq 0 ]] ; then
		
		echo "it's a right sided lesion"

		fslmaths $T1_in_MNI_SyN_brain_flip -mas $MNI_right $T1_fake_right_hemi

		fslmaths $T1_fake_right_hemi -bin $T1_fake_right_hemi_mask

		fslmaths $T1_in_MNI_brain -mas $MNI_right $T1_real_right_hemi

		fslmaths $T1_real_right_hemi -bin -mas $L_mask_in_MNI_dil_binv $T1_real_right_hemi_mask_min_lesion

		antsRegistrationSyN.sh -d 3 -f $T1_real_right_hemi -m $T1_fake_right_hemi -t s \
			-x $T1_real_right_hemi_mask_min_lesion,$T1_fake_right_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}T1_right_fake_2_real_hemi -n ${ncpu}

		fslmaths $T1_right_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $T1_lesion_fill

		fslmaths $T1_in_MNI -mas $L_mask_in_MNI_dil_binv -add $T1_lesion_fill $T1_with_lesion_fill_MNI

		antsRegistrationSyN.sh -d 3 -f $T1_orig -m $T1_with_lesion_fill_MNI -t s \
			-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_T1_sim_lesionless_ -n ${ncpu}

		if [[ $wf -eq 1 ]] ; then

			fslmaths $T2_in_MNI_SyN_brain_flip -mas $MNI_right $T2_fake_right_hemi

			fslmaths $T2_fake_right_hemi -bin $T2_fake_right_hemi_mask

			fslmaths $T2_in_MNI_brain -mas $MNI_right $T2_real_right_hemi

			fslmaths $T2_real_right_hemi -bin -mas $L_mask_in_MNI_dil_binv $T2_real_right_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $T2_real_right_hemi -m $T2_fake_right_hemi -t s \
				-x $T2_real_right_hemi_mask_min_lesion,$T2_fake_right_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}T2_right_fake_2_real_hemi -n ${ncpu}

			fslmaths $T2_right_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $T2_lesion_fill

			fslmaths $T2_in_MNI -mas $L_mask_in_MNI_dil_binv -add $T2_lesion_fill $T2_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $T2_in_T1 -m $T2_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_T2_sim_lesionless_ -n ${ncpu}

		elif [[ $wf -eq 2 ]] ; then

			fslmaths $FLAIR_in_MNI_SyN_brain_flip -mas $MNI_right $FLAIR_fake_right_hemi

			fslmaths $FLAIR_fake_right_hemi -bin $FLAIR_fake_right_hemi_mask

			fslmaths $FLAIR_in_MNI_brain -mas $MNI_right $FLAIR_real_right_hemi

			fslmaths $FLAIR_real_right_hemi -bin -mas $L_mask_in_MNI_dil_binv $FLAIR_real_right_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $FLAIR_real_right_hemi -m $FLAIR_fake_right_hemi -t s \
				-x $FLAIR_real_right_hemi_mask_min_lesion,$FLAIR_fake_right_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}FLAIR_right_fake_2_real_hemi -n ${ncpu}

			fslmaths $FLAIR_right_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $FLAIR_lesion_fill

			fslmaths $FLAIR_in_MNI -mas $L_mask_in_MNI_dil_binv -add $FLAIR_lesion_fill $FLAIR_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $FLAIR_in_T1 -m $FLAIR_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_FLAIR_sim_lesionless_ -n ${ncpu}

		elif [[ $wf -eq 3 ]] ; then

			fslmaths $T2_in_MNI_SyN_brain_flip -mas $MNI_right $T2_fake_right_hemi

			fslmaths $T2_fake_right_hemi -bin $T2_fake_right_hemi_mask

			fslmaths $T2_in_MNI_brain -mas $MNI_right $T2_real_right_hemi

			fslmaths $T2_real_right_hemi -bin -mas $L_mask_in_MNI_dil_binv $T2_real_right_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $T2_real_right_hemi -m $T2_fake_right_hemi -t s \
				-x $T2_real_right_hemi_mask_min_lesion,$T2_fake_right_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}T2_right_fake_2_real_hemi -n ${ncpu}

			fslmaths $T2_right_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $T2_lesion_fill

			fslmaths $T2_in_MNI -mas $L_mask_in_MNI_dil_binv -add $T2_lesion_fill $T2_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $T2_in_T1 -m $T2_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_T2_sim_lesionless_ -n ${ncpu}

			fslmaths $FLAIR_in_MNI_SyN_brain_flip -mas $MNI_right $FLAIR_fake_right_hemi

			fslmaths $FLAIR_fake_right_hemi -bin $FLAIR_fake_right_hemi_mask

			fslmaths $FLAIR_in_MNI_brain -mas $MNI_right $FLAIR_real_right_hemi

			fslmaths $FLAIR_real_right_hemi -bin -mas $L_mask_in_MNI_dil_binv $FLAIR_real_right_hemi_mask_min_lesion

			antsRegistrationSyN.sh -d 3 -f $FLAIR_real_right_hemi -m $FLAIR_fake_right_hemi -t s \
				-x $FLAIR_real_right_hemi_mask_min_lesion,$FLAIR_fake_right_hemi_mask -o ${preproc}/sub-${subj}_${ses_long}FLAIR_right_fake_2_real_hemi -n ${ncpu}

			fslmaths $FLAIR_right_fake_2_real_hemi_SyN -mas $L_mask_in_MNI_dil $FLAIR_lesion_fill

			fslmaths $FLAIR_in_MNI -mas $L_mask_in_MNI_dil_binv -add $FLAIR_lesion_fill $FLAIR_with_lesion_fill_MNI

			antsRegistrationSyN.sh -d 3 -f $FLAIR_in_T1 -m $FLAIR_with_lesion_fill_MNI -t s \
				-x $MNI_brain_mask_in_nat_min_lesion,$MNI_mask_min_lesion -o ${output}/sub-${subj}_${ses_long}T1nat_FLAIR_sim_lesionless_ -n ${ncpu}

		fi
		
	fi
	
	echo " part 3 done " >> $wf_mark3
		
else
		
echo " third part already done, skipping to recon-all "
		

fi

	# recon all

	# need to add -T2pial and -FLAIRpial !

	recall_scripts=${output}/sub-${subj}/scripts;
	search_wf_mark4=($(find ${recall_scripts} -type f | grep recon-all.done));
	
if [[ ! ${search_wf_mark1} ]] ; then
	

	if [[ $wf -eq 1 ]]; then

		recon-all -i $T1_with_lesion_fill_T1nat -s sub-${subj} -sd $output -T2 $T2_with_lesion_fill_T1nat -T2pial -openmp ${ncpu} -parallel -all

	elif [[ $wf -eq 2 ]] ; then

		recon-all -i $T1_with_lesion_fill_T1nat -s sub-${subj} -sd $output -FLAIR $FLAIR_with_lesion_fill_T1nat -FLAIRpial -openmp ${ncpu} -parallel -all

	elif [[ $wf -eq 3 ]] ; then

		recon-all -i $T1_with_lesion_fill_T1nat -s sub-${subj} -sd $output -T2 $T2_with_lesion_fill_T1nat -T2pial -FLAIR  $FLAIR_with_lesion_fill_T1nat -openmp ${ncpu} -parallel -all

	fi
	
	recon_all_pid=$!
	
	echo ${recon_all_pid}

	sleep 2
	
else
	
	echo " recon-all already done, skipping. "
	
fi


	## After recon-all is finished we need to calculate percent lesion/lobe overlap
	# need to make labels array

	lesion_lobes_report=${output}/percent_lobes_lesion_overlap_report.txt

	touch $lesion_lobes_report

	echo " Percent overlap between lesion and each lobe " >> $lesion_lobes_report

	echo " each lobe mask voxel count and volume in cmm is reported " >> $lesion_lobes_report

	echo " overlap in voxels and volume (cmm) are reported " >> $lesion_lobes_report

	declare -a labels=("RT_Frontal"  "LT_Frontal"  "RT_Temporal"  "LT_Temporal"  "RT_Parietal"  "LT_Parietal" \
	"RT_Occipital"  "LT_Occipital"  "RT_Cingulate"  "LT_Cingulate"  "RT_Insula"  "LT_Insula"  "RT_Putamen"  "LT_Putamen" \
	"RT_Caudate"  "LT_Caudate"  "RT_Thalamus"  "LT_Thalamus" "RT_Pallidum"  "LT_Pallidum"  "RT_Accumbens"  "LT_Accumbens"  "RT_Amygdala"  "LT_Amygdala" \
	"RT_Hippocampus"  "LT_Hippocampus"  "RT_PWM"  "LT_PWM");

	declare -a wm=("4001"  "3001"  "4005"  "3005"  "4006"  "3006" "4004"  "3004"  "4003"  "3003"  "4007"  "3007" \ 
	"0"  "0"  "0"  "0"  "0"  "0"  "0" "0"  "0"  "0"  "0"  "0"  "0"  "0"  "5002"  "5001");

	declare -a gm=("2001"  "1001"  "2005"  "1005"  "2006"  "1006"  "2004"  "1004"  "2003"  "1003"  "2007"  "1007" \ 
	"51"  "12"  "50"  "11"  "49"  "10"  "52"  "13"  "58"  "26"  "54"  "18"  "53"  "17"  "0"  "0");


	fs_lobes_mgz=${output}/sub-${subj}/mri/lobes_ctx_wm_fs.mgz

	fs_lobes_nii=${output}/sub-${subj}/mri/lobes_ctx_wm_fs.nii

	 labelslength=${#labels[@]}

	 wmslength=${#wm[@]}

	 gmslength=${#gm[@]}

 	fs_lobes_mark=${fs_lobes_nii}
 	search_wf_mark5=($(find ${fs_lobes_nii} -type f | grep lobes_ctx_wm_fs.nii));
	
if [[ ! $search_wf_mark5 ]]; then

	 # quick sanity check

	 if [[ ["${labelslength}" -eq "${wmslength}"] && ["${gmslength}" -eq "${wmslength}"] ]]; then

 		 echo "we're doing okay captain! " ${labelslength} " " ${wmslength} " " ${gmslength}

 	 else

 		 echo "we have a problem captain! " ${labelslength} " " ${wmslength} " " ${gmslength}

 	 fi

	 mri_annotation2label --subject sub-${subj} --sd ${output} --hemi rh --lobesStrict ${output}/sub-${subj}/label/rh.lobesStrict

	 mri_annotation2label --subject sub-${subj} --sd ${output} --hemi lh --lobesStrict ${output}/sub-${subj}/label/lh.lobesStrict

	 mri_aparc2aseg --s sub-${subj} --sd ${output} --labelwm --hypo-as-wm --rip-unknown --volmask --annot lobesStrict --o ${fs_lobes_mgz}

	 mri_convert -rl $T1_with_lesion_fill_T1nat -rt nearest $fs_lobes_mgz $fs_lobes_nii
	 
	 
 else
	 
	 echo " lobes fs image already done, skipping. "
	 
	 
 fi

	 # use for loop to read all values and indexes
	 

  	search_wf_mark6=($(find ${ROIs} -type f | grep LT_PWM_bin.nii.gz));
	
 if [[ ! $search_wf_mark6 ]]; then
	for i in {0..11}; do

	    echo "Now working on ${labels[$i]}"

	    fslmaths ${fs_lobes_nii} -thr ${gm[$i]} -uthr ${gm[$i]} ${ROIs}/${labels[$i]}_gm.nii.gz

	    fslmaths ${fs_lobes_nii} -thr ${wm[$i]} -uthr ${wm[$i]} ${ROIs}/${labels[$i]}_wm.nii.gz

	    fslmaths ${dir_rois} ${ROIs}/${labels[$i]}_gm.nii.gz -add ${ROIs}/${labels[$i]}_wm.nii.gz -bin ${ROIs}/${labels[$i]}_bin.nii.gz

	done

	i=""

	for i in {12..25}; do

	 	echo "Now working on ${labels[$i]}"

		fslmaths ${fs_lobes_nii} -thr ${gm[$i]} -uthr ${gm[$i]} -bin ${ROIs}/${labels[$i]}_bin.nii.gz

 	done

	i=""

	for i in {26..27}; do

		echo "Now working on ${labels[$i]}"

		fslmaths ${fs_lobes_nii} -thr ${wm[$i]} -uthr ${wm[$i]} -bin ${ROIs}/${labels[$i]}_bin.nii.gz

	done
	
else
	
	echo " isolating lobe labels already done, skipping to lesion overlap check"
	
fi

	i=""

	   # Now to check overlap and quantify existing overlaps
	   # we also need to calculate volume and no. of vox for each lobe out of FS
	   # also lesion volume
	   
	   # the way I'm using bc now isn't working, it just prints 0 percent, or empty space.. need to firgure this out.

	   l_vol=($(fslstats $L_mask_in_T1 -V))
	   echo " * The lesion occupies " ${l_vol[0]} " voxels in total with " ${l_vol[0]} " cmm volume. " >> $lesion_lobes_report

	   for (( i=0; i<${labelslength}; i++ )); do

		   fslmaths ${ROIs}/${labels[$i]}_bin.nii.gz -mas $L_mask_in_T1  ${overlap}/${labels[$i]}_intersect_L_mask.nii.gz

		   b=($(fslstats ${overlap}/${labels[$i]}_intersect_L_mask.nii.gz -V))
		   a=($( echo ${b[0]} | cut -c1-1))

		   vol_lobe=($(fslstats $ROIs/${labels[$i]}_bin.nii.gz -V))

		   echo " - The " ${labels[$i]} " label is " ${vol_lobe[0]} " voxels in total, with a volume of " ${vol_lobe[1]} " cmm volume. " \
			   >> $lesion_lobes_report

		   if [[ $a -ne 0 ]]; then

			   vol_ov=($(fslstats ${overlap}/${labels[$i]}_intersect_L_mask.nii.gz -V))
			   
			   ov_perc=($(echo "scale=4; (${vol_ov[1]}/${vol_lobe[1]})*100" | bc ))

			   echo " ** The lesion overlaps with the " ${labels[$i]} " in " ${vol_ov[1]} \
				   " cmm " ${ov_perc} " percent of total lobe volume " >> $lesion_lobes_report

		   else

			   echo " No overlap between the lesion and " ${labels[$i]} " lobe. " >> $lesion_lobes_report

		   fi

	   done
	   
	   end=$(date +%s)
	   
	   echo $start
	   echo $end

done


