#!/bin/bash -e
# set -x
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
#
# v0.1 - dd 09/11/2018 - alpha version
version="v1.2 - dd 23/11/2021"

# To Do
#  - fod calc msmt-5tt in stead of dhollander


### TESTING
new_style=1

# -----------------------------------  MAIN  ---------------------------------------------
# this script defines a few functions:
#  - Usage (for information to the novice user)
#  - kul_e2cl (for logging)
#
# this script will skip if some steps are already processed

kul_main_dir=`dirname "$0"`
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` performs dMRI preprocessing.

Usage:

  `basename $0` -p participant <OPT_ARGS>

Example:

  `basename $0` -p pat001 -n 6 -d "tax dhollander"

Required arguments:

	 -p:  praticipant (BIDS name of the participant)


Optional arguments:

	 -d:  dwiprep options: can be dhollander, tax and/or tournier (default = dhollander) e.g. "tax dhollander"
	 -s:  session (BIDS session)
	 -n:  number of cpu for parallelisation (default 6)
	 -b:  use Synb0-DISCO instead of topup (requires docker)
	 -e:  options to pass to eddy (default "--slm=linear --repol")
	 -r:  use reverse phase data only for topup and not for further processing
	 -m:  specify the dwi2mask method (1=hdbet, 2=b02template-ants, 3=legacy; 1=default)
	 -v:  show output from mrtrix commands (0=silent, 1=normal, 2=verbose; default=1)

Documentation:

	This script preprocesses dMRI data using MRtrix3.
	It uses input data organised in the BIDS format.
	It perfoms the following:
		1/ contactenation of different dMRI datasets accounting for differential intensity scaling using dwicat
		2/ dwidenoise
		3/ mrdegibs
		4/ dwifslpreproc, either using topup or synb0-disco
		5/ dwibiascorrect
		6/ upsampling to an isotropic resolution of 1.3 mm 
		7/ creation of dwi_mask
		8/ response estimation
		9: outputs an ADC, FA en DEC image in the for quality assurance purpose


USAGE

	exit 1
}

# set +x

# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
ncpu=6 # default if option -n is not given
silent=1 # default if option -v is not given
# Specify additional options for FSL eddy
eddy_options="--slm=linear --repol"
dwipreproc_options="dhollander"
dwi2mask_method=1

# Set required options
p_flag=0
s_flag=0
rev_only_topup=0
synb0=0
#fmap=0
verbose_level=1

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:n:d:t:e:m:v:rbv" OPT; do

		case $OPT in
		p) #participant
			p_flag=1
			participant=$OPTARG
		;;
		s) #session
			s_flag=1
			ses=$OPTARG
		;;
		n) #ncpu
			ncpu=$OPTARG
		;;
		d) #dwipreproc_options
			dwipreproc_options=$OPTARG
		;;
		e) #eddy_options
			eddy_options=$OPTARG
		;;
		m) #dwi2mask options
			dwi2mask_method=$OPTARG
		;;
		r) #rev_only topup
			rev_only_topup=1
		;;
		b) #use SynB0-DISCO
			synb0=1
		;;
		v) #verbose
			silent=$OPTARG
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


# MRTRIX verbose or not?
if [ $silent -eq 0 ] ; then
	export MRTRIX_QUIET=1
	verbose_level=0
elif [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
	verbose_level=1
elif [ $silent -eq 2 ] ; then
	export MRTRIX_QUIET=0
	verbose_level=2
fi

# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Some parallelisation
FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

d=$(date "+%Y-%m-%d_%H-%M-%S")
log=log/log_${d}.txt


# --- MAIN ----------------

# Check mrtrix3 version
if [ $mrtrix_version_revision_major -eq 2 ]; then
	mrtrix3new=0
	echo " you are using an older version of MRTrix3 $mrtrix_version_revision_major"
	echo " this is not supported. Exitting"
	exit 1
elif [ $mrtrix_version_revision_major -eq 3 ] && [ $mrtrix_version_revision_minor -lt 100 ]; then
	mrtrix3new=1
	echo " you are using a new version of MRTrix3 $mrtrix_version_revision_major $mrtrix_version_revision_minor but not the latest"
elif [ $mrtrix_version_revision_major -eq 3 ] && [ $mrtrix_version_revision_minor -gt 100 ]; then
	mrtrix3new=2
	echo " you are using the newest version of MRTrix3 $mrtrix_version_revision_major $mrtrix_version_revision_minor"
else 
	echo "cannot find correct mrtrix versions - exitting"
	exit 1
fi

# start
bids_participant=BIDS/sub-${participant}

# Either a session is given on the command line
# If not the session(s) need to be determined.
if [ $s_flag -eq 1 ]; then

	# session is given on the command line
	search_sessions=BIDS/sub-${participant}/ses-${ses}

else

	# search if any sessions exist
	search_sessions=($(find BIDS/sub-${participant} -type d | grep dwi))

fi

num_sessions=${#search_sessions[@]}

echo "  Number of BIDS sessions: $num_sessions"
echo "    notably: ${search_sessions[@]}"


# ---- BIG LOOP for processing each session
for current_session in `seq 0 $(($num_sessions-1))`; do

# set up directories
#cd $cwd
long_bids_participant=${search_sessions[$current_session]}
#echo $long_bids_participant
bids_participant=${long_bids_participant%dwi}

# Create the Directory to write preprocessed data in
preproc=${cwd}/dwiprep/sub-${participant}/$(basename $bids_participant)
#echo $preproc

# Directory to put raw mif data in
raw=${preproc}/raw

# set up preprocessing & logdirectory
mkdir -p ${preproc}/raw
mkdir -p ${preproc}/log

KUL_LOG_DIR="${preproc}/log"
mkdir -p $KUL_LOG_DIR
echo $KUL_LOG_DIR

kul_e2cl " Start processing $bids_participant" ${preproc}/${log}


# STEP 1 - CONVERSION of BIDS to MIF ---------------------------------------------

function KUL_dwiprep_convert {
# test if conversion has been done
declare -a dwi_pes
declare -a fdwi_pes
declare -a pedirs
declare -a peds
if [ ! -f ${preproc}/dwi_orig.mif ]; then

	kul_e2cl " Preparing datasets from BIDS directory..." ${preproc}/${log}

	# convert dwi
	bids_dwi_search="$bids_participant/dwi/sub-*_dwi.nii.gz"
	bids_dwi_found=$(ls $bids_dwi_search)

	number_of_bids_dwi_found=$(echo $bids_dwi_found | wc -w)

	if [ $number_of_bids_dwi_found -eq 1 ]; then #FLAG, if comparing dMRI sequences, they should not be catted

		kul_e2cl "   only 1 dwi dataset, scaling not necessary" ${preproc}/${log}
		dwi_base=${bids_dwi_found%%.*}
		mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
		-json_import ${dwi_base}.json -strides 1:3 -force \
		-clear_property comments -nthreads $ncpu ${preproc}/dwi_orig.mif
		
		mrinfo -force ${preproc}/dwi_orig.mif -export_pe_table ${raw}/dwi_orig_petable.txt
		mapfile -t dwi_pes < ${raw}/dwi_orig_petable.txt
		fdwi_pes=(${dwi_pes[0]})
		# we assume each part has only 1 pe
		pedir="${fdwi_pes[0]},${fdwi_pes[1]},${fdwi_pes[2]}"

	else

		kul_e2cl "   found $number_of_bids_dwi_found dwi datasets, checking number_of slices (and adjusting), scaling & catting" ${preproc}/${log}

		### find number of slices of the multiple datasets and correct if necessary
		dwi_i=1
		for dwi_file in $bids_dwi_found; do
			dwi_base=${dwi_file%%.*}

			# read the number of slices
			ns_dwi[dwi_i]=$(mrinfo ${dwi_base}.nii.gz -size | awk '{print $(NF-1)}')
			kul_e2cl "   dataset p${dwi_i} has ${ns_dwi[dwi_i]} as number of slices" ${preproc}/${log}
			 
			((dwi_i++))

		done

		max=10000000
		for i in "${ns_dwi[@]}"
			do
			# Update max if applicable
			if [[ "$i" -lt "$max" ]]; then
				max="$i"
			fi

		done

		# echo "Max is: $max"
		((max--))


		# convert each file and make a mean b0
		dwi_i=1
		for dwi_file in $bids_dwi_found; do
			dwi_base=${dwi_file%%.*}
			declare -a ped_${dwi_i}

			mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
			-json_import ${dwi_base}.json ${raw}/dwi_p${dwi_i}.mif -strides 1:3 -coord 2 0:${max} -force -clear_property comments -nthreads $ncpu

			#dwiextract -quiet -bzero ${raw}/dwi_p${dwi_i}.mif - | mrmath -axis 3 - mean ${raw}/b0s_p${dwi_i}.mif -force

			((dwi_i++))
		done

		# dwicat the files together
		dwicat ${raw}/dwi_p?.mif ${preproc}/dwi_orig.mif #-nocleanup 

	fi

else

	echo " Conversion has been done already... skipping to next step"

fi
}

task_in="KUL_dwiprep_convert"
KUL_task_exec $verbose_level "kul_dwiprep part 1: convert data to mif" "$KUL_LOG_DIR/1_convert_mif"

# Only keep the desired part of the dMRI
if [ $rev_only_topup -eq 1 ]; then
	

	#TO DO

		#mrinfo ${raw}/dwi_p${dwi_i}.mif -export_pe_table ${raw}/dwi_p${dwi_i}_petable.txt
		#peds[$((dwi_i-1))]=$(mrinfo ${raw}/dwi_p${dwi_i}.mif -size)
		#mapfile -t ${dwi_pes} < ${raw}/dwi_p${dwi_i}_petable.txt
		#echo $dwi_pes
		#fdwi_pes=(${dwi_pes[0]})
		#echo ${fdwi_pes[@]}
		#pedirs[$((dwi_i-1))]="${fdwi_pes[0]},${fdwi_pes[1]},${fdwi_pes[2]}"
		#echo ${pedirs[$((dwi_i-1))]}

		#for xx in ${!pedirs[@]}; do
		#
		#	if [[ ${pedirs[$xx]} == ${pedirs[0]} ]] && [[ ${peds[$xx]} -ge ${peds[0]} ]]; then
		#		echo "same pe as first part"
		#		pedir=${pedirs[$xx]}
		#		# pedinv=${ped[$xx]}
		#	elif [[ ! ${pedirs[$xx]} == ${pedirs[0]} ]] && [[ ${peds[$xx]} -lt ${peds[0]} ]]; then
		#		echo "different pe than first part"
		#		pedir=${ped[0]}
		#		pedinv=${pedirs[$xx]}
		#	elif [[ ! ${pedirs[$xx]} == ${pedirs[0]} ]] && [[ ${peds[$xx]} -gt ${peds[0]} ]]; then
		#		pedir=${ped[$xx]}
		#		pedinv=${pedirs[0]}
		#	elif [[ ${pedirs[$xx]} == ${pedirs[0]} ]] && [[ ${peds[$xx]} == ${peds[0]} ]]; then
		#		pedir=${ped[0]}
		#		pedinv=${pedirs[$xx]}
		#	fi
		#
		#done

	pedir="0,-1,0"
	dwiextract ${preproc}/dwi_orig.mif -pe ${pedir} ${preproc}/dwi_orig_norev.mif -force
	dwi_orig=dwi_orig_norev.mif

else

	dwi_orig=dwi_orig.mif

fi



# STEP 2 - DWI Preprocessing ---------------------------------------------

#echo ${preproc}
cd ${preproc}
mkdir -p dwi

# Make a descent initial mask
if [ ! -f dwi_orig_mask.nii.gz ]; then
	kul_e2cl "   Making an initial brain mask..." ${log}

	task_in1="dwiextract ${dwi_orig} dwi/bzeros.mif -bzero -force \
	&& dwiextract ${dwi_orig} dwi/nonbzeros.mif -no_bzero -force \
	&& mrcat -force -nthreads ${ncpu} dwi/bzeros.mif dwi/nonbzeros.mif dwi/rearranged_dwis.mif"
	
	if [ $mrtrix3new -eq 2 ]; then
		if [ $dwi2mask_method -eq 1 ];then
			task_in2="dwi2mask hdbet \
				dwi/rearranged_dwis.mif dwi_orig_mask.nii.gz -nthreads $ncpu -force"
		elif [ $dwi2mask_method -eq 2 ];then
			task_in2="dwi2mask b02template -software antsfull -template ${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod.nii.gz \
				${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod_brain_mask.nii.gz \
				dwi/rearranged_dwis.mif dwi_orig_mask.nii.gz -nthreads $ncpu -force"
		elif [ $dwi2mask_method -eq 3 ];then
			task_in2="dwi2mask legacy \
				dwi/rearranged_dwis.mif dwi_orig_mask.nii.gz -nthreads $ncpu -force"
		fi
		task_in="$task_in1; $task_in2"
		KUL_task_exec $verbose_level "kul_dwiprep: make an initial mask" "$KUL_LOG_DIR/1_convert_mif"
	else
		dwi2mask \
				dwi/rearranged_dwis.mif dwi_orig_mask.nii.gz -nthreads $ncpu -force
	fi
fi

# Do some qa: make FA/ADC of unprocessed images
mkdir -p qa

if [ ! -f qa/adc_orig.nii.gz ]; then

	kul_e2cl "   Calculating FA/ADC/dec..." ${log}
	dwi2tensor $dwi_orig dwi_orig_dt.mif -mask dwi_orig_mask.nii.gz -force
	tensor2metric dwi_orig_dt.mif -fa qa/fa_orig.nii.gz -force
	tensor2metric dwi_orig_dt.mif -adc qa/adc_orig.nii.gz -force

fi

# check if first 2 steps of dwi preprocessing are done
if [ ! -f dwi/degibbs.mif ] && [ ! -f dwi_preproced.mif ]; then

	kul_e2cl " Start part 1 of preprocessing: dwidenoise & mrdegibbs" ${log}

	# dwidenoise
	# kul_e2cl "   dwidenoise..." ${log}
	# STEFAN: pretty sure the mask gives a masked topup_in in kul_dwifslpreproc,... \
	#    and this get's converted to b0.nii for synb0. No eyes in there anymore and this makes the fieldmap look really bad
	# dwidenoise $dwi_orig dwi/denoise.mif -noise dwi/noiselevel.mif -mask dwi_orig_mask.nii.gz -nthreads $ncpu -force
	task_in="dwidenoise $dwi_orig dwi/denoise.mif -noise dwi/noiselevel.mif -mask dwi_orig_mask.nii.gz -nthreads $ncpu -force"
	KUL_task_exec $verbose_level "kul_dwiprep part 2: dwidenoise" "$KUL_LOG_DIR/2_dwidenoise"

	# mrdegibbs
	# kul_e2cl "   mrdegibbs..." ${log}
	task_in="mrdegibbs dwi/denoise.mif dwi/degibbs.mif -nthreads $ncpu -force; \
		rm dwi/denoise.mif"
	KUL_task_exec $verbose_level "kul_dwiprep part 2: mrdegibbs" "$KUL_LOG_DIR/2_mrdegibbs"

else

	echo "   part 1 of preprocessing has been done already... skipping to next step"

fi


# check if step 3 of dwi preprocessing is done (dwipreproc, i.e. motion and distortion correction takes very long)
if [ ! -f dwi/geomcorr.mif ]  && [ ! -f dwi_preproced.mif ]; then

	# motion and distortion correction using rpe_header
	kul_e2cl "   Start part 2 of preprocessing: dwipreproc (this takes time!)..." ${log}

	# Make the directory for the output of eddy_qc
	mkdir -p eddy_qc/raw

	# prepare for Synb0-disco
	if [ $synb0 -eq 1 ]; then


		if [ $new_style -eq 1 ];then
			
			# run KUL_synb0 first and then use the new dwifslpreproc -topupfield option
			KUL_synb0.sh -p $participant

		else 

			# find T1
			bids_T1_found=($(find $cwd/$bids_participant/anat -type f -name "*T1w.nii.gz" ! -name "*gadolinium*")) 
			number_of_bids_T1_found=${#bids_T1_found[@]}
			if [ $number_of_bids_T1_found -gt 1 ]; then
				kul_e2cl "   more than 1 T1 dataset, using first only for Synb0-disco" ${preproc}/${log}
			fi
			#echo $bids_T1_found
			Synb0_T1=${bids_T1_found[0]}
			echo "The used T1 for sSynb0-disco is $Synb0_T1"

		
		fi
	fi

	# prepare for fmap
	#if [ $fmap -eq 1 ]; then
	#	# find fmap
	#	bids_fmap_search="$cwd/$bids_participant/fmap/sub-*_epi.nii.gz"
	#	echo $bids_fmap_search
	#	bids_fmap_found=$(ls $bids_fmap_search)
	#	number_of_bids_fmap_found=$(echo $bids_fmap_found | wc -w)
	#	if [ $number_of_bids_fmap_found -lt 1 ]; then
	#		kul_e2cl "   more than 1 fmap dataset, using first only for topup" ${preproc}/${log}
	#	fi
	#	#echo $bids_fmap_found
	#	fmap_epi=${bids_fmap_found[0]}
	#	echo "The used fmap for topup is $fmap_epi"
	#	fmap_base=${fmap_epi%%.*}
	#	mrconvert ${fmap_base}.nii.gz -json_import ${fmap_base}.json raw/dwi_reverse_phase.mif -strides 1:3 -force -clear_property comments -nthreads $ncpu
	#fi

	# prepare eddy_options
	#
	echo "eddy_options: $eddy_options"
	full_eddy_options="--cnr_maps --residuals "${eddy_options}
	echo "full_eddy_options: $full_eddy_options"


	# PREPARE FOR TOPUP and EDDY, but allways use rpe_header, since we start from BIDS format
	regular_dwipreproc=1

	# read the pe table of the b0s of dwi_orig.mif
	IFS=$'\n'
	pe=($(dwiextract dwi_orig.mif -bzero - | mrinfo -petable -))
	#echo $pe
	# count how many b0s there are
	n_pe=$(echo ${#pe[@]})
	#echo $n_pe


	# in case there is only 1 b0max
	if [ $n_pe -gt 4 ]; then
		echo "Prepare for topup: more than 4 b0s, now checking pe_scheme"

		# extract first b0
		dwiextract dwi_orig.mif -bzero - | mrconvert - -coord 3 0 raw/b0s_pe0.mif -force
		# get the pe_scheme of the first b0
		previous_pe=$(echo ${pe[0]})

		# read over the following b0s, and only keep 1 with a new b0 scheme

		for i in `seq 1 $(($n_pe-1))`; do

			current_pe=$(echo ${pe[$i]})

			if [ $previous_pe = $current_pe ]; then
				echo previous_pe=$previous_pe, current_pe=$current_pe
				echo "same pe scheme, skip"

			else

				echo "Prepare for topup: more than 5 b0s, but some have different pe_scheme"
				regular_dwipreproc=0

				echo previous_pe=$previous_pe, current_pe=$current_pe
				echo "new pe_scheme, convert"
				dwiextract dwi_orig.mif -bzero - | mrconvert - -coord 3 $i raw/b0s_pe${i}.mif -force
				break 

			fi

			previous_pe=$current_pe

		done

		mrcat raw/b0s_pe*.mif dwi/se_epi_for_topup.mif -force

	fi

	echo "regular_dwipreproc: $regular_dwipreproc"
	echo "mrtrix3new: $mrtrix3new"
	echo "synb0: $synb0"
	echo "fmap: $fmap"
	echo "rev_only_topup: $rev_only_topup"

	if [ $new_style -eq 1 ];then

		if [ $synb0 -eq 0 ]; then

			# note: maybe add -eddy_mask
			dwifslpreproc_option="-se_epi dwi/se_epi_for_topup.mif -align_seepi"

		else
			
			dwifslpreproc_option="-topup_files ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/synb0/topup_"

		fi

		task_in="dwifslpreproc ${dwifslpreproc_option} -rpe_header \
				-eddyqc_all eddy_qc/raw -eddy_options \"${full_eddy_options} \" -force -nthreads $ncpu -nocleanup \
				dwi/degibbs.mif dwi/geomcorr.mif"
		echo $task_in
		KUL_task_exec $verbose_level "kul_dwiprep part 3: motion and distortion correction" "$KUL_LOG_DIR/3_dwifslpreproc"
		wait

	else

		if [ $regular_dwipreproc -eq 1 ]; then


			if [ $synb0 -eq 0 ]; then
				if [ $fmap -eq 0 ]; then
					kul_dwifslpreproc dwi/degibbs.mif dwi/geomcorr.mif -rpe_header \
						-eddyqc_all eddy_qc/raw -eddy_options "${full_eddy_options} " -force -nthreads $ncpu -nocleanup
				else
					mrcat raw/b0s_pe*.mif raw/se_epi_for_topup.mif -force
					kul_dwifslpreproc -se_epi raw/se_epi_for_topup.mif -align_seepi dwi/degibbs.mif dwi/geomcorr.mif -rpe_header \
						-eddyqc_all eddy_qc/raw -eddy_options "${full_eddy_options} " -force -nthreads $ncpu -nocleanup
				fi
			else
				kul_dwifslpreproc dwi/degibbs.mif dwi/geomcorr.mif \
					-synb0_disco_T1 "$Synb0_T1" \
					-rpe_header -eddyqc_all eddy_qc/raw -eddy_options "${full_eddy_options} " -force -nthreads $ncpu -nocleanup
			fi


		else

			if [ $synb0 -eq 1 ]; then
				echo "EXPERIMENTAL: Not tested well; contact Stefan to fix scipt KUL_dwiprep.sh synb0=1 case, but regular_dwiprproc=0, lines 485 on"
				kul_dwifslpreproc dwi/degibbs.mif dwi/geomcorr.mif \
				-synb0_disco_T1 "$Synb0_T1" \
				-rpe_header -eddyqc_all eddy_qc/raw -eddy_options "${full_eddy_options} " -force -nthreads $ncpu -nocleanup
			
			
			
			else

				# concat all b0 with different pe_schemes
				mrcat raw/b0s_pe*.mif raw/se_epi_for_topup.mif -force

				echo $rev_only_topup

				if [ $rev_only_topup -eq 0 ]; then

					kul_dwifslpreproc -se_epi raw/se_epi_for_topup.mif -align_seepi dwi/degibbs.mif dwi/geomcorr.mif -rpe_header \
					-eddyqc_all eddy_qc/raw -eddy_options "${full_eddy_options} " -force -nthreads $ncpu -nocleanup

				else

					#kul_dwifslpreproc -se_epi raw/se_epi_for_topup.mif -align_seepi dwi/degibbs.mif dwi/geomcorr.mif -rpe_pair \
					#-eddyqc_all eddy_qc/raw -eddy_options "${full_eddy_options} " -force -nthreads $ncpu -nocleanup
					kul_dwifslpreproc -se_epi raw/se_epi_for_topup.mif -align_seepi dwi/degibbs.mif dwi/geomcorr.mif -rpe_header \
					-eddyqc_all eddy_qc/raw -eddy_options "${full_eddy_options} " -force -nthreads $ncpu -nocleanup
				fi

			fi
		
		fi

	fi		


	temp_dir=$(ls -d *dwifslpreproc*)


	# create an intermediate mask of the dwi data
	kul_e2cl "    creating intermediate mask of the dwi data..." ${log}

	if [ $mrtrix3new -eq 2 ]; then
      dwiextract dwi/geomcorr.mif dwi/geomcorr_bzeros.mif -bzero -force \
	    && dwiextract dwi/geomcorr.mif dwi/geomcorr_nonbzeros.mif -no_bzero -force \
	    && mrcat -force -nthreads ${ncpu} dwi/geomcorr_bzeros.mif dwi/geomcorr_nonbzeros.mif dwi/rearranged_geomcorr_dwis.mif
		if [ $dwi2mask_method -eq 1 ];then
			dwi2mask hdbet \
				dwi/rearranged_geomcorr_dwis.mif dwi/dwi_intermediate_mask.nii.gz -nthreads $ncpu -force
		elif [ $dwi2mask_method -eq 2 ];then
			dwi2mask b02template -software antsfull -template ${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod.nii.gz \
				${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod_brain_mask.nii.gz \
				dwi/rearranged_geomcorr_dwis.mif dwi/dwi_intermediate_mask.nii.gz -nthreads $ncpu -force
		elif [ $dwi2mask_method -eq 3 ];then
			dwi2mask legacy \
				dwi/rearranged_geomcorr_dwis.mif dwi/dwi_intermediate_mask.nii.gz -nthreads $ncpu -force
		fi
	else
		dwi2mask dwi/rearranged_geomcorr_dwis.mif dwi_intermediate_mask.nii.gz -nthreads $ncpu -force
	fi

	# check id eddy_quad is available
	#machine_type=$(uname)
	#echo $machine_type
	test_eddy_quad=$(command -v eddy_quad)
	echo $test_eddy_quad

	if [ $test_eddy_quad = "" ]; then

		echo "Eddy_quad skipped (which was not found in the path)"

	else

		# Let's run eddy_quad
		rm -rf eddy_qc/quad
		kul_e2cl "   running eddy_quad..." ${log}
		eddy_quad $temp_dir/dwi_post_eddy --eddyIdx $temp_dir/eddy_indices.txt \
			--eddyParams $temp_dir/dwi_post_eddy.eddy_parameters --mask $temp_dir/eddy_mask.nii \
			--bvals $temp_dir/bvals --bvecs $temp_dir/bvecs --output-dir eddy_qc/quad --verbose

		# make an mriqc/fmriprep style report (i.e. just link qc.pdf into main dwiprep directory)
		echo $cwd
		echo $preproc
		if [[ -z ${ses} ]]; then
			ln -s $cwd/${preproc}/eddy_qc/quad/qc.pdf $cwd/${preproc}/../sub-${participant}.pdf &
		else
			ln -s $cwd/${preproc}/eddy_qc/quad/qc.pdf $cwd/${preproc}/../sub-${participant}_ses-${ses}.pdf &
		fi

	fi

	# clean-up the above dwipreproc temporary directory
	# rm -rf $temp_dir

else

	echo "   part 2 of preprocessing has been done already... skipping to next step"

fi


# check if next 4 steps of dwi preprocessing are done
if [ ! -f dwi_preproced.mif ]; then

	kul_e2cl " Start part 3 of preprocessing: dwibiascorrect, upsampling & creation of a final dwi_mask" ${log}

	# bias field correction
	kul_e2cl "    dwibiascorrect" ${log}
	dwibiascorrect ants dwi/geomcorr.mif dwi/biascorr.mif -bias dwi/biasfield.mif -nthreads $ncpu -force -mask dwi/dwi_intermediate_mask.nii.gz


	# upsample the images
	kul_e2cl "    upsampling resolution..." ${log}
	mrgrid -nthreads $ncpu -force -axis 1 5,5 dwi/biascorr.mif crop - | mrgrid -axis 1 5,5 -force - pad - | mrgrid -voxel 1.3 -force - regrid dwi/upsampled.mif 

	rm dwi/biascorr.mif

	# copy to main directory for subsequent processing
	kul_e2cl "    saving..." ${log}
	mrconvert dwi/upsampled.mif dwi_preproced.mif -set_property comments "Preprocessed dMRI data." -nthreads $ncpu -force
	rm dwi/upsampled.mif

	# create a final mask of the dwi data
	kul_e2cl "    creating mask of the dwi data..." ${log}
	if [ $mrtrix3new -eq 2 ]; then
		dwiextract dwi_preproced.mif dwi/preproced_bzeros.mif -bzero -force \
		&& dwiextract dwi_preproced.mif dwi/preproced_nonbzeros.mif -no_bzero -force \
		&& mrcat -force -nthreads ${ncpu} dwi/preproced_bzeros.mif dwi/preproced_nonbzeros.mif dwi/rearranged_preproced_dwis.mif
		if [ $dwi2mask_method -eq 1 ];then
			dwi2mask hdbet \
				dwi/rearranged_preproced_dwis.mif dwi_mask.nii.gz -nthreads $ncpu -force
		elif [ $dwi2mask_method -eq 2 ];then
			dwi2mask b02template -software antsfull -template ${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod.nii.gz \
				${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod_brain_mask.nii.gz \
				dwi/rearranged_preproced_dwis.mif dwi_mask.nii.gz -nthreads $ncpu -force
		elif [ $dwi2mask_method -eq 3 ];then
			dwi2mask legacy \
				dwi/rearranged_preproced_dwis.mif dwi_mask.nii.gz -nthreads $ncpu -force
		fi
	else
		dwi2mask dwi/rearranged_preproced_dwis.mif dwi_mask.nii.gz -nthreads $ncpu -force
	fi

	# create mean b0 of the dwi data
	kul_e2cl "    creating mean b0 of the dwi data ..." ${log}
	dwiextract -quiet dwi_preproced.mif -bzero - | mrmath -axis 3 - mean dwi_b0.nii.gz -force

else

	echo "   part 3 of preprocessing has been done already... skipping to next step"

fi



# STEP 3 - RESPONSE ---------------------------------------------
mkdir -p response
# response function estimation (we compute following algorithms: tournier, tax and dhollander)

echo " Using dwipreproc_options: $dwipreproc_options"

if [[ $dwipreproc_options == *"dhollander"* ]]; then

	if [ ! -f response/wm_response.txt ]; then
		kul_e2cl "   Calculating dhollander dwi2response..." ${log}
		dwi2response dhollander dwi_preproced.mif response/dhollander_wm_response.txt -mask dwi_mask.nii.gz \
		response/dhollander_gm_response.txt response/dhollander_csf_response.txt -nthreads $ncpu -force

	else

		echo " dwi2response dhollander already done, skipping..."

	fi

	if [ ! -f response/dhollander_wmfod.mif ]; then
		kul_e2cl "   Calculating dhollander dwi2fod & normalising it..." ${log}

		dwi2fod msmt_csd dwi_preproced.mif response/dhollander_wm_response.txt response/dhollander_wmfod.mif \
		response/dhollander_gm_response.txt response/dhollander_gm.mif \
		response/dhollander_csf_response.txt response/dhollander_csf.mif -mask dwi_mask.nii.gz -force -nthreads $ncpu

		dwi2fod msmt_csd dwi_preproced.mif response/dhollander_wm_response.txt response/dhollander_wmfod_noGM.mif \
		response/dhollander_csf_response.txt response/dhollander_csf_noGM.mif -mask dwi_mask.nii.gz -force -nthreads $ncpu

		#mtnormalise response/dhollander_wmfod.mif response/dhollander_wmfod_norm.mif \
		#response/dhollander_gm.mif response/dhollander_gm_norm.mif \
		#response/dhollander_csf.mif response/dhollander_csf_norm.mif -mask dwi_mask.nii.gz -force -nthreads $ncpu

		#mtnormalise response/dhollander_wmfod_noGM.mif response/dhollander_wmfod_norm_noGM.mif \
		#response/dhollander_csf_noGM.mif response/dhollander_csf_norm_noGM.mif -mask dwi_mask.nii.gz -force -nthreads $ncpu

	else

		echo " dwi2fod dhollander already done, skipping..."

	fi

fi

if [[ $dwipreproc_options == *"tax"* ]]; then

	if [ ! -f response/tax_response.txt ]; then
		kul_e2cl "   Calculating tax dwi2response..." ${log}
		dwi2response tax dwi_preproced.mif response/tax_response.txt -nthreads $ncpu -force

	else

		echo " dwi2response tax already done, skipping..."

	fi

	if [ ! -f response/tax_wmfod.mif ]; then
		kul_e2cl "   Calculating tax dwi2fod..." ${log}
		dwi2fod csd dwi_preproced.mif response/tax_response.txt response/tax_wmfod.mif  \
		-mask dwi_mask.nii.gz -force -nthreads $ncpu -mask dwi_mask.nii.gz

	else

		echo " dwi2fod tax already done, skipping..."

	fi

fi

if [[ $dwipreproc_options == *"tournier"* ]]; then

	if [ ! -f response/tournier_response.txt ]; then
		kul_e2cl "   Calculating tournier dwi2response..." ${log}
		dwi2response tournier dwi_preproced.mif response/tournier_response.txt -nthreads $ncpu -force

	else

		echo " dwi2response already done, skipping..."

	fi

	if [ ! -f response/tournier_wmfod.mif ]; then
		kul_e2cl "   Calculating tournier dwi2fod..." ${log}
		dwi2fod csd dwi_preproced.mif response/tournier_response.txt response/tournier_wmfod.mif  \
		-mask dwi_mask.nii.gz -force -nthreads $ncpu

	else

		echo " dwi2fod already done, skipping..."

	fi

fi

# STEP 4 - DO QA ---------------------------------------------
# Make an FA/dec image

mkdir -p qa

if [ ! -f qa/dec.mif ]; then

	kul_e2cl "   Calculating FA/ADC/dec..." ${log}
	dwi2tensor dwi_preproced.mif dwi_dt.mif -force -mask dwi_mask.nii.gz
	tensor2metric dwi_dt.mif -fa qa/fa.nii.gz -mask dwi_mask.nii.gz -force
	tensor2metric dwi_dt.mif -adc qa/adc.nii.gz -mask dwi_mask.nii.gz -force

	if [[ $dwipreproc_options == *"tournier"* ]]; then

		fod2dec response/tournier_wmfod.mif qa/tournier_dec.mif -force -mask dwi_mask.nii.gz
	fi

	if [[ $dwipreproc_options == *"tax"* ]]; then
		fod2dec response/tax_wmfod.mif qa/tax_dec.mif -force -mask dwi_mask.nii.gz
	fi

	if [[ $dwipreproc_options == *"dhollander"* ]]; then
		fod2dec response/dhollander_wmfod.mif qa/dhollander_dec.mif -force -mask dwi_mask.nii.gz
		#fod2dec response/dhollander_wmfod_norm.mif qa/dhollander_dec_norm.mif -force
	fi

	#mrconvert dwi/noiselevel.mif qa/noiselevel.nii.gz

fi


# We finished processing current session
# write a "done" log file for this session
echo "done" > log/${current_session}_done.log


# STEP 5 - CLEANUP - here we clean up (large) temporary files
#rm -fr dwi/degibbs.mif
#rm -rf dwi/geomcorr.mif
#rm -rf raw

echo " Finished processing session $bids_participant"


# ---- END of the BIG loop over sessions
done

# write a file to indicate that dwiprep runned succesfully
#   his file will be checked by KUL_preproc_all
#   dwiprep_file_to_check=dwiprep/sub-${BIDS_participant}/dwiprep_is_done.log

echo "done" > ../dwiprep_is_done.log


kul_e2cl "Finished " ${log}
