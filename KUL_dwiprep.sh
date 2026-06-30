#!/bin/bash
# Bash shell script to process diffusion & structural 3D-T1w MRI data
#
# Requires Mrtrix3, FSL, ants
#
# @ Stefan Sunaert - UZ/KUL - prof.sunaert@gmail.com
# @ Ahmed Radwan - KUL - radwanphd@gmail.com
#
# v0.1 - dd 09/11/2018 - created
version="v1.3 - dd 27/11/2021"

# temp var
mb=2

# To Do
#  - fod calc msmt-5tt in stead of dhollander

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd, mrtrix3new & $log_dir is made in main_functions

# --------------------------------------------------------------------------------------------------------
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
	 -c:  use shard-recon instead of eddy
	 -e:  options to pass to eddy (default "--slm=linear --repol")
	 -v:  show output from mrtrix commands (0=silent, 1=normal, 2=verbose; default=1)
	 -r:  use reverse phase data only for topup and not for further processing
	 -m:  specify the dwi2mask method (1=hdbet, 2=b02template-ants, 3=legacy; 3=default)


Documentation:

	This script preprocesses dMRI data using MRtrix3.
	It uses input data organised in the BIDS format.
	It perfoms the following:
		1/ contactenation of different dMRI datasets accounting for differential intensity scaling using dwicat
		2/ dwidenoise
		3/ mrdegibs
		4/ dwifslpreproc or shard-recon, either using topup or synb0-disco
		5/ dwibiascorrects
		6/ upsampling to an isotropic resolution of 1.3 mm 
		7/ creation of dwi_mask
		8/ response estimation
		9: outputs an ADC, FA en DEC image in the for quality assurance purpose


USAGE

	exit 1
}

# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
ncpu=6 # default if option -n is not given
verbose_level=1
local_verbose_level=1 # default if option -v is not given
rev_only_topup=0
synb0=0
total_errorcount=0
dwipreproc_options="dhollander"
dwi2mask_method=3
shard=0

# Specify additional options for FSL eddy
eddy_options="--slm=linear --repol"


# Set required options
p_flag=0
s_flag=0
v_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:n:d:e:m:v:rbc" OPT; do

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
		v) #verbose
			local_verbose_level=$OPTARG
		;;
		r) #rev_only topup
			rev_only_topup=1
		;;
		b) #use SynB0-DISCO
			synb0=1
		;;
		c) #shard
			shard=1
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

# MRTRIX and others verbose or not?
if [ $local_verbose_level -eq 0 ] ; then
	export MRTRIX_QUIET=1
	verbose_level=0
elif [ $local_verbose_level -eq 1 ] ; then
	export MRTRIX_QUIET=1
	verbose_level=1
elif [ $local_verbose_level -eq 2 ] ; then
	export MRTRIX_QUIET=0
	verbose_level=2
fi


# REST OF SETTINGS ---

# timestamp
start=$(date +%s)

# Some parallelisation
FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

#d=$(date "+%Y-%m-%d_%H-%M-%S")
#log=log/log_${d}.txt

if [[ $synb0 -eq 1 ]] && [[ $mrtrix3new -lt 2 ]]; then
	kul_echo "Synb0 usage needs a newer version of MRTrix3 than the one you have installed. Please update to 3.0.3-xxx > 100"
	exit 1
fi


### Functions ---------------------------------------
function KUL_dwiprep_convert {
# test if conversion has been done
if [ ! -f ${preproc}/dwi_orig.mif ]; then

	kul_echo " Preparing datasets from BIDS directory..."
	#echo "preproc: $preproc"

	if [ $number_of_bids_dwi_found -eq 1 ]; then #FLAG, if comparing dMRI sequences, they should not be catted

		kul_echo "   only 1 dwi dataset, scaling not necessary"
		dwi_base=${bids_dwi_found%%.*}
		mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
		-json_import ${dwi_base}.json -strides 1:3 -force \
		-clear_property comments -nthreads $ncpu ${preproc}/dwi_orig.mif
		
	else

		kul_echo "   found $number_of_bids_dwi_found dwi datasets, checking number_of slices (and adjusting), scaling & catting"

		### find number of slices of the multiple datasets and correct if necessary
		dwi_i=1
		for dwi_file in $bids_dwi_found; do
			dwi_base=${dwi_file%%.*}

			# read the number of slices
			ns_dwi[dwi_i]=$(mrinfo ${dwi_base}.nii.gz -size | awk '{print $(NF-1)}')
			kul_echo "   dataset p${dwi_i} has ${ns_dwi[dwi_i]} as number of slices"
			 
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
			#echo "dwi_base: $dwi_base"
			
			kul_echo "   converting ${dwi_base}.nii.gz to mif"
			mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
			-json_import ${dwi_base}.json ${raw}/dwi_p${dwi_i}.mif -strides 1:3 -coord 2 0:${max} \
			-force -clear_property comments -nthreads $ncpu

			((dwi_i++))
		done

		# dwicat the files together
		kul_echo "   performing dwicat"
		dwicat ${raw}/dwi_p?.mif ${preproc}/dwi_orig.mif #-nocleanup 

	fi

else

	kul_echo " Conversion has been done already... skipping to next step"

fi
}

function kul_dwi2mask {

	task_in1="dwiextract ${dwi2mask_image_in} dwi/bzeros.mif -bzero -force \
	&& dwiextract ${dwi2mask_image_in} dwi/nonbzeros.mif -no_bzero -force \
	&& mrcat -force -nthreads ${ncpu} dwi/bzeros.mif dwi/nonbzeros.mif dwi/rearranged_dwis.mif"
	
	if [ $mrtrix3new -eq 2 ]; then
		if [ $dwi2mask_method -eq 1 ];then
			task_in2="dwi2mask hdbet \
				dwi/rearranged_dwis.mif ${dwi2mask_mask_out} -nthreads $ncpu -force"
		elif [ $dwi2mask_method -eq 2 ];then
			task_in2="dwi2mask b02template -software antsfull -template ${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod.nii.gz \
				${kul_main_dir}/atlasses/Temp_4_KUL_dwiprep/UKBB_fMRI_mod_brain_mask.nii.gz \
				dwi/rearranged_dwis.mif ${dwi2mask_mask_out} -nthreads $ncpu -force"
		elif [ $dwi2mask_method -eq 3 ];then
			task_in2="dwi2mask legacy \
				dwi/rearranged_dwis.mif ${dwi2mask_mask_out} -nthreads $ncpu -force"
		fi
	else
		task_in2="dwi2mask \
				dwi/rearranged_dwis.mif ${dwi2mask_mask_out} -nthreads $ncpu -force"
	fi
	task_in="$task_in1; $task_in2"
	KUL_task_exec $verbose_level "${dwi2mask_message}" "${dwi2mask_logfile}"
	rm -f dwi/rearranged_dwis.mi
	
}

function kul_mrview_figure {

	capture="-capture.folder $capture_dir -capture.prefix $capture_prefix -noannotations -capture.grab"
	cmd="mrview -load $underlay 
		$mode_plane \
		$overlay \
		$capture \
		-force \
		$mrview_exit"
	echo $cmd
	eval $cmd

}

# --- MAIN ----------------
# start
bids_participant=BIDS/sub-${participant}

# Either a session is given on the command line
# If not the session(s) need to be determined.
if [ $s_flag -eq 1 ]; then

	# session is given on the command line
	search_sessions=${cwd}/BIDS/sub-${participant}/ses-${ses}

else

	# search if any sessions exist
	search_sessions=($(find ${cwd}/BIDS/sub-${participant} -type d | grep dwi))

fi

num_sessions=${#search_sessions[@]}

kul_echo "  Number of BIDS sessions: $num_sessions"
kul_echo "    notably: ${search_sessions[@]}"


# ---- BIG LOOP for processing each session
for current_session in `seq 0 $(($num_sessions-1))`; do

	# set up directories
	long_bids_participant=${search_sessions[$current_session]}
	#echo $long_bids_participant
	bids_participant=${long_bids_participant%dwi}

	echo "long_bids_participant: $long_bids_participant"
	echo "bids_participant: $bids_participant"

	# Create the Directory to write preprocessed data in
	preproc=${cwd}/dwiprep/sub-${participant}/$(basename $bids_participant)
	echo "preproc: $preproc"

	# Directory to put raw mif data in
	raw=${preproc}/raw

	# set up preprocessing & logdirectory
	mkdir -p ${preproc}/raw
	mkdir -p ${preproc}/log

	KUL_LOG_DIR="${preproc}/log"
	mkdir -p $KUL_LOG_DIR
	#echo $KUL_LOG_DIR
	d=$(date "+%Y-%m-%d_%H-%M-%S")
	#log=log/log_${d}.txt
	log=${KUL_LOG_DIR}/${script}_${d}.log

	kul_echo " Start processing $bids_participant"
	cd ${preproc}

	# STEP 1 - CONVERSION of BIDS to MIF ---------------------------------------------
	bids_dwi_search="$bids_participant/dwi/sub-*_dwi.nii.gz"
	bids_dwi_found=$(ls $bids_dwi_search)
	number_of_bids_dwi_found=$(echo $bids_dwi_found | wc -w)
	task_in="KUL_dwiprep_convert"
	KUL_task_exec $verbose_level "kul_dwiprep part 1: convert data to mif" "1_convert_mif"

	echo "number_of_bids_dwi_found: $number_of_bids_dwi_found"

	# Only keep the desired part of the dMRI
	# do this if rev_only_topup -eq 1, but not if there is only 1 dwi file (which means there is no rev phase)
	if [ $rev_only_topup -eq 1 ] && [ $number_of_bids_dwi_found -gt 1 ]; then

		if [ ! -f dwi_orig_norev.mif ]; then

			# Try to detect the phase encoding (pe) direction of the reverse phase images automatically
			# Assumed: this is the pe direction with the least images
			pe_directions_present=($(mrinfo ${preproc}/dwi_orig.mif -petable | sort | uniq -c))
			pe_first_count=${pe_directions_present[0]}
			pe_second_count=${pe_directions_present[5]}
			echo "preproc: $preproc"
			echo "pe_directions_present: ${pe_directions_present[@]}"
			echo "pe_first_count: $pe_first_count"
			echo "pe_second_count: $pe_second_count"
			if [ $pe_first_count -gt $pe_second_count ]; then
				pe_direction_to_keep="${pe_directions_present[1]},${pe_directions_present[2]},${pe_directions_present[3]}"
			else
				pe_direction_to_keep="${pe_directions_present[6]},${pe_directions_present[7]},${pe_directions_present[8]}"
			fi
			kul_echo "keeping $pe_direction_to_keep phase encoding direction"

			dwiextract ${preproc}/dwi_orig.mif -pe ${pe_direction_to_keep} ${preproc}/dwi_orig_norev.mif -force
			dwi_orig=dwi_orig_norev.mif

		fi 

	else

		dwi_orig=dwi_orig.mif

	fi


	# STEP 2 - DWI Preprocessing ---------------------------------------------

	mkdir -p dwi

	# Make a descent initial mask
	if [ ! -f dwi/dwi_orig_mask.nii.gz ]; then
		kul_echo "Making an initial brain mask..."
		dwi2mask_image_in=${dwi_orig}
		dwi2mask_mask_out="dwi/dwi_orig_mask.nii.gz"
		dwi2mask_message="kul_dwiprep part1: make an initial mask"
		dwi2mask_logfile="1_convert_mif"
		kul_dwi2mask
	fi

	# Do some qa: make FA/ADC of unprocessed images
	mkdir -p qa

	if [ ! -f qa/adc_orig.nii.gz ]; then

		kul_echo "Calculating initial FA/ADC/dec..."
		task_in="dwi2tensor $dwi_orig dwi_orig_dt.mif -mask dwi/dwi_orig_mask.nii.gz -force; \
			tensor2metric dwi_orig_dt.mif -fa qa/fa_orig.nii.gz -force; \
			tensor2metric dwi_orig_dt.mif -adc qa/adc_orig.nii.gz -force"
		KUL_task_exec $verbose_level "kul_dwiprep part 1: initial /ADC/DEC" "1_qa"
	fi

	# check if first 2 steps of dwi preprocessing are done
	if [ ! -f dwi/degibbs.mif ] && [ ! -f dwi_preproced.mif ]; then

		kul_echo "Start part 1 of preprocessing: dwidenoise & mrdegibbs"

		# dwidenoise
		kul_echo "dwidenoise..."
		# STEFAN: pretty sure the mask gives a masked topup_in in kul_dwifslpreproc,... \
		#    and this get's converted to b0.nii for synb0. No eyes in there anymore and this makes the fieldmap look really bad
		# dwidenoise $dwi_orig dwi/denoise.mif -noise dwi/noiselevel.mif -mask dwi/dwi_orig_mask.nii.gz -nthreads $ncpu -force
		task_in="dwidenoise $dwi_orig dwi/denoise.mif -noise dwi/noiselevel.mif -nthreads $ncpu -force"
		KUL_task_exec $verbose_level "kul_dwiprep part 2: dwidenoise" "2_dwidenoise"

		# mrdegibbs
		kul_echo "mrdegibbs..."
		task_in="mrdegibbs dwi/denoise.mif dwi/degibbs.mif -nthreads $ncpu -force; \
			rm dwi/denoise.mif"
		KUL_task_exec $verbose_level "kul_dwiprep part 2: mrdegibbs" "2_mrdegibbs"

	else

		kul_echo "part 1 of preprocessing has been done already... skipping to next step"

	fi


	# check if step 3 of dwi preprocessing is done (dwipreproc, i.e. motion and distortion correction takes very long)
	if [ ! -f dwi/geomcorr.mif ]  && [ ! -f dwi_preproced.mif ]; then

		# motion and distortion correction using rpe_header
		kul_echo "Start part 2 of preprocessing: dwipreproc (this takes time!)..."

		# prepare for Synb0-disco
		if [ $synb0 -eq 1 ]; then
				
			# run KUL_synb0 first and then use the new dwifslpreproc -topupfield option
			cd ${cwd}
			if [ $s_flag -eq 1 ]; then
				synb0_ses="-s $ses"
			else
				synb0_ses=""
			fi
			task_in="KUL_synb0.sh -p $participant $synb0_ses"
			KUL_task_exec $verbose_level "kul_dwiprep part 3: synb0" "3_synb0"
			cd ${preproc}
		fi

		# prepare eddy_options
		#echo "eddy_options: $eddy_options"
		full_eddy_options="--cnr_maps --residuals "${eddy_options}
		#echo "full_eddy_options: $full_eddy_options"


		# PREPARE FOR TOPUP and EDDY, but allways use rpe_header, since we start from BIDS format

		# read the pe table of the b0s of dwi_orig.mif
		IFS=$'\n'
		pe=($(dwiextract dwi_orig.mif -bzero - | mrinfo -petable -))
		echo "pe: $pe"
		# count how many b0s there are
		n_pe=$(echo ${#pe[@]})
		echo "n_pe: $n_pe"
		

		# in case there is a reverse phase information
		if [ $n_pe -gt 1 ]; then
			kul_echo "Prepare for topup: checking pe_scheme"

			# extract first b0
			dwiextract dwi_orig.mif -bzero - | mrconvert - -coord 3 0 raw/b0s_pe0.mif -force
			# get the pe_scheme of the first b0
			previous_pe=$(echo ${pe[0]})
			
			# read over the following b0s, and only keep 1 with a new b0 scheme

			for i in `seq 1 $(($n_pe-1))`; do

				current_pe=$(echo ${pe[$i]})
				echo "current_pe: $curre"
				if [ $previous_pe = $current_pe ]; then
					kul_echo "previous_pe=$previous_pe, current_pe=$current_pe"
					kul_echo "same pe scheme, skip"
				else
					kul_echo "previous_pe=$previous_pe, current_pe=$current_pe"
					kul_echo "new pe_scheme found, converting"
					dwiextract dwi_orig.mif -bzero - | mrconvert - -coord 3 $i raw/b0s_pe${i}.mif -force
					break 
				fi

				previous_pe=$current_pe

			done

			mrcat raw/b0s_pe*.mif dwi/se_epi_for_topup.mif -force

			# if shard-recon is the motion-correction algo, and there is reverse phase info, run topup first
			if [ $shard -eq 1 ] && [ $synb0 -eq 0 ]; then

				kul_echo "running topup for shard-recon"
				mkdir -p shard
				mrconvert dwi/se_epi_for_topup.mif shard/topup_in.nii -strides -1,+2,+3,+4 -export_pe_table shard/topup_datain.txt -force

				# Check if the input has an uneven number fo slices, than use another topup config.
				# see https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=FSL;899f842a.2008
				topupin_num_slices=$(mrinfo -size shard/topup_in.nii | cut -d' ' -f 3)
				kul_echo "Number of slices in topupin data: $topupin_num_slices"

				if [ $((topupin_num_slices%2)) -eq 1 ]; then
				
					topup_cfg=b02b0_1.cnf

				else

					topup_cfg=b02b0.cnf

				fi

				task_in="topup --imain=shard/topup_in.nii \
					--datain=shard/topup_datain.txt \
					--out=shard/field \
					--fout=shard/fieldmap.nii.gz \
					--config=$topup_cfg --verbose"
				#	--subsamp=1,1,1,1,1,1,1,1,1 \
            	#	--miter=10,10,10,10,10,20,20,30,30 \
            	#	--lambda=0.00033,0.000067,0.0000067,0.000001,0.00000033,0.000000033,0.0000000033,0.000000000033,0.00000000000067"

				KUL_task_exec $verbose_level "kul_dwiprep part 3: topup for shard-recon" "3_topup_shard"


			elif [ $shard -eq 1 ] && [ $synb0 -eq 1 ]; then
				
				#  in case of synb0
				if [ $s_flag -eq 1 ];then
					dwifslprep_ses="ses-${ses}/"
				else
					dwifslprep_ses=""
				fi

				#echo ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/${dwifslprep_ses}synb0/topup_fieldmap.nii.gz
				cp ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/${dwifslprep_ses}synb0/topup_fieldmap.nii.gz \
					shard/fieldmap.nii.gz

			fi

		fi

		kul_echo "mrtrix3new: $mrtrix3new"
		kul_echo "synb0: $synb0"
		kul_echo "shard: $shard"
		kul_echo "rev_only_topup: $rev_only_topup"
		kul_echo "n_pe: $n_pe"


		# Set the options for dwifslpreproc
		if [ $synb0 -eq 1 ] && [ $shard -eq 0 ]; then

			#  in case of synb0
			if [ $s_flag -eq 1 ];then
				dwifslprep_ses="ses-${ses}/"
			else
				dwifslprep_ses=""
			fi
			dwifslpreproc_option="-topup_files ${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}/${dwifslprep_ses}synb0/topup_"
		
		elif [ $synb0 -eq 0 ] && [ $shard -eq 0 ] && [ $n_pe -gt 1 ]; then
		
			#  in case a se_epi is given for tupop (a large number of b0)
			dwifslpreproc_option="-se_epi dwi/se_epi_for_topup.mif -align_seepi"
		
		elif [ $shard -eq 1 ] && [ $n_pe -gt 1 ]; then

			shard_fieldmap="-fieldmap shard/fieldmap.nii.gz"

		elif [ $shard -eq 1 ] && [ $n_pe -eq 1 ]; then	

			shard_fieldmap=""

		else

			# otherwise run the standard
			dwifslpreproc_option=""
		
		fi

		# Now run dwifslpreproc of shard-recon
		# note: maybe add -eddy_mask

		if [ $shard -eq 1 ]; then
			
			mkdir -p shard

			task_in1="mrinfo dwi/degibbs.mif -export_grad_mrtrix shard/grad.b -force"
			num_shells_tmp=($(mrinfo dwi/degibbs.mif -shell_bvalues))
			num_shells=${#num_shells_tmp[@]}

			if [ $num_shells -eq 2 ]; then
				lmax=""
				rlmax="-rlmax 2,2"
			elif [ $num_shells -eq 3 ]; then
				lmax=""
				rlmax="-rlmax 4,2,0"
			elif [ $num_shells -eq 6 ]; then
				lmax="-lmax 0,4,4,6,8,8"
				rlmax=""
			fi

			task_in2="dwimotioncorrect dwi/degibbs.mif shard/postmc-mssh.mif $shard_fieldmap \
				-mask dwi/dwi_orig_mask.nii.gz \
				$lmax $rlmax -mb ${mb} -sorder 1,0 -export_motion shard/motion.txt -export_weights shard/sliceweights.txt \
				-force -nocleanup -nthreads $ncpu"
			

			task_in3="mssh2amp shard/postmc-mssh.mif shard/grad.b dwi/geomcorr.mif -nonnegative"

			task_in="$task_in1;$task_in2;$task_in3"
			KUL_task_exec $verbose_level "kul_dwiprep part 3: shard-recon" "3_shard-recon"

			#motionstats motion.txt sliceweights.txt -sorder 1,0 -plot
		
		else
			
			# Make the directory for the output of eddy_qc
			mkdir -p eddy_qc/raw
			
			# for some old studies, the header info is not available
			# in that case one can set an environment variable 
			# e.g. export KUL_dwiprep_custom_dwifslpreproc="-rpe_none -pe_dir ap"
			# it is entirely up to the user to set the correct parameters
			if [ -z $KUL_dwiprep_custom_dwifslpreproc ]; then
				default_dwifslpreproc="-rpe_header"
			else
				default_dwifslpreproc=$KUL_dwiprep_custom_dwifslpreproc
				dwifslpreproc_option=""
			fi

			task_in="dwifslpreproc ${dwifslpreproc_option} $default_dwifslpreproc \
					-eddyqc_all eddy_qc/raw -eddy_options \"${full_eddy_options} \" -force -nthreads $ncpu -nocleanup \
					dwi/degibbs.mif dwi/geomcorr.mif"
			KUL_task_exec $verbose_level "kul_dwiprep part 3: motion and distortion correction" "3_dwifslpreproc"
		
		fi

		# create an intermediate mask of the dwi data
		kul_echo "creating intermediate mask of the dwi data..."
		dwi2mask_image_in="dwi/geomcorr.mif"
		dwi2mask_mask_out="dwi/dwi_intermediate_mask.nii.gz"
		dwi2mask_message="kul_dwiprep: make an intermediate mask"
		dwi2mask_logfile="3_mask"
		kul_dwi2mask

		temp_dir=$(ls -d *dwifslpreproc*)

		if [ ! command -v eddy_quad &> /dev/null ]; then

			kul_echo "Eddy_quad skipped (which was not found in the path)"

		elif [ $shard -eq 1 ]; then

			kul_echo "shard-recon qa is coming soon"

		else

			# Let's run eddy_quad
			rm -rf eddy_qc/quad

			kul_echo "   running eddy_quad..."
			task_in="eddy_quad $temp_dir/dwi_post_eddy --eddyIdx $temp_dir/eddy_indices.txt \
				--eddyParams $temp_dir/dwi_post_eddy.eddy_parameters --mask $temp_dir/eddy_mask.nii \
				--bvals $temp_dir/bvals --bvecs $temp_dir/bvecs --output-dir eddy_qc/quad" # --verbose
			KUL_task_exec $verbose_level "kul_dwiprep part 3: eddy_quad" "3_eddy_quad"

		fi

		# clean-up the above dwipreproc temporary directory
		# rm -rf $temp_dir

	else

		kul_echo "   part 2 of preprocessing has been done already... skipping to next step"

	fi


	# check if next 4 steps of dwi preprocessing are done
	if [ ! -f dwi_preproced.mif ]; then

		kul_echo "Start part 3 of preprocessing: dwibiascorrect, upsampling & creation of a final dwi_mask"

		# bias field correction
		kul_echo "    dwibiascorrect"
		task_in="dwibiascorrect ants dwi/geomcorr.mif dwi/biascorr.mif \
			-bias dwi/biasfield.mif -nthreads $ncpu -force -mask dwi/dwi_intermediate_mask.nii.gz"
		KUL_task_exec $verbose_level "kul_dwiprep part 4: dwibiascorrect" "4_dwibiascorrect"

		
		# upsample the images
		kul_echo "    upsampling resolution..."
		task_in="mrgrid -nthreads $ncpu -force -axis 1 5,5 dwi/biascorr.mif crop - | mrgrid -axis 1 5,5 -force - pad - | mrgrid -voxel 1.3 -force - regrid dwi/upsampled.mif"
		KUL_task_exec $verbose_level "kul_dwiprep part 5: upsampling resolution" "5_upsample"
		rm dwi/biascorr.mif
	

		# copy to main directory for subsequent processing
		kul_echo "    saving..."
		task_in="mrconvert dwi/upsampled.mif dwi_preproced.mif -set_property comments \"Preprocessed dMRI data.\" -nthreads $ncpu -force"
		KUL_task_exec $verbose_level "kul_dwiprep part 5: saving" "5_saving"
		rm dwi/upsampled.mif


		# create a final mask of the dwi data
		kul_echo "    creating mask of the dwi data..."
		dwi2mask_image_in="dwi_preproced.mif"
		dwi2mask_mask_out="dwi_mask.nii.gz"
		dwi2mask_message="kul_dwiprep part 6: make a final mask"
		dwi2mask_logfile="6_mask"
		kul_dwi2mask


		# create mean b0 of the dwi data
		kul_echo "    creating mean b0 of the dwi data ..."
		dwiextract -quiet dwi_preproced.mif -bzero - | mrmath -axis 3 - mean dwi_b0.nii.gz -force

	else

		kul_echo "part 3 of preprocessing has been done already... skipping to next step"

	fi



	# STEP 3 - RESPONSE ---------------------------------------------
	mkdir -p response
	# response function estimation (we compute following algorithms: tournier, tax and dhollander)

	kul_echo "Using dwipreproc_options: $dwipreproc_options"

	if [[ $dwipreproc_options == *"dhollander"* ]]; then
		if [ ! -f response/dhollander_wm_response.txt ]; then
			kul_echo "Calculating dhollander dwi2response..."
			task_in="dwi2response dhollander dwi_preproced.mif response/dhollander_wm_response.txt -mask dwi_mask.nii.gz \
			response/dhollander_gm_response.txt response/dhollander_csf_response.txt -nthreads $ncpu -force"
			KUL_task_exec $verbose_level "kul_dwiprep part 6: estimate response dhollander" "6_response_dhollander"
		else
			kul_echo " dwi2response dhollander already done, skipping..."
		fi

		if [ ! -f response/dhollander_wmfod.mif ]; then
			kul_echo "Calculating dhollander dwi2fod & normalising it..."

			task_in1="dwi2fod msmt_csd dwi_preproced.mif response/dhollander_wm_response.txt response/dhollander_wmfod.mif \
			response/dhollander_gm_response.txt response/dhollander_gm.mif \
			response/dhollander_csf_response.txt response/dhollander_csf.mif -mask dwi_mask.nii.gz -force -nthreads $ncpu"

			task_in2="dwi2fod msmt_csd dwi_preproced.mif response/dhollander_wm_response.txt response/dhollander_wmfod_noGM.mif \
			response/dhollander_csf_response.txt response/dhollander_csf_noGM.mif -mask dwi_mask.nii.gz -force -nthreads $ncpu"
			task_in="$task_in1;$task_in2"
			KUL_task_exec $verbose_level "kul_dwiprep part 7: dwi2fod dhollander" "6_dwi2fod_dhollander"
		else
			kul_echo " dwi2fod dhollander already done, skipping..."
		fi
	fi

	if [[ $dwipreproc_options == *"tax"* ]]; then
		if [ ! -f response/tax_response.txt ]; then
			kul_echo "Calculating tax dwi2response..."
			task_in1="dwi2response tax dwi_preproced.mif response/tax_response.txt -nthreads $ncpu -force"
			KUL_task_exec $verbose_level "kul_dwiprep part 6: estimate response tax" "6_response_tax"
		else
			kul_echo " dwi2response tax already done, skipping..."
		fi

		if [ ! -f response/tax_wmfod.mif ]; then
			kul_echo "Calculating tax dwi2fod..."
			task_in="dwi2fod csd dwi_preproced.mif response/tax_response.txt response/tax_wmfod.mif  \
			-mask dwi_mask.nii.gz -force -nthreads $ncpu -mask dwi_mask.nii.gz"
			KUL_task_exec $verbose_level "kul_dwiprep part 7: dwi2fod tax" "6_dwi2fod_tax"
		else
			kul_echo " dwi2fod tax already done, skipping..."
		fi
	fi

	if [[ $dwipreproc_options == *"tournier"* ]]; then
		if [ ! -f response/tournier_response.txt ]; then
			kul_echo "Calculating tournier dwi2response..."
			task_in="dwi2response tournier dwi_preproced.mif response/tournier_response.txt -nthreads $ncpu -force"
			KUL_task_exec $verbose_level "kul_dwiprep part 6: estimate response tournier" "6_response_tournier"
		else
			kul_echo " dwi2response already done, skipping..."
		fi

		if [ ! -f response/tournier_wmfod.mif ]; then
			kul_echo "Calculating tournier dwi2fod..."
			task_in="dwi2fod csd dwi_preproced.mif response/tournier_response.txt response/tournier_wmfod.mif  \
			-mask dwi_mask.nii.gz -force -nthreads $ncpu"
			KUL_task_exec $verbose_level "kul_dwiprep part 7: dwi2fod tax" "6_dwi2fod_tournier"
		else
			kul_echo " dwi2fod already done, skipping..."
		fi
	fi

	# STEP 4 - DO QA ---------------------------------------------
	# Make an FA/dec image

	mkdir -p qa

	if [ ! -f qa/noiselevel.nii.gz ]; then

		kul_echo "Calculating final FA/ADC/dec..."
		task_in="dwi2tensor dwi_preproced.mif dwi_dt.mif -force -mask dwi_mask.nii.gz; \
			tensor2metric dwi_dt.mif -fa qa/fa.nii.gz -mask dwi_mask.nii.gz -force; \
			tensor2metric dwi_dt.mif -adc qa/adc.nii.gz -mask dwi_mask.nii.gz -force"
		KUL_task_exec $verbose_level "kul_dwiprep part 8: final FA/ADC/DEC" "8_qa"

		if [[ $dwipreproc_options == *"tournier"* ]]; then
			task_in="fod2dec response/tournier_wmfod.mif qa/tournier_dec.mif -force -mask dwi_mask.nii.gz"
			KUL_task_exec $verbose_level "kul_dwiprep part 8: final FA/ADC/DEC" "8_qa"
		fi

		if [[ $dwipreproc_options == *"tax"* ]]; then
			task_in="fod2dec response/tax_wmfod.mif qa/tax_dec.mif -force -mask dwi_mask.nii.gz"
			KUL_task_exec $verbose_level "kul_dwiprep part 8: final FA/ADC/DEC" "8_qa"
		fi

		if [[ $dwipreproc_options == *"dhollander"* ]]; then
			task_in="fod2dec response/dhollander_wmfod.mif qa/dhollander_dec.mif -force -mask dwi_mask.nii.gz"
			KUL_task_exec $verbose_level "kul_dwiprep part 8: final FA/ADC/DEC" "8_qa"
		fi

		mrconvert dwi/biasfield.mif qa/biasfield.nii.gz
		mrconvert dwi/noiselevel.mif qa/noiselevel.nii.gz

	fi

	# make some QA figures
	if [ ! -f qa/fa.png ]; then
		mode_plane="-mode 2"
		overlay=""
		capture_dir="qa"
		mrview_exit="-exit"
		
		underlays="raw/b0s_pe*.mif"
		for f in $underlays; do
			underlay=$f
			if [ -f $underlay ]; then
				capture_prefix=$(basename -s .mif $f)
				KUL_mrview_figure.sh -p ${participant} -u $underlay -d "$(pwd)/$capture_dir" -f $capture_prefix -t 2 
				#kul_mrview_figure
			fi
		done
		
		underlay="qa/fa_orig.nii.gz"
		capture_prefix="fa_orig_"
		KUL_mrview_figure.sh -p ${participant} -u $underlay -d "$(pwd)/$capture_dir" -f $capture_prefix -t 2
		#kul_mrview_figure

		underlay="qa/fa.nii.gz"
		capture_prefix="fa_"
		KUL_mrview_figure.sh -p ${participant} -u $underlay -d "$(pwd)/$capture_dir" -f $capture_prefix -t 2
		#kul_mrview_figure

	fi

	# We finished processing current session
	# write a "done" log file for this session
	if [ $total_errorcount -eq 0 ]; then
		touch log/${current_session}_done.log
	else
		touch log/${current_session}_done_with_errors.log
	fi

	# STEP 5 - CLEANUP - here we clean up (large) temporary files
	#rm -fr dwi/degibbs.mif
	#rm -rf dwi/geomcorr.mif
	#rm -rf raw
	#rm -rf $temp_dir

	kul_echo "Finished processing session $bids_participant"


	# ---- END of the BIG loop over sessions
done

# write a file to indicate that dwiprep runned succesfully (or not)
if [ $total_errorcount -eq 0 ]; then
	touch ../dwiprep_is_done.log
	echo "Finished"
else
	touch ../dwiprep_is_done_with_errors.log
	echo "Finished with errors. Look at ${log}"
fi
