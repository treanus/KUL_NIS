#!/bin/bash -e
# Bash shell script to run synb0 from BIDS and store the output of topup in the BIDS derivatives
#
# Requires docker
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 06/09/2021
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` runs synb0 tuned for KUL/UZLeuven data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -s precovid -n 6

Required arguments:

     -p:  participant name

Optional arguments:
     
     -s:  session
	 -c:  cleanup the topup_fieldmap (remove signal in air anterior from eyes) in development
     -n:  number of cpu to use (default 15)
     -v:  show output from commands

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ncpu=15
session=""
sdc=1
cleanup=0
dir_epi="PA"

# Set required options
p_flag=0
s_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:n:cv" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
		s) #session
			session=$OPTARG
			s_flag=1
		;;
		c) #cleanup
			cleanup=1
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
        v) #verbose
			silent=0
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
	echo "Option -p (or -a) is required: give the BIDS name of the participant." >&2
	echo
	exit 2
fi

# ----  MAIN  ------------------------

# Either a session is given on the command line
# If not the session(s) need to be determined.
if [ $s_flag -eq 1 ]; then

    # session is given on the command line
    search_sessions=BIDS/sub-${participant}/ses-${session}

else

    # search if any sessions exist
    search_sessions=($(find BIDS/sub-${participant} -type d | grep dwi))

fi    
 
num_sessions=${#search_sessions[@]}
    
echo "  Number of BIDS sessions: $num_sessions"
echo "    notably: ${search_sessions[@]}"


# ---- BIG LOOP for processing each session
for i in `seq 0 $(($num_sessions-1))`; do

	# set up 
	long_bids_subj=${search_sessions[$i]}
	bids_subj=${long_bids_subj%dwi}
	#echo $bids_subj

	if [[ $bids_subj == *"ses-"* ]];then
		ses=${bids_subj#*ses-}
		ses=${ses%*/}
		sessuf1="/ses-${ses}"
		sessuf2="_ses-${ses}"
	else
		ses=""
		sessuf1=""
		sessuf2=""
	fi
	#echo $ses
	#echo $sessuf1
	#echo $sessuf2
	

	# run synb0
	# prepare for Synb0-disco
	if [ $sdc -eq 1 ]; then

		synb0_scratch="${cwd}/synb0_work_${bids_subj}"
		bids_dmri_found=($(find $cwd/$bids_subj/dwi -type f -name "*dwi.nii.gz")) 
		number_of_bids_dmri_found=${#bids_dmri_found[@]}

		test_file=$synb0_scratch/OUTPUTS/topup_fieldcoef.nii.gz

		if [ ! -f $test_file ];then

			# clean
			rm -fr $synb0_scratch
			
			# setup
			mkdir -p $synb0_scratch/INPUTS
			mkdir -p $synb0_scratch/OUTPUTS


			# find T1
			bids_T1_found=($(find $cwd/$bids_subj/anat -type f -name "*T1w.nii.gz" ! -name "*gadolinium*")) 
			number_of_bids_T1_found=${#bids_T1_found[@]}
			if [ $number_of_bids_T1_found -gt 1 ]; then
				kul_e2cl "   more than 1 T1 dataset, using first only for Synb0-disco" ${preproc}/${log}
			fi
			#echo $bids_T1_found
			Synb0_T1=${bids_T1_found[0]}
			echo "The used T1 for synb0-disco is $Synb0_T1"

			cp $Synb0_T1 $synb0_scratch/INPUTS/T1_full.nii.gz


			# extract the B0
			# find dMRI
			if [ $number_of_bids_dmri_found -gt 1 ]; then
				kul_e2cl "   more than 1 dMRI b0 dataset, using first only for Synb0-disco" ${preproc}/${log}
			fi
			#echo $bids_dmri_found
			Synb0_dmri=${bids_dmri_found[0]}
			echo "The used dMRI for synb0-disco is $Synb0_dmri"
			dwi_base=${Synb0_dmri%%.*}

			mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
				-json_import ${dwi_base}.json $synb0_scratch/dwi_p1.mif -strides 1:3 -force -clear_property comments -nthreads $ncpu

			dwiextract -quiet -bzero $synb0_scratch/dwi_p1.mif $synb0_scratch/dwi_p1_b0s.mif -force
			mrconvert $synb0_scratch/dwi_p1_b0s.mif -coord 3 0 $synb0_scratch/INPUTS/b0.nii.gz -strides -1,+2,+3,+4 -export_pe_table $synb0_scratch/topup_datain.txt


			# adjust the FOV of the T1 to match the b0
			mrgrid $synb0_scratch/INPUTS/b0.nii.gz regrid $synb0_scratch/INPUTS/b0_as_T1.nii.gz \
				-template $synb0_scratch/INPUTS/T1_full.nii.gz
			mrgrid -mask $synb0_scratch/INPUTS/b0_as_T1.nii.gz $synb0_scratch/INPUTS/T1_full.nii.gz crop \
				$synb0_scratch/INPUTS/T1_crop.nii.gz
			mrgrid $synb0_scratch/INPUTS/T1_crop.nii.gz crop -axis 1 10,10 $synb0_scratch/INPUTS/T1.nii.gz
			rm $synb0_scratch/INPUTS/T1_crop.nii.gz
			rm $synb0_scratch/INPUTS/T1_full.nii.gz
			rm $synb0_scratch/INPUTS/b0_as_T1.nii.gz
		

			# read topup_datain.txt and add line
			topup_data=($(cat $synb0_scratch/topup_datain.txt))
			echo "${topup_data[0]} ${topup_data[1]} ${topup_data[2]} ${topup_data[3]}" > $synb0_scratch/INPUTS/acqparams.txt
			echo "${topup_data[0]} ${topup_data[1]} ${topup_data[2]} 0.000" >> $synb0_scratch/INPUTS/acqparams.txt


			# run synb0
			cmd="docker run -u $(id -u) --rm \
				-v $synb0_scratch/INPUTS/:/INPUTS/ \
				-v $synb0_scratch/OUTPUTS:/OUTPUTS/ \
				-v $FS_LICENSE:/extra/freesurfer/license.txt \
				--user $(id -u):$(id -g) \
				hansencb/synb0" 

			echo "  we run synb0 using command: $cmd"
			eval $cmd


			# make a json
			json_file=${synb0_scratch}/sub-${participant}${sessuf2}_dir-${dir_epi}_epi.json
			echo "{" > $json_file
			echo "\"PhaseEncodingDirection\": \"j\"," >> $json_file
			echo "\"TotalReadoutTime\": 0.000," >> $json_file
			# need to reprogram this, so that there is a synthetic synb0 in fmap for each dwi-dataset 
			echo "\"IntendedFor\": " >> $json_file
			for i in `seq 0 $(($number_of_bids_dmri_found-1))`; do
				if [ $i -eq $(($number_of_bids_dmri_found-1)) ];then
					comma=""
				else
					comma=", "
				fi
				echo "\"${bids_dmri_found[i]#*$participant/}\"$comma" >> $json_file
			done
			#echo "]" >> $json_file
			echo "}" >> $json_file


			# add these to the BIDS derivatives		
			mkdir -p ${cwd}/BIDS/derivatives/synb0/sub-${participant}${sessuf1}/topup
			if [ $cleanup -eq 1 ];then
				mrgrid $synb0_scratch/OUTPUTS/topup_fieldcoef.nii.gz crop -axis 1 0,5 - | \
					mrgrid - pad -axis 1 0,5 $synb0_scratch/OUTPUTS/topup_fieldcoef_clean.nii.gz
				fieldmap2copy=$synb0_scratch/OUTPUTS/topup_fieldcoef_clean.nii.gz
			else
				fieldmap2copy=$synb0_scratch/OUTPUTS/topup_fieldcoef.nii.gz
			fi
			cp $fieldmap2copy \
				${cwd}/BIDS/derivatives/synb0/sub-${participant}${sessuf1}/topup/topup_fieldcoef.nii.gz
				#${cwd}/BIDS/derivatives/synb0/sub-${participant}${sessuf1}/topup/sub-${participant}${sessuf2}_topup_fieldcoef.nii.gz
			cp $synb0_scratch/OUTPUTS/topup_movpar.txt \
				${cwd}/BIDS/derivatives/synb0/sub-${participant}${sessuf1}/topup/topup_movpar.txt
				#${cwd}/BIDS/derivatives/synb0/sub-${participant}${sessuf1}/topup/sub-${participant}${sessuf2}_topup_movpar.txt
		
		else

			echo "  $bids_subj already done"

		fi


	elif [ $sdc -eq 2 ]; then

		# TOPUP CASE

		topup_scratch="${cwd}/topup_work_${participant}"
		test_file=$qsi_data/sub-${participant}/fmap/sub-${participant}_dir-${dir_epi}_epi.nii.gz

		if [ ! -f $test_file ];then

			# clean
			rm -fr $topup_scratch
			
			# setup
			mkdir -p $topup_scratch

			# find the dMRI to be used for topup
			bids_dmri_found=($(find $cwd/$bids_subj/dwi -type f -name "*dwi.nii.gz")) 
			#echo ${bids_dmri_found[@]}
			number_of_bids_dmri_found=${#bids_dmri_found[@]}
			#echo $number_of_bids_dmri_found

			i=0
			intended_for=""
			for dmri in ${bids_dmri_found[@]}; do
				dwi_base=${dmri%%.*}
				pe=$(grep PhaseEncodingDirection\" $dwi_base.json | awk  'BEGIN{FS="\""}{print $4}')
				if [ $pe = $topup_pe_dir ]; then
					echo "The $topup_pe_dir direction is file $dwi_base"
					p=$i
				else
					intended_for[i]=${dmri[i]}
				fi
				i=$((i++))
				#echo "i= $i"
			done

			#echo "p= $p"

			dwi_base=${bids_dmri_found[$p]%%.*}
			echo $dwi_base
			mrconvert ${dwi_base}.nii.gz -fslgrad ${dwi_base}.bvec ${dwi_base}.bval \
					-json_import ${dwi_base}.json $topup_scratch/dwi_p1.mif -strides 1:3 -force -clear_property comments -nthreads $ncpu
			dwiextract -quiet -bzero $topup_scratch/dwi_p1.mif $topup_scratch/dwi_p1_b0s.mif -force
			mrconvert $topup_scratch/dwi_p1_b0s.mif -coord 3 0 $topup_scratch/b0.nii.gz -strides -1,+2,+3,+4 -export_pe_table $topup_scratch/topup_datain.txt
			topup_data=($(cat $topup_scratch/topup_datain.txt))
			trt=${topup_data[3]}

			# make a json
			json_file=${topup_scratch}/sub-${participant}_dir-${dir_epi}_epi.json
			number_of_intended_for=${#intended_for[@]}
			#echo $number_of_intended_for
			echo "{" > $json_file
			echo "\"PhaseEncodingDirection\": \"$topup_pe_dir\"," >> $json_file
			echo "\"TotalReadoutTime\": $trt," >> $json_file
			echo "\"IntendedFor\": [" >> $json_file
			for i in `seq 0 $(($number_of_intended_for-1))`; do
				#echo $i
				if [ $i -eq $(($number_of_intended_for-1)) ];then
					comma=""
				else
					comma=", "
				fi
				echo "\"${intended_for[i]#*$participant/}\"$comma" >> $json_file
			done
			echo "]" >> $json_file
			echo "}" >> $json_file


			# add these to the BIDS
			mkdir -p $qsi_data/sub-${participant}/fmap
			cp $json_file $qsi_data/sub-${participant}/fmap
			mrconvert $topup_scratch/dwi_p1_b0s.mif $qsi_data/sub-${participant}/fmap/sub-${participant}_dir-${dir_epi}_epi.nii.gz -force
		
		fi

	elif [ $sdc -eq 4 ]; then

		echo "Not doing any SDC"

	fi

done
