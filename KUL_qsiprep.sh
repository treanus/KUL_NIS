#!/bin/bash -e
# Bash shell script to run qsiprep or mrtrix_connectome
#
# Requires docker
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 19/03/2021
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` runs qsiprep tuned for KUL/UZLeuven data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe 

Required arguments:

     -p:  participant name

Optional arguments:
     
     -s:  session
	 -w:  workflow (1=mrtrix_connectome, 2=qsiprep; default:1) 
     -t:  susceptibility correction model (1=synb0, 2=topup j, 3=topup j-; 4=none default:1)
     -m:  hmc_model when using qsiprep workflow (1=none,2=eddy,3=3dSHORE; default:2)
     -g:  use gpu (does not work an MacOs)
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
hmc=2
gpu=0
sdc=1
wfl=1
session=""

# Set required options
p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "w:p:s:n:m:t:gv" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
		s) #session
			session=$OPTARG
		;;
        w) #workflow
			wfl=$OPTARG
		;;
		m) #hmc
			hmc=$OPTARG
		;;
		t) #sdc
			sdc=$OPTARG
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
		g) #gpu
			gpu=1
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
	echo "Option -p is required: give the BIDS name of the participant." >&2
	echo
	exit 2
fi

if [ $hmc -eq 1 ]; then
    hmc_type="none"
elif [ $hmc -eq 2 ]; then
    hmc_type="eddy"
elif [ $hmc -eq 3 ]; then
    hmc_type="3dSHORE"
else
    echo "Wrong hmc type; exitting"
    exit
fi

if [ $sdc -eq 1 ]; then
    sdc_type="synb0"
	dir_epi="PA"
elif [ $sdc -eq 2 ]; then
    sdc_type="topup j"
	topup_pe_dir="j"
	dir_epi="PA"
elif [ $sdc -eq 3 ]; then
    sdc_type="topup j-"
	topup_pe_dir="j-"
	dir_epi="AP"	
elif [ $sdc -eq 4 ]; then
    sdc_type="none"
else
    echo "Wrong sdc type; exitting"
    exit
fi

# set the version of qsiprep
version="0.13.1"
version="latest"

# ----  MAIN  ------------------------

#echo $participant
if [ -z $session ];then
	sessuf=""
else
	sessuf="/ses-"
fi
bids_subj=BIDS/sub-${participant}${sessuf}${session}
echo $bids_subj

qsi_data="${cwd}/BIDS"
qsi_scratch="${cwd}/qsiprep_work_${participant}"


if [ $gpu -eq 1 ]; then
	gpu_cmd1="--gpus all"
	gpu_cmd2="--eddy-config /data/derivatives/eddy_params.json"
	eddy_config_settings="{
  	\"flm\": \"linear\",
	\"slm\": \"linear\",
	\"fep\": false,
	\"interp\": \"spline\",
	\"nvoxhp\": 1000,
	\"fudge_factor\": 10,
	\"dont_sep_offs_move\": false,
	\"dont_peas\": false,
	\"niter\": 5,
	\"method\": \"jac\",
	\"repol\": true,
	\"num_threads\": 1,
	\"is_shelled\": true,
	\"use_cuda\": true,
	\"cnr_maps\": true,
	\"residuals\": false,
	\"output_type\": \"NIFTI_GZ\",
	\"args\": \"\"
	}"
	echo $eddy_config_settings > $qsi_data/derivatives/eddy_params.json
else
	gpu_cmd1=""
	gpu_cmd2="--eddy-config /data/derivatives/eddy_params.json"
	eddy_config_settings="{
  	\"flm\": \"linear\",
	\"slm\": \"linear\",
	\"fep\": false,
	\"interp\": \"spline\",
	\"nvoxhp\": 1000,
	\"fudge_factor\": 10,
	\"dont_sep_offs_move\": false,
	\"dont_peas\": false,
	\"niter\": 5,
	\"method\": \"jac\",
	\"repol\": true,
	\"num_threads\": $ncpu,
	\"is_shelled\": true,
	\"use_cuda\": false,
	\"cnr_maps\": true,
	\"residuals\": false,
	\"output_type\": \"NIFTI_GZ\",
	\"args\": \"\"
	}"
	echo $eddy_config_settings > $qsi_data/derivatives/eddy_params.json
fi

# run synb0
# prepare for Synb0-disco
if [ $sdc -eq 1 ]; then

	synb0_scratch="${cwd}/synb0_work_${bids_subj}"
	bids_dmri_found=($(find $cwd/$bids_subj/dwi -type f -name "*dwi.nii.gz")) 
	number_of_bids_fmri_found=${#bids_dmri_found[@]}

	test_file=$synb0_scratch/OUTPUTS/b0_u.nii.gz

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

		cp $Synb0_T1 $synb0_scratch/INPUTS/T1.nii.gz


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
			hansencb/synb0 \
			--notopup"
		echo $cmd
		eval $cmd

		synb0out=2

		if [ -z $session ];then
			sessuf2=""
		else
			sessuf2="_ses-${session}"
		fi

		if [ $synb0out -eq 1 ];then 
			# make a json
			json_file=${synb0_scratch}/sub-${participant}_fieldmap.json
			echo "{" > $json_file
			echo "\"Units\": \"Hz\"," >> $json_file
			echo "\"IntendedFor\": [" >> $json_file
			for i in `seq 0 $(($number_of_bids_fmri_found-1))`; do
				if [ $i -eq $(($number_of_bids_fmri_found-1)) ];then
					comma=""
				else
					comma=", "
				fi
				echo "\"${bids_dmri_found[i]#*$participant/}\"$comma" >> $json_file
			done
			echo "]" >> $json_file
			echo "}" >> $json_file


			# add these to the BIDS
			mkdir -p $qsi_data/sub-${participant}/fmap
			cp $json_file $qsi_data/sub-${participant}/fmap
			cp $synb0_scratch/OUTPUTS/topup_fieldcoef.nii.gz $qsi_data/sub-${participant}/fmap/sub-${participant}_fieldmap.nii.gz
		
		else

			# make a json

			sessuf2="ses-$session"
			json_file=${synb0_scratch}/sub-${participant}${sessuf2}_dir-${dir_epi}_epi.json
			echo "{" > $json_file
			echo "\"PhaseEncodingDirection\": \"j\"," >> $json_file
			echo "\"TotalReadoutTime\": 0.000," >> $json_file
			echo "\"IntendedFor\": [" >> $json_file
			for i in `seq 0 $(($number_of_bids_fmri_found-1))`; do
				if [ $i -eq $(($number_of_bids_fmri_found-1)) ];then
					comma=""
				else
					comma=", "
				fi
				echo "\"${bids_dmri_found[i]#*$participant/}\"$comma" >> $json_file
			done
			echo "]" >> $json_file
			echo "}" >> $json_file


			# add these to the BIDS
			mkdir -p ${cwd}/${bids_subj}/fmap
			cp $json_file ${cwd}/${bids_subj}/fmap
			cp $synb0_scratch/OUTPUTS/b0_u.nii.gz ${cwd}/${bids_subj}/fmap/sub-${participant}${sessuf2}_dir-${dir_epi}_epi.nii.gz

		fi
	
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


if [ $wfl -eq 1 ]; then
	
	outputdir="$cwd/MRtrix3_connectome"
	scratchdir="$cwd/MRtrix3_connectome_sub-${participant}"

	test_file="$cwd/MRtrix3_connectome/MRtrix3_connectome-preproc/sub-${participant}/dwi/sub-${participant}_desc-preproc_dwi.nii.gz"
	
	if [ ! -f $test_file ];then
		my_cmd="docker run -i --rm \
			-v $cwd/BIDS:/bids_dataset \
			-v $outputdir:/output \
			$gpu_cmd1 \
			bids/mrtrix3_connectome \
			/bids_dataset /output preproc --participant_label $participant \
			--output_verbosity 4 \
			--template_reg ants "
	else
		my_cmd="echo Already preprocessed"
	fi

elif [ $wfl -eq 2 ]; then

	outputdir="$cwd/qsiprep"

	my_cmd="docker run --rm -it \
		-v $FS_LICENSE:/opt/freesurfer/license.txt:ro \
		-v $qsi_data:/data:ro \
		-v $qsi_out:/out \
		-v $qsi_scratch:/scratch \
		$gpu_cmd1 \
		pennbbl/qsiprep:$version \
		/data /out participant \
		-w /scratch \
		--run-uuid $(id -u) \
		--acquisition_type main \
		--output-resolution 1.2 \
		--hmc_model $hmc_type \
		$gpu_cmd2 \
		--prefer-dedicated-fmaps \
		--participant_label $participant \
		--recon_spec mrtrix_multishell_msmt_noACT "
fi

echo $my_cmd

mkdir -p $outputdir

eval $my_cmd


