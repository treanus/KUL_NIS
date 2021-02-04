#!/bin/bash
# Bash shell script to prepare fMRI/DTI results for Brainlab Elements Server
#
# Requires Mrtrix3, Karawun
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/12/2020
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` prepares fMRI/DTI data for upload to Brainlab Elements Server.

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -s -v 

Required arguments:

	 -p:  participant name


Optional arguments:

	 -s:  simple colors (red for fMRI, blue for dti)
	 -v:  show output from mrtrix commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
simple=0 # default if option -s is not given
silent=1 # default if option -v is not given

# Set required options
#p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:sv" OPT; do

		case $OPT in
		s) #simple
			simple=1
		;;
		p) #participant
			participant=$OPTARG
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
#if [ $p_flag -eq 0 ] ; then
#	echo
#	echo "Option -p is required: give the BIDS name of the participant." >&2
#	echo
#	exit 2
#fi


# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then

	export MRTRIX_QUIET=1

fi

# --- MAIN ---
printf "\n\n\n"

# Check if all necessary data exists
if [ -f "DICOM_classic/DICOM/IM_0002" ]; then
    exist_dicom=1
else
    exist_dicom=0
    echo ""
    echo "Error: No classic dicom found (DICOM_classic/DICOM/IM_0002)"
    exit
fi
if test -n "$(find 16bit_Elements -maxdepth 1 -name '*rBO.img' -print -quit)"; then
    exist_B0=1
else
    exist_B0=0
    echo ""
    echo "Error: No B0.img in folder 16bit_Elements"
    exit
fi
if test -n "$(find 16bit_Anat -maxdepth 1 -name '*anat.img' -print -quit)"; then
    exist_T1=1
else
    exist_T1=0
    echo ""
    echo "Error: No *anat.img found in folder 16bit_Anat)"
    exit
fi

# Prepare the nii files for the tracts
mkdir -p for_elements_$participant/dti

# Convert the B0
echo "Now converting the B0"
mrconvert 16bit_Elements/*rBO.img for_elements_$participant/B0.nii

# convert each tract
cd 16bit_Elements
files=($(find -name "*_r_*.img" -type f -printf '%P\n'))

color=3
for t in "${files[@]}"; do 
    #echo $t
    printf "\n"
    d1=${t#*_r_}
    def=${d1%.*}
    #echo $def
    PS3="File ${t} corresponds to tract: "

    select opt in Other Ignore CST_L CST_R SLF_L SLF_R IFOF_L IFOF_R \
            ILF_L ILF_R SMAPT_L SMAPT_R FAT_L FAT_R CINGc_L \
            CINGc_R CINGt_L CINGt_R UNC_L UNC_R $def; do 
        #echo "you have selected $REPLY"
        #echo "this is $opt"
        if [ $simple -eq 0 ]; then
            color=$REPLY
        fi
        if [ $REPLY -eq 1 ]; then
            read -p "Describe the tract: " tract
        elif [ $REPLY -eq 2 ]; then
            break
        else
            tract=$opt
        fi    
        mrcalc $t 1 -gt $color -mult ../for_elements_$participant/dti/${tract}.nii
        break
    done
done

cd ..
#fi

# Prepare the nii for the fMRI act
mkdir -p for_elements_$participant/fmri

# Convert the T1w anat
printf "\n\n\n"
echo "Now converting the T1"
mrconvert 16bit_Anat/*anat.img for_elements_$participant/T1.nii

cd 16bit_Elements
files=( $(find -name "*_act_*.img" -type f -printf '%P\n') )

color=1
for fmri in "${files[@]}"; do 
    printf "\n"
    #echo $fmri
    d1=${fmri#*_act_}
    def=${d1%.*}
    #echo $def
    PS3="File ${fmri} corresponds to: "

    select opt in Other Ignore HAND LIP FOOT LANGUAGE $def; do
        if [ $REPLY -eq 1 ]; then
            read -p "Describe the task: " task
        elif [ $REPLY -eq 2 ]; then
            break
        else
            task=$opt
        fi
        if [ $simple -eq 0 ]; then
            color=$REPLY
        fi
        #echo $task
        #echo $color
        mrcalc $fmri 2 -gt $color -mult ../for_elements_$participant/fmri/${task}.nii
        break
    done
done
cd ..

# Call Karawun
printf "\n\n\n"
echo "Now preparing for Karawun DICOM conversion"

echo "Converting tracts now... wait..."
importTractography -d DICOM_classic/DICOM/IM_0002 -o for_elements_$participant/Upload_dti \
    -n for_elements_$participant/B0.nii -l for_elements_$participant/dti/*.nii

echo "Converting fMRI now... wait..."
importTractography -d DICOM_classic/DICOM/IM_0002 -o for_elements_$participant/Upload_fmri \
    -n for_elements_$participant/T1.nii -l for_elements_$participant/fmri/*.nii

echo "Finished"