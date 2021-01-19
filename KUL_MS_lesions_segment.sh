#!/bin/bash -e 
# Sarah Cappelle & Stefan Sunaert
# 19/01/2021
# This script is the first part of Sarah's Study1
# This script computes a MS lesion map using freesurfer samseg
# 

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` computes a lesion mask in MS subjects and the median of the MTR, T1T2 and T1FLAIR ratio's.

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -v 

Required arguments:

	 -p:  participant name


Optional arguments:

     -s:  session of the participant
     -a:  automatic mode (just work on all images in the BIDS folder)
     -n:  number of cpu to use
     -v:  show output from commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
auto=0 # default if option -s is not given
silent=1 # default if option -v is not given
outputdir="Study_MSlesions"
ncpu=6

# Set required options
#p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:n:av" OPT; do

		case $OPT in
		a) #automatic mode
			auto=1
		;;
		p) #participant
			participant=$OPTARG
		;;
        s) #session
			session=$OPTARG
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
#if [ $p_flag -eq 0 ] ; then
#	echo
#	echo "Option -p is required: give the BIDS name of the participant." >&2
#	echo
#	exit 2
#fi

ants_verbose=1
fs_silent=""
# verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
    ants_verbose=0
    fs_silent=" > /dev/null 2>&1"
fi

# --- FUNCTIONS ---
# --- MAIN ---
printf "\n\n\n"

# here we give the data
if [ $auto -eq 0 ]; then
    datadir="$cwd/BIDS/sub-${participant}/ses-$session/anat"
    T1w_all=("$datadir/sub-${participant}_ses-${session}_T1w.nii.gz")
    FLAIR=("$datadir/sub-${participant}_ses-${session}_FLAIR.nii.gz")
else
    T1w_all=($(find BIDS -type f -name "*T1w.nii.gz" | sort ))
fi

mkdir -p $outputdir/compute
touch $outputdir/mijn_resultaten.csv

echo "participant_and_session, \
    NAWM_mtr, NAWM_t1t2, NAWM_t1flair, \
    NAGM_mtr, NAGM_t1t2, NAGM_t1flair, \
    Thal_rh_mtr, Thal_rh_t1t2, Thal_rh_t1flair, \
    Thal_lh_mtr, Thal_lh_t1t2, Thal_lh_t1flair, \
    CSF_lateral_mtr, CSF_lateral_t1t2, CSF_lateral_t1flair, \
    CSF_3rd_mtr, CSF_3rd_t1t2, CSF_3rd_t1flair, \
    CSF_4th_mtr, CSF_4th_t1t2, CSF_4th_t1flair, \
    MSlesion_mtr, MSlesion_t1t2, MSlesion_t1flair \
    Volume_NAWM, Volume_NAGM, Volume_CSF, Volume_MSlesions" > $outputdir/mijn_resultaten.csv

flair=0
for test_T1w in ${T1w_all[@]}; do

    base0=${test_T1w##*/};base=${base0%_T1w*}
    check_done="$outputdir/compute/${base}.done"

    if [ ! -f $check_done ];then

        # Test whether FLAIR also exist
        test_FLAIR="${test_T1w%_T1w*}_FLAIR.nii.gz"
        if [ -f $test_FLAIR ];then
            #echo "The FLAIR exists"
            flair=1
            T1w=$test_T1w
            FLAIR=$test_FLAIR
        fi

        if [ $flair -eq 1 ];then 

            participant_and_session=$base

            T1w_iso="$cwd/$outputdir/compute/${participant_and_session}_T1w_iso.nii.gz"
            FLAIR_reg2T1w="$cwd/$outputdir/compute/${participant_and_session}_FLAIR_reg2T1w.nii.gz"

            # make the T1 iso 1mm 
            echo " starting on $T1w: make an 1x1x1mm isotropic"
            mrgrid $T1w regrid -voxel 1 $T1w_iso -force

            # coregister the flair to the T1w
            echo " coregistering FLAIR to T1w_iso"
            my_cmd="mri_coreg --mov $FLAIR --ref $T1w_iso \
            --reg $cwd/$outputdir/compute/${participant_and_session}_flair2T1.lta $fs_silent"
            eval $my_cmd
            my_cmd="mri_vol2vol --mov $FLAIR --reg $cwd/$outputdir/compute/${participant_and_session}_flair2T1.lta \
                --o $FLAIR_reg2T1w --targ $T1w_iso $fs_silent"
            eval $my_cmd

            # run samseg
            echo " running samseg (takes about 20 minutes)"
            my_cmd="run_samseg --input $T1w_iso $FLAIR_reg2T1w --pallidum-separate \
            --lesion --lesion-mask-pattern 0 1 --output $cwd/$outputdir/compute/samsegOutput_$participant_and_session \
            --threads $ncpu $fs_silent"
            eval $my_cmd

            # compute masks
            SamSeg="$cwd/$outputdir/compute/samsegOutput_${participant_and_session}/seg.mgz"
            MTR="T1T2FLAIRMTR_ratio/${participant_and_session}_MTC_ratio.nii.gz"
            T1T2="T1T2FLAIRMTR_ratio/${participant_and_session}_T1T2w_ratio.nii.gz"
            T1FLAIR="T1T2FLAIRMTR_ratio/${participant_and_session}_T1FLAIR_ratio.nii.gz"

            NAWM_lh="$cwd/$outputdir/${participant_and_session}_NAWM_lh.nii.gz"
            NAWM_rh="$cwd/$outputdir/${participant_and_session}_NAWM_rh.nii.gz"
            NAWM="$cwd/$outputdir/${participant_and_session}_NAWM.nii.gz"
            NAGM_lh="$cwd/$outputdir/${participant_and_session}_NAGM_lh.nii.gz"
            NAGM_rh="$cwd/$outputdir/${participant_and_session}_NAGM_rh.nii.gz"
            NAGM="$cwd/$outputdir/${participant_and_session}_NAGM.nii.gz"
            Thal_lh="$cwd/$outputdir/${participant_and_session}_Thal_lh.nii.gz"
            Thal_rh="$cwd/$outputdir/${participant_and_session}_Thal_rh.nii.gz"
            CSF_lateral_lh="$cwd/$outputdir/${participant_and_session}_CSF_lateral_lh.nii.gz"
            CSF_lateral_rh="$cwd/$outputdir/${participant_and_session}_CSF_lateral_rh.nii.gz"
            CSF_lateral="$cwd/$outputdir/${participant_and_session}_CSF_lateral.nii.gz"
            CSF_3rd="$cwd/$outputdir/${participant_and_session}_CSF_3rd.nii.gz"
            CSF_4th="$cwd/$outputdir/${participant_and_session}_CSF_4th.nii.gz"
            CSF="$cwd/$outputdir/${participant_and_session}_CSF.nii.gz"
            MSlesion="$cwd/$outputdir/${participant_and_session}_MSLesion.nii.gz"

            mrcalc $SamSeg 2 -eq $NAWM_lh -force
            mrcalc $SamSeg 41 -eq $NAWM_rh -force
            mrcalc $NAWM_lh $NAWM_rh -add $NAWM -force
            mrcalc $SamSeg 3 -eq $NAGM_lh -force
            mrcalc $SamSeg 42 -eq $NAGM_rh -force
            mrcalc $NAGM_lh $NAGM_rh -add $NAGM -force
            mrcalc $SamSeg 10 -eq $Thal_lh -force
            mrcalc $SamSeg 49 -eq $Thal_rh -force
            mrcalc $SamSeg 4 -eq $CSF_lateral_lh -force
            mrcalc $SamSeg 43 -eq $CSF_lateral_rh -force
            mrcalc $CSF_lateral_lh $CSF_lateral_rh -add $CSF_lateral -force
            mrcalc $SamSeg 14 -eq $CSF_3rd -force
            mrcalc $SamSeg 15 -eq $CSF_4th -force
            mrcalc $CSF_lateral $CSF_3rd $CSF_4th -add $CSF -force
            mrcalc $SamSeg 99 -eq $MSlesion -force

            NAWM_mtr=$(mrstats -mask $NAWM -output median $MTR)
            NAWM_t1t2=$(mrstats -mask $NAWM -output median $T1T2)
            NAWM_t1flair=$(mrstats -mask $NAWM -output median $T1FLAIR)
            NAGM_mtr=$(mrstats -mask $NAGM -output median $MTR)
            NAGM_t1t2=$(mrstats -mask $NAGM -output median $T1T2)
            NAGM_t1flair=$(mrstats -mask $NAGM -output median $T1FLAIR)
            Thal_rh_mtr=$(mrstats -mask $Thal_rh -output median $MTR)
            Thal_rh_t1t2=$(mrstats -mask $Thal_rh -output median $T1T2)
            Thal_rh_t1flair=$(mrstats -mask $Thal_rh -output median $T1FLAIR)
            Thal_lh_mtr=$(mrstats -mask $Thal_lh -output median $MTR)
            Thal_lh_t1t2=$(mrstats -mask $Thal_lh -output median $T1T2)
            Thal_lh_t1flair=$(mrstats -mask $Thal_lh -output median $T1FLAIR)
            CSF_lateral_mtr=$(mrstats -mask $CSF_lateral -output median $MTR)
            CSF_lateral_t1t2=$(mrstats -mask $CSF_lateral -output median $T1T2)
            CSF_lateral_t1flair=$(mrstats -mask $CSF_lateral -output median $T1FLAIR)
            CSF_3rd_mtr=$(mrstats -mask $CSF_3rd -output median $MTR)
            CSF_3rd_t1t2=$(mrstats -mask $CSF_3rd -output median $T1T2)
            CSF_3rd_t1flair=$(mrstats -mask $CSF_3rd -output median $T1FLAIR)
            CSF_4th_mtr=$(mrstats -mask $CSF_4th -output median $MTR)
            CSF_4th_t1t2=$(mrstats -mask $CSF_4th -output median $T1T2)
            CSF_4th_t1flair=$(mrstats -mask $CSF_4th -output median $T1FLAIR)
            
            MSlesion_mtr=$(mrstats -mask $MSlesion -output median $MTR)
            MSlesion_t1t2=$(mrstats -mask $MSlesion -output median $T1T2)
            MSlesion_t1flair=$(mrstats -mask $MSlesion -output median $T1FLAIR)

            Volume_NAWM=$(mrstats -ignorezero -output count $NAWM)
            Volume_NAGM=$(mrstats -ignorezero -output count $NAGM)
            Volume_CSF=$(mrstats -ignorezero -output count $CSF)
            Volume_MSlesions=$(mrstats -ignorezero -output count $MSlesion)

            echo "$participant_and_session, \
                $NAWM_mtr, $NAWM_t1t2, $NAWM_t1flair, \
                $NAGM_mtr, $NAGM_t1t2, $NAGM_t1flair, \
                $Thal_rh_mtr, $Thal_rh_t1t2, $Thal_rh_t1flair, \
                $Thal_lh_mtr, $Thal_lh_t1t2, $Thal_lh_t1flair, \
                $CSF_lateral_mtr, $CSF_lateral_t1t2, $CSF_lateral_t1flair, \
                $CSF_3rd_mtr, $CSF_3rd_t1t2, $CSF_3rd_t1flair, \
                $CSF_4th_mtr, $CSF_4th_t1t2, $CSF_4th_t1flair, \
                $MSlesion_mtr, $MSlesion_t1t2, $MSlesion_t1flair \
                $Volume_NAWM, $Volume_NAGM, $Volume_CSF, $Volume_MSlesions" >> $outputdir/mijn_resultaten.csv
        
            touch $check_done

        else

            echo " $participant_and_session does not have a flair - nothing to do here"
        
        fi
    
    else
        echo " $participant_and_session  already done".
    fi

done


