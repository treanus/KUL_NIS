#!/bin/bash -e
# Sarah Cappelle & Stefan Sunaert
# 19/01/2021
# This script is the first part of Sarah's Study1
# This script computes a MS lesion map using freesurfer samseg
# 
v="0.9"

kul_main_dir=`dirname "$0"`
script="$0"
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` computes statistics on the MSlesion images as well as on MTR, T1T2 and T1FLAIR ratio's.

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -v 

Required arguments:

	 -p:  participant name


Optional arguments:

     -s:  session of the participant
     -a:  automatic mode (just work on all images in the T1T2FLAIRMTR folder)
     -t:  type (1=nocalib, 2=lincalib, 3=nonlincalib, 4=nonlincalib2)
     -v:  show output from commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
auto=0 # default if option -s is not given
silent=1 # default if option -v is not given
outputdir="T1T2FLAIRMTR_ratio"
type_sel=1;type=""

# Set required options
#p_flag=0

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:s:n:t:av" OPT; do

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
        t) #session
			type_sel=$OPTARG
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


# verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
fi

# --- FUNCTIONS ---
function KUL_create_results_file {
    if [ ! -f $outdir/stats/$my_results_file ];then
        touch $outdir/stats/$my_results_file
        echo "participant_and_session, \
        NAWM_mtr, NAWM_t1t2, NAWM_t1flair, \
        NAGM_mtr, NAGM_t1t2, NAGM_t1flair, \
        Thal_rh_mtr, Thal_rh_t1t2, Thal_rh_t1flair, \
        Thal_lh_mtr, Thal_lh_t1t2, Thal_lh_t1flair, \
        CSF_lateral_mtr, CSF_lateral_t1t2, CSF_lateral_t1flair, \
        CSF_3rd_mtr, CSF_3rd_t1t2, CSF_3rd_t1flair, \
        CSF_4th_mtr, CSF_4th_t1t2, CSF_4th_t1flair, \
        MSlesion_mtr, MSlesion_t1t2, MSlesion_t1flair, \
        Volume_NAWM, Volume_NAGM, Volume_CSF, Volume_MSlesions" > $outdir/stats/$my_results_file
    fi
}

function KUL_compute_stats {
    # define the input images
    SamSeg="$cwd/$outdir/samseg/seg.mgz"
    MTR="$outdir/${participant_and_session}_ratio-MTC.nii.gz"
    T1T2="$outdir/${participant_and_session}_ratio-T1T2w_${type}.nii.gz"
    T1FLAIR="$outdir/${participant_and_session}_ratio-T1FLAIR_${type}.nii.gz"
    echo $MTR
    echo $T1T2

    # define the output images
    NAWM_lh="$cwd/$outdir/rois/${participant_and_session}_NAWM_lh.nii.gz"
    NAWM_rh="$cwd/$outdir/rois/${participant_and_session}_NAWM_rh.nii.gz"
    NAWM="$cwd/$outdir/rois/${participant_and_session}_NAWM.nii.gz"
    NAGM_lh="$cwd/$outdir/rois/${participant_and_session}_NAGM_lh.nii.gz"
    NAGM_rh="$cwd/$outdir/rois/${participant_and_session}_NAGM_rh.nii.gz"
    NAGM="$cwd/$outdir/rois/${participant_and_session}_NAGM.nii.gz"
    Thal_lh="$cwd/$outdir/rois/${participant_and_session}_Thal_lh.nii.gz"
    Thal_rh="$cwd/$outdir/rois/${participant_and_session}_Thal_rh.nii.gz"
    CSF_lateral_lh="$cwd/$outdir/rois/${participant_and_session}_CSF_lateral_lh.nii.gz"
    CSF_lateral_rh="$cwd/$outdir/rois/${participant_and_session}_CSF_lateral_rh.nii.gz"
    CSF_lateral="$cwd/$outdir/rois/${participant_and_session}_CSF_lateral.nii.gz"
    CSF_3rd="$cwd/$outdir/rois/${participant_and_session}_CSF_3rd.nii.gz"
    CSF_4th="$cwd/$outdir/rois/${participant_and_session}_CSF_4th.nii.gz"
    CSF="$cwd/$outdir/rois/${participant_and_session}_CSF.nii.gz"
    MSlesion="$cwd/$outdir/rois/${participant_and_session}_MSLesion.nii.gz"

    echo " making VOIs"
    # do the computation of the masks
    if [ ! -f $CSF ]; then
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
        mrcalc $CSF_lateral $CSF_3rd -add $CSF_4th -add $CSF -force
        #mrcalc $SamSeg 99 -eq $MSlesion -force
    fi

    echo " computing stats"
    # do the stats
    NAGM_mtr="NA"
    NAWM_mtr="NA"
    Thal_rh_mtr="NA"
    Thal_lh_mtr="NA"
    CSF_lateral_mtr="NA"
    CSF_3rd_mtr="NA"
    CSF_4th_mtr="NA"
    MSlesion_mtr="NA"
    if [ -f $MTR ]; then
        NAGM_mtr=$(mrstats -mask $NAGM -output median $MTR)
        NAWM_mtr=$(mrstats -mask $NAWM -output median $MTR)
        Thal_rh_mtr=$(mrstats -mask $Thal_rh -output median $MTR)
        Thal_lh_mtr=$(mrstats -mask $Thal_lh -output median $MTR)
        CSF_lateral_mtr=$(mrstats -mask $CSF_lateral -output median $MTR)
        CSF_3rd_mtr=$(mrstats -mask $CSF_3rd -output median $MTR)
        CSF_4th_mtr=$(mrstats -mask $CSF_4th -output median $MTR)
        MSlesion_mtr=$(mrstats -mask $MSlesion -output median $MTR)
    fi

    NAGM_t1t2="NA"
    NAWM_t1t2="NA"
    Thal_rh_t1t2="NA"
    Thal_lh_t1t2="NA"
    CSF_lateral_t1t2="NA"
    CSF_3rd_t1t2="NA"
    CSF_4th_t1t2="NA"
    MSlesion_t1t2="NA"
    if [ -f $T1T2 ]; then
        NAWM_t1t2=$(mrstats -mask $NAWM -output median $T1T2)
        NAGM_t1t2=$(mrstats -mask $NAGM -output median $T1T2)
        Thal_rh_t1t2=$(mrstats -mask $Thal_rh -output median $T1T2)
        Thal_lh_t1t2=$(mrstats -mask $Thal_lh -output median $T1T2)
        CSF_lateral_t1t2=$(mrstats -mask $CSF_lateral -output median $T1T2)
        CSF_3rd_t1t2=$(mrstats -mask $CSF_3rd -output median $T1T2)
        CSF_4th_t1t2=$(mrstats -mask $CSF_4th -output median $T1T2)
        MSlesion_t1t2=$(mrstats -mask $MSlesion -output median $T1T2)
    fi

    NAGM_t1flair="NA"
    NAWM_t1flair="NA"
    Thal_rh_t1flair2="NA"
    Thal_lh_t1flair="NA"
    CSF_lateral_t1flair="NA"
    CSF_3rd_t1flair="NA"
    CSF_4th_t1flair="NA"
    MSlesion_t1flair="NA"
    if [ -f $T1FLAIR ]; then
        NAWM_t1flair=$(mrstats -mask $NAWM -output median $T1FLAIR)
        NAGM_t1flair=$(mrstats -mask $NAGM -output median $T1FLAIR)
        Thal_rh_t1flair=$(mrstats -mask $Thal_rh -output median $T1FLAIR)
        Thal_lh_t1flair=$(mrstats -mask $Thal_lh -output median $T1FLAIR)
        CSF_lateral_t1flair=$(mrstats -mask $CSF_lateral -output median $T1FLAIR)
        CSF_3rd_t1flair=$(mrstats -mask $CSF_3rd -output median $T1FLAIR)
        CSF_4th_t1flair=$(mrstats -mask $CSF_4th -output median $T1FLAIR)
        MSlesion_t1flair=$(mrstats -mask $MSlesion -output median $T1FLAIR)
    fi

    Volume_NAWM=$(mrstats -ignorezero -output count $NAWM)
    Volume_NAGM=$(mrstats -ignorezero -output count $NAGM)
    Volume_CSF=$(mrstats -ignorezero -output count $CSF)
    Volume_MSlesions=$(mrstats -ignorezero -output count $MSlesion)

    # write the stats to a .csv file
    echo " saving"
    echo "$participant_and_session, \
        $NAWM_mtr, $NAWM_t1t2, $NAWM_t1flair, \
        $NAGM_mtr, $NAGM_t1t2, $NAGM_t1flair, \
        $Thal_rh_mtr, $Thal_rh_t1t2, $Thal_rh_t1flair, \
        $Thal_lh_mtr, $Thal_lh_t1t2, $Thal_lh_t1flair, \
        $CSF_lateral_mtr, $CSF_lateral_t1t2, $CSF_lateral_t1flair, \
        $CSF_3rd_mtr, $CSF_3rd_t1t2, $CSF_3rd_t1flair, \
        $CSF_4th_mtr, $CSF_4th_t1t2, $CSF_4th_t1flair, \
        $MSlesion_mtr, $MSlesion_t1t2, $MSlesion_t1flair, \
        $Volume_NAWM, $Volume_NAGM, $Volume_CSF, $Volume_MSlesions" >> $outdir/stats/$my_results_file

}

# --- MAIN ---
printf "\n\n\n"

# here we give the data
if [ $auto -eq 0 ]; then
    if [ -z "$session" ]; then
        fullsession1=""
        fullsession2=""
    else
        fullsession1="ses-${session}/"
        fullsession2="ses-${session}_"
    fi
    datadir="$cwd/BIDS/sub-${participant}/${fullsession1}anat"
    T1w_all=("$datadir/sub-${participant}_${fullsession2}T1w.nii.gz")
else
    T1w_all=($(find BIDS -type f -name "*T1w.nii.gz" | sort ))
fi

#echo $session
#echo $fullsession1
#echo $fullsession2
#echo $datadir
#echo $T1w_all



if [ $type_sel -eq 1 ]; then
    type="calib-none"
elif [ $type_sel -eq 2 ]; then
    type="calib-lin"
elif  [ $type_sel -eq 3 ]; then
    type="calib-nonlin"
elif  [ $type_sel -eq 4 ]; then
    type="calib-nonlin2"
elif  [ $type_sel -eq 5 ]; then
    type="calib-nonlin2b"
elif  [ $type_sel -eq 6 ]; then
    type="calib-nonlin3"
else
    exit
fi


for test_T1w in ${T1w_all[@]}; do

    base0=${test_T1w##*/};base=${base0%_T1w*}
    local_participant=${base%_ses*}
    local_session="ses-${base##*ses-}"
    outdir=$outputdir/$local_participant/$local_session
    mkdir -p $outdir/stats
    check_done="$outdir/stats/${base}_stats.done"

    #if [ ! -f $check_done ];then

        participant_and_session=$base
        echo "Processing $participant_and_session"

        my_results_file="${participant_and_session}_${type}_results.csv"
        KUL_create_results_file

        KUL_compute_stats
    
        touch $check_done
    
    #else
    #    echo " $base already done".
    #fi

done


