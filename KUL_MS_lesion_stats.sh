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
     -n:  number of threads to use
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
ncpu=15

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


# verbose or not?
if [ $silent -eq 1 ] ; then
	export MRTRIX_QUIET=1
fi

# --- FUNCTIONS ---
function KUL_create_results_file {
    if [ ! -f $outdir/stats/$my_results_file ];then
        touch $outdir/stats/$my_results_file
        echo "participant_and_session, \
        NAWM_lh_mtr, NAWM_lh_t1t2, NAWM_lh_t1flair, \
        NAWM_rh_mtr, NAWM_rh_t1t2, NAWM_rh_t1flair, \
        NAGM_lh_mtr, NAGM_lh_t1t2, NAGM_lh_t1flair, \
        NAGM_rh_mtr, NAGM_rh_t1t2, NAGM_rh_t1flair, \
        Thal_lh_mtr, Thal_lh_t1t2, Thal_lh_t1flair, \
        Thal_rh_mtr, Thal_rh_t1t2, Thal_rh_t1flair, \
        CSF_lateral_lh_mtr, CSF_lateral_lh_t1t2, CSF_lateral_lh_t1flair, \
        CSF_lateral_rh_mtr, CSF_lateral_rh_t1t2, CSF_lateral_rh_t1flair, \
        CSF_3rd_mtr, CSF_3rd_t1t2, CSF_3rd_t1flair, \
        CSF_4th_mtr, CSF_4th_t1t2, CSF_4th_t1flair, \
        MSlesion_mtr, MSlesion_t1t2, MSlesion_t1flair, \
        CC_Posterior_mtr, CC_Posterior_t1t2, CC_Posterior_t1flair, \
        CC_Mid_Posterior_mtr, CC_Mid_Posterior_t1t2, CC_Mid_Posterior_t1flair, \
        CC_Central_mtr, CC_Central_t1t2, CC_Central_t1flair, \
        CC_Mid_Anterior_mtr, CC_Mid_Anterior_t1t2, CC_Mid_Anterior_t1flair, \
        CC_Anterior_mtr, CC_Anterior_t1t2, CC_Anterior_t1flair, \
        Volume_NAWM_lh, Volume_NAWM_rh, \
        Volume_NAGM_lh, Volume_NAGM_rh, \
        Volume_Thalamus_lh, Volume_Thalamus_rh, \
        Volume_CSF, Volume_CSF_lateral_lh, Volume_CSF_lateral_rh, \
        Volume_CSF_3rd, Volume_CSF_4th, \
        Volume_MSLesions, \
        Volume_CC_Posterior, Volume_CC_Mid_Posterior, Volume_CC_Central, \
        Volume_CC_Mid_Anterior, Volume_CC_Anterior, \
        Volume_TIV" > $outdir/stats/$my_results_file
    fi
}

function KUL_compute_stats {
    # define the input images
    SamSeg="$cwd/$outdir/samseg/seg.mgz"
    fastfs="$cwd/$outdir/fs/${participant_and_session}/mri/aparc.DKTatlas+aseg.deep.withCC.mgz"
    MTR="$outdir/${participant_and_session}_ratio-MTC.nii.gz"
    T1T2="$outdir/${participant_and_session}_ratio-T1T2w_${type}.nii.gz"
    T1FLAIR="$outdir/${participant_and_session}_ratio-T1FLAIR_${type}.nii.gz"
    #echo $MTR
    #echo $T1T2

    # define the output images
    NAWM_lh="$cwd/$outdir/rois/${participant_and_session}_NAWM_lh.nii.gz"
    NAWM_rh="$cwd/$outdir/rois/${participant_and_session}_NAWM_rh.nii.gz"
    NAGM_lh="$cwd/$outdir/rois/${participant_and_session}_NAGM_lh.nii.gz"
    NAGM_rh="$cwd/$outdir/rois/${participant_and_session}_NAGM_rh.nii.gz"
    Thal_lh="$cwd/$outdir/rois/${participant_and_session}_Thal_lh.nii.gz"
    Thal_rh="$cwd/$outdir/rois/${participant_and_session}_Thal_rh.nii.gz"
    CC_Posterior="$cwd/$outdir/rois/${participant_and_session}_CC_Posterior.nii.gz"
    CC_Mid_Posterior="$cwd/$outdir/rois/${participant_and_session}_CC_Mid_Posterior.nii.gz"
    CC_Central="$cwd/$outdir/rois/${participant_and_session}_CC_Central.nii.gz"
    CC_Mid_Anterior="$cwd/$outdir/rois/${participant_and_session}_CC_Mid_Anterior.nii.gz"
    CC_Anterior="$cwd/$outdir/rois/${participant_and_session}_CC_Anterior.nii.gz"
    CSF_lateral_lh="$cwd/$outdir/rois/${participant_and_session}_CSF_lateral_lh.nii.gz"
    CSF_lateral_rh="$cwd/$outdir/rois/${participant_and_session}_CSF_lateral_rh.nii.gz"
    CSF_3rd="$cwd/$outdir/rois/${participant_and_session}_CSF_3rd.nii.gz"
    CSF_4th="$cwd/$outdir/rois/${participant_and_session}_CSF_4th.nii.gz"
    MSlesion="$cwd/$outdir/rois/${participant_and_session}_MSLesion.nii.gz"

    echo " making VOIs"
    # do the computation of the masks
    if [ ! -f $CC_Anterior ]; then
        #mrcalc $SamSeg 99 -eq $MSlesion -force # already computed in T1T2FLAIRMTR script
        mrcalc $SamSeg 2 -eq $NAWM_lh -force -nthreads $ncpu
        mrcalc $SamSeg 41 -eq $NAWM_rh -force -nthreads $ncpu
        mrcalc $SamSeg 3 -eq $NAGM_lh -force -nthreads $ncpu
        mrcalc $SamSeg 42 -eq $NAGM_rh -force -nthreads $ncpu
        mrcalc $SamSeg 10 -eq $Thal_lh -force -nthreads $ncpu
        mrcalc $SamSeg 49 -eq $Thal_rh -force -nthreads $ncpu
        mrcalc $SamSeg 4 -eq $CSF_lateral_lh -force -nthreads $ncpu
        mrcalc $SamSeg 43 -eq $CSF_lateral_rh -force -nthreads $ncpu
        mrcalc $SamSeg 14 -eq $CSF_3rd -force -nthreads $ncpu
        mrcalc $SamSeg 15 -eq $CSF_4th -force -nthreads $ncpu
        
        # note fastsurfer images have a different FOV. Regridding to samseg.
        mrcalc $fastfs 251 -eq - | mrgrid -template $NAWM_lh - regrid - | mrcalc - 0.9 -gt $CC_Posterior -force -nthreads $ncpu
        mrcalc $fastfs 252 -eq - | mrgrid -template $NAWM_lh - regrid - | mrcalc - 0.9 -gt $CC_Mid_Posterior -force -nthreads $ncpu
        mrcalc $fastfs 253 -eq - | mrgrid -template $NAWM_lh - regrid - | mrcalc - 0.9 -gt $CC_Central -force -nthreads $ncpu
        mrcalc $fastfs 254 -eq - | mrgrid -template $NAWM_lh - regrid - | mrcalc - 0.9 -gt $CC_Mid_Anterior -force -nthreads $ncpu
        mrcalc $fastfs 255 -eq - | mrgrid -template $NAWM_lh - regrid - | mrcalc - 0.9 -gt $CC_Anterior -force -nthreads $ncpu
    fi

    echo " computing stats"
    # do the stats
    if [ -f $MTR ]; then
        NAGM_lh_mtr=$(mrstats -mask $NAGM_lh -output median $MTR -nthreads $ncpu)
        NAGM_rh_mtr=$(mrstats -mask $NAGM_rh -output median $MTR -nthreads $ncpu)
        NAWM_lh_mtr=$(mrstats -mask $NAWM_lh -output median $MTR -nthreads $ncpu)
        NAWM_rh_mtr=$(mrstats -mask $NAWM_rh -output median $MTR -nthreads $ncpu)
        Thal_rh_mtr=$(mrstats -mask $Thal_rh -output median $MTR -nthreads $ncpu)
        Thal_lh_mtr=$(mrstats -mask $Thal_lh -output median $MTR -nthreads $ncpu)
        CC_Posterior_mtr=$(mrstats -mask $CC_Posterior -output median $MTR -nthreads $ncpu)
        CC_Mid_Posterior_mtr=$(mrstats -mask $CC_Mid_Posterior -output median $MTR -nthreads $ncpu)
        CC_Central_mtr=$(mrstats -mask $CC_Central -output median $MTR -nthreads $ncpu)
        CC_Mid_Anterior_mtr=$(mrstats -mask $CC_Mid_Anterior -output median $MTR -nthreads $ncpu)
        CC_Anterior_mtr=$(mrstats -mask $CC_Anterior -output median $MTR -nthreads $ncpu)
        CSF_lateral_lh_mtr=$(mrstats -mask $CSF_lateral_lh -output median $MTR -nthreads $ncpu)
        CSF_lateral_rh_mtr=$(mrstats -mask $CSF_lateral_rh -output median $MTR -nthreads $ncpu)
        CSF_3rd_mtr=$(mrstats -mask $CSF_3rd -output median $MTR -nthreads $ncpu)
        CSF_4th_mtr=$(mrstats -mask $CSF_4th -output median $MTR -nthreads $ncpu)
        MSlesion_mtr=$(mrstats -mask $MSlesion -output median $MTR -nthreads $ncpu)
    fi

    if [ -f $T1T2 ]; then
        NAGM_lh_t1t2=$(mrstats -mask $NAGM_lh -output median $T1T2 -nthreads $ncpu)
        NAGM_rh_t1t2=$(mrstats -mask $NAGM_rh -output median $T1T2 -nthreads $ncpu)
        NAWM_lh_t1t2=$(mrstats -mask $NAWM_lh -output median $T1T2 -nthreads $ncpu)
        NAWM_rh_t1t2=$(mrstats -mask $NAWM_rh -output median $T1T2 -nthreads $ncpu)
        Thal_rh_t1t2=$(mrstats -mask $Thal_rh -output median $T1T2 -nthreads $ncpu)
        Thal_lh_t1t2=$(mrstats -mask $Thal_lh -output median $T1T2 -nthreads $ncpu)
        CC_Posterior_t1t2=$(mrstats -mask $CC_Posterior -output median $T1T2 -nthreads $ncpu)
        CC_Mid_Posterior_t1t2=$(mrstats -mask $CC_Mid_Posterior -output median $T1T2 -nthreads $ncpu)
        CC_Central_t1t2=$(mrstats -mask $CC_Central -output median $T1T2 -nthreads $ncpu)
        CC_Mid_Anterior_t1t2=$(mrstats -mask $CC_Mid_Anterior -output median $T1T2 -nthreads $ncpu)
        CC_Anterior_t1t2=$(mrstats -mask $CC_Anterior -output median $T1T2 -nthreads $ncpu)
        CSF_lateral_lh_t1t2=$(mrstats -mask $CSF_lateral_lh -output median $T1T2 -nthreads $ncpu)
        CSF_lateral_rh_t1t2=$(mrstats -mask $CSF_lateral_rh -output median $T1T2 -nthreads $ncpu)
        CSF_3rd_t1t2=$(mrstats -mask $CSF_3rd -output median $T1T2 -nthreads $ncpu)
        CSF_4th_t1t2=$(mrstats -mask $CSF_4th -output median $T1T2 -nthreads $ncpu)
        MSlesion_t1t2=$(mrstats -mask $MSlesion -output median $T1T2 -nthreads $ncpu)
    fi

    if [ -f $T1FLAIR ]; then
        NAGM_lh_t1flair=$(mrstats -mask $NAGM_lh -output median $T1FLAIR -nthreads $ncpu)
        NAGM_rh_t1flair=$(mrstats -mask $NAGM_rh -output median $T1FLAIR -nthreads $ncpu)
        NAWM_lh_t1flair=$(mrstats -mask $NAWM_lh -output median $T1FLAIR -nthreads $ncpu)
        NAWM_rh_t1flair=$(mrstats -mask $NAWM_rh -output median $T1FLAIR -nthreads $ncpu)
        Thal_rh_t1flair=$(mrstats -mask $Thal_rh -output median $T1FLAIR -nthreads $ncpu)
        Thal_lh_t1flair=$(mrstats -mask $Thal_lh -output median $T1FLAIR -nthreads $ncpu)
        CC_Posterior_t1flair=$(mrstats -mask $CC_Posterior -output median $T1FLAIR -nthreads $ncpu)
        CC_Mid_Posterior_t1flair=$(mrstats -mask $CC_Mid_Posterior -output median $T1FLAIR -nthreads $ncpu)
        CC_Central_t1flair=$(mrstats -mask $CC_Central -output median $T1FLAIR -nthreads $ncpu)
        CC_Mid_Anterior_t1flair=$(mrstats -mask $CC_Mid_Anterior -output median $T1FLAIR -nthreads $ncpu)
        CC_Anterior_t1flair=$(mrstats -mask $CC_Anterior -output median $T1FLAIR -nthreads $ncpu)
        CSF_lateral_lh_t1flair=$(mrstats -mask $CSF_lateral_lh -output median $T1FLAIR -nthreads $ncpu)
        CSF_lateral_rh_t1flair=$(mrstats -mask $CSF_lateral_rh -output median $T1FLAIR -nthreads $ncpu)
        CSF_3rd_t1flair=$(mrstats -mask $CSF_3rd -output median $T1FLAIR -nthreads $ncpu)
        CSF_4th_t1flair=$(mrstats -mask $CSF_4th -output median $T1FLAIR -nthreads $ncpu)
        MSlesion_t1flair=$(mrstats -mask $MSlesion -output median $T1FLAIR -nthreads $ncpu)
    fi

    # read some volumes from fastsurfer
    fastfs_stats="$cwd/$outdir/fs/${participant_and_session}/stats/aparc.DKTatlas+aseg.deep.volume.stats"
    s_CC_Posterior="CC_Posterior"
    t=$(grep $s_CC_Posterior $fastfs_stats) 
    Volume_CC_Posterior=$(echo $t | cut -d " " -f 4 )
    s_CC_Mid_Posterior="CC_Mid_Posterior"
    t=$(grep $s_CC_Mid_Posterior $fastfs_stats) 
    Volume_CC_Mid_Posterior=$(echo $t | cut -d " " -f 4 )
    s_CC_Central="CC_Central"
    t=$(grep $s_CC_Central $fastfs_stats) 
    Volume_CC_Central=$(echo $t | cut -d " " -f 4 )
    s_CC_Mid_Anterior="CC_Mid_Anterior"
    t=$(grep $s_CC_Mid_Anterior $fastfs_stats) 
    Volume_CC_Mid_Anterior=$(echo $t | cut -d " " -f 4 )
    s_CC_Anterior="CC_Anterior"
    t=$(grep $s_CC_Anterior $fastfs_stats) 
    Volume_CC_Anterior=$(echo $t | cut -d " " -f 4 )
    
    # read more volumes, now from samseg
    samseg_stats="$cwd/$outdir/samseg/samseg.stats"
    s_CSF="CSF"
    t=$(grep $s_CSF $samseg_stats)
    Volume_CSF=$(echo $t | cut -d "," -f 2 )
    s_CSF_4th="4th-Ventricle"
    t=$(grep $s_CSF_4th $samseg_stats)
    Volume_CSF_4th=$(echo $t | cut -d "," -f 2 )
    s_CSF_3rd="3rd-Ventricle"
    t=$(grep $s_CSF_3rd $samseg_stats)
    Volume_CSF_3rd=$(echo $t | cut -d "," -f 2 )
    s_CSF_lateral_lh="Left-Lateral-Ventricle"
    t=$(grep $s_CSF_lateral_lh $samseg_stats)
    Volume_CSF_lateral_lh=$(echo $t | cut -d "," -f 2 )
    s_CSF_lateral_rh="Right-Lateral-Ventricle"
    t=$(grep $s_CSF_lateral_rh $samseg_stats)
    Volume_CSF_lateral_rh=$(echo $t | cut -d "," -f 2 )
    s_NAGM_lh="Left-Cerebral-Cortex"
    t=$(grep $s_NAGM_lh $samseg_stats)
    Volume_NAGM_lh=$(echo $t | cut -d "," -f 2 )
    s_NAGM_rh="Right-Cerebral-Cortex"
    t=$(grep $s_NAGM_rh $samseg_stats)
    Volume_NAGM_rh=$(echo $t | cut -d "," -f 2 )
    s_NAWM_lh="Left-Cerebral-White-Matter"
    t=$(grep $s_NAWM_lh $samseg_stats) 
    Volume_NAWM_lh=$(echo $t | cut -d "," -f 2 )
    s_NAWM_rh="Right-Cerebral-White-Matter"
    t=$(grep $s_NAWM_rh $samseg_stats) 
    Volume_NAWM_rh=$(echo $t | cut -d "," -f 2 )
    s_MSLesions="Lesions"
    t=$(grep $s_MSLesions $samseg_stats) 
    Volume_MSLesions=$(echo $t | cut -d "," -f 2 )
    s_Thalamus_lh="Left-Thalamus"
    t=$(grep $s_Thalamus_lh $samseg_stats) 
    Volume_Thalamus_lh=$(echo $t | cut -d "," -f 2 )
    s_Thalamus_rh="Right-Thalamus"
    t=$(grep $s_Thalamus_rh $samseg_stats) 
    Volume_Thalamus_rh=$(echo $t | cut -d "," -f 2 )

    # read more volumes, now from samseg tiv
    samseg_stats2="$cwd/$outdir/samseg/sbtiv.stats"
    s_TIV="Intra-Cranial"
    t=$(grep $s_TIV $samseg_stats2) 
    Volume_TIV=$(echo $t | cut -d "," -f 2 )

    # write the stats to a .csv file
    echo " saving"
    echo "$participant_and_session, \
        $NAWM_lh_mtr, $NAWM_lh_t1t2, $NAWM_lh_t1flair, \
        $NAWM_rh_mtr, $NAWM_rh_t1t2, $NAWM_rh_t1flair, \
        $NAGM_lh_mtr, $NAGM_lh_t1t2, $NAGM_lh_t1flair, \
        $NAGM_rh_mtr, $NAGM_rh_t1t2, $NAGM_rh_t1flair, \
        $Thal_lh_mtr, $Thal_lh_t1t2, $Thal_lh_t1flair, \
        $Thal_rh_mtr, $Thal_rh_t1t2, $Thal_rh_t1flair, \
        $CSF_lateral_lh_mtr, $CSF_lateral_lh_t1t2, $CSF_lateral_lh_t1flair, \
        $CSF_lateral_rh_mtr, $CSF_lateral_rh_t1t2, $CSF_lateral_rh_t1flair, \
        $CSF_3rd_mtr, $CSF_3rd_t1t2, $CSF_3rd_t1flair, \
        $CSF_4th_mtr, $CSF_4th_t1t2, $CSF_4th_t1flair, \
        $MSlesion_mtr, $MSlesion_t1t2, $MSlesion_t1flair, \
        $CC_Posterior_mtr, $CC_Posterior_t1t2, $CC_Posterior_t1flair, \
        $CC_Mid_Posterior_mtr, $CC_Mid_Posterior_t1t2, $CC_Mid_Posterior_t1flair, \
        $CC_Central_mtr, $CC_Central_t1t2, $CC_Central_t1flair, \
        $CC_Mid_Anterior_mtr, $CC_Mid_Anterior_t1t2, $CC_Mid_Anterior_t1flair, \
        $CC_Anterior_mtr, $CC_Anterior_t1t2, $CC_Anterior_t1flair, \
        $Volume_NAWM_lh, $Volume_NAWM_rh, \
        $Volume_NAGM_lh, $Volume_NAGM_rh, \
        $Volume_Thalamus_lh, $Volume_Thalamus_rh, \
        $Volume_CSF, $Volume_CSF_lateral_lh, $Volume_CSF_lateral_rh, \
        $Volume_CSF_3rd, $Volume_CSF_4th, \
        $Volume_MSLesions, \
        $Volume_CC_Posterior, $Volume_CC_Mid_Posterior, $Volume_CC_Central, \
        $Volume_CC_Mid_Anterior, $Volume_CC_Anterior, \
        $Volume_TIV" >> $outdir/stats/$my_results_file

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
    type="calib-nonlin3"
else
    echo "Error: wrong type"
    exit
fi

for test_T1w in ${T1w_all[@]}; do

    base0=${test_T1w##*/};base=${base0%_T1w*}
    local_participant=${base%_ses*}
    local_session="ses-${base##*ses-}"
    outdir=$outputdir/$local_participant/$local_session
    mkdir -p $outdir/stats
    check_done="$outdir/stats/${base}_stats.done"

    if [ ! -f $check_done ];then

        participant_and_session=$base
        echo "Processing $participant_and_session"

        my_results_file="${participant_and_session}_${type}_results.csv"
        KUL_create_results_file

        KUL_compute_stats
    
        touch $check_done
    
    else
        echo " $base already done".
    fi

done
