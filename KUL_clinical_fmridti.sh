#!/bin/bash -e
# Bash shell script to analyse clinical fMRI/DTI
#
# Requires matlab fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 16/12/2020
version="0.1"

kul_main_dir=`dirname "$0"`
source $kul_main_dir/KUL_main_functions.sh
cwd=$(pwd)

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` is a batch analysis of clinical fMRI/DTI data

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -d DICOM/JohnDoe.zip

Required arguments:

     -p:  participant name

Optional arguments:

     -d:  dicom zip file (or directory)
     -v:  show output from commands


USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given

# Set required options
p_flag=0
d_flag=0 

if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:d:v" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        d) #dicomzip
			dicomzip=$OPTARG
            d_flag=1
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

# MRTRIX verbose or not?
if [ $silent -eq 1 ] ; then

	export MRTRIX_QUIET=1

fi

# --- functions ---

function KUL_antsApply_Transform {

    antsApplyTransforms -d 3 --float 1 \
        --verbose 1 \
        -i $input \
        -o $output \
        -r $reference \
        -t $transform \
        -n Linear
}

function KUL_convert2bids {
    # convert the DICOM to BIDS
    if [ ! -d "BIDS/sub-${participant}" ];then
        KUL_dcm2bids_new.sh -d $dicomzip -p ${participant} -c study_config/sequences.txt -e
    else
        echo "BIDS conversion already done"
    fi
}

function KUL_run_fmriprep {
    if [ ! -f fmriprep/sub-${participant}.html ]; then
        cp study_config/run_fmriprep.txt KUL_LOG/sub-${participant}_run_fmriprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_fmriprep.txt 
        rm -fr fmriprep_work_${participant}
    else
        echo "fmriprep already done"
    fi
}

function KUL_run_dwiprep {
    if [ ! -f dwiprep/sub-${participant}/dwiprep_is_done.log ]; then
        cp study_config/run_dwiprep.txt KUL_LOG/sub-${participant}_run_dwiprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_dwiprep.txt
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_dwiprep.txt 
    else
        echo "dwiprep already done"
    fi
}

function KUL_run_freesurfer {
    if [ ! -f freesurfer/sub-${participant}.done ]; then
        cp study_config/run_freesurfer.txt KUL_LOG/sub-${participant}_run_freesurfer.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_freesurfer.txt
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_freesurfer.txt 
    else
        echo "freesurfer already done"
    fi
}

function KUL_compute_SPM_matlab {
    fmriresults="$computedir/RESULTS/stats_$shorttask"
    mkdir -p $fmriresults
    pcf="${scriptsdir}/stats_${shorttask}.m" #participant config file
    pjf="${scriptsdir}/stats_${shorttask}_job.m" #participant job file
    # get rid of - in filename, since this breaks -r in matlab
    pcf=${pcf/run-/run}
    pjf=${pjf/run-/run}
    #echo "$pcf -- $pjf"
    cp $tcf $pcf
    cp $tjf $pjf
    sed -i.bck "s|###JOBFILE###|$pjf|" $pcf
    sed -i.bck "s|###FMRIDIR###|$fmridatadir|" $pjf
    sed -i.bck "s|###FMRIFILE###|$fmrifile|" $pjf
    sed -i.bck "s|###FMRIRESULTS###|$fmriresults|" $pjf
    $matlab_exe -nodisplay -nosplash -nodesktop -r "run('$pcf');exit;"
            
    result=$computedir/RESULTS/MNI/${shorttask}_space-MNI152NLin6Asym.nii
    cp $fmriresults/spmT_0001.nii $result
            
    result_global=$cwd/RESULTS/sub-$participant/SPM_${shorttask}.nii
            
    # since SPM analysis was in MNI space, we transform back in native space
    input=$result
    output=$result_global
    transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5
    reference=$result
    #echo "input=$input"
    #echo "output=$output"
    #echo "transform=$transform"
    #echo "reference=$reference"
    KUL_antsApply_Transform

    gm_mask="$fmriprepdir/../anat/sub-${participant}_label-GM_probseg.nii.gz"
    gm_mask2=$computedir/RESULTS/gm_mask_${shorttask}.nii.gz
    gm_result_global=$cwd/RESULTS/sub-$participant/SPM_${shorttask}_gm.nii
    mrgrid $gm_mask regrid -template $result_global $gm_mask2
    mrcalc $result_global $gm_mask2 0.3 -gt -mul $gm_result_global

} 

function KUL_compute_SPM {
    #  setup variables
    computedir="$cwd/compute/SPM/sub-$participant"
    fmridatadir="$computedir/fmridata"
    scriptsdir="$computedir/scripts"
    fmriprepdir="fmriprep/sub-$participant/func"
    searchtask="_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii"
    matlab_exe=$(which matlab)
    #  the template files in KNT for SPM analysis
    tcf="$kul_main_dir/share/spm12/spm12_fmri_stats_1run.m" #template config file
    tjf="$kul_main_dir/share/spm12/spm12_fmri_stats_1run_job.m" #template job file

    mkdir -p $fmridatadir
    mkdir -p $scriptsdir
    mkdir -p $computedir/RESULTS/MNI


    # Provide the anatomy
    cp -f $fmriprepdir/../anat/sub-${participant}_desc-preproc_T1w.nii.gz $globalresultsdir/T1w.nii.gz
    gunzip -f $globalresultsdir/T1w.nii.gz

    if [ ! -f KUL_LOG/${participant}_SPM.done ]; then
        echo "Preparing for SPM"
        tasks=( $(find $fmriprepdir -name "*${searchtask}.gz" -type f) )
        #echo ${tasks[@]}

        # we loop over the found tasks
        for task in ${tasks[@]}; do
            d1=${task#*_task-}
            shorttask=${d1%_space*}
            #echo "$task -- $shorttask"
            if [[ ! "$shorttask" = *"rest"* ]]; then
                echo " Analysing task $shorttask"
                fmrifile="${shorttask}${searchtask}"
                cp $fmriprepdir/*$fmrifile.gz $fmridatadir
                gunzip -f $fmridatadir/*$fmrifile.gz
                KUL_compute_SPM_matlab

                # do the combined analysis
                if [[ "$shorttask" == *"run-2" ]]; then
                    echo "this is run2, we run full analysis now" 
                    shorttask=${shorttask%_run-2}
                    #echo $shorttask
                    tcf="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs.m" #template config file
                    tjf="$kul_main_dir/share/spm12/spm12_fmri_stats_2runs_job.m" #template job file
                    fmrifile="${shorttask}"
                    KUL_compute_SPM_matlab
                fi
            
            fi
        done
        echo "Done" > KUL_LOG/${participant}_SPM.done
    else
        echo "SPM analysis already done"
    fi
}

function KUL_segment_tumor {
    bidsdir="BIDS/sub-$participant"
    T1w=($(find $bidsdir -name "*T1w.nii.gz" ! -name "*gadolinium*" -type f ))
    nT1w=${#T1w[@]}
    echo "number of non-contrast T1w: $nT1w"
    cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
    ncT1w=${#cT1w[@]}
    echo "number of contrast enhanced T1w: $ncT1w"
    FLAIR=($(find $bidsdir -name "*FLAIR.nii.gz" -type f ))
    nFLAIR=${#FLAIR[@]}
    echo "number of FLAIR: $nFLAIR"
    T2w=($(find $bidsdir -name "*T2w.nii.gz" -type f ))
    nT2w=${#T2w[@]}
    echo "number of T2w: $nT2w"

    mkdir -p $globalresultsdir
    
    # this will segment the lesion automatically
    hdgliodir="compute/hdglio/input"
    if [ $hdglio -eq 1 ]; then

        if [ $nT1w -eq 1 ] && [ $ncT1w -eq 1 ] && [ $nFLAIR -eq 1 ] && [ $nT2w -eq 1 ];then
            mkdir -p $hdgliodir
            mkdir -p compute/hdglio/output
            if [ ! -f "compute/hdglio/output/volumes.txt" ]; then
                cp $T1w $hdgliodir/T1.nii.gz
                cp $cT1w $hdgliodir/CT1.nii.gz
                cp $FLAIR $hdgliodir/FLAIR.nii.gz
                cp $T2w $hdgliodir/T2.nii.gz
                echo "Running HD-GLIO-AUTO using docker"
                docker run --gpus all --mount type=bind,source=${cwd}/compute/hdglio/input,target=/input \
                 --mount type=bind,source=${cwd}/compute/hdglio/output,target=/output \
                 jenspetersen/hd-glio-auto
                mrcalc compute/hdglio/output/segmentation.nii.gz 1 -ge $globalresultsdir/lesion.nii.gz
            else
                echo "HD-GLIO-AUTO already done"
            fi
        else
            echo "Not possible to run HD-GLIO-AUTO"
        fi 
    fi
}

function KUL_run_VBG {
    if [ $vbg -eq 1 ]; then
        vgb_test="lesion_wf/output_LWF/sub-${participant}/sub-${participant}_aparc+aseg.nii.gz"
        if [[ ! -f $vbg_test ]]; then
            KUL_VBG.sh -p ${participant} -l $globalresultsdir/lesion.nii -z T1 -b -B 1 -t -F -n 12 -v
            mkdir -p freesurfer
            ln -s ${cwd}/lesion_wf/output_LWF/sub-${participant}/sub-${participant}_FS_output/sub-${participant}/ freesurfer
            echo "done" > freesurfer/sub-${participant}_freesurfer_is.done
        else
            echo "KUL_VBG has already run"
        fi
    fi
}


# --- MAIN ---
hdglio=1
vbg=1
ncpu=12
globalresultsdir=$cwd/RESULTS/sub-$participant

# STEP 1 - BIDS conversion
KUL_convert2bids

read -p "Are you happy? (y/n) " answ
if [[ ! "$answ" == "y" ]]; then
    exit 1
fi

# STEP 2 - run fmriprep/dwiprep and continue
KUL_run_fmriprep &
KUL_run_dwiprep &

# STEP 3 - run HD-GLIO-AUTO
KUL_segment_tumor

# STEP 4 - run VBG+freesurfer or freesurfer only
if [ $vbg -eq 1 ];then
    KUL_run_VBG &
else
    KUL_run_freesurfer &
fi

# WAIT FOR ALL TO FINISH
wait

KUL_compute_SPM &

#KUL_dwiprep_anat.sh -p $participant -n $ncpu
#KUL_dwiprep_fibertract.sh -p $participant -n $ncpu -c study_config/tracto_tracts.csv -r study_config/tracto_rois.csv -w dhollander_wmfod_reg2T1w -v

# run FSL Melodic
computedir="$cwd/compute/FSL/melodic/sub-$participant"
fmridatadir="$computedir/fmridata"
scriptsdir="$computedir/scripts"
fmriprepdir="fmriprep/sub-$participant/func"
globalresultsdir=$cwd/RESULTS/sub-$participant
searchtask="_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii"

mkdir -p $fmridatadir
mkdir -p $computedir/RESULTS
mkdir -p $globalresultsdir

if [ ! -f KUL_LOG/${participant}_melodic.done ]; then
    echo "Preparing for Melodic"
    tasks=( $(find $fmriprepdir -name "*${searchtask}.gz" -type f -printf '%P\n') )
    # we loop over the found tasks
    for task in ${tasks[@]}; do
        d1=${task#*_task-}
        shorttask=${d1%_space*}
        #echo "$task -- $shorttask"
        echo " Analysing task $shorttask"
        fmrifile="${shorttask}${searchtask}"
        cp $fmriprepdir/*$fmrifile.gz $fmridatadir
        gunzip $fmridatadir/*$fmrifile.gz
        fmriresults="$computedir/stats_$shorttask"
        mkdir -p $fmriresults
        melodic_in="$fmridatadir/sub-${participant}_task-$fmrifile"
        # find the TR
        tr=$(mrinfo $melodic_in -spacing | cut -d " " -f 4)
        # make model and contrast
        dyn=$(mrinfo $melodic_in -size | cut -d " " -f 4)
        t_glm_con="$kul_main_dir/share/FSL/fsl_glm.con"
        t_glm_mat="$kul_main_dir/share/FSL/fsl_glm_${dyn}dyn.mat"        
        #melodic -i Melodic/sub-Croes/fmridata/sub-Croes_task-LIP_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii -o test/ --report --Tdes=glm.mat --Tcon=glm.con
        melodic -i $melodic_in -o $fmriresults --report --tr=$tr --Tdes=$t_glm_mat --Tcon=$t_glm_con --Oall
        # now we compare to known networks
        mkdir -p $fmriresults/kul
        fslcc --noabs -p 3 -t .204 $kul_main_dir/atlasses/Yeo2011_rsfMRI_in_FSL_Space/yeo2011_7_liberal_combined.nii.gz \
         $fmriresults/melodic_IC.nii.gz > $fmriresults/kul/kul_networks.txt

        while IFS=$' ' read network ic stat; do
            echo $network
            network_name=$(sed "${network}q;d" $kul_main_dir/atlasses/Yeo2011_rsfMRI_in_FSL_Space/yeo2011_7_liberal_combined_networks.txt)
            echo $network_name
            icfile="$fmriresults/stats/thresh_zstat${ic}.nii.gz"
            network_file="$fmriresults/kul/melodic_${network_name}_ic${ic}.nii.gz"
            echo $icfile
            echo $network_file
            mrcalc $icfile 2 -gt $icfile -mul $network_file
            # since Melodic analysis was in MNI space, we transform back in native space
            input=$network_file
            output=$globalresultsdir/melodic_${shorttask}_${network_name}_ic${ic}.nii
            transform=${cwd}/fmriprep/sub-${participant}/anat/sub-${participant}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5
            reference=$network_file
            echo "input=$input"
            echo "output=$output"
            echo "transform=$transform"
            echo "reference=$reference"
            KUL_antsApply_Transform
        done < $fmriresults/kul/kul_networks.txt
    done
    echo "Done" > KUL_LOG/${participant}_melodic.done
else
    echo "Melodic analysis already done"
fi
