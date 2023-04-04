#!/bin/bash
# Bash shell script to analyse clinical fMRI/DTI
#
# Requires matlab fmriprep
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 07/12/2021
version="0.9"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

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

     -t:  processing type
        1: (DEFAULT) do automatic tumor intra-axial segmentation and vbg (tumor with T1w, cT1w, T2w and FLAIR)
        2: do automatic tumor extra-axial segmentation and vbg (tumor with T1w, cT1w, T2w and FLAIR)
        3: do vbg with manual mask (tumor but missing one of T1w, cT1w, T2w and FLAIR; 
                    put lesion.nii.gz in RESULTS/sub-{participant}/Lesion)
        4: clinical fMRI without glioma (cavernoma, epilepsy, etc... cT1w)
        5: clinical dMRI for DBS of essential tremor (DRT tract)
        6: cliniacl dMRI for DBS of Parkinson's disease (CSHD pathway)
     -d:  dicom zip file (or directory)
     -s:  scaffold (make a default DICOM and study_config)
     -B:  make a backup and cleanup 
     -r:  redo certain steps (program will ask)
     -R:  make results ready
        1: use cT1w as underlay
        2: use FLAIR as underlay
        3: use SWI as underlay
        4: Use T1w as underlay
     -n:  number of threads to use (default 48)
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
ants_verbose=1
ncpu=48
bc=0 
type=1
redo=0
results=0 
verbose_level=1
dbs=0
scaffold=0 

# Set required options
p_flag=0
d_flag=0


if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:t:d:n:v:R:Brs" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        t) #type
			type=$OPTARG
		;;
        d) #dicomzip
			dicomzip=$OPTARG
            d_flag=1
		;;
        n) #ncpu
			ncpu=$OPTARG
		;;
		B) #backup&clean
			bc=1
		;;
        r) #redo
			redo=1
		;;
        s) #scaffold
			scaffold=1
		;;
        R) #make results
			results=$OPTARG
		;;
        v) #verbose
            verbose_level=$OPTARG
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

KUL_LOG_DIR="KUL_LOG/${script}/sub-${participant}"
mkdir -p $KUL_LOG_DIR

# MRTRIX and others verbose or not?
if [ $verbose_level -lt 2 ] ; then
	export MRTRIX_QUIET=1
    silent=1
    str_silent=" > /dev/null 2>&1" 
    ants_verbose=0
elif [ $verbose_level -eq 2 ] ; then
    silent=0
    str_silent="" 
    ants_verbose=1
fi

# Determine what to process depending on patient lesion type
if [ $type -eq 1 ]; then
    hdglio=1
    vbg=1
elif [ $type -eq 2 ]; then
    hdglio=1
    vbg=2
elif [ $type -eq 3 ]; then
    hdglio=0
    vbg=3
elif [ $type -eq 4 ]; then
    hdglio=0
    vbg=0
elif [ $type -eq 5 ]; then
    hdglio=0
    vbg=0
    dbs=1
elif [ $type -eq 6 ]; then
    hdglio=0
    vbg=0
    dbs=2
fi


# GLOBAL defs
globalresultsdir=$cwd/RESULTS/sub-$participant
derivativesdir=${cwd}/BIDS/derivatives/KUL_compute/sub-${participant}

function KUL_scaffold {

    echo "Making scaffold for clinical_sub-${participant}_type${type}"
    mkdir -p $cwd/clinical_sub-${participant}_type${type}/DICOM
    mkdir -p $cwd/clinical_sub-${participant}_type${type}/study_config
    rm -fr $cwd/clinical_sub-${participant}_type${type}/study_config/*
    if [ $type -lt 5 ]; then
        echo "Setting up for a tumor/epilepsy/... patient (type: $type)"
        cp ${kul_main_dir}/study_config/clinical_fmri_dmri/* $cwd/clinical_sub-${participant}_type${type}/study_config
    else
        echo "Setting up for a DBS patient (type: $type)"
        cp ${kul_main_dir}/study_config/clinical_dmri_dbs/* $cwd/clinical_sub-${participant}_type${type}/study_config
    fi

    exit 0

}

# Scaffold
if [ $scaffold -eq 1 ]; then
    KUL_scaffold
fi

# The BACKUP and clean option
if [ $bc -eq 1 ]; then
    # clean some stuff
    clean_dwiprep="./dwiprep/sub-${participant}/sub-${participant}/*dwifsl*tmp* \
        ./dwiprep/sub-${participant}/sub-${participant}/raw \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi_orig* \
        ./dwiprep/sub-${participant}/sub-${participant}/dwi_preproced.mif"
    clean_other="./fmriprep_work* \
        ./BIDS/tmp_dcm2bids"

    rm -fr $clean_dwiprep $clean_other

    while true; do
        read -s -p "Give a password to encrypt the backup with: " password
        echo
        read -s -p "Password (again): " password2
        echo
        [ "$password" = "$password2" ] && break
        echo "Please try again"
    done    
    today=$(date +%Y_%m_%d)
    7z a -bd -y -p${password} -mhe=on -mx=7 ../Finished_${today}_sub-${participant}.7z * 

    exit 0
fi

# The make RESULTS option
if [ $results -gt 0 ];then

    if [ $type -lt 5 ];then

        mrview_tracts[0]="Tract-csd_CST_LT"
        mrview_rgb[0]="0.678,0.847,0.902"
        mrview_tracts[1]="Tract-csd_CST_RT"
        mrview_rgb[1]="0,0,1"
        mrview_tracts[2]="Tract-csd_AF_all_LT"
        mrview_rgb[2]="1,0,0"
        mrview_tracts[3]="Tract-csd_AF_all_RT"
        mrview_rgb[3]="0,1,0"
        mrview_tracts[4]="Tract-csd_CCing_LT"
        mrview_rgb[4]="1,1,0"
        mrview_tracts[5]="Tract-csd_CCing_RT"
        mrview_rgb[5]="1,0.647,0"
        mrview_tracts[6]="Tract-csd_TCing_LT"
        mrview_rgb[6]="1,1,0"
        mrview_tracts[7]="Tract-csd_TCing_RT"
        mrview_rgb[7]="1,0.646,0"
        mrview_tracts[8]="Tract-csd_FAT_LT"
        mrview_rgb[8]="1,0.647,0"
        mrview_tracts[9]="Tract-csd_FAT_RT"
        mrview_rgb[9]="1,1,0"
        mrview_tracts[10]="Tract-csd_ILF_LT"
        mrview_rgb[10]="0,0,1"
        mrview_tracts[11]="Tract-csd_ILF_RT"
        mrview_rgb[11]="0.678,0.847,0.784"
        mrview_tracts[12]="Tract-csd_IFOF_LT"
        mrview_rgb[12]="0.75,0.25,0.75"
        mrview_tracts[13]="Tract-csd_IFOF_RT"
        mrview_rgb[13]="1,0.75,0.8"
        mrview_tracts[14]="Tract-csd_UF_LT"
        mrview_rgb[14]="0,0.784,0"
        mrview_tracts[15]="Tract-csd_UF_RT"
        mrview_rgb[15]="0.784,0,0"
        mrview_tracts[16]="Tract-csd_OR_occlobe_LT"
        mrview_rgb[16]="0.2,0.784,0.4"
        mrview_tracts[17]="Tract-csd_OR_occlobe_RT"
        mrview_rgb[17]="0.784,0.2,0.4"
        mrview_tracts[18]="Tract-csd_MdLF_LT"
        mrview_rgb[18]="0.6,0.784,0.04"
        mrview_tracts[19]="Tract-csd_MdLF_RT"
        mrview_rgb[19]="0.784,0.6,0.04"
        mrview_tracts[20]="Tract-csd_ML_LT"
        mrview_rgb[20]="0,0.5,0.04"
        mrview_tracts[21]="Tract-csd_ML_RT"
        mrview_rgb[21]="0,0.5,0.5"
        ntracts=22

    elif [ $type -eq 5 ];then

        mrview_tracts[0]="Tract-csd_CST_LT"
        mrview_rgb[0]="0.678,0.847,0.902"
        mrview_tracts[1]="Tract-csd_CST_RT"
        mrview_rgb[1]="0,0,1"
        mrview_tracts[2]="Tract-csd_DRT_LT"
        mrview_rgb[2]="1,0,0.23"
        mrview_tracts[3]="Tract-csd_DRT_RT"
        mrview_rgb[3]="0.23,1,0"
        ntracts=4

    elif [ $type -eq 6 ];then

        mrview_tracts[0]="Tract-csd_CST_LT"
        mrview_rgb[0]="0.678,0.847,0.902"
        mrview_tracts[1]="Tract-csd_CST_RT"
        mrview_rgb[1]="0,0,1"
        mrview_tracts[2]="Tract-csd_CSHDP_LT"
        mrview_rgb[2]="0.23,0.12,0"
        mrview_tracts[3]="Tract-csd_CSHDP_RT"
        mrview_rgb[3]="1,0.12,0.20"
        ntracts=4

    fi

    #echo "ntracts: $ntracts"
    result_type=0

    if [ $results -eq 1 ]; then

        underlay=$globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_cT1w"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_cT1w"

    elif [ $results -eq 2 ]; then

        underlay=$globalresultsdir/Anat/FLAIR_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_FLAIR"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_FLAIR"

    elif [ $results -eq 3 ]; then

        underlay=$globalresultsdir/Anat/SWI_reg2_T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_SWI"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_SWI"

    else

        underlay=$globalresultsdir/Anat/T1w.nii.gz
        resultsdir_png="$globalresultsdir/Tracto_figures_T1w"
        resultsdir_dcm="$globalresultsdir/PACS/Tracto_T1w"

    fi

    mrview_resolution=256

    rm -fr $resultsdir_png
    mkdir -p $resultsdir_png
    mkdir -p $resultsdir_dcm
    

    for tract_set_i in $(eval echo "{0..$ntracts..2}"); do

        if [ $tract_set_i -lt $ntracts ]; then
            tract_set=(${mrview_tracts[@]:$tract_set_i:2})
            tractname=${tract_set[0]:0:-3}
            tract_i=tract_set_i
            #echo $tractname 
        else
            tract_set=(${mrview_tracts[@]})
            tractname="Tract-csd_ALL"
            tract_i=0
        fi

        mrview_tck=""
        
        for tract in ${tract_set[@]}; do 
            #echo $tract_set_i
            #echo "$tract ${mrview_rgb[$tract_i]} on $underlay"
            if [ -f $globalresultsdir/Tracto/${tract}.tck ]; then
                mrview_tck="$mrview_tck -tractography.load $globalresultsdir/Tracto/${tract}.tck -tractography.colour ${mrview_rgb[$tract_i]}"
            fi
            tract_i=$(($tract_i+1))
        done
        
        #echo $mrview_tck


        ori[0]="TRA"
        ori[1]="SAG"
        ori[2]="COR"

        for orient in ${ori[@]}; do

            if [[ "$orient" == "TRA" ]]; then
                underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF)}')
            elif [[ "$orient" == "SAG" ]]; then
                underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF-2)}')
            else
                underlay_slices=$(mrinfo $underlay -size | awk '{print $(NF-1)}')
            fi
        

            if [ $result_type -eq 0 ]; then
                i=0
                echo "Making ${tractname}_${orient} on $(basename $underlay)"
                mkdir -p $resultsdir_png/${tractname}_${orient}
                voxel_index="-capture.folder $resultsdir_png/${tractname}_${orient} \
                    -capture.prefix ${tractname}_${orient} -noannotations -orientationlabel 1"
                while [ $i -lt $underlay_slices ]
                do
                    #echo Number: $i
                    if [[ "$orient" == "TRA" ]]; then
                        voxel_index="$voxel_index -voxel 0,0,$i -capture.grab"
                        plane=2
                    elif [[ "$orient" == "SAG" ]]; then
                        voxel_index="$voxel_index -voxel $i,0,0 -capture.grab"
                        plane=0
                    else
                        voxel_index="$voxel_index -voxel 0,$i,0 -capture.grab"
                        plane=1
                    fi    
                    let "i+=1" 
                done
                mode_plane="-mode 1 -plane $plane"
                mrview_exit="-exit"
            else
                voxel_index=""
                mode_plane="-mode 2"
                mrview_exit=""
            fi
            #echo $voxel_index

        

            cmd="mrview -size $mrview_resolution,$mrview_resolution
                -load $underlay \
                $mode_plane \
                -tractography.lighting 1 \
                -tractography.slab 1.5 \
                -tractography.thickness 0.3 \
                $mrview_tck \
                $voxel_index \
                -force \
                $mrview_exit"
            #echo $cmd
            eval $cmd

            if [[ "$mrview_exit" = "-exit" ]];then
                cmd="convert $resultsdir_png/${tractname}_${orient}/${tractname}_${orient}*.png \
                    $resultsdir_png/${tractname}_${orient}/${tractname}_${orient}.tiff"
                #echo $cmd
                eval $cmd

                if [ -f DICOM/smartbrain.dcm ]; then
                    dcmdir="$resultsdir_dcm/${tractname}_${orient}"
                    echo "Making dicoms in $dcmdir"
                    mkdir -p $dcmdir
                    cmd="KUL_nii2dcm.py -s ${tractname}_${orient} \
                        $resultsdir_png/${tractname}_${orient}/${tractname}_${orient}.tiff \
                        DICOM/smartbrain.dcm \
                        $dcmdir"
                    echo $cmd
                    eval $cmd
                fi
            fi
        done
    
    done

    exit

fi

# --- functions ---
function KUL_check_redo {
    if [ $redo -eq 1 ];then
        if [ $type -lt 3 ];then
            read -p "Redo: tumor segmentation? (y/n) " answ
            if [[ "$answ" == "y" ]]; then
                #echo $answ
                #echo "rm ${cwd}/KUL_LOG/sub-${participant}_anat_*.done"
                rm -rf ${cwd}/KUL_LOG/sub-${participant}_anat_*.done >/dev/null 2>&1
                rm -rf $derivativesdir/KUL_anat_biascorrect >/dev/null 2>&1
                rm -rf $derivativesdir/KUL_anat_register_rigid >/dev/null 2>&1
                rm -rf $derivativesdir/KUL_anat_segment_tumor >/dev/null 2>&1
                rm -rf $globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz
            fi
        fi
        read -p "Redo: fmriprep? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_run_dwiprep.txt"
            rm -rf ${cwd}/fmriprep/sub-${participant} >/dev/null 2>&1
            rm -rf ${cwd}/fmriprep_work >/dev/null 2>&1
            rm -f ${cwd}/fmriprep/sub-${participant}.html >/dev/null 2>&1
        fi
        read -p "Redo: KUL_dwiprep? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_run_dwiprep.txt"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_run_dwiprep.txt >/dev/null 2>&1
            rm -rf ${cwd}/dwiprep/sub-${participant} >/dev/null 2>&1
            rm -rf $derivativesdir/synb0 >/dev/null 2>&1
        fi
        read -p "Redo: SPM? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_SPM.done"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_SPM.done >/dev/null 2>&1
            rm -f ${cwd}/RESULTS/sub-${participant}/Melodic/* >/dev/null 2>&1
            rm -fr $derivativesdir/SPM/* >/dev/null 2>&1
        fi
        read -p "Redo: Melodic? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_melodic.done"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_melodic.done >/dev/null 2>&1
            rm -f ${cwd}/RESULTS/sub-${participant}/SPM/* >/dev/null 2>&1
            rm -fr $derivativesdir/FSL_melodic/* >/dev/null 2>&1
        fi
        read -p "Redo: KUL_VBG? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_VBG.log"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_VBG.log >/dev/null 2>&1
            rm -fr $derivativesdir/KUL_VBG/* >/dev/null 2>&1
        fi
        read -p "Redo: msbp? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_MSBP.done"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_MSBP.done >/dev/null 2>&1
            rm -fr ${cwd}/BIDS/derivatives/cmp/sub-${participant} >/dev/null 2>&1
            rm -fr ${cwd}/BIDS/derivatives/nipype/sub-${participant} >/dev/null 2>&1
            rm -f ${cwd}/BIDS/derivatives/sub-${participant}_anatomical_config.ini >/dev/null 2>&1
        fi
        read -p "Redo: KUL_FWT? (y/n) " answ
        if [[ "$answ" == "y" ]]; then
            echo $answ
            echo "rm ${cwd}/KUL_LOG/sub-${participant}_FWT.done"
            rm -f ${cwd}/KUL_LOG/sub-${participant}_FWT.done >/dev/null 2>&1
            rm -fr $derivativesdir/KUL_FWT/* >/dev/null 2>&1
        fi
    fi
}

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
        if [ $d_flag -eq 1 ]; then
            KUL_dcm2bids.sh -d $dicomzip -p ${participant} -c study_config/sequences.txt -e -v
        else
            echo "Error: no dicom zip file given."
        fi
    else
        echo "BIDS conversion already done"
    fi
}

function KUL_check_data {
    
    mkdir -p $globalresultsdir
    echo -e "\n\nAn overview of the BIDS data:"
    bidsdir="BIDS/sub-$participant"
    T1w=($(find $bidsdir -name "*T1w.nii.gz" ! -name "*gadolinium*" -type f ))
    nT1w=${#T1w[@]}
    echo "  number of non-contrast T1w: $nT1w"

    cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
    ncT1w=${#cT1w[@]}
    # check if it has been moved away
    if [ $ncT1w -eq 0 ]; then
        if [ -f "$kulderivativesdir/sub-${participant}/cT1w/sub-${participant}_ce-gadolinium_T1w.nii.gz" ]; then
            ncT1w=-1
        fi
    fi
    echo "  number of contrast enhanced T1w: $ncT1w"


    FLAIR=($(find $bidsdir -name "*FLAIR.nii.gz" -type f ))
    nFLAIR=${#FLAIR[@]}
    echo "  number of FLAIR: $nFLAIR"
    FGATIR=($(find $bidsdir -name "*FGATIR.nii.gz" -type f ))
    nFGATIR=${#FGATIR[@]}
    echo "  number of FGATIR: $nFGATIR"
    T2w=($(find $bidsdir -name "*T2w.nii.gz" -type f ))
    nT2w=${#T2w[@]}
    echo "  number of T2w: $nT2w"
    SWI=($(find $bidsdir -name "*run-01_SWI.nii.gz" -type f ))
    nSWI=${#SWI[@]}
    SWIp=($(find $bidsdir -name "*run-02_SWI.nii.gz" -type f ))
    nSWIp=${#SWIp[@]}
    echo "  number of SWI magnitude: $nSWI"
    echo "  number of SWI phase: $nSWIp"

    # check the T1w
    if [ $nT1w -eq 0 ]; then
        echo "No T1w (without Gd) found. Fmriprep will not run."
        echo " Is the BIDS dataset correct?"
        read -p "Are you sure you want to continue? (y/n)? " answ
        if [[ "$answ" == "n" ]]; then
            exit 1
        fi
    fi 

    # check the cT1w
    if [ $ncT1w -eq -1 ]; then
        echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
        echo " At least one is missing."
        echo " The contrast T1w has been moved to the derivatives folder, due to previous processing."
        read -p "Do you want to restore it? (y/n)? " answ
        if [[ "$answ" == "y" ]]; then
            cp -f $kulderivativesdir/sub-${participant}/cT1w/sub-${participant}_ce-gadolinium_T1w.* $bidsdir/anat
            cT1w=($(find $bidsdir -name "*T1w.nii.gz" -name "*gadolinium*" -type f ))
            ncT1w=${#cT1w[@]}
        fi
    fi

    # check hd-glio-auto requirements
    if [ $hdglio -eq 1 ]; then
        if [ $nT1w -lt 1 ] || [ $ncT1w -lt 1 ] || [ $nT2w -lt 1 ] || [ $nT1w -lt 1 ]; then
            echo "For running hd-glio-auto a T1w, cT1w, T2w and FLAIR are required."
            echo " At least one is missing. Is the BIDS dataset correct?"
            read -p "Are you sure you want to continue? (y/n)? " answ
            if [[ "$answ" == "n" ]]; then
                exit 1
            fi
        fi
    fi 

    # check the BIDS
    find_fmri=($(find ${cwd}/BIDS/sub-${participant} -name "*_bold.nii.gz"))
    n_fMRI=${#find_fmri[@]}
    if [ $n_fMRI -eq 0 ]; then
        echo "WARNING: no fMRI data"
    fi

    find_dwi=($(find ${cwd}/BIDS/sub-${participant} -name "*_dwi.nii.gz"))
    n_dwi=${#find_dwi[@]}
    if [ $n_dwi -eq 0 ]; then
        echo "WARNING: no dwi data"
    fi
    echo -e "\n\n"

}

function KUL_rigid_register {
    warp_field="${registeroutputdir}/${source_mri_label}_reg2_T1w"
    output_mri="${globalresultsdir}/Anat/${source_mri_label}_reg2_T1w.nii.gz"
    #echo "Rigidly registering $source_mri to $target_mri"
    antsRegistration --verbose $ants_verbose --dimensionality 3 \
    --output [$warp_field,$output_mri] \
    --interpolation BSpline \
    --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform [$target_mri,$source_mri,1] \
    --transform Rigid[0.1] \
    --metric MI[$target_mri,$source_mri,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    #echo "Done rigidly registering $source_mri to $target_mri"
}

function KUL_run_fmriprep {
    if [ ! -f fmriprep/sub-${participant}.html ]; then
        
        # preparing for fmriprep
        cp study_config/run_fmriprep.txt KUL_LOG/sub-${participant}_run_fmriprep.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        if [ $n_fMRI -gt 0 ]; then
            fmriprep_options="--fs-no-reconall --use-aroma --use-syn-sdc "
        else
            fmriprep_options="--fs-no-reconall --anat-only "
        fi
        sed -i.bck "s/fmriprep_options: /fmriprep_options: ${fmriprep_options}/" KUL_LOG/sub-${participant}_run_fmriprep.txt
        rm -f KUL_LOG/sub-${participant}_run_fmriprep.txt.bck
        
        # running fmriprep
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_fmriprep.txt 
        
        # cleaning the working directory
        rm -fr fmriprep_work_${participant}
        
        # copying the result to the global results dir
        cp -f fmriprep/sub-$participant/anat/sub-${participant}_desc-preproc_T1w.nii.gz $globalresultsdir/Anat/T1w_fmriprep.nii.gz
        #gunzip -f $globalresultsdir/Anat/T1w.nii.gz
        
        # create a GM mask in the global results dir
        mrcalc fmriprep/sub-$participant/anat/sub-${participant}_dseg.nii.gz 1 -eq \
            fmriprep/sub-$participant/anat/sub-${participant}_dseg.nii.gz -mul - | \
            maskfilter - median - | \
            maskfilter - dilate $globalresultsdir/Anat/T1w_GM.nii.gz

        # add to the report
        if [ -f fmriprep/sub-${participant}.html ]; then
            ln -s ${cwd}/fmriprep/sub-${participant}.html ${cwd}/REPORT/sub-${participant}_03_fmriprep.html
        fi

    else
        echo "Fmriprep already done"
    fi
}

function KUL_run_dwiprep {
    if [ $n_dwi -gt 0 ];then
        if [ ! -f dwiprep/sub-${participant}/dwiprep_is_done.log ]; then
            cp study_config/run_dwiprep.txt KUL_LOG/sub-${participant}_run_dwiprep.txt
            sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_dwiprep.txt
            rm -f KUL_LOG/sub-${participant}_run_dwiprep.txt.bck
            
            KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_dwiprep.txt 
            
        else
            echo "Dwiprep already done"
        fi
    fi
}

function KUL_run_freesurfer {
    if [ ! -f BIDS/derivatives/freesurfer/${participant}_freesurfer_is.done ]; then
        cp study_config/run_freesurfer.txt KUL_LOG/sub-${participant}_run_freesurfer.txt
        sed -i.bck "s/BIDS_participants: /BIDS_participants: ${participant}/" KUL_LOG/sub-${participant}_run_freesurfer.txt
        rm -f KUL_LOG/sub-${participant}_run_freesurfer.txt.bck
        KUL_preproc_all.sh -e -c KUL_LOG/sub-${participant}_run_freesurfer.txt 
    else
        echo "Freesurfer already done"
    fi
}


function KUL_run_fastsurfer {
if [ ! -f KUL_LOG/sub-${participant}_FastSurfer.done ]; then
    #echo "Hybrid parcellation flag is set, now starting FastSurfer/FreeSurfer hybrid recon-all based part of VBG"

    # make your log file
    #prep_log="KUL_LOG/sub-${participant}_run_fastsurfer.txt" 
    #if [[ ! -f ${prep_log} ]] ; then
    #    touch ${prep_log}
    #else
    #    echo "${prep_log} already created"
    #fi
    kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"

    fs_output="${cwd}/BIDS/derivatives/freesurfer"
    #output_d="$derivativesdir/FastSurfer"
    #str_op="${output_d}/${participant}"
    #fasu_output="${str_op}fastsurfer"
    fasu_output="$derivativesdir/FastSurfer"
    #T1_4_parc="${str_op}_T1_nat_4parc.mgz"
    if [ $vbg -eq 1 ];then
        T1_4_parc=$derivativesdir/KUL_VBG/output_VBG/sub-${participant}/sub-Casier_T1_nat_4parc.mgz
    else
        T1_4_parc="${cwd}/$T1w"
    fi
    #echo $T1_4_parc
    #ls -l $T1_4_parc

    #recall_scripts="${fs_output}/sub-${participant}/scripts"
    #echo $recall_scripts
 
    mkdir -p ${fs_output} >/dev/null 2>&1
    mkdir -p ${fasu_output} >/dev/null 2>&1

    # Run recon-all and convert the T1 to .mgz for display
    # running with -noskulltrip and using brain only inputs
    # for recon-all
    # if we can run up to skull strip, break, fix with hd-bet result then continue it would be much better
    # if we can switch to fast-surf, would be great also
    # another possiblity is using recon-all -skullstrip -clean-bm -gcut -subjid <subject name>
    
    #echo "starting recon-all stage 1"
    task_in="recon-all -i ${T1_4_parc} -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -autorecon1 -no-isrunning"
    KUL_task_exec $verbose_level "FastSurfer part 1: recon-all stage 1" "$KUL_LOG_DIR/FastSurfer"
    #echo "done recon-all stage 1"



    #task_in="mri_convert -rl ${fs_output}/${participant}/mri/brainmask.mgz ${T1_BM_4_FS} ${clean_BM_mgz}"
    #task_exec

    #task_in="mri_mask ${FS_brain} ${T1_BM_4_FS} ${new_brain} && mv ${new_brain} ${fs_output}/${participant}/mri/brainmask.mgz && cp \
    #${fs_output}/${participant}/mri/brainmask.mgz ${fs_output}/${participant}/mri/brainmask.auto.mgz"
    #task_exec


    FaSu_loc=$(which run_fastsurfer.sh)
    #nvd_cu=$(nvcc --version)
    user_id_str=$(id -u $(whoami))
    T1_4_FaSu=$(basename ${T1_4_parc})
    nvram=$(echo $(nvidia-smi --query-gpu=memory.free --format=csv) | rev | cut -d " " -f2 | rev)
    if [[ ! -z ${nvram} ]]; then
        if [[ ${nvram} -lt 6000 ]]; then
            batch_fasu="4"
        elif [[ ${nvram} -gt 6500 ]] && [[ ${nvram} -lt 7000 ]]; then
            batch_fasu="6"
        elif [[ ${nvram} -gt 7000 ]]; then
            batch_fasu="8"
        fi
    else
        batch_fasu="2"
    fi


    if [[ ! -z ${FaSu_loc} ]]; then

        if [ ${nvram} -lt 5500 ]; then
            FaSu_cpu=" --no_cuda "
            FaSu_mode="cpu"
            #echo " Running FastSurfer without CUDA " | tee -a ${prep_log}
        else
            FaSu_cpu=""
            FaSu_mode="cuda-gpu"
            #echo " Running FastSurfer with CUDA " | tee -a ${prep_log}
        fi

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        task_in="run_fastsurfer.sh --t1 ${T1_4_parc} \
        --sid sub-${participant} --sd ${fasu_output} --fsaparc --parallel --threads ${ncpu} \
        --fs_license $FS_LICENSE --py python ${FaSu_cpu} --ignore_fs_version --batch ${batch_fasu}"
        kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"
        KUL_task_exec $verbose_level "FastSurfer part 2: Fastsurfer itself (script & $FaSu_mode mode)" "$KUL_LOG_DIR/FastSurfer"

    else

        # it's a good idea to run autorecon1 first anyway
        # then use the orig from that to feed to FaSu

        echo "Local FastSurfer not found, switching to Docker version" | tee -a ${prep_log}
        T1_4_FaSu=$(basename ${T1_4_parc})
        dir_4_FaSu=$(dirname ${T1_4_parc})

        if [ ${nvram} -lt 5500 ]; then
            FaSu_v="cpu"
        else
            FaSu_v="gpu"
        fi

        task_in="docker run -v ${dir_4_FaSu}:/data -v ${fasu_output}:/output \
        -v $FREESURFER_HOME:/fs60 --rm --user ${user_id_str} fastsurfer:${FaSu_v} \
        --fs_license /fs60/$(basename $FS_LICENSE) --sid sub-${participant} \
        --sd /output/ --t1 /data/${T1_4_FaSu} \
        --parallel --threads ${ncpu}"
        kul_log_file="KUL_LOG/sub-${participant}_run_fastsurfer.txt"
        KUL_task_exec $verbose_level "FastSurfer part 2: Fastsurfer itself (docker & $FaSu_mode mode)" "$KUL_LOG_DIR/FastSurfer"

    fi

    #fs_output="${cwd}/BIDS/derivatives/freesurfer"
    #fasu_output="$derivativesdir/FastSurfer"
    

    # time to copy the surfaces and labels from FaSu to FS dirtask_exec
    # here we run FastSurfer first and 

    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/surf ${output_d}/${participant}_FS_output/${participant}/
    #cp -rf ${output_d}/${participant}fastsurfer/${participant}/label ${output_d}/${participant}_FS_output/${participant}/
    
    #cp -rf ${output_d}/sub-${participant}/surf/* $fs_output/sub-${participant}/surf/
    #cp -rf ${output_d}/sub-${participant}/label/* $fs_output/sub-${participant}/surf/label/ 

    rsync -azv ${fasu_output}/sub-${participant}/ ${fs_output}/sub-${participant}/

    #task_in="recon-all -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} -parallel -all -noskullstrip"
    #task_exec

    task_in="recon-all -s sub-${participant} -sd ${fs_output} -openmp ${ncpu} \
        -parallel -no-isrunning -make all"
    KUL_task_exec $verbose_level "FastSurfer part 3: recon-all -make-all" "FastSurfer"

    #exit

    #task_in="mri_convert -rl ${fs_output}/${participant}/mri/brain.mgz ${T1_brain_clean} ${fs_output}/${participant}/mri/real_T1.mgz"

    #task_exec

    #task_in="mri_convert -rl ${fs_output}/${participant}/mri/brain.mgz -rt nearest ${Lmask_o} ${fs_output}/${participant}/mri/Lmask_T1_bin.mgz"

    #task_exec

    #fs_parc_mgz="${fs_output}/${participant}/mri/aparc+aseg.mgz"
    touch KUL_LOG/sub-${participant}_FastSurfer.done

else
    echo "Already done FastSurfer"
fi
}


function KUL_segment_tumor {
    
    # Segmentation of the tumor
    
    # check if it needs to be performed
    if [ $hdglio -eq 1 ];then

        if [ ! -f "$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" ]; then
            KUL_anat_segment_tumor.sh -p $participant -v $verbose_level
            lesion_png=RESULTS/sub-${participant}/Lesion/sub-${participant}_tumor_segment.png
            if [ -f $lesion_png ]; then
                cp -f $lesion_png REPORT/sub-${participant}_01_tumor_segment.png
            fi
        else
            echo "Tumor segmentation already done"
        fi

    fi

}

function KUL_run_VBG {
    if [ $vbg -gt 0 ]; then
        if [ $vbg -eq 1 ]; then
            vbg_extra_axial=""
            vbg_lesion="$derivativesdir/KUL_anat_segment_tumor/sub-${participant}_lesion_and_cavity.nii.gz"
        elif [ $vbg -eq 2 ]; then
            vbg_extra_axial="-E"
            vbg_lesion="$derivativesdir/KUL_anat_segment_tumor/sub-${participant}_lesion_and_cavity.nii.gz"
        elif [ $vbg -eq 3 ]; then
            vbg_extra_axial=""
            vbg_lesion="${cwd}/RESULTS/sub-${participant}/Lesion/lesion.nii.gz"
        fi

        vbg_test="$derivativesdir/KUL_VBG/output_VBG/sub-${participant}/sub-${participant}_T1_nat_filled.nii.gz"
        if [[ ! -f $vbg_test ]]; then
            echo "Computing KUL_VBG"
            mkdir -p ${cwd}/BIDS/derivatives/freesurfer/sub-${participant}
            mkdir -p $derivativesdir/KUL_VBG
            
            # See to it that freesurfer 6 is used
            export FREESURFER_HOME=/usr/local/KUL_apps/freesurfer_6.0.0
            export SUBJECTS_DIR=$FREESURFER_HOME/subjects
            export FS_LICENSE=$FREESURFER_HOME/license.txt
            source $FREESURFER_HOME/SetUpFreeSurfer.sh
            fs_v=$(recon-all --version)
            kul_echo "Using $fs_v for VBG"

            task_in="KUL_VBG.sh -S ${participant} \
                -l $vbg_lesion \
                -o $derivativesdir/KUL_VBG \
                -m $derivativesdir/KUL_VBG \
                $vbg_extra_axial \
                -z T1 -b -B 1 -t -P 1 -n $ncpu"
            KUL_task_exec $verbose_level "KUL_VBG" "7_VBG"

            # copy the output of VBG to the derivatives freesurfer directory
            cp -r $derivativesdir/KUL_VBG/output_VBG/sub-${participant}/sub-${participant}_FS_output/sub-${participant} \
                BIDS/derivatives/freesurfer/

            # add to the report
            KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz \
                -t 2 -d REPORT -f 04_VBG_input
            KUL_mrview_figure.sh -p ${participant} -u BIDS/derivatives/KUL_compute/sub-${participant}/KUL_VBG/output_VBG/sub-${participant}/sub-${participant}_T1_nat_filled.nii.gz \
                -t 2 -d REPORT -f 04_VBG_output
            #convert label:"Results of VBG:\nTop: Original T1w\nBottom: Filled T1w" REPORT/VBG_temp_caption.png
            montage REPORT/sub-${participant}_04_VBG_input.png REPORT/sub-${participant}_04_VBG_output.png -tile 1x2 -geometry +0+0 \
                -mode Concatenate REPORT/sub-${participant}_04_VBG.png   
            rm REPORT/sub-${participant}_04_VBG_input.png REPORT/sub-${participant}_04_VBG_output.png  

            echo "Done computing KUL_VBG"

        else
            echo "KUL_VBG has already run"
        fi

    fi
}

function KUL_run_msbp {

    if [ ! -f KUL_LOG/sub-${participant}_MSBP.done ]; then

        echo "Running MSBP"

        # there seems tpo be a problem with docker if the fsaverage dir is a soft link; so we delete the link and hardcopy it
        rm -fr $cwd/BIDS/derivatives/freesurfer/fsaverage
        cp -r $FREESURFER_HOME/subjects/fsaverage $cwd/BIDS/derivatives/freesurfer/fsaverage

        task_in="docker run --rm -u $(id -u) -v $cwd/BIDS:/bids_dir \
         -v $cwd/BIDS/derivatives:/output_dir \
         -v $FS_LICENSE:/opt/freesurfer/license.txt \
         sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
         --participant_label $participant --isotropic_resolution 1.0 --thalamic_nuclei \
         --brainstem_structures --skip_bids_validator --fs_number_of_cores $ncpu \
         --multiproc_number_of_cores $ncpu"
        KUL_task_exec $verbose_level "MSBP" "10_msbp"

        echo "Done MSBP"
        touch KUL_LOG/sub-${participant}_MSBP.done
        
    else
        echo "MSBP already done"
    fi
}

function KUL_run_FWT {
    if [ $n_dwi -gt 0 ];then
        config="tracks_list.txt"
        if [ ! -f KUL_LOG/sub-${participant}_FWT.done ]; then

            task_in="KUL_FWT_make_VOIs.sh -p ${participant} \
            -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
            -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
            -c $cwd/study_config/${config} \
            -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
            -o $kulderivativesdir/sub-${participant}/FWT \
            -n $ncpu"
            KUL_task_exec $verbose_level "KUL_FWT voi generation" "12_FWTvoi"


            task_in="KUL_FWT_make_TCKs.sh -p ${participant} \
            -F $cwd/BIDS/derivatives/freesurfer/sub-${participant}/mri/aparc+aseg.mgz \
            -M $cwd/BIDS/derivatives/cmp/sub-${participant}/anat/sub-${participant}_label-L2018_desc-scale3_atlas.nii.gz \
            -c $cwd/study_config/${config} \
            -d $cwd/dwiprep/sub-${participant}/sub-${participant} \
            -o $kulderivativesdir/sub-${participant}/FWT \
            -T 1 -a iFOD2 \
            -f 1 \
            -Q -S \
            -n $ncpu"
            KUL_task_exec $verbose_level "KUL_FWT tract generation" "12_FWTtck"

            rm -fr $globalresultsdir/Tracto/*
            mcp -o "$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*_fin_map_BT_iFOD2.nii.gz" \
                "$globalresultsdir/Tracto/Tract-csd_#2.nii.gz"
            mcp -o "$kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*/*_fin_BT_iFOD2.tck" \
                "$globalresultsdir/Tracto/Tract-csd_#2.tck"
            pdfunite $kulderivativesdir/sub-${participant}/FWT/sub-${participant}_TCKs_output/*_output/Screenshots/*fin_BT_iFOD2_inMNI_screenshot2_niGB.pdf $globalresultsdir/Tracto/Tracts_Summary.pdf
            
            # add to report
            cp -f RESULTS/sub-${participant}/Tracto/Tracts_Summary.pdf REPORT/sub-${participant}_06_Tract_Summary.pdf

            touch KUL_LOG/sub-${participant}_FWT.done
            
        else
            echo "FWT already done"
        fi
    fi
}

function KUL_anatomical_biascorrect {
    check="KUL_LOG/sub-${participant}_anat_biascorrect.done"
    if [ ! -f $check ]; then

        KUL_anat_biascorrect.sh -p $participant -v $verbose_level

        touch $check
    else 
        echo "Anatomical bias correction already done"
    fi
}

function KUL_register_anatomical_images {
    check="KUL_LOG/sub-${participant}_anat_reg.done"
    if [ ! -f $check ]; then

        KUL_anat_register.sh -p $participant -c -v $verbose_level
        cp  $cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_anat_register_rigid/*reg2_T1w.nii.gz $globalresultsdir/Anat/
        cp  $cwd/BIDS/derivatives/KUL_compute/sub-${participant}/KUL_anat_register_rigid/T1w.nii.gz $globalresultsdir/Anat/
        touch $check
    
    else 
        echo "Anatomical registration already done"
    fi
}

function KUL_clear_cT1w {
    
    # a funtion to remove the cT1w (gadolinium enhanced T1w) away since it conflicts during msbp
    clear_cT1w_outputdir="$kulderivativesdir/sub-${participant}/cT1w"
    mkdir -p $clear_cT1w_outputdir

    if [ $ncT1w -gt 0 ]; then
        source_mri="${cT1w%*.nii.gz}*"
        if [ -f $cT1w ]; then
            mv $source_mri $clear_cT1w_outputdir 
        fi
    fi

}

function KUL_fmriproc {

    if [ $n_fMRI -gt 0 ];then

        if [ ! -f ${cwd}/KUL_LOG/sub-${participant}_SPM.done ]; then
            task_in="KUL_fmriproc_spm.sh -p $participant"
            KUL_task_exec $verbose_level "KUL_fmriproc_spm" "7_fmriproc_spm"

            # add to report
            for spm in RESULTS/sub-${participant}/SPM/*.nii; do
                #echo $spm
                max_T=$(mrstats -output max $spm)
                #echo $max_T
                thresh=$(awk "BEGIN {print $max_T/3}")
                task=$(basename $spm)
                mrcalc $spm $thresh -gt REPORT/spm_tmp_$task
                KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz -o REPORT/spm_tmp_$task \
                    -t 2 -d REPORT -f 05_afMRI_${task}_Thr_${thresh}
                rm -f REPORT/spm_tmp_$task
            done
        fi

        if [ ! -f ${cwd}/KUL_LOG/sub-${participant}_melodic.done ]; then
            task_in="KUL_fmriproc_conn.sh -p $participant"
            KUL_task_exec $verbose_level "KUL_fmriproc_conn" "8_fmriproc_conn"

            # add to report
            for spm in RESULTS/sub-${participant}/Melodic/*.nii; do
                #echo $spm
                max_T=$(mrstats -output max $spm)
                #echo $max_T
                thresh=$(awk "BEGIN {print $max_T/3}")
                task=$(basename $spm)
                mrcalc $spm $thresh -gt REPORT/spm_tmp_$task
                KUL_mrview_figure.sh -p ${participant} -u RESULTS/sub-${participant}/Anat/T1w.nii.gz -o REPORT/spm_tmp_$task \
                    -t 2 -d REPORT -f 05_rsfMRI_${task}_Thr_${thresh}
                rm -f REPORT/spm_tmp_$task
            done
        fi
    fi

}

function KUL_run_dwiprep_anat {

    dwi_anat_check=${cwd}/KUL_LOG/sub-${participant}_dwiprep_anat.done
    if [ ! -f $dwi_anat_check ]; then
        task_in="KUL_dwiprep_anat.sh -p $participant -n $ncpu"
        KUL_task_exec $verbose_level "KUL_dwiprep_anat" "11_dwiprep_anat"

        if [ -f dwiprep/sub-${participant}/sub-${participant}/qa/sub-${participant}_T1w_with_fa.png ]; then
            cp -f dwiprep/sub-${participant}/sub-${participant}/qa/sub-${participant}_T1w_with_fa.png \
                REPORT/sub-${participant}_02_T1w_with_fa.png
            cp -f dwiprep/sub-${participant}/sub-${participant}/qa/sub-${participant}_T1w_brain_with_fa.png \
                REPORT/sub-${participant}_02_T1w_brain_with_fa.png
            cp -f dwiprep/sub-${participant}/sub-${participant}/eddy_qc/quad/qc.pdf REPORT/sub-${participant}_02_eddy_qc.pdf
        fi

        touch $dwi_anat_check
    fi

}

function KUL_run_cT1w_subtraction {

    cT1w_subt_check=${cwd}/KUL_LOG/sub-${participant}_cT1w_subtraction.done
    if [ ! -f $cT1w_subt_check ]; then
        if [ $nT1w -gt 0 ] && [ $ncT1w -gt 0 ]; then
            echo "Creating T1w subtraction image"
            mkdir -p $derivativesdir/cT1w_subtraction
            # get the mask
            cp -f fmriprep/sub-$participant/anat/sub-${participant}_desc-brain_mask.nii.gz $derivativesdir/cT1w_subtraction/brain_mask.nii.gz
            # get the lesion & exclude the tumor if any
            #echo $hdglio
            if [ $hdglio -eq 1 ];then
                if [  -f "$globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz" ]; then
                    mrcalc $derivativesdir/cT1w_subtraction/brain_mask.nii.gz $globalresultsdir/Lesion/sub-${participant}_lesion_and_cavity.nii.gz \
                        -sub $derivativesdir/cT1w_subtraction/brain_final_mask.nii.gz -force
                    tomatch=$derivativesdir/cT1w_subtraction/brain_final_mask.nii.gz
                else
                    tomatch=$derivativesdir/cT1w_subtraction/brain_mask.nii.gz
                fi
            else
                tomatch=$derivativesdir/cT1w_subtraction/brain_mask.nii.gz
            fi
            # historgam match the cT1w and T1w
            mrhistmatch linear \
                -mask_input $tomatch \
                -mask_target $tomatch \
                $globalresultsdir/Anat/T1w.nii.gz $globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz \
                $globalresultsdir/Anat/T1w_matched2_cT1w.nii.gz -force
            # subtract them    
            mrcalc $globalresultsdir/Anat/cT1w_reg2_T1w.nii.gz $globalresultsdir/Anat/T1w_matched2_cT1w.nii.gz -sub \
                $globalresultsdir/Anat/cT1w_T1w_subtracted.nii.gz -force
            touch $cT1w_subt_check
        fi
    else
        echo "T1w subtraction already done"    
    fi

}


# --- MAIN ---
# Prepare for sending to pacs.
cp -f $kul_main_dir/share/KUL_mevislab_* .

# Check if scaffolding has happend
if [ ! -d $cwd/study_config ]; then
    KUL_scaffold
fi


# STEP 1 - BIDS conversion
KUL_convert2bids

KUL_check_participant

kulderivativesdir=$cwd/BIDS/derivatives/KUL_compute
mkdir -p $kulderivativesdir
mkdir -p $globalresultsdir/Anat
mkdir -p $globalresultsdir/SPM
mkdir -p $globalresultsdir/Melodic
mkdir -p $globalresultsdir/Tracto
mkdir -p $globalresultsdir/PACS/fMRI
mkdir -p $cwd/RESULTS
mkdir -p $cwd/REPORT

if [ $KUL_DEBUG -gt 0 ]; then 
    echo "kulderivativesdir: $kulderivativesdir"
    echo "globalresultsdir: $globalresultsdir"
fi

# Run BIDS validation
check_in=${KUL_LOG_DIR}/1_bidscheck.done
if [ ! -f $check_in ]; then

    docker run -ti --rm -v ${cwd}/BIDS:/data:ro bids/validator /data

    read -p "Are you happy? (y/n) " answ
    if [[ ! "$answ" == "y" ]]; then
        exit 1
    else
        touch $check_in
    fi
fi


# Check if fMRI and/or dwi data are present and/or to redo some processing
echo "Starting KUL_clinical_fmridti"
KUL_check_data
KUL_check_redo


# STEP 1 - run bias correction
KUL_anatomical_biascorrect


# STEP 2 - register all anatomical other data to the T1w without contrast
KUL_register_anatomical_images


# STEP 3 - run tumor segmentation 
KUL_segment_tumor
    

# STEP 4 - run fmriprep and continue
KUL_run_fmriprep &


# STEP 5 - run dwiprep and continue
KUL_run_dwiprep &


wait

# STEP 6 - generate Gd contrast T1w subtraction
KUL_run_cT1w_subtraction


# STEP 7 - get rid of the Gadolinium T1w image since it may conflict with msbp
KUL_clear_cT1w


# STEP 8 - run SPM & melodic
KUL_fmriproc


# STEP 9 - run VBG
KUL_run_VBG 
wait


# STEP 8 - dev
#KUL_run_fastsurfer
#KUL_run_freesurfer # let msbp also do FS
#wait 


# STEP 10 - run SPM/melodic/msbp
KUL_run_msbp
wait



# STEP 11 run dwiprep_anat
KUL_run_dwiprep_anat


# STEP 12 - run Fun With Tracts
KUL_run_FWT



# STEP 13 - call yourself to make tractography figures
fig_check=${cwd}/KUL_LOG/sub-${participant}_figures.done
if [ ! -f $fig_check ]; then
    #echo $type
    #echo $ncT1w
    if [ $ncT1w -gt 0 ]  || [ $ncT1w -eq -1 ]; then 
        KUL_clinical_fmridti_new.sh -p $participant -t $type -R 1 
    fi
    #echo $nFLAIR
    if [ $nFLAIR -gt 0 ]; then 
        KUL_clinical_fmridti_new.sh -p $participant -t $type -R 2
    fi
    if [ $nSWI -gt 0 ]; then 
        KUL_clinical_fmridti_new.sh -p $participant -t $type -R 3
    fi
    if [ $nT1w -gt 0 ] && [ $ncT1w -lt 1 ] && [ $nFLAIR -eq 0 ]; then 
        KUL_clinical_fmridti_new.sh -p $participant -t $type -R 4
    fi
    touch $fig_check
else 
    echo "Figures already done"
fi



# STEP 14 - run Karawun
karawun_check=${cwd}/KUL_LOG/sub-${participant}_karawun.done
if [ ! -f $karawun_check ]; then
    if [ $type -lt 5 ]; then
        kul_echo "Running KUL_karawun for a tumor patient"
        KUL_karawun_prepare.sh -p ${participant} -t 1 -r 3
    elif [ $type -eq 5 ]; then
        kul_echo "Running KUL_karawun for a DBS ET patient"
        KUL_karawun_prepare.sh -p ${participant} -t 2 -r 10 
    elif [ $type -eq 6 ]; then
        kul_echo "Running KUL_karawun for a DBS Parkinson patient"
        KUL_karawun_prepare.sh -p ${participant} -t 3 -r 10 
    fi
    touch $karawun_check
fi

echo "Finished"
