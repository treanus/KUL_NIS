#!/bin/bash
# Bash shell script to visualise (f)MRI/dMRI results
#
#
# @ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be
# 11/10/2022
version="0.1"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")
source $kul_main_dir/KUL_main_functions.sh
# $cwd & $log_dir is made in main_functions

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` visualises (f)MRI/dMRI results

Usage:

  `basename $0` <OPT_ARGS>

Example:

  `basename $0` -p JohnDoe -t 1

Required arguments:

     -p:  participant name

Optional arguments:

     -u:  underlay
            specify a file 
        or use (from RESULTS/Anat directory)
            1: use cT1w as underlay
            2: use FLAIR as underlay
            3: use SWI as underlay
            4: Use T1w as underlay (default)
     -o:  overlay
            specify 1 file to use as overlay
        or
            specify a text file with many overlays and settings
     -d:  output directory (default=RESULTS/sub-participant/View)
     -f:  output_filename (default=underlay_with_overlay)
     -t:  view type
        0: or not given = open mrview to view
        1: produce TRA/SAG/COR png files of every slice
        2: produce a png montage with TRA sections of every couple of slices
     -v:  show output from commands (0=silent, 1=normal, 2=verbose; default=1)

USAGE

	exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
#
# Set defaults
silent=1 # default if option -v is not given
underlay=""
overay=""
view=0
output_dir=""
d_given=0
output_file=""
verbose_level=1

# Set required options
p_flag=0


if [ "$#" -lt 1 ]; then
	Usage >&2
	exit 1

else

	while getopts "p:u:o:d:f:t:v:" OPT; do

		case $OPT in
		p) #participant
			participant=$OPTARG
            p_flag=1
		;;
        u) #underlay
			underlay=$OPTARG
		;;
        o) #overlay
			overlay=$OPTARG
		;;
        t) #view type
			view=$OPTARG
		;;
		d) #output_dir
			output_dir=$OPTARG
            d_given=1
		;;
        f) #output_file
			output_file=$OPTARG
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

function KUL_mrview {

    cmd="mrview \
        $mrview_underlay \
        $mrview_mode \
        $mrview_plane \
        $mrview_overlay \
        $mrview_annotations \
        $mrview_capture \
        $mrview_exit"
    #echo $cmd
    eval $cmd

}

#### MAIN

# Check the underlay or set default
echo "underlay: $underlay"
if [ -z "$underlay" ];then
    mrview_underlay="$cwd/RESULTS/sub-$participant/Anat/T1w.nii.gz"
else
    if [ -n "$underlay" ] && [ "$underlay" -eq "$underlay" ] 2>/dev/null; then
        if [ $underlay -eq 4 ]; then 
            mrview_underlay="$cwd/RESULTS/sub-$participant/Anat/T1w.nii.gz"
        elif [ $underlay -eq 1 ]; then 
            mrview_underlay="$cwd/RESULTS/sub-$participant/Anat/cT1w_reg2_T1w.nii.gz"
        elif [ $underlay -eq 2 ]; then 
            mrview_underlay="$cwd/RESULTS/sub-$participant/Anat/FLAIR_reg2_T1w.nii.gz"
        elif [ $underlay -eq 3 ]; then 
            mrview_underlay="$cwd/RESULTS/sub-$participant/Anat/SWI_reg2_T1w.nii.gz"
        fi
    else
        mrview_underlay="$underlay"
    fi
fi

# Check once more that the underlay exists on disk
if [ ! -f "$underlay" ];then
    echo "The underlay does not exist. Sorry quitting."
    exit 1
fi

# Check the output dir or set default
if [ -z "$output_dir" ]; then 

    underlay_name_tmp=$(basename $mrview_underlay)
    #echo "underlay_name_tmp: $underlay_name_tmp"
    underlay_name=${underlay_name_tmp%%.nii.gz}
    output_dir=$cwd/RESULTS/sub-$participant/View/$underlay_name

fi
#echo "output_dir: $output_dir"


# Check overlay 
if [ -z "$overlay" ]; then

    mrview_overlay=""
    base_overlay=""

elif [[ "$overlay" == *".txt" ]]; then

    #echo "found txt file"
    base_overlay="$overlay"
    mrview_overlay=$(cat $overlay)
    
else

    mrview_overlay="-overlay.load $overlay -overlay.opacity 0.4"
    base_overlay_tmp=$(basename $overlay)
    base_overlay=${base_overlay_tmp%%.nii.gz}
    
fi
#echo "mrview_overlay: $mrview_overlay"
#echo "base_overlay: $base_overlay"


# Check out file name
if [ -z $output_file ]; then
    final_output_dir="${output_dir}_with_${base_overlay}"
    final_output_file=$(basename $final_output_dir)
else
    # strip .png if any
    output_file=${output_file%%.png}
    if [ -z $base_overlay ]; then
        final_output_dir="${output_dir}"
    elif [ $d_given -eq 1 ]; then
        final_output_dir="${output_dir}"
    else
        final_output_dir="${output_dir}_with_${base_overlay}"
    fi
    final_output_file="sub-${participant}_$output_file"
fi
mkdir -p $final_output_dir
#echo "Setting final_output_dir: $final_output_dir"
#echo "Final output filename: $final_output_file"

# Check the view mode or set default
#echo "view_type: $view"
if [ $view -eq 0 ]; then

    mrview_mode="-mode 2"
    mrview_exit="&"
    mrview_capture=""
    mrview_plane=""
    mrview_annotations=""
    KUL_mrview

elif [ $view -eq 1 ]; then

    mrview_mode="-mode 1"
    mrview_exit="-exit"
    mrview_annotations="-noannotations"

    ori[0]="TRA"
    ori[1]="SAG"
    ori[2]="COR"

    for orient in ${ori[@]}; do

        if [[ "$orient" == "TRA" ]]; then
            underlay_slices=$(mrinfo $mrview_underlay -size | awk '{print $(NF)}')
        elif [[ "$orient" == "SAG" ]]; then
            underlay_slices=$(mrinfo $mrview_underlay -size | awk '{print $(NF-2)}')
        else
            underlay_slices=$(mrinfo $mrview_underlay -size | awk '{print $(NF-1)}')
        fi
    
        i=0
        #echo ${orient}
        mkdir -p $final_output_dir/${orient}
        voxel_index=""
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
        mrview_capture="-capture.folder $final_output_dir/${orient} -capture.prefix tmp $voxel_index"
        mrview_plane="-plane $plane"
        KUL_mrview

    done


elif [ $view -eq 2 ]; then

    mrview_mode="-mode 1"
    mrview_exit="-exit"
    mrview_plane="-plane 2"
    mrview_annotations="-noannotations"

    underlay_slices=$(mrinfo $mrview_underlay -size | awk '{print $(NF)}')
    #echo "underlay_slices: $underlay_slices"
    voxel_index=""
    i=0
    while [ $i -lt $underlay_slices ]
    do
        #echo Number: $i
        voxel_index="$voxel_index -voxel 0,0,$i -capture.grab"
        let "i+=7" 
    done

    mrview_capture="-capture.folder $final_output_dir -capture.prefix tmp $voxel_index"
    KUL_mrview

    montage $final_output_dir/tmp*.png \
        -mode Concatenate \
        $final_output_dir/${final_output_file}.png
    rm -f $final_output_dir/tmp*.png

fi



