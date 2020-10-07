    #!/bin/bash
# @ Stefan Sunaert & Ahmed Radwan- UZ/KUL - stefan.sunaert@uzleuven.be
#
v="v0.1 - dd 07/10/2020"

# This is the main script of the KUL_NeuroImaging_Toools
#


pbs_cpu=$(grep pbs_cpu $conf | grep -v \# | sed 's/[^0-9]//g')
        pbs_mem=$(grep pbs_mem $conf | grep -v \# | sed 's/[^0-9]//g')
        pbs_lp=$(grep pbs_lp $conf | grep -v \# | cut -d':' -f 2 | tr -d '\r')
        pbs_email=$(grep pbs_email $conf | grep -v \# | cut -d':' -f 2 | tr -d '\r')
        pbs_walltime=$(grep pbs_walltime $conf | grep -v \# | cut -d':' -f 2- | tr -d '\r')
        pbs_singularity_mriqc=$(grep pbs_singularity_mriqc $conf | grep -v \# | cut -d':' -f 2 | tr -d '\r')
        pbs_singularity_fmriprep=$(grep pbs_singularity_fmriprep $conf | grep -v \# | cut -d':' -f 2 | tr -d '\r')

        if [ $silent -eq 0 ]; then

            echo "  pbs_cpu: $pbs_cpu"
            echo "  pbs_mem: $pbs_mem"
            echo "  pbs_lp: $pbs_lp"
            echo "  pbs_email: $pbs_email"
            echo "  pbs_walltime: $pbs_walltime"
            echo "  pbs_singularity_mriqc: ${pbs_singularity_mriqc}"
            echo "  pbs_singularity_fmriprep: $pbs_singularity_fmriprep"

        fi
cp $kul_main_dir/VSC/master_dwiprep.pbs VSC/run_mriqc.pbs



    perl  -pi -e "s/##LP##/${pbs_lp}/g" VSC/run_mriqc.pbs
    perl  -pi -e "s/##CPU##/${pbs_cpu}/g" VSC/run_mriqc.pbs
    perl  -pi -e "s/##GPU##/${pbs_mem}/g" VSC/run_mriqc.pbs
    esc_pbs_email=$(echo $pbs_email | sed 's#\([]\!\(\)\#\%\@\*\$\/&\-\=[]\)#\\\1#g')
    perl  -pi -e "s/##EMAIL##/${esc_pbs_email}/g" VSC/run_mriqc.pbs
    esc_pbs_walltime=$(echo $pbs_walltime | sed 's#\([]\!\(\)\#\%\@\*\$\/&\-\=[]\)#\\\1#g')
    perl  -pi -e "s/##WALLTIME##/${esc_pbs_walltime}/g" VSC/run_mriqc.pbs
    
    #esc_pbs_singularity_fmriprep=$(echo $pbs_singularity_fmriprep | sed 's#\([]\!\(\)\#\%\@\*\$\/&\-\=[]\)#\\\1#g')
    #perl  -pi -e "s/##FMRIPREP##/${esc_pbs_singularity_fmriprep}/g" VSC/run_mriqc.pbs
    #esc_pbs_singularity_mriqc=$(echo $pbs_singularity_mriqc | sed 's#\([]\!\(\)\#\%\@\*\$\/&\-\=[]\)#\\\1#g')
    #perl  -pi -e "s/##MRIQC##/${esc_pbs_singularity_mriqc}/g" VSC/run_mriqc.pbs
    
    esc_task_command=$(echo $task_command | sed 's#\([]\!\(\)\#\%\@\*\$\/&\-\=[]\)#\\\1#g')
    perl  -pi -e "s/##COMMAND##/${esc_task_command}/g" VSC/run_mriqc.pbs

