#!/bin/bash

dicom="DICOM/Croes_Rudi_Joseph/"
participant="Croes"


# convert the DICOM to BIDS
KUL_dcm2bids.sh -d $dicom -p $participant -c study_config/sequences.txt -e

# run fmriprep
cp study_config/run_fmriprep.txt KUL_LOG/run_fmriprep_$participant.txt
sed -i "s/BIDS_participants: /BIDS_participants: $participant/" KUL_LOG/run_fmriprep_$participant.txt
KUL_preproc_all.sh -e -c KUL_LOG/run_fmriprep_Croes.txt 

# run SPM12
tcf="/DATA/test/study_config/test.m" #template config file
tjf="/DATA/test/study_config/test_job.m" #template job file
pcf="/DATA/test/study_config/test_$participant.m" #participant config file
pjf="/DATA/test/study_config/test_job_$participant.m" #participant job file
cp $tcf $pcf
cp $tjf $pjf
sed -i "s|###JOBFILE###|$pjf|" $pcf

fmridir="/DATA/test/SPM"
fmrifile="HAND_bold.nii"
fmriresults="/DATA/test/SPM/RESULTS"
sed -i "s|###FMRIDIR###|$fmridir|" $pjf
sed -i "s|###FMRIFILE###|$fmrifile|" $pjf
sed -i "s|###FMRIRESULTS###|$fmriresults|" $pjf

/usr/local/MATLAB/R2018a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('$pcf');exit ; "

