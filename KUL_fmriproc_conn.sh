mkdir -p CONN/input
cd CONN/input

ln -s  ../../fmriprep/sub-*/*/func/sub-*_ses-*_task-rest_run-02_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz .
ln -s  ../../fmriprep/sub-*/*/func/sub-*_ses-*task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz .
ln -s ../../fmriprep/sub-*/anat/sub-*_space-MNI152NLin6Asym_desc-preproc_T1w.nii.gz .
ln -s ../../fmriprep/sub-*/anat/sub-*_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz .
ln -s ../../fmriprep/sub-*/anat/sub-*_space-MNI152NLin6Asym_label*.nii.gz .
ln -s ../../fmriprep/sub-*/anat/sub-*_space-MNI152NLin2009cAsym_label*.nii.gz .

cp  fmriprep/sub-*/*/func/sub-*_ses-*_task-rest_run-02_*confounds_regressors.tsv .
mmv "*.tsv" "#1.txt"


