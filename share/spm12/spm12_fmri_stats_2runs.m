% List of open inputs
addpath /home/stefan/KUL_apps/spm12
nrun = 1; % enter the number of runs here
jobfile = {'###JOBFILE###'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
