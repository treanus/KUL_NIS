% List of open inputs

spm_path = [getenv('KUL_apps_DIR') filesep 'spm12'];
addpath(spm_path)
nrun = 1; % enter the number of runs here
jobfile = {'###JOBFILE###'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
