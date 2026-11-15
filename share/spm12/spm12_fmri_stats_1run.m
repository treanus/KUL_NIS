% List of open inputs

% SPM12 location. KUL_MATLAB_APPS is exported by the installer's environment
% block ($SOFTWARE_ROOT/src/matlab_apps); KUL_apps_DIR is the legacy name from
% the old /usr/local/KUL_apps layout, kept as a fallback for existing installs.
% Both are checked explicitly because the previous one-liner addpath()'d
% '/spm12' without complaint whenever the variable was unset, and the run only
% failed later and obscurely, at the first spm() call.
matlab_apps = getenv('KUL_MATLAB_APPS');
if isempty(matlab_apps)
    matlab_apps = getenv('KUL_apps_DIR');
end
if isempty(matlab_apps)
    error('KUL:noMatlabApps', ['Neither KUL_MATLAB_APPS nor KUL_apps_DIR is set. ' ...
        'Source the KUL environment block, or run the installer spm12 section.']);
end
spm_path = fullfile(matlab_apps, 'spm12');
if exist(fullfile(spm_path, 'spm.m'), 'file') ~= 2
    error('KUL:noSPM', 'SPM12 not found at %s (no spm.m there).', spm_path);
end
addpath(spm_path)
nrun = 1; % enter the number of runs here
jobfile = {'###JOBFILE###'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
