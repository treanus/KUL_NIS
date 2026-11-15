% SPM thresholding: Bizzi/Roessler presurgical fMRI method
% Reference: Bizzi et al. (2008) NeuroImage 41:1211-1223
%            Roessler et al. (2005) Clin Neurophysiol 116:2309-2315
%
% Primary map : uncorrected p < 0.001, cluster extent k >= 10 voxels
% Secondary map: FWE-corrected p < 0.05 (Random Field Theory, no cluster extent)
%
% Substitutions performed by KUL_fmriproc_spm.sh before calling:
%   ###FMRIRESULTS###  -> absolute path to the SPM results directory

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

fmriresults_dir = '###FMRIRESULTS###';
cd(fmriresults_dir);

load('SPM.mat');
df = SPM.xX.erdf;

p_unc  = 0.001;
k_min  = 50;
p_FWE  = ###PFWE###;

% Height thresholds
t_thresh_unc = tinv(1 - p_unc, df);
t_thresh_FWE = spm_uc(p_FWE, [1 df], 'T', SPM.xVol.R, 1, SPM.xVol.S);

fprintf('\n=== Bizzi/Roessler Thresholding ===\n');
fprintf('Residual df            = %.1f\n', df);
fprintf('Uncorrected p<%.3f    ->  t > %.3f\n', p_unc, t_thresh_unc);
fprintf('FWE-corrected p<%.2f  ->  t > %.3f\n', p_FWE, t_thresh_FWE);

Vt = spm_vol('spmT_0001.nii');
Yt = spm_read_vols(Vt);
Yt(isnan(Yt)) = 0;
Yt(isinf(Yt)) = 0;

% --- Primary: uncorrected p<0.001, k>=10 voxels ---
Y_unc = Yt;
Y_unc(Y_unc < t_thresh_unc) = 0;
[L, num] = spm_bwlabel(double(Y_unc > 0), 26);
n_removed = 0;
for c = 1:num
    idx = find(L == c);
    if numel(idx) < k_min
        Y_unc(idx) = 0;
        n_removed = n_removed + 1;
    end
end
n_survived = num - n_removed;
fprintf('p<0.001 unc: %d clusters survived (k>=%d), %d sub-threshold clusters removed\n', ...
    n_survived, k_min, n_removed);

Vout = Vt;
Vout.fname = fullfile(fmriresults_dir, sprintf('spmT_0001_p001unc_k%d.nii', k_min));
spm_write_vol(Vout, Y_unc);

% --- Secondary: FWE p<###PFWE###, cluster extent k>=k_min ---
Y_FWE = Yt;
Y_FWE(Y_FWE < t_thresh_FWE) = 0;
[L_FWE, num_FWE] = spm_bwlabel(double(Y_FWE > 0), 26);
n_FWE_removed = 0;
for c = 1:num_FWE
    idx = find(L_FWE == c);
    if numel(idx) < k_min
        Y_FWE(idx) = 0;
        n_FWE_removed = n_FWE_removed + 1;
    end
end
n_FWE_survived = num_FWE - n_FWE_removed;
fprintf('FWE p<%.3f: %d clusters survived (k>=%d), %d sub-threshold removed\n', ...
    p_FWE, n_FWE_survived, k_min, n_FWE_removed);

% build output filename from p_FWE value (e.g. 0.01 -> FWE01, 0.005 -> FWE005)
pfwe_str = strrep(sprintf('%.3f', p_FWE), '0.', '');
pfwe_str = regexprep(pfwe_str, '0+$', '');  % trim trailing zeros
fwe_fname = fullfile(fmriresults_dir, sprintf('spmT_0001_FWE%s_k%d.nii', pfwe_str, k_min));

Vout2 = Vt;
Vout2.fname = fwe_fname;
spm_write_vol(Vout2, Y_FWE);

fprintf('Written: %s\n', Vout.fname);
fprintf('Written: %s\n', Vout2.fname);
