%-----------------------------------------------------------------------
% Job saved on 18-Dec-2020 11:29:19 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'###FMRIDIR###'};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '###FMRIFILE###_run-1';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'###FMRIDIR###'};
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '###FMRIFILE###_run-2';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{3}.spm.stats.fmri_spec.dir = {'###FMRIRESULTS###'};
matlabbatch{3}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{3}.spm.stats.fmri_spec.timing.RT = 1.5;
matlabbatch{3}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{3}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{3}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('File Selector (Batch Mode): Selected Files (###FMRIFILE###_run-1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.stats.fmri_spec.sess(1).cond.name = 'ON';
matlabbatch{3}.spm.stats.fmri_spec.sess(1).cond.onset = [30
                                                         90
                                                         150
                                                         210
                                                         270];
matlabbatch{3}.spm.stats.fmri_spec.sess(1).cond.duration = 30;
matlabbatch{3}.spm.stats.fmri_spec.sess(1).cond.tmod = 0;
matlabbatch{3}.spm.stats.fmri_spec.sess(1).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{3}.spm.stats.fmri_spec.sess(1).cond.orth = 1;
matlabbatch{3}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{3}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{3}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
matlabbatch{3}.spm.stats.fmri_spec.sess(1).hpf = 128;
matlabbatch{3}.spm.stats.fmri_spec.sess(2).scans(1) = cfg_dep('File Selector (Batch Mode): Selected Files (###FMRIFILE###_run-2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.stats.fmri_spec.sess(2).cond.name = 'ON';
matlabbatch{3}.spm.stats.fmri_spec.sess(2).cond.onset = [30
                                                         90
                                                         150
                                                         210
                                                         270];
matlabbatch{3}.spm.stats.fmri_spec.sess(2).cond.duration = 30;
matlabbatch{3}.spm.stats.fmri_spec.sess(2).cond.tmod = 0;
matlabbatch{3}.spm.stats.fmri_spec.sess(2).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{3}.spm.stats.fmri_spec.sess(2).cond.orth = 1;
matlabbatch{3}.spm.stats.fmri_spec.sess(2).multi = {''};
matlabbatch{3}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch{3}.spm.stats.fmri_spec.sess(2).multi_reg = {''};
matlabbatch{3}.spm.stats.fmri_spec.sess(2).hpf = 128;
matlabbatch{3}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{3}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{3}.spm.stats.fmri_spec.volt = 1;
matlabbatch{3}.spm.stats.fmri_spec.global = 'None';
matlabbatch{3}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{3}.spm.stats.fmri_spec.mask = {''};
matlabbatch{3}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{4}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{4}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{5}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.con.consess{1}.tcon.name = 'ON';
matlabbatch{5}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{5}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{5}.spm.stats.con.delete = 1;
