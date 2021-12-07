%-----------------------------------------------------------------------
% Job saved on 06-Dec-2021 17:00:40 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'###FMRIDIR###'};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '###FMRIFILE###_run-1';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'###FMRIDIR###'};
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '###FMRIFILE###_run-2';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'###FMRIDIR###'};
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '###FMRIFILE###_run-3';
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{4}.spm.stats.fmri_spec.dir = {'###FMRIRESULTS###'};
matlabbatch{4}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{4}.spm.stats.fmri_spec.timing.RT = ###TR###;
matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{4}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('File Selector (Batch Mode): Selected Files (TAAL_run-1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.spm.stats.fmri_spec.sess(1).cond.name = 'ON';
matlabbatch{4}.spm.stats.fmri_spec.sess(1).cond.onset = [30
                                                         90
                                                         150
                                                         210
                                                         270];
matlabbatch{4}.spm.stats.fmri_spec.sess(1).cond.duration = 30;
matlabbatch{4}.spm.stats.fmri_spec.sess(1).cond.tmod = 0;
matlabbatch{4}.spm.stats.fmri_spec.sess(1).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{4}.spm.stats.fmri_spec.sess(1).cond.orth = 1;
matlabbatch{4}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{4}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{4}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
matlabbatch{4}.spm.stats.fmri_spec.sess(1).hpf = 128;
matlabbatch{4}.spm.stats.fmri_spec.sess(2).scans(1) = cfg_dep('File Selector (Batch Mode): Selected Files (TAAL_run-2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.spm.stats.fmri_spec.sess(2).cond.name = 'ON';
matlabbatch{4}.spm.stats.fmri_spec.sess(2).cond.onset = [30
                                                         90
                                                         150
                                                         210
                                                         270];
matlabbatch{4}.spm.stats.fmri_spec.sess(2).cond.duration = 30;
matlabbatch{4}.spm.stats.fmri_spec.sess(2).cond.tmod = 0;
matlabbatch{4}.spm.stats.fmri_spec.sess(2).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{4}.spm.stats.fmri_spec.sess(2).cond.orth = 1;
matlabbatch{4}.spm.stats.fmri_spec.sess(2).multi = {''};
matlabbatch{4}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch{4}.spm.stats.fmri_spec.sess(2).multi_reg = {''};
matlabbatch{4}.spm.stats.fmri_spec.sess(2).hpf = 128;
matlabbatch{4}.spm.stats.fmri_spec.sess(3).scans(1) = cfg_dep('File Selector (Batch Mode): Selected Files (TAAL_run-3)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.spm.stats.fmri_spec.sess(3).cond.name = 'ON';
matlabbatch{4}.spm.stats.fmri_spec.sess(3).cond.onset = [30
                                                         90
                                                         150
                                                         210
                                                         270];
matlabbatch{4}.spm.stats.fmri_spec.sess(3).cond.duration = 30;
matlabbatch{4}.spm.stats.fmri_spec.sess(3).cond.tmod = 0;
matlabbatch{4}.spm.stats.fmri_spec.sess(3).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{4}.spm.stats.fmri_spec.sess(3).cond.orth = 1;
matlabbatch{4}.spm.stats.fmri_spec.sess(3).multi = {''};
matlabbatch{4}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
matlabbatch{4}.spm.stats.fmri_spec.sess(3).multi_reg = {''};
matlabbatch{4}.spm.stats.fmri_spec.sess(3).hpf = 128;
matlabbatch{4}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{4}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{4}.spm.stats.fmri_spec.volt = 1;
matlabbatch{4}.spm.stats.fmri_spec.global = 'None';
matlabbatch{4}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{4}.spm.stats.fmri_spec.mask = {''};
matlabbatch{4}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{5}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{5}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{6}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.con.consess{1}.tcon.name = 'ON';
matlabbatch{6}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{6}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{6}.spm.stats.con.delete = 1;
