%-----------------------------------------------------------------------
% Job saved on 16-Dec-2020 18:48:24 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'###FMRIDIR###'};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '###FMRIFILE###';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{2}.spm.stats.fmri_spec.dir = {'###FMRIRESULTS###'};
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 1.5;
matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('File Selector (Batch Mode): Selected Files (###FMRIFILE###)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.name = 'ON';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.onset = [30
                                                      90
                                                      150
                                                      210
                                                      270];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.duration = 30;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.orth = 1;
matlabbatch{2}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{2}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = {''};
matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{2}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{2}.spm.stats.fmri_spec.volt = 1;
matlabbatch{2}.spm.stats.fmri_spec.global = 'None';
matlabbatch{2}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{2}.spm.stats.fmri_spec.mask = {''};
matlabbatch{2}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{4}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'ON';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.delete = 1;
