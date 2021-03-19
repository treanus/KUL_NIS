% type
%analysis_type = 'calib-none'; % no calibration
%analysis_type = 'calib-lin'; % Ganzetti eye/musle 
analysis_type = 'calib-nonlin'; % Ganzetti newer non-linear
%analysis_type = 'calib-nonlin2'; % Cappelle & Sunaert - mimics Ganzetti, but no MNI needed
%analysis_type = 'calib-nonlin3'; % Cappelle & Sunaert - use brain tissue, without lesions

% Seems to be 1 real outlier in the data, notably sub-172_ses20161028, this
% is row number 294
remove_outlier = 1;


% define the results file
results_file = ['ALL_' analysis_type '.csv'];

% read the results file into a variable t
t = readtable(results_file,'Delimiter',',','TreatAsEmpty',{'.','NA'});

% here we remove the outlier
if remove_outlier == 1
    display("removing outliner")
    t(294,:) = [];
end

% number of scans
s = size(t,1);

% close all figures
close all

% make a line graph with measures from all subjects
figure(1)
plot(1:s,t.NAWM_lh_mtr,'r','LineWidth',2)
hold on
plot(1:s,t.NAWM_lh_t1t2,'g','LineWidth',2)
plot(1:s,t.NAWM_lh_t1flair,'b','LineWidth',2)
legend('MTR','T1T2','T1FLAIR')
title(['Line Graph of NAWM_L (' analysis_type ')'])
hold off

% Boxplots of volumes
figure(2)
Vol_WM = t.Volume_NAWM_lh+t.Volume_NAWM_rh;
Vol_GM = t.Volume_NAGM_lh+t.Volume_NAGM_rh;
l = {'TIV','WM','GM','CSF'};
boxplot([t.Volume_TIV Vol_GM Vol_WM t.Volume_CSF],'Labels',l)
title(['Boxplot of volumes (in mm^3)'])

figure(3)
Vol_CC = t.Volume_CC_Anterior + t.Volume_CC_Mid_Anterior + t.Volume_CC_Central + t.Volume_CC_Mid_Posterior + t.Volume_CC_Posterior;
l = {'MS lesions','Thalamus_L','Thalamus_R','CC'};
boxplot([t.Volume_MSLesions t.Volume_Thalamus_lh t.Volume_Thalamus_rh Vol_CC],'Labels',l)
title(['Boxplot of volumes (in mm^3)'])

% Boxplots of measurements
figure(4)
l = {'NAWM_L','NAWM_R','NAGM_L','NAGM_R','CSF_lat_L','CSF_lat_R','MSlesion'};
h = boxplot([t.NAWM_lh_mtr t.NAWM_rh_mtr t.NAGM_lh_mtr t.NAGM_rh_mtr t.CSF_lateral_lh_mtr t.CSF_lateral_rh_mtr t.MSlesion_mtr],'Notch','on','Labels',l);
ylabel('MTR')

figure(5)
boxplot([t.NAWM_lh_t1t2 t.NAWM_rh_t1t2 t.NAGM_lh_t1t2 t.NAGM_rh_t1t2 t.CSF_lateral_lh_t1t2 t.CSF_lateral_rh_t1t2 t.MSlesion_t1t2],'Notch','on','Labels',l)
ylabel('T1T2')

figure(6)
boxplot([t.NAWM_lh_t1flair t.NAWM_rh_t1flair t.NAGM_lh_t1flair t.NAGM_rh_t1flair t.CSF_lateral_lh_t1flair t.CSF_lateral_rh_t1flair t.MSlesion_t1flair],'Notch','on','Labels',l)
ylabel('T1FLAIR')

% Plot how several parts of the CC correlate in size
figure(7)
l = {'Anterior', 'Mid-Anterior', 'Central', 'Mid-Posterior', 'Posterior'};
corrplot([t.Volume_CC_Anterior t.Volume_CC_Mid_Anterior t.Volume_CC_Central t.Volume_CC_Mid_Posterior t.Volume_CC_Posterior],'varNames',l)

% Make regression plots (and stats) between Lesion Volume and Lesion
% ratio's
figure(8)
reg_mtr = fitlm(t.Volume_MSLesions, t.MSlesion_mtr)
subplot(3,1,1)
plot(reg_mtr)
xl='MS lesion volume (mm^3)';
yl='MTR';
xlabel(xl)
ylabel(yl)
title ([xl ' vs. ' yl])

reg_t1t2 = fitlm(t.Volume_MSLesions, t.MSlesion_t1t2)
subplot(3,1,2)
plot(reg_t1t2)
yl='T1T2';
xlabel(xl)
ylabel(yl)
title ([xl ' vs. ' yl])

reg_flair = fitlm(t.Volume_MSLesions, t.MSlesion_t1flair)
subplot(3,1,3)
plot(reg_flair)
yl='T1FLAIR';
xlabel(xl)
ylabel(yl)
title ([xl ' vs. ' yl])


% Make regression plots (and stats) between MTR and T1/T2, T1/error closing GZ fileFLAIR in the
% lesion
figure(9)
reg_mtt1 = fitlm(t.MSlesion_mtr,t.MSlesion_t1t2)
subplot(3,1,1)
plot(reg_mtt1)
xl='median MTR';
yl='median T1T2';
xlabel(xl)
ylabel(yl)
title ([xl ' vs. ' yl ' in MS lesions'])

reg_mtfl = fitlm(t.MSlesion_mtr,t.MSlesion_t1flair)
subplot(3,1,2)
plot(reg_mtfl)
xl='median MTR';
yl='median T1FLAIR';
xlabel(xl)
ylabel(yl)
title ([xl ' vs. ' yl ' in MS lesions'])

reg_t2fl = fitlm(t.MSlesion_t1t2,t.MSlesion_t1flair)
subplot(3,1,3)
plot(reg_t2fl)
xl='median T1T2';
yl='median T1FLAIR';
xlabel(xl)
ylabel(yl)
title ([xl ' vs. ' yl ' in MS lesions'])

% Make the figure of D. Pareto et al 2020
mean_mtr_NAWM=mean(mean([t.NAWM_lh_mtr t.NAWM_rh_mtr]));
mean_mtr_NAGM=mean(mean([t.NAGM_lh_mtr t.NAGM_rh_mtr]));
mean_mtr_CSF=mean(mean([t.CSF_lateral_lh_mtr t.CSF_lateral_rh_mtr]));

mean_t1t2_NAWM=mean(mean([t.NAWM_lh_t1t2 t.NAWM_rh_t1t2]));
mean_t1t2_NAGM=mean(mean([t.NAGM_lh_t1t2 t.NAGM_rh_t1t2]));
mean_t1t2_CSF=mean(mean([t.CSF_lateral_lh_t1t2 t.CSF_lateral_rh_t1t2]));

mean_t1fl_NAWM=mean(mean([t.NAWM_lh_t1flair t.NAWM_rh_t1flair]));
mean_t1fl_NAGM=mean(mean([t.NAGM_lh_t1flair t.NAGM_rh_t1flair]));
mean_t1fl_CSF=mean(mean([t.CSF_lateral_lh_t1flair t.CSF_lateral_rh_t1flair]));

figure(9)
subplot(3,1,1)
hold on
scatter(mean_mtr_NAWM,mean_t1t2_NAWM,'ok','LineWidth',2)
text(mean_mtr_NAWM,mean_t1t2_NAWM,'  NAWM')
scatter(mean_mtr_NAGM,mean_t1t2_NAGM,'ok','LineWidth',2)
text(mean_mtr_NAGM,mean_t1t2_NAGM,'  NAGM')
scatter(mean_mtr_CSF,mean_t1t2_CSF,'ok','LineWidth',2)
text(mean_mtr_CSF,mean_t1t2_CSF,'  CSF')

subplot(3,1,2)
hold on
scatter(mean_mtr_NAWM,mean_t1t2_NAWM,'ok','LineWidth',2)
text(mean_mtr_NAWM,mean_t1t2_NAWM,'  NAWM')
scatter(mean_mtr_NAGM,mean_t1t2_NAGM,'ok','LineWidth',2)
text(mean_mtr_NAGM,mean_t1t2_NAGM,'  NAGM')
scatter(mean_mtr_CSF,mean_t1t2_CSF,'ok','LineWidth',2)
text(mean_mtr_CSF,mean_t1t2_CSF,'  CSF')

subplot(3,1,3)
hold on
scatter(mean_t1t2_NAWM,mean_t1fl_NAWM,'ok','LineWidth',2)
text(mean_t1t2_NAWM,mean_t1fl_NAWM,'  NAWM')
scatter(mean_t1t2_NAGM,mean_t1fl_NAGM,'ok','LineWidth',2)
text(mean_t1t2_NAGM,mean_t1fl_NAGM,'  NAGM')
scatter(mean_t1t2_CSF,mean_t1fl_CSF,'ok','LineWidth',2)
text(mean_t1t2_CSF,mean_t1fl_CSF,'  CSF')


