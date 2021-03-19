c1 = dir('ALL_CSV/sub-P*_ses-*_T1w_iso_biascorrected_histogram.csv');
c2 = dir('ALL_CSV/sub-P*_ses-*_T1w_iso_biascorrected_calib-lin_histogram.csv');
c3 = dir('ALL_CSV/sub-P*_ses-*_T1w_iso_biascorrected_calib-nonlin_histogram.csv');
c4 = dir('ALL_CSV/sub-P*_ses-*_T1w_iso_biascorrected_calib-nonlin2_histogram.csv');
c4b = dir('ALL_CSV/sub-P*_ses-*_T1w_iso_biascorrected_calib-nonlin2b_histogram.csv');
c4c = dir('ALL_CSV/sub-P*_ses-*_T1w_iso_biascorrected_calib-nonlin3_histogram.csv');
c5 = dir('ALL_CSV/sub-P*_ses-*_T2w_iso_biascorrected_reg2T1w_histogram.csv');
c6 = dir('ALL_CSV/sub-P*_ses-*_T2w_iso_biascorrected_calib-lin_reg2T1w_histogram.csv');
c7 = dir('ALL_CSV/sub-P*_ses-*_T2w_iso_biascorrected_calib-nonlin_reg2T1w_histogram.csv');
c8 = dir('ALL_CSV/sub-P*_ses-*_T2w_iso_biascorrected_calib-nonlin2_reg2T1w_histogram.csv');
c8b = dir('ALL_CSV/sub-P*_ses-*_T2w_iso_biascorrected_calib-nonlin2b_reg2T1w_histogram.csv');
c8c = dir('ALL_CSV/sub-P*_ses-*_T2w_iso_biascorrected_calib-nonlin3_reg2T1w_histogram.csv');

%c9 = dir('ALL_CSV/sub-P*_ses-*_FLAIR_iso_biascorrected_reg2T1w_histogram.csv');
%c10 = dir('ALL_CSV/sub-P*_ses-*_FLAIR_iso_biascorrected_calib-lin_reg2T1w_histogram.csv');
%c11 = dir('ALL_CSV/sub-P*_ses-*_FLAIR_iso_biascorrected_calib-nonlin_reg2T1w_histogram.csv');
%c12 = dir('ALL_CSV/sub-P*_ses-*_FLAIR_iso_biascorrected_calib-nonlin2_reg2T1w_histogram.csv');

close all
hold off
b=1;
e=100;
m=1000000;
xmax=1200;
n=size(c1,1)-1
for i = 1:n
    figure(1)
    f = ['ALL_CSV/' c1(i).name];
    h = csvread(f,1,0);
    subplot(6,1,1)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax])
    hold on

    f = ['ALL_CSV/' c2(i).name];
    h = csvread(f,1,0);
    subplot(6,1,2)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax/6])
    hold on
    
    f = ['ALL_CSV/' c3(i).name];
    h = csvread(f,1,0);
    subplot(6,1,3)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax/6])
    hold on
    
    f = ['ALL_CSV/' c4(i).name];
    h = csvread(f,1,0);
    subplot(6,1,4)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax])
    hold on
    
    f = ['ALL_CSV/' c4b(i).name];
    h = csvread(f,1,0);
    subplot(6,1,5)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax])
    hold on
    
    f = ['ALL_CSV/' c4c(i).name];
    h = csvread(f,1,0);
    subplot(6,1,6)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax/6])
    hold on
    
    figure(2)
    f = ['ALL_CSV/' c5(i).name];
    h = csvread(f,1,0);
    subplot(6,1,1)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax])
    hold on

    f = ['ALL_CSV/' c6(i).name];
    h = csvread(f,1,0);
    subplot(6,1,2)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax])
    hold on
    
    f = ['ALL_CSV/' c7(i).name];
    h = csvread(f,1,0);
    subplot(6,1,3)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax])
    hold on
    
    f = ['ALL_CSV/' c8(i).name];
    h = csvread(f,1,0);
    subplot(6,1,4)
    plot(h(1,b:e),h(2,b:e))
    ylim([0 m]); xlim([0 xmax])
    hold on
    
%     f = ['ALL_CSV/' c8b(i).name];
%     h = csvread(f,1,0);
%     subplot(6,1,5)
%     plot(h(1,b:e),h(2,b:e))
%     ylim([0 m]); xlim([-100 xmax])
%     hold on
%     
%     f = ['ALL_CSV/' c8c(i).name];
%     h = csvread(f,1,0);
%     subplot(6,1,6)
%     plot(h(1,b:e),h(2,b:e))
%     ylim([0 m]); xlim([-100 xmax])
%     hold on
    
%     figure(3)
%     f = ['ALL_CSV/' c9(i).name];
%     h = csvread(f,1,0);
%     subplot(4,1,1)
%     plot(h(1,b:e),h(2,b:e))
%     ylim([0 m]); xlim([-100 xmax])
%     hold on
% 
%     f = ['ALL_CSV/' c10(i).name];
%     h = csvread(f,1,0);
%     subplot(4,1,2)
%     plot(h(1,b:e),h(2,b:e))
%     ylim([0 m]); xlim([-100 xmax])
%     hold on
%     
%     f = ['ALL_CSV/' c11(i).name];
%     h = csvread(f,1,0);
%     subplot(4,1,3)
%     plot(h(1,b:e),h(2,b:e))
%     ylim([0 m]); xlim([-100 xmax])
%     hold on
%     
%     f = ['ALL_CSV/' c12(i).name];
%     h = csvread(f,1,0);
%     subplot(4,1,4)
%     plot(h(1,b:e),h(2,b:e))
%     ylim([0 m]); xlim([-100 xmax])
%     hold on
    
end
