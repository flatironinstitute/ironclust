vcDir0 = '/mnt/home/jjun/ceph/groundtruth/paired_recordings/english/';
S = dir(vcDir0);
csFiles = {S.name}; csFiles=csFiles([S.isdir]);
csFiles = setdiff(csFiles, {'.','..'});
csFiles = cellfun(@(x)fullfile(vcDir0, x), csFiles, 'UniformOutput', 0);

rec_duration_sec = [7.0262e+02   1.1190e+03   1.2490e+03   5.7872e+02   4.6000e+02   7.5000e+02   2.8790e+03];
    
%% run ironclust
for iFile = 1:numel(csFiles)
    try
        irc2(csFiles{iFile});
    catch
        fprintf(2, '%s: %s\n', csFiles{iFile}, lasterr());
    end
end

% SNR, FP, FN, Accuracy (step_sec_drift=10s, 300s). ksort2 uses 1s batch size
mr_snr_fp_fn_acc_irc2 = [2.9 5.8 11.4 84; 7.3 8.7 15.5 78.2; 2 74.8 33.7 22.3; 5.7 13.3 35.2 59; .6 9 46.6 50.7; 4.5 9.2 11.7 81; 3.3 56.2 18 40];
runtime_irc2 = [46.2, 41.8, 58.4, 33.9 14.9, 33.7, 78.4];

% SNR, FP, FN, Accuracy (step_sec_drift=20s, 300s)
% mr_snr_fp_fn_acc_irc2 = [2.9 6.6 21.5 74.4; 7.3 6.1 26 70.6; 2 25 38.1 51.3; 5.7 18.7 61.4 35.5; .6 14.9 46.7 48.7; 4.5 9.5 11.6 80.9; 3.3 80.6 13.9 18.8];

% SNR, FP, FN, Accuracy (step_sec_drift=2s, 300s). ksort2 uses 1s batch size
% mr_snr_fp_fn_acc_irc2 = [78.8; 76.1; 15.5; 41.9; 50.6; 81; 18.8];
% runtime_irc2 = [56.5 81.2 124 52.6 22.3 56.3 151.3];

% SNR, FP, FN, Accuracy (step_sec_drift=8s, 240s). ksort2 uses 1s batch size
% mr_snr_fp_fn_acc_irc2 = [2.9 5.4 22.8 74; 7.3 8.4 16.4 77.6; 2 71.0 34.5 25.2; 5.7 13.5 34.9 59.1; 0.6 5.9 47.6 50.8; 4.5 9.4 12 80.6; 3.3 55.9 18.1 40.2];
% runtime_irc2 = [44.7, 89.7, 64.7, 40.9, 18.6, 44.7, 89.7];

%% run kilosort
% kilosort2 CUDA code needs to be recompiled.
% Before running matlab, run
% `module load cuda gcc`
%
% in matlab, run
% `setenv('MW_NVCC_PATH', '/cm/shared/sw/pkg/devel/cuda/10.0.130_410.48/bin/')`
for iFile = 1:numel(csFiles)
    try
        run_ksort2(csFiles{iFile});
    catch
        fprintf(2, '%s: %s\n', csFiles{iFile}, lasterr());
    end
end

% SNR, FP, FN, Accuracy
mr_snr_fp_fn_acc_ksort2 = [2.9 3.7 2.8 93.7; 7.3 .7 16.9 82.5; 2 23.3 33.0 55.7; 5.7 3.8 39.8 58.8; .6 25.2 35.9 52.7; 4.5 8.8 7.4 85; 3.3 10.9 50.3 46.9];
runtime_ksort2 = [143.4 171.8 208.8 87.8 61.4 116.8 192.5];

%% comparison plot
DETECTION_THRESH = 4;
text_ = @(x,y,c)text(x,y,c,'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
text1_ = @(x)text(4, x, 'threshold (=\theta)', 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', 'r');
text2_ = @(x,y)text(x, y, {'o    IronClust2', 'x    KiloSort2'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
line01_ = @()line([0 1], [0 1], 'Color','k');
legend01_ = @()legend({'SNR<\theta', 'SNR\geq\theta'}, 'Location', 'NorthWest');
legend_snr_ = @(x)legend({'IronClust2<\theta', 'IronClust2\geq\theta', 'KiloSort2<\theta', 'KiloSort2\geq\theta'}, 'Location', x);

figure('Color','w'); 
vrSnr_gt = mr_snr_fp_fn_acc_irc2(:,1);
vl_subthresh = vrSnr_gt<DETECTION_THRESH;

% accuracy
subplot 241; cla; hold on;
plot(vrSnr_gt(vl_subthresh), mr_snr_fp_fn_acc_irc2(vl_subthresh,end)/100, 'ro');
plot(vrSnr_gt(~vl_subthresh), mr_snr_fp_fn_acc_irc2(~vl_subthresh,end)/100, 'ko');
plot(vrSnr_gt(vl_subthresh), mr_snr_fp_fn_acc_ksort2(vl_subthresh,end)/100, 'rx');
plot(vrSnr_gt(~vl_subthresh), mr_snr_fp_fn_acc_ksort2(~vl_subthresh,end)/100, 'kx');
line([DETECTION_THRESH DETECTION_THRESH], [0 1], 'Color', 'r'); text1_(.7);
axis([0 10 0 1]); grid on; title('Paired recording groundtruth (32chan.)');
xlabel('SNR (Vmin/Vrms)'); ylabel('Accuracy'); legend_snr_('SouthEast');

subplot 245; cla; hold on;
plot(mr_snr_fp_fn_acc_irc2(vl_subthresh,end)/100, mr_snr_fp_fn_acc_ksort2(vl_subthresh,end)/100, 'ro');
plot(mr_snr_fp_fn_acc_irc2(~vl_subthresh,end)/100, mr_snr_fp_fn_acc_ksort2(~vl_subthresh,end)/100, 'ko');
axis([0 1 0 1]); grid on; line01_(); 
xlabel('IronClust2'); ylabel('KiloSort2'); title('Accuracy'); legend01_();

% false positive
subplot 242; cla; hold on;
plot(vrSnr_gt(vl_subthresh), mr_snr_fp_fn_acc_irc2(vl_subthresh,2)/100, 'ro');
plot(vrSnr_gt(~vl_subthresh), mr_snr_fp_fn_acc_irc2(~vl_subthresh,2)/100, 'ko');
plot(vrSnr_gt(vl_subthresh), mr_snr_fp_fn_acc_ksort2(vl_subthresh,2)/100, 'rx');
plot(vrSnr_gt(~vl_subthresh), mr_snr_fp_fn_acc_ksort2(~vl_subthresh,2)/100, 'kx');
axis([0 10 0 1]); grid on;
line([DETECTION_THRESH DETECTION_THRESH], [0 1], 'Color', 'r'); text1_(.7);
xlabel('SNR (Vmin/Vrms)'); ylabel('False Positives'); legend_snr_('NorthEast');

subplot 246; cla; hold on;
plot(mr_snr_fp_fn_acc_irc2(vl_subthresh,2)/100, mr_snr_fp_fn_acc_ksort2(vl_subthresh,2)/100, 'ro');
plot(mr_snr_fp_fn_acc_irc2(~vl_subthresh,2)/100, mr_snr_fp_fn_acc_ksort2(~vl_subthresh,2)/100, 'ko');
axis([0 1 0 1]); grid on; line01_(); 
xlabel('IronClust2'); ylabel('KiloSort2'); title('False Positives'); legend01_();

% false negative
subplot 243; cla; hold on;
plot(vrSnr_gt(vl_subthresh), mr_snr_fp_fn_acc_irc2(vl_subthresh,3)/100, 'ro');
plot(vrSnr_gt(~vl_subthresh), mr_snr_fp_fn_acc_irc2(~vl_subthresh,3)/100, 'ko');
plot(vrSnr_gt(vl_subthresh), mr_snr_fp_fn_acc_ksort2(vl_subthresh,3)/100, 'rx');
plot(vrSnr_gt(~vl_subthresh), mr_snr_fp_fn_acc_ksort2(~vl_subthresh,3)/100, 'kx');
axis([0 10 0 1]); grid on;
line([DETECTION_THRESH DETECTION_THRESH], [0 1], 'Color', 'r'); text1_(.7);
xlabel('SNR (Vmin/Vrms)'); ylabel('False Negatives'); legend_snr_('NorthEast');

subplot 247; hold on;
plot(mr_snr_fp_fn_acc_irc2(vl_subthresh,3)/100, mr_snr_fp_fn_acc_ksort2(vl_subthresh,3)/100, 'ro');
plot(mr_snr_fp_fn_acc_irc2(~vl_subthresh,3)/100, mr_snr_fp_fn_acc_ksort2(~vl_subthresh,3)/100, 'ko');
grid on; axis([0 1 0 1]); line01_(); 
xlabel('IronClust2'); ylabel('KiloSort2'); title('False Negatives'); legend01_();

% speed
subplot 244; cla; hold on;
plot(rec_duration_sec, runtime_irc2, 'ko');
plot(rec_duration_sec, runtime_ksort2, 'kx');
xlabel('Recording duration (sec)'); ylabel('Runtime (sec)');
legend({'IronClust2', 'KiloSort2'}, 'Location', 'SouthEast');
axis ([0 3000 0 250]); grid on;
title('IronClust2: GPU + 12 Cores; KiloSort2: GPU');

subplot 248; cla; hold on;
switch 2
    case 1
        plot(rec_duration_sec./runtime_irc2, rec_duration_sec./runtime_ksort2, 'ko');
        title('Speed (x Realtime)');
        axis([0 40 0 40]); grid on;
        plot([0 40], [0 40], 'k-');
        plot([0 40], [0 40/3], 'k--'); text_(30, 10, '3x');
    case 2
        plot(runtime_irc2, runtime_ksort2, 'ko');
        title('Runtime (s)');
        axis([0 250 0 250]); grid on;
        plot([0 250], [0 250], 'k-');
        plot([0 250], [0 250*3], 'k--'); text_(225/3, 225, '3x');
end
xlabel('IronClust2'); ylabel('KiloSort2');