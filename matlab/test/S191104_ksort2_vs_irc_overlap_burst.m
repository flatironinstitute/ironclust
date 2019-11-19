%% OVERLAP ANALYSIS 

SNR_THRESH = 4;

vcFile_mat_irc = 'C:\tmp\irc2\hybrid_synth\static_siprobe\rec_64c_1200s_11\raw_geom_score.mat';
vcFile_mat_ksort2 = 'C:\tmp\ksort2\hybrid_synth\static_siprobe\rec_64c_1200s_11\raw_geom_score.mat';
S_irc = load(vcFile_mat_irc);
S_ksort2 = load(vcFile_mat_ksort2);
% mnCnt_gt = cell2mat(cellfun(@(x)histcounts(x, 1:6)', S_irc.S_overlap.cvnOverlap_gt, 'UniformOutput', 0));
% mrFrac_gt = mnCnt_gt./sum(mnCnt_gt,1); mrFrac_gt=mrFrac_gt(:,vl_gt)';
% 
% % total error due to misassignment. probablity 
% vl_gt = S_irc.vrSnr_gt >= SNR_THRESH;
% mr_err_irc = (1 - S_irc.S_overlap.mpHit_overlap_gt(vl_gt,:)) .* mrFrac_gt;
% mr_err_ksort2 = (1 - S_ksort2.S_overlap.mpHit_overlap_gt(vl_gt,:)) .* mrFrac_gt;
% 
% fh_mean = @(x)nanmean(x(:));
% fh_bar = @(x)[fh_mean(x(:,1)), fh_mean(x(:,1:2)), fh_mean(x)];
% figure; hold on;
% bar([fh_bar(mr_err_irc); fh_bar(mr_err_ksort2)]);

% mean(S_irc.S_score_clu.vrAccuracy(vl_gt))
% mean(S_ksort2.S_score_clu.vrAccuracy(vl_gt))

mrBar = [];
for i=1:2
    switch i
        case 1, S1=S_irc;
        case 2, S1=S_ksort2;
    end    
    viClu_gt = find(S1.vrSnr_gt >= SNR_THRESH);
%     vl_gtspk = ismember(S1.S_gt.viClu, viClu_gt);
    vl = cell2mat(S1.S_score_clu.cvlHit_gt(viClu_gt)');
    vn = cell2mat(S1.S_overlap.cvnOverlap_gt(viClu_gt)');
%     vl = S1.S_overlap.vlHit_gtspk(vl_gtspk);
%     vn = S1.S_overlap.vnOverlap_gt(vl_gtspk);
    score1 = (1-[mean(vl(vn==1)), mean(vl(vn<=2)), mean(vl(vn<=3)), mean(vl(vn<=inf))]);
    mrBar = [mrBar; score1];
    vn_data = cumsum(histcounts(vn, 1:5));
end
figure('Color','w'); bar(mrBar',1);
csXlabel0 = {'0 (%d)', '<=1 (%d)', '<=2 (%d)', 'any (%d)'};
csXlabel = arrayfun(@(x)sprintf(csXlabel0{x}, vn_data(x)), 1:numel(vn_data), 'UniformOutput', 0);
set(gca,'XTickLabel', csXlabel);
xlabel('# Overlapping spikes');
grid on;
ylabel('Error rate');
legend({'IronClust','KiloSort2'}, 'location', 'NE');
ylim([0 .08]);
title(sprintf('Hybrid static (SNR>=%d)',SNR_THRESH));

%% BURST analysis using bionet
SNR_THRESH = 8;

vcFile_mat_irc = 'C:\tmp\irc2\bionet\bionet_static\static_8x_A_2A\raw_geom_score.mat';
vcFile_mat_ksort2 = 'C:\tmp\ksort2\bionet\bionet_static\static_8x_A_2A\raw_geom_score.mat';
S_irc = load(vcFile_mat_irc);
S_ksort2 = load(vcFile_mat_ksort2);

mrBar = [];
for i=1:2
    switch i
        case 1, S1=S_irc;
        case 2, S1=S_ksort2;
    end    
    viClu_gt = find(S1.vrSnr_gt >= SNR_THRESH);
%     vl_gtspk = ismember(S1.S_gt.viClu, viClu_gt);
    vl = cell2mat(S1.S_score_clu.cvlHit_gt(viClu_gt)');
    vn = cell2mat(S1.S_gt.cvnBurst_clu(viClu_gt));
%     vl = S1.S_overlap.vlHit_gtspk(vl_gtspk);
%     vn = S1.S_overlap.vnOverlap_gt(vl_gtspk);
    score1 = 1-[mean(vl(vn==0)), mean(vl(vn<=1)), mean(vl(vn<=2)), mean(vl(vn<=inf))];
    mrBar = [mrBar; score1];
    vn_data = cumsum(histcounts(vn, 0:4));
end
figure('Color','w'); bar(mrBar',1);
csXlabel0 = {'0 (%d)', '<=1 (%d)', '<=2 (%d)', 'any (%d)'};
csXlabel = arrayfun(@(x)sprintf(csXlabel0{x}, vn_data(x)), 1:numel(vn_data), 'UniformOutput', 0);
set(gca,'XTickLabel', csXlabel);
xlabel('Burst index');
grid on;
ylabel('Error rate');
legend({'IronClust','KiloSort2'}, 'Location', 'NE');
ylim([0 .1]);
title(sprintf('BIONET static (SNR>=%d)',SNR_THRESH));


% plot(S_irc.vrSnr_gt, S_irc.S_overlap.mpHit_overlap_gt(:,3), 's');
% figure; hold on;
% plot(S_irc.vrSnr_gt, S_irc.S_score_clu.vrAccuracy, 'o');
% 
% figure; hold on;
% plot(S_ksort2.vrSnr_gt, S_ksort2.S_overlap.mpHit_overlap_gt(:,1), 'o');
% plot(S_ksort2.vrSnr_gt, S_ksort2.S_overlap.mpHit_overlap_gt(:,2), 'x');
% % plot(S_ksort2.vrSnr_gt, S_ksort2.S_overlap.mpHit_overlap_gt(:,3), 's');
% 
% figure; hold on;
% plot(S_ksort2.vrSnr_gt, S_ksort2.S_score_clu.vrAccuracy, 'o');

%% bursting comparison