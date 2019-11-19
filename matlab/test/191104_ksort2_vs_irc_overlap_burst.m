SNR_THRESH = 8;

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
    vl_gtspk = ismember(S1.S_gt.viClu, viClu_gt);
    vl = S1.S_overlap.vlHit_gtspk(vl_gtspk);
    vn = S1.S_overlap.vnOverlap_gt(vl_gtspk);
    score1 = (1-[mean(vl(vn==1)), mean(vl(vn<=2)), mean(vl)]);
    mrBar = [mrBar; score1];
end
figure; bar(mrBar');


%%
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