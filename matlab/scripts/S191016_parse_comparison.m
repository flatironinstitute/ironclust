% Example
S_score1 = load('~/ceph/DanEnglish/juxta_cell_output/m14_190326_155432/raw_geom_score.mat')
%   struct with fields:
% 
%                 P: [1×1 struct]
%           S_burst: []
%              S_gt: [1×1 struct]
%         S_overlap: []
%       S_score_clu: [1×1 struct]
%     S_score_ksort: []
%         cviSpk_gt: {[3556×1 int32]}
%       cviTime_clu: {}
%        miCluMatch: [66×1 double]
%              mrFp: [66×1 single]
%            mrMiss: [66×1 single]
%          trWav_gt: [61×32 single]
%         viSite_gt: 10
%           vnCluGt: 3556
%         vnSite_gt: 0
%          vrSnr_gt: 2.9096
%      vrSnr_min_gt: 2.9096
%       vrSnr_sd_gt: 0.5934
%        vrVmin_clu: 77.6700
%         vrVmin_gt: 77.6700
%         vrVpp_clu: 102.8782
%       vrVrms_site: []

%%
% set the ground truth unit index
iGt1 = 1;

% Find the peak channel (which has the most negative peak)
iChan_max = S_score1.viSite_gt(iGt1); 

% find the SNR of the groundtruth unit 
snr_gt1 = S_score1.vrSnr_gt(iGt1);

% find the best matching sorted unit for the given ground truth unit
iClu1 = S_score1.miCluMatch(1,iGt1);

% find various scores for the best matching unit
accuracy1 = S_score1.S_score_clu.vrAccuracy(iGt1); % numel(match) / numel(union(sorted_unit, true_unit))
false_positive1 = S_score1.S_score_clu.vrFp(iGt1);
false_negative1 = S_score1.S_score_clu.vrMiss(iGt1);

%% plot the groundtruth waveform
figure; 
subplot 121; 
title(sprintf('Groundtruth unit %d', iGt1)); grid on;
vx_plot = S_score1.vrSnr_gt
plot(vx_plot, S_score1.trWav_gt(:,:,iGt1));

