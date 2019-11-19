%%
run_ksort2 C:\tmp\groundtruth\hybrid_synth\static_siprobe\rec_64c_1200s_11;
run_ksort2 C:\tmp\groundtruth\hybrid_synth\drift_siprobe\rec_64c_1200s_11;
irc2 clear; irc2 test-static;
irc2 clear; irc2 test-drift;

%%

csFiles = {...
    'C:\tmp\irc2\hybrid_synth\static_siprobe\rec_64c_1200s_11\raw_geom_score.mat', 
    'C:\tmp\ksort2\hybrid_synth\static_siprobe\rec_64c_1200s_11\raw_geom_score.mat',
    'C:\tmp\irc2\hybrid_synth\drift_siprobe\rec_64c_1200s_11\raw_geom_score.mat',
    'C:\tmp\ksort2\hybrid_synth\drift_siprobe\rec_64c_1200s_11\raw_geom_score.mat'};

S_static_irc2 = load(csFiles{1});
S_static_ksort2 = load(csFiles{2});
S_drift_irc2 = load(csFiles{3});
S_drift_ksort2 = load(csFiles{4});

%%
SNR_THRESH = 6;
fPlot_subthresh = 0;

vlSnr = S_static_irc2.vrSnr_gt>SNR_THRESH;
figure('Color','W');
for iPlot=1:6
    subplot(2,3,iPlot); hold on;
    switch iPlot
        case 1, [vrX, vrY] = deal(S_static_irc2.S_score_clu.vrAccuracy, S_static_ksort2.S_score_clu.vrAccuracy); vcTitle='Accuracy (static)';
        case 4, [vrX, vrY] = deal(S_drift_irc2.S_score_clu.vrAccuracy, S_drift_ksort2.S_score_clu.vrAccuracy); vcTitle='Accuracy (drift)';
        case 2, [vrX, vrY] = deal(1-S_static_irc2.S_score_clu.vrFp, 1-S_static_ksort2.S_score_clu.vrFp); vcTitle='Precision (static)';
        case 5, [vrX, vrY] = deal(1-S_drift_irc2.S_score_clu.vrFp, 1-S_drift_ksort2.S_score_clu.vrFp); vcTitle='Precision (drift)';
        case 3, [vrX, vrY] = deal(1-S_static_irc2.S_score_clu.vrMiss, 1-S_static_ksort2.S_score_clu.vrMiss); vcTitle='Recall (static)';
        case 6, [vrX, vrY] = deal(1-S_drift_irc2.S_score_clu.vrMiss, 1-S_drift_ksort2.S_score_clu.vrMiss); vcTitle='Recall (drift)';
    end
    plot(vrX(vlSnr), vrY(vlSnr), 'k.'); 
    if fPlot_subthresh, plot(vrX(~vlSnr), vrY(~vlSnr), 'r.'); end
    xlabel('IronClust'); ylabel('KiloSort2'); grid on; axis([.1 1 .1 1]); plot([.1 1], [.1 1], 'k-');
    set(gca,'XScale','log','YScale','log');
    title(vcTitle); axis square;
end

%% runtime comparison plot
t_irc2_static = [51.4+7.8, 11.2];
t_irc2_drift = [47.7+6.1, 15.7];
t_ksort2_static = [114.905239, 255.6-114.905239];
t_ksort2_drirt = [121.4, 262.2-121.4];
figure('Color','w'); 
bar([t_irc2_static; t_irc2_drift; t_ksort2_static; t_ksort2_drirt], 'stacked');
grid on;
ylabel('Runtime (s)');
set(gca, 'XTickLabel', {'IronClust-Static', 'IronClust-Drift', 'KiloSort2-Static', 'KiloSort2-Drift'});
xtickangle(30);
legend({'detect-sort', 'auto-curation'}, 'location', 'NW');
set(gca, 'XLim', [.5 4.5]);
title('Hybrid groundtruth (64ch, 1200s)');