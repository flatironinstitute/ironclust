addpath jsonlab-1.5
vcDir = 'test';
%vcFile1 = 'average_runtime_2019-1015.csv';
vcFile_json = 'analysis_2019-1118.json';
%mr1 = readmatrix(fullfile(vcDir, vcFile1), 'NumHeaderLines', 1);
%tbl1 = readtable(fullfile(vcDir, vcFile1));
% csStudy1 = tbl1.study;
% csSorter = tbl1.Properties.VariableNames(3:end);
if ~exist('S_json')
    if ~exist('spikeforest_analysis.mat')
        S_json = loadjson(fullfile(vcDir, vcFile_json));
        save spikeforest_analysis
    else
        load spikeforest_analysis
    end
end
title_ = @(x)title(x,'Interpreter','none');
figure_ = @()figure('Color','w');

%% compare runtime
vcStudyName1 = 'hybrid_static_siprobe';
% vcStudyName1 = 'hybrid_drift_siprobe';
% vcStudyName1 = 'hybrid_static_tetrode';
% vcStudyName1 = 'hybrid_drift_tetrode';
% vcStudyName1 = {'paired_boyden32c', 'paired_crcns', 'paired_mea64c', 'paired_kampff'};

vS_sar = cell2mat(S_json.StudyAnalysisResults);
if ischar(vcStudyName1)
    vi_sar1 = strcmpi({vS_sar.studyName}, vcStudyName1);
else
    vi_sar1 = ismember({vS_sar.studyName}, vcStudyName1);
end
S_sar1 = vS_sar(vi_sar1);

% selected study 1
csRecording1 = cellfun(@(x)char(x), [S_sar1.recordingNames], 'UniformOutput',0);
vS_sr1 = cell2mat([S_sar1.sortingResults]);
csSorter = {vS_sr1.sorterName};
vrSnr_unit1 = S_sar1.trueSnrs;
viRecording_unit1 = S_sar1.trueRecordingIndices;
n_unit = numel(viRecording_unit1);
n_sorter = numel(csSorter);
n_recording = numel(csRecording1);
cvrAcc_sorter1 = {vS_sr1.accuracies};
cvrPrc_sorter1 = {vS_sr1.precisions};
cvrRcl_sorter1 = {vS_sr1.recalls};
getnum_ = @(x)str2num(x(x>='0' & x<='9'));

try
    vnChans_recording1 = cellfun(@(x)str2num(x(5:6)), csRecording1);
    vrDuration_recording1 = cellfun(@(x)getnum_(cell2mat(regexp(x,'c_\d+s_','match'))), csRecording1);
catch
    vnChans_recording1 = [];
    vrDuration_recording1 = [];
end
% collect runtime
switch 1 
    case 1
        cvr_recording_sorter1 = {vS_sr1.cpuTimesSec};
        mrRun_recording_sorter1 = nan(n_recording, n_sorter);
        for iSorter = 1:numel(cvr_recording_sorter1)
            vr1 = cvr_recording_sorter1{iSorter};
            if ~isempty(vr1)
                mrRun_recording_sorter1(:,iSorter) = vr1;
            end
        end
    case 2 % get it from /SortingResults
        mrRun_recording_sorter1 = zeros(n_recording, n_sorter);
        for iSorter = 1:numel(cvrRun_sorter1)
            vcSorter2 = csSorter{iSorter};
            vl2 = cellfun(@(x)strcmpi(x.sorterName, vcSorter2) & strcmpi(x.studyName, vcStudyName1), S_json.SortingResults);
            if any(vl2)
                vS_sr2 = cell2mat(S_json.SortingResults(vl2));                  
%                 [vl12, vi12] = ismember(csRecordingName1, {vS_sr2.recordingName});            
                mrRun_recording_sorter1(:,iSorter) = [vS_sr2.cpuTimeSec];
            end
        end
end
mrRecSnr_unit1 = [viRecording_unit1(:), vrSnr_unit1(:)];

%
SNR_THRESH = 8;
toVec_ = @(x)x(:);
cell2vec_ = @(x,i)toVec_(x{i});
cmr_AccPrcRcl_unit_sorter1 = cell(size(csSorter));
mr_AccPrcRcl_unit_sorter0 = nan(n_unit, 3);
[mrAcc_recording_sorter1, mrPrc_recording_sorter1, mrRcl_recording_sorter1] = deal(nan(n_recording, n_sorter));
for iSorter = 1:numel(csSorter)
    [c1,c2,c3] = deal(cell2vec_(cvrAcc_sorter1, iSorter), cell2vec_(cvrPrc_sorter1, iSorter), cell2vec_(cvrRcl_sorter1, iSorter));
    if iscell(c1)
        cmr_AccPrcRcl_unit_sorter1{iSorter} = mr_AccPrcRcl_unit_sorter0;
    else
        vl1 = vrSnr_unit1 > SNR_THRESH;
        cmr_AccPrcRcl_unit_sorter1{iSorter} = [c1,c2,c3];
        mrAcc_recording_sorter1(:,iSorter) = grpstats(c1(vl1), viRecording_unit1(vl1));
        mrPrc_recording_sorter1(:,iSorter) = grpstats(c2(vl1), viRecording_unit1(vl1));
        mrRcl_recording_sorter1(:,iSorter) = grpstats(c3(vl1), viRecording_unit1(vl1));
    end       
end %for

% table output
disp(['acc_', vcStudyName1]); disp(array2table(mrAcc_recording_sorter1, 'RowNames', csRecording1, 'VariableNames', csSorter));
disp(['prc_', vcStudyName1]); disp(array2table(mrPrc_recording_sorter1, 'RowNames', csRecording1, 'VariableNames', csSorter));
disp(['rcl_', vcStudyName1]); disp(array2table(mrRcl_recording_sorter1, 'RowNames', csRecording1, 'VariableNames', csSorter));
disp(['run_', vcStudyName1]); disp(array2table(mrRun_recording_sorter1, 'RowNames', csRecording1, 'VariableNames', csSorter));

%% plot by SNR
% figure; hold on; 
% plot(vrSnr_unit1, cell2vec_(cvrAcc_sorter1, 2), 'k.'); 
% plot(vrSnr_unit1, cell2vec_(cvrAcc_sorter1, 4), 'r.');
% ylabel('Accuracy'); legend({'IronClust2', 'KiloSort2'}); grid on;
% 
% figure; hold on; 
% plot(vrSnr_unit1, cell2vec_(cvrPrc_sorter1, 2), 'k.'); 
% plot(vrSnr_unit1, cell2vec_(cvrPrc_sorter1, 4), 'r.');
% ylabel('Precision'); legend({'IronClust2', 'KiloSort2'}); grid on;
% 
% figure; hold on; 
% plot(vrSnr_unit1, cell2vec_(cvrRcl_sorter1, 2), 'k.'); 
% plot(vrSnr_unit1, cell2vec_(cvrRcl_sorter1, 4), 'r.');
% ylabel('Recall'); legend({'IronClust2', 'KiloSort2'}); grid on;

mrAcc_Prc_Rcl_irc2_ks2 = [...
    cell2vec_(cvrAcc_sorter1, 2), cell2vec_(cvrAcc_sorter1, 5), ...
    cell2vec_(cvrPrc_sorter1, 2), cell2vec_(cvrPrc_sorter1, 5), ...
    cell2vec_(cvrRcl_sorter1, 2), cell2vec_(cvrRcl_sorter1, 5)];
cviMode = {[1,2],[3,4],[5,6]};
csMode = {'Accuracy','Precision', 'Recall'};
for iMode = 1
    vrX_plot = unique(floor(vrSnr_unit1/2))*2+1;
    [vr_mu1, vr_sd1] = grpstats(mrAcc_Prc_Rcl_irc2_ks2(:,cviMode{iMode}), floor(vrSnr_unit1/2), {'mean','std'});
    figure_(); vhBar=errorbar(repmat(vrX_plot(:),1,2), vr_mu1, vr_sd1, '-'); grid on;
    vhBar(1).Color='k'; vhBar(2).Color='r';
    ylabel(csMode{iMode}); xlabel('SNR'); axis([0 20 0 1.1]);
    legend({'IronClust','KiloSort2'},'Location','SE');
end

%% Plot accuracy, precision, recall
figure_(); hold on;
vcColor = 'kbr';
iSorter = 5;
vrAcc1 = cell2vec_(cvrAcc_sorter1, iSorter);
hold on; plot(vrSnr_unit1, vrAcc1, 'k.');
for iMode = 1:3
    vrX_plot = unique(floor(vrSnr_unit1/2))*2+1;
    vrPlot1 = mrAcc_Prc_Rcl_irc2_ks2(:,(iMode-1)*2+1);
    [vr_mu1, vr_sd1] = grpstats(vrPlot1, floor(vrSnr_unit1/2), {'mean','std'});
    vhBar=errorbar(vrX_plot+(iMode-2)/4, vr_mu1, vr_sd1, '.'); 
    vhBar(1).Color=vcColor(iMode);
end
ylabel('Score'); xlabel('SNR'); axis([0 20 0 1]); %set(gca,'XTickLabel',(get(gca,'XTick')-1)*2+1);
set(gca,'XTick',1:2:18);
legend({'Accuracy', 'Accuracy', 'Precision', 'Recall'},'Location','SE'); grid on;
title(csSorter{iSorter});

figure_(); hold on; 
plot(vrSnr_unit1, cell2vec_(cvrAcc_sorter1, 2), 'k.');
plot(vrSnr_unit1, cell2vec_(cvrAcc_sorter1, 5), 'r.');
grid on;
axis([0 20 0 1]);

%% siprobe drift accuracy comparison
acc_hybrid_drift_rec_64c_1200s_11 = [0.50002         0.92918     0.51656    0.79644      0.90934      NaN          0.4942          0.68278         0.47172        NaN];
acc_hybrid_static_rec_64c_1200s_11 = [0.88507         0.96705     0.96133    0.94571      0.96724      NaN         0.88433          0.91591         0.50815        NaN   ];
x = mean(acc_hybrid_static_rec_64c_1200s_11,1);
y = mean(acc_hybrid_drift_rec_64c_1200s_11,1);
y(isnan(y)) = 0;
[~,viSorter] = sort(y, 'descend');
figure_(); hBar1=bar([x(viSorter); y(viSorter)]',1);
hBar1(1).FaceColor='k'; hBar1(2).FaceColor='r'; legend({'static', 'drift'});
title_('rec_64c_1200s_11');
set(gca,'XTickLabel', csSorter(viSorter)); xtickangle(45); grid on; colormap([1 1 1;1 0 0]); ylabel('Accuracy');


%% compare tetrode drift comparison
acc_hybrid_static_tetrode = [
    0.55425         0.80825       0.723    0.81375        0.814      NaN         0.91725            0.815         0.56975        NaN   
    0.4255         0.96117     0.72983    0.98617       0.7585      NaN         0.81367          0.94317         0.55283        NaN   
    0.42538           0.962     0.77087    0.98275      0.70462      NaN         0.85488          0.76388         0.61413        NaN   
];
acc_hybrid_drift_tetrode = [
    0.506         0.81375       0.738       0.76      0.71625      NaN           0.736          0.79675           0.704        NaN   
    0.279         0.67675       0.372    0.71125      0.57975      NaN           0.629          0.76525         0.59125        NaN   
    0.333         0.95112     0.58588    0.95375      0.56375      NaN         0.60262          0.59575           0.564        NaN   
];
x = mean(acc_hybrid_static_tetrode,1);
y = mean(acc_hybrid_drift_tetrode,1);
y(isnan(y)) = 0;
[~,viSorter] = sort(y, 'descend');
figure_(); hBar1=bar([x(viSorter); y(viSorter)]',1);
hBar1(1).FaceColor='k'; hBar1(2).FaceColor='r'; legend({'static', 'drift'});
title_('rec_4c_1200s');
set(gca,'XTickLabel', csSorter(viSorter)); xtickangle(45); grid on; colormap([1 1 1;1 0 0]); ylabel('Accuracy');


%% compare runtime. matlab processors used the actual runtime to take out the installation time
run_hybrid_static_siprobe_rec_64c_1200s_11 = [
    264.8         95.6      114.36      175.79      133.42       NaN         788.35           1522.1          1546.9         NaN
    ];
run_hybrid_static_siprobe_rec_32c_1200s_11 = [
    154.3         41.6      49       152.81      163.27       NaN         742.45           1144.9          850.85         NaN   
    150.45        46.9       58.24       154.84      162.83       NaN         725.49           1271.8          960.64         NaN   
    152.39         44.5      54       151.15      141.24       NaN         747.63           1286.7           908.7         NaN   
    ];
run_hybrid_static_siprobe_rec_16c_1200s_11 = [
    86.946         26.7      33.57       108.82      65.02       NaN         297.83           655.02          531.38         NaN   
    90.578         29.6      35.81      109.04      67.54       NaN         302.32           764.01          603.52         NaN   
    90.06         33.3      36.47      109.23      93.93       NaN         291.54           769.77          599.99         NaN   
    ];
x = mean(run_hybrid_static_siprobe_rec_16c_1200s_11,1);
y = mean(run_hybrid_static_siprobe_rec_32c_1200s_11,1);
z = mean(run_hybrid_static_siprobe_rec_64c_1200s_11,1);

z(isnan(z)) = inf; [~,viSorter] = sort(z, 'ascend');
figure_(); hBar1=bar([x(viSorter); y(viSorter); z(viSorter)]',1);
hBar1(1).FaceColor='b'; hBar1(2).FaceColor='r'; hBar1(3).FaceColor='k'; 
title_('run_hybrid_static_siprobe');
set(gca,'XTickLabel', csSorter(viSorter)); xtickangle(45); grid on; ylabel('Runtime (s)');
hold on; hPlot1 = line(get(gca,'XLim'), [1200,1200], 'Color','k');
legend({'16ch', '32ch', '64ch'},'Location','NW');


%% compare runtime. matlab processors used the actual runtime to take out the installation time
run_hybrid_static_siprobe_rec_64c_1200s_11 = [
    264.8         95.6      114.36      175.79      133.42       NaN         788.35           1522.1          1546.9         NaN
    ];
run_hybrid_static_siprobe_rec_64c_600s_11 = [
    154.19         50.9      108      82.66      82.58       NaN         426.99           1273.4          960.73         NaN   
    ];

x = mean(run_hybrid_static_siprobe_rec_64c_600s_11,1);
y = mean(run_hybrid_static_siprobe_rec_64c_1200s_11,1);

y(isnan(y)) = inf; [~,viSorter] = sort(y, 'ascend');
figure_(); hBar1=bar([x(viSorter); y(viSorter)]',1);
hBar1(1).FaceColor='g'; hBar1(2).FaceColor='k'; 
title_('run_hybrid_static_siprobe');
set(gca,'XTickLabel', csSorter(viSorter)); xtickangle(45); grid on; ylabel('Runtime (s)');
hold on; hPlot1 = line(get(gca,'XLim'), [1200,1200], 'Color','k');
hold on; hPlot1 = line(get(gca,'XLim'), [600,600], 'Color','g');
legend({'600s', '1200s'},'Location','NW');


%% plot by technology. do a weighted average by number of clusters
SNR_THRESH = 8;

% collect from json
[tbl_acc0, tbl_prc0, tbl_rcl0, tbl_cpu, csSorter0, vnUnits_acc1] = ...
    sf_json2tbl_study(S_json, SNR_THRESH);

% tbl_acc0 = readtable(fullfile(vcDir, 'average_accuracy_2019-1015.csv')); 
% tbl_prc0 = readtable(fullfile(vcDir, 'average_precision_2019-1015.csv'));
% tbl_rcl0 = readtable(fullfile(vcDir, 'average_recall_2019-1015.csv'));
% tbl_cnt0 = readtable(fullfile(vcDir, 'count_accuracy_2019-1015.csv')); disp('tbl_cnt1'); disp(tbl_cnt0);

title_ = @(x)title(x,'Interpreter','none');
figure_ = @()figure('Color','w');
% csSorter0 = tbl_acc0.Properties.VariableNames(3:end);

% vnUnits = cellfun(@(x)sum(x.trueSnrs>SNR_THRESH), S_json.StudyAnalysisResults)';
% [tbl_acc0.nUnits,tbl_prc0.nUnits, tbl_rcl0.nUnits] = deal(vnUnits);

disp('tbl_acc0'); disp(tbl_acc0); 
disp('tbl_prc0'); disp(tbl_prc0); 
disp('tbl_rcl0'); disp(tbl_rcl0); 

%tbl2mr_ = @(x,vi)table2array(x(vi,3:end));
tbl2mr_ = @(x,vi)x(vi,:);
disp_ = @(x)fprintf('%s\n', sprintf('%0.3f\t', x));
disp_nunits_snr_ = @(x)disp_([nansum(x(:,1:end-1).*x(:,end),1) ./ nansum(x(:,end)), nansum(x(:,end))]);
mean_nunits_snr_ = @(x)[nansum(x(:,1:end-1).*x(:,end),1) ./ nansum(x(:,end)), nansum(x(:,end))];

disp(csSorter0)
vr_acc_paired1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [0:3]+1));
vr_acc_hybrid_siprobe1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [40]+1));
vr_acc_synth_siprobe = mean_nunits_snr_(tbl2mr_(tbl_acc0, [14,28:33]+1));
vr_acc_synth_mea1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [38]+1));

vr_acc_hybrid_tetrode1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [43]));
vr_acc_synth_tetrode1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [17:24,34:37]+1));
vr_acc_manual_tetrode1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [25:27]+1));

vr_acc_paired_monotrode1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [4:7]+1));
vr_acc_synth_monotrode1 = mean_nunits_snr_(tbl2mr_(tbl_acc0, [8:13]+1));


mr_acc1 = [vr_acc_paired1; vr_acc_hybrid_siprobe1; vr_acc_synth_siprobe; vr_acc_synth_mea1; 
          vr_acc_hybrid_tetrode1; vr_acc_synth_tetrode1; vr_acc_manual_tetrode1; 
          vr_acc_paired_monotrode1; vr_acc_synth_monotrode1]; 
% vnUnits_acc1 = mr_acc1(:,end);
cs_dataset_acc1 = {'Paired-assorted', 'SiProbe-hybrid', 'SiProbe-synth', 'MEA-synth', 'Tetrode-hybrid', 'Tetrode-synth', 'Tetrode-manual', 'Monotrode-paired', 'Monotrode-synth'};
% cs_xlabel = arrayfun(@(i)[cs_dataset_acc1{i}, sprintf(' (%d)',vnUnits_acc1(i))], 1:numel(vnUnits_acc1), 'UniformOutput',0);
disp(mr_acc1);

[vrOrder_acc1, viOrder_acc1] = sort(mean_nunits_snr_(mr_acc1), 'descend'); viOrder_acc1(1)=[]; vrOrder_acc1(1)=[];
figure_(); bar(mr_acc1(:,viOrder_acc1),1); grid on; ylabel('Accuracy'); grid on;
set(gca,'XTickLabel', cs_xlabel); xtickangle(45);
legend(csSorter0(viOrder_acc1));


%% collect static SNR for paired recordings



