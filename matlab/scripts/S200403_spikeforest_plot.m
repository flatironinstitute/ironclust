addpath jsonlab-1.5
switch 3
    case 1 % Flatiron
        vcDir = '/mnt/home/jjun/src/spikeforest2/working/runs/';
        vcFile_json = '2020_04_14a/output.json';
    case 2 % linux
        vcDir = '/home/jjun/src/ironclust/matlab/scripts/';
        vcFile_json = 'spikeforest2_output.json';
    case 3 % MAC
        vcDir = '/Users/jamesjun/src/ironclust/matlab/scripts/';
        vcFile_json = 'spikeforest2_output.json';
end
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

%% create a struct

vS_sar = cell2mat(S_json.StudyAnalysisResults);
csStudySetName = unique({vS_sar.studySetName});
csStudySetName_sar = {vS_sar.studySetName};
csStudyName_sar = {vS_sar.studyName};
S_result = struct();
csFieldName = {'accuracies', 'precisions', 'recalls', 'cpuTimesSec'};
for iStudySet = 1:numel(csStudySetName)
    vcStudySetName = csStudySetName{iStudySet};
    S_studyset1 = struct();
    vi_sar1 = find(strcmpi(csStudySetName_sar, vcStudySetName));
    for iStudy = 1:numel(vi_sar1)
        S_sar1 = vS_sar(vi_sar1(iStudy));
        cS_sorter1 = S_sar1.sortingResults;
        vrSnr = S_sar1.trueSnrs;
        viRec = S_sar1.trueRecordingIndices;
        if isempty(vrSnr), continue; end
        S_study1 = struct();
        for iSorter = 1:numel(cS_sorter1)
            try
                S_ = cS_sorter1{iSorter};
                S1 = struct();
                S1.viRec = viRec;
                S1.accuracies = S_.accuracies;
                S1.vrSnr = vrSnr;
                S1.cpuTimesSec = S_.cpuTimesSec;
                S_study1.(S_.sorterName) = S1;
            catch
                S_study1.(S_.sorterName) = [];
            end
        end
        S_studyset1.(S_sar1.studyName) = S_study1;
    end
    S_result.(vcStudySetName) = S_studyset1;
end

%% select study and plot accuracy vs length

vcStudySetName = 'LONG_STATIC';
%vcStudySetName = 'LONG_DRIFT';

vcMetric = 'accuracy'; % select from {'accuracy', 'count', 'cpuTimesSec'}
P = struct('snr_thresh', 8, 'accuracy_thresh', .8);

vrTimeDur = [300, 600, 1200, 2400, 4800];
vrChan = [8,16];
csSorterName = {'HerdingSpikes2', 'IronClust', 'JRClust', 'KiloSort', 'KiloSort2', 'Klusta', 'MountainSort4', 'SpykingCircus', 'Tridesclous'};
csMetric = {'accuracy_mean', 'count_accuracy', 'cpuTimesSec_mean'};

study_name_ = @(t,c)sprintf('%s_%ds_%dc', vcStudySetName, t, c);
[trMean_time_chan_sorter, trSEM_time_chan_sorter] = ...
    deal(nan(numel(vrTimeDur), numel(vrChan), numel(csSorterName)));
nansem_ = @(x)nanstd(x)/sqrt(numel(x));
for iTimeDur = 1:numel(vrTimeDur)
    timeDur = vrTimeDur(iTimeDur);
    for iChan = 1:numel(vrChan)
        nChan = vrChan(iChan);
        S_study = S_result.(vcStudySetName).(study_name_(timeDur, nChan));
        for iSorter = 1:numel(csSorterName)
            try
                S_ = S_study.(csSorterName{iSorter});
                viRec_ = S_.viRec;
                switch vcMetric
                    case 'accuracy'
                        viRec_(S_.vrSnr<P.snr_thresh) = nan;
                        means_ = grpstats(S_.accuracies, viRec_);
                        mean_ = nanmean(means_);
                        sem_ = nansem_(means_);
                    case 'count'
                        viRec_(S_.accuracies<P.accuracy_thresh) = nan;
                        mean_ = grpstats(S_.accuracies, viRec_, {'count'});
                        sem_ = nan;
                    case 'cpuTimesSec'
                        means_ = grpstats(S_.cpuTimesSec, viRec_);
                        mean_ = nanmean(means_);
                        sem_ = nansem_(means_);
                    otherwise
                        error('invalid metric');
                end
                trMean_time_chan_sorter(iTimeDur, iChan, iSorter) = mean_;
                trSEM_time_chan_sorter(iTimeDur, iChan, iSorter) = sem_;
            catch
                ; 
            end
        end
    end
end

% plot
nTime = size(trMean_time_chan_sorter,1);
mrColor = flipud([1,1,1] .* linspace(.25,.75,nTime)');
for iPlot=2
    iChan = iPlot;
    figure('Color','w');
    title(sprintf('%s-%dch', vcStudySetName, vrChan(iChan)),'Interpreter','none');    
    mrMean_dur_sorter = squeeze(trMean_time_chan_sorter(:,iChan,:));
    mrSEM_dur_sorter = squeeze(trSEM_time_chan_sorter(:,iChan,:));
    if false
        [~, viSorter] = sort(nanmean(mrMean_dur_sorter), 'ascend');
    else
        viSorter = 1:numel(csSorterName);
    end
    %subplot(2,1,iPlot); 
    hold on;
    vhBar = bar(mrMean_dur_sorter(:,viSorter)',1); drawnow;
    cmrErrX = arrayfun(@(x)x.getSingleBarExtentsArray(1), vhBar, 'UniformOutput',0);
    mrErrX = cell2mat(cellfun(@(x)x(2,:)', cmrErrX, 'UniformOutput',0));
    
    set(gca,'XTickLabel', csSorterName, 'XTick', viSorter); xtickangle(gca, -20);
    ylabel('Accuracy (SNR>8)');
    mrZero_dur_sorter = zeros(size(mrMean_dur_sorter));
    errorbar_ = @(c,yd,yu)errorbar(mrErrX, mrMean_dur_sorter(:,viSorter)', yd, yu, ...
        'LineStyle','none','Color',c,'CapSize',0,'LineWidth',1);
%     switch iPlot
%         case 1
%             arrayfun(@(x)set(x,'FaceColor','none','EdgeColor','k','LineWidth',1), vhBar);
%             errorbar_('k', mrSEM_dur_sorter(:,viSorter)', mrZero_dur_sorter');
%         case 2
            arrayfun(@(x,y)set(x,'EdgeColor','none', 'FaceColor', mrColor(y,:)), vhBar, 1:nTime);
            errorbar_('k', mrZero_dur_sorter', mrSEM_dur_sorter(:,viSorter)');
%     end
    legend(arrayfun(@(x)sprintf('%dmin',x), vrTimeDur/60, 'UniformOutput',0));
    set(gca,'YTick',0:.1:1,'YGrid','on');
    axis([.5, 11 0 1]);
end

