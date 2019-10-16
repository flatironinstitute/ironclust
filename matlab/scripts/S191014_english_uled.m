% define constants

csDir = {'m90_190915_125303', 'm90_190916_105341', 'm90_190917_112752', 'm90_190918_105939'};
vcFile_xml0 = '/mnt/ceph/users/jjun/DanEnglish/PV_ChR2_Chronic/%s/amplifier.xml';
% 8.9, 6.6, 5.7, 3.4 hours
fPlot = 0;

P = struct('nChans_ain', 8, 'viChans_ain', [6,7], 'nSamples', [], ...
    'nSamples_last', 1e8, 'thresh', 2^15*.75, 'spkLim', [-300, 3300], ...
    'subsample', 1, 'uV_per_bit', .195);
P_filt = struct('sRateHz', [], 'freqLim', [300,6000], 'freqLim_width', [100,1000], 'fGpu', 1, 'fft_thresh', 8);
title_ = @(x)title(x, 'Interpreter', 'none');

% Total duration
% m90_190915_125303
% # stim
%            3        3331        3320
% m90_190917_112752
% # stim
%            0        1062        1060
% m90_190916_105341
% # stim
%                     5266        5240
% m90_190918_105939
% # stim
%            7        3550        3533

% last 1 hour (1e8 samples) duration
% m90_190918_105939
% # stim
%         1495        1496
% m90_190916_105341
% # stim
%         1675        1675
% m90_190917_112752
% # stim
%    999   999
% m90_190918_105939
% # stim
%         1495        1496


%% load and save analysis


% iFile = 3;
% for iFile = 1:numel(csDir)
for iFile = 1:1
    vcDataID1 = csDir{iFile};
    disp(vcDataID1); % show dataID
    vcFile_xml1 = sprintf(vcFile_xml0, csDir{iFile});
    S_xml1 = irc('call', 'load_xml_neuroscope_', {vcFile_xml1});

    
    % load analog in (led control)
    vcFile_ain1 = fullfile(fileparts(vcFile_xml1), 'analogin.dat');
    mnWav_ain1 = irc('call', 'load_bin', {vcFile_ain1, 'uint16', [P.nChans_ain, P.nSamples]});    
    if isempty(P.nSamples_last)
        mnWav_ain1 = mnWav_ain1(P.viChans_ain,:)';
    else
        mnWav_ain1 = mnWav_ain1(P.viChans_ain, [end-P.nSamples_last+1:end])';
    end    
    % disp(std(diff(single(mnWav_ain))));
    if fPlot
        figure; plot(mnWav_ain1(1:100:end,:)); title_(vcDataID1);
    end
    cviTime_on_stim1 = arrayfun(@(i)find(diff(mnWav_ain1(:,i)>=P.thresh)>0)+1, 1:size(mnWav_ain1,2), 'UniformOutput',0);
    cviTime_off_stim1 = arrayfun(@(i)find(diff(mnWav_ain1(:,i)<P.thresh)>0), 1:size(mnWav_ain1,2), 'UniformOutput',0);
    disp('# stim');
    disp(cellfun(@numel, cviTime_on_stim1));
    if 0
        iStim = 1;
        figure; plot(mnWav_ain1(:,iStim)); title_(vcDataID1);
        hold on; plot(cviTime_on_stim1{iStim}, mnWav_ain1(cviTime_on_stim1{iStim},iStim), 'r.');
        hold on; plot(cviTime_off_stim1{iStim}, mnWav_ain1(cviTime_off_stim1{iStim},iStim), 'k.');
    end
    clear mnWav_ain


    % load waveform recording and analyze last hour
    mnWav_T1 = irc('call', 'load_bin', {S_xml1.vcFile_dat, 'int16', [S_xml1.nChans, P.nSamples]});
    cviSite2chan = {S_xml1.viSite2chan1, S_xml1.viSite2chan2, S_xml1.viSite2chan3, S_xml1.viSite2chan4};
    P_filt.sRateHz = S_xml1.sRateHz/P.subsample;
    if ~isempty(P.nSamples_last)
        mnWav_T1 = mnWav_T1(:, [end-P.nSamples_last+1:end]);
    end

    vS_shank1 = cell(size(cviSite2chan)); % to be converted to vector of struct
    for iShank = 1:numel(cviSite2chan)
        mrWav1 = single(mnWav_T1(cviSite2chan{iShank},1:P.subsample:end)'); % temporal downlsample
    %     mrWav1 = mrWav1 - mean(mrWav1,2);
        ctnWav1 = cellfun(@(x)irc('call', 'mr2tr', {mrWav1, P.spkLim, x/P.subsample}), cviTime_on_stim1, 'UniformOutput', 0);
        cmrWav_raw_mu1 = cellfun(@(x)squeeze(mean(x,2)), ctnWav1, 'UniformOutput', 0);
        cmrWav_raw_sd1 = cellfun(@(x)squeeze(std(x,[],2)), ctnWav1, 'UniformOutput', 0);

        % filter
        mrWav1 = fft_filter(mrWav1, P_filt); % filter the shank    
        ctnWav1 = cellfun(@(x)irc('call', 'mr2tr', {mrWav1, P.spkLim, x/P.subsample}), cviTime_on_stim1, 'UniformOutput', 0);
        cmrWav_mu1 = cellfun(@(x)squeeze(mean(x,2)), ctnWav1, 'UniformOutput', 0);
        cmrWav_sd1 = cellfun(@(x)squeeze(std(x,[],2)), ctnWav1, 'UniformOutput', 0);

        vS_shank1{iShank} = struct('cviTime_on_stim1', cviTime_on_stim1, 'cviTime_off_stim1', cviTime_off_stim1, ...
            'cmrWav_mu1', cmrWav_mu1, 'cmrWav_sd1', cmrWav_sd1, ...
            'cmrWav_raw_mu1', cmrWav_raw_mu1, 'cmrWav_raw_sd1', cmrWav_raw_sd1);
        fprintf('.');
    end
    clear mrWav1 mnWav_T1
    vS_shank1 = cell2mat(vS_shank1);
    save(fullfile(fileparts(vcFile_xml1), 'analysis.mat'), 'vS_shank1', 'S_xml1', 'P_filt', 'P', 'vcDataID1');
end %for iFile


%% plot
% iFile = 1;
% vcFile_xml1 = sprintf(vcFile_xml0, csDir{iFile});
% load(fullfile(fileparts(vcFile_xml1), 'analysis.mat')

trWav_mu1_shank = {vS_shank1.cmrWav_mu1}; trWav_mu1_shank = cat(3,trWav_mu1_shank{:});
trWav_sd1_shank = {vS_shank1.cmrWav_sd1}; trWav_sd1_shank = cat(3,trWav_sd1_shank{:});
trWav_raw_mu1_shank = {vS_shank1.cmrWav_raw_mu1}; trWav_raw_mu1_shank = cat(3,trWav_raw_mu1_shank{:});
trWav_raw_sd1_shank = {vS_shank1.cmrWav_raw_sd1}; trWav_raw_sd1_shank = cat(3,trWav_raw_sd1_shank{:});

for iMode_plot = 1:2
    figure('Color','w','Name', vcDataID1); 
    ax=[];
    vrT_plot = (P.spkLim(1):P.spkLim(2))/S_xml1.sRateHz*1000;
    for iPlot=1:8
        ax(iPlot)=subplot(4,2,iPlot);
        switch iMode_plot
            case 2
                plot(vrT_plot, trWav_sd1_shank(:,:,iPlot) * P.uV_per_bit);
                ylabel('Filtered SD (uV)'); ylim([0 50]);
            case 1
                plot(vrT_plot, trWav_raw_mu1_shank(:,:,iPlot) * P.uV_per_bit);
                ylabel('Raw mean (uV)'); ylim([-1e4, 1e4]);
        end
        title(sprintf('LED%d, Shank%d', mod(iPlot-1,2)+1, ceil(iPlot/2)));        
        grid on; xlabel('Time since LED pulse (ms)');    
    end
    linkaxes(ax,'xy');
    set(gca,'XLim', vrT_plot([1,end])); grid on;
end
