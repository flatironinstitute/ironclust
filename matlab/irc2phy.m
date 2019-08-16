function irc2phy(vcFile_prm, vcDir_out)
% J. James Jun 2019 Aug 16
% Code cited from
%   https://github.com/JaneliaSciComp/JRCLUST/blob/9e77b8422c9a16dc65231c60977ba6af2b52fe91/%2Bjrclust/%2Bexport/phy.m
% phy data format:
%   https://phy.readthedocs.io/en/latest/visualization/

% Input
% ----
% EXPORT_MODE: 1 for loading all recordings to memory, 2 for paged read

if nargin<2, vcDir_out=''; end
if isempty(vcDir_out)
    vcDir_out = fullfile(fileparts(vcFile_prm), 'phy'); 
end
mkdir_(vcDir_out);

[S0, P] = irc('call', 'load_cached_', {vcFile_prm});

[vnAmp_spk, viTime_spk, viSite_spk, S_clu] = deal(S0.vrAmp_spk, S0.viTime_spk, S0.viSite_spk, get_(S0, 'S_clu'));
if ~isempty(S_clu)
    viClu_spk = S_clu.viClu;
    viClu_spk(viClu_spk<0) = 0;
end

writeNPY_ = @(x,y)writeNPY(x, fullfile(vcDir_out, y));

writeNPY_(uint64(vnAmp_spk), 'amplitudes.npy');
writeNPY_(uint64(viTime_spk), 'spike_times.npy');
writeNPY_(int32(viSite_spk) - 1, 'spike_sites.npy');
writeNPY_(int32(P.viSite2Chan) - 1, 'channel_map.npy');
writeNPY_(P.mrSiteXY, 'channel_positions.npy'); % dimension?
if ~isempty(S_clu)
    writeNPY_(uint32(viClu_spk) - 1, 'spike_clusters.npy'); % -1 for zero indexing
end

% read feature file and write to it
writeNPY_(pcFeatures, 'pc_features.npy');

    [nSites, nSpikes] = deal(res.featuresShape(1), res.featuresShape(3));

    % take just primary peak
    spikeFeatures = squeeze(spikeFeatures(:, 1, :))'; % nSpikes x nFeatures
    if hCfg.nPCsPerSite == 1
        pcFeatures = reshape(spikeFeatures, size(spikeFeatures, 1), 1, []);
    else
        nSites = nSites/hCfg.nPCsPerSite;
        pcFeatures = zeros(nSpikes, hCfg.nPCsPerSite, nSites, 'single');
        for i = 1:hCfg.nPCsPerSite
            pcFeatures(:, i, :) = spikeFeatures(:, ((i-1)*nSites+1):i*nSites);
        end
    end

    featuresFile = fullfile(hCfg.outputDir, 'pc_features.npy');
    writeNPY(pcFeatures, featuresFile);
    hCfg.updateLog('phy', sprintf('Saved spikeFeatures to %s', featuresFile), 0, 0);

    indFile = fullfile(hCfg.outputDir, 'pc_feature_ind.npy');
    writeNPY(uint32(hCfg.siteNeighbors(1:nSites, res.spikeSites)') - 1, indFile); % -1 for zero indexing
    hCfg.updateLog('phy', sprintf('Saved spikeFeature indices to %s', indFile), 0, 0);

% param file
if exist(fullfile(hCfg.outputDir, 'params.py'), 'file') ~= 2
    paramFile = fullfile(hCfg.outputDir, 'params.py');
    fid = fopen(paramFile, 'w');

    rawRecordings = cellfun(@(x) sprintf('r''%s''', x), hCfg.rawRecordings, 'UniformOutput', 0);
    rawRecordings = ['[' strjoin(rawRecordings, ', ') ']'];
    fprintf(fid, 'dat_path = %s\n', rawRecordings);

    fprintf(fid, 'n_channels_dat = %i\n', hCfg.nChans);

    fprintf(fid, 'dtype = ''%s''\n', dtype2NPY(hCfg.dataType));

    fprintf(fid, 'offset = %d\n', hCfg.headerOffset);

    if hCfg.sampleRate ~= floor(hCfg.sampleRate)
        fprintf(fid,'sample_rate = %i\n', hCfg.sampleRate);
    else
        fprintf(fid,'sample_rate = %i.\n', hCfg.sampleRate);
    end

    fprintf(fid,'hp_filtered = False');
    fclose(fid);

    hCfg.updateLog('phy', sprintf('Saved params to %s', paramFile), 0, 0);
end
end % func


%% LOCAL FUNCTIONS
function dtype = dtype2NPY(dtype)
    switch dtype
        case 'single'
            dtype = 'float32';
        case 'double'
            dtype = 'float64';
    end
end

