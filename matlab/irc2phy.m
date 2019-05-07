%-------------------------------------------------------------------------
% Converts irc format to phy format and calls 
% assume all sites are used (mda format)
function irc2phy(vcFile_prm, savePath)

% Load irc output
P = loadParam_(vcFile_prm);
S0 = load(strrep(vcFile_prm, '.prm', '_jrc.mat'));
S_clu = S0.S_clu;

% rez struct (KiloSort)
%           ops: [1×1 struct]
%            xc: [1×128 double]
%            yc: [1×128 double]
%       xcoords: [1×128 double]
%       ycoords: [1×128 double]
%     connected: [128×1 logical]
%          Wrot: [128×128 double]
%          temp: [1×1 struct]
%        errall: [1×39 double]
%           dWU: [32×128×1024 single]
%       nspikes: [1024×2639 double]
%           st3: [2301962×4 double]
%         cProj: [2301962×16 single]
%      simScore: [1024×1024 single]
%        iNeigh: [16×1024 double]
%       cProjPC: [2301962×3×16 single]
%      iNeighPC: [16×1024 double]
%             W: [44×1024×3 single]
%             U: [128×1024×3 single]
%            mu: [1024×1 single]
%           t2p: [1024×2 double]
%         nbins: [1025×1 double]
%          ypos: [1×1024 double]
%       WrotInv: [128×128 double]
%          Wraw: [128×44×1024 single]

% https://github.com/kwikteam/phy-contrib/blob/master/docs/template-gui.md
% {st3(:,1:5), connected, xcoords, ycoords, W, U, cProj, iNeigh, cProjPC, iNeighPC, Wrot, simScore, ops:{Nchan, chanMap, Nfilt, fbinary, NchanTOT, fs}}
[nSpikes, nSites, nClu] = deal(numel(S0.viTime_spk), numel(P.viSite2Chan), S_clu.nClu);
[xcoords, ycoords, kcoords] = deal(P.mrSiteXY(:,1), P.mrSiteXY(:,2), P.viShank_site);
connected = true(1, nSites);

st3 = zeros(nSpikes, 4, 'double'); 
st3(:,1) = S0.viTime_spk;
st3(:,2) = S0.S_clu.viClu;
st3(:,3) = abs(S0.vrAmp_spk);

% compute Wrot by whitening the filtered traces

% cProjPC
% pc_features.npy - [nSpikes, nFeaturesPerChannel, nPCFeatures] single

% iNeighPC
% pc_feature_ind.npy - [nTemplates, nPCFeatures] uint32 matrix specifying which pcFeatures are included in the pc_features matrix.

% build a cluster template. use code from zack's waveform visualizer
templates = permute(S_clu.tmrWav_spk_clu, [3,1,2]); % make sure it spans full channels, don't be confused with trWav_clu which stores a subset

rez = makeStruct_(st3, connected, xcoords, ycoords, kcoords, templates, cProj, iNeigh, cProjPC, iNeighPC);
rez.ops = struct('Nchan', nSites, 'chanMap', P.viSite2Chan, 'Nfilt', nClu, 'fbinary', P.vcFile, 'NchanTOT', P.nChans, 'fs', P.sRateHz);

mkdir_(savePath);
rezToPhy(rez, savePath);

end %func


%==========================================================================
function out1 = makeStruct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = loadParam_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mkdir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

%-------------------------------------------------------------------------
% copied from KiloSort/finalPass/rezToPhy.m
% 4/30/2019
function [spikeTimes, clusterIDs, amplitudes, templates, templateFeatures, ...
    templateFeatureInds, pcFeatures, pcFeatureInds] = rezToPhy(rez, savePath)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)
%
% spikeTimes will be in samples, not seconds
[fSave_template, fSave_whitening] = deal(0, 1); % set constants

outputs = {'amplitudes.npy', 'channel_map.npy', 'channel_positions.npy', 'pc_features.npy', ...
           'pc_feature_ind.npy', 'similar_templates.npy', 'spike_clusters.npy', 'spike_templates.npy', ...
           'spike_times.npy', 'templates.npy', 'templates_ind.npy', 'template_features.npy', ...
           'template_feature_ind.npy', 'whitening_mat.npy', 'whitening_mat_inv.npy'};

fs = dir(fullfile(savePath, '*.npy'));
for i = 1:length(fs)
    fname = fs(i).name;
    % don't delete .npy files which have nothing to do with us
    if find(strcmp(fname, outputs))
        delete(fullfile(savePath, fname));
    end
end
if exist(fullfile(savePath, '.phy'), 'dir')
    rmdir(fullfile(savePath, '.phy'), 's');
end

spikeTimes = uint64(rez.st3(:,1));
% [spikeTimes, ii] = sort(spikeTimes);
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4
    spikeClusters = uint32(1+rez.st3(:,5));
end
amplitudes = rez.st3(:,3);

Nchan = rez.ops.Nchan;

% try
%     load(rez.ops.chanMap);
% catch
%    chanMap0ind  = [0:Nchan-1]';
%    connected    = ones(Nchan, 1);
%    xcoords      = ones(Nchan, 1);
%    ycoords      = (1:Nchan)';
% end
% chanMap0 = chanMap(connected>1e-6);

connected   = rez.connected(:);
xcoords     = rez.xcoords(:);
ycoords     = rez.ycoords(:);
kcoords     = rez.kcoords(:);
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;

if isfield(rez, 'templates')
    templates = rez.templates;
else
    nt0 = size(rez.W,1);
    U = rez.U;
    W = rez.W;

    % for i = 1:length(chanMap0)
    %     chanMap0(i) = chanMap0(i) - sum(chanMap0(i) > chanMap(connected<1e-6));
    % end
    % [~, invchanMap0] = sort(chanMap0);

    templates = zeros(Nchan, nt0, rez.ops.Nfilt, 'single');
    for iNN = 1:rez.ops.Nfilt
       templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
    end
    templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
end %switch
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatures = rez.cProj;
templateFeatureInds = uint32(rez.iNeigh);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

if ~isempty(savePath)
    
    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
    writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    if size(rez.st3,2)>4
        writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
    writeNPY(templates, fullfile(savePath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
%     Fs = rez.ops.fs;
    conn        = logical(connected);
    chanMap0ind = int32(chanMap0ind);
    
    writeNPY(chanMap0ind(conn), fullfile(savePath, 'channel_map.npy'));
    %writeNPY(connected, fullfile(savePath, 'connected.npy'));
%     writeNPY(Fs, fullfile(savePath, 'Fs.npy'));
    writeNPY([xcoords(conn) ycoords(conn)], fullfile(savePath, 'channel_positions.npy'));
    writeNPY(kcoords(conn), fullfile(savePath, 'channel_shanks.npy'));
    
    if fSave_template
        writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
        writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
    end
    writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
    writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    if fSave_whitening % just store an identity matrix
        whiteningMatrix = eye(Nchan);
%         whiteningMatrix = rez.Wrot/200;
        whiteningMatrixInv = whiteningMatrix^-1;
        writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
        writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
    end
    
%     if isfield(rez, 'simScore')
%         similarTemplates = rez.simScore;
%         writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
%     end
    
     %make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        fid = fopen(fullfile(savePath,'params.py'), 'w');
        
        [~, fname, ext] = fileparts(rez.ops.fbinary);
        
        fprintf(fid,['dat_path = ''',fname ext '''\n']);
        fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        if mod(rez.ops.fs,1)
            fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
        else
            fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
        end
        fprintf(fid,'hp_filtered = False');
        fclose(fid);
    end
end
end %func