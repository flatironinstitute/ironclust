function irc2phy(vcFile_prm, vcDir_out)
% J. James Jun 2019 Aug 16
% Code cited from
%   https://github.com/JaneliaSciComp/JRCLUST/blob/9e77b8422c9a16dc65231c60977ba6af2b52fe91/%2Bjrclust/%2Bexport/phy.m
%   https://github.com/petersenpeter/KilosortWrapper/blob/master/rezToPhy_KSW.m
% phy data format:
%   https://github.com/kwikteam/phy-contrib/blob/master/docs/template-gui.md
%   https://phy.readthedocs.io/en/latest/visualization/
% [10/28/2019]
% - Supports irc2 format

% Input
% ----
% EXPORT_MODE: 1 for loading all recordings to memory, 2 for paged read

if nargin<2, vcDir_out=''; end
if isempty(vcDir_out)
    vcDir_out = fullfile(fileparts(vcFile_prm), 'phy'); 
end
mkdir_(vcDir_out);

vcFile_mat = strrep(vcFile_prm, '.prm', '_irc.mat');
if exist_file_(vcFile_mat)
    S0 = load(vcFile_mat);
    P = S0.P;
else
    [S0, P] = load_cached_(vcFile_prm);
end

[vrAmp_spk, viTime_spk, viSite_spk, S_clu] = deal(S0.vrAmp_spk, S0.viTime_spk, S0.viSite_spk, get_(S0, 'S_clu'));
nSpikes = numel(viTime_spk);
if ~isempty(S_clu)
    viClu_spk = S_clu.viClu;
    viClu_spk(viClu_spk<0) = 0;
end
nSites = numel(P.viSite2Chan);

writeNPY_(uint64(abs(vrAmp_spk)), fullfile(vcDir_out, 'amplitudes.npy'));
writeNPY_(uint64(viTime_spk), fullfile(vcDir_out, 'spike_times.npy'));
writeNPY_(int32(viSite_spk) - 1, fullfile(vcDir_out, 'spike_sites.npy'));
writeNPY_(int32(P.viSite2Chan) - 1, fullfile(vcDir_out, 'channel_map.npy'));
writeNPY_(P.mrSiteXY, fullfile(vcDir_out, 'channel_positions.npy')); % dimension?
if ~isempty(S_clu)
    writeNPY_(uint32(viClu_spk)-1, fullfile(vcDir_out, 'spike_templates.npy'));
    writeNPY_(single(S_clu.mrWavCor), fullfile(vcDir_out, 'similar_templates.npy'));
end

% read feature file and write to it
vcFile_fet = strrep(vcFile_prm, '.prm', '_fet.irc');
if exist_file_(vcFile_fet) % irc vers2 format    
    trFet_spk = load_bin_(vcFile_fet, S0.type_fet, S0.dimm_fet);
    nSites_fet = size(trFet_spk, 2);
    writeNPY_(permute(trFet_spk, [3,1,2]), fullfile(vcDir_out, 'pc_features.npy'));    
else
    % irc ver1 format
    trFet_spk = load_spkfet_(S0, P);
    trFet_spk = reshape(trFet_spk(:,1,:), [], P.nPcPerChan, nSpikes);
    nSites_fet = size(trFet_spk, 1);
    writeNPY_(permute(trFet_spk, [3,2,1]), fullfile(vcDir_out, 'pc_features.npy'));
end

% write locations of features
writeNPY_(uint32(P.miSites(1:nSites_fet, S_clu.viSite_clu)') - 1, fullfile(vcDir_out, 'pc_feature_ind.npy')); % -1 for zero indexing
    
% Templates file
switch 1
    case 2
        [~, nSites, nClu] = size(S_clu.tmrWav_spk_clu);
        writeNPY_(permute(S_clu.tmrWav_spk_clu, [3,1,2]), fullfile(vcDir_out, 'templates.npy'));
        writeNPY_(repmat(0:nSites, nClu, 1), fullfile(vcDir_out, 'templates_ind.npy'));        
    case 1
        if isfield(S_clu, 'trWav_spk_clu')
            trWav_clu = S_clu.trWav_spk_clu(:,1:nSites_fet,:);
        elseif isfield(S_clu, 'trWav_clu')
            trWav_clu = S_clu.trWav_clu(:,1:nSites_fet,:);
        else
            error('irc2phy: `trWav_spk_clu` or `trWav_clu` not found');
        end
        writeNPY_(permute(trWav_clu, [3,1,2]), fullfile(vcDir_out, 'templates.npy'));
        writeNPY_(uint32(P.miSites(1:nSites_fet, S_clu.viSite_clu)') - 1, fullfile(vcDir_out, 'template_ind.npy'));
end %switch
writeNPY_(eye(nSites), fullfile(vcDir_out, 'whitening_mat.npy'));
writeNPY_(eye(nSites), fullfile(vcDir_out, 'whitening_mat_inv.npy'));

% param file
write_params_(vcDir_out, P);
end %func


%--------------------------------------------------------------------------
function write_params_(vcDir_out, P)
fOverwrite = 1;

vcFile_params = fullfile(vcDir_out, 'params.py');
if ~exist_file_(vcFile_params) || fOverwrite
    fid = fopen(vcFile_params, 'w');
%     csFile_merge = get_raw_files_(P);
    [vcDir, vcFile, vcExt] = fileparts(P.vcFile);
    if isempty(vcDir), vcDir = pwd(); end
    rawRecordings = sprintf('''%s''', fullfile(vcDir, [vcFile, vcExt]));
%     rawRecordings = cellfun(@(x) sprintf('r''%s''', x), csFile_merge, 'UniformOutput', 0);
%     rawRecordings = ['[' strjoin(rawRecordings, ', ') ']'];
    fprintf(fid, 'dat_path = %s\n', rawRecordings);
    fprintf(fid, 'n_channels_dat = %i\n', P.nChans);
    fprintf(fid, 'dtype = ''%s''\n', dtype2NPY_(P.vcDataType));
    fprintf(fid, 'offset = %d\n', P.header_offset);
    if P.sRateHz == floor(P.sRateHz)
        fprintf(fid,'sample_rate = %i\n', P.sRateHz);
    else
        fprintf(fid,'sample_rate = %i.\n', P.sRateHz);
    end
    fprintf(fid,'hp_filtered = False');
    fclose(fid);
    fprintf('Wrote to %s\n', vcFile_params);
end
end % func


%--------------------------------------------------------------------------
function y = writeNPY_(x, vcFile)
try
    if ~exist('writeNPY', 'file')
        addpath('npy-matlab');
    end
    writeNPY(x, vcFile);
    fprintf('Wrote to %s\n', vcFile);
catch
    fprintf(2, 'Error writing to %s\n', vcFile);
end
end %func


%--------------------------------------------------------------------------
function dtype = dtype2NPY_(dtype)
    switch dtype
        case 'single'
            dtype = 'float32';
        case 'double'
            dtype = 'float64';
    end
end % func


%--------------------------------------------------------------------------
function out1 = exist_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mkdir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_raw_files_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_spkfet_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_bin_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

function [out1, out2] = load_cached_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end


