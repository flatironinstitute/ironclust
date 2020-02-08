function vcFile_params = irc2phy(vcFile_prm, vcDir_out)
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
assert(exist_file_(vcFile_prm), ['File does not exist: ', vcFile_prm]);
if isempty(vcDir_out)
    vcDir_out = fullfile(fileparts(vcFile_prm), 'phy'); 
end
mkdir_(vcDir_out);
delete(fullfile(vcDir_out, '*'));

[S0, S_auto, trPc_spk, trPc2_spk, P] = load_irc2_(vcFile_prm);

nSites = numel(P.viSite2Chan);

writeNPY_(uint64(abs(S0.vrAmp_spk)), fullfile(vcDir_out, 'amplitudes.npy'));
writeNPY_(uint64(S0.viTime_spk), fullfile(vcDir_out, 'spike_times.npy'));
writeNPY_(int32(S0.viSite_spk) - 1, fullfile(vcDir_out, 'spike_sites.npy'));
writeNPY_(int32(P.viSite2Chan) - 1, fullfile(vcDir_out, 'channel_map.npy'));
writeNPY_(double(P.mrSiteXY), fullfile(vcDir_out, 'channel_positions.npy')); % dimension?
writeNPY_(uint32(S_auto.viClu)-1, fullfile(vcDir_out, 'spike_templates.npy'));

% read feature file and write to it 
[nPc, nSites_fet, ~] = size(trPc_spk);
writeNPY_(permute(trPc_spk(1:min(3,nPc),:,:), [3,1,2]), fullfile(vcDir_out, 'pc_features.npy'));    
    
% write locations of features
writeNPY_(uint32(P.miSites(1:nSites_fet, S_auto.viSite_clu)') - 1, ...
    fullfile(vcDir_out, 'pc_feature_ind.npy')); % -1 for zero indexing
    
% Templates file: compute telmpates using most popular sites and compute
% waveform similarity
[trWav_clu, viSite_clu, mrWavCor] = clu_wav_(S0, S_auto, trPc_spk, trPc2_spk);
writeNPY_(permute(trWav_clu, [3,1,2]), fullfile(vcDir_out, 'templates.npy'));
writeNPY_(uint32(P.miSites(1:nSites_fet, viSite_clu)') - 1, fullfile(vcDir_out, 'template_ind.npy'));
writeNPY_(single(mrWavCor), fullfile(vcDir_out, 'similar_templates.npy'));

writeNPY_(eye(nSites), fullfile(vcDir_out, 'whitening_mat.npy'));
writeNPY_(eye(nSites), fullfile(vcDir_out, 'whitening_mat_inv.npy'));

% param file
vcFile_params = write_params_(vcDir_out, P);
end %func


%--------------------------------------------------------------------------
function vcFile_params = write_params_(vcDir_out, P)
fOverwrite = 1;

vcFile_params = fullfile(vcDir_out, 'params.py');
if ~exist_file_(vcFile_params) || fOverwrite
    fid = fopen(vcFile_params, 'w');
    [vcDir, vcFile, vcExt] = fileparts(P.vcFile);
    if isempty(vcDir), vcDir = pwd(); end
    switch lower(vcExt)
        case '.bin'
            rawRecordings = fullfile(vcDir, [vcFile, vcExt]);
        case '.mda'
            rawRecordings = fullfile(vcDir, [vcFile, '.bin']);
            mda2bin_(P.vcFile, rawRecordings);
        otherwise, error('Unsupported raw format: %s', P.vcFile);
    end
    if ispc()
        rawRecordings = strrep(rawRecordings, filesep(), [filesep(),filesep()]);
    end
    fprintf(fid, 'dat_path = ''%s''\n', rawRecordings);
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
    if exist_file_(vcFile), delete(vcFile); end
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
% Call from irc.m
function cell_out = call_irc_(dbstack1, cell_input, nargout)
vcFunc = dbstack1(1).name;
try
    switch nargout
        case 0, cell_out{1} = []; irc('call', vcFunc, cell_input);
        case 1, cell_out{1} = irc('call', vcFunc, cell_input);
        case 2, [cell_out{1}, cell_out{2}] = irc('call', vcFunc, cell_input);
        case 3, [cell_out{1}, cell_out{2}, cell_out{3}] = irc('call', vcFunc, cell_input);
        case 4, [cell_out{1}, cell_out{2}, cell_out{3}, cell_out{4}] = irc('call', vcFunc, cell_input);
        otherwise, error('call_irc_: undefined func: %s', vcFunc);
    end
catch ME
    fprintf(2, 'call_irc_: %s\n', ME.message);
    rethrow ME;
end
end %func


%--------------------------------------------------------------------------
% Call from irc2.m
function cout = call_irc2_(dbstack1, cin, nargout)
vcFunc = dbstack1(1).name;
try
    switch nargout
        case 0, cout{1} = []; irc2('call', vcFunc, cin);
        case 1, cout{1} = irc2('call', vcFunc, cin);
        case 2, [cout{1}, cout{2}] = irc2('call', vcFunc, cin);
        case 3, [cout{1}, cout{2}, cout{3}] = irc2('call', vcFunc, cin);
        case 4, [cout{1}, cout{2}, cout{3}, cout{4}] = irc2('call', vcFunc, cin);
        case 5, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}] = irc2('call', vcFunc, cin);
        case 6, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}, cout{6}] = irc2('call', vcFunc, cin);
        case 7, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}, cout{6}, cout{7}] = irc2('call', vcFunc, cin);
        case 8, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}, cout{6}, cout{7}, cout{8}] = irc2('call', vcFunc, cin);
        otherwise, error('call_irc2_: too many output');
    end
catch ME
    fprintf(2, 'call_irc2_: %s\n', ME.message);
    rethrow ME;
end
end %func


function varargout = load_fet_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = exist_file_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = mkdir_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = clu_wav_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = mda2bin_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_irc2_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
