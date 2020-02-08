function [hCfg, res] = irc(prmFile)
% [hCfg, res] = irc(filename)
% import ironclust results
NUM_PC = 3;

% delete previous 
% delete(fullfile(vcDir_out, '*'), 'f');
[hCfg, res] = deal([]);
[vcDir_,vcFile_,vcExt_] = fileparts(prmFile);
filename_ = fullfile(vcDir_, vcFile_);

try
    vcFile_prm = strrep(prmFile, '_jrclust.prm', '.prm');
    [S0, S_auto, trPc_spk, trPc2_spk, P] = irc2('call', 'load_irc2_', {vcFile_prm});
    [trWav_clu, viSite_clu, mrWavCor] = irc2('call', 'clu_wav_', {S0, S_auto, trPc_spk, trPc2_spk});
    S_clu = S0.S_clu;
catch ME
    error('Failed to load ''%s'': %s', prmFile, ME.message);
end

% new parameter file prmFile
S0.P.vcFile_prm = prmFile;
hCfg = jrclust.Config(prmFile);

% try to import params directly
oldPrms = fieldnames(S0.P);

% handle special cases in P
if isfield(S0.P, 'spkLim_raw_ms') && isempty(S0.P.spkLim_raw_ms)
    if isfield(S0.P, 'spkLim_ms') && isfield(S0.P, 'spkLim_raw_factor')
        S0.P.spkLim_raw_ms = S0.P.spkLim_ms * S0.P.spkLim_raw_factor;
    end
end

for i = 1:numel(oldPrms)
    propname = oldPrms{i};
    if isfield(hCfg.oldParamSet, propname) && ~isempty(S0.P.(propname))
         % these will be converted automatically
         try
            hCfg.(propname) = S0.P.(propname);
         catch % use default
         end
    end
end

% construct res
res = struct();
res.spikeTimes = S0.viTime_spk;
res.spikeSites = S0.viSite_spk;
res.spikeAmps = S0.vrAmp_spk;
res.spikePositions = S0.mrPos_spk;
nSites = numel(S0.P.viSite2Chan);
res.spikesBySite = vi2cell_(S0.viSite_spk, nSites, 0);
if ~isempty(S0.viSite2_spk)
    res.spikesBySite2 = vi2cell_(S0.viSite2_spk, nSites, 0);
end
res.spikeSites2 = S0.viSite2_spk;
res.siteThresh = S0.vrThresh_site;

% Save local waveforms
trWav_spk = pc2wav_(trPc_spk, S0.mrPv_global); % full waveform
single2int16_ = @(x)int16(x/P.uV_per_bit);
dimm_spk = size(trWav_spk);
nPc = min(size(trPc_spk,1), NUM_PC);
res.rawShape = dimm_spk;
res.filtShape = dimm_spk;
res.featuresShape = [size(trPc_spk,2)*nPc, 2, size(trPc_spk,3)];
res.spikesRaw = single2int16_(trWav_spk); % for now use filtered for raw for backward compatibility
write_bin_([filename_, '_raw.jrc'], res.spikesRaw, 1);
res.spikesFilt = single2int16_(trWav_spk);
write_bin_([filename_, '_filt.jrc'], res.spikesFilt, 1);
trPc2mr_ = @(x)reshape(permute(x(1:nPc,:,:),[2,1,3]), size(x,2)*nPc, 1, size(x,3)); % trim to have 
res.spikeFeatures = cat(2, trPc2mr_(trPc_spk), trPc2mr_(trPc2_spk));
write_bin_([filename_, '_features.jrc'], res.spikeFeatures, 1);
trWav_spk = [];

% copy cluster info
sRes.spikeRho = S_clu.rho;
sRes.spikeDelta = S_clu.delta;
sRes.spikeNeigh = S_clu.nneigh;
sRes.ordRho = S_clu.ordrho;

sRes.spikeClusters = S_auto.viClu;
cviSpk_clu = vi2cell_(S_auto.viClu, S_auto.nClu, 1);
sRes.spikesByCluster = cviSpk_clu;
if isfield(S_clu, 'csNote_clu')
    sRes.clusterNotes = S_clu.csNote_clu;
end

% keep most dense cluster centers from icl post merge
clusterCenters = zeros(max(sRes.spikeClusters), 1);
for iCluster = 1:numel(clusterCenters)
    iSpikes = sRes.spikesByCluster{iCluster};
    iCenters = intersect(S_auto.icl, iSpikes);
    if isempty(iCenters)
        [~, densest] = max(sRes.spikeRho(iSpikes));
        clusterCenters(iCluster) = iSpikes(densest);
    else
        [~, densest] = max(sRes.spikeRho(iCenters));
        clusterCenters(iCluster) = iCenters(densest);
    end
end
sRes.clusterCenters = clusterCenters;
    
% save cluster mean waveforms
trPc_full_spk = expand_pc_(trPc_spk, trPc2_spk, S0, P);
trPc_full_clu = cellfun(@(x)mean(trPc_full_spk(:,:,x),3), cviSpk_clu, 'UniformOutput', 0);
trPc_full_spk = [];
trPc_full_clu = cat(3, trPc_full_clu{:});
trWav_full_spk = pc2wav_(trPc_full_clu, S0.mrPv_global); % full waveform
sRes.meanWfGlobal = trWav_full_spk;
sRes.meanWfGlobalRaw = trWav_full_spk;
trWav_full_spk = [];

sRes.waveformSim = mrWavCor;
sRes.meanWfLocal = trWav_clu;
sRes.meanWfLocalRaw = trWav_clu;
sRes.meanWfRawHigh = trWav_clu;
sRes.meanWfRawLow = trWav_clu;
sRes.clusterSites = viSite_clu;
sRes.initialClustering = S_auto.viClu;

if isfield(S_clu, 'viSite_min_clu')
    sRes.unitPeakSites = S_clu.viSite_min_clu;
end
if isfield(S_clu, 'vnSite_clu')
    sRes.nSitesOverThresh = S_clu.vnSite_clu;
end
sRes.unitCount = S_auto.vnSpk_clu;

if isfield(S_clu, 'vrDc2_site')
    sRes.rhoCutSite = S_clu.vrDc2_site;
end
if isfield(S_clu, 'vrIsiRatio_clu')
    sRes.unitISIRatio = S_clu.vrIsiRatio_clu;
end
if isfield(S_clu, 'vrIsoDist_clu')
    sRes.unitIsoDist = S_clu.vrIsoDist_clu;
end
if isfield(S_clu, 'vrLRatio_clu')
    sRes.unitLRatio = S_clu.vrLRatio_clu;
end
if isfield(S_clu, 'vrPosX_clu') && isfield(S_clu, 'vrPosY_clu')
    sRes.clusterCentroids = [S_clu.vrPosX_clu(:), S_clu.vrPosY_clu(:)];
end
if isfield(S_clu, 'vrSnr_clu')
    sRes.unitSNR = S_clu.vrSnr_clu;
end
if isfield(S_clu, 'vrVmin_clu')
    sRes.unitPeaks = S_clu.vrVmin_clu;
end
if isfield(S_clu, 'vrVmin_clu')
    sRes.unitPeaksRaw = S_clu.vrVmin_clu;
end
if isfield(S_clu, 'vrVpp_clu')
    sRes.unitVpp = S_clu.vrVpp_clu;
end
if isfield(S_clu, 'vrVpp_uv_clu')
    sRes.unitVppRaw = S_clu.vrVpp_uv_clu;
end
if isfield(S_clu, 'vrVrms_site')
    sRes.siteRMS = S_clu.vrVrms_site;
end

hClust = jrclust.sort.DensityPeakClustering(hCfg, sRes, res);
msgs = hClust.inconsistentFields();
assert(isempty(msgs), strjoin(msgs, ', '));

% remove quality scores/initial clustering from sRes
sRes = rmfield_(sRes, {'clusterCentroids', 'clusterNotes', 'initialClustering', ...
                      'meanWfGlobal', 'meanWfGlobalRaw', 'meanWfLocal', 'meanWfLocalRaw', ...
                      'meanWfRawHigh', 'meanWfRawLow', 'nSitesOverThresh', 'waveformSim', ...
                      'siteRMS', 'unitISIRatio', 'unitIsoDist', 'unitLRatio', 'unitPeakSites', ...
                      'unitPeaks', 'unitPeaksRaw', 'unitSNR', 'unitVpp', 'unitVppRaw'});
res = jrclust.utils.mergeStructs(res, sRes);
res.hClust = hClust;

delete_(strrep(prmFile, '.prm', '_res.mat'));

end %func


%--------------------------------------------------------------------------
function S = rmfield_(S, csName)
if ischar(csName), csName={csName}; end
for iField=1:numel(csName)
    if isfield(S, csName{iField})
        S = rmfield(S, csName{iField});
    end
end
end %func


%--------------------------------------------------------------------------
function renameFile(oldfile, newfile)
    try
        movefile(oldfile, newfile);
    catch ME
        warning('failed to rename %s to %s: %s (try to move it manually?)', oldfile, newfile, ME.message);
    end
end


%--------------------------------------------------------------------------
function trPc_spk1 = expand_pc_(trPc1_spk, trPc2_spk, S0, P)

miSites = P.miSites;
[nPc, nSites_spk, nSpk] = size(trPc1_spk);
nSites = numel(P.viSite2Chan);
trPc_spk1 = zeros([nPc, nSites, nSpk], 'single');
if isempty(trPc2_spk)
    cviSite_spk = {S0.viSite_spk};
    ctrPc_spk = {trPc1_spk};
else
    cviSite_spk = {S0.viSite_spk, S0.viSite2_spk};
    ctrPc_spk = {trPc1_spk, trPc2_spk};
end
for iFet = numel(cviSite_spk):-1:1
    cviSpk_site = vi2cell_(cviSite_spk{iFet}, nSites);
    trPc_ = ctrPc_spk{iFet};
    for iSite = 1:nSites
       viSpk1 = cviSpk_site{iSite};
       if isempty(viSpk1), continue; end
       viSite1 = miSites(:,iSite);       
       trPc_spk1(:,viSite1,viSpk1) = trPc_(:,:,viSpk1);
    end
end
end %func


%--------------------------------------------------------------------------
function trWav_spk = pc2wav_(trPc_spk, mrPv)
dimm_pc = size(trPc_spk);
dimm_wav = [size(mrPv,1), dimm_pc(2), dimm_pc(3)];
trWav_spk = reshape(mrPv * reshape(trPc_spk, dimm_pc(1), []), dimm_wav);
end %func


%--------------------------------------------------------------------------
% 4/23/2019 JJJ: delete either cell of files or multiple arguments
function delete_(varargin)
for iArg = 1:nargin
    csFiles = varargin{iArg};
    if ischar(csFiles), csFiles = {csFiles}; end
    for i=1:numel(csFiles)
        try
            if iscell(csFiles)
                delete(csFiles{i});
            else
                delete(csFiles(i));
            end
        catch
    %         disperr_();
        end
    end
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

function varargout = vi2cell_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = file2struct_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = exist_dir_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = write_bin_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
