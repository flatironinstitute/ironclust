function irc2klusters(vcFile_prm, vcDir_out)
% J. James Jun 2019 Jun 27
% modified from https://github.com/brendonw1/KilosortWrapper/blob/master/Kilosort2Neurosuite.m
% Original author: 
%   By Peter Petersen 2018
%   petersen.peter@gmail.com
% klusters file format:
%   http://klusters.sourceforge.net/UserManual/data-files.html#cluster-file

[fOverwrite, fZeroBase] = deal(0, 1);

t1 = tic;
if nargin<2, vcDir_out=''; end
if isempty(vcDir_out)
    vcDir_out = fullfile(fileparts(vcFile_prm), 'klusters'); 
end
mkdir_(vcDir_out);

[S0, P] = irc('call', 'load_cached_', {vcFile_prm});
assert(isfield(S0, 'viTime_spk') && isfield(S0, 'S_clu'), 'incorrect format');
S_clu = S0.S_clu;
nClu = S0.S_clu.nClu;

viTime_spk = uint64(S0.viTime_spk(:));
viClu_spk = uint32(S_clu.viClu(:));
viShank_site = P.viShank_site(:);
[~, vcFile_base] = fileparts(vcFile_prm);
% [nChans, samples] = deal(P.nChans, sum(S0.vnSamples_file));


% compute template
fprintf('Computing templates per unit\n\t'); t_template = tic;
[mnWav_T, trWav_clu] = deal([]);
P1 = setfield(P, 'fWav_raw_show', 0); % return filered
P1.spkLim = max(abs(P.spkLim)) * [-1,1]; % must be symmetric
nSamples_spk = diff(P1.spkLim) + 1;
viSite_clu = zeros(nClu,1);
for iClu = 1:nClu
    viTime_spk1 = S0.viTime_spk(S_clu.cviSpk_clu{iClu});
    [mrWav_clu1, mnWav_T] = irc('call','load_wav_med_', {P1, viTime_spk1, mnWav_T}); % not filtering   
    if isempty(trWav_clu)
        trWav_clu = zeros(size(mrWav_clu1,1), size(mrWav_clu1,2), nClu, 'single');
    else
        trWav_clu(:,:,iClu) = mrWav_clu1;
    end
    [~, viSite_clu(iClu)] = max(range(mrWav_clu1));
    fprintf('.');
end %for
trWav_clu = permute(trWav_clu, [2,1,3]);
mnWav_T = int16(mnWav_T);
fprintf('\n\ttook %0.1fs\n', toc(t_template));


% Shank loop
fprintf('Exporting shanks\n\t'); t_shank = tic;
viShank_clu = viShank_site(viSite_clu);
viShank_unique = unique(viShank_clu);
[cviSpk_shank, cviSite_shank] = deal(cell(size(viShank_unique)));
for iiShank = 1:numel(viShank_unique)    
    iShank1 = viShank_unique(iiShank);
    fprintf('Loading data for shank %d:\n', iShank1);
    viClu_shank1 = find(viShank_clu == iShank1);
    viSpk_shank1 = find(ismember(viClu_spk,viClu_shank1));
    cviSpk_shank{iiShank} = viSpk_shank1;    
    nSpk1 = numel(viSpk_shank1);
    viClu_spk1 = viClu_spk(viSpk_shank1);
    nClu1 = numel(viClu_shank1);
    
    vcFile_clu1 = fullfile(vcDir_out, sprintf('%s.clu.%d', vcFile_base, iShank1));
    write_file_(vcFile_clu1, [nClu1; viClu_spk(viSpk_shank1)-fZeroBase], '%.0f\n'); % zero-base

    vcFile_res1 = fullfile(vcDir_out, sprintf('%s.res.%d', vcFile_base, iShank1));
    viTime_spk1 = viTime_spk(viSpk_shank1);
    write_file_(vcFile_res1, viTime_spk1-fZeroBase, '%.0f\n'); % zero-base
    
    vcFile_spk1 = fullfile(vcDir_out, sprintf('%s.spk.%d', vcFile_base, iShank1));    
    if exist(vcFile_spk1, 'file') ~= 2 || fOverwrite
        mnWav1 = mnWav_T(viShank_site == iShank1, :)';
        mnWav1 = int16(irc('call', 'filt_car', {single(mnWav1), P1}));
        tnWav_spk1 = irc('call', 'mr2tr', {mnWav1, P1.spkLim, viTime_spk1});
        mnWav1 = []; % clear memory
        tnWav_spk1 = permute(tnWav_spk1, [3,1,2]); % dimm: [nSites1, nSamples_spk, nSpk1]    
        write_file_(vcFile_spk1, tnWav_spk1);
    end
    
    vcFile_fet1 = fullfile(vcDir_out, sprintf('%s.fet.%d', vcFile_base, iShank1));
    if exist(vcFile_fet1, 'file') ~= 2 || fOverwrite
        write_fet_(vcFile_fet1, tnWav_spk1, viTime_spk1);        
    end
    tnWav_spk1 = []; % clear memory
    
    vcFile_par1 = fullfile(vcDir_out, sprintf('%s.par.%d', vcFile_base, iShank1));
    viSite_shank1 = find(viShank_site == iShank1);
    cviSite_shank{iiShank} = viSite_shank1;
    
    fprintf('.');
end %for
fprintf('\n\ttook %0.1fs\n', toc(t_shank));


% write parameter
vcFile_par = fullfile(vcDir_out, sprintf('%s.par', vcFile_base));
write_par_(vcFile_par, cviSite_shank, P1);

end %func


%--------------------------------------------------------------------------
function write_file_(vcFile, vnData, vcFormat)
% mode 1: write formatted numbers
%  write_file_(vcFile, vnData, vcFormat)
% mode 2: write binary
%  write_file_(vcFile, vnData)
% mode 3: write text (separated by new lines)
%  write_file_(vcFile, csData)

if nargin<3, vcFormat = []; end

t_write = tic;
fid=fopen(vcFile, 'w');
if ~isempty(vcFormat)
    fprintf(fid, vcFormat, vnData);
else
    if iscell(vnData)
        csLines = vnData;
        for iLine = 1:numel(csLines)
            fprintf(fid, '%s\n', csLines{iLine});
        end %for
    else
        fwrite(fid, vnData, class(vnData));
    end
end
fclose(fid);
fprintf('Wrote to %s (took %0.1fs)\n', vcFile, toc(t_write));
end %func


%--------------------------------------------------------------------------
function write_fet_(vcFile_fet1, tnWav_spk1, viTime_spk1)
% waveforms: int16
% capture three PCA
% tnWav_spk1: [nSites1 x nSamples_spk x nSpk1]

% constants
[fZeroBase, nPc] = deal(1, 3);

[nSites1, nSamples1, nSpk1] = size(tnWav_spk1);
mnWav_spk1 = reshape(tnWav_spk1,[],nSpk1);
wranges = int64(range(mnWav_spk1,1));
wpowers = int64(sum(mnWav_spk1.^2,1)/size(mnWav_spk1,1)/100);
mnWav_spk1 = [];

fprintf('Computing PCA\n\t');
t_pca = tic;
tnWav_spk2 = permute(tnWav_spk1, [2,3,1]);
trFet_spk1 = zeros([nSpk1, nPc, nSites1], 'single');
for iSite = 1:nSites1
    mrFet1 = zscore(single(tnWav_spk2(:,:,iSite)),[],2);
    trFet_spk1(:,:,iSite) = pca(mrFet1, 'NumComponents',nPc, 'Centered', false);
    fprintf('.');
end
tnWav_spk2 = []; 
mnFet_spk1 = reshape(permute(trFet_spk1, [2,3,1]), [], nSpk1);
factor = (2^15)./max(abs(mnFet_spk1'));
mnFet_spk1 = int64(mnFet_spk1 .* factor');
fprintf('\n\ttook %0.1fs\n', toc(t_pca));

Fet = [mnFet_spk1; wranges; wpowers; int64(viTime_spk1(:)'-fZeroBase)];
nFeatures = size(Fet, 1); 
vcFormat = ['%d', repmat('\t%d', 1, nFeatures-1), '\n'];

% write to file
t_write = tic;
fid = fopen(vcFile_fet1, 'w');
fprintf(fid, '%d\n', nFeatures);
fprintf(fid, vcFormat, Fet);
fclose(fid);
fprintf('Wrote to %s (took %0.1fs)\n', vcFile_fet1, toc(t_write));
end %func


%--------------------------------------------------------------------------
% write par files for each shanks (easier than writing an xml file)
function write_par_(vcFile_par, cviSite_shank, P)

[nBitsPerSample, fZeroBase, nPc] = deal(16, 1, 3);
highPass = round(P.freqLim(1));
sample_interval_us = round(1e6/P.sRateHz);
nShanks = numel(cviSite_shank);
nSites = sum(cellfun(@sum, cviSite_shank));

csLines_ = {};
csLines_{end+1} = sprintf('%d %d', nSites, nBitsPerSample);
csLines_{end+1} = sprintf('%d %d', sample_interval_us, highPass);
csLines_{end+1} = sprintf('%d', nShanks);
for iiShank = 1:nShanks
    viSites1 = cviSite_shank{iiShank} - fZeroBase;
    csLines_{end+1} = sprintf('%d, %s', numel(viSites1), sprintf('%d ', viSites1));
end
write_file_(vcFile_par, csLines_);


for iiShank = 1:nShanks
    viSite_shank1 = cviSite_shank{iiShank} - fZeroBase;
    nSites1 = numel(viSite_shank1);
    nSamples_spk = diff(P.spkLim)+1;
    
    csLines_ = {};
    csLines_{end+1} = sprintf('%d %d %d', nSites1, nSites1, sample_interval_us);
    csLines_{end+1} = sprintf('%d ', viSite_shank1);
    csLines_{end+1} = '10 2'; %  # refractory sample index after detection, RMS integration window length
    csLines_{end+1} = '90'; %  # approximate firing frequency in Hz
    csLines_{end+1} = sprintf('%d %d', nSamples_spk, 1-P.spkLim(1)); % # number of samples in each waveform, sample index of the peak
    csLines_{end+1} = '12 6'; % # window length to realign the spikes, sample index of the peak (detection program)
    csLines_{end+1} = '4 4'; % # number of samples (before and after the peak) to use for reconstruction and features
    csLines_{end+1} = sprintf('%d %d', nPc, nSamples_spk); % # number of principal components (features) per electrode, number of samples used for the PCA
    csLines_{end+1} = sprintf('%d', highPass);
    
    vcFile_par1 = strrep(vcFile_par, '.par', sprintf('.par.%d', iiShank));
    write_file_(vcFile_par1, csLines_);
end % for
end %func


%--------------------------------------------------------------------------
function mkdir_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
%function out1 = title_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end