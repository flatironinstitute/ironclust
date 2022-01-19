function irc2klusters(vcFile_prm, vcDir_out)
% J. James Jun 2019 Aug 16
% Code cited from
%   https://github.com/brendonw1/KilosortWrapper/blob/master/Kilosort2Neurosuite.m
% klusters file format:
%   http://klusters.sourceforge.net/UserManual/data-files.html#cluster-file

% Input
% ----
% EXPORT_MODE: 1 for loading all recordings to memory, 2 for paged read


[fOverwrite, fZeroBase, EXPORT_MODE] = deal(1, 1, 2);

t1 = tic;
if nargin<2, vcDir_out=''; end
if isempty(vcDir_out)
    vcDir_out = fullfile(fileparts(vcFile_prm), 'klusters'); 
end
mkdir_(vcDir_out);

[S0, P] = load_cached_(vcFile_prm);
assert(isfield(S0, 'viTime_spk') && isfield(S0, 'S_clu'), 'incorrect format');
S_clu = S0.S_clu;
nClu = S0.S_clu.nClu;

viTime_spk = uint64(S0.viTime_spk(:));
viClu_spk = uint32(S_clu.viClu(:));
viShank_site = P.viShank_site(:);
[~, vcFile_base] = fileparts(vcFile_prm);
% [nChans, samples] = deal(P.nChans, sum(S0.vnSamples_file));


% compute template
%fprintf('Computing templates per unit\n\t'); t_template = tic;
[mnWav_T, trWav_clu] = deal([]);
P1 = setfield(P, 'fWav_raw_show', 0); % return filered
P1.spkLim = max(abs(P.spkLim)) * [-1,1]; % must be symmetric
nSamples_spk = diff(P1.spkLim) + 1;

switch EXPORT_MODE
    case 2
        viSite_clu = S_clu.viSite_clu;
    case 1
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
end


% Shank loop
fprintf('Exporting shanks\n'); t_shank = tic;
viShank_clu = viShank_site(viSite_clu);
viShank_unique = unique(viShank_clu);
[cviSpk_shank, cviSite_shank, csFile_spk_shank, csFile_fet_shank, cviTime_spk_shank] = deal(cell(size(viShank_unique)));
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
    vcFile_fet1 = fullfile(vcDir_out, sprintf('%s.fet.%d', vcFile_base, iShank1));
    if EXPORT_MODE == 1
        if exist(vcFile_spk1, 'file') ~= 2 || fOverwrite
            mnWav1 = mnWav_T(viShank_site == iShank1, :)';
            mnWav1 = int16(irc('call', 'filt_car', {single(mnWav1), P1}));
            tnWav_spk1 = irc('call', 'mr2tr', {mnWav1, P1.spkLim, viTime_spk1});
            mnWav1 = []; % clear memory
            tnWav_spk1 = permute(tnWav_spk1, [3,1,2]); % dimm: [nSites1, nSamples_spk, nSpk1]    
            write_file_(vcFile_spk1, tnWav_spk1);
            
            if exist(vcFile_fet1, 'file') ~= 2 || fOverwrite
                write_fet_(vcFile_fet1, tnWav_spk1, viTime_spk1);        
            end
            tnWav_spk1 = []; % clear memory
            fprintf('.');
        end
    else
        csFile_spk_shank{iiShank} = vcFile_spk1;
        csFile_fet_shank{iiShank} = vcFile_fet1;
        cviTime_spk_shank{iiShank} = viTime_spk1;        
    end    
    vcFile_par1 = fullfile(vcDir_out, sprintf('%s.par.%d', vcFile_base, iShank1));
    viSite_shank1 = find(viShank_site == iShank1);
    cviSite_shank{iiShank} = viSite_shank1;
end %for

% save spk waveforms by shank
if EXPORT_MODE == 2
    S_shank = makeStruct_(cviSpk_shank, cviSite_shank, csFile_spk_shank, csFile_fet_shank, cviTime_spk_shank);
    save_spk_klusters_(vcFile_prm, P1, S_shank);
end %switch

% write parameter files
vcFile_par = fullfile(vcDir_out, sprintf('%s.par', vcFile_base));
write_par_(vcFile_par, cviSite_shank, P1);
end %func


%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
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
fprintf('\tWrote to %s (took %0.1fs)\n', vcFile, toc(t_write));
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
% function mnWav1 = read_chan_(P, viChan)
% % Get around fread bug (matlab) where built-in fread resize doesn't work
% 
% MAX_MEMORY = 1e9;   % 1GB buffer memory
% 
% if nargin<2, viChan = []; end
% [vcFile, nChans, vcDataType] = deal(P.vcFile, P.nChans, P.vcDataType);
% 
% nLoad = P.nSamples
% 
% if isempty(dimm_wav)
%     mnWav1 = fread(fid_bin, inf, ['*', vcDataType]);
% else
%     if numel(dimm_wav)==1, dimm_wav = [dimm_wav, 1]; end
%     mnWav1 = fread(fid_bin, prod(dimm_wav), ['*', vcDataType]);
%     if numel(mnWav1) == prod(dimm_wav)
%         mnWav1 = reshape(mnWav1, dimm_wav);
%     else
%         dimm2 = floor(numel(mnWav1) / dimm_wav(1));
%         if dimm2 >= 1
%             nSamples1 = dimm_wav(1) * dimm2;
%             mnWav1 = reshape(mnWav1(1:nSamples1), dimm_wav(1), dimm2);
%         else
%             mnWav1 = [];
%         end
%     end
% end
% end %func


%--------------------------------------------------------------------------
% save filtered waveform in a transposed format (_wav.jrc)
% supports single file only and non-transposed
function save_spk_klusters_(vcFile_prm, P1, S_shank)
% ouput
% -----
% S_shank: {cviSpk_shank, cviSite_shank, csFile_spk_shank, csFile_fet_shank}
%   csFile_spk_shank: output file names to store spike waveforms by shank
%   csFile_fet_shank: output file names to store spike features by shank
%   cviSite_shank: site numbers per shank


if nargin<2, P1 = []; end
if nargin<3, S_shank = []; end

NUM_PC = 3;

fprintf('Exporting waveforms and features for all shanks\n\t'); t_spk = tic;
    
[S0, P] = load_cached_(vcFile_prm);
if isempty(P1), P1 = P; end

[nPad_filt, nChans, vcDataType] = deal(P.nPad_filt, P.nChans, P.vcDataType);
vcFile_wav = strrep(P.vcFile_prm, '.prm', '_wav.jrc');
fid_w = [];

% plan loading
[fid1, nBytes_file, header_offset] = fopen_(P.vcFile, 'r');
[nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file, P);
nSamples_file = floor(nBytes_file / bytesPerSample_(vcDataType) / nChans);

% load files in steps with overlap
mnWav_pre1 = [];
iSample_file1 = 0; % zero base
v_fid_spk_shank = [];
cell_pv_shank = cell(size(S_shank.cviSite_shank));

for iLoad1 = 1:nLoad1
    nSamples1 = ifeq_(iLoad1 == nLoad1, nSamples_last1, nSamples_load1);
    mnWav1 = load_file_(fid1, nSamples1, P);
    if iLoad1 < nLoad1
        mnWav_post1 = load_file_preview_(fid1, P);
    else
        mnWav_post1 = [];
    end
    ilim_load1 = iSample_file1 + [1, nSamples1];
    [mnWav1, ilim_wav1] = filter_(mnWav1, mnWav_pre1, mnWav_post1, P);
    switch 2
        case 2 % save waveforms and fet, don't write the waveforms            
            if isempty(v_fid_spk_shank) % open files
                v_fid_spk_shank = cellfun(@(x)fopen(x, 'w'), S_shank.csFile_spk_shank, 'UniformOutput', 1);
                v_fid_fet_shank = cellfun(@(x)fopen(x, 'w'), S_shank.csFile_fet_shank, 'UniformOutput', 1);
            end
            % write to files
            for iShank = 1:numel(S_shank.cviSite_shank)       
                viTime_shank1 = S_shank.cviTime_spk_shank{iShank};
                viSite_shank1 = S_shank.cviSite_shank{iShank};
                [trWav11, viTime11] = mr2tr_sub_(mnWav1(:,viSite_shank1), ilim_wav1, ilim_load1, P1.spkLim, viTime_shank1);
                trWav11 = meanSubt_(trWav11);
                if isempty(cell_pv_shank{iShank})
                    cell_pv_shank{iShank} = calc_pv_(trWav11, NUM_PC);
                end
                mnPc11 = calc_pc_(trWav11, cell_pv_shank{iShank}); % project feature
                
                % write waveform and compute waveform ranges and powers
                mrWav12 = reshape(permute(trWav11, [3,1,2]), [], size(trWav11,2)); % [nSites11, nSamples_spk, nSpk11]
                if ~strcmpi(class_(mrWav12), 'int16')
                    mrWav12 = int16(gather_(mrWav12) * get_set_(P, 'scale_filter', 200)); % make sure dynamic range is properly represented
                end
                fwrite(v_fid_spk_shank(iShank), mrWav12, 'int16'); % [nSites11, nSamples_spk, nSpk11]                
                vn_range11 = int64(range(mrWav12,1));
                vn_power11 = int64(sum(mrWav12.^2,1)/size(mrWav12,1)/100);                                
                clear mrWav12
                
                % save feature vector                                              
                mrFet11 = [mnPc11; vn_range11; vn_power11; int64(viTime11(:)'-1)];
                vcFormat = ['%d', repmat('\t%d', 1, size(mrFet11,1)-1), '\n'];
                if iLoad1 == 1
                    fprintf(v_fid_fet_shank(iShank), '%d\n', size(mrFet11,1));
                end
                fprintf(v_fid_fet_shank(iShank), vcFormat, mrFet11);
            end
            
        case 1 % save filtered waveforms to a file
            mnWav1 = mnWav1(ilim_wav1(1):ilim_wav1(end), :); % trim
            if isempty(fid_w), fid_w = fopen(vcFile_wav, 'w'); end
            fwrite_sub_(fid_w, gather_(mnWav1), iSample_file1, nSamples_file); % provide offset
    end
    % end
    if iLoad1 < nLoad1, mnWav_pre1 = mnWav1(end-P.nPad_filt+1:end, :); end
    clear mnWav1
    iSample_file1 = iSample_file1 + nSamples1;
    fprintf('.');
end %for
fprintf('\n');

% close files
fclose_(v_fid_spk_shank);
fclose_(v_fid_fet_shank);
fclose_(fid_w);
cellfun(@(x)fprintf('\tWrote to %s\n', x), S_shank.csFile_spk_shank);
cellfun(@(x)fprintf('\tWrote to %s\n', x), S_shank.csFile_fet_shank);
fprintf('\ttook %0.1fs\n', toc(t_spk));
end %func


%--------------------------------------------------------------------------
function trPv1 = calc_pv_(trWav1, nPc)
% trWav1: [nSamples_spk, nSpk1, nSites1]
% trPv1: [nSamples_spk, nPc, nSites1]

MAX_SAMPLE = 10000;

[nSamples_spk, nSpk1, nSites1] = size(trWav1);
trPv1 = zeros(nSamples_spk, nPc, nSites1, 'like', trWav1);
iMid = round(nSamples_spk/2);

viSpk_sub = subsample_vr_(1:size(trWav1,2), MAX_SAMPLE);
trWav2 = trWav1(:,viSpk_sub,:);

for iSite=1:nSites1
    mrWav11 = trWav2(:,:,iSite);
    [mrPv11, vrD11] = eig(mrWav11 * mrWav11');
    mrPv11 = zscore_(mrPv11(:,end:-1:end-nPc+1)); % sort largest first
    
    % spike center should be negative
%     vrSign11 = (median((mrPv11' * mrWav11)') > 0)*2-1;
    vrSign11 = (mrPv11(iMid,:) < 0) * 2 - 1; %1 or -1 depending on the sign
    trPv1(:,:,iSite) = bsxfun(@times, mrPv11, vrSign11);    
end %for
end %func


%--------------------------------------------------------------------------
function mnFet_spk1 = calc_pc_(trWav1, trPv)
% waveforms: int16
% capture three PCA
% trWav11: [nSamples_spk, nSpk11, nSites11]   
% trPv: [nSamples_spk, nPc, nSites11]

% constants
[~, nSpk1, nSites1] = size(trWav1);
[~, nPc, ~] = size(trPv);

trFet1 = zeros(nPc, nSpk1, nSites1, 'like', trWav1);
for iSite1 = 1:size(trPv,3)
    trFet1(:,:,iSite1) = trPv(:,:,iSite1)' * trWav1(:,:,iSite1); % [nPc, nSpk11]    
end
trFet1 = gather_(trFet1);
mnFet_spk1 = reshape(permute(trFet1, [1,3,2]), [], nSpk1); % [nPc, nSites, nSpk]
factor = (2^15)./max(abs(mnFet_spk1'));
mnFet_spk1 = int64(mnFet_spk1 .* factor');
end %func


%--------------------------------------------------------------------------
function [trWav_spk1, viTime_spk1] = mr2tr_sub_(mnWav1, ilim1, ilim_load1, spkLim, viTime_spk)
ilim_spk = [find(viTime_spk >= ilim_load1(1), 1, 'first'), find(viTime_spk <= ilim_load1(end), 1, 'last')];
viSpk1 = ilim_spk(1):ilim_spk(end);
viTime1 = viTime_spk(viSpk1) - ilim_load1(1) + ilim1(1);
trWav_spk1 = mr2tr_(mnWav1, spkLim, viTime1);
viTime_spk1 = viTime_spk(viSpk1);
end %func


%--------------------------------------------------------------------------
function [mnWav1, ilim_trim1] = filter_(mnWav_center1, mnWav_pre1, mnWav_post1, P)
[nPre1, nPost1] = deal(size(mnWav_pre1, 1), size(mnWav_post1, 1));
% mnWav1 = cat(1, mnWav_pre1, mnWav_center1, mnWav_post1);
mnWav1 = [mnWav_pre1; mnWav_center1; mnWav_post1];

% try filtering in GPU
fGpu = 0;
if fGpu
    try
        mnWav1 = filt_car_(single(gpuArray_(mnWav1)), P);
    catch
        fGpu = 0;
    end
end
if ~fGpu
    mnWav1 = filt_car_(single(mnWav1), P);
end
ilim_trim1 = [nPre1+1, size(mnWav1,1)-nPost1];
end %func


%--------------------------------------------------------------------------
function fwrite_sub_(fid_w, mnWav1, iSample_file1, nSamples_file)
% input
% -----
% mnWav1: nSamples1 x nChans
[nSamples1, nChans] = size(mnWav1);
nBytesPerSample1 = bytesPerSample_(mnWav1);
offset_pre = 0;
for iChan = 1:nChans
    offset1 = (nSamples_file * (iChan-1) + iSample_file1) * nBytesPerSample1;
    switch 2
        case 2, fseek(fid_w, offset1 - offset_pre, 'cof');
        case 1, fseek(fid_w, offset1, 'bof');
    end
    fwrite_(fid_w, mnWav1(:,iChan));
    offset_pre = offset1;
end
end %func


%--------------------------------------------------------------------------
% call irc
function out1 = mkdir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = fwrite_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = read_chan_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = bytesPerSample_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_file_preview_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = filt_car_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ifeq_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2tr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = fclose_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = zscore_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gpuArray_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gather_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsample_vr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = class_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

function [out1, out2] = load_cached_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end

function [out1, out2, out3] = fopen_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = plan_load_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function varargout = meanSubt_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
