%--------------------------------------------------------------------------
% 4/29/2019 JJJ
% called from `irc convert-mda-xxx`
function convert_mda(vcMode, varargin)

switch vcMode
    case 'reliability', convert_mda_reliability_(varargin{:});
    case 'samplerate', convert_mda_samplerate_(varargin{:});
        
    case 'monotrode', convert_mda_monotrode_(varargin{:});
    case 'gnode', convert_mda_gnode_(varargin{:});
    case 'waveclus', convert_mda_waveclus_(varargin{:});
    case 'buzsaki', convert_mda_buzsaki_(varargin{:}); 
    case 'yass', convert_mda_yass_(varargin{:});
    case 'kampff1', convert_mda_kampff1_(varargin{:});
    case 'kampff2', convert_mda_kampff2_(varargin{:});
    case 'mea', convert_mda_mea_(varargin{:});
    case 'manual', convert_mda_manual_(varargin{:});
    case 'boyden', convert_mda_boyden_(varargin{:});
    case 'mearec100', convert_mda_mearec100_(varargin{:});
        
    case 'emouse', convert_mda_emouse_(varargin{:});
    case 'extract-mda', extract_mda_(varargin{:});
end %switch
end %func


%--------------------------------------------------------------------------
% interpolate mda files
function convert_mda_reliability_(vcDir_in, vcDir_out)
% extract top channels only

if isempty(vcDir_in)
    vcDir_in = '/mnt/ceph/users/jjun/groundtruth/magland_synth/datasets_noise20_K20_C8/';
end
if isempty(vcDir_out)
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/test_reliability';
end

vr_randseed = 1:30; % repeat 10x by adding noise
rand_scale1 = 20/4;

S_json = struct('scale_factor', .01, 'spike_sign', -1);
vS_Files_raw = dir(fullfile(vcDir_in, '**/raw.mda'));
csFiles_raw = cellfun(@(x,y)fullfile(x,y), {vS_Files_raw.folder}, {vS_Files_raw.name}, 'UniformOutput', 0);

% look for raw files in it
for iFile = 1:numel(csFiles_raw)
    vcFile_raw1 = csFiles_raw{iFile};
    vcDir_in1 = fileparts(vcFile_raw1);    
    try        
        % determine the peak channel
        mrGt1 = readmda_(fullfile(vcDir_in1, 'firings_true.mda'));
        mrWav1 = readmda_(fullfile(vcDir_in1, 'raw.mda'))';        
        mrSiteXY1 = csvread(fullfile(vcDir_in1, 'geom.csv'));
        S_json1 = loadjson_(fullfile(vcDir_in1, 'params.json'));
        [nSamples1, nChans1] = size(mrWav1);
        
        for iRand = 1:numel(vr_randseed)
            rng(vr_randseed(iRand));
            % add random number to the array  
            mrWav11 = mrWav1;
            for iChan = 1:nChans1
                mrWav11(:,iChan) = mrWav1(:,iChan) + randn(nSamples1,1,'single') * rand_scale1;
            end
            vcDir_out11 = fullfile(vcDir_out, strrep(vcDir_in1, vcDir_in, ''), sprintf('rand%03d', iRand));
            export_spikeforest_(vcDir_out11, mrWav11', mrSiteXY1, S_json1, mrGt1);        
        end
        fprintf('%d: %s converted\n', iFile, vcDir_in1);
    catch     
        fprintf(2, '%d: error processing %s\n', iFile, vcDir_in1);
    end
end %for iFile
end %func


%--------------------------------------------------------------------------
% interpolate mda files
function convert_mda_samplerate_(vcDir_in, vcDir_out)
% extract top channels only

if isempty(vcDir_in)
    vcDir_in = '/mnt/ceph/users/jjun/groundtruth/magland_synth/datasets_noise20_K20_C8';
end
if isempty(vcDir_out)
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/test_samplerate';
end

vr_samplerate = [15000, 20000, 30000, 40000, 60000];
S_json = struct('scale_factor', .01, 'spike_sign', -1);
vS_Files_raw = dir(fullfile(vcDir_in, '**/raw.mda'));
csFiles_raw = cellfun(@(x,y)fullfile(x,y), {vS_Files_raw.folder}, {vS_Files_raw.name}, 'UniformOutput', 0);

% look for raw files in it
for iFile = 1:numel(csFiles_raw)
    vcFile_raw1 = csFiles_raw{iFile};
    vcDir_in1 = fileparts(vcFile_raw1);    
    try        
        % determine the peak channel
        mrGt1 = readmda_(fullfile(vcDir_in1, 'firings_true.mda'));
        viTime_gt1 = mrGt1(2,:); % one base
        mrWav1 = readmda_(fullfile(vcDir_in1, 'raw.mda'));        
        mrSiteXY1 = csvread(fullfile(vcDir_in1, 'geom.csv'));
        S_json1 = loadjson_(fullfile(vcDir_in1, 'params.json'));
        samplerate1 = S_json1.samplerate;        
        mrFft1 = fft(single(mrWav1)');
        
        for iSampleRate = 1:numel(vr_samplerate)
            samplerate11 = vr_samplerate(iSampleRate);
            S_json11 = setfield(S_json1, 'samplerate', samplerate11);                       
            % scale waveforms using inverse fft
            mrWav11 = resample_ifft_(mrFft1, samplerate1, samplerate11);
            mrWav11 = mrWav11';
            % scale mrGt time index
            viTime_gt11 = round((viTime_gt1-1) * samplerate11/samplerate1) + 1; % scale
            mrGt11 = mrGt1;
            mrGt11(2,:) = viTime_gt11;     
            
            vcDir_out11 = fullfile(vcDir_out, sprintf('sr%d', samplerate11), strrep(vcDir_in1, vcDir_in, ''));
            export_spikeforest_(vcDir_out11, mrWav11, mrSiteXY1, S_json11, mrGt11);        
        end
        fprintf('%d: %s converted\n', iFile, vcDir_in1);
    catch     
        fprintf(2, '%d: error processing %s\n', iFile, vcDir_in1);
    end
end %for iFile
end %func


%--------------------------------------------------------------------------
function mrWav2 = resample_ifft_(mrFft1, samplerate1, samplerate2)
[n1, nChans] = size(mrFft1);
n2 = ceil((n1-1) * samplerate2 / samplerate1) + 1;

mrWav2 = zeros(n2, nChans, 'like', mrFft1);
vrFft2 = zeros(n2, 1, 'like', mrFft1);
if samplerate2 >= samplerate1
    % increase sample rate
    vrFilt2 = fft_lowpass_(n2, samplerate1/2*.8, 1000, samplerate2); 
    n1a = ceil(n1/2);
    n1b = n1-n1a;
    vi1 = [1:n1a, n2-n1b+1:n2];
    for iChan = 1:nChans    
        vrFft2(:) = 0;
        vrFft2(vi1) = mrFft1(:,iChan);
        mrWav2(:,iChan) = real(ifft(vrFft2 .* vrFilt2));
    end
else
    % reduce sample rate
    vrFilt2 = fft_lowpass_(n2, samplerate2/2*.8, 1000, samplerate2); 
    n2a = ceil(n2/2);
    n2b = n2-n2a;
    vi2 = [1:n2a, n1-n2b+1:n1];
    for iChan = 1:nChans    
%         vrFft2(:) = 0;
        vrFft2 = mrFft1(vi2,iChan);
        mrWav2(:,iChan) = real(ifft(vrFft2 .* vrFilt2));
    end    
end
end %func


%--------------------------------------------------------------------------
function filt = fft_lowpass_(N, fhi, fwid_hi, sRateHz)
[n1, n2, f] = get_freq_(N, sRateHz);
filt = sqrt((1-erf((abs(f)-fhi)/fwid_hi))/2);
end %func


%--------------------------------------------------------------------------
function [n1, n2, f] = get_freq_(N, sRateHz)
% [n1, n2, f] = get_freq_(N, sRateHz)
if mod(N,2)==0
    if nargin>=2
        df = sRateHz/N;
        f = df * [0:N/2, -N/2+1:-1]';
    end
    n1 = N/2+1;
else
    if nargin>=2
        df = sRateHz/N;
        f = df * [0:(N-1)/2, -(N-1)/2:-1]'; 
    end
    n1 = (N-1)/2+1;
end
n2 = N-n1;
end %func


%--------------------------------------------------------------------------
% load mda single type and save as int16
function convert_mda_monotrode_(vcDir_in, vcDir_out)
% extract top channels only

if isempty(vcDir_in)
    vcDir_in = '/mnt/ceph/users/jjun/groundtruth/paired_recordings';
end
if isempty(vcDir_out)
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/paired_monotrode';
end

csDir = {'mea64c', 'boyden32c', 'crcns', 'kampff'};

% load timestamps and parse out non-noise clusters
S_json = struct('scale_factor', .01, 'spike_sign', -1);
nChans_out = 4; % number of channels to export

for iDir = 1:numel(csDir)
    vcDir_in1 = fullfile(vcDir_in, csDir{iDir});
    vcDir_out1 = fullfile(vcDir_out, csDir{iDir});    
    try
        export_monotrode_(vcDir_in1, vcDir_out1, nChans_out)
    catch
        disp(lasterr());
    end  
end %for

end %func


%--------------------------------------------------------------------------
% Pierre Yger dataset
function export_monotrode_(vcDir_in, vcDir_out, nChans_out)
spkLim = [-30,30];

% look for .raw files
vS_Files_raw = dir(fullfile(vcDir_in, '**/raw.mda'));
csFiles_raw = cellfun(@(x,y)fullfile(x,y), {vS_Files_raw.folder}, {vS_Files_raw.name}, 'UniformOutput', 0);

% look for raw files in it
for iFile = 1:numel(csFiles_raw)
    vcFile_raw1 = csFiles_raw{iFile};
    vcDir_in1 = fileparts(vcFile_raw1);    
    try        
        % determine the peak channel
        mrGt1 = readmda_(fullfile(vcDir_in1, 'firings_true.mda'));
        viTime_gt1 = mrGt1(2,:);
        mrWav1 = readmda_(fullfile(vcDir_in1, 'raw.mda'));        
        mrWav_mean1 = mr2tr3_(mrWav1', spkLim, viTime_gt1);
        mrWav_mean1 = squeeze(mean(single(mrWav_mean1),2));
        [~, iChan0_in1] = max(max(mrWav_mean1) - min(mrWav_mean1));
        mrSiteXY1 = csvread(fullfile(vcDir_in1, 'geom.csv'));
        [~, viSrt_site_] = sort(pdist2(mrSiteXY1, mrSiteXY1(iChan0_in1,:)), 'ascend');
        nChans_out1 = min(nChans_out, size(mrSiteXY1,1));
        viChan_out1 = sort(viSrt_site_(1:nChans_out1));
        
        % get n-nearest sites
        for iChan_out1 = 1:nChans_out1
            vcDir_out11 = strrep(sprintf('%s_ch%d',vcDir_in1,iChan_out1), vcDir_in, vcDir_out);
            mkdir_(vcDir_out11);            
            copyfile(fullfile(vcDir_in1, 'params.json'), fullfile(vcDir_out11, 'params.json'));
            copyfile(fullfile(vcDir_in1, 'firings_true.mda'), fullfile(vcDir_out11, 'firings_true.mda'));
            iChan_out = viChan_out1(iChan_out1);
            csvwrite(fullfile(vcDir_out11, 'geom.csv'), mrSiteXY1(iChan_out,:));
            writemda_(mrWav1(iChan_out,:), fullfile(vcDir_out11, 'raw.mda'));
        end
        fprintf('%d: %s converted\n', iFile, vcDir_in1);
    catch     
        fprintf(2, '%d: error processing %s\n', iFile, vcDir_in1);
    end
end %for iFile
end %func


%--------------------------------------------------------------------------
% load mda single type and save as int16
function convert_mda_gnode_(vcDir_in, vcDir_out)

if isempty(vcDir_in)
    vcDir_in = '/mnt/home/jjun/ceph/recordings/gnode';
end
if isempty(vcDir_out)
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/waveclus_synth';
end

[csFiles_gt, csFiles_h5] = dir_set_(vcDir_in, '_gt.gdf', '_X.h5');

% load timestamps and parse out non-noise clusters
S_json = struct('samplerate', 24000, 'scale_factor', .01, 'spike_sign', -1);
mrSiteXY = [0, 0];
nChans = 1;
fFlipPolarity = 1; % mean polarity is positive on the waveform
fOverwrite = 0;
for iFile = 1:numel(csFiles_h5)
    try
        % determine input and output directories
        vcFile_gt1 = csFiles_gt{iFile};
        [vcDir1, vcDir2] = fileparts(vcFile_gt1);
        [~,vcDir1] = fileparts(vcDir1);
        vcDir2 = strrep(vcDir2, '_gt','');        
        vcDir_out1 = fullfile(vcDir_out, vcDir1, vcDir2);
        
        % firings_true.mda
        vcFile_h5_ = csFiles_h5{iFile};
        sRateHz1 = h5read(vcFile_h5_, '/srate');
        vrWav1 = single(h5read(vcFile_h5_, '/X'));
        
        mr_gt1 = dlmread(vcFile_gt1);
        S_json1 = setfield(S_json, 'samplerate', sRateHz1);
        nSpk1 = size(mr_gt1,1);        
        mrGt1 = zeros(3, nSpk1, 'double');
        mrGt1(2,:) = mr_gt1(:,2); % factor 4 needed for dif
        mrGt1(3,:) = mr_gt1(:,1); % cluster label, starts with 1
          
        fprintf('\n%d/%d: %s, #ch=%d, sRateHz=%0.1f, duration=%0.1fs, nSpikes=%d\n', ...
            iFile, numel(csFiles_h5), vcDir1, nChans, sRateHz1, mr_gt1(end,2)/sRateHz1, nSpk1);        
            
        % raw.mda
        if ~exist_file_(fullfile(vcDir_out1, 'raw.mda')) || fOverwrite
            vrWav1 = single(h5read(vcFile_h5_, '/X'));
            if fFlipPolarity
                vrWav1 = -vrWav1;
            end
        else
            vrWav1 = [];
        end
        export_spikeforest_(vcDir_out1, vrWav1, mrSiteXY, S_json1, mrGt1);        
    catch
        disp(lasterr());
    end  
end %for

end %func


%--------------------------------------------------------------------------
% load mda single type and save as int16
function convert_mda_waveclus_(vcDir_in, vcDir_out)

if isempty(vcDir_in)
    vcDir_in = '/mnt/ceph/users/jjun/recordings/waveclus/sim2';
end
if isempty(vcDir_out)
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/waveclus_synth/sim2_c1';
end
csFiles_mat = dir_(fullfile(vcDir_in, 'simulation_*.mat'));
S_gt = load(fullfile(vcDir_in, 'ground_truth.mat'));

% load timestamps and parse out non-noise clusters
S_json = struct('samplerate', 24000, 'scale_factor', .01, 'spike_sign', -1);
mrSiteXY = [0, 0];
nChans = 1;
fFlipPolarity = 1; % mean polarity is positive on the waveform
for iFile = 1:numel(csFiles_mat)
    try
        % firings_true.mda
        viClu_spk1 = S_gt.spike_classes{iFile};
        viSpk_time1 = S_gt.spike_first_sample{iFile}; %sampled at 24 KHz
        vi_valid1 = find(viClu_spk1>0);
        [viClu_spk1, viSpk_time1] = deal(viClu_spk1(vi_valid1), viSpk_time1(vi_valid1));
        mrWav_clu1 = S_gt.su_waveforms{iFile}'; %sampled at 96 KHz
        [~,iPeak1] = max(mean(abs(mrWav_clu1),2));
        mrGt1 = zeros(3, numel(viSpk_time1), 'double');
        mrGt1(2,:) = round(viSpk_time1 + iPeak1/4); % factor 4 needed for dif
        mrGt1(3,:) = viClu_spk1; % cluster label, starts with 1
        
        % output directory
        [~, vcDir1] = fileparts(csFiles_mat{iFile});
        vcDir_out1 = fullfile(vcDir_out, vcDir1);
          
        fprintf('\n%d/%d: %s, #ch=%d, sRateHz=%0.1f\n', ...
            iFile, numel(csFiles_mat), vcDir1, nChans, S_json.samplerate);        
            
        % raw.mda
%         if ~exist_file_(fullfile(vcDir_out1, 'raw.mda'))
            vrWav1 = load(csFiles_mat{iFile});
            vrWav1 = single(vrWav1.data);
            if fFlipPolarity
                vrWav1 = -vrWav1;
            end
%         else
%             vrWav1 = [];
%         end
        export_spikeforest_(vcDir_out1, vrWav1, mrSiteXY, S_json, mrGt1);        
    catch
        disp(lasterr());
    end  
end %for

end %func


%--------------------------------------------------------------------------
% 4/29/2019 JJJ: Buzsaki format
function convert_mda_buzsaki_(vcDir_in, vcDir_out)
if nargin<1, vcDir_in = []; end
if nargin<2, vcDir_out = []; end
if isempty(vcDir_in)
    vcDir_in = '/mnt/ceph/users/jjun/recordings/Peter_MS21_180719_155941_concat/Kilosort_2018-12-07_151933/Peter_MS21_180719_155941_concat.spikes.cellinfo.mat';
end
if isempty(vcDir_out)
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/manual_sortings/buzsaki_petersen';
end
S_spikes = load(vcDir_in); S_spikes = S_spikes.spikes;
[viUID_clu, cvrTime_clu, viShankID_clu, viSite_clu] = get_(S_spikes, 'UID', 'times', 'shankID', 'maxWaveformCh1');
cell_channel = cellfun(@(x)readNPY(fullfile(fileparts(vcDir_in), x)), ...
        {'channel_map.npy', 'channel_positions.npy', 'channel_shanks.npy'}, 'UniformOutput',0);
[viSite2Chan, mrPosXY_site, viShank] = deal(cell_channel{1}+1, cell_channel{2}, cell_channel{3});
S_json = struct('samplerate', 20000, 'scale_factor', .195, 'spike_sign', -1);
vcFile_dat = fullfile(fileparts(fileparts(vcDir_in)), 'Peter_MS21_180719_155941_concat.dat');
% mnWav = load_bin_(vcFile, vcDataType, dimm, header)
% create `nShanks` number of ground truth
S_recording = struct('vcDataType', 'int16', 'nChans', 128);
viShank_uniq = unique(viShankID_clu);
viShank_uniq = setdiff(viShank_uniq, [7, 14]); % explicitly exclude shanks
nShanks = numel(viShank_uniq);
mnWav_T = [];
for iShank1 = 1:nShanks
    iShank = viShank_uniq(iShank1);
    viSite1 = find(viShank==iShank);
%     if numel(viSite1) < S_recording.nSites_shank, continue; end
    % reverse lookup for channel_map
    viChan1 = irc('call', 'reverse_lookup', {viSite1, viSite2Chan});
    mrSiteXY1 = mrPosXY_site(viSite1,:);
    vcDir_out1 = fullfile(vcDir_out, sprintf('shank%d', iShank1));
    cvrTime_clu1 = cvrTime_clu(viShankID_clu==iShank);
    viClu1 = cell2mat(arrayfun(@(x)repmat(x, 1, numel(cvrTime_clu1{x})), 1:numel(cvrTime_clu1), 'UniformOutput', 0));
    viTime1 = cell2mat(cvrTime_clu1') * S_json.samplerate;
    [viTime1, viSrt] = sort(viTime1); viClu1 = viClu1(viSrt);    
    mrGt1 = [ones(1, numel(viClu1)); viTime1(:)'; viClu1(:)'];    
    if isempty(mnWav_T)
        if ~exist_file_(fullfile(vcDir_out1, 'raw.mda'))
            mnWav_T = load_bin_(vcFile_dat, S_recording.vcDataType, S_recording.nChans);
        end
    end
    if isempty(mnWav_T)
        mnWav_T1 = [];
    else
        mnWav_T1 = mnWav_T(viChan1,:);
    end
    export_spikeforest_(vcDir_out1, mnWav_T1, mrSiteXY1, S_json, mrGt1);
end
end %func


%--------------------------------------------------------------------------
function convert_mda_yass_(vcDir_in, vcDir_out)

if isempty(vcDir_in)
    vcDir_in = 'K:\PeterLee';
end
if isempty(vcDir_out)
    vcDir_out = 'K:\spikeforest\groundtruth\visapy_mea';
end
[csFiles_gt, csFiles_h5] = dir_set_(vcDir_in, 'ViSAPy_ground_truth.gdf', 'ViSAPy_nonfiltered.h5');

% load timestamps and parse out non-noise clusters
S_json = struct('spike_sign', -1);
sample_offset = 16001;

for iFile = 1:numel(csFiles_gt)
    try
        % firings_true.mda
        gt = textread(csFiles_gt{iFile},'%u');
        gt = reshape(gt,2,[])';
        vi_gt1 = find(gt(:,2) > sample_offset);
        mrGt1 = zeros(3, size(gt,1), 'double');
        mrGt1(2,:) = gt(vi_gt1,2) - sample_offset; % time
        mrGt1(3,:) = gt(vi_gt1,1); % cluster label
        
        % params.json
        sRateHz1 = double(str2num_(h5read(csFiles_h5{iFile}, '/srate')));
        vx = h5read(csFiles_h5{iFile}, '/electrode/x');
        vy = h5read(csFiles_h5{iFile}, '/electrode/y');
        vz = h5read(csFiles_h5{iFile}, '/electrode/z');
        mrSiteXY1 = [vx(:), vz(:)];        
        nChans1 = numel(vx);
        scale_factor1 = 1e-4;
        S_json1 = struct_set_(S_json, 'samplerate', sRateHz1, 'scale_factor', scale_factor1);
        
        % output directory
        [~, vcDir12] = fileparts(fileparts(csFiles_gt{iFile}));
        vcDir_out1 = fullfile(vcDir_out, vcDir12);
          
        fprintf('\n%d/%d: %s, #ch=%d, sRateHz=%0.1f\n', ...
            iFile, numel(csFiles_gt), vcDir12, nChans1, sRateHz1);        
            
        % raw.mda
        if ~exist_file_(fullfile(vcDir_out1, 'raw.mda'))
            mnWav1 = h5read(csFiles_h5{iFile},'/data');
            mnWav1 = mnWav1(:,sample_offset+1:end);        
        else
            mnWav1 = [];
        end
        export_spikeforest_(vcDir_out1, mnWav1, mrSiteXY1, S_json1, mrGt1);        
    catch
        disp(lasterr());
    end  
end %for

end %func


%--------------------------------------------------------------------------
function convert_mda_kampff1_(vcDir_in, vcDir_out)
if isempty(vcDir_in)
    vcDir_in = 'K:\AdamKampff\kampff2\';
end
if isempty(vcDir_out)
    vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\neuropix32c\';
end
viChan_A = [1 3 6 8]' + [0:8:63]; viChan_A = viChan_A(:);
viChan_B = [2 4 5 7]' + [0:8:63]; viChan_B = viChan_B(:);
cviChan = {viChan_A, viChan_B};
csDir = {'2015_09_03_Pair_9_0A', '2015_09_03_Pair_9_0B'};

vcFile_raw = fullfile(vcDir_in, 'raw.mda');
vcFile_geom = fullfile(vcDir_in, 'geom.csv');
mrWav = readmda_(vcFile_raw);
mrSiteXY = csvread(vcFile_geom);

for iDir = 1:numel(csDir)
    vcDir_out1 = fullfile(vcDir_out, csDir{iDir});
    mkdir_(vcDir_out1);
    
    % write geom.csv
    viChan1 = cviChan{iDir};
    csvwrite(fullfile(vcDir_out1, 'geom.csv'), mrSiteXY(viChan1,:));
    
    % copy params.json and gorund truth
    copyfile(fullfile(vcDir_in, 'params.json'), fullfile(vcDir_out1, 'params.json'));
    copyfile(fullfile(vcDir_in, 'firings_true.mda'), fullfile(vcDir_out1, 'firings_true.mda'));
    
    % write raw.mda
    writemda_(mrWav(viChan1,:), fullfile(vcDir_out1, 'raw.mda'));
end %for

end %func


%--------------------------------------------------------------------------
function convert_mda_kampff2_(vcDir_in, vcDir_out)
if isempty(vcDir_in)
    vcDir_in = 'K:\AdamKampff\Neuropixels';
end
if isempty(vcDir_out)
    vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\neuropix32c';
end

S_npx = struct('sRateHz', 30000, 'nChans', 384, 'nChans_out', 32, 'vcDataType', 'int16', 'spike_sign', -1); % 320 um span
S_npx.viChan_max = [235 190 219 182 179 249 184 230 162 176 110 190];
S_npx.csDir_in = {'c28', 'c21', 'c24', 'c46', 'c45', 'c26', 'c19', 'c27', 'c14', 'c44', 'c42', 'c16'};

convert_mda_kampff_(vcDir_in, vcDir_out, S_npx);
end %func


%--------------------------------------------------------------------------
% Pierre Yger dataset
function convert_mda_mea_(vcDir_in, vcDir_out, vcFile_prm)
if isempty(vcDir_in)
    vcDir_in = 'K:\PierreYger\';
end
if isempty(vcDir_out)
    vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\mea64c';
end

% load the probe layout
S_map = load(fullfile(vcDir_in, 'chanMap.mat'));
mrSiteXY_in = [S_map.xcoords(:), S_map.ycoords(:)];
if ~isempty(get_(S_map, 'exclude_chan'))
    mrSiteXY_in(S_map.exclude_chan,:) = inf; % exclude channel will not get selected
end

S_mea = struct('nChans_out', 64, 'spike_sign', -1, 'nChans', 256, 'vcDataType', 'uint16');
S_params = struct('spike_sign', -1, 'samplerate', 20000); 
bytes_per_sample = bytesPerSample_(S_mea.vcDataType); % @TODO: use S_npx.vcDataType

% look for .raw files
vS_Files_raw = dir(fullfile(vcDir_in, '*/*.raw'));
csFiles_raw = cellfun(@(x,y)fullfile(x,y), {vS_Files_raw.folder}, {vS_Files_raw.name}, 'UniformOutput', 0);
csFiles_raw = csFiles_raw(~contains(csFiles_raw, 'juxta'));

% look for raw files in it
for iFile = 1:numel(csFiles_raw)
    vcFile_raw1 = csFiles_raw{iFile};
    [~,vcDir1] = fileparts(fileparts(vcFile_raw1));
    vcFile_txt1 = strrep(vcFile_raw1, '.raw', '.txt');
    vcFile_npy1 = strrep(vcFile_raw1, '.raw', '.triggers.npy');  
%     if ~exist_file_(vcFile_npy1), vcFile_npy1 = strrep(vcFile_raw1, '.raw', '.snippets.npy'); end
    try
        assert(exist_file_(vcFile_npy1) && exist_file_(vcFile_txt1), sprintf('%d: %s error\n', iFile, vcDir1));
        S_txt1 = meta2struct_(vcFile_txt1);
        
        % output files
        vcDir_out1 = fullfile(vcDir_out, vcDir1);
        mkdir_(vcDir_out1);

        % write geom.csv
        iChan0_in1 = S_txt1.channel+1;
        [~, viSrt_site_] = sort(pdist2(mrSiteXY_in, mrSiteXY_in(iChan0_in1,:)), 'ascend');
        viChan_out1 = sort(viSrt_site_(1:S_mea.nChans_out));
        iSite_gt1 = find(viChan_out1 == iChan0_in1);        
        csvwrite(fullfile(vcDir_out1, 'geom.csv'), mrSiteXY_in(viChan_out1,:));

        % write raw.mda
        S_dir1 = dir(vcFile_raw1);
        nBytes = S_dir1.bytes - S_txt1.padding;
        nSamples1 = floor(nBytes / S_mea.nChans / bytes_per_sample);
        mrWav_ = load_bin_(vcFile_raw1, S_mea.vcDataType, [S_mea.nChans, nSamples1], S_txt1.padding);
        mrWav_ = mrWav_(viChan_out1,:);
        writemda(mrWav_, fullfile(vcDir_out1, 'raw.mda'), S_mea.vcDataType);
        
        % write param.json
        struct2json_(S_params, fullfile(vcDir_out1, 'params.json'));

        % export firings_true.mda
        viTime_gt1 = readNPY(vcFile_npy1);
        nSpikes = numel(viTime_gt1);
        mrFirings_out = ones(3, nSpikes, 'double');
        mrFirings_out(1,:) = iSite_gt1;
        mrFirings_out(2,:) = viTime_gt1;
        writemda(mrFirings_out, fullfile(vcDir_out1, 'firings_true.mda'), 'float64');

        fprintf('%d: %s converted\n', iFile, vcDir1);
    catch        
        disperr_();
    end
end %for iFile
end %func


%--------------------------------------------------------------------------
% extract time and put it in the folders
% divide files into 600 s chunk and evaluate accuracy
function convert_mda_manual_(vcDir_in, vcDir_out)
if isempty(vcDir_in)
    vcDir_in = 'K:\JasonChung\tetrodes_manual\';
end
if isempty(vcDir_out)
    vcDir_out = 'K:\spikeforest\groundtruth\manual_sortings\tetrode_1200s\'; % cratet sorter 1,2,3
end
S_out = struct('t_dur_sec', 1200); % export 10 minute recordings

% look for .raw files
vcFile_raw = '20160426_DL15_02_r1.nt16.mda';
csFiles_gt = {'manclust1_firings.mda', 'manclust2_firings.mda', 'manclust3_firings.mda'};
csDir_out = {'sorter1', 'sorter2', 'sorter3'};

mrWav1 = readmda_(fullfile(vcDir_in, vcFile_raw));
nSamples = size(mrWav1,2);
vcFile_param = fullfile(vcDir_in, 'params.json');
S_param = loadjson_(vcFile_param);
nSamples_per_file = S_out.t_dur_sec * S_param.samplerate;
nTime = floor(nSamples / nSamples_per_file);
        
for iFile = 1:numel(csFiles_gt)   
    mrGt1 = readmda_(fullfile(vcDir_in, csFiles_gt{iFile})); 
    viTime_gt1 = mrGt1(2,:);    
    for iTime = 1:nTime
        vcDir_out11 = fullfile(vcDir_out, sprintf('%s_%d', csDir_out{iFile}, iTime));                    
        try
            mkdir_(vcDir_out11);

            % copy params.json and geom.csv
            copyfile(vcFile_param, fullfile(vcDir_out11, 'params.json'));
            copyfile(fullfile(vcDir_in, 'geom.csv'), fullfile(vcDir_out11, 'geom.csv'));

            % write raw.mda
            viTime_out11 = (1:nSamples_per_file) + (iTime-1)*nSamples_per_file;        
            writemda_(mrWav1(:,viTime_out11), fullfile(vcDir_out11, 'raw.mda'));

            % write firings_true.mda
            mrGt11 = mrGt1(:,viTime_gt1 >= viTime_out11(1) & viTime_gt1 <= viTime_out11(end));
            mrGt11(2,:) = mrGt11(2,:) - (viTime_out11(1) - 1);
            writemda_(mrGt11, fullfile(vcDir_out11, 'firings_true.mda'));

            fprintf('%d_%d: converted to %s\n', iFile, iTime, vcDir_out11);
        catch
            fprintf(2, '%d_%d: error creating %s\n', iFile, iTime, vcDir_out11);
        end
    end %for    
end %for
end %func


%--------------------------------------------------------------------------
function convert_mda_boyden_(vcDir_in, vcDir_out)
if isempty(vcDir_in)
    vcDir_in = 'K:\spikeforest\groundtruth\paired_recordings\boyden\';
end
if isempty(vcDir_out)
    %vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\boyden32c\';
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/paired_recordings\boyden32c';
end
S_mea = struct('nChans_out', 32);

% look for .raw files
vS_Files_raw = dir(fullfile(vcDir_in, '*/raw.mda'));
csFiles_raw = cellfun(@(x,y)fullfile(x,y), {vS_Files_raw.folder}, {vS_Files_raw.name}, 'UniformOutput', 0);

% look for raw files in it
for iFile = 1:numel(csFiles_raw)
    vcFile_raw1 = csFiles_raw{iFile};
    vcDir_in1 = fileparts(vcFile_raw1);
    [~, vcDir1] = fileparts(vcDir_in1);
    vcDir_out1 = fullfile(vcDir_out, vcDir1);
    try
%         if exist_dir_(vcDir_out1), continue; end
        mkdir_(vcDir_out1);
        vcFile_out_raw1 = fullfile(vcDir_out1, 'raw.mda');
        if exist_file_(vcFile_out_raw1), continue; end
        
        copyfile(fullfile(vcDir_in1, 'params.json'), fullfile(vcDir_out1, 'params.json'));
        copyfile(fullfile(vcDir_in1, 'firings_true.mda'), fullfile(vcDir_out1, 'firings_true.mda'));
        
        % determine the peak channel
        mrGt1 = readmda_(fullfile(vcDir_in1, 'firings_true.mda'));
        viTime_gt1 = mrGt1(2,:);
        mrWav1 = readmda_(fullfile(vcDir_in1, 'raw.mda'));
        mrWav_mean1 = mr2tr3_(mrWav1', [-30,30], viTime_gt1);
        mrWav_mean1 = squeeze(mean(mrWav_mean1,2));
        [~, iChan0_in1] = max(max(mrWav_mean1)-min(mrWav_mean1));
        
        % write geom.csv
        mrSiteXY1 = csvread(fullfile(vcDir_in1, 'geom.csv'));
        [~, viSrt_site_] = sort(pdist2(mrSiteXY1, mrSiteXY1(iChan0_in1,:)), 'ascend');
        viChan_out1 = sort(viSrt_site_(1:S_mea.nChans_out));
        iSite_gt1 = find(viChan_out1 == iChan0_in1);        
        csvwrite(fullfile(vcDir_out1, 'geom.csv'), mrSiteXY1(viChan_out1,:));

        % write raw.mda
        writemda_(mrWav1(viChan_out1,:), vcFile_out_raw1);       

        fprintf('%d: %s converted\n', iFile, vcDir1);
    catch     
        fprintf(2, '%d: error processing %s\n', iFile, vcDir_in1);
    end
end %for iFile
end %func


%--------------------------------------------------------------------------
% Apr 28, 2019 JJJ
% convet mearec 100chan format to 64 channel format
function convert_mda_mearec100_(vcDir_in, vcDir_out)
if nargin<1, vcDir_in = []; end
if nargin<2, vcDir_out = []; end

if isempty(vcDir_in)
    vcDir_in = '/mnt/ceph/users/jjun/groundtruth/mearec_synth/sqmea';
end
if isempty(vcDir_out)
    vcDir_out = '/mnt/ceph/users/jjun/groundtruth/mearec_synth/sqmea64c';
end
S_settings = struct('nChans_out', 64, 'spkLim', [-30,30], 'vcDataType', 'single');

% look for .raw files
vS_Files_raw = dir(fullfile(vcDir_in, '**/raw.mda'));
csFiles_raw = cellfun(@(x,y)fullfile(x,y), {vS_Files_raw.folder}, {vS_Files_raw.name}, 'UniformOutput', 0);

% look for raw files in it
for iFile = 1:numel(csFiles_raw)
    vcFile_raw1 = csFiles_raw{iFile};
    vcDir_in1 = fileparts(vcFile_raw1);
    vcDir_out1 = strrep(strrep(vcDir_in1, vcDir_in, vcDir_out), '_C100', '_C64');
    try
%         if exist_dir_(vcDir_out1), continue; end
        mkdir_(vcDir_out1);
        vcFile_out_raw1 = fullfile(vcDir_out1, 'raw.mda');
        if exist_file_(vcFile_out_raw1), continue; end
        
        copyfile(fullfile(vcDir_in1, 'params.json'), fullfile(vcDir_out1, 'params.json'));
        copyfile(fullfile(vcDir_in1, 'firings_true.mda'), fullfile(vcDir_out1, 'firings_true.mda'));
        
        % determine the peak channel
        mrGt1 = readmda_(fullfile(vcDir_in1, 'firings_true.mda'));
        viTime_gt1 = mrGt1(2,:);
        mrWav1 = readmda_(fullfile(vcDir_in1, 'raw.mda'));
        mrWav1 = cast(mrWav1, S_settings.vcDataType);
        mrWav_mean1 = mr2tr3_(mrWav1', S_settings.spkLim, viTime_gt1);
        mrWav_mean1 = squeeze(mean(mrWav_mean1,2));
        [~, iChan0_in1] = max(max(mrWav_mean1)-min(mrWav_mean1));
        mrSiteXY1 = csvread(fullfile(vcDir_in1, 'geom.csv'));
        % write geom.csv
        switch 2
            case 2 % nearest square site arrangement
                xy0 = mrSiteXY1(iChan0_in1,:);                
                vrX0 = unique(mrSiteXY1(:,1));
                vrY0 = unique(mrSiteXY1(:,2));                
                [~,viX_] = sort(abs(vrX0-xy0(1)));
                [~,viY_] = sort(abs(vrY0-xy0(2)));
                vrX_ = vrX0(viX_(1:sqrt(S_settings.nChans_out)));                
                vrY_ = vrY0(viY_(1:sqrt(S_settings.nChans_out)));                                                
                viChan_out1 = find(ismember(mrSiteXY1(:,1), vrX_) & ismember(mrSiteXY1(:,2), vrY_));
                
            case 1 % nearest n sites                
                [~, viSrt_site_] = sort(pdist2(mrSiteXY1, mrSiteXY1(iChan0_in1,:)), 'ascend');
                viChan_out1 = sort(viSrt_site_(1:S_settings.nChans_out));                
        end
        iSite_gt1 = find(viChan_out1 == iChan0_in1);   % not used     
        csvwrite(fullfile(vcDir_out1, 'geom.csv'), mrSiteXY1(viChan_out1,:));

        % write raw.mda
        writemda_(mrWav1(viChan_out1,:), vcFile_out_raw1);       
        fprintf('%d: %s converted\n', iFile, vcDir_out1);
    catch     
        fprintf(2, '%d: error processing %s\n', iFile, vcDir_in1);
    end
end %for iFile
end %func


%--------------------------------------------------------------------------
function convert_mda_emouse_(vcDir_in, vcDir_out)

if nargin<1, vcDir_in=''; end
if nargin<2, vcDir_out=''; end
assert(~isempty(vcDir_in), 'convert_mda: convert_mda_emouse_: `vcDir_in` must be provided');
% if isempty(vcDir_in), vcDir_in='/mnt/ceph/users/jjun/groundtruth/hybrid_synth/linear_drift/'; end
if isempty(vcDir_out), vcDir_out=vcDir_in; end

% make directory
mkdir_(vcDir_out);

% write geom.csv
S_prb = load(fullfile(vcDir_in, 'chanMap_3B_64sites.mat'));
vcFile_geom = fullfile(vcDir_out, 'geom.csv');
csvwrite(vcFile_geom, [S_prb.xcoords(:), S_prb.ycoords(:)]);
fprintf('Wrote to %s\n', vcFile_geom);
nChans = numel(S_prb.xcoords);
sRateHz = S_prb.fs;

% write raw.mda
vcFile_bin = fullfile(vcDir_in, 'sim_binary.imec.ap.bin');
if exist_file_(vcFile_bin)
    vcFile_raw = fullfile(vcDir_out, 'raw.mda');
    append_bin2mda_({vcFile_bin}, vcFile_raw, nChans, 'int16');
    fprintf('Wrote to %s\n', vcFile_raw);
else
    fprintf(2, '%s does not exist\n', vcFile_bin);
end

% Write to params.json
vcFile_json = fullfile(vcDir_out, 'params.json');
struct2json_(struct('spike_sign', -1, 'samplerate', sRateHz), vcFile_json);
% fprintf('Wrote to %s\n', vcFile_json);

% write firings_true.mda
vcFile_true = fullfile(vcDir_out, 'firings_true.mda');
S_gt = load(fullfile(vcDir_in, 'eMouseGroundTruth.mat'));
mrGt = zeros(numel(S_gt.gtRes), 3);
mrGt(:,2) = S_gt.gtRes(:);
mrGt(:,3) = S_gt.gtClu(:);
writemda_(vcFile_true, mrGt');
fprintf('Wrote to %s\n', vcFile_true);
end %func


%--------------------------------------------------------------------------
function extract_mda_(vcDir_in, vcDir_out, viSite, tLim)

fOverwrite = 1;

if nargin<1, vcDir_in=''; end
if nargin<2, vcDir_out=''; end
assert(~isempty(vcDir_in), 'convert_mda: extract_mda_: `vcDir_in` must be provided');
if isempty(vcDir_out), vcDir_out=vcDir_in; end
if ischar(viSite), viSite = str2num(viSite); end
if ischar(tLim), tLim=str2num(tLim); end

% make directory
mkdir_(vcDir_out);

% write geom.csv
mrSiteXY = csvread(fullfile(vcDir_in, 'geom.csv'));
vcFile_geom = fullfile(vcDir_out, 'geom.csv');
csvwrite(vcFile_geom, mrSiteXY(viSite,:));
fprintf('Wrote to %s\n', vcFile_geom);

% write raw.mda
vcFile_raw = fullfile(vcDir_out, 'raw.mda');
if ~exist_file_(vcFile_raw) || fOverwrite
    [S_mda, fid_r] = readmda_header_(fullfile(vcDir_in, 'raw.mda'));
    [vcDataType, dimm] = get_(S_mda, 'vcDataType', 'dimm');
    nChans = dimm(1);
    nBytes_skip = (round(tLim(1)*dimm(2))-1) * bytesPerSample_(vcDataType) * nChans;
    if nBytes_skip>0, fseek(fid_r, nBytes_skip, 'cof'); end
    nSamples_copy = round(diff(tLim)*dimm(2));
    mr_ = fread_(fid_r, [nChans, nSamples_copy], vcDataType);
    fclose(fid_r);
    writemda_(vcFile_raw, mr_(viSite,:));
    fprintf('Wrote to %s\n', vcFile_raw);
else
    S_mda = readmda_header_(fullfile(vcDir_in, 'raw.mda'));
end
nLim = round(S_mda.dimm(2) * tLim);
if nLim(1)<1, nLim(1)=1; end

% Write to params.json
vcFile_json = fullfile(vcDir_out, 'params.json');
copyfile(fullfile(vcDir_in, 'params.json'), vcFile_json, 'f');
fprintf('Wrote to %s\n', vcFile_json);

% write firings_true.mda
vcFile_true = fullfile(vcDir_out, 'firings_true.mda');
mrGt = readmda_(fullfile(vcDir_in, 'firings_true.mda'))';
mrGt = mrGt(mrGt(:,2) >= nLim(1) & mrGt(:,2) <= nLim(2), :);
mrGt(:,2) = mrGt(:,2) - nLim(1) + 1;
writemda_(vcFile_true, mrGt');
fprintf('Wrote to %s\n', vcFile_true);
end %func


%--------------------------------------------------------------------------
% Call from irc.m
function cout = call_irc2_(dbstack1, cell_input, nargout)
vcFunc = dbstack1(1).name;
try
    switch nargout
        case 0, cout{1} = []; irc2('call', vcFunc, cell_input);
        case 1, cout{1} = irc2('call', vcFunc, cell_input);
        case 2, [cout{1}, cout{2}] = irc2('call', vcFunc, cell_input);
        case 3, [cout{1}, cout{2}, cout{3}] = irc2('call', vcFunc, cell_input);
        case 4, [cout{1}, cout{2}, cout{3}, cout{4}] = irc2('call', vcFunc, cell_input);
        otherwise, error('call_irc2_: undefined func: %s', vcFunc);
    end
catch ME
    fprintf(2, 'call_irc_: %s\n', ME.message);
    rethrow ME;
end
end %func


%--------------------------------------------------------------------------
% 1/31/2019 JJJ: get the field(s) of a struct or index of an array or cell
function varargout = get_(varargin)
% same as struct_get_ function
% retrieve a field. if not exist then return empty
% [val1, val2] = get_(S, field1, field2, ...)
% [val] = get_(cell, index)

if nargin==0, varargout{1} = []; return; end
S = varargin{1};
if isempty(S), varargout{1} = []; return; end

if isstruct(S)
    for i=2:nargin
        vcField = varargin{i};
        try
            varargout{i-1} = S.(vcField);
        catch
            varargout{i-1} = [];
        end
    end
elseif iscell(S)
    try    
        varargout{1} = S{varargin{2:end}};
    catch
        varargout{1} = [];
    end
else
    try    
        varargout{1} = S(varargin{2:end});
    catch
        varargout{1} = [];
    end
end
end %func


%--------------------------------------------------------------------------
% irc2.m
function varargout = export_spikeforest_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = write_bin_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = append_bin2mda_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = writemda_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = loadjson_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct2json_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = mkdir_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = exist_file_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = ifeq_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = str2num_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = readmda_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = dir_set_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = dir_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = prb2geom_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = disperr_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = convert_mda_kampff_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end

function varargout = mr2tr3_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = chan2site_prb_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_bin_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = file_dimm_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = gt2mda_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = loadParam_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_batch_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_gt_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = bytesPerSample_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = meta2struct_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_set_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end

function varargout = readmda_header_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = trim_wav_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = fread_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
