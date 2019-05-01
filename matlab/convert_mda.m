%--------------------------------------------------------------------------
% 4/29/2019 JJJ
% called from `irc convert-mda-xxx`
function convert_mda(vcMode, varargin)

switch vcMode
    case 'buzsaki', convert_mda_buzsaki_(varargin{:}); 
    case 'yass', convert_mda_yass_(varargin{:});
    case 'kampff1', convert_mda_kampff1_(varargin{:});
    case 'kampff2', convert_mda_kampff2_(varargin{:});
    case 'mea', convert_mda_mea_(varargin{:});
    case 'manual', convert_mda_manual_(varargin{:});
    case 'boyden', convert_mda_boyden_(varargin{:});
    case 'mearec100', convert_mda_mearec100_(varargin{:});
end %switch
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
S_recording = struct('vcDataType', 'int16', 'nChans', 128, 'nSites_shank', 10);
viShank_uniq = unique(viShankID_clu);
nShanks = numel(viShank_uniq);
mnWav_T = [];
for iShank1 = 1:nShanks
    iShank = viShank_uniq(iShank1);
    viSite1 = find(viShank==iShank);
    if numel(viSite1) < S_recording.nSites_shank, continue; end
    % reverse lookup for channel_map
    viChan1 = irc('call', 'reverse_lookup', {viSite1, viSite2Chan});
    mrSiteXY1 = mrPosXY_site(viSite1,:);
    vcDir_out1 = fullfile(vcDir_out, sprintf('shank%d', iShank));
    cvrTime_clu1 = cvrTime_clu(viShankID_clu==iShank);
    viClu1 = cell2mat(arrayfun(@(x)repmat(x, 1, numel(cvrTime_clu1{x})), 1:numel(cvrTime_clu1), 'UniformOutput', 0));
    viTime1 = cell2mat(cvrTime_clu1');
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
    vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\boyden32c\';
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


%==========================================================================
% call irc.m

%--------------------------------------------------------------------------
% 0 output
function writemda_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function prb2geom_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function struct2json_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function export_spikeforest_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function disperr_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function convert_mda_kampff_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function mkdir_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end

%--------------------------------------------------------------------------
% 1 output
function out1 = exist_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = loadjson_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2tr3_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = chan2site_prb_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_bin_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file_dimm_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gt2mda_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ifeq_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = loadParam_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_batch_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_gt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = str2num_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = bytesPerSample_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = readmda_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = meta2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

%--------------------------------------------------------------------------
% 2 outputs
function [out1, out2] = dir_set_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end

%--------------------------------------------------------------------------
% varargout

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