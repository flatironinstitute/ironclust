% works on mda format only for now
% no manual interface now

%--------------------------------------------------------------------------
function irc2(vcDir_in, vcDir_out, vcFile_arg)
if nargin<1, vcDir_in = get_test_data_(); end
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end

if isempty(vcDir_out)
    vcDir_out = strrep(vcDir_in, 'groundtruth', 'irc');
end

P = makeParam_(vcDir_in, vcDir_out, vcFile_arg);
S0 = detect_(P);
S0 = sort_(S0, P);
S0 = auto_(S0, P);
save_(S0);
validate_(S0);
end %func


%--------------------------------------------------------------------------
function S0 = sort_(S0, P)

% drift processing

% cluster, channel loop
% rho
% delta

end %func


%--------------------------------------------------------------------------
function S0 = detect_(P)
% keep all features in the memory, no disk storage

% parfor loop
if get_set_(P, 'fParfor', 1)
    gcp_ = gcp();
else
    gcp_ = [];
end

% P.fGpu = 1;
% else
%     gcp_ = [];
% end

% fParfor=0, fGpu=0: 90s
% fParfor=0, fGpu=1: 26s **GPU is fastest**
% fParfor=1, fGpu=0: 30s (local, 9 workers)
% fParfor=1, fGpu=1: 35s (local, 9 workers)
t1=tic;

% load one
S_paged = readmda_paged_(P); % initialize
[nLoads, viOffset_load] = deal(S_paged.nLoads, S_paged.viOffset_load);
[mrWav1, nlim_wav1, fDone] = readmda_paged_(); % process first part
cS_detect = cell(nLoads, 1);
cS_detect{1} = detect_paged_(mrWav1, P, makeStruct_(nlim_wav1)); % process the first part
[vrThresh_site, mrPv_global] = struct_get_(cS_detect{1}, 'vrThresh_site', 'mrPv_global');
S_cache1 = makeStruct_(vrThresh_site, mrPv_global);

for iLoad = 2:nLoads          
    [mrWav1, nlim_wav1, fDone] = readmda_paged_(); % process first part
    S_cache1.nlim_wav1 = nlim_wav1; % trim waveform
    if isempty(gcp_)            
        cS_detect{iLoad} = detect_paged_(mrWav1, P, S_cache1);            
    else
        vS_out(iLoad-1) = parfeval(gcp_, @(x,y)detect_paged_(x,P,y), 1, mrWav1, S_cache1);
    end    
end
if ~isempty(gcp_)
    for iLoad = 2:nLoads
        [completedIdx, S_] = fetchNext(vS_out);
        cS_detect{completedIdx+1} = S_;
    end
end
t_detect = toc(t1);

% Save output
[viSite_spk, viTime_spk, vrAmp_spk, trPc_spk] = detect_merge_(cS_detect, viOffset_load);
S0 = makeStruct_(viSite_spk, viTime_spk, vrAmp_spk, trPc_spk, vrThresh_site, mrPv_global, t_detect);
fprintf('detect_: took %0.1fs (fParfor=%d, fGpu=%d)\n', t_detect, P.fParfor, P.fGpu);
end %func


%--------------------------------------------------------------------------
function [viSite_spk, viTime_spk, vrAmp_spk, trPc_spk] = detect_merge_(cS_detect, viOffset_load)

vnSpk_load = cellfun(@(x)numel(x.viSite_spk), cS_detect);
miSpk_load = [0; cumsum(vnSpk_load)];
miSpk_load = [miSpk_load(1:end-1)+1, miSpk_load(2:end)];

nSpk = sum(vnSpk_load);
[viSite_spk, viTime_spk, vrAmp_spk] = ...
    deal(zeros(nSpk, 1, 'int32'), zeros(nSpk, 1, 'int64'), zeros(nSpk, 1, 'single'));
viOffset_load = int64(viOffset_load);
trPc_spk = [];

for iLoad = 1:numel(cS_detect)
    S1 = cS_detect{iLoad};
    viSpk1 = miSpk_load(iLoad,1):miSpk_load(iLoad,2);
    viSite_spk(viSpk1) = S1.viSite_spk;
    viTime_spk(viSpk1) = int64(S1.viTime_spk) + viOffset_load(iLoad);
    vrAmp_spk(viSpk1) = S1.vrAmp_spk;
    if isempty(trPc_spk)
        trPc_spk = zeros(size(S1.trPc_spk,1), size(S1.trPc_spk,2), nSpk, 'single');
    end
    trPc_spk(:,:,viSpk1) = S1.trPc_spk;
end
end %func


%--------------------------------------------------------------------------
function S_detect = detect_paged_(mrWav_T, P, S_cache)
if nargin<3, S_cache = []; end

% filter and trim 
% nlim_wav1 = struct_get_(S_cache, 'nlim_wav1');
mrWav2 = filter_(mrWav_T', P);
S_detect = get_spikes_(mrWav2, P, S_cache);

end %func


%--------------------------------------------------------------------------
function S_detect = get_spikes_(mrWav2, P, S_cache)
[vrThresh_site, nlim_wav1, mrPv_global] = ...
    struct_get_(S_cache, 'vrThresh_site', 'nlim_wav1', 'mrPv_global');

if isempty(nlim_wav1)
    [nPad_pre, nPad_post] = deal(0);
else
    nPad_pre = nlim_wav1(1)-1;
    nPad_post = size(mrWav2,1) - nlim_wav1(2);
end

% common mode rejection
if P.blank_thresh > 0
    if isempty(vnWav11)
        vnWav11 = mr2ref_(mrWav2, P.vcCommonRef, P.viSiteZero); %vrWav_mean1(:);    
    end
    vlKeep_ref = car_reject_(vnWav11(:), P);
    fprintf('Rejecting %0.3f %% of time due to motion\n', (1-mean(vlKeep_ref))*100 );
else
    vlKeep_ref = [];
end

% detect spikes or use the one passed from the input (importing)
if isempty(vrThresh_site), [vrThresh_site, fGpu] = mr2thresh_(mrWav2, P); end
[viTime_spk, vrAmp_spk, viSite_spk] = detect_spikes_(mrWav2, vrThresh_site, vlKeep_ref, P);
[viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)gather_(x), viTime_spk, vrAmp_spk, viSite_spk);    

% reject spikes within the overlap region
if ~isempty(nlim_wav1)
    viKeep_spk = find(viTime_spk >= nlim_wav1(1) & viTime_spk <= nlim_wav1(2));
    [viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)x(viKeep_spk), viTime_spk, vrAmp_spk, viSite_spk);    
end%if

% extract spike waveforms
trWav_spk = get_spkwav_(mrWav2, viSite_spk, viTime_spk, P);

% extract spike feaures
if isempty(mrPv_global)
    [mrPv_global, vrD_global] = get_pv_(trWav_spk, P); 
end
trPc_spk = gather_(project_pc_(trWav_spk, mrPv_global));
%     trWav_spk1 = pc2spkwav_(trPc_spk, mrPv_global); % recover waveform

% return struct
if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end
S_detect = makeStruct_(trPc_spk, mrPv_global, viTime_spk, vrAmp_spk, viSite_spk);
end %func


%--------------------------------------------------------------------------
function trWav_spk1 = pc2spkwav_(trPc_spk, mrPv_global)
[nPc_spk, nSites_spk, nSpk] = size(trPc_spk);
nSamples_spk = size(mrPv_global,1);
dimm_spk = [nSamples_spk, nSites_spk, nSpk];
trWav_spk1 = reshape(mrPv_global * reshape(trPc_spk, size(trPc_spk,1), []), dimm_spk);
end %func


%--------------------------------------------------------------------------
function trPc_spk = project_pc_(trWav_spk, mrPv)
[nSamples, nSites, nSpk] = size(trWav_spk);
nPc = size(mrPv,2);
mr = reshape(trWav_spk, size(trWav_spk,1), []);
trPc_spk = reshape((mrPv' * mr), nPc, nSites, nSpk);
end %func


%--------------------------------------------------------------------------
function [mrPv1, vrD1, fh_proj] = get_pv_(tr, P)
%tr: nSamples x nSpikes x nChans
MAX_SAMPLE = 10000;        

viSpk_sub = subsample_vr_(1:size(tr,3), MAX_SAMPLE);
switch 2
    case 2, mr1 = reshape(tr(:, :, viSpk_sub), size(tr,1), []); 
    case 1, mr1 = squeeze(tr(:, 1, viSpk_sub)); 
end
nPc_spk = get_set_(P, 'nPc_spk', 6); % # components to compress spike waveforms

switch 2
    case 1
        % mrSpkWav1 = meanSubt_(mrSpkWav1);
        [mrPv1, vrD1] = eig(mr1 * mr1');
        mrPv1 = zscore_(fliplr(mrPv1)); % sort largest first
        vrD1 = flipud(diag(vrD1));
    case 2
        [mrPv1, ~, vrD1] = pca(gather_(mr1)','Centered',false, 'NumComponents', nPc_spk);
%         fh_proj = @(x)(mrPv' * x)';
%         mr2 = mrPv * mrPc'; % confirm the agreement
%         mr3 = mrPv * (mrPv' * mr1);
end
% spike center should be negative
% iMid = 1-P.spkLim(1);
% vrSign = (mrPv1(iMid,:) < 0) * 2 - 1; %1 or -1 depending on the sign
% mrPv = bsxfun(@times, mrPv1, vrSign);
end %func


%--------------------------------------------------------------------------
function trWav_spk = get_spkwav_(mrWav, viSite_spk, viTime_spk, P)
nSpks = numel(viTime_spk);
nSites = numel(P.viSite2Chan);
spkLim_wav = P.spkLim;
% spkLim_raw = P.spkLim_raw;
nSites_spk = size(P.miSites,1); %(P.maxSite * 2) + 1;

% Realignment parameters
fRealign_spk = get_set_(P, 'fRealign_spk', 0); %0,1,2
if isempty(viSite_spk)
    trWav_spk = permute(mr2tr_(mrWav, spkLim_wav, viTime_spk), [1,3,2]);
else
    trWav_spk = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'like', mrWav);
    for iSite = 1:nSites
        viiSpk11 = find(viSite_spk == iSite);
        if isempty(viiSpk11), continue; end
        viTime_spk11 = viTime_spk(viiSpk11); %already sorted by time
        viSite11 = P.miSites(:,iSite);
        try
            tnWav_spk1 = mr2tr_(mrWav, spkLim_wav, viTime_spk11, viSite11);
            trWav_spk(:,:,viiSpk11) = permute(tnWav_spk1, [1,3,2]);
        catch % GPU failure
            disperr_('mn2tn_wav_: GPU failed'); 
        end
    end
end
% if 0 %10/19/2018 JJJ
%     trWav_spk = meanSubt_spk_(trWav_spk);
% end
end %func


%--------------------------------------------------------------------------
function mnWav2 = filter_(mnWav1, P)
%-----
% Filter
fprintf('\tFiltering spikes...\n\t'); t_filter = tic;
if get_set_(P, 'fSmooth_spatial', 0)
    mnWav1 = spatial_smooth_(mnWav1, P);
end
vcDataType_filter = get_set_(P, 'vcDataType_filter', 'single');
try    
    mnWav1_ = cast_(mnWav1, vcDataType_filter);
    [mnWav1_, P.fGpu] = gpuArray_(mnWav1_, P.fGpu);
    [mnWav2, vnWav11] = filt_car_(mnWav1_, P);    
catch % GPU failure
    P.fGpu = 0;
    mnWav1_ = cast_(mnWav1, vcDataType_filter);
    [mnWav2, vnWav11] = filt_car_(mnWav1_, P);    
end
end %func


%--------------------------------------------------------------------------
function [out, nlim, fDone] = readmda_paged_(P)
% manage overlap
persistent S_mda fid iLoad nLoads nSamples_load nSamples_last nSamples_pad vcFile

if nargin==1
    vcFile = P.vcFile;
    [S_mda, fid] = readmda_header_(vcFile);
    S_mda.dimm = S_mda.dimm(:)'; % formatting
    assert(S_mda.nBytes_missing<1, 'readmda_paged_: mda file is incomplete');
    iLoad = 0;
    [nLoads, nSamples_load, nSamples_last] = plan_load_(S_mda.nBytes_data, P);
    nSamples_pad = get_set_(P, 'nPad_filt', 100);
    viOffset_load = (0:nLoads-1) * nSamples_load; % return offset
    S_page = makeStruct_(nLoads, nSamples_load, nSamples_last, nSamples_pad, viOffset_load);
    [out, nlim, fDone] = deal(S_page, [], 0);    
else
    if iLoad <= nLoads
        iLoad = iLoad + 1;
        if nLoads == 1 % no padding if single load
            dimm1 = S_mda.dimm;
            nlim = [1, dimm1(end)]; 
        elseif iLoad == 1 % first
            dimm1 = [S_mda.dimm(1:end-1), nSamples_load + nSamples_pad];
            nlim = [1, nSamples_load];
        elseif iLoad == nLoads % last
            dimm1 = [S_mda.dimm(1:end-1), nSamples_last + nSamples_pad];
            nlim = [1, nSamples_last] + nSamples_pad;
        else % middle (double padded)
            dimm1 = [S_mda.dimm(1:end-1), nSamples_load + nSamples_pad*2];
            nlim = [1, nSamples_load] + nSamples_pad;
        end
        t1=tic;
        out = fread_(fid, dimm1, S_mda.vcDataType);  % transpose MDA
        t_dur1 = toc(t1);
        mb_loaded1 = prod(dimm1) * bytesPerSample_(S_mda.vcDataType) / 1e6;
        fprintf('Read %s (%d/%d), took %0.1fs (%0.1f MB/s, %0.1f MB)\n', ...
            vcFile, iLoad, nLoads, t_dur1, mb_loaded1/t_dur1, mb_loaded1);
        if iLoad == nLoads
            fclose(fid);
            fDone = 1;
        else
            fDone = 0;
            frewind_(fid, [S_mda.dimm(1:end-1), nSamples_pad*2], S_mda.vcDataType);
        end
    else
        [out, nlim] = deal([]);
        fDone = 1;
    end
end

end %func


%--------------------------------------------------------------------------
function vcDir_in = get_test_data_()
if ispc()
    vcDir_in = 'C:\tmp\groundtruth\hybrid_synth\static_siprobe\rec_64c_1200s_11'; 
elseif isunix()
    vcDir_in = '~/ceph/groundtruth/hybrid_synth/static_siprobe/rec_64c_1200s_11'; 
end
end %func


%--------------------------------------------------------------------------
function P = makeParam_(vcDir_in, vcDir_out, vcFile_arg)
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end

% assume there is raw.mda, geom.csv, params.json, firings_true.mda
P = file2struct_(ircpath_(read_cfg_('default_prm', 0)));  %P = defaultParam();
P.vcFile = fullfile(vcDir_in, 'raw.mda');
P.vcDir_out = vcDir_out;

% probe file
P.mrSiteXY = csvread(fullfile(vcDir_in, 'geom.csv'));
P.viSite2Chan = 1:size(P.mrSiteXY,1);

% load json file
S_json = loadjson_(fullfile(vcDir_in, 'params.json'));
P.sRateHz = get_set_(S_json, 'samplerate', P.sRateHz);
P.fInverse_file = get_set_(S_json, 'spike_sign', -1) == -1;

% read param
if isstruct(vcFile_arg)
    S_arg = vcFile_arg;
elseif exist_file_(vcFile_arg)
    S_arg = meta2struct_(vcFile_arg);
else
    S_arg = [];
end
P = struct_merge_(P, S_arg);

% derived fields
P = fill_param_(P);
end %func


%--------------------------------------------------------------------------
function P = fill_param_(P)
P.spkRefrac = round(P.spkRefrac_ms * P.sRateHz / 1000);
P.spkLim = round(P.spkLim_ms * P.sRateHz / 1000);
% P.spkLim_raw = calc_spkLim_raw_(P);
if ~isempty(get_(P, 'viChanZero')) && isempty(P.viSiteZero)
    [~, viSiteZero] = ismember(P.viChanZero, P.viSite2Chan);
    P.viSiteZero = viSiteZero(viSiteZero>0);
end
if ~isempty(get_(P, 'viSiteZero')), P.viSiteZero(P.viSiteZero > numel(P.viSite2Chan)) = []; end
if ~isfield(P, 'viShank_site'), P.viShank_site = []; end
[P.miSites, P.nSites_fet] = findNearSites_(P.mrSiteXY, P);
% LFP sampling rate
if ~isempty(get_(P, 'nSkip_lfp'))
    P.sRateHz_lfp = P.sRateHz / P.nSkip_lfp;
else
    P.sRateHz_lfp = get_set_(P, 'sRateHz_lfp', 2500);
    P.nSkip_lfp = round(P.sRateHz / P.sRateHz_lfp);
end
if isempty(get_(P, 'nTime_clu'))
    try
        batch_sec_drift = get_set_(P, 'batch_sec_drift', 300);
        P.nTime_clu = max(round(recording_duration_(P) / batch_sec_drift), 1);
        P.nTime_clu = min(P.nTime_clu, get_set_(P, 'nBatch_max_drift', 32));
        if fEditFile
            fprintf('\tnTime_clu = %d (batch_sec_drift = %0.1f s)\n', P.nTime_clu, batch_sec_drift);
        end
    catch
        P.nTime_clu = 1;
    end
end
if isempty(get_(P, 'nTime_drift'))
    try
        step_sec_drift = get_set_(P, 'step_sec_drift', 10);
        P.nTime_drift = max(round(recording_duration_(P) / step_sec_drift), 1);
        if fEditFile
            fprintf('\tnTime_drift = %d (step_sec_drift = %0.1f s)\n', P.nTime_drift, step_sec_drift);
        end
    catch
        P.nTime_drift = 64;
    end
end
end %func


%--------------------------------------------------------------------------
% from irc.m (v4.9.11)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 8/22/18 JJJ: changed from the cell output to varargout
% 9/26/17 JJJ: Created and tested
function varargout = struct_get_(varargin)
% Obtain a member of struct
% cvr = cell(size(varargin));
% if varargin is given as cell output is also cell
S = varargin{1};
for iArg=1:nargout
    vcName = varargin{iArg+1};
    if iscell(vcName)
        csName_ = vcName;
        cell_ = cell(size(csName_));
        for iCell = 1:numel(csName_)
            vcName_ = csName_{iCell};
            if isfield(S, vcName_)
                cell_{iCell} = S.(vcName_);
            end
        end %for
        varargout{iArg} = cell_;
    elseif ischar(vcName)
        if isfield(S, vcName)
            varargout{iArg} = S.(vcName);
        else
            varargout{iArg} = [];
        end
    else
        varargout{iArg} = [];
    end
end %for
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
function varargout = multifun_(hFun, varargin)
% apply same function to the input, unary function only

if nargout ~= numel(varargin), error('n arg mismatch'); end
for i=1:nargout
    try
        varargout{i} = hFun(varargin{i});
    catch
        varargout{i} = varargin{i};
    end
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
function frewind_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function disperr_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end

function out1 = meta2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_merge_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ircpath_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = read_cfg_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = loadjson_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = filesize_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = bytesPerSample_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = fread_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gather_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2ref_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = car_reject_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_copy_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cast_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2tr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsample_vr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = zscore_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

function [out1, out2] = readmda_header_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = mr2thresh_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = gpuArray_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = filt_car_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = findNearSites_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end

function [out1, out2, out3] = plan_load_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = detect_spikes_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end

