% works on mda format only for now
% no manual interface now
% analyzes single mda file format
% todo: multiple mda files
% todo: multiple formats

%--------------------------------------------------------------------------
function varargout = irc2(vcDir_in, vcDir_out, vcFile_arg, vcArg3)
% irc2(vcDir_in, vcDir_out, vcFile_arg)
% irc2(vcCmd, vcArg1, vcArg2)

if nargin<1, vcDir_in = ''; end
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end
if nargin<4, vcArg3 = ''; end

persistent vcFile_prm_

% batch processing. it uses default param for now
if iscell(vcDir_in) && iscell(vcDir_out)
    [csDir_in, csDir_out, csFile_arg] = deal(vcDir_in, vcDir_out, vcFile_arg);
    parfor iFile = 1:numel(csDir_in)
        try
            fprintf('irc2 batch-processing %s (%d/%d)\n', csDir_in{iFile}, iFile, numel(csDir_in));
            irc2(csDir_in{iFile}, csDir_out{iFile}, csFile_arg{iFile});
        catch
            disp(lasterr());
        end
    end 
    return;
end

[vcCmd, vcArg1, vcArg2] = deal(vcDir_in, vcDir_out, vcFile_arg); 
if isempty(vcFile_arg)
    vcFile_arg = file2struct_(read_cfg_('default2_prm'));
end

P = []; 
switch lower(vcCmd)
    case 'juxta'
        convert_mda_ui('english'); return;
    case 'version'
        if nargout==0, version_(); 
        else, varargout{1} = version_(); 
        end
        return;
    case 'scoreboard', irc2_scoreboard(); return;
    case {'detect-sort', 'sort', 'auto', '', 'describe', 'verify', 'manual'}
        if isempty(vcDir_out)
            vcFile_prm = vcFile_prm_;
            fprintf('irc2 (%s) opening %s\n', version_(), vcFile_prm);
        else
            vcFile_prm = vcDir_out;
            vcFile_prm_ = vcFile_prm;
        end
        if isempty(vcFile_prm), fprintf(2, 'provide .prm file.\n'); return; end
        if ~exist_file_(vcFile_prm)
            fprintf(2, 'File does not exist: %s\n', vcFile_prm); return;
        end
        P = file2struct_(vcFile_prm);
        S0 = load0_(P.vcFile_prm);
        fPlot_gt = [];
        switch lower(vcCmd)
            case 'detect-sort', clear_();
            case 'sort', clear_('sort');
            case 'describe', describe_(S0); return;
            case {'verify', 'validate'}, validate_(P, fPlot_gt); return;
            case 'manual', manual_(P); return;
        end
    case 'benchmark'
        if nargout==0, benchmark_(vcArg1, vcArg2, vcArg3); 
        else, varargout{1} = benchmark_(vcArg1, vcArg2, vcArg3); 
        end
        return;
    case 'plot', irc('plot', vcArg1, vcArg2); return;
    case 'clear', clear_(); vcFile_prm_=[]; return;
    case 'clear-sort', clear_('sort'); return;        
    case {'test-static', 'test-drift', 'test-tetrode', 'test-tetrode2', 'test-tetrode3', ...
            'test-bionet', 'test-bionet1', 'test-monotrode', ...
            'test-monotrode1', 'test-monotrode2', 'test-monotrode3'}
        vcDir_in = get_test_data_(strsplit_get_(vcCmd,'-',2));
        fPlot_gt = [];
    case 'export', irc('export', vcArg1); return;
    otherwise
        fPlot_gt=0;
        clear_();
end

fprintf('Running irc2.m (%s)\n', version_());

if isempty(P)
    P = makeParam_(vcDir_in, vcDir_out, vcFile_arg);
end
vcFile_prm_ = P.vcFile_prm;
S0 = get(0, 'UserData');

if isempty(struct_get_(S0, 'trPc_spk'))
    S0 = detect_(P); 
    set(0, 'UserData', S0);
end
if isempty(struct_get_(S0, 'S_clu'))
    S0 = sort_(S0, P);
    set(0, 'UserData', S0);
end
S0 = auto_(S0, P);
set(0, 'UserData', S0);
describe_(S0);

% output
save0_(S0);
% struct_save_(S0, strrep(P.vcFile_prm, '.prm', '_irc.mat'));
vcFile_firings_mda = fullfile(P.vcDir_out, 'firings.mda');
save_firings_mda_(S0, vcFile_firings_mda);

% Validate
if get_set_(P, 'fValidate', 1) && ~strcmpi(vcCmd, 'auto')
    validate_(P, fPlot_gt);
end
end %func


%--------------------------------------------------------------------------
function manual_(P)
global fDebug_ui


if ~is_sorted_(P)
    fprintf(2, 'File must to be sorted first (run "irc spikesort %s")\n', P.vcFile_prm); 
    return; 
end
S0 = load0_(P.vcFile_prm);
P.fParfor = 0;
fDebug_ui = false;

hMsg = msgbox_('Plotting... (this closes automatically)'); t1=tic;
S0 = figures_manual_(P); %create figures for manual interface
clear mouse_figure;
clear get_fig_cache_ get_tag_ %clear persistent figure handles

% Set fields
S0 = struct_merge_(S0, ...
    struct('iCluCopy', 1, 'iCluPaste', [], 'hCopy', [], 'hPaste', [], 'nSites', numel(P.viSite2Chan)));
set(0, 'UserData', S0);

% hFigRD
try
    S0.S_clu = plot_FigRD_(S0.S_clu, P); % ask user before doing so
    set(0, 'UserData', S0);
catch
    ;
end

% Set initial amplitudes
plot_FigWavCor_(S0);  % need `S_clu.mrWavCor`
S0 = plot_FigWav_(S0); % hFigWav %do this after for ordering

close_(get_fig_('FigTrial')); %close previous FigTrial figure
close_(get_fig_('FigTrial_b')); %close previous FigTrial figure
S0 = button_CluWav_simulate_(1, [], S0); %select first clu
auto_scale_proj_time_(S0);
S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
S_log = load_(strrep(P.vcFile_prm, '.prm', '_log.mat'), [], 0);
if ~isempty(S_log), S0.cS_log = {S_log}; end
save_log_('start', S0); %crash proof log

% Finish up
close_(hMsg);
fprintf('UI creation took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function flag = is_sorted_(P)
% return true if already detected. .spkwav file must exist

vcFile_mat = subsFileExt_(P.vcFile_prm, '_irc.mat');
if exist_file_(vcFile_mat)
    csNames = whos('-file', vcFile_mat);
    flag = ismember('S_clu', {csNames.name});    
else
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
function S_bench = benchmark_(vcDir_in1, vcDir_out1, vcParam1)
if nargin<3, vcParam1 = ''; end

% specify the output folder
vcFile_mat1 = dir_(fullfile(vcDir_out1, '*_irc.mat'));
if isempty(vcFile_mat1)
    % process the data
    fprintf('Running benchmark: ''%s'' using ''%s'': ', vcDir_in1, vcParam1); t1=tic;
    [~, vcConsoleOut] = system(sprintf('./run_irc %s %s %s', vcDir_in1, vcDir_out1, vcParam1));
    fprintf('took %0.1fs\n', toc(t1));
else
    % load already processed data
    vcFile_prm1 = strrep(vcFile_mat1{1}, '_irc.mat', '.prm');
    vcConsoleOut = evalc(sprintf('irc2(''describe'', ''%s'');', vcFile_prm1));
    fprintf('Loaded benchmark from %s\n', vcFile_mat1{1});
end
% parse the output
[memory_gb, vcFile_prm, runtime_sec, runtime_detect_sec, runtime_sort_sec, runtime_merge_sec] = ...
    parse_console_out_(vcConsoleOut, ...
        'memory usage (GiB):', 'Parameter file:', 'Total runtime (s):', ...
        'Detect + feature (s):', 'Cluster runtime (s):', 'merge runtime (s):');
% str2num_strip_ = @(x)str2double(x((x>='0' & x<='9') | x=='.'));
S_bench = makeStruct_func_(@(x)str2num_(x), memory_gb, runtime_sec, runtime_detect_sec, runtime_sort_sec, runtime_merge_sec);
S_bench = struct_add_(S_bench, vcFile_prm, vcConsoleOut);

end %func


%--------------------------------------------------------------------------
function S = makeStruct_func_(varargin)
% pass the function handle and add the field to the struct after applying a
% function

S = struct();
fh = varargin{1};
for i = 2:nargin
    S.(inputname(i)) =  fh(varargin{i}); 
end
end %func


%--------------------------------------------------------------------------
% 11/28/2018 JJJ: string to number
function varargout = str2num_(varargin)
for iArg = 1:nargin
    val1 = varargin{iArg};
    if ischar(val1)
        % strip non numeric characters
        val1 = val1((val1>='0' & val1<='9') | val1=='.');
        varargout{iArg} = str2num(val1);
    else
        varargout{iArg} = nan;
    end
end
end %func


%--------------------------------------------------------------------------
% returns string for a given argument
function varargout = parse_console_out_(varargin)
% [val1, val2, ...] = parse_system_out_(system_out, name1, name2, ...)

cs1 = strsplit(varargin{1}, '\n');
for iArg_out = 1:nargout
    iLine1 = find(contains(cs1, varargin{iArg_out+1}), 1, 'first');
    if ~isempty(iLine1)
        cs2 = strsplit(cs1{iLine1}, ' ');
        varargout{iArg_out} = cs2{end};
    else
        varargout{iArg_out} = [];
    end
end %for
end %func


%--------------------------------------------------------------------------
% 11/6/18 JJJ: Displaying the version number of the program and what's used. #Tested
function [vcVer, vcDate, vcHash] = version_()
vcVer = 'v5.0.4';
vcDate = '10/9/2019';
vcHash = file2hash_();

if nargout==0
    fprintf('%s (%s) installed, MD5: %s\n', vcVer, vcDate, vcHash);
    return;
end
end %func


%--------------------------------------------------------------------------
function validate_(P, fPlot_gt)
if nargin<2, fPlot_gt = []; end
if isempty(fPlot_gt), fPlot_gt = read_cfg_('fPlot_gt'); end

vcFile_gt_mda = get_(P, 'vcFile_gt');
if ~exist_file_(vcFile_gt_mda), return; end
vcFile_firings_mda = fullfile(P.vcDir_out, 'firings.mda');
S_score = irc('validate-mda', P.vcFile_gt, vcFile_firings_mda, P.vcFile, fPlot_gt); % assume that groundtruth file exists
struct_save_(S_score, fullfile(P.vcDir_out, 'raw_geom_score.mat'), 1);

end %func


%--------------------------------------------------------------------------
function save0_(S0)
if nargin<1, S0 = []; end
if isempty(S0), S0 = get(0, 'UserData'); end
vcFile_mat = strrep(S0.P.vcFile_prm, '.prm', '_irc.mat');
vcFile_fet = strrep(S0.P.vcFile_prm, '.prm', '_fet.irc');
vcFile_fet2 = strrep(S0.P.vcFile_prm, '.prm', '_fet2.irc');
vcFile_knn = strrep(S0.P.vcFile_prm, '.prm', '_knn.irc');

trPc_spk = gather_(get_(S0, 'trPc_spk'));
if ~isempty(trPc_spk)
    S0.trPc_spk = [];
    S0.dimm_fet = size(trPc_spk);
    S0.type_fet = class(trPc_spk);
    write_bin_(vcFile_fet, trPc_spk);
end

trPc2_spk = gather_(get_(S0, 'trPc2_spk'));
if ~isempty(trPc2_spk)
    S0.trPc2_spk = [];
    write_bin_(vcFile_fet2, trPc2_spk);
end

S_clu = get_(S0, 'S_clu');
miKnn = gather_(get_(S_clu, 'miKnn'));
if ~isempty(miKnn)
    S0.S_clu.miKnn = [];
    S0.dimm_knn = size(miKnn);
    S0.type_knn = class(miKnn);
    write_bin_(vcFile_knn, miKnn);
end

struct_save_(S0, vcFile_mat);
end %


%--------------------------------------------------------------------------
function S0 = load0_(vcFile_prm)
vcFile_mat = strrep(vcFile_prm, '.prm', '_irc.mat');
vcFile_fet = strrep(vcFile_prm, '.prm', '_fet.irc');
vcFile_fet2 = strrep(vcFile_prm, '.prm', '_fet2.irc');
vcFile_knn = strrep(vcFile_prm, '.prm', '_knn.irc');

S0 = load(vcFile_mat);
if isempty(get_(S0, 'trPc_spk'))
    if exist_file_(vcFile_fet)
        S0.trPc_spk = load_bin_(vcFile_fet, S0.type_fet, S0.dimm_fet);
    end
end
if isempty(get_(S0, 'trPc2_spk'))
    if exist_file_(vcFile_fet2)
        S0.trPc2_spk = load_bin_(vcFile_fet2, S0.type_fet, S0.dimm_fet);
    end
end

S_clu = get_(S0, 'S_clu');
if isempty(get_(S_clu, 'miKnn'))
    if exist_file_(vcFile_knn)
        S0.S_clu.miKnn = load_bin_(vcFile_knn, S0.type_knn, S0.dimm_knn);
    end
end
set(0, 'UserData', S0);
end %func


%--------------------------------------------------------------------------
% negative index means from the end, 0 index means end index
function vc1 = strsplit_get_(vc,delim,idx)
cs = strsplit(vc, delim);
idx = mod(idx-1,numel(cs))+1;
vc1 = cs{idx};
end %func


%--------------------------------------------------------------------------
function S0 = sort_(S0, P)
switch get_set_(P, 'sort_mode', 1)
    case 1, S0 = sort1_(S0, P);
    case 2, S0 = sort2_(S0, P);    
end %switch
end %func


%--------------------------------------------------------------------------
function vcDir_in = clear_(vcMode)
if nargin<1, vcMode = 'all'; end

switch vcMode
    case 'all'
        set(0, 'UserData', []);
        disp('irc2: cleared all');
        vcDir_in = '';
    case 'sort'
        S0 = get(0, 'UserData');
        S0.S_clu = [];
        set(0,'UserData',S0);
        disp('irc2: cleared sort');
        try
            vcDir_in = fileparts(S0.P.vcFile);
        catch
            vcDir_in = '';
        end
end
end

%--------------------------------------------------------------------------
function csDesc = describe_(S0)
P=S0.P;

runtime_total = S0.runtime_detect + S0.runtime_sort + S0.runtime_automerge;
tDur = recording_duration_(S0.P, S0); 
memory_sort = S0.memory_sort - S0.memory_init;
memory_detect = S0.memory_detect - S0.memory_init;
nSites = numel(P.viSite2Chan);
nShanks = max(P.viShank_site);
nSitesPerEvent = size(P.miSites,1);
nSpk = numel(S0.viTime_spk);
nFeatures = P.nSites_fet * P.nPcPerChan;

csDesc = {};
try
    csDesc = {};
    csDesc{end+1} = sprintf('Recording format');
    csDesc{end+1} = sprintf('    Recording file:         %s', P.vcFile);
    csDesc{end+1} = sprintf('    Probe file:             %s', P.probe_file);
    csDesc{end+1} = sprintf('    Recording Duration:     %0.1fs', tDur);
    csDesc{end+1} = sprintf('    Data Type:              %s', P.vcDataType);
    csDesc{end+1} = sprintf('    #Channels in file:      %d', P.nChans);
    csDesc{end+1} = sprintf('    #Sites:                 %d', nSites);
    csDesc{end+1} = sprintf('    #Shanks:                %d', nShanks);
    csDesc{end+1} = sprintf('Pre-processing');
    csDesc{end+1} = sprintf('    Filter type:            %s', P.vcFilter);
    csDesc{end+1} = sprintf('    Filter range (Hz):      %0.1f-%0.1f', P.freqLim);
    csDesc{end+1} = sprintf('    Common ref:             %s', P.vcCommonRef);
    csDesc{end+1} = sprintf('    FFT threshold:          %d', get_set_(P, 'fft_thresh', 0));
    csDesc{end+1} = sprintf('Events');
    csDesc{end+1} = sprintf('    #Spikes:                %d', nSpk);
    csDesc{end+1} = sprintf('    Feature extracted:      %s', P.vcFet);
    csDesc{end+1} = sprintf('    #Sites/event:           %d', nSitesPerEvent);
    csDesc{end+1} = sprintf('    maxDist_site_um:        %0.0f', P.maxDist_site_um);    
    csDesc{end+1} = sprintf('    maxDist_site_spk_um:    %0.0f', P.maxDist_site_spk_um);    
    csDesc{end+1} = sprintf('    #Features/event:        %d', nFeatures);    
catch
    ;
end
if isfield(S0, 'S_clu')
    S_clu = S0.S_clu;
    csDesc{end+1} = sprintf('Cluster');    
    csDesc{end+1} = sprintf('    #Clusters:              %d', S_clu.nClu);
    csDesc{end+1} = sprintf('    #Unique events:         %d', sum(S_clu.viClu>0));
    csDesc{end+1} = sprintf('    min. spk/clu:           %d', P.min_count);
    csDesc{end+1} = sprintf('    Cluster method:         %s', P.vcCluster);
    csDesc{end+1} = sprintf('    knn:                    %d', P.knn);
    csDesc{end+1} = sprintf('    nTime_clu:              %d', P.nTime_clu);
    csDesc{end+1} = sprintf('    nTime_drift:            %d', P.nTime_drift);
    csDesc{end+1} = sprintf('    fSpatialMask_clu:       %d', P.fSpatialMask_clu);
    csDesc{end+1} = sprintf('Auto-merge');   
    csDesc{end+1} = sprintf('    delta_cut:              %0.3f', get_set_(P, 'delta_cut', 1));
    csDesc{end+1} = sprintf('    maxWavCor:              %0.3f', P.maxWavCor);
end
try
    csDesc{end+1} = sprintf('Runtime (s)');
    csDesc{end+1} = sprintf('    Detect + feature (s):   %0.1fs', S0.runtime_detect);    
    csDesc{end+1} = sprintf('    Cluster runtime (s):    %0.1fs', S0.runtime_sort);
    csDesc{end+1} = sprintf('    merge runtime (s):      %0.1fs', S0.runtime_automerge);
    csDesc{end+1} = sprintf('    Total runtime (s):      %0.1fs', runtime_total);
    csDesc{end+1} = sprintf('    Runtime speed:          x%0.1f realtime', tDur / runtime_total);    

    csDesc{end+1} = sprintf('memory usage (GiB):         %0.3f', max(memory_detect,memory_sort)/2^30);
    csDesc{end+1} = sprintf('    detect(GiB):            %0.3f', memory_detect/2^30);
    csDesc{end+1} = sprintf('    sort(GiB):              %0.3f', memory_sort/2^30);

    csDesc{end+1} = sprintf('Execution');
    csDesc{end+1} = sprintf('    fGpu (GPU use):         %d', P.fGpu);
    csDesc{end+1} = sprintf('    fParfor (parfor use):   %d', P.fParfor);
    csDesc{end+1} = sprintf('    Parameter file:         %s', P.vcFile_prm);
catch
    ;
end
if nargout==0
    cellfun(@(x)disp(x), csDesc);
end
end %func


%--------------------------------------------------------------------------
function save_firings_mda_(S0, vcFile_firings_mda)
S_gt1 = struct('viTime', S0.viTime_spk, 'viSite', S0.viSite_spk, 'viClu', S0.S_clu.viClu);
gt2mda_(S_gt1, vcFile_firings_mda);
end %func
            

%--------------------------------------------------------------------------
% auto merge
function S0 = auto_(S0, P)
if nargin<2, P=[]; end
if nargin<1, S0=[]; end
if isempty(S0)
    S0 = get(0, 'UserData');
elseif isempty(P)
    P = S0.P;
end

fprintf('\nauto-merging...\n'); runtime_automerge=tic;


% refresh clu, start with fundamentals
S0.S_clu = struct_copy_(S0.S_clu, ...
    'rho', 'delta', 'ordrho', 'nneigh', 'P', 'miKnn', 'S_drift');
S0.S_clu = postCluster_(S0.S_clu, P); % peak merging

maxWavCor = get_set_(P, 'maxWavCor', .99);
if maxWavCor<1
    post_merge_mode = get_set_(P, 'post_merge_mode', 1);
    switch post_merge_mode 
        case 1, mrDist_clu = wave_similarity_clu1_(S0, P);
        case 2, mrDist_clu = calc_dist_ccg(S0, P);
        otherwise, mrDist_clu = calc_dist_clu_(S0, P, post_merge_mode);
    end
    S0.S_clu = templateMatch_post_(S0.S_clu, P, mrDist_clu);
end

% compute SNR per cluster and remove small SNR
S0.S_clu = S_clu_sort_(S0.S_clu, 'viSite_clu');
S0.S_clu = S_clu_refrac_(S0.S_clu, P); % refractory violation removal
S0.S_clu.mrCC = correlogram_(S0.S_clu, get0('viTime_spk'), P);
S0.S_clu = calc_clu_wav_(S0, P);

S0.runtime_automerge = toc(runtime_automerge);
fprintf('\n\tauto-merging took %0.1fs (fGpu=%d, fParfor=%d)\n', ...
    S0.runtime_automerge, P.fGpu, P.fParfor);
end %func


%--------------------------------------------------------------------------
function [S_clu, viClu_delete] = calc_clu_wav_(S0, P)

S_clu = S0.S_clu;
S_clu.nClu = max(S_clu.viClu);
cviSpk_clu = arrayfun_(@(x)find(S_clu.viClu==x), (1:S_clu.nClu)');
viSite_clu = cellfun(@(x)mode(S0.viSite_spk(x)), cviSpk_clu);
nT = size(S0.mrPv_global,1);
nSites_spk = size(S0.trPc_spk,2);
nSites = numel(P.viSite2Chan);
trWav_clu = zeros(nT, nSites_spk, S_clu.nClu, 'single');
tmrWav_clu = zeros(nT, nSites, S_clu.nClu, 'single'); %global waveform
for iClu = 1:S_clu.nClu
    viSpk1 = cviSpk_clu{iClu};
    if isempty(viSpk1), continue; end
    viSite1 = S0.viSite_spk(viSpk1);
    iSite1 = viSite_clu(iClu);
    trPc1 = S0.trPc_spk(:,:,viSpk1);
    mrPc1 = mean(trPc1(:,:,viSite1==iSite1),3);
    mrWav1 = S0.mrPv_global * mrPc1;
    trWav_clu(:,:,iClu) = mrWav1;
    viSite_clu1 = P.miSites(:,iSite1);
    tmrWav_clu(:,viSite_clu1,iClu) = mrWav1;
end
[trWav_raw_clu, tmrWav_raw_clu] = deal(trWav_clu, tmrWav_clu); % for now raw is not saved
S_clu = struct_add_(S_clu, trWav_clu, trWav_raw_clu, tmrWav_clu, viSite_clu, tmrWav_raw_clu);

% compute SNR
mrWav1_clu = squeeze_(trWav_clu(:,1,:),2); 
vrVmax_clu = max(mrWav1_clu)';
vrVmin_clu = min(mrWav1_clu)';
vrVpp_clu = vrVmax_clu - vrVmin_clu;
vrRms_site = S0.vrThresh_site(:) / S0.P.qqFactor;
S_clu.vrSnr_clu = abs(vrVmin_clu(:)) ./ vrRms_site(viSite_clu);
S_clu.vrSnr2_clu = abs(vrVpp_clu(:)) ./ vrRms_site(viSite_clu);

if ~isempty(get_(P, 'min_snr_clu'))
    viClu_delete = find(S_clu.vrSnr_clu < P.min_snr_clu);
    if ~isempty(viClu_delete)
        nClu_pre = S_clu.nClu;
        S_clu = delete_clu_(S_clu, viClu_delete);
        nClu_post = S_clu.nClu;
        fprintf('calc_clu_wav_: %d->%d clusters, %d removed below SNR=%0.1f\n', ...
            nClu_pre, nClu_post, numel(viClu_delete), P.min_snr_clu);
    end
else
    viClu_delete = [];
end

% update similarity
S0.S_clu = S_clu;
mrWavCor = wave_similarity_clu1_(S0, P);
S_clu.mrWavCor = mrWavCor + mrWavCor'; % make it symmetric

S_clu = S_clu_position_(S_clu);
S_clu = S_clu_refresh_(S_clu);
if ~isfield(S_clu, 'csNote_clu'), S_clu.csNote_clu = cell(S_clu.nClu, 1); end
end %func


%--------------------------------------------------------------------------
function mrDist_clu = calc_dist_clu__(S0, P, fMode)
% output
% -----
% mrDist_clu: 0..1, 1 being most similar, 0 being most dissimilar
S_clu = S0.S_clu;
MIN_COUNT = P.min_count;

fprintf('calc_dist_clu_: \n\t'); t1=tic;

% compute average waveforms by clusters
nClu = max(S_clu.viClu);
nSites = max(S0.viSite_spk);
if 0
    trPc_spk = S0.trPc_spk(:,1:P.nSites_fet,:);
else
    trPc_spk = S0.trPc_spk;
end
nPc = size(trPc_spk,1);
% compute average clu waveforms (sort by centers)
[viClu, vrPos_spk, mrPos_spk] = deal(S_clu.viClu, S0.vrPow_spk, S0.mrPos_spk);
[nPc_spk, nSites_spk, nSpk] = size(S0.trPc_spk);
[cviSpk_drift, mlDrift, viDrift_spk] = struct_get_(S0.S_clu.S_drift, 'cviSpk_drift', 'mlDrift', 'viDrift_spk');
nDrift = numel(cviSpk_drift);

if isempty(cviSpk_drift)
    mrDist_clu = pos_dist_clu_(viClu, mrPos_spk, vrPos_spk, nClu, MIN_COUNT);
    return;
end


% compute PC per cluster
normalize_ = @(x)x ./ sqrt(sum(x.^2));
qrPc_clu_drift = zeros(nPc, nSites, nClu, nDrift, 'single');
miSite_clu_drift = zeros(nClu, nDrift);
for iDrift = 1:nDrift
    viSpk1 = find(viDrift_spk==iDrift);
    [viClu1, viSite1, trPc1] = deal(viClu(viSpk1), S0.viSite_spk(viSpk1), trPc_spk(:,:,viSpk1));     
    for iClu = 1:nClu
        vl11 = viClu1==iClu;
        if sum(vl11) < MIN_COUNT, continue; end
        viSite11 = viSite1(vl11);
        iSite11 = mode(viSite11);      
        trPc11 = trPc1(:,:,vl11);
        miSite_clu_drift(iClu, iDrift) = iSite11;
        switch 2
            case 1       
                vl11 = vl11 & viSite1==iSite11;
                if sum(vl11) < MIN_COUNT, continue; end
                viSite_pc_11 = P.miSites(:,iSite11);
                qrPc_clu_drift(:, viSite_pc_11, iClu, iDrift) = median3_(trPc11,3);                
            case 2
                trPc_full11 = trPc_full_(trPc11, viSite11, P.miSites, MIN_COUNT);
                qrPc_clu_drift(:, :, iClu, iDrift) = nanmedian(trPc_full11,3);
        end
    end
end
qrPc_clu_drift(isnan(qrPc_clu_drift)) = 0;
% qrPc_clu_drift = reshape(qrPc_clu_drift, [], nClu, nDrift);


mrDist_clu = zeros(nClu, 'single');
for iDrift = 1:nDrift
    viDrift1 = find(mlDrift(:,iDrift));
    qrPc_clu2 = qrPc_clu_drift(:,:,:,viDrift1);
    trPc_clu1 = qrPc_clu_drift(:,:,:,iDrift);
    for iClu1 = 1:nClu
        iSite11 = miSite_clu_drift(iClu1,iDrift);
        if iSite11==0, continue; end
        viSite11 = P.miSites(:,iSite11);
        vrPc11 = trPc_clu1(:,viSite11,iClu1);  vrPc11 = vrPc11(:);
        if all(vrPc11==0), continue; end        
        
        dimm12 = size(qrPc_clu2);
        trPc12 = reshape(qrPc_clu2(:,viSite11,:,:), [], dimm12(3), dimm12(4));        
        vrDist1 = min(sum((trPc12 - vrPc11).^2), [], 3) / sum(vrPc11.^2);
        vrDist1 = 1 - min(vrDist1(:), 1);
        mrDist_clu(:,iClu1) = max(mrDist_clu(:,iClu1), vrDist1);
    end
end %for

fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function trPc_full = trPc_full_(trPc, viSite, miSites, MIN_COUNT)

nSites = size(miSites,2);
[viSite1_uniq, ~, cviSpk1_uniq] = unique_count_(viSite);

trPc_full = nan(size(trPc,1), nSites, size(trPc,3), 'single');
for iUniq=1:numel(cviSpk1_uniq)
    viSpk1 = cviSpk1_uniq{iUniq};
    if numel(viSpk1) >= MIN_COUNT
        viSite1 = miSites(:, viSite1_uniq(iUniq));
        trPc_full(:,viSite1,viSpk1) = trPc(:,:,viSpk1);
    end
end
end %func


%--------------------------------------------------------------------------
function mrDist_clu = calc_dist_clu_(S0, P, fMode)
% output
% -----
% mrDist_clu: 0..1, 1 being most similar, 0 being most dissimilar
S_clu = S0.S_clu;
MIN_COUNT = P.min_count;

% compute average waveforms by clusters
nClu = max(S_clu.viClu);
nSites = max(S0.viSite_spk);
if 0
    trPc_spk = S0.trPc_spk(:,1:P.nSites_fet,:);
else
    trPc_spk = S0.trPc_spk;
end
% compute average clu waveforms (sort by centers)
[viClu, vrPos_spk, mrPos_spk] = deal(S_clu.viClu, S0.vrPow_spk, S0.mrPos_spk);
[nPc_spk, nSites_spk, nSpk] = size(S0.trPc_spk);
[cviSpk_drift, mlDrift] = struct_get_(S0.S_clu.S_drift, 'cviSpk_drift', 'mlDrift');
unique_vec_ = @(x)unique(x(:));
if P.fParfor
    gcp_ = gcp();
else
    gcp_ = [];
end

fprintf('calc_dist_clu_ (fParfor=%d): \n\t', P.fParfor); t1=tic;
if ~isempty(cviSpk_drift)
    mrDist_clu = zeros(nClu, 'single');
    cmrDist_clu = cell(size(cviSpk_drift));
    for iDrift = 1:numel(cviSpk_drift)
        switch 1
            case 5
                viSpk1 = S0.S_clu.miKnn(:,cviSpk_drift{iDrift});
                viSpk1 = unique([viSpk1(:); cell2mat_(cviSpk_drift(find(mlDrift(:, iDrift))))]);
            case 4
                viSpk1 = cell2mat_(cviSpk_drift(find(mlDrift(:, iDrift))));
                viSpk1 = unique_vec_(S0.S_clu.miKnn(:,viSpk1));
            case 3, viSpk1 = unique_vec_(S0.S_clu.miKnn(:,cviSpk_drift{iDrift}));
            case 2, viSpk1 = cviSpk_drift{iDrift};
            case 1, viSpk1 = cell2mat_(cviSpk_drift(find(mlDrift(:, iDrift))));
        end
        [viClu1, viSite1, trPc1] = deal(viClu(viSpk1), S0.viSite_spk(viSpk1), trPc_spk(:,:,viSpk1)); 
        if fMode==1, fMode = S0.mrPv_global; end
        if isempty(gcp_)
            cmrDist_clu{iDrift} = wav_dist_clu_(viClu1, viSite1, trPc1, fMode, nClu, P);
        else
            vS_out(iDrift) = parfeval(gcp_, ...
                @(x,y,z)wav_dist_clu_(x,y,z, fMode, nClu, P), 1, viClu1, viSite1, trPc1);
        end
        fprintf('.');
    end %for
    if ~isempty(gcp_)
        for iDrift1 = 1:numel(cviSpk_drift)
            [iDrift, mrDist_clu1] = fetchNext(vS_out);
            cmrDist_clu{iDrift} = mrDist_clu1;
        end
    end
    mrDist_clu = max(cat(3,cmrDist_clu{:}), [], 3);
else
    mrDist_clu = pos_dist_clu_(viClu, mrPos_spk, vrPos_spk, nClu, MIN_COUNT);
end
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function mrDist_clu = pos_dist_clu_(viClu, mrPos_spk, vrPow_spk, nClu, MIN_COUNT)

switch 2
    case 3, fh_pow_ = @(x)sqrt(x);
    case 2, fh_pow_ = @(x)log(x);
    case 1, fh_pow_ = @(x)x;
end

% compute average pc per cluster per site
cviSpk_clu = arrayfun(@(x)find(viClu==x), 1:nClu, 'UniformOutput', 0);
[mrPos_med_clu, mrPos_sd_clu] = deal(nan(size(mrPos_spk,2)+1, nClu));
for iClu = 1:nClu
    viSpk1 = cviSpk_clu{iClu};
    if numel(viSpk1) >= MIN_COUNT
        mrPos_spk1 = [mrPos_spk(viSpk1,:), fh_pow_(vrPow_spk(viSpk1))];
        [mrPos_sd_clu(:,iClu), mrPos_med_clu(:,iClu)] = sd_mad_(mrPos_spk1);
    end
end

% compare cluster by cluster
sd_global = nanmedian(mrPos_sd_clu,2);
mrDist_clu = zeros(nClu);
for iClu1 = 1:nClu    
    switch 2
        case 2, vrCorr12 = (mrPos_med_clu - mrPos_med_clu(:,iClu1)) ./ sd_global;  
        case 1
            [mu1, sd1] = deal(mrPos_med_clu(:,iClu1), mrPos_sd_clu(:,iClu1));
            vrCorr12 = (mrPos_med_clu - mu1) ./ sd1;  
    end
    mrDist_clu(:,iClu1) = sqrt(sum(vrCorr12.^2))';
end
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wav_dist_clu_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P)
switch 1
    case 2
        % this uses new method of site referencing
        mrDist_clu = wav_dist_clu2_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P);
    case 1
        mrDist_clu = wav_dist_clu1_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P);
end %switch
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wav_dist_clu2_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P)

% compute average waveforms by clusters
MIN_COUNT = P.min_count;
nSites = size(P.miSites,2);
mrDist_clu = zeros(nClu, 'single');

% site loop
cviSpk_site = arrayfun(@(x)find(viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
for iSite = 1:nSites
    viSpk1 = find(viSite_spk==iSite);
    [mrFet12, viSpk12, ~, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, [], P, iSite);
    if isempty(viSpk12), continue; end
    viClu12 = viClu(viSpk12);
    [viClu_uniq1, vnClu_uniq1, cviClu_uniq1] = unique_count_(viClu12);    
    [viClu_uniq1, vnClu_uniq1, cviClu_uniq1] = ...
        multifun_(@(x)x(vnClu_uniq1 >= MIN_COUNT & viClu_uniq1 ~= 0), viClu_uniq1, vnClu_uniq1, cviClu_uniq1);
    if numel(viClu_uniq1)<2, continue; end
    mrPc12 = cell2matfun_(@(x)median2_(mrFet12(:,x),2), cviClu_uniq1');
%     mrPc12 = cell2matfun_(@(x)mean(mrFet12(:,x),2), cviClu_uniq1');
    switch 2
        case 2
            mrDist1 = zeros(size(mrPc12,2), 'single');
            for iClu1 = 1:size(mrPc12,2)
                mrDist1(:,iClu1) = sum((mrPc12 - mrPc12(:,iClu1)).^2) ./ sum(mrPc12(:,iClu1).^2);
            end
            mrDist1 = 1 - min(mrDist1,1);
        case 1, mrDist1 = corr(mrPc12);
    end
    mrDist_clu(viClu_uniq1,viClu_uniq1) = max(mrDist1, mrDist_clu(viClu_uniq1,viClu_uniq1));
end %for
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wav_dist_clu1_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P)

% compute average waveforms by clusters
MIN_COUNT = P.min_count;
[ctrPc_mu_clu, ctrPc_sd_clu, cviSite_clu] = deal(cell(nClu,1));
dimm_pc = size(trPc_spk);

% compute average pc per cluster per site
for iClu = 1:nClu
    viSpk1 = find(viClu==iClu);
    if isempty(viSpk1), continue; end
    viSite1 = viSite_spk(viSpk1);   
    
    [viSite1_uniq1, vn_uniq1, cviSpk1_uniq1] = unique_count_(viSite1);    
    [viSite1_uniq1, cviSpk1_uniq1] = ...
        multifun_(@(x)x(vn_uniq1>=MIN_COUNT), viSite1_uniq1, cviSpk1_uniq1);
    if isempty(cviSpk1_uniq1), continue; end    
    cviSite_clu{iClu} = viSite1_uniq1;
    switch 2
        case 2
            ctrPc_mu_clu{iClu} = mean_tr_(trPc_spk(:,:,viSpk1), cviSpk1_uniq1);
        case 1
            [ctrPc_mu_clu{iClu}, ctrPc_sd_clu{iClu}] = med_rms_tr_(trPc_spk(:,:,viSpk1), cviSpk1_uniq1);
    end
end
% construct arrays and sites
[viSite_all, viClu_all] = cvr2vr_vi_(cviSite_clu);
trPc_all = cat(3, ctrPc_mu_clu{:});

% compare cluster by cluster
mrDist_clu = zeros(nClu);
for iClu1 = 1:nClu
    [trPc_mu1, trPc_sd1, viSite1] = deal(ctrPc_mu_clu{iClu1}, ctrPc_sd_clu{iClu1}, cviSite_clu{iClu1});
    vrDist_clu1 = zeros(nClu,1);
    for iiSite1 = 1:numel(viSite1)        
        vi_all2 = find(viSite_all == viSite1(iiSite1) & viClu_all ~= iClu1);
        if isempty(vi_all2), continue; end
        mrPc_mu1 = trPc_mu1(:,:,iiSite1);
        if ~isempty(trPc_sd1)
            mrPc_sd1 = trPc_sd1(:,:,iiSite1);
        else
            mrPc_sd1 = [];
        end
        viClu2 = viClu_all(vi_all2);
        trPc2 = trPc_all(:,:,vi_all2);
        vrCorr12 = wavcor_pc_(mrPc_mu1, trPc2, mrPv_global, mrPc_sd1);
        vrDist_clu1(viClu2) = max(vrDist_clu1(viClu2), vrCorr12(:));
    end
    mrDist_clu(:,iClu1) = vrDist_clu1;
end
end %func


%--------------------------------------------------------------------------
function trPc_site = mean_tr_(trPc, cviSpk_site)
if numel(cviSpk_site) == 1
    trPc_site = mean(trPc,3);
else        
    trPc1 = permute(trPc, [3,1,2]);
    trPc_site = cell2matfun_(@(x)mean(trPc1(x,:,:)), cviSpk_site);    
    trPc_site = permute(trPc_site, [2,3,1]);
end   
end %func


%--------------------------------------------------------------------------
function [trPc_med_site, trPc_sd_site] = med_rms_tr_(trPc, cviSpk_site)
if numel(cviSpk_site) == 1
    trPc_med_site = median3_(trPc,3);
    if nargout>=2
        trPc_sd_site = std(trPc,[],3);
    end
else        
    trPc1 = permute(trPc, [3,1,2]);
    trPc_med_site = cell2matfun_(@(x)median3_(trPc1(x,:,:)), cviSpk_site);    
    trPc_med_site = permute(trPc_med_site, [2,3,1]);
    if nargout>=2
        trPc_sd_site = cell2matfun_(@(x)std(trPc1(x,:,:)), cviSpk_site);
        trPc_sd_site = permute(trPc_sd_site, [2,3,1]);
    end
end   
end %func


%--------------------------------------------------------------------------
function B = median3_(A,idimm)
if nargin<2, idimm=1; end

A = sort(A,idimm);
imid = ceil(size(A,idimm)/2);
if idimm==1
    B = A(imid,:,:);
elseif idimm==2
    B = A(:,imid,:);
else
    B = A(:,:,imid);
end
end %func


%--------------------------------------------------------------------------
function mat = cell2matfun_(varargin)
mat = cell2mat(cellfun(varargin{1}, varargin{2:end}, 'UniformOutput', 0));
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wave_similarity_clu_(S0, P, mode_sim)
S_clu = S0.S_clu;
MIN_COUNT = P.min_count;
if nargin<3, mode_sim = []; end
if isempty(mode_sim), mode_sim = 4; end

% compute average waveforms by clusters
nClu = max(S_clu.viClu);
cviSpk_clu = arrayfun(@(x)find(S_clu.viClu==x), 1:nClu, 'UniformOutput', 0);
% compute average clu waveforms (sort by centers)
[ctrPc_clu, ctrPos_clu, cviSite_clu] = deal(cell(nClu,1));
dimm_pc = size(S0.trPc_spk); 
dimm_pc(2) = min(dimm_pc(2), P.nSites_fet);
dimm_pos = size(S0.mrPos_spk);
switch 2
    case 2, fh_pow_ = @(x)log(x);
    case 1, fh_pow_ = @(x)x;
end
% compute average pc per cluster per site
% find centroid of the cluster
for iClu = 1:nClu
    viSpk1 = cviSpk_clu{iClu};
    viSite1 = S0.viSite_spk(viSpk1);
    [viSite1_uniq1, vn_uniq1, cviSpk1_uniq1] = unique_count_(viSite1);
    trPc_spk1 = S0.trPc_spk(1:dimm_pc(1), 1:dimm_pc(2), viSpk1);
    mrPos_spk1 = [S0.mrPos_spk(viSpk1,:), fh_pow_(S0.vrPow_spk(viSpk1))];
    [viSite1_uniq1, cviSpk1_uniq1] = ...
        multifun_(@(x)x(vn_uniq1>=MIN_COUNT), viSite1_uniq1, cviSpk1_uniq1);
    trPc1 = zeros(dimm_pc(1), dimm_pc(2), numel(cviSpk1_uniq1), 'single');
    trPos1 = zeros(dimm_pos(2)+1, 2, numel(cviSpk1_uniq1), 'single');    
    for iiSite1 = 1:numel(cviSpk1_uniq1)
        % compute average waveforms by peak sites
        viSpk11 = cviSpk1_uniq1{iiSite1};
        trPc1(:,:,iiSite1) = mean(trPc_spk1(:,:,viSpk11), 3);
        [sd_, med_] = sd_mad_(mrPos_spk1(viSpk11,:));
        trPos1(:,:,iiSite1) = [med_(:), sd_(:)];
    end
    [ctrPc_clu{iClu}, ctrPos_clu{iClu}, cviSite_clu{iClu}] = deal(trPc1, trPos1, viSite1_uniq1);
end
% construct arrays and sites
[viSite_all, viClu_all] = cvr2vr_vi_(cviSite_clu);
[trPc_all, trPos_all] = deal(cat(3, ctrPc_clu{:}), cat(3, ctrPos_clu{:}));

% compare cluster by cluster
mrPv_global = S0.mrPv_global(:,1:dimm_pc(1));
mr2vr_ = @(x)x(:);
tr2mr_ = @(x)reshape(x,[],size(x,3));
mrDist_clu = zeros(nClu);
norm_mean_ = @(x)sqrt(mean(x.^2));
norm_sum_ = @(x)sqrt(sum(x.^2));
normalize_mean_ = @(x)x ./ sqrt(mean(x.^2));
normalize_sum_ = @(x)x ./ sqrt(sum(x.^2));
for iClu1 = 1:nClu
    [trPc1, viSite1, trPos1] = deal(ctrPc_clu{iClu1}, cviSite_clu{iClu1}, ctrPos_clu{iClu1});
    vrDist_clu1 = zeros(nClu,1);
    for iiSite1 = 1:numel(viSite1)
        [mrPc1, mrPos1] = deal(trPc1(:,:,iiSite1), trPos1(:,:,iiSite1));
        vi_all2 = find(viSite_all == viSite1(iiSite1) & viClu_all ~= iClu1);
        if isempty(vi_all2), continue; end
        viClu2 = viClu_all(vi_all2);
        [trPc2, trPos2] = deal(trPc_all(:,:,vi_all2), trPos_all(:,:,vi_all2));
        switch mode_sim
            case 7 % centroid merging
                [mu1, sd1] = deal(mrPos1(:,1), mrPos1(:,2));
                mu2 = squeeze_(trPos2(:,1,:), 2);
                vrCorr12 = norm_mean_((mu1-mu2) ./ sd1);
                vrCorr12 = 1 - min(vrCorr12/10,1);
            case 6
                a = normalize_sum_(mr2vr_(pc2wav_(mrPv_global, mrPc1)));
                b = normalize_sum_(tr2mr_(pc2wav_(mrPv_global, trPc2)));
                c = a'*a;
                vrCorr12 = 1-abs(c-a'*b) ./ c;
            case 5 % rms after normalization
                a = normalize_sum_(mr2vr_(pc2wav_(mrPv_global, mrPc1)));
                b = normalize_sum_(tr2mr_(pc2wav_(mrPv_global, trPc2)));
                vrCorr12 = 1 - std(a-b, 1); % standardized rms
            case 4 % rms after normalization, .97 threshold
                a = normalize_sum_(mrPc1(:));
                b = normalize_sum_(reshape(trPc2, [], size(trPc2,3)));
                vrCorr12 = 1 - std(a-b, 1); % standardized rms
            case 3
                vrCorr12 = normalize_sum_(mrPc1(:))' * normalize_sum_(tr2mr_(trPc2));
            case 2
                vrCorr12 = 1-pdist2(mrPc1(:)', reshape(trPc2, [], size(trPc2,3))', 'cosine');
            case 1
                vrCorr12 = wavcor_(pc2wav_(mrPv_global, mrPc1), pc2wav_(mrPv_global, trPc2));
            otherwise, error('wave_similarity_clu_: invalid mode');
        end
        vrDist_clu1(viClu2) = max(vrDist_clu1(viClu2), vrCorr12(:));
    end
    mrDist_clu(:,iClu1) = vrDist_clu1;
end
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wave_similarity_clu1_(S0, P)
S_clu = S0.S_clu;

[KNN, MAX_SAMPLE, NUM_PC, fUse_raw] = deal(16, 4000, 3, 0);
nAve_knn = min(KNN, get_set_(P, 'knn', 30));
MIN_COUNT = get_set_(P, 'min_count');
trPc_spk = S0.trPc_spk(:,1:P.nSites_fet,:);

fprintf('Automated merging (post-hoc)\n'); t1=tic;
[viClu_spk, miKnn, vrRho_spk] = struct_get_(S_clu, 'viClu', 'miKnn', 'rho');
miKnn = miKnn(1:nAve_knn,:);
[viSite_spk, viTime_spk, mrPv] = struct_get_(S0, 'viSite_spk', 'viTime_spk', 'mrPv_global');
frac_thresh = get_set_(P, 'thresh_core_knn', .75);
nShift_max = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
viShift = -nShift_max:nShift_max;
dimm_spk = [size(trPc_spk,1), size(trPc_spk,2), size(S0.trPc_spk,3)];

% create template (nTemplate per cluster)
[ctrPc_clu, cviSite_clu] = deal(cell(S_clu.nClu, 1));
nSites = max(viSite_spk);
nDrift = get_set_(P, 'nTime_drift', 64);
nSpk_min = get_set_(P, 'knn', 30);
fprintf('\tComputing template\n\t'); t_template = tic;

% mrMean_site = zeros(nDrift, S_clu.nClu);
for iClu = 1:S_clu.nClu
    viSpk1 = find(S_clu.viClu == iClu);
    viSpk1 = viSpk1(:);
    viiSpk1 = round(linspace(1, numel(viSpk1), nDrift+1));
    [vlKeep_clu1, viSite_clu1] = deal(true(nDrift, 1), zeros(nDrift,1));
    trPc_drift1 = zeros(dimm_spk(1), dimm_spk(2), nDrift, 'single');
    [miKnn1, vrRho1] = deal(miKnn(:,viSpk1), vrRho_spk(viSpk1)');
    for iDrift = 1:nDrift
        vii1 = viiSpk1(iDrift):viiSpk1(iDrift+1);
        [viSpk11, vrRho11, miKnn11] = deal(viSpk1(vii1), vrRho1(vii1), miKnn1(:,vii1));
        
        switch 1 %3 % expand selection using miKnn
            case 1 % only connect to neighbors with higher density
                viSpk11 = miKnn11(vrRho_spk(miKnn11) >= vrRho11);
                viSpk11 = unique(viSpk11);
                iSite11 = mode(viSite_spk(viSpk11));
                vl_ = viSite_spk(viSpk11) == iSite11;
                viSpk11 = viSpk11(vl_);   
            case 2 % only connect to neighbors with higher density
                viSpk11 = unique(miKnn11(:));
                iSite11 = mode(viSite_spk(viSpk11));
                vl_ = viSite_spk(viSpk11) == iSite11;
                viSpk11 = viSpk11(vl_);                  
%                 mrMean_site(iDrift, iClu) = mean(vl_);
            case 4
                viSpk11 = viSpk11(vrRho11 > median(vrRho11));
                viSpk11 = miKnn(:, viSpk11);
                viSpk11 = viSpk11(:); % density based selection                
                iSite11 = mode(viSite_spk(viSpk11));
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11);   
            case 3     
                viSpk11 = miKnn(:, viSpk11);
                viSpk11 = viSpk11(:); % density based selection                
                iSite11 = mode(viSite_spk(viSpk11));
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11);   
            case 6               
                viSpk11 = miKnn(:, viSpk11);
                viSpk11 = viSpk11(:);
                viSpk11 = viSpk11(vrRho11 > quantile(vrRho11, frac_thresh));
                viSpk11 = viSpk11(viClu_spk(viSpk11) == iClu);
                iSite11 = mode(viSite_spk(viSpk11));
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11);
            case 5 % find local peak
                vrRho1 = S_clu.rho(viSpk1);
                [~,ii11] = max(vrRho1(vii1)); % local peak
                iSpk11 = viSpk1(vii1(ii11)); % center spike index
                iSite11 = viSite_spk(iSpk11);
                % find neighbors around the local peak
                viSpk11 = miKnn(:, miKnn(:,iSpk11));
                viSpk11 = viSpk11(:);
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11 & viClu_spk(viSpk11) == iClu);     
        end
        viSite_clu1(iDrift) = iSite11;
        if numel(viSpk11) < nSpk_min
            vlKeep_clu1(iDrift) = false;
            continue;
        end
        trPc_drift1(:,:,iDrift) = mean(trPc_spk(:,:,viSpk11),3);
    end
    ctrPc_clu{iClu} = trPc_drift1(:,:,vlKeep_clu1);
    cviSite_clu{iClu} = viSite_clu1(vlKeep_clu1);
    fprintf('.');
end
fprintf('\n\ttook %0.1fs\n', toc(t_template));


% merge the templates: todo, faster code
normalize_ = @(x)bsxfun(@rdivide, x, sqrt(sum(x.^2)));
fh_norm_tr = @(x)normalize_(reshape(pc2wav_(mrPv,x), [], size(x,3)));

mrDist_clu = nan(S_clu.nClu, 'single');
for iClu1 = 1:(S_clu.nClu-1)
    viSite_clu1 = cviSite_clu{iClu1};
    if isempty(viSite_clu1), continue; end
    mr1_ = fh_norm_tr(shift_trWav_(ctrPc_clu{iClu1}, viShift));
    viSite_clu1 = repmat(viSite_clu1(:), numel(viShift), 1);
    vrDist_clu1 = zeros(S_clu.nClu, 1, 'single');
    for iClu2 = (iClu1+1):S_clu.nClu
        viSite2 = cviSite_clu{iClu2};
        if ~any(ismember(viSite_clu1, viSite2)), continue; end
        viSite12 = intersect(viSite_clu1, viSite2);
%         if isempty(viSite12), continue; end
        mr2_ = fh_norm_tr(ctrPc_clu{iClu2});   
        for iSite12_ = 1:numel(viSite12)
            iSite12 = viSite12(iSite12_);
            mrDist12 = mr2_(:, viSite2==iSite12)' * mr1_(:, viSite_clu1==iSite12);
            vrDist_clu1(iClu2) = max(vrDist_clu1(iClu2), max(mrDist12(:)));
        end
    end
    mrDist_clu(:, iClu1) = vrDist_clu1;
end %for
end %func


%--------------------------------------------------------------------------
function [vr, vr_med] = sd_mad_(mr)
vr_med = median2_(mr);
vr = median2_(abs(mr - vr_med)) / .6745;
end % func


%--------------------------------------------------------------------------
function B = median2_(A,idimm)
if nargin<2, idimm=1; end

A = sort(A,idimm);
imid = ceil(size(A,idimm)/2);
if idimm==1
    B = A(imid,:);
else
    B = A(:,imid);
end
end %func


%--------------------------------------------------------------------------
function vrCorr12 = wavcor_pc_(mrPc1, trPc2, fMode, mrPc_sd1)
if numel(fMode)>1
    mrPv_global = fMode;
    fMode = 1;
end
if nargin<4, mrPc_sd1 = []; end

switch fMode 
    case 1 % correlation distance of waveform
        mr1 = mrPv_global * mrPc1;
        dimm1 = size(mr1);
        mr2 = reshape(mrPv_global * reshape(trPc2, size(trPc2,1), []), prod(dimm1), []);
        vr1 = mr1(:);
        vrCorr12 = sum((mr2-vr1).^2) ./ sum(vr1.^2);
        vrCorr12 = 1 - min(vrCorr12,1);
        
    case 100 % correlation distance of waveform
        mrWav1 = mrPv_global * mrPc1;                
        a_ = zscore_(mrWav1(:), 1);
        mrWav2 = mrPv_global * reshape(trPc2, size(trPc2,1), []);
        b_ = zscore_(reshape(mrWav2, [], size(trPc2,3)), 1);
        vrCorr12 = (a_'*b_) / size(mrWav1,1) / size(mrWav1,2);
        
    case 2 % RMS distance of PC
        vr1 = mrPc1(:);
        mrD = vr1 - reshape(trPc2,[],size(trPc2,3));
        vrCorr12 = sqrt(sum(mrD.^2) ./ sum(vr1.^2));
        vrCorr12 = 1 - min(vrCorr12/10,1);
        
    case 3.1 % RMS distance of PC
        vr1 = mrPc1(:);
        mrD = vr1 - reshape(trPc2,[],size(trPc2,3));
        vrCorr12 = 1 - min(sum(mrD.^2) ./ sum(vr1.^2), 1);
        
    case 3.3 % both RMS distance and normalized difference power
        vr1 = mrPc1(:);        
        mr2 = reshape(trPc2,[],size(trPc2,3));
        vr12 = sum((mr2-vr1).^2) ./ sum(vr1.^2);        
        vr1 = vr1 / sqrt(sum(vr1.^2));
        mr2 = mr2 ./ sqrt(sum(mr2.^2));
        vr22 = sum((mr2-vr1).^2);        
        vrCorr12 = 1 - min(min(vr12, vr22), 1);    
        
    case 3
        vr1 = mrPc1(:);        
        mr2 = reshape(trPc2,[],size(trPc2,3));        
        vr1 = vr1 / sqrt(sum(vr1.^2));
        mr2 = mr2 ./ sqrt(sum(mr2.^2));
        vrCorr12 = vr1' * mr2;
        
    case 4 % RMS distance of PC. same result as 2
        mrWav1 = mrPv_global * mrPc1;
        vr1 = mrWav1(:);
        mrWav2 = reshape(mrPv_global * reshape(trPc2, size(trPc2,1), []), [], size(trPc2,3));
        mrD = vr1 - mrWav2;
        vrCorr12 = sqrt(sum(mrD.^2) ./ sum(vr1.^2));
        vrCorr12 = 1 - min(vrCorr12,1);
        
    case 5
        mrWav1 = mrPv_global * mrPc1;
        a_ = zscore_(mrWav1(:), 1);
        mrWav2 = mrPv_global * reshape(trPc2, size(trPc2,1), []);
        b_ = zscore_(reshape(mrWav2, [], size(trPc2,3)), 1);
        vrCorr12 = (a_'*b_) / size(mrWav1,1) / size(mrWav1,2);    
        
    case 6 % normalized difference
        vr1 = mrPc1(:);        
        mr2 = reshape(trPc2,[],size(trPc2,3));        
        vr1 = vr1 / sqrt(sum(vr1.^2));
        mr2 = mr2 ./ sqrt(sum(mr2.^2));
        vrCorr12 = 1 - min(sum((vr1 - mr2).^2), 1);               
        
    case 7
        vr1 = mrPc1(:);
        mr1 = reshape(trPc2,[],size(trPc2,3));
        vrCorr1 = 1 - min(sum((vr1 - mr1).^2) ./ sum(vr1.^2), 1);
        vr2 = vr1/sqrt(sum(vr1.^2));
        mr2 = mr1./sqrt(sum(mr1.^2));
        vrCorr2 = vr2' * mr2;  
        vrCorr12 = (vrCorr1 + vrCorr2)/2;
        
    case 8
        vr1 = mrPc1(:);
        mrD = vr1 - reshape(trPc2,[],size(trPc2,3));
        vrCorr12 = 1 - min(sum(mrD.^2) ./ sum(mrPc_sd1(:).^2), 1);        
        
end %switch
end %func


%--------------------------------------------------------------------------
function vrCorr12 = wavcor_(mrWav1, trWav2, nShift)
% mrWav1 = meanSubt_(mrWav1);
% trWav2 = meanSubt_(trWav2);

if nargin<3, nShift = 2; end
% try time delay and match
switch 1
    case 1
        if nShift == 0
            a_ = zscore_(mrWav1(:));
            b_ = zscore_(reshape(trWav2, [], size(trWav2,3)),1);
            vrCorr12 = (a_'*b_) / size(mrWav1,1) / size(mrWav1,2);
        else
            vi0 = (nShift+1):(size(mrWav1,1)-nShift);
            b_ = zscore_(reshape(trWav2(vi0,:,:), [], size(trWav2,3)),1);
            mrCorr12 = zeros(nShift*2+1, size(trWav2,3));
            for iShift1 = 1:(nShift*2+1)
                a_ = mrWav1(vi0 + iShift1-1-nShift,:);
                mrCorr12(iShift1,:) = zscore_(a_(:))' * b_;
            end
            vrCorr12 = max(mrCorr12)' / numel(vi0) / size(mrWav1,2);
        end
    case 2
        nT = size(mrWav1,1);
        mrCorr12 = zeros(nShift*2+1, size(trWav2,3));
        [cvi1, cvi2] = shift_range_(nT, nShift);
        for iShift1 = 1:numel(cvi1)
            [vi1,vi2] = deal(cvi1{iShift1}, cvi2{iShift1});
            a_ = mrWav1(vi1,:);
            b_ = zscore_(reshape(trWav2(vi2,:,:), [], size(trWav2,3)),1);
            mrCorr12(iShift1,:) = zscore_(a_(:))' * b_ / numel(vi1);
        end
        vrCorr12 = max(mrCorr12)' / size(mrWav1,2);
end
end %func


%--------------------------------------------------------------------------
function out = pc2wav_(mrPv, mrPc)
switch ndims(mrPc)
    case 2, out = mrPv * mrPc;
    case 3        
        [nPc, nSites, nWav] = size(mrPc);
        out = reshape(mrPv * reshape(mrPc, nPc,[]), [], nSites, nWav);
end
end %func


%--------------------------------------------------------------------------
function [vr, vi] = cvr2vr_vi_(cvr)
vr = cell2mat_(cvr);
vn1 = cellfun(@(x)size(x,1), cvr);
vi = cell2mat_(arrayfun(@(x)repmat(x, vn1(x),1), 1:numel(cvr), 'UniformOutput', 0)');
end %func


%--------------------------------------------------------------------------
% todo: delayed execution, use parfeval
function S0 = sort1_(S0, P)

% drift processing
fprintf('Clustering\n'); 
runtime_sort = tic;

[nSites_fet, miSites, knn, nPcPerChan] = struct_get_(P, 'nSites_fet', 'miSites', 'knn', 'nPcPerChan');

% determine feature dimension
% if ~isempty(nFet_max)
%     nPcPerChan = round(nFet_max / nSites_fet);
%     nPcPerChan = min(max(nPcPerChan, 1), size(S0.trPc_spk,1));
% end
switch 1
    case 1
        get_pc_ = @(x)reshape(x(1:nPcPerChan,1:nSites_fet,:), nPcPerChan*nSites_fet, 1, size(x,3));
    case 2
        get_pc_ = @(x)reshape(x(1:nPcPerChan,:,:), nPcPerChan*size(x,2), 1, size(x,3));
    case 3
        get_pc_ = @(x)reshape(x, size(x,1)*size(x,2), 1, size(x,3));
end
if isempty(get_(S0, 'trPc2_spk'))
    trPc_spk = get_pc_(S0.trPc_spk);
else
    trPc_spk = cat(2, get_pc_(S0.trPc_spk), get_pc_(S0.trPc2_spk));
end

nSites = size(P.miSites,2);
cviSpk_site = arrayfun(@(x)find(S0.viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
cviSpk2_site = arrayfun(@(x)find(S0.viSite2_spk==x), 1:nSites, 'UniformOutput', 0)';
% vnSpk_site = cellfun(@numel, cviSpk_site);
nSpk = numel(S0.viSite_spk);

% parfor loop
if get_set_(P, 'fParfor', 1)
    gcp_ = gcp();
else
    gcp_ = [];
end
if 0
    [S0.mrPos_spk, S0.vrPow_spk] = calc_pos_spk_(S0.trPc_spk, S0.viSite_spk, P);
end
S_drift = calc_drift_(S0, P);
[viDrift_spk, mlDrift] = struct_get_(S_drift, 'viDrift_spk', 'mlDrift');

% Calculate Rho
P_sort = struct_set_(P, 'mlDrift', mlDrift, 'fSort_drift', 1);
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
miKnn = zeros(knn, nSpk, 'int32');
fprintf('Calculating Rho\n\t'); t1=tic;
[cvrRho, cmiKnn, cvrDelta, cviNneigh] = deal(cell(nSites,1));
% send jobs
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site2_(trPc_spk, cviSpk_site, cviSpk2_site, viDrift_spk, iSite, P_sort);
    if isempty(mrFet12), continue; end
    if ~isempty(gcp_)
        vS_out(iSite) = parfeval(gcp_, ...
            @(x,y,z,a)rho_knn_(x,y,z,a, P_sort), 2, mrFet12, viSpk12, viDrift12, n1);
    else
        [cvrRho{iSite}, cmiKnn{iSite}, P_sort.fGpu] = ...
            rho_knn_(mrFet12, viSpk12, viDrift12, n1, P_sort);      
    end
end %for
% collect jobs
if ~isempty(gcp_)
    for iSite1 = 1:nSites
        [iSite, vrRho1, miKnn1] = fetchNext(vS_out);
        [cvrRho{iSite}, cmiKnn{iSite}] = deal(vrRho1, miKnn1);
    end
end
% assemble jobs
for iSite = 1:nSites     
    viSpk1 = cviSpk_site{iSite};
    if isempty(viSpk1), continue; end
    [vrRho(viSpk1), miKnn(:,viSpk1)] = deal(cvrRho{iSite}, cmiKnn{iSite});
end %for
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), P_sort.fGpu, P_sort.fParfor);


% Calculate Delta
fprintf('Calculating Delta\n\t'); t2=tic;
% send jobs
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site2_(trPc_spk, cviSpk_site, cviSpk2_site, viDrift_spk, iSite, P_sort);
    if isempty(mrFet12), continue; end
    vrRho12 = vrRho(viSpk12);
    if ~isempty(gcp_)
        vS_out(iSite) = parfeval(gcp_, ...
            @(x,y,z,a,b)delta_knn_(x,y,z,a,b, P_sort), 2, mrFet12, viSpk12, vrRho12, viDrift12, n1);
    else
        [cvrDelta{iSite}, cviNneigh{iSite}, P_sort.fGpu] = ...
            delta_knn_(mrFet12, viSpk12, vrRho12, viDrift12, n1, P_sort);
    end
end %for
% collect jobs
if ~isempty(gcp_)
    for iSite1 = 1:nSites
        [iSite, vrDelta1, viNneigh1] = fetchNext(vS_out);
        [cvrDelta{iSite}, cviNneigh{iSite}] = deal(vrDelta1, viNneigh1);
    end
end
% assemble jobs
for iSite = 1:nSites     
    viSpk1 = cviSpk_site{iSite};
    if isempty(viSpk1), continue; end
    [vrDelta(viSpk1), viNneigh(viSpk1)] = deal(cvrDelta{iSite}, cviNneigh{iSite});
end %for
vrRho = vrRho / max(vrRho) / 10;     % divide by 10 to be compatible with previous version displays
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t2), P_sort.fGpu, P_sort.fParfor);

% output
runtime_sort = toc(runtime_sort);
[~, ordrho] = sort(vrRho, 'descend');
memory_sort = memory_matlab_();
S0.S_clu = struct('rho', vrRho, 'delta', vrDelta, 'ordrho', ordrho, 'nneigh', viNneigh, ...
    'P', P, 'miKnn', miKnn, 'S_drift', S_drift);
S0.runtime_sort = runtime_sort;
S0.memory_sort = memory_sort;
end %func


%--------------------------------------------------------------------------
function [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site2_(trPc_spk, cviSpk_site, cviSpk2_site, viDrift_spk, iSite, P)
% decide whether to use 1, 2, or 3 features

if isempty(cviSpk_site{iSite})
    [mrFet12, viSpk12, viDrift12] = deal([]);
    [n1, n2] = deal(0); 
    return;
end
% [nSites_fet, miSites] = struct_get_(P, 'nSites_fet', 'miSites');

[viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
if isempty(viSpk2)
    mrFet12 = squeeze_(trPc_spk(:,1,viSpk1),2);
else
    mrFet12 = [squeeze_(trPc_spk(:,1,viSpk1),2), squeeze_(trPc_spk(:,2,viSpk2),2)];
end
viSpk12 = [viSpk1; viSpk2];
[n1, n2] = deal(numel(viSpk1), numel(viSpk2));
        
if isempty(viDrift_spk)
    viDrift12 = [];
else
    viDrift12 = viDrift_spk(viSpk12);
end    

if get_set_(P, 'nFet_max', 0)
    if size(mrFet12,1) > nFet_max
        [~, mrFet12] = pca(mrFet12', 'NumComponents', nFet_max);
        mrFet12 = mrFet12';
    end
end
if get_set_(P, 'fWhiten_fet', 0)
	mrFet12 = whiten_(mrFet12);
end
end %func



%--------------------------------------------------------------------------
% todo: delayed execution, use parfeval
function S0 = sort2_(S0, P)

% drift processing
fprintf('Clustering\n'); 
runtime_sort = tic;

[nSites_fet, miSites, knn, nFet_max, nPcPerChan] = struct_get_(P, 'nSites_fet', 'miSites', 'knn', 'nFet_max', 'nPcPerChan');

% determine feature dimension
if ~isempty(nFet_max)
    nPcPerChan = floor(nFet_max / nSites_fet);
    nPcPerChan = min(max(nPcPerChan, 1), size(S0.trPc_spk,1));
end
trPc_spk = S0.trPc_spk(1:nPcPerChan,:,:);

nSites = size(P.miSites,2);
cviSpk_site = arrayfun(@(x)find(S0.viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
vnSpk_site = cellfun(@numel, cviSpk_site);
nSpk = numel(S0.viSite_spk);

% parfor loop
if get_set_(P, 'fParfor', 1)
    gcp_ = gcp();
else
    gcp_ = [];
end
if 0
    [S0.mrPos_spk, S0.vrPow_spk] = calc_pos_spk_(S0.trPc_spk, S0.viSite_spk, P);
end
S_drift = calc_drift_(S0, P);
[viDrift_spk, mlDrift] = struct_get_(S_drift, 'viDrift_spk', 'mlDrift');

% Calculate Rho
P_sort = struct_set_(P, 'mlDrift', mlDrift, 'fSort_drift', 1);
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
miKnn = zeros(knn, nSpk, 'int32');
fprintf('Calculating Rho\n\t'); t1=tic;
[cvrRho, cmiKnn, cvrDelta, cviNneigh] = deal(cell(nSites,1));
% send jobs
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, viDrift_spk, P_sort, iSite);
    if isempty(mrFet12), continue; end
    if ~isempty(gcp_)
        vS_out(iSite) = parfeval(gcp_, ...
            @(x,y,z,a)rho_knn_(x,y,z,a, P_sort), 2, mrFet12, viSpk12, viDrift12, n1);
    else
        [cvrRho{iSite}, cmiKnn{iSite}, P_sort.fGpu] = ...
            rho_knn_(mrFet12, viSpk12, viDrift12, n1, P_sort);      
    end
end %for
% collect jobs
if ~isempty(gcp_)
    for iSite1 = 1:nSites
        [iSite, vrRho1, miKnn1] = fetchNext(vS_out);
        [cvrRho{iSite}, cmiKnn{iSite}] = deal(vrRho1, miKnn1);
    end
end
% assemble jobs
for iSite = 1:nSites     
    viSpk1 = cviSpk_site{iSite};
    if isempty(viSpk1), continue; end
    [vrRho(viSpk1), miKnn(:,viSpk1)] = deal(cvrRho{iSite}, cmiKnn{iSite});
end %for
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), P_sort.fGpu, P_sort.fParfor);


% Calculate Delta
fprintf('Calculating Delta\n\t'); t2=tic;
% send jobs
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, viDrift_spk, P_sort, iSite);
    if isempty(mrFet12), continue; end
    vrRho12 = vrRho(viSpk12);
    if ~isempty(gcp_)
        vS_out(iSite) = parfeval(gcp_, ...
            @(x,y,z,a,b)delta_knn_(x,y,z,a,b, P_sort), 2, mrFet12, viSpk12, vrRho12, viDrift12, n1);
    else
        [cvrDelta{iSite}, cviNneigh{iSite}, P_sort.fGpu] = ...
            delta_knn_(mrFet12, viSpk12, vrRho12, viDrift12, n1, P_sort);
    end
end %for
% collect jobs
if ~isempty(gcp_)
    for iSite1 = 1:nSites
        [iSite, vrDelta1, viNneigh1] = fetchNext(vS_out);
        [cvrDelta{iSite}, cviNneigh{iSite}] = deal(vrDelta1, viNneigh1);
    end
end
% assemble jobs
for iSite = 1:nSites     
    viSpk1 = cviSpk_site{iSite};
    if isempty(viSpk1), continue; end
    [vrDelta(viSpk1), viNneigh(viSpk1)] = deal(cvrDelta{iSite}, cviNneigh{iSite});
end %for
vrRho = vrRho / max(vrRho) / 10;     % divide by 10 to be compatible with previous version displays
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t2), P_sort.fGpu, P_sort.fParfor);

% output
runtime_sort = toc(runtime_sort);
[~, ordrho] = sort(vrRho, 'descend');
memory_sort = memory_matlab_();
S0.S_clu = struct('rho', vrRho, 'delta', vrDelta, 'ordrho', ordrho, 'nneigh', viNneigh, ...
    'P', P, 'miKnn', miKnn, 'S_drift', S_drift);
S0.runtime_sort = runtime_sort;
S0.memory_sort = memory_sort;
end %func


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1, fGpu] = delta_knn_(mrFet12, viSpk12, vrRho12, viDrift12, n1, P)
if isempty(mrFet12)
    [vrDelta1, viNneigh1] = deal([]);
    return;
end
switch get_set_(P, 'fSort_drift', 1)
    case 1
        [vrDelta1, viNneigh1, fGpu] = delta_drift_knn_(mrFet12, viDrift12, P.mlDrift, vrRho12, n1, P);
        viNneigh1 = viSpk12(viNneigh1);
    case 0
        SINGLE_INF = 3.402E+38;
        vrRho12 = vrRho12(:);
        [vi12, vi1] = deal(1:size(mrFet12,2), 1:n1);

        % do cuda
        [vrDelta1, viNneigh1, fGpu] = cuda_delta_knn_(mrFet12, vrRho12, vi12, vi1, P);

        vrDelta1 = gather_(vrDelta1) .* vrRho12(vi1);
        viNneigh1 = viSpk12(gather_(viNneigh1));
        viNan = find(isnan(vrDelta1) | isinf(vrDelta1));
        viNneigh1(viNan) = viNan;
        vrDelta1(viNan) = sqrt(SINGLE_INF);
end %switch
end %func


%--------------------------------------------------------------------------
function [vrRho1, miKnn1, fGpu] = rho_knn_(mrFet12, viSpk12, viDrift12, n1, P)
fGpu = P.fGpu;
if isempty(mrFet12)
    [vrRho1, miKnn1] = deal([]);
    return;
end
switch get_set_(P, 'fSort_drift', 1)
    case 1
        [vrRho1, fGpu, miKnn1] = rho_drift_knn_(mrFet12, viDrift12, P.mlDrift, n1, P);  
        miKnn1 = viSpk12(miKnn1);
    case 0
        knn = get_set_(P, 'knn', 30);
        [vi12, vi1] = deal(1:size(mrFet12,2), 1:n1);
        [vrKnn1, fGpu, miKnn1] = cuda_knn_(mrFet12, vi12, vi1, P);
        vrRho1 = gather_(1./vrKnn1);

        n_ = size(miKnn1,1);
        if n_ == knn
            miKnn1(:,:) = gather_(miKnn1);
        else        
            miKnn1(1:n_,:) = gather_(miKnn1);
            miKnn1(n_+1:end,:) = repmat(gather_(miKnn1(end,:)), [knn-n_, 1]);
        end
        miKnn1 = viSpk12(miKnn1);    
end %switch
end %func


%--------------------------------------------------------------------------
function [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, viDrift_spk, P, iSite)

if isempty(cviSpk_site{iSite})
    [mrFet12, viSpk12, viDrift12] = deal([]);
    [n1, n2] = deal(0); 
    return;
end
[nSites_fet, miSites] = struct_get_(P, 'nSites_fet', 'miSites');

mrFet1 = trPc_spk(:, 1:nSites_fet, cviSpk_site{iSite});
mrFet1 = reshape(mrFet1, [], size(mrFet1,3));
viSites1 = miSites(1:nSites_fet,iSite);

% find neighboring sites
[~, mi_] = ismember(miSites, viSites1);
mi_(:,iSite)=0; %exclude seif
[~,mi2] = sort(mi_);
viSites2 = find(sum(mi_>0) == nSites_fet); %neighbor sites
miSites2 = mi2(end-nSites_fet+1:end, viSites2); % sites to index

% extract fet and build distance matrix
switch 1
    case 2
        mrFet12 = mrFet1;
    case 1
        cmrFet2 = cell(1, numel(viSites2));
        for iiSite2 = 1:numel(viSites2)
            trFet2_ = trPc_spk(:, miSites2(:,iiSite2), cviSpk_site{viSites2(iiSite2)});
            cmrFet2{iiSite2} = reshape(trFet2_, [], size(trFet2_,3));
        end
        mrFet12 = [mrFet1, cell2mat_(cmrFet2)];
end
n1 = numel(cviSpk_site{iSite});
n2 = size(mrFet12,2) - n1;
viSpk12 = cell2mat(cviSpk_site([iSite, viSites2]));
if isempty(viDrift_spk)
    viDrift12 = [];
else
    viDrift12 = viDrift_spk(viSpk12);
end

% spatial mask
if get_set_(P, 'fSpatialMask_clu', 1) && nSites_fet >= get_set_(P, 'nChans_min_car', 8)
    vrSpatialMask = spatialMask_(P, iSite, nSites_fet, P.maxDist_site_um);
    vrSpatialMask = repmat(vrSpatialMask(:)', [P.nPcPerChan, 1]);
    mrFet12 = bsxfun(@times, mrFet12, vrSpatialMask(:));
end
if 0
	mrFet12 = whiten_(mrFet12);
end
end %func


%--------------------------------------------------------------------------
% from file exchange
function [mr_out, whMat, invMat, vr_mean] = whiten_(mr_in, epsilon)
if nargin<2, epsilon=[]; end
if isempty(epsilon), epsilon = 0.0001; end
vr_mean = mean(mr_in,2);
X = mr_in - vr_mean;
A = X*X';
[V,D,~] = svd(A);
whMat = sqrt(size(mr_in,2)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
mr_out = whMat' * X;
if nargout>=3
    invMat = pinv(whMat);
end
end


%--------------------------------------------------------------------------
% uses parfor and temp file (memory efficient)
% write a temp file and delete
function S0 = detect_(P)
% keep all features in the memory, no disk storage

% parfor loop
if get_set_(P, 'fParfor', 1)
    gcp_ = gcp();
else
    gcp_ = [];
end

runtime_detect = tic; 
memory_init = memory_matlab_();

% load one
S_paged = readmda_paged_(P); % initialize
[nLoads, viOffset_load] = deal(S_paged.nLoads, S_paged.viOffset_load);
[mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part
cS_detect = cell(nLoads, 1);
cS_detect{1} = detect_paged_(mrWav_T1, P, makeStruct_(nlim_wav1)); % process the first part
mrWav_T1 = [];
[vrThresh_site, mrPv_global] = struct_get_(cS_detect{1}, 'vrThresh_site', 'mrPv_global');
S_cache = makeStruct_(vrThresh_site, mrPv_global);
fprintf('Memory use: %0.3f GiB\n', memory_matlab_()/2^30);
if ~isempty(gcp_)
    % write to temp files for each worker
    csFile_wav = arrayfun_(@(x)strrep(P.vcFile_prm, '.prm', sprintf('_wav%d.irc', x)), 1:nLoads-1);
    [cell_vnlim_wav, cell_dimm_wav] = deal(cell(size(csFile_wav)));   
    fprintf('detect_: writing %d files to disk\n', nLoads-1); t_write=tic;
    for iLoad = 1:nLoads-1
        [mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part
        write_bin_(csFile_wav{iLoad}, mrWav_T1);
        cell_dimm_wav{iLoad} = size(mrWav_T1);
        cell_vnlim_wav{iLoad} = nlim_wav1;
        mrWav_T1 = [];
        S_cache1 = S_cache;
        S_cache1.nlim_wav1 = cell_vnlim_wav{iLoad};
        S_cache1.dimm_wav1 = cell_dimm_wav{iLoad};
        S_cache1.vcFile_wav1 = csFile_wav{iLoad};
        vS_out(iLoad) = parfeval(gcp_, @(x)detect_paged_([],P,x), 1, S_cache1);
        fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
    end
    fprintf('\n\tWriting to disk took %0.1fs\n', toc(t_write));   
    fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
    for iLoad = 1:nLoads-1
        [completedIdx, S_] = fetchNext(vS_out);
        cS_detect{completedIdx+1} = S_;
    end
else
    for iLoad = 2:nLoads
        [mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part    
        S_cache.nlim_wav1 = nlim_wav1; % trim waveform
        cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache);        
        fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
        mrWav_T1 = [];
    end
end
S0 = detect_merge_(cS_detect, viOffset_load);
switch 1
    case 2
        trPc_spk1 = reshape(S0.mrVp_spk, 1, size(S0.mrVp_spk,1), size(S0.mrVp_spk,2));
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(S0.trPc_spk1, S0.viSite_spk, P);
    case 1
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(S0.trPc_spk, S0.viSite_spk, P);
end %switch
% Save output
runtime_detect = toc(runtime_detect);
memory_detect = memory_matlab_();
S0 = struct_add_(S0, vrPow_spk, vrThresh_site, mrPv_global, runtime_detect, P, memory_detect, memory_init, mrPos_spk);
fprintf('detect_: took %0.1fs (fParfor=%d, fGpu=%d)\n', runtime_detect, P.fParfor, P.fGpu);
end %func


%--------------------------------------------------------------------------
% uses parfor and temp file (memory efficient)
% write a temp file and delete
function S0 = detect2_(P)
% keep all features in the memory, no disk storage

% parfor loop
if get_set_(P, 'fParfor', 1)
    gcp_ = gcp();
else
    gcp_ = [];
end

runtime_detect = tic; 
memory_init = memory_matlab_();

% load one
S_paged = readmda_paged_(P); % initialize
[nLoads, viOffset_load] = deal(S_paged.nLoads, S_paged.viOffset_load);
[mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part
cS_detect = cell(nLoads, 1);
cS_detect{1} = detect_paged_(mrWav_T1, P, makeStruct_(nlim_wav1)); % process the first part
mrWav_T1 = [];
[vrThresh_site, mrPv_global] = struct_get_(cS_detect{1}, 'vrThresh_site', 'mrPv_global');
S_cache = makeStruct_(vrThresh_site, mrPv_global);
fprintf('Memory use: %0.3f GiB\n', memory_matlab_()/2^30);
if ~isempty(gcp_)
    % write to temp files for each worker
    csFile_wav = arrayfun_(@(x)strrep(P.vcFile_prm, '.prm', sprintf('_wav%d.irc', x)), 1:nLoads-1);
    [cell_vnlim_wav, cell_dimm_wav] = deal(cell(size(csFile_wav)));   
    fprintf('detect_: writing %d files to disk\n', nLoads-1); t_write=tic;
    for iLoad = 1:nLoads-1
        [mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part
        write_bin_(csFile_wav{iLoad}, mrWav_T1);
        cell_dimm_wav{iLoad} = size(mrWav_T1);
        cell_vnlim_wav{iLoad} = nlim_wav1;
        fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
        mrWav_T1 = [];
    end
    fprintf('\n\tWriting to disk took %0.1fs\n', toc(t_write));
    % process and collect work
    parfor iLoad = 1:nLoads-1
        S_cache1 = S_cache;
        S_cache1.nlim_wav1 = cell_vnlim_wav{iLoad};
        S_cache1.dimm_wav1 = cell_dimm_wav{iLoad};
        S_cache1.vcFile_wav1 = csFile_wav{iLoad};
        cS_detect{iLoad+1} = detect_paged_([], P, S_cache1); 
       % delete_(csFile_wav{iLoad}); % delete tmp file after done
    end    
else
    for iLoad = 2:nLoads
        [mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part    
        S_cache.nlim_wav1 = nlim_wav1; % trim waveform
        cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache);        
        fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
        mrWav_T1 = [];
    end
end
S0 = detect_merge_(cS_detect, viOffset_load);
switch 1
    case 2
        trPc_spk1 = reshape(S0.mrVp_spk, 1, size(S0.mrVp_spk,1), size(S0.mrVp_spk,2));
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(S0.trPc_spk1, S0.viSite_spk, P);
    case 1
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(S0.trPc_spk, S0.viSite_spk, P);
end %switch
% Save output
runtime_detect = toc(runtime_detect);
memory_detect = memory_matlab_();
S0 = struct_add_(S0, vrPow_spk, vrThresh_site, mrPv_global, runtime_detect, P, memory_detect, memory_init, mrPos_spk);
fprintf('detect_: took %0.1fs (fParfor=%d, fGpu=%d)\n', runtime_detect, P.fParfor, P.fGpu);
end %func


%--------------------------------------------------------------------------
% uses parfeval
function S0 = detect1_(P)
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
runtime_detect = tic; 
memory_init = memory_matlab_();

% load one
S_paged = readmda_paged_(P); % initialize
[nLoads, viOffset_load] = deal(S_paged.nLoads, S_paged.viOffset_load);
[mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part
cS_detect = cell(nLoads, 1);
cS_detect{1} = detect_paged_(mrWav_T1, P, makeStruct_(nlim_wav1)); % process the first part
mrWav_T1 = [];
[vrThresh_site, mrPv_global] = struct_get_(cS_detect{1}, 'vrThresh_site', 'mrPv_global');
S_cache1 = makeStruct_(vrThresh_site, mrPv_global);
fprintf('Memory use: %0.3f GiB\n', memory_matlab_()/2^30);
for iLoad = 2:nLoads          
    [mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); % process first part
    S_cache1.nlim_wav1 = nlim_wav1; % trim waveform
    if isempty(gcp_)            
        cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache1);            
    else
        vS_out(iLoad-1) = parfeval(gcp_, @(x,y)detect_paged_(x,P,y), 1, mrWav_T1, S_cache1);
    end    
    mrWav_T1 = [];
    fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
end
if ~isempty(gcp_)
    for iLoad = 2:nLoads
        [completedIdx, S_] = fetchNext(vS_out);
        cS_detect{completedIdx+1} = S_;
    end
end
S0 = detect_merge_(cS_detect, viOffset_load);
switch 1
    case 2
        trPc_spk1 = reshape(S0.mrVp_spk, 1, size(S0.mrVp_spk,1), size(S0.mrVp_spk,2));
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(S0.trPc_spk1, S0.viSite_spk, P);
    case 1
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(S0.trPc_spk, S0.viSite_spk, P);
end %switch
% Save output
runtime_detect = toc(runtime_detect);
memory_detect = memory_matlab_();
S0 = struct_add_(S0, vrPow_spk, vrThresh_site, mrPv_global, runtime_detect, P, memory_detect, memory_init, mrPos_spk);
fprintf('detect_: took %0.1fs (fParfor=%d, fGpu=%d)\n', runtime_detect, P.fParfor, P.fGpu);
end %func


%--------------------------------------------------------------------------
function [mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk, viSite_spk, P)
[~, nSites_spk, nSpk] = size(trPc_spk);
switch 1
    case 4, mrW = squeeze_(abs(trPc_spk(1,:,:)),1);
    case 3, mrW = squeeze_(sum(trPc_spk(1:P.nPcPerChan,:,:).^2),1);
    case 2, mrW = squeeze_(sum(trPc_spk.^2),1);
    case 1, mrW = squeeze_(trPc_spk(1,:,:),1) .^2;
end
vrPow_spk = sum(mrW,1)';

fprintf('Calculating spike positions\n\t');
t1=tic;
toVec_ = @(x)x(:);
switch 1        
    case 1 % centroid       
        trSiteXY_spk = reshape(P.mrSiteXY(P.miSites(:, viSite_spk),:), nSites_spk, nSpk, 2);
        mrPos_spk = [sum(mrW .* trSiteXY_spk(:,:,1))', sum(mrW .* trSiteXY_spk(:,:,2))'] ./ vrPow_spk;
%         figure; plot(mrPos_spk(:,1), mrPos_spk(:,2), '.', 'MarkerSize',1); axis([10 60 10 210])
        % figure; plot3(mrPos_spk(:,1), mrPos_spk(:,2), log(mean(mrW)'), '.');
        
    case 2 % 2d peak (non-uniform kernel)
        mrPos_spk = zeros(nSpk, 2, 'single');
        nSites = max(viSite_spk);
        cviSpk_site = arrayfun(@(x)find(viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
        for iSite = 1:nSites
            [viSpk1, viSite1] = deal(cviSpk_site{iSite}, P.miSites(:,iSite));
            [mrSiteXY1, mrW1] = deal(P.mrSiteXY(viSite1,:), mrW(:,viSpk1));
            mrPos_spk(viSpk1,:) = find_peak_2d_(mrSiteXY1, mrW1);
            fprintf('.');
        end
        
    case 3 % centroid by pc components
        [nPc, nSites_spk, nSpk] = size(trPc_spk);
        trPos_pc_spk = zeros(nSpk, 2, nPc, 'single');
        trSiteXY_spk = reshape(P.mrSiteXY(P.miSites(:, viSite_spk),:), nSites_spk, nSpk, 2);
        for iPc = 1:size(trPc_spk,1)
            mrW1 = squeeze_(trPc_spk(iPc,:,:),1).^2;
            mrPos1 = [sum(mrW1 .* trSiteXY_spk(:,:,1))', sum(mrW1 .* trSiteXY_spk(:,:,2))'] ./ sum(mrW1)';
            trPos_pc_spk(:,:,iPc) = mrPos1;
        end %for
        trPos_pc_spk = permute(trPos_pc_spk, [1,3,2]);        
end %switch
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function mrXY_spk1 = find_peak_2d_(mrXY1, mrW1)
% mrW1=double(mrW1);
% F = scatteredInterpolant(mrXY1(:,1), mrXY1(:,2), mrW1(:,1), 'natural');

% define grid
d0 = pdist(mrXY1); d0 = min(d0(d0>0));
dxy = 2; %4x interp
[xy_min, xy_max] = deal(min(mrXY1), max(mrXY1));
[xlim1, ylim1] = deal(...
        [xy_min(1)-d0, xy_max(1)+d0], ...
        [xy_min(2)-d0, xy_max(2)+d0]);
[xx,yy] = meshgrid(xlim1(1):dxy:xlim1(2), ylim1(1):dxy:ylim1(2));

% convolve and evaluate
mrXY2 = single([xx(:), yy(:)]);
A = exp(-(pdist2(mrXY1, mrXY2) / (d0)).^2); % 4 neighbors half
[~,viMax] = max(single(A)' * mrW1);
mrXY_spk1 = mrXY2(viMax,:);
mrXY_spk1 = mrXY_spk1 + (rand(size(mrXY_spk1))*2-1);

if 0
    mrW2 = reshape(mrW1(:,1)'*A, size(xx));
    figure; mesh(xx,yy,mrW2/10);
    hold on; plot3(mrXY1(:,1), mrXY1(:,2), mrW1(:,1), 'ro');
    hold on; plot3(mrXY_spk1(1,1), mrXY_spk1(1,2), mrW2(viMax(1)), 'r*');
    view(2)
end
end %func


%--------------------------------------------------------------------------
function S0 = detect_merge_(cS_detect, viOffset_load)

vnSpk_load = cellfun(@(x)numel(x.viSite_spk), cS_detect);
miSpk_load = [0; cumsum(vnSpk_load)];
miSpk_load = [miSpk_load(1:end-1)+1, miSpk_load(2:end)];

nSpk = sum(vnSpk_load);
[viSite_spk, viTime_spk, vrAmp_spk] = ...
    deal(zeros(nSpk, 1, 'int32'), zeros(nSpk, 1, 'int64'), zeros(nSpk, 1, 'single'));
viOffset_load = int64(viOffset_load);
[trPc_spk, mrVp_spk, trPc2_spk, viSite2_spk] = deal([]);

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
    if ~isempty(S1.trPc2_spk)
        if isempty(trPc2_spk)
            trPc2_spk = zeros(size(trPc_spk), 'single'); 
            viSite2_spk = zeros(size(viSite_spk), 'int32');
        end
        viSite2_spk(viSpk1) = S1.viSite2_spk;
        trPc2_spk(:,:,viSpk1) = S1.trPc2_spk;
    end
    if 0
        if isempty(mrVp_spk)
            mrVp_spk = zeros(size(S1.mrVp_spk,1), nSpk, 'single');
        end
        mrVp_spk(:,viSpk1) = S1.mrVp_spk;
    end
end

S0 = makeStruct_(viSite_spk, viTime_spk, vrAmp_spk, trPc_spk, mrVp_spk, viSite2_spk, trPc2_spk); 
end %func


%--------------------------------------------------------------------------
function S_detect = detect_paged_(mrWav_T, P, S_cache)
% S_detect = detect_paged_(mrWav_T, P, S_cache)
% S_detect = detect_paged_(file_name, P, S_cache)
 
if nargin<3, S_cache = []; end
if isempty(mrWav_T)
    mrWav_T = load_bin_(S_cache.vcFile_wav1, P.vcDataType, S_cache.dimm_wav1);
    delete_(S_cache.vcFile_wav1); % delete file 
end
% if ~P.fParfor && P.fGpu
%     mrWav_T = gpuArray_(mrWav_T);
% end
% filter and trim 
% nlim_wav1 = struct_get_(S_cache, 'nlim_wav1');
mrWav2 = filter_transpose_(mrWav_T, P);
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

% if get_set_(P, 'nSites_global', 8) >= 1
%     nSites = size(mrWav2,2);
%     if nSites <= get_set_(P, 'nSites_global', 8)
%         viSite_spk=ones(size(viSite_spk), 'like', viSite_spk); 
%     end
% end

% reject spikes within the overlap region
if ~isempty(nlim_wav1)
    viKeep_spk = find(viTime_spk >= nlim_wav1(1) & viTime_spk <= nlim_wav1(2));
    [viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)x(viKeep_spk), viTime_spk, vrAmp_spk, viSite_spk);    
end%if

% extract spike waveforms
trWav_spk = get_spkwav_(mrWav2, viSite_spk, viTime_spk, P);
if 0
    lim1 = round(size(trWav_spk,1) * [1/4,3/4]);
    mrVp_spk = min(trWav_spk(lim1(1):lim1(2),:,:));  
    mrVp_spk = gather_(squeeze_(mrVp_spk,1));
else
    mrVp_spk = [];
end
if 0
    dimm_spk = size(trWav_spk);
    trWav_spk = reshape(meanSubt_(reshape(trWav_spk, [], size(trWav_spk,3))), dimm_spk);
end

% extract spike feaures
if isempty(mrPv_global)
    [mrPv_global, vrD_global] = get_pv_(trWav_spk, P); 
end
trPc_spk = gather_(project_pc_(trWav_spk, mrPv_global, P));

if get_set_(P, 'sort_mode', 1) == 1
    viSite2_spk = find_site_spk23_(trWav_spk, viSite_spk, P);
    trWav_spk = []; %clear memory
    trWav2_spk = mn2tn_wav_spk2_(mrWav2, viSite2_spk, viTime_spk, P);
    trPc2_spk = gather_(project_pc_(trWav2_spk, mrPv_global, P));
else
    [viSite2_spk, trPc2_spk] = deal([]);
end

% return struct
if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end
S_detect = makeStruct_(trPc_spk, mrPv_global, viTime_spk, vrAmp_spk, viSite_spk, ...
    mrVp_spk, trPc2_spk, viSite2_spk, vrThresh_site);
end %func


%--------------------------------------------------------------------------
function trWav_spk1 = pc2spkwav_(trPc_spk, mrPv_global)
[nPc_spk, nSites_spk, nSpk] = size(trPc_spk);
nSamples_spk = size(mrPv_global,1);
dimm_spk = [nSamples_spk, nSites_spk, nSpk];
trWav_spk1 = reshape(mrPv_global * reshape(trPc_spk, size(trPc_spk,1), []), dimm_spk);
end %func


%--------------------------------------------------------------------------
function trPc_spk = project_pc1_(trWav_spk, mrPv)
[nSamples, nSites, nSpk] = size(trWav_spk);
nPc = size(mrPv,2);
mr = reshape(trWav_spk, size(trWav_spk,1), []);
trPc_spk = reshape((mrPv' * mr), nPc, nSites, nSpk);
end %func


%--------------------------------------------------------------------------
% 5/28/2019 JJJ: returns a cell of matrix
function trPc_spk = project_pc_(trWav_spk, mrPv, P)
% [mrFet1, mrFet2, mrFet3] = deal([]);
project_ = @(x,y)reshape(x' * reshape(y, size(y,1), []), size(x,2), size(y,2), size(y,3));
trPc_spk = project_(mrPv, trWav_spk); % project 0-shift

if get_set_(P, 'fInterp_fet', 0) == 0, return; end
    
% interpolate and shift
nPc = size(mrPv,2);
switch 1
    case 3
        viShift = [0, -1,-.5,.5,1,1.5,-1.5];     
    case 2
        viShift = [0, -1,-.5,.5,1,1.5,-1.5,2,-2]; 
    case 1
        viShift = [0, -1,-.5,.5,1]; 
end
trPv_shift = arrayfun(@(x)vr2mr_shift_(mrPv(:,x), viShift), 1:nPc, 'UniformOutput', 0);
trPv_shift = permute(cat(3, trPv_shift{:}),[1,3,2]);

% find best alignment using first channel
mrWav1 = squeeze_(trWav_spk(:,1,:), 2); % first channel
mrPv1_shift = squeeze_(trPv_shift(:,1,:),2); %first pv shift
[~, viMax_spk] = max(abs(mrPv1_shift' * mrWav1));  % project chan1 shift   

for iPc = 1:nPc
    for iShift=2:numel(viShift)
        viSpk1 = find(viMax_spk == iShift);
        if isempty(viSpk1), continue; end
        trPc_spk(:,:,viSpk1) = project_(trPv_shift(:,:,iShift), trWav_spk(:,:,viSpk1));
    end %for
end %for
end %func


%--------------------------------------------------------------------------
function mrPv1 = vr2mr_shift_(vr1, viShift)
vi0 = (1:numel(vr1))';
mrPv1 = zeros(numel(vr1), numel(viShift), 'like', vr1);
mrPv1(:,1) = vr1;
for iShift = 2:numel(viShift)
    mrPv1(:,iShift) = interp1(vi0, vr1, vi0+viShift(iShift), 'pchip', 'extrap');
end
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
% trFft = fft(tr);
% trFft = abs(trFft(2:end/2,:,:));
% figure; imagesc(std(trFft,1,3))
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
function [mrWav2, vrWav_mean2] = filter_transpose_(mnWav_T, P)
%-----
% Filter
fprintf('\tFiltering spikes...'); t_filter = tic;
if get_set_(P, 'fSmooth_spatial', 0)
    mnWav_T = spatial_smooth_(mnWav_T, P);
end
vcDataType_filter = get_set_(P, 'vcDataType_filter', 'single');
switch 2
    case 3
        mrWav2 = filt_car_(single(mnWav_T'), P);
        fCar=0;
    case 2
        mrWav2 = fft_filter_transpose(single(mnWav_T), P);
        fCar=1;
    case 1
        mrWav2 = fft_filter(single(mnWav_T'), P);
        fCar=1;
end
if fCar
    %global subtraction before 
    nChans = size(mrWav2, 2);
    if nChans >= get_set_(P, 'nChans_min_car', 8) %&& ndims(mnWav2) >= 3
        [mrWav2, vrWav_mean2] = wav_car_(mrWav2, P); 
    else
      vrWav_mean2 = [];
    end
    fprintf(' took %0.1fs\n', toc(t_filter));
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
%         out = out'; % transpose
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
function vcDir_in = get_test_data_(vcMode)
if nargin<1, vcMode = 'static'; end
switch vcMode
    case 'drift', vcDir_in = 'groundtruth/hybrid_synth/drift_siprobe/rec_64c_1200s_11'; 
    case 'static', vcDir_in = 'groundtruth/hybrid_synth/static_siprobe/rec_64c_1200s_11'; 
    case 'tetrode', vcDir_in = 'groundtruth/hybrid_synth/static_tetrode/rec_4c_1200s_11'; 
    case 'tetrode2', vcDir_in = 'groundtruth/hybrid_synth/static_tetrode/rec_4c_1200s_21'; 
    case 'tetrode3', vcDir_in = 'groundtruth/hybrid_synth/static_tetrode/rec_4c_1200s_31'; 
    case 'bionet', vcDir_in = 'groundtruth/bionet/bionet_static/static_8x_A_2A';
    case 'bionet1', vcDir_in = 'groundtruth/bionet/bionet_drift/drift_8x_A_2A';     
    case 'monotrode', vcDir_in = 'groundtruth/waveclus_synth/quiroga_difficult1/C_Difficult1_noise005';
    case 'monotrode1', vcDir_in = 'groundtruth/waveclus_synth/quiroga_difficult1/C_Difficult1_noise01';
    case 'monotrode2', vcDir_in = 'groundtruth/waveclus_synth/quiroga_difficult1/C_Difficult1_noise02';
    case 'monotrode3', vcDir_in = 'groundtruth/waveclus_synth/sim2_2K10/simulation_94';
end
if ispc()
    vcDir_in = strrep(vcDir_in, '/', '\');    
    vcDir_in = fullfile('C:\tmp', vcDir_in);
elseif isunix()
    vcDir_in = fullfile('~/ceph', vcDir_in);
end
end %func


%--------------------------------------------------------------------------
function P = makeParam_(vcDir_in, vcDir_out, vcFile_arg)
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end

if isempty(vcDir_out), vcDir_out = strrep(vcDir_in, 'groundtruth', 'irc2'); end
if ~exist_dir_(vcDir_out), mkdir(vcDir_out); end

% assume there is raw.mda, geom.csv, params.json, firings_true.mda
P = file2struct_(ircpath_(read_cfg_('default_prm', 0)));
P2 = file2struct_(ircpath_(read_cfg_('default2_prm', 0)));
P = struct_merge_(P, P2);

% now only supporting .mda file
P.vcFile = fullfile(vcDir_in, 'raw.mda');
P.vcDir_out = vcDir_out;
S_mda = readmda_header_(P.vcFile);
P.nChans = S_mda.dimm(1);
P.vcDataType = S_mda.vcDataType;
P.vnSamples_file = S_mda.dimm(2); % scalar
P.probe_file = fullfile(vcDir_in, 'geom.csv');

% probe file (mda format)
P.mrSiteXY = csvread(P.probe_file);
P.viSite2Chan = 1:size(P.mrSiteXY,1);
P.viShank_site = ones(size(P.viSite2Chan));
P.vrSiteHW = [12,12];

% load json file
S_json = loadjson_(fullfile(vcDir_in, 'params.json'));
P.sRateHz = get_set_(S_json, 'samplerate', P.sRateHz);
P.fInverse_file = get_set_(S_json, 'spike_sign', -1) == -1;

% read param
S_arg = [];
if isstruct(vcFile_arg)
    S_arg = vcFile_arg;
elseif exist_file_(vcFile_arg)
    if matchFileExt_(vcFile_arg, '.prm')
        S_arg = file2struct_(vcFile_arg);
    elseif matchFileExt_(vcFile_arg, '.txt')
        S_arg = import_spikeforest_args_(vcFile_arg);
    end
end
P = struct_merge_(P, S_arg);

% derived fields
P = fill_param_(P);
P.vcFile_prm = fullfile(vcDir_out, 'raw_geom.prm');

vcFile_gt = fullfile(vcDir_in, 'firings_true.mda');
if exist_file_(vcFile_gt)
    P.vcFile_gt = vcFile_gt;
end

edit_prm_file_(P, P.vcFile_prm);
fprintf('Created %s\n', P.vcFile_prm);
end %func


%--------------------------------------------------------------------------
function P = import_spikeforest_args_(vcArg_txt)
P = struct();

if exist_file_(vcArg_txt)    
    S_txt = meta2struct_(vcArg_txt);
else    
    return;
end
        
P.sRateHz = get_set_(S_txt, 'samplerate', 30000);
if isfield(S_txt, 'detect_sign')
    P.fInverse_file = ifeq_(S_txt.detect_sign>0, 1, 0);
elseif isfield(S_txt, 'spike_sign')
    P.fInverse_file = ifeq_(S_txt.spike_sign>0, 1, 0);
end

mask_out_artifacts = get_set_(S_txt, 'mask_out_artifacts', []);
if strcmpi(mask_out_artifacts, 'true')
    P.blank_thresh = 10; 
else
    P.blank_thresh = []; 
end    

% set adjacency radius
P.maxDist_site_um = get_set_(S_txt, 'adjacency_radius', 50);
P.maxDist_site_spk_um = get_set_(S_txt, 'adjacency_radius_out', 75);
if P.maxDist_site_spk_um<=0, P.maxDist_site_spk_um = inf; end
%P.maxDist_site_merge_um = P.maxDist_site_spk_um * 0.4667;   

% set frequency
freq_min = get_set_(S_txt, 'freq_min', []);
freq_max = get_set_(S_txt, 'freq_max', []);
if ~isempty(freq_min) && ~isempty(freq_max)
    P.freqLim = [freq_min, freq_max]; 
end

% vcCommonRef and whiten
fWhiten = strcmpi(get_(S_txt, 'whiten'), 'True');
vcCommonRef = get_set_(S_txt, 'common_ref_type', 'none');
if fWhiten
    if isempty(vcCommonRef) || strcmpi(vcCommonRef, 'none')
        P.vcCommonRef = 'whiten';
    else
        P.vcCommonRef = ['whiten-', vcCommonRef];
    end
else
    P.vcCommonRef = vcCommonRef;
end

% copy fields (don't copy empty fields)
P = struct_copyas_(P, S_txt, ...
    {'detect_threshold', 'pc_per_chan', 'merge_thresh', 'scale_factor'}, ...
    {'qqFactor', 'nPcPerChan', 'maxWavCor', 'uV_per_bit'}, 1);

% String parameters
P = struct_copyas_(P, S_txt, {'filter_type', 'feature_type'}, {'vcFilter', 'vcFet'});

% same name
P = struct_copyas_(P, S_txt, ...
    {'knn', 'batch_sec_drift', 'step_sec_drift', 'min_count', 'nSites_whiten', ...
    'fft_thresh', 'delta_cut', 'fft_thresh_low', 'post_merge_mode', 'sort_mode'});

% set boolean
P = set_bool_(P, 'fGpu', S_txt);
P = set_bool_(P, 'fSave_spkwav', S_txt);

% override settings
P.fParfor = 0; %disable parfor when running spikeforest platform
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
        fprintf('\tnTime_clu = %d (batch_sec_drift = %0.1f s)\n', P.nTime_clu, batch_sec_drift);
    catch
        P.nTime_clu = 1;
    end
end
if isempty(get_(P, 'nTime_drift'))
    try
        step_sec_drift = get_set_(P, 'step_sec_drift', 10);
        P.nTime_drift = max(round(recording_duration_(P) / step_sec_drift), 1);
        fprintf('\tnTime_drift = %d (step_sec_drift = %0.1f s)\n', P.nTime_drift, step_sec_drift);
    catch
        P.nTime_drift = 64;
    end
end

P.bytesPerSample = bytesPerSample_(P.vcDataType);
P = struct_default_(P, 'vcFile_prm', subsFileExt_(P.vcFile, '.prm'));
% if ~isempty(get_(P, 'gain_boost')), P.uV_per_bit = P.uV_per_bit / P.gain_boost; end
P.spkThresh = P.spkThresh_uV / P.uV_per_bit;
P = struct_default_(P, 'cvrDepth_drift', {});
% P = struct_default_(P, {'maxSite_fet', 'maxSite_detect', 'maxSite_sort','maxSite_pix', 'maxSite_dip', 'maxSite_merge', 'maxSite_show'}, P.maxSite);
P = struct_default_(P, 'mrColor_proj', [.75 .75 .75; 0 0 0; 1 0 0]);
P.mrColor_proj = reshape(P.mrColor_proj(:), [], 3); %backward compatible
P = struct_default_(P, {'blank_thresh', 'thresh_corr_bad_site', 'tlim_load'}, []);
if numel(P.tlim)==1, P.tlim = [0, P.tlim]; end
if isfield(P, 'rejectSpk_mean_thresh'), P.blank_thresh = P.rejectSpk_mean_thresh; end
P.vcFilter = get_filter_(P);
if isempty(get_(P, 'vcFilter_show'))
    P.vcFilter_show = P.vcFilter;
end
P.nC_max = read_cfg_('nC_max'); % override nC_max (gpu parameter)
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
% 8/2/17 JJJ: added '.' if dir is empty
% 7/31/17 JJJ: Substitute file extension
function varargout = subsFileExt_(vcFile, varargin)
% Substitute the extension part of the file
% [out1, out2, ..] = subsFileExt_(filename, ext1, ext2, ...)
if isempty(vcFile), varargout{1}=[]; return; end
[vcDir_, vcFile_, ~] = fileparts(vcFile);
if isempty(vcDir_), vcDir_ = '.'; end
for i=1:numel(varargin)
    vcExt_ = varargin{i};    
    varargout{i} = [vcDir_, filesep(), vcFile_, vcExt_];
end
end %func


%--------------------------------------------------------------------------
function [vi_uniq, vn_uniq, cvi_uniq] = unique_count_(vi)
% count number of unique elements and sort by vn_uniq
vi_uniq = unique(vi);
cvi_uniq = arrayfun(@(x)find(vi==x), vi_uniq, 'UniformOutput', 0);
vn_uniq = cellfun(@numel, cvi_uniq);
[vn_uniq, ix] = sort(vn_uniq, 'descend');
vi_uniq = vi_uniq(ix);
cvi_uniq = cvi_uniq(ix);
end %func


%--------------------------------------------------------------------------
% viSite1_uniq1, vn_uniq1, cviSpk1_uniq1
function [viUnique, vnUnique, cvi_uniq] = unique_count__(vi)
% vn = [1 1 1 1 4 4 4 8];
% vn = [1 1 1 1];
[viUnique, vnUnique, cvi_uniq] = deal([], [], {});
if isempty(vi), return; end
vi1 = sort(vi(:)');
vi_ = find(diff(sort(vi1)) > 0);
if isempty(vi_)
    vnUnique = 1;
    viUnique = vi1(1);
else
    vnUnique = diff([0, vi_, numel(vi)]);
    viUnique = [vi1(1), vi1(vi_+1)];
end
if nargout >= 3
    cvi_uniq = arrayfun(@(x)find(vi==x), viUnique, 'UniformOutput', 0);
end
end %func


%--------------------------------------------------------------------------
function z = zscore_(x, flag, dim)
if isempty(x), z=[]; return; end
if nargin < 2, flag = 0; end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
mu = mean(x,dim);
sigma = std(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
end %func

%--------------------------------------------------------------------------
function S = struct_add_(varargin)
% S = struct_add_(S, var1, var2, ...)
% output
% S.var1=var1; S.var2=var2; ...

S = varargin{1};
for i=2:numel(varargin)
    try
        S.(inputname(i)) = varargin{i};
    catch
        disperr_();
    end
end
end %func


%--------------------------------------------------------------------------
function frewind_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function disperr_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function struct_save_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function edit_prm_file_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function write_bin_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function delete_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function plot_FigWavCor_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function close_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function auto_scale_proj_time_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function save_log_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end

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
function out1 = cell2mat_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_default_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_filter_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = set0_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = postCluster_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gt2mda_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_dir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = meanSubt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
% function out1 = zscore_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = templateMatch_post_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = memory_matlab_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = calc_drift_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = recording_duration_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = squeeze_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refresh_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = find_site_spk23_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mn2tn_wav_spk2_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = matchFileExt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_bin_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = arrayfun_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = delete_clu_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = shift_trWav_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_copyas_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = set_bool_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ifeq_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = spatial_smooth_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2hash_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = dir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = msgbox_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = figures_manual_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = plot_FigRD_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = plot_FigWav_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_sort_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refrac_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = correlogram_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = keyPressFcn_cell_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = button_CluWav_simulate_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_fig_cache_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_position_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

function [out1, out2] = readmda_header_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = mr2thresh_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = gpuArray_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = filt_car_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = findNearSites_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = spatialMask_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = shift_range_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = wav_car_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = get_fig_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end

function [out1, out2, out3] = plan_load_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = detect_spikes_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = cuda_delta_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = cuda_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = unique_count_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = rho_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = delta_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end


