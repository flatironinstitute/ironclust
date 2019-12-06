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
fDetect = read_cfg_('fForceRerun');

if iscell(vcDir_in) && iscell(vcDir_out)    
    [csDir_in, csDir_out, csFile_arg] = deal(vcDir_in, vcDir_out, vcFile_arg);
    parfor iFile = 1:numel(csDir_in)
        try
            if exist_file_(fullfile(csDir_out{iFile}, 'firings.mda')) && ~fDetect
                continue;
            end
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
[fDetect, fSort] = deal(exist_file_(vcDir_in) || exist_dir_(vcDir_in)); % cmd mode
[P, S0, fPlot_gt, fValidate] = deal([]); 
switch lower(vcCmd)
    case 'compile-deploy', compile_cuda_([], '1'); return
    case {'compile', 'install'}, compile_cuda_([], '0'); return
    case 'readmda_header', varargout{1} = readmda_header_(vcArg1); return;
    case 'mcc', irc('mcc'); return; 
    case {'join-mda', 'join_mda', 'joinmda'}
        join_mda_(vcArg1, vcArg2); return;
    case {'readmda', 'read-mda', 'read_mda'}
        [varargout{1}, varargout{2}] = readmda_(vcArg1); return;
    case 'import-clip'
        [S0, P] = import_clip_(vcArg1); 
    case 'edit', edit_(vcArg1); return;
    case 'juxta'
        convert_mda_ui('english'); return;
    case 'version'
        if nargout==0, version_(); 
        else, varargout{1} = version_(); 
        end
        return;
    case 'scoreboard', irc2_scoreboard(); return;
    case {'spikesort', 'detect-sort', 'sort', 'auto', '', 'describe', 'manual', ...
            'verify', 'auto-verify', 'sort-verify', 'spikesort-verify'}
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
        switch lower(vcCmd)
            case {'detect-sort', 'spikesort', 'spikesort-verify'}, clear_(); fDetect = 1; fSort = 1;
            case {'sort', 'sort-verify'}, clear_('sort'); fDetect = 0; fSort = 1;          
            case 'describe', describe_(vcFile_prm); return;
            case {'verify', 'validate'}, validate_(P, fPlot_gt); return;
            case 'manual', manual_(P); return;
        end
        fValidate = contains(vcCmd, {'auto-verify', 'sort-verify', 'spikesort-verify'});
    case 'benchmark'
        if nargout==0, benchmark_(vcArg1, vcArg2, vcArg3); 
        else, varargout{1} = benchmark_(vcArg1, vcArg2, vcArg3); 
        end
        return;
    case 'plot', irc('plot', vcArg1, vcArg2); return;
    case 'clear', clear_(); vcFile_prm_=[]; return;
    case 'clear-sort', clear_('sort'); return;     
    case {'test-mcc', 'test_mcc', 'testmcc'}, test_mcc_(vcArg1); return;
    case {'test-static', 'test-drift', 'test-tetrode', 'test-tetrode2', 'test-tetrode3', ...
            'test-bionet', 'test-bionet1', 'test-monotrode', ...
            'test-monotrode1', 'test-monotrode2', 'test-monotrode3'}
        vcDir_in = get_test_data_(strsplit_get_(vcCmd,'-',2));
        [fDetect, fSort, fValidate] = deal(1, 1, 1);
    case 'export', irc('export', vcArg1); return;
    case {'export-phy', 'phy'}, irc2phy(vcArg1, vcArg2); return;
    case {'export-klusters', 'klusters', 'neurosuite'}, irc2klusters_v2(vcArg1, vcArg2); return;
    otherwise, clear_();
end

fprintf('Running irc2.m (%s)\n', version_());

if isempty(P)
    P = makeParam_(vcDir_in, vcDir_out, vcFile_arg);
end
vcFile_prm_ = P.vcFile_prm;
vcFile_mat = strrep(P.vcFile_prm, '.prm', '_irc.mat');
if exist_file_(vcFile_mat) && ~fDetect
    S0 = load0_(P.vcFile_prm);
else
    clear_(P.vcFile_prm);
    S0 = detect_(P); 
    save0_(S0);
end
% set(0, 'UserData', S0);

if isempty(struct_get_(S0, 'S_clu')) || fSort
    S0.S_clu = sort_(S0, P);
%     set(0, 'UserData', S0);
end
% set(0, 'UserData', S0);
S0 = auto_(S0, P);
describe_(S0);
% set(0, 'UserData', S0);

% save clu
save_clu_(S0.S_clu, P);

% output
vcFile_firings_mda = fullfile(P.vcDir_out, 'firings.mda');
save_firings_mda_(S0, vcFile_firings_mda);

% Validate
if fValidate, validate_(P, fPlot_gt); end
end %func


%--------------------------------------------------------------------------
function save_clu_(S_clu, P)
if isempty(S_clu), return; end

vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');

miKnn = get_(S_clu, 'miKnn');
if ~isempty(miKnn)
    S_clu.miKnn = [];
    S_clu.dimm_knn = size(miKnn);
    S_clu.type_knn = class(miKnn);
    write_bin_([vcFile_prm_, '_knn.irc'], miKnn);
end

struct_save_(S_clu, [vcFile_prm_, '_clu_irc.mat'], 1);
end %func


%--------------------------------------------------------------------------
function [S0, P] = import_clip_(vcFile_mat)
% import monotrode clips

S_mat = load(vcFile_mat);
P = file2struct_(ircpath_(read_cfg_('default_prm', 0)));
P2 = file2struct_(ircpath_(read_cfg_('default2_prm', 0)));
P = struct_merge_(P, P2);
P_ = S_mat.par;
P.sRateHz = double(P_.sr);
P.qqFactor = double(P_.stdmin);
P.spkLim = [-double(P_.w_pre), double(P_.w_post)];
P.miSites = 1;
P.viSite2Chan = 1;
P.nSites_fet = 1;
P.fParfor = 0;
P.mrSiteXY = [0,0];
P.nTime_clu = 1;
P.nTime_drift = 1;
P.nC_max = read_cfg_('nC_max'); % override nC_max (gpu parameter)
P.viShank_site = 1;
P.vcDir_out = fileparts(vcFile_mat);
P.fPlot_gt = 0;
P.vcFile_prm = fullfile(P.vcDir_out, 'raw_geom.prm');
edit_prm_file_(P, P.vcFile_prm);

nSpikes = numel(S_mat.index);
S0 = struct('viTime_spk', int64(S_mat.index), ...
    'viSite_spk', ones(nSpikes,1,'int32'), 'P', P);
S0.mrPos_spk = zeros(nSpikes, 2, 'single');
S0.vrAmp_spk = max(abs(S_mat.spikes),[],2);
[mrPv, mrPc_spk] = pca(S_mat.spikes, ...
    'NumComponents', P.nPc_spk, 'Centered', 'off');
S0.trPc_spk = reshape(single(mrPc_spk)', [P.nPc_spk, 1, nSpikes]);
S0.mrPv_global = single(mrPv);
S0.runtime_detect = nan;
S0.memory_init = memory_matlab_();
S0.memory_detect = nan;
end %func


%--------------------------------------------------------------------------
function trPc_spk = load_fet_(S0, P, iFet)
if nargin<3, iFet = 1; end

% return empty if trPc_spk should be loaded partially
if get_set_(P, 'fLargeRecording', 0)
    trPc_spk=[]; 
    return; 
end

vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');
fprintf('Loading feature %d...', iFet); t1=tic;
if iFet==1
    if ~isempty(get_(S0, 'trPc_spk'))
        trPc_spk = S0.trPc_spk;
        return;
    end
    vcFile_fet = [vcFile_prm_, '_fet.irc'];
    if exist_file_(vcFile_fet)
        trPc_spk = load_bin_(vcFile_fet, S0.type_fet, S0.dimm_fet);
    else
        csFiles_in = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet_%d.irc',x)], 1:S0.nLoads);
        trPc_spk = load_bin_merge_(csFiles_in, S0.type_fet, S0.dimm_fet);
    end
elseif iFet==2
    if isempty(get_(S0, 'viSite2_spk'))
        trPc_spk = [];
        return;
    end
    if ~isempty(get_(S0, 'trPc2_spk'))
        trPc_spk = S0.trPc2_spk;
        return;
    end    
    vcFile_fet = [vcFile_prm_, '_fet2.irc'];
    if exist_file_(vcFile_fet)
        trPc_spk = load_bin_(vcFile_fet, S0.type_fet, S0.dimm_fet);
    else
        csFiles_in = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet2_%d.irc',x)], 1:S0.nLoads);
        trPc_spk = load_bin_merge_(csFiles_in, S0.type_fet, S0.dimm_fet);
    end
else
    error('invalid feature id');
end
fprintf(' took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function trPc_spk = load_bin_merge_(csFiles_in, type_fet, dimm_fet)
try
    trPc_spk = zeros(dimm_fet, type_fet);
catch
    trPc_spk = []; return; % memory out
end
iOffset = 0;
for iFile = 1:numel(csFiles_in)
    nBytes_file1 = filesize_(csFiles_in{iFile});
    dimm_fet1 = dimm_fet;
    if numel(dimm_fet) == 1
        dimm_fet1 = nBytes_file1 / bytesPerSample_(type_fet);        
    else
        dimm_fet1(end) = nBytes_file1 / bytesPerSample_(type_fet) / prod(dimm_fet(1:end-1));
    end
    tr_ = load_bin_(csFiles_in{iFile}, type_fet, dimm_fet1);
    vi1 = (1:dimm_fet1(end)) + iOffset;
    trPc_spk(:,:,vi1) = tr_;
    iOffset = vi1(end);
end %for
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

% keep trPc_spk loaded for manual
S0.trPc_spk = load_fet_(S0, P, 1);
if isempty(get_(S0.S_clu, 'trWav_spk_clu'))
    S0.S_clu = calc_clu_wav_(S0, P);    
    S0.S_clu = S_clu_quality_(S0.S_clu, P);
end
set(0, 'UserData', S0);

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
if ~isempty(vcFile_mat1) && ~read_cfg_('fForceRerun')
    % load already processed data
    vcFile_prm1 = strrep(vcFile_mat1{1}, '_irc.mat', '.prm');
    vcConsoleOut = evalc(sprintf('irc2(''describe'', ''%s'');', vcFile_prm1));
    fprintf('Loaded benchmark from %s\n', vcFile_mat1{1});
else
    % process the data
    fprintf('Running benchmark: ''%s'' using ''%s'': ', vcDir_in1, vcParam1); t1=tic;
    vcCmd = ifeq_(ispc(), 'run_irc', './run_irc');
    [~, vcConsoleOut] = system(sprintf('%s %s %s %s', vcCmd, vcDir_in1, vcDir_out1, vcParam1));
    fprintf('took %0.1fs\n', toc(t1));    
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
vcVer = 'v5.2.1';
vcDate = '12/6/2019';
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

vcFile_prm_ = strrep(S0.P.vcFile_prm, '.prm', '');

trPc_spk = gather_(get_(S0, 'trPc_spk'));
S0.trPc_spk = [];
fSave_fet = get_set_(S0.P, 'fSave_fet', 1);
if ~isempty(trPc_spk) && fSave_fet
    S0.trPc_spk = [];
    S0.dimm_fet = size(trPc_spk);
    S0.type_fet = class(trPc_spk);
    write_bin_([vcFile_prm_, '_fet.irc'], trPc_spk);
end

trPc2_spk = gather_(get_(S0, 'trPc2_spk'));
S0.trPc2_spk = [];
if ~isempty(trPc2_spk) && fSave_fet
    S0.trPc2_spk = [];
    write_bin_([vcFile_prm_, '_fet2.irc'], trPc2_spk);
end

S_clu = get_(S0, 'S_clu');
S0.S_clu = [];
struct_save_(S0, [vcFile_prm_, '_irc.mat'], 1);

if ~isempty(S_clu), save_clu_(S_clu, P); end
end %


%--------------------------------------------------------------------------
function S0 = load0_(vcFile_prm)
fprintf('Loading %s... ', vcFile_prm); t1=tic;
% set(0, 'UserData', []);
vcFile_mat = strrep(vcFile_prm, '.prm', '_irc.mat');
vcFile_clu_mat = strrep(vcFile_prm, '.prm', '_clu_irc.mat');
vcFile_knn = strrep(vcFile_prm, '.prm', '_knn.irc');

S0 = load(vcFile_mat);
if isempty(get_(S0, 'S_clu')) && exist_file_(vcFile_clu_mat)
    S0.S_clu = load(vcFile_clu_mat);
end

if isempty(get_(S0.S_clu, 'miKnn'))
    if exist_file_(vcFile_knn)
        S0.S_clu.miKnn = load_bin_(vcFile_knn, S0.S_clu.type_knn, S0.S_clu.dimm_knn);
    end
end
% set(0, 'UserData', S0);
fprintf('took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
% negative index means from the end, 0 index means end index
function vc1 = strsplit_get_(vc,delim,idx)
cs = strsplit(vc, delim);
idx = mod(idx-1,numel(cs))+1;
vc1 = cs{idx};
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
    otherwise
        if exist_file_(vcMode)
            vcFile_prm_ = strrep(vcMode, '.prm', '');
            delete([vcFile_prm_, '*.irc']);
            delete([vcFile_prm_, '*_irc.mat']);
        else
            error('pass a parameter file');
        end
end
end

%--------------------------------------------------------------------------
function csDesc = describe_(S0)
% usage
% ----
% describe_(S0)
% describe_(vcFile_prm)
% csDesc = describe_(S0)

if ischar(S0)
    vcFile_prm = S0;
    S0 = load0_(vcFile_prm);
end 
P=S0.P;

try
    [runtime_detect, runtime_sort, runtime_automerge] = ...
        deal(S0.runtime_detect, S0.S_clu.runtime_sort, S0.S_clu.runtime_automerge);
    runtime_total = runtime_detect + runtime_sort + runtime_automerge;
catch
    runtime_total = nan;
end
try
    tDur = recording_duration_(S0.P, S0); 
catch
    tDur = nan;
end
memory_sort = S0.S_clu.memory_sort - S0.memory_init;
memory_detect = S0.memory_detect - S0.memory_init;
nSites = numel(P.viSite2Chan);
nShanks = max(P.viShank_site);
nSites_spk = size(P.miSites,1);
nSpk = numel(S0.viTime_spk);
try
    nFeatures = S0.S_clu.nFeatures;
    nPcPerChan = nFeatures / nSites_spk;
catch
    nFeatures = P.nSites_fet * P.nPcPerChan;
    nPcPerChan = P.nPcPerChan;
end

csDesc = {};
try
    csDesc = {};
    csDesc{end+1} = sprintf('');    
    csDesc{end+1} = sprintf('------------------------------');    
    csDesc{end+1} = sprintf('Summary of %s', P.vcFile_prm);
    csDesc{end+1} = sprintf('------------------------------');    
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
    csDesc{end+1} = sprintf('    Filter range (Hz):      [%0.1f, %0.1f]', P.freqLim);
    csDesc{end+1} = sprintf('    Common ref:             %s', P.vcCommonRef);
    csDesc{end+1} = sprintf('    FFT threshold:          %d', get_set_(P, 'fft_thresh', 0));
    csDesc{end+1} = sprintf('Events');
    csDesc{end+1} = sprintf('    #Spikes:                %d', nSpk);
    csDesc{end+1} = sprintf('    Feature extracted:      %s', P.vcFet);    
    csDesc{end+1} = sprintf('    #Sites/event:           %d', nSites_spk);
    csDesc{end+1} = sprintf('    maxDist_site_um:        %0.0f', P.maxDist_site_um);    
    csDesc{end+1} = sprintf('    maxDist_site_spk_um:    %0.0f', P.maxDist_site_spk_um);
    csDesc{end+1} = sprintf('    spkLim_ms:              [%0.3f, %0.3f]', P.spkLim_ms);
    csDesc{end+1} = sprintf('    #Features/event:        %d', nFeatures);    
    csDesc{end+1} = sprintf('    #PC/chan:               %d', nPcPerChan);
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
    csDesc{end+1} = sprintf('Auto-merge');   
    csDesc{end+1} = sprintf('    delta_cut:              %0.3f', get_set_(P, 'delta_cut', 1));
    csDesc{end+1} = sprintf('    maxWavCor:              %0.3f', P.maxWavCor);
end
try
    csDesc{end+1} = sprintf('Runtime (s)');
    csDesc{end+1} = sprintf('    Detect + feature (s):   %0.1fs', runtime_detect);    
    csDesc{end+1} = sprintf('    Cluster runtime (s):    %0.1fs', runtime_sort);
    csDesc{end+1} = sprintf('    merge runtime (s):      %0.1fs', runtime_automerge);
    csDesc{end+1} = sprintf('    Total runtime (s):      %0.1fs', runtime_total);
    csDesc{end+1} = sprintf('    Runtime speed:          x%0.1f realtime', tDur / runtime_total);    

    csDesc{end+1} = sprintf('memory usage (GiB):         %0.3f', max(memory_detect, memory_sort)/2^30);
    csDesc{end+1} = sprintf('    detect(GiB):            %0.3f', memory_detect/2^30);
    csDesc{end+1} = sprintf('    sort(GiB):              %0.3f', memory_sort/2^30);

    csDesc{end+1} = sprintf('Execution');
    csDesc{end+1} = sprintf('    fGpu (GPU use):         %d', P.fGpu);
    csDesc{end+1} = sprintf('    fParfor (parfor use):   %d', P.fParfor);
    csDesc{end+1} = sprintf('    fSave_fet:              %d', get_set_(P, 'fSave_fet', 1));
    csDesc{end+1} = sprintf('    fLargeRecording:        %d', get_set_(P, 'fLargeRecording', 0));    
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

fprintf('\nauto-merging...\n'); runtime_automerge = tic;


% refresh clu, start with fundamentals
S0.S_clu = struct_copy_(S0.S_clu, ...
    'rho', 'delta', 'ordrho', 'nneigh', 'P', 'miKnn', 'S_drift', ...
    'nFeatures', 'runtime_sort', 'memory_sort');
S0.S_clu = postCluster_(S0.S_clu, P, S0.viSite_spk); % peak merging

maxWavCor = get_set_(P, 'maxWavCor', .99);
viClu_delete = [];
if maxWavCor<1
    post_merge_mode = get_set_(P, 'post_merge_mode', 1);    
    switch post_merge_mode 
        case 1, [mrDist_clu, viClu_delete] = wave_similarity_clu_(S0, P);
        case 2, mrDist_clu = calc_dist_ccg(S0, P);
    end
    S0.S_clu = templateMatch_post_(S0.S_clu, P, mrDist_clu);
end
% delete low snr clu
if ~isempty(viClu_delete)
    nClu_pre = S0.S_clu.nClu;
    S0.S_clu = delete_clu_(S0.S_clu, viClu_delete);
    fprintf('calc_clu_wav_: %d->%d clusters, %d removed below SNR=%0.1f\n', ...
        nClu_pre, S0.S_clu.nClu, numel(viClu_delete), P.min_snr_clu);
end
% compute SNR per cluster and remove small SNR
S0.S_clu = S_clu_refresh_(S0.S_clu, 1, S0.viSite_spk);
S0.S_clu = S_clu_sort_(S0.S_clu, 'viSite_clu');
S0.S_clu = S_clu_refrac_(S0.S_clu, P, [], S0.viTime_spk); % refractory violation removal

S0.S_clu.runtime_automerge = toc(runtime_automerge);
fprintf('\n\tauto-merging took %0.1fs (fGpu=%d, fParfor=%d)\n', ...
    S0.S_clu.runtime_automerge, P.fGpu, P.fParfor);
end %func


%--------------------------------------------------------------------------
function [S_clu, viClu_delete] = calc_clu_wav_(S0, P)

S_clu = S0.S_clu;
S_clu.nClu = max(S_clu.viClu);
cviSpk_clu = arrayfun_(@(x)find(S_clu.viClu==x), (1:S_clu.nClu)');
viSite_clu = cellfun(@(x)mode(S0.viSite_spk(x)), cviSpk_clu);
nT = size(S0.mrPv_global,1);
if isempty(get_(S0, 'trPc_spk'))
    trPc_spk = load_fet_(S0, P, 1);
else
    trPc_spk = S0.trPc_spk;
end
nSites_spk = size(trPc_spk,2);
nSites = numel(P.viSite2Chan);
trWav_clu = zeros(nT, nSites_spk, S_clu.nClu, 'single');
tmrWav_clu = zeros(nT, nSites, S_clu.nClu, 'single'); %global waveform
for iClu = 1:S_clu.nClu
    viSpk1 = cviSpk_clu{iClu};
    if isempty(viSpk1), continue; end
    viSite1 = S0.viSite_spk(viSpk1);
    iSite1 = viSite_clu(iClu);
    trPc1 = trPc_spk(:,:,viSpk1);
    mrPc1 = mean(trPc1(:,:,viSite1==iSite1),3);
    mrWav1 = S0.mrPv_global * mrPc1;
    trWav_clu(:,:,iClu) = mrWav1;
    viSite_clu1 = P.miSites(:,iSite1);
    tmrWav_clu(:,viSite_clu1,iClu) = mrWav1;
end
[trWav_raw_clu, tmrWav_raw_clu] = deal(trWav_clu, tmrWav_clu); % for now raw is not saved
[trWav_spk_clu, tmrWav_spk_clu] = deal(trWav_clu, tmrWav_clu); % for now raw is not saved
S_clu = struct_add_(S_clu, trWav_clu, trWav_raw_clu, tmrWav_clu, viSite_clu, tmrWav_raw_clu, trWav_spk_clu, tmrWav_spk_clu);

% compute SNR
mrWav1_clu = squeeze_(trWav_clu(:,1,:),2); 
vrVmax_clu = max(mrWav1_clu)';
vrVmin_clu = min(mrWav1_clu)';
vrVpp_clu = vrVmax_clu - vrVmin_clu;
vrRms_site = S0.vrThresh_site(:) / S0.P.qqFactor;
S_clu.vrSnr_clu = abs(vrVmin_clu(:)) ./ vrRms_site(viSite_clu);
S_clu.vrSnr2_clu = abs(vrVpp_clu(:)) ./ vrRms_site(viSite_clu);

viClu_delete = [];

% update similarity
S0.S_clu = S_clu;
mrWavCor = wave_similarity_clu_(S0, P);
S_clu.mrWavCor = mrWavCor + mrWavCor'; % make it symmetric

S_clu = S_clu_position_(S_clu);
S_clu = S_clu_refresh_(S_clu, 1, S0.viSite_spk);
if ~isfield(S_clu, 'csNote_clu'), S_clu.csNote_clu = cell(S_clu.nClu, 1); end
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
% keeps trPc in the main memory, not sent out to workers
function [mrDist_clu, viClu_remove] = wave_similarity_clu_(S0, P)
S_clu = S0.S_clu;

fprintf('Automated merging (post-hoc)\n'); t1=tic;
[viClu_spk, miKnn, vrRho_spk] = struct_get_(S_clu, 'viClu', 'miKnn', 'rho');
if false
    KNN_MAX = 16;
    nAve_knn = min(KNN_MAX, get_set_(P, 'knn', 30));
    miKnn = miKnn(1:nAve_knn,:);
end
[viSite_spk, viSite2_spk, mrPv] = ...
    struct_get_(S0, 'viSite_spk', 'viSite2_spk', 'mrPv_global');
nShift_max = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
viShift = -nShift_max:nShift_max;
dimm_spk = [S0.dimm_fet(1), P.nSites_fet, S0.dimm_fet(3)];

% create template (nTemplate per cluster)
cS_clu = cell(S_clu.nClu, 1);
nDrift = get_set_(P, 'nTime_drift', 64);
nSpk_min = get_set_(P, 'knn', 30);
nClu = S_clu.nClu;

fprintf('\tComputing template\n'); t_template = tic;
cviSpk_clu = arrayfun_(@(x)find(S_clu.viClu==x), (1:S_clu.nClu)');
S_pre = makeStruct_(vrRho_spk, viSite_spk, miKnn, nDrift, dimm_spk, ...
    nSpk_min, viSite2_spk, P);
fParfor = P.fParfor;
if fParfor
    try
        parfor iClu = 1:S_clu.nClu
            cS_clu{iClu} = wave_similarity_clu_pre_(cviSpk_clu{iClu}, S_pre);
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iClu = 1:S_clu.nClu
        cS_clu{iClu} = wave_similarity_clu_pre_(cviSpk_clu{iClu}, S_pre);
    end
end
[ctrPc_clu, cviSite_clu, viClu_remove] = merge_clu_pre_(cS_clu, S0, P);
fprintf('\n\ttook %0.1fs\n', toc(t_template));

% merge the templates: todo, faster code
mrDist_clu = nan(S_clu.nClu, 'single');
S_post = makeStruct_(cviSite_clu, ctrPc_clu, nClu, viShift, mrPv);
if fParfor
    try
        parfor iClu1 = 1:(nClu-1)
            mrDist_clu(:, iClu1) = wav_similarity_clu_post_(iClu1, S_post);
        end %for
    catch
        fParfor = 0;
    end 
end
if ~fParfor
    for iClu1 = 1:(nClu-1)
        mrDist_clu(:, iClu1) = wav_similarity_clu_post_(iClu1, S_post);
    end %for
end
end %func


%--------------------------------------------------------------------------
% multiply triangular weighing factor
function [ctrPc_clu, cviSite_clu, viClu_remove] = merge_clu_pre_(cS_pre3, S0, P)
[ctrPc_clu, cviSite_clu] = deal(cell(size(cS_pre3)));

% P.fLargeRecording = 0; % debug

trPc_spk = load_fet_(S0, P, 1);
if isempty(trPc_spk)
    [ctrPc_clu, cviSite_clu] = merge_clu_pre_load_(cS_pre3, S0, P);
else
    switch get_set_(P, 'mode_merge_mean', 2)
        case 1, mean_ = @(y,x)mean(y(:,:,x),3);
        case 2
            nSites_fet = min(S0.dimm_fet(2), P.nPc_spk);
            mean_ = @(y,x)mean(y(:,1:nSites_fet,x),3);
        case 3, mean_ = @(y,x)median(y(:,:,x),3);
    end

    for iClu = 1:numel(cS_pre3)
        S_pre3 = cS_pre3{iClu};
        cviSite_clu{iClu} = [S_pre3.viSite1(:); S_pre3.viSite2(:)];
        trPc_clu1 = cellfun_(@(x)mean_(trPc_spk,x), S_pre3.cviSpk1);
        ctrPc_clu{iClu} = cat(3,trPc_clu1{:}); % faster than pre-allocating
    end
    trPc_spk = []; % clear from memory

    % Second peak average
    trPc2_spk = load_fet_(S0, P, 2);
    if isempty(trPc2_spk), return; end
    for iClu = 1:numel(cS_pre3)
        S_pre3 = cS_pre3{iClu};
        trPc_clu2 = cellfun_(@(x)mean_(trPc2_spk,x), S_pre3.cviSpk2);
        trPc_clu2 = cat(3,trPc_clu2{:});
        ctrPc_clu{iClu} = cat(3, ctrPc_clu{iClu}, trPc_clu2);
    end
end
% remove low-snr peaks
viClu_remove = find_low_snr_clu_(S0, P, ctrPc_clu, cviSite_clu);
end %func


%--------------------------------------------------------------------------
function viClu_remove = find_low_snr_clu_(S0, P, ctrPc_clu, cviSite_clu)
[nRemoved, nTotal] = deal(0);
for iClu = 1:numel(ctrPc_clu)
    [trPc_clu1, viSite_clu1] = deal(ctrPc_clu{iClu}, cviSite_clu{iClu});
    if isempty(trPc_clu1), continue; end
    vrVmin_clu1 = min(S0.mrPv_global * squeeze_(trPc_clu1(:,1,:),2),[],1);
    vrThresh1 = S0.vrThresh_site(viSite_clu1);
    vrSnr_clu1 = abs(vrVmin_clu1(:)) ./ vrThresh1(:) * P.qqFactor;
    vlKeep1 = vrSnr_clu1 >= P.qqFactor;    
    [nRemoved, nTotal] = deal(nRemoved + sum(~vlKeep1), nTotal + numel(vlKeep1));
    [ctrPc_clu{iClu}, cviSite_clu{iClu}] = deal(trPc_clu1(:,:,vlKeep1), viSite_clu1(vlKeep1));
end %for
viClu_remove = find(cellfun(@isempty, ctrPc_clu));
fprintf('merge: removed %d/%d templates below SNR=%0.3f\n', nRemoved, nTotal, P.min_snr_clu);    
end %func


%--------------------------------------------------------------------------
% long recording version
function [ctrPc_clu, cviSite_clu] = merge_clu_pre_load_(cS_pre3, S0, P)
[ctrPc_clu, cviSite_clu] = deal(cell(size(cS_pre3)));

% P.fLargeRecording=0; % debug
[cviFet_clu, ccviSpk_clu] = deal(cell(size(cS_pre3)));
for iClu = 1:numel(cS_pre3)
    S_pre3 = cS_pre3{iClu};
    cviSite_clu{iClu} = [S_pre3.viSite1(:); S_pre3.viSite2(:)];    
    n1 = numel(S_pre3.viSite1);
    n2 = numel(S_pre3.viSite2);
    cviFet_clu{iClu} = [repmat(1,n1,1); repmat(2,n2,1)];
    if n2>0
        ccviSpk_clu{iClu} = [cS_pre3{iClu}.cviSpk1, cS_pre3{iClu}.cviSpk2];
    else
        ccviSpk_clu{iClu} = cS_pre3{iClu}.cviSpk1;
    end
end
ctrPc_clu = mean_fet_site_load_(cviSite_clu, cviFet_clu, ccviSpk_clu, S0, P);
% viClu_remove = find_low_snr_clu_(S0, P, ctrPc_clu, cviSite_clu);
end %func


%--------------------------------------------------------------------------
function ctrPc_clu = mean_fet_site_load_(cviSite_clu, cviFet_clu, ccviSpk_clu, S0, P)
fprintf('mean_fet_load_... '); t1=tic;
viSite_all = cell2mat_(cviSite_clu);
viFet_all = cell2mat_(cviFet_clu);
cviSpk_all = [ccviSpk_clu{:}];
ctrPc_all = cell(size(cviSpk_all));
switch get_set_(P, 'mode_merge_mean', 2)
    case 1, mean_ = @(y,x)mean(y(:,:,x),3);
    case 2
        nSites_fet = min(S0.dimm_fet(2), P.nPc_spk);
        mean_ = @(y,x)mean(y(:,1:nSites_fet,x),3);
    case 3, mean_ = @(y,x)median(y(:,:,x),3);
end

% read the whole site and trim
[cviSpk_site, cviSpk2_site] = calc_cviSpk_site_(S0, P);
[type_fet, dimm_fet] = deal(S0.type_fet, S0.dimm_fet);
nSites = size(P.miSites, 2);
for iSite = 1:nSites
    vi_all1 = find(viSite_all == iSite & viFet_all == 1);
    [trPc_site1, trPc_site2] = load_fet_drift_(S0, iSite);
    
    [~,cviSpk_all1] = cellfun_(@(x)ismember(x,cviSpk_site{iSite}), cviSpk_all(vi_all1));
    ctrPc_all(vi_all1) = cellfun_(@(x)mean_(trPc_site1,x), cviSpk_all1);
    if ~isempty(trPc_site2)
        vi_all2 = find(viSite_all == iSite & viFet_all == 2);
        [~,cviSpk_all2] = cellfun_(@(x)ismember(x,cviSpk2_site{iSite}), cviSpk_all(vi_all2));
        ctrPc_all(vi_all2) = cellfun_(@(x)mean_(trPc_site2,x), cviSpk_all2);
    end
end

% repackage to clu index
vnSpk_clu = cellfun(@numel, cviSite_clu);
ctrPc_clu = arrayfun_(@(x,y)cat(3, ctrPc_all{x:y}), cumsum(vnSpk_clu)-vnSpk_clu+1, cumsum(vnSpk_clu));
fprintf('took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function mr = mask_mr_(mr, ml)
mr(~ml) = 0;
end %func


%--------------------------------------------------------------------------
function vrDist_clu1 = wav_similarity_clu_post_(iClu1, S3)
[cviSite_clu, ctrPc_clu, nClu, viShift, mrPv] = ...
    struct_get_(S3, 'cviSite_clu', 'ctrPc_clu', 'nClu', 'viShift', 'mrPv');
vrDist_clu1 = zeros(nClu, 1, 'single');
viSite_clu1 = cviSite_clu{iClu1};
if isempty(viSite_clu1), return; end
mr1 = normalize_tr_(shift_trWav_(ctrPc_clu{iClu1}, viShift), mrPv);
viSite_clu1 = repmat(viSite_clu1(:), numel(viShift), 1);
for iClu2 = (iClu1+1):nClu
    viSite2 = cviSite_clu{iClu2};
    if ~any(ismember(viSite_clu1, viSite2)), continue; end
    viSite12 = intersect(viSite_clu1, viSite2);
    mr2 = normalize_tr_(ctrPc_clu{iClu2}, mrPv);   
    for iSite12_ = 1:numel(viSite12)
        iSite12 = viSite12(iSite12_);
        mrDist12 = mr2(:, viSite2==iSite12)' * mr1(:, viSite_clu1==iSite12);
        vrDist_clu1(iClu2) = max(vrDist_clu1(iClu2), max(mrDist12(:)));
    end
end
end %func


%--------------------------------------------------------------------------
function mrWav1 = normalize_tr_(mrPc, mrPv)
mrWav = reshape(pc2wav_(mrPv,mrPc), [], size(mrPc,3));
mrWav1 = bsxfun(@rdivide, mrWav, sqrt(sum(mrWav.^2)));
end %func


%--------------------------------------------------------------------------
function S_pre3 = wave_similarity_clu_pre_(viSpk1, S)
switch 2
    case 1, S_pre3 = wave_similarity_clu_pre1_(viSpk1, S);
    case 2, S_pre3 = wave_similarity_clu_pre2_(viSpk1, S);
end
end %func


%--------------------------------------------------------------------------
% iterate by drift, not sites
function S_pre3 = wave_similarity_clu_pre2_(viSpk1, S)

[vrRho_spk, viSite_spk, viSite2_spk, miKnn, nDrift, nSpk_min, P] = ...
    struct_get_(S, 'vrRho_spk', 'viSite_spk', 'viSite2_spk', 'miKnn', 'nDrift', 'nSpk_min', 'P');
viSpk1 = viSpk1(:);
viiSpk1 = round(linspace(1, numel(viSpk1), nDrift+1));
vrRho1 = vrRho_spk(viSpk1)';
viSite1 = viSite_spk(viSpk1)';
[cviSite1, cviSite2, cviSpk1, cviSpk2] = deal({});
for iDrift = 1:nDrift
    vii1 = viiSpk1(iDrift):viiSpk1(iDrift+1);
    viSite11 = viSite1(vii1);
    iSite11 = mode(viSite11);
    vii1 = vii1(viSite11 == iSite11);
    if ~isempty(miKnn)
        miKnn11 = miKnn(:, viSpk1(vii1)); % todo: load from disk
    else
        miKnn11 = load_miKnn_site_(P, iSite11, viSite_spk, viSpk1(vii1));
    end
    viSpk11 = unique(miKnn11(vrRho_spk(miKnn11) >= vrRho1(vii1)));
    
    viSite_spk111 = viSite_spk(viSpk11);
    iSite111 = mode(viSite_spk111);
    viSpk111 = viSpk11(viSite_spk111 == iSite111);   
    if numel(viSpk111) >= nSpk_min
        cviSpk1{end+1} = viSpk111;
        cviSite1{end+1} = iSite111;
    end
    
    if ~isempty(viSite2_spk)
        viSite_spk112 = viSite2_spk(viSpk11);
        iSite112 = mode(viSite_spk112);
        viSpk112 = viSpk11(viSite_spk112 == iSite112);   
        if numel(viSpk112) >= nSpk_min
            cviSpk2{end+1} = viSpk112;
            cviSite2{end+1} = iSite112;
        end    
    end
end
[viSite1, viSite2] = deal(cell2mat(cviSite1), cell2mat(cviSite2)); 
S_pre3 = makeStruct_(viSite1, viSite2, cviSpk1, cviSpk2);
end %func


%--------------------------------------------------------------------------
% iterate by drift, not sites
function S_pre3 = wave_similarity_clu_pre1_(viSpk1, S)

[vrRho_spk, viSite_spk, viSite2_spk, miKnn, nDrift, dimm_spk, nSpk_min] = ...
    struct_get_(S, 'vrRho_spk', 'viSite_spk', 'viSite2_spk', 'miKnn', 'nDrift', 'dimm_spk', 'nSpk_min');
viSpk1 = viSpk1(:);
viiSpk1 = round(linspace(1, numel(viSpk1), nDrift+1));
[miKnn1, vrRho1] = deal(miKnn(:,viSpk1), vrRho_spk(viSpk1)');
[viSite1, viSite2, cviSpk1, cviSpk2] = deal({});
for iDrift = 1:nDrift
    vii1 = viiSpk1(iDrift):viiSpk1(iDrift+1);
    miKnn11 = miKnn1(:,vii1);
    viSpk11 = unique(miKnn11(vrRho_spk(miKnn11) >= vrRho1(vii1)));
    
    viSite_spk111 = viSite_spk(viSpk11);
    iSite111 = mode(viSite_spk111);
    viSpk111 = viSpk11(viSite_spk111 == iSite111);   
    if numel(viSpk111) >= nSpk_min
        cviSpk1{end+1} = viSpk111;
        viSite1{end+1} = iSite111;
    end
    
    if ~isempty(viSite2_spk)
        viSite_spk112 = viSite2_spk(viSpk11);
        iSite112 = mode(viSite_spk112);
        viSpk112 = viSpk11(viSite_spk112 == iSite112);   
        if numel(viSpk112) >= nSpk_min
            cviSpk2{end+1} = viSpk112;
            viSite2{end+1} = iSite112;
        end    
    end
end
[viSite1, viSite2] = deal(cell2mat(viSite1), cell2mat(viSite2)); 
S_pre3 = makeStruct_(viSite1, viSite2, cviSpk1, cviSpk2);
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
function [cviSpk_site, cviSpk2_site] = calc_cviSpk_site_(S0, P)
nSites = size(P.miSites,2);
cviSpk_site = vi2cell_(S0.viSite_spk, nSites);
% cviSpk_site = arrayfun(@(x)find(S0.viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
viSite2_spk = get_(S0, 'viSite2_spk');
if ~isempty(viSite2_spk)
    cviSpk2_site = vi2cell_(S0.viSite2_spk, nSites);
%     cviSpk2_site = arrayfun(@(x)find(S0.viSite2_spk==x), 1:nSites, 'UniformOutput', 0)';
else
    cviSpk2_site = cell(size(cviSpk_site));
end
end %func


%--------------------------------------------------------------------------
% todo: delayed execution, use parfeval
function S_clu = sort_(S0, P)

% drift processing
fprintf('Clustering\n'); 
runtime_sort = tic;
[trPc_spk, nFeatures] = get_pc_sort_(S0, P);
S_drift = calc_drift_(S0, P);

if isempty(trPc_spk)
    [vrRho, vrDelta, miKnn, viNneigh, memory_sort] = sort_disk_(S0, P, S_drift);
else
    [vrRho, vrDelta, miKnn, viNneigh, memory_sort] = sort_ram_(S0, P, trPc_spk, S_drift);
end

% output
vrRho = vrRho / max(vrRho) / 10;     % divide by 10 to be compatible with previous version displays
[~, ordrho] = sort(vrRho, 'descend');
S_clu = struct('rho', vrRho, 'delta', vrDelta, 'ordrho', ordrho, 'nneigh', viNneigh, ...
    'P', P, 'miKnn', miKnn, 'S_drift', S_drift, 'nFeatures', nFeatures);
S_clu.runtime_sort = toc(runtime_sort);
S_clu.memory_sort = memory_sort;
end %func


%--------------------------------------------------------------------------
function [vrRho, vrDelta, miKnn, viNneigh, memory_sort] = sort_disk_(S0, P, S_drift)

[cviSpk_site, cviSpk2_site] = calc_cviSpk_site_(S0, P);
% viDrift_spk = get_viDrift_spk_(S_drift);
[viLim_drift, mlDrift] = get_(S_drift, 'viLim_drift', 'mlDrift');
[viSite_spk] = struct_get_(S0, 'viSite_spk');
nSpk = numel(viSite_spk);
nSites = size(P.miSites,2);
fParfor = get_set_(P, 'fParfor', true) && nSites>1;
nPcPerChan = get_set_(P, 'nPcPerChan', 0);
miKnn = [];

% Calculate Rho
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
viNneigh = zeros(nSpk, 1, 'int64');
mlPc = calc_mlPc_(nPcPerChan, S0.dimm_fet);
S_fet = struct_add_(P, mlDrift, fParfor, mlPc, viLim_drift);
S_fet = struct_merge_(S_fet, ...
    struct_copy_(S0, 'type_fet', 'dimm_fet', 'nLoads', 'ccviSpk_site_load', 'ccviSpk_site2_load', 'vnSpk_load'));

fprintf('Calculating Rho (sort_disk_)\n\t'); t1=tic;
for iSite = 1:nSites         
    [viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
    if isempty(viSpk1), continue; end
    vrRho(viSpk1) = rho_knn_disk_(iSite, viSpk1, viSpk2, S_fet);    
end
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), S_fet.fGpu, fParfor);

fprintf('Calculating Delta (sort_disk_)\n\t'); t1=tic;
for iSite = 1:nSites 
    [viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
    if isempty(viSpk1), continue; end
    [vrDelta(viSpk1), viNneigh(viSpk1)] = delta_knn_disk_(iSite, viSpk1, viSpk2, vrRho, S_fet);
end
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), S_fet.fGpu, fParfor);
memory_sort = memory_matlab_();
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function n = bytesPerSample_(vcDataType)
% Return number of bytes per data type
% syntax
% -----
% n = bytesPerSample_(vcDataType)   % pass string
% n = bytesPerSample_(variable) % pass variable

if ~ischar(vcDataType), vcDataType = class_(vcDataType); end
    
switch lower(vcDataType)
    case {'char', 'byte', 'int8', 'uint8', 'uchar'}
        n = 1;    
    case {'int16', 'uint16'}
        n = 2;
    case {'single', 'float', 'int32', 'uint32', 'float32'}
        n = 4;
    case {'double', 'int64', 'uint64'}
        n = 8;
    otherwise
        n = [];
        fprintf(2, 'Unsupported data type: %s\n', vcDataType);
end
end %func


%--------------------------------------------------------------------------
function [vc, fGpu] = class_(vr)
% Return the class for GPU or CPU arrays 
if isempty(vr)
    vc = class(gather_(vr));
else
    vc = class(gather_(vr(1)));
end
if nargout>=2, fGpu = isGpu_(vr); end
end %func


%--------------------------------------------------------------------------
function [vrRho, vrDelta, miKnn, viNneigh, memory_sort] = sort_disk1_(S0, P, S_drift)

[cviSpk_site, cviSpk2_site] = calc_cviSpk_site_(S0, P);
% viDrift_spk = get_viDrift_spk_(S_drift);
[viLim_drift, mlDrift] = get_(S_drift, 'viLim_drift', 'mlDrift');
vcFile_prm = P.vcFile_prm;
[viSite_spk, vnSpk_load, nLoads] = struct_get_(S0, 'viSite_spk', 'vnSpk_load', 'nLoads');
nSpk = numel(viSite_spk);
nSites = size(P.miSites,2);
knn = get_set_(P, 'knn', 30);
fParfor = get_set_(P, 'fParfor', true) && nSites>1;
nPcPerChan = get_set_(P, 'nPcPerChan', 0);
miKnn = [];

% Calculate Rho
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
viNneigh = zeros(nSpk, 1, 'int64');
mlPc = calc_mlPc_(nPcPerChan, S0.dimm_fet);
P_sort = struct_add_(P, mlDrift, fParfor, mlPc);
P_sort = struct_merge_(P_sort, ...
    struct_copy_(S0, 'type_fet', 'dimm_fet', 'nLoads', 'ccviSpk_site_load', 'ccviSpk_site2_load'));

fprintf('Calculating Rho (sort_disk_)\n\t'); t1=tic;
[cvrRho_site, cvrDelta_site, cviNneigh_site] = deal(cell(nSites,1));
if fParfor   
    try
        parfor iSite = 1:nSites 
            try
                [viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
                if isempty(viSpk1), continue; end
                mrFet12 = load_fet_drift_(P_sort, iSite, P_sort.mlPc);
                viSpk12 = int64([viSpk1; viSpk2]); 
                [cvrRho_site{iSite}, miKnn_] = rho_knn_(mrFet12, viSpk12, ...
                    discretize(viSpk12, viLim_drift), numel(viSpk1), P_sort);
                save_miKnn_site_(vcFile_prm, iSite, miKnn_);
            catch
                fParfor = 0;
            end
        end
    catch
        fParfor = 0;
    end
end
for iSite = 1:nSites         
    [viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
    if isempty(viSpk1), continue; end
    if isempty(cvrRho_site{iSite})
        mrFet12 = load_fet_drift_(P_sort, iSite, P_sort.mlPc); 
        viSpk12 = int64([viSpk1; viSpk2]); 
        [cvrRho_site{iSite}, miKnn_] = rho_knn_(mrFet12, viSpk12, ...
            discretize(viSpk12, viLim_drift), numel(viSpk1), P_sort);
        save_miKnn_site_(vcFile_prm, iSite, miKnn_);
    end
    vrRho(viSpk1) = cvrRho_site{iSite};
end
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), P_sort.fGpu, fParfor);

fprintf('Calculating Delta (sort_disk_)\n\t'); t1=tic;
if fParfor
    try
        parfor iSite = 1:nSites 
            try
                [viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
                if isempty(viSpk1), continue; end
                mrFet12 = load_fet_drift_(P_sort, iSite, P_sort.mlPc);  
                viSpk12 = int64([viSpk1; viSpk2]); 
                [cvrDelta_site{iSite}, cviNneigh_site{iSite}] = ...
                    delta_knn_(mrFet12, viSpk12, [vrRho(viSpk1); vrRho(viSpk2)], ...
                        discretize(viSpk12, viLim_drift), numel(viSpk1), P_sort);
            catch
                fParfor = 0;
            end
        end
    catch
        fParfor = 0;
    end
end
for iSite = 1:nSites 
    [viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
    if isempty(viSpk1), continue; end
    if isempty(cvrDelta_site{iSite}) && isempty(cviNneigh_site{iSite})
        mrFet12 = load_fet_drift_(P_sort, iSite, P_sort.mlPc);  
        viSpk12 = int64([viSpk1; viSpk2]); 
        [cvrDelta_site{iSite}, cviNneigh_site{iSite}] = ...
            delta_knn_(mrFet12, viSpk12, [vrRho(viSpk1); vrRho(viSpk2)], ...
                discretize(viSpk12, viLim_drift), numel(viSpk1), P_sort);
    end
    [vrDelta(viSpk1), viNneigh(viSpk1)] = deal(cvrDelta_site{iSite}, cviNneigh_site{iSite});
end
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), P_sort.fGpu, fParfor);
memory_sort = memory_matlab_();
end %func


%--------------------------------------------------------------------------
function save_miKnn_site_(vcFile_prm, iSite1, miKnn1)
vcFile_knn_site1 = strrep(vcFile_prm, '.prm', sprintf('_knn_site_%d.irc', iSite1));
write_bin_(vcFile_knn_site1, miKnn1);
end %func


%--------------------------------------------------------------------------
function [miKnn1, vl] = load_miKnn_site_(P, iSite1, viSite_spk, viSpk1)
if nargin<4, viSpk1 = []; end

vcFile_knn_site1 = strrep(P.vcFile_prm, '.prm', sprintf('_knn_site_%d.irc', iSite1));
viSpk_site1 = find(viSite_spk==iSite1);
dimm_knn1 = [get_set_(P, 'knn', 30), numel(viSpk_site1)];
miKnn1 = load_bin_(vcFile_knn_site1, 'int64', dimm_knn1);
if ~isempty(viSpk1)
    [vl, vi] = ismember(viSpk1, viSpk_site1);
    miKnn1 = miKnn1(:, vi(vl));
end
end %func


%--------------------------------------------------------------------------
% load miKnn for given list of spikes
function miKnn1 = load_miKnn_spk_(P, viSite_spk, viSpk1)
% usage
% ----
% miKnn1 = load_miKnn_spk_(P, viSite_spk, viSpk1)
% cmiKnn1 = load_miKnn_spk_(P, viSite_spk, cviSpk1)

fCell = iscell(viSpk1);
if fCell
    vn_cell = cellfun(@numel, viSpk1);
    viSpk1 = cell2mat(viSpk1);
end
viSite1 = viSite_spk(viSpk1);
nSites = max(viSite1);
miKnn1 = zeros(get_set_(P, 'knn', 30), numel(viSpk1), 'int64');
for iSite = 1:nSites
    viiSpk1 = find(viSite1 == iSite);
    if isempty(viiSpk1), continue; end
    miKnn1(:,viiSpk1) = load_miKnn_site_(P, iSite, viSite_spk, viSpk1(viiSpk1));
end
if fCell
    cmiKnn1 = cell(numel(vn_cell), 1);
    vi_lim = [0; cumsum(vn_cell)]+1;
    for iCell = 1:numel(vn_cell)
        vi1 = vi_lim(iCell):vi_lim(iCell+1)-1;
        cmiKnn1{iCell} = miKnn1(:, vi1);
    end
    miKnn1 = cmiKnn1;
end
end %func


%--------------------------------------------------------------------------
% load fet from files
function [out1, out2] = load_fet_drift_(S_in, iSite, mlPc)
if nargin<3, mlPc = []; end

[type_fet, dimm_fet, nLoads] = struct_get_(S_in, 'type_fet', 'dimm_fet', 'nLoads');
vcFile_prm = get_(S_in, 'vcFile_prm');
if isempty(vcFile_prm)
    try
        vcFile_prm = S_in.P.vcFile_prm;
    catch
        error('load_fet_drift_: P.vcFile_prm not defined');
    end
end
vcFile_prm_ = strrep(vcFile_prm, '.prm', '');
fh_sum = @(x,y)sum(cellfun(@numel, x(1:y-1)));
vnSpk1_load = cellfun(@(x)numel(x{iSite}), S_in.ccviSpk_site_load);
vnSpk2_load = cellfun(@(x)numel(x{iSite}), S_in.ccviSpk_site2_load);
bytes_per_spk = bytesPerSample_(type_fet) * dimm_fet(1) * dimm_fet(2);
vnBytes_offset1_load = cellfun(@(x)fh_sum(x, iSite), S_in.ccviSpk_site_load) * bytes_per_spk;
vnBytes_offset2_load = cellfun(@(x)fh_sum(x, iSite), S_in.ccviSpk_site2_load) * bytes_per_spk;
csFiles_in1 = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet_site_%d.irc',x)], 1:nLoads);
csFiles_in2 = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet2_site_%d.irc',x)], 1:nLoads);
if ~isempty(mlPc)
    nFeatures = sum(mlPc(:));
    tr2mr_ = @(x)reshape(x, [], size(x,3));
    fet2pc_ = @(x)reshape(x(mlPc(:),:), [nFeatures, size(x,2)]); 
    if nargout==1
        mrFet12 = zeros(nFeatures, sum(vnSpk1_load)+sum(vnSpk2_load), type_fet);
    else
        mrFet1 = zeros(nFeatures, sum(vnSpk1_load), type_fet);
        mrFet2 = zeros(nFeatures, sum(vnSpk2_load), type_fet);
        mrFet12 = [];
    end
else
    trPc1 = zeros(dimm_fet(1), dimm_fet(2), sum(vnSpk1_load), type_fet);
    trPc2 = zeros(dimm_fet(1), dimm_fet(2), sum(vnSpk2_load), type_fet);
end
for iLoad = 1:numel(csFiles_in1)
    viSpk1 = (1:vnSpk1_load(iLoad)) + sum(vnSpk1_load(1:iLoad-1));
    viSpk2 = (1:vnSpk2_load(iLoad)) + sum(vnSpk2_load(1:iLoad-1));
    dimm_fet1 = [dimm_fet(1), dimm_fet(2), vnSpk1_load(iLoad)];
    dimm_fet2 = [dimm_fet(1), dimm_fet(2), vnSpk2_load(iLoad)];
    trPc_load1 = load_bin_(csFiles_in1{iLoad}, type_fet, dimm_fet1, vnBytes_offset1_load(iLoad));
    trPc_load2 = load_bin_(csFiles_in2{iLoad}, type_fet, dimm_fet2, vnBytes_offset2_load(iLoad));    
    if ~isempty(mlPc)
        if isempty(mrFet12)
            mrFet1(:,viSpk1) = fet2pc_(tr2mr_(trPc_load1));
            if ~isempty(trPc_load2), mrFet2(:,viSpk2) = fet2pc_(tr2mr_(trPc_load2)); end            
        else
            mrFet12(:,viSpk1) = fet2pc_(tr2mr_(trPc_load1));
            if ~isempty(trPc_load2)
                mrFet12(:,sum(vnSpk1_load)+viSpk2) = fet2pc_(tr2mr_(trPc_load2)); 
            end
        end
    else
        trPc1(:,:,viSpk1) = trPc_load1;
        if ~isempty(trPc_load2), trPc2(:,:,viSpk2) = trPc_load2; end
    end
end %for
if ~isempty(mlPc)
    if isempty(mrFet12)
        [out1, out2] = deal(mrFet1, mrFet2);
    else
        out1 = mrFet12;
    end
else
    [out1, out2] = deal(trPc1, trPc2);
end
end %func


%--------------------------------------------------------------------------
function [vrRho, vrDelta, miKnn, viNneigh, memory_sort] = sort_ram_(S0, P, trPc_spk, S_drift)

viDrift_spk = get_viDrift_spk_(S_drift);
mlDrift = S_drift.mlDrift;

nSpk = numel(S0.viSite_spk);
nSites = size(P.miSites,2);
knn = get_set_(P, 'knn', 30);
[cviSpk_site, cviSpk2_site] = calc_cviSpk_site_(S0, P);
fParfor = get_set_(P, 'fParfor', true) && nSites>1;

% Calculate Rho
P_sort = struct_set_(P, 'mlDrift', mlDrift, 'fSort_drift', 1, 'fParfor', 0);
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
viNneigh = zeros(nSpk, 1, 'int64');
miKnn = zeros(knn, nSpk, 'int64');
[cvrRho, cmiKnn, cvrDelta, cviNneigh] = deal(cell(nSites,1));
% send jobs

if fParfor 
    gcp_ = gcp();
else
    gcp_ = [];
end

fprintf('Calculating Rho (sort_ram_)\n\t'); t1=tic;
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1] = ...
        pc2fet_site2_(trPc_spk, cviSpk_site, cviSpk2_site, viDrift_spk, iSite);
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
[cvrRho, cmiKnn] = deal([]); % clear memory
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), P_sort.fGpu, fParfor);


% Calculate Delta
fprintf('Calculating Delta (sort_ram_)\n\t'); t2=tic;
% send jobs
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1] = ...
        pc2fet_site2_(trPc_spk, cviSpk_site, cviSpk2_site, viDrift_spk, iSite);
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
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t2), P_sort.fGpu, fParfor);
memory_sort = memory_matlab_();
end %func


%--------------------------------------------------------------------------
function [trPc, nFeatures] = get_pc_sort_(S0, P)

% dimm_fet = size(trPc);
[type_fet, dimm_fet, nLoads, vnSpk_load] = ...
    struct_get_(S0, 'type_fet', 'dimm_fet', 'nLoads', 'vnSpk_load');
[nSites_spk, nSpk] = deal(dimm_fet(2), dimm_fet(3));
nPcPerChan = get_set_(P, 'nPcPerChan', 0);
mlPc = calc_mlPc_(nPcPerChan, S0.dimm_fet);
nC_max = read_cfg_('nC_max', 0);
nFeatures = sum(mlPc(:));
if nFeatures > nC_max
    fprintf('get_pc_sort_: feature trimmed %d->%d\n', nFeatures, nC_max);
    nFeatures = nC_max;
    vi_ = find(mlPc(:));
    mlPc(vi_(nC_max+1:end)) = false;
end

% return empty if trPc_spk should be loaded partially
if get_set_(P, 'fLargeRecording', 0)
    trPc=[]; 
    return; 
end
vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');
csFiles_in1 = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet_%d.irc',x)], 1:nLoads);
try
    if isempty(get_(S0, 'viSite2_spk'))
        trPc = zeros([nFeatures, 1, nSpk], type_fet);
        csFiles_in2 = [];
    else
        trPc = zeros([nFeatures, 2, nSpk], type_fet);
        csFiles_in2 = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet2_%d.irc',x)], 1:nLoads);
    end
catch
    fprintf('get_pc_sort_: memory out, switched to `fLargeRecording=1`\n');
    trPc=[]; 
    return; 
end

fprintf('get_pc_sort_...'); t1=tic;

tr2mr_ = @(x)reshape(x, [], size(x,3));
pc2fet_ = @(x)reshape(x(mlPc(:),:), [sum(mlPc(:)), 1, size(x,2)]);
loadfet_ = @(csFiles, iLoad)pc2fet_(tr2mr_(...
    load_bin_(csFiles{iLoad}, type_fet, ...
        [dimm_fet(1), dimm_fet(2), vnSpk_load(iLoad)])));

iSpk_offset = 0;
for iLoad = 1:nLoads
    nSpk1 = vnSpk_load(iLoad);
    viSpk1 = (1:nSpk1) + iSpk_offset;
    trPc(:,1,viSpk1) = loadfet_(csFiles_in1, iLoad);
    if ~isempty(csFiles_in2)
        trPc(:,2,viSpk1) = loadfet_(csFiles_in2, iLoad);
    end
    iSpk_offset = viSpk1(end);
end
fprintf(' took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [mlPc, vnPc_site] = calc_mlPc_(nPcPerChan, dimm_fet)
[nPc_spk, nSites_spk] = deal(dimm_fet(1), dimm_fet(2));

if isempty(nPcPerChan), nPcPerChan = 0; end
if numel(nPcPerChan) == 1
    if nPcPerChan == 0
        vnPc_site = nPc_spk:-1:1;
    else
        vnPc_site = repmat(nPcPerChan, [1, nSites_spk]);
    end
else
    vnPc_site = nPcPerChan;
end

% determine miPc
mlPc = false(nPc_spk, nSites_spk);
nSites1 = min(nSites_spk, numel(vnPc_site));
for iSite = 1:nSites1
    nPc1 = vnPc_site(iSite);
    mlPc(1:nPc1, iSite) = true;
end
% mlPc(1,:) = true;
end %func


%--------------------------------------------------------------------------
function mrFet = load_pc2fet_site1_(fid_r, cviSpk_site, iSite, P_)
% skip, load, and return 
[mlPc, dimm_fet, type_fet] = struct_get_(P_, 'mlPc', 'dimm_fet', 'type_fet');
vnSpk_site = cellfun(@numel, cviSpk_site);
if iSite==1
    iOffset = 0;    
else
    iOffset = sum(vnSpk_site(1:iSite-1));
end
nSpk1 = vnSpk_site(iSite);
nBytes_skip = iOffset * bytesPerSample_(type_fet) * prod(dimm_fet(1:2));
fseek(fid_r, nBytes_skip, 'bof');
trPc = load_bin_(fid_r, type_fet, [dimm_fet(1)*dimm_fet(2), nSpk1]);
mrFet = reshape(trPc(mlPc(:),:), [sum(mlPc(:)), nSpk1]); 
end %func


%--------------------------------------------------------------------------
function [mrFet12, viSpk12, viDrift12, n1] = pc2fet_site2_(trPc_spk, cviSpk_site, cviSpk2_site, viDrift_spk, iSite)
% decide whether to use 1, 2, or 3 features

if isempty(cviSpk_site{iSite})
    [mrFet12, viSpk12, viDrift12] = deal([]);
    n1 = 0;
    return;
end
% [nSites_fet, miSites] = struct_get_(P, 'nSites_fet', 'miSites');

[viSpk1, viSpk2] = deal(cviSpk_site{iSite}, cviSpk2_site{iSite});
viSpk12 = [viSpk1; viSpk2];
n1 = numel(viSpk1);

if isempty(trPc_spk)
    mrFet12 = [];
elseif isstruct(trPc_spk)
    mrFet12 = trPc_spk;
    mrFet12.iSite = iSite;
elseif isempty(viSpk2)
    mrFet12 = squeeze_(trPc_spk(:,1,viSpk1),2);
else
    mrFet12 = [squeeze_(trPc_spk(:,1,viSpk1),2), squeeze_(trPc_spk(:,2,viSpk2),2)];
end
        
if isempty(viDrift_spk)
    viDrift12 = [];
else
    viDrift12 = viDrift_spk(viSpk12);
end    
end %func


%--------------------------------------------------------------------------
function vrRho1 = rho_knn_disk_(iSite, viSpk1, viSpk2, S_fet)
% S_fet contains {type_fet, dimm_fet, nLoads, vcFile_prm, mlPc, mlDrift, viLim_drift}
    
[mlDrift, viLim_drift, vcFile_prm] = ...
    struct_get_(S_fet, 'mlDrift', 'viLim_drift', 'vcFile_prm');
n1 = numel(viSpk1);
nDrift = size(mlDrift,1);
knn = get_set_(S_fet, 'knn', 30);
[cvrRho_drift1, cmiKnn_drift1] = deal(cell(1, nDrift));
cviiSpk1_drift = vi2cell_(discretize(viSpk1, viLim_drift), nDrift);
cviiSpk2_drift = vi2cell_(discretize(viSpk2, viLim_drift), nDrift);

fParfor = get_set_(S_fet, 'fParfor', 0);
if fParfor
    try
        for iDrift = 1:nDrift    
            [cvrRho_drift1{iDrift}, cmiKnn_drift1{iDrift}] = ...
                rho_knn_disk_aux_(iSite, viSpk1, viSpk2, ...
                    cviiSpk1_drift{iDrift}, cviiSpk2_drift{iDrift}, S_fet);
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iDrift = 1:nDrift    
        [cvrRho_drift1{iDrift}, cmiKnn_drift1{iDrift}] = ...
            rho_knn_disk_aux_(iSite, viSpk1, viSpk2, ...
                cviiSpk1_drift{iDrift}, cviiSpk2_drift{iDrift}, S_fet);
    end    
end

% Combine miKnn1 and save
vrRho1 = zeros(n1, 1, 'single');
miKnn1 = zeros(knn, n1, 'int64');
for iDrift = 1:nDrift
    vii1_ = cviiSpk1_drift{iDrift};
    [vrRho1(vii1_), miKnn1(:,vii1_)] = ...
        deal(cvrRho_drift1{iDrift}, cmiKnn_drift1{iDrift});    
end
save_miKnn_site_(vcFile_prm, iSite, miKnn1);
end %func


%--------------------------------------------------------------------------
function [vrRho1, miKnn1] = rho_knn_disk_aux_(iSite, viSpk1, viSpk2, vii1, vii2, S_fet)
% [vii1_, vii2_] = deal(cviiSpk1_drift{iDrift}, cviiSpk2_drift{iDrift});
if isempty(vii1), [vrRho1, miKnn1] = deal([]); return; end    
knn = get_set_(S_fet, 'knn', 30);
mrFet1_ = load_fet2pc_site_(S_fet, 1, iSite, viSpk1(vii1));    
mrFet2_ = load_fet2pc_site_(S_fet, 2, iSite, viSpk2(vii2));    
[vrRho1, ~, miKnn1] = cuda_knn_([mrFet1_, mrFet2_], mrFet1_, S_fet);
vrRho1 = gather_(1./vrRho1);

% save miKnn
viSpk12_ = int64([viSpk1(vii1); viSpk2(vii2)]);
miKnn1 = viSpk12_(gather_(miKnn1));    
n2 = size(miKnn1,1);
if n2 < knn       
    miKnn1 = [miKnn1; repmat(miKnn1(end,:), [knn-n2, 1])];
end
end %func


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1] = delta_knn_disk_(iSite, viSpk1, viSpk2, vrRho, S_fet)
% S_fet contains {type_fet, dimm_fet, nLoads, vcFile_prm, mlPc, mlDrift, viLim_drift}
    
[mlDrift, viLim_drift, vcFile_prm] = ...
    struct_get_(S_fet, 'mlDrift', 'viLim_drift', 'vcFile_prm');
n1 = numel(viSpk1);
nDrift = size(mlDrift,1);
[cvrDelta_drift1, cviNneigh_drift1] = deal(cell(1, nDrift));
cviiSpk1_drift = vi2cell_(discretize(viSpk1, viLim_drift), nDrift);
cviiSpk2_drift = vi2cell_(discretize(viSpk2, viLim_drift), nDrift);
[vrRho1, vrRho2] = deal(vrRho(viSpk1), vrRho(viSpk2));
fParfor = get_set_(S_fet, 'fParfor', 0);
if fParfor
    try
        for iDrift = 1:nDrift    
            [cvrDelta_drift1{iDrift}, cviNneigh_drift1{iDrift}] = ...
                delta_knn_disk_aux_(iSite, viSpk1, viSpk2, ...
                    cviiSpk1_drift{iDrift}, cviiSpk2_drift{iDrift}, vrRho1, vrRho2, S_fet);
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iDrift = 1:nDrift    
        [cvrDelta_drift1{iDrift}, cviNneigh_drift1{iDrift}] = ...
            delta_knn_disk_aux_(iSite, viSpk1, viSpk2, ...
                cviiSpk1_drift{iDrift}, cviiSpk2_drift{iDrift}, vrRho1, vrRho2, S_fet);
    end    
end

% Combine miKnn1 and save
vrDelta1 = zeros(n1, 1, 'single');
viNneigh1 = zeros(n1, 1, 'int64');
for iDrift = 1:nDrift
    vii1_ = cviiSpk1_drift{iDrift};
    [vrDelta1(vii1_), viNneigh1(vii1_)] = ...
        deal(cvrDelta_drift1{iDrift}, cviNneigh_drift1{iDrift});    
end
end %func


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1] = delta_knn_disk_aux_(iSite, viSpk1, viSpk2, vii1, vii2, vrRho1, vrRho2, S_fet)
% [vii1_, vii2_] = deal(cviiSpk1_drift{iDrift}, cviiSpk2_drift{iDrift});
if isempty(vii1), [vrDelta1, viNneigh1] = deal([]); return; end    

mrFet1_ = load_fet2pc_site_(S_fet, 1, iSite, viSpk1(vii1));    
mrFet2_ = load_fet2pc_site_(S_fet, 2, iSite, viSpk2(vii2));    
vrRho12 = [vrRho1(vii1); vrRho2(vii2)];
[vrDelta1, viNneigh1] = cuda_delta_knn_([mrFet1_, mrFet2_], vrRho12, [], numel(vii1), S_fet);
vrDelta1 = gather_(vrDelta1);

% save miKnn
viSpk12_ = int64([viSpk1(vii1); viSpk2(vii2)]);
viNneigh1 = viSpk12_(gather_(viNneigh1));    
end %func


%--------------------------------------------------------------------------
% load specific site and drift
% todo: speed up by using cache
function out1 = load_fet2pc_site_(S_fet, iFet, iSite, viSpk)

[type_fet, dimm_fet, nLoads, vcFile_prm, mlPc, vnSpk_load] = ...
    struct_get_(S_fet, ...
    'type_fet', 'dimm_fet', 'nLoads', 'vcFile_prm', 'mlPc', 'vnSpk_load');
vcFile_prm_ = strrep(vcFile_prm, '.prm', '');
bytes_per_spk = bytesPerSample_(type_fet) * dimm_fet(1) * dimm_fet(2);
fh_sum = @(x,y)sum(cellfun(@numel, x(1:y-1)));
viSpk_offset_load = cumsum([0; vnSpk_load]);
switch iFet
    case 1
        vnSpk_load = cellfun(@(x)numel(x{iSite}), S_fet.ccviSpk_site_load);
        vnBytes_offset_load = cellfun(@(x)fh_sum(x, iSite), S_fet.ccviSpk_site_load) * bytes_per_spk;
        csFiles_fet = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet_site_%d.irc',x)], 1:nLoads);
        cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site_load);
    case 2
        vnSpk_load = cellfun(@(x)numel(x{iSite}), S_fet.ccviSpk_site2_load);
        csFiles_fet = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet2_site_%d.irc',x)], 1:nLoads);
        vnBytes_offset_load = cellfun(@(x)fh_sum(x, iSite), S_fet.ccviSpk_site2_load) * bytes_per_spk;
        cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site2_load);
end
if ~isempty(mlPc)
    nFeatures = sum(mlPc(:));
    tr2mr_ = @(x)reshape(x, [], size(x,3));
    fet2pc_ = @(x)reshape(x(mlPc(:),:), [nFeatures, size(x,2)]); 
    mrFet1 = zeros(nFeatures, numel(viSpk), type_fet);
else
    trPc1 = zeros(dimm_fet(1), dimm_fet(2), numel(viSpk), type_fet);
end
iOffset_spk = 0;
for iLoad = 1:numel(csFiles_fet)
    [vl_, vi_] = ismember(viSpk, cviSpk_load{iLoad} + viSpk_offset_load(iLoad));
    if ~any(vl_), continue; end
    viiSpk_ = vi_(vl_);
    viSpk1 = (1:numel(viiSpk_)) + iOffset_spk;
    dimm_fet1 = [dimm_fet(1), dimm_fet(2), vnSpk_load(iLoad)];
    trPc_load1 = load_bin_(csFiles_fet{iLoad}, type_fet, dimm_fet1, vnBytes_offset_load(iLoad));    
    if ~isempty(mlPc)
        mrFet1(:,viSpk1) = fet2pc_(tr2mr_(trPc_load1(:,:,viiSpk_)));
    else
        trPc1(:,:,viSpk1) = trPc_load1(:,:,viiSpk_);
    end
    iOffset_spk = viSpk1(end);
end
if ~isempty(mlPc)
    out1 = mrFet1;
else
    out1 = trPc1;
end
end %func
   
    
%--------------------------------------------------------------------------
function [vrRho1, miKnn1, fGpu] = rho_knn_(mrFet12, viSpk12, viDrift12, n1, P)
% mrFet12: matrix or struct

[fGpu, mlDrift] = struct_get_(P, 'fGpu', 'mlDrift');
if isempty(mrFet12)
    [vrRho1, miKnn1] = deal([]);
    return;
end

nT_drift = size(mlDrift,1);
viDrift_spk1 = viDrift12(1:n1);
knn = get_set_(P, 'knn', 30);
cvi1 = vi2cell_(viDrift_spk1, nT_drift);
vrRho1 = zeros(n1, 1, 'single');
miKnn1 = zeros(knn, n1, 'int64');
for iDrift = 1:nT_drift
    vi1 = cvi1{iDrift};
    if isempty(vi1), continue; end    
    vi2 = find(mlDrift(viDrift12, iDrift));
    [vrRho_, fGpu, miKnn_] = cuda_knn_(mrFet12(:,vi2), mrFet12(:,vi1), P);
    vrRho1(vi1) = gather_(1./vrRho_);
    
    miKnn_ = vi2(gather_(miKnn_));    
    n2 = size(miKnn_,1);
    if n2 < knn       
        miKnn_ = [miKnn_; repmat(miKnn_(end,:), [knn-n2, 1])];
    end
    miKnn1(:,vi1) = viSpk12(miKnn_);
end
end %func


%--------------------------------------------------------------------------
function [vrKnn, fGpu, miKnn] = cuda_knn_(mrFet2, mrFet1, P)
% it computes the index of KNN

persistent CK

knn = get_set_(P, 'knn', 30);
[CHUNK, nC_max, nThreads] = deal(4, P.nC_max, 512); % tied to cuda_knn_index.cu
[n2, n1, nC] = deal(size(mrFet2,2), size(mrFet1,2), size(mrFet1,1));
fGpu = P.fGpu && nC <= nC_max;
knn = min(knn, n2);

if fGpu    
    for iRetry = 1:2
        try
            gmrFet2 = gpuArray_(mrFet2, P.fGpu);
            gmrFet1 = gpuArray_(mrFet1, P.fGpu);
            if isempty(CK)
                CK = parallel.gpu.CUDAKernel('cuda_knn_index.ptx','cuda_knn_index.cu'); % auto-compile if ptx doesn't exist
                CK.ThreadBlockSize = [nThreads, 1];          
                CK.SharedMemorySize = 4 * CHUNK * (nC_max + nThreads*2); % @TODO: update the size
            end
            CK.GridSize = [ceil(n1 / CHUNK / CHUNK), CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]    
            vrKnn = zeros([n1, 1], 'single', 'gpuArray');
            vnConst = int32([n2, n1, nC, knn]);            
            miKnn = zeros([knn, n1], 'int32', 'gpuArray'); 
            [vrKnn, miKnn] = feval(CK, vrKnn, miKnn, gmrFet2, gmrFet1, vnConst);
            return;
        catch % use CPU, fail-safe
            CK = [];
            fGpu = 0; % using CPU
        end
    end
end
if ~fGpu    
    [vrKnn, miKnn] = knn_cpu_(mrFet2, mrFet1, knn);
end
end %func


%--------------------------------------------------------------------------
% 10/16/2018 JJJ: Using native function in matlab (pdist2 function)
% 9/20/2018 JJJ: Memory-optimized knn
function [vrKnn, miKnn] = knn_cpu_(mrFet2, mrFet1, knn)

nStep_knn = 1000;
mrFet2_T = mrFet2';
n1 = size(mrFet1,2);
miKnn = int32(zeros([knn, n1], 'like', mrFet1));
vrKnn = zeros([n1, 1], 'like', mrFet1);
for i1 = 1:nStep_knn:n1
    vi1_ = i1:min(i1+nStep_knn-1, n1);
    [mrKnn1, miKnn(:,vi1_)] = pdist2(mrFet2_T, mrFet1(:,vi1_)', 'euclidean', 'Smallest', knn);
    vrKnn(vi1_) = mrKnn1(end,:);
end %for
end %func


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1, fGpu] = delta_knn_(mrFet12, viSpk12, vrRho12, viDrift12, n1, P)
if isempty(mrFet12)
    [vrDelta1, viNneigh1] = deal([]);
    return;
end
mlDrift = P.mlDrift;
% fGpu = isGpu_(mrFet12); 
knn = get_set_(P, 'knn', 30);
SINGLE_INF = 3.402E+38;
fGpu = P.fGpu;
nT_drift = size(mlDrift,1);
vrRho12 = vrRho12(:);
[viDrift_spk1, vrRho1] = deal(viDrift12(1:n1)', vrRho12(1:n1)');
[vrDelta1, viNneigh1] = deal([]);
for iDrift = 1:nT_drift
    vi1_ = find(viDrift_spk1==iDrift);
    if isempty(vi1_), continue; end
    vi2_ = find(mlDrift(viDrift12, iDrift));
    [vrDelta1_, viNneigh1_, fGpu] = cuda_delta_knn_(mrFet12, vrRho12, vi2_, vi1_, P);
    if isempty(vrDelta1)
        vrDelta1 = zeros([1, n1], 'like', vrDelta1_);
        viNneigh1 = zeros([1, n1], 'like', viNneigh1_);
    elseif ~isGpu_(vrDelta1) && fGpu
        vrDelta1_ = gather_(vrDelta1_);
        viNneigh1_ = gather_(viNneigh1_);
    end
    vrDelta1(vi1_) = vrDelta1_;
    viNneigh1(vi1_) = viNneigh1_;
end

vrDelta1 = gather_(vrDelta1) .* vrRho1;
viNneigh1 = gather_(viNneigh1);
viNan = find(isnan(vrDelta1) | isinf(vrDelta1));
viNneigh1(viNan) = viNan;
vrDelta1(viNan) = sqrt(SINGLE_INF);
viNneigh1 = viSpk12(viNneigh1);
end %func


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1, fGpu] = cuda_delta_knn_(mrFet, vrRho, vi2, vi1, P)

persistent CK
CHUNK = get_set_(P, 'CHUNK', 16);
nC_max = get_set_(P, 'nC_max', 60);
nThreads = get_set_(P, 'nThreads', 128);
if isempty(vi2) && numel(vi1)==1
    vi2 = 1:size(mrFet,2);
    vi1 = 1:vi1;
end
[n2, n1, nC] = deal(numel(vi2), numel(vi1), size(mrFet,1));
fGpu = P.fGpu && nC <= nC_max;

if fGpu    
    for iRetry = 1:2
        try
            [gmrFet2, gmrFet1, gvrRho2, gvrRho1] = gpuArray_deal_(mrFet(:,vi2), mrFet(:,vi1), vrRho(vi2), vrRho(vi1));
            if isempty(CK)
                CK = parallel.gpu.CUDAKernel('cuda_delta_knn.ptx','cuda_delta_knn.cu'); % auto-compile if ptx doesn't exist
                CK.ThreadBlockSize = [nThreads, 1];          
                CK.SharedMemorySize = 4 * CHUNK * (nC_max + nThreads*2 + 1); % @TODO: update the size
            end
            CK.GridSize = [ceil(n1 / CHUNK / CHUNK), CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]    
            vrDelta1 = zeros([n1, 1], 'single', 'gpuArray'); 
            viNneigh1 = zeros([n1, 1], 'uint32', 'gpuArray'); 
            vnConst = int32([n2, n1, nC]);
            [vrDelta1, viNneigh1] = feval(CK, vrDelta1, viNneigh1, gmrFet2, gmrFet1, gvrRho2, gvrRho1, vnConst);
            viNneigh1 = uint32(vi2(viNneigh1));
            return;
        catch % use CPU, fail-safe
            CK = [];
            fGpu = 0;
        end
    end
end
if ~fGpu
    [vrDelta1, viNneigh1] = delta_knn_cpu_(mrFet, vrRho, vi2, vi1);    
end
end %func


%--------------------------------------------------------------------------
% 9/20/2018 JJJ: Memory-optimized knn for computing delta (dpclus)
function [vrDelta1, viNneigh1] = delta_knn_cpu_(mrFet, vrRho, vi2, vi1)
nStep_knn = 1000;
n1 = numel(vi1);
vrDelta1 = zeros([n1, 1], 'single');
viNneigh1 = zeros([n1, 1], 'uint32');
[mrFet2, mrFet1] = deal(mrFet(:,vi2), mrFet(:,vi1));
[vrRho2, vrRho1] = deal(vrRho(vi2), vrRho(vi1));
fh_dist_ = @(y)bsxfun(@plus, sum(y.^2), bsxfun(@minus, sum(mrFet2.^2)', 2*mrFet2'*y));
%fh_dist_ = @(y)pdist2(y', mrFet2', 'cosine')';
for i1 = 1:nStep_knn:n1
    vi1_ = i1:min(i1+nStep_knn-1, n1);
    mrD_ = fh_dist_(mrFet1(:,vi1_));
    mrD_(bsxfun(@le, vrRho2, vrRho1(vi1_)')) = inf;
%     if true
%         mrD_(vi2(:) > vi1(:)') = inf;
%     end
    [vrDelta1(vi1_), viNneigh1(vi1_)] = min(mrD_);
end
vrDelta1 = sqrt(abs(vrDelta1));
viNneigh1 = uint32(vi2(viNneigh1));
end %func


%--------------------------------------------------------------------------
function varargout = gpuArray_deal_(varargin)
if nargin > nargout
    fGpu = varargin{end}; 
else
    fGpu = 1;
end

for iArg = 1:nargout
    vr_ = varargin{iArg};
    if fGpu
        try
            if ~isa(vr_, 'gpuArray'), vr_ = gpuArray(vr_); end
        catch
            fGpu = 0;
        end
    else
        if isa(vr_, 'gpuArray'), vr_ = gather_(vr_); end
    end
    varargout{iArg} = vr_;
end %for
end %func


%--------------------------------------------------------------------------
function flag = isGpu_(vr)
try
    flag = isa(vr, 'gpuArray');
catch
    flag = 0;
end
end


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
fParfor = get_set_(P, 'fParfor', false);

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
fprintf('Memory use: %0.3f GiB\n', (memory_matlab_()-memory_init)/2^30);
delete_file_fet_(P); % clear fet
[vcFile, vS_load] = readmda_paged_('close'); % close the file
cS_detect{1} = detect_paged_save_(cS_detect{1}, P, 1);    
viSite2Chan = get_(P, 'viSite2Chan');
if ~fDone
    if fParfor
        parfor iLoad = 2:nLoads  % change to for loop for debugging
            S_load1 = vS_load(iLoad);
            mrWav_T1 = load_file_part_(vcFile, S_load1, viSite2Chan);
            S_cache1 = setfield(S_cache, 'nlim_wav1', S_load1.nlim);
            cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache1);  mrWav_T1 = [];
            cS_detect{iLoad} = detect_paged_save_(cS_detect{iLoad}, P, iLoad);        
        end
    else
        for iLoad = 2:nLoads  % change to for loop for debugging
            S_load1 = vS_load(iLoad);
            mrWav_T1 = load_file_part_(vcFile, S_load1, viSite2Chan);
            S_cache1 = setfield(S_cache, 'nlim_wav1', S_load1.nlim);
            cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache1);  mrWav_T1 = [];
            cS_detect{iLoad} = detect_paged_save_(cS_detect{iLoad}, P, iLoad);        
        end
    end
end
S0 = detect_merge_(cS_detect, viOffset_load, P);
fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
runtime_detect = toc(runtime_detect);
memory_detect = memory_matlab_();
S0 = struct_add_(S0, vrThresh_site, mrPv_global, runtime_detect, P, memory_detect, memory_init, nLoads);
fprintf('Detection took %0.1fs and used %0.3f GiB (fParfor=%d, fGpu=%d)\n', ...
    runtime_detect, (memory_detect-memory_init)/2^30, P.fParfor, P.fGpu);
end %func


%--------------------------------------------------------------------------
function delete_file_fet_(P)
vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');
delete([vcFile_prm_, '_fet*.irc']);
end %func


%--------------------------------------------------------------------------
function mrWav_T1 = load_file_part_(vcFile, S_load1, viSite2Chan)
if nargin<3, viSite2Chan=[]; end
fid1 = fopen(vcFile, 'r');
fseek(fid1, S_load1.offset, 'bof');
mrWav_T1 = fread_(fid1, S_load1.dimm, S_load1.vcDataType);      
if ~isempty(viSite2Chan), mrWav_T1 = mrWav_T1(viSite2Chan,:); end
fclose_(fid1);
end %func


%--------------------------------------------------------------------------
function [mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk, viSite_spk, P)
[~, nSites_spk, nSpk] = size(trPc_spk);
switch 1
    case 4, mrW = squeeze_(abs(trPc_spk(1,:,:)),1);
    case 3, mrW = squeeze_(sum(trPc_spk(1:P.nPcPerChan,:,:).^2),1);
    case 2, mrW = squeeze_(sum(trPc_spk.^2),1);
    case 1, mrW = squeeze_(trPc_spk(1,1:P.nSites_fet,:),1) .^2;
end
vrPow_spk = sum(mrW,1)';

% fprintf('Calculating spike positions\n\t'); t1=tic;
toVec_ = @(x)x(:);
switch 1        
    case 1 % centroid       
        trSiteXY_spk = reshape(P.mrSiteXY(P.miSites(1:P.nSites_fet, viSite_spk),:), P.nSites_fet, nSpk, 2);
        mrPos_spk = [sum(mrW .* trSiteXY_spk(:,:,1))', sum(mrW .* trSiteXY_spk(:,:,2))'] ./ vrPow_spk;
%         figure; plot(mrPos_spk(:,1), mrPos_spk(:,2), '.', 'MarkerSize',1); axis([10 60 10 210])
        % figure; plot3(mrPos_spk(:,1), mrPos_spk(:,2), log(mean(mrW)'), '.');
        
    case 2 % 2d peak (non-uniform kernel)
        mrPos_spk = zeros(nSpk, 2, 'single');
        nSites = max(viSite_spk);
        cviSpk_site = vi2cell_(viSite_spk, nSites);
%         cviSpk_site = arrayfun(@(x)find(viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
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
% fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [cviSpk_site, nSites] = vi2cell_(viSite_spk, nSites)
if nargin<2, nSites = []; end
if isempty(nSites), nSites = max(viSite_spk); end
cviSpk_site = arrayfun(@(x)find(viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
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
% ALSO COMPUTE POSITION AND POWER
% save _fet.irc and _fet2.irc
function S0 = detect_merge_(cS_detect, viOffset_load, P)
%    [mrPos_spk, vrPow_spk] = calc_pos_spk_(S0.trPc_spk, S0.viSite_spk, P);
fSave_fet = get_set_(P, 'fSave_fet', 1);
vnSpk_load = cellfun(@(x)numel(x.viSite_spk), cS_detect);
miSpk_load = [0; cumsum(vnSpk_load)];
miSpk_load = [miSpk_load(1:end-1)+1, miSpk_load(2:end)];

nSpk = sum(vnSpk_load);
[viSite_spk, viTime_spk, mrPos_spk] = ...
    deal(zeros(nSpk, 1, 'int32'), zeros(nSpk, 1, 'int64'), zeros(nSpk, 2, 'single'));
[vrPow_spk, vrAmp_spk] = deal(zeros(nSpk, 1, 'single'));
viOffset_load = int64(viOffset_load);
[mrVp_spk, fid_fet, fid_fet2, type_fet, dimm_fet, viSite2_spk, trPc_spk, trPc2_spk] = deal([]);
[ccviSpk_site_load, ccviSpk_site2_load] = deal(cell(size(cS_detect)));
for iLoad = 1:numel(cS_detect)
    S1 = cS_detect{iLoad};
    if isempty(fid_fet) && ~isempty(S1.trPc_spk) && fSave_fet
        fid_fet = fopen(strrep(P.vcFile_prm, '.prm', '_fet.irc'), 'w');
        type_fet = class(S1.trPc_spk);
        dimm_fet = [size(S1.trPc_spk,1), size(S1.trPc_spk,2), nSpk];
    elseif isempty(type_fet) && isempty(dimm_fet)
        type_fet = S1.type_fet;
        dimm_fet = [S1.dimm_fet(1), S1.dimm_fet(2), nSpk];
    end
    viSpk1 = miSpk_load(iLoad,1):miSpk_load(iLoad,2);
    viSite_spk(viSpk1) = S1.viSite_spk;
    viTime_spk(viSpk1) = int64(S1.viTime_spk) + viOffset_load(iLoad);
    vrAmp_spk(viSpk1) = S1.vrAmp_spk;
    if isfield(S1, 'mrPos_spk') && isfield(S1, 'vrPow_spk')
        [mrPos_spk(viSpk1,:), vrPow_spk(viSpk1)] = deal(S1.mrPos_spk, S1.vrPow_spk);    
    else
        [mrPos_spk(viSpk1,:), vrPow_spk(viSpk1)] = calc_pos_spk_(S1.trPc_spk, S1.viSite_spk, P);            
    end
    if ~isempty(fid_fet)
        write_bin_(fid_fet, S1.trPc_spk); 
    elseif ~fSave_fet
        if isempty(trPc_spk), trPc_spk = zeros(dimm_fet, 'single'); end
        trPc_spk(:,:,viSpk1) = S1.trPc_spk;
    end
    ccviSpk_site_load{iLoad} = get_(S1, 'cviSpk_site');
    
    % secondary peak 
    if isempty(fid_fet2) && ~isempty(S1.trPc2_spk) && fSave_fet
        fid_fet2 = fopen(strrep(P.vcFile_prm, '.prm', '_fet2.irc'), 'w');        
    end    
    if ~isempty(get_(S1, 'viSite2_spk'))
        if isempty(viSite2_spk)
            viSite2_spk = zeros(nSpk, 1, 'int32');
        end
        viSite2_spk(viSpk1) = S1.viSite2_spk;
    end
    if ~isempty(fid_fet2)
        write_bin_(fid_fet2, S1.trPc2_spk); 
    elseif ~fSave_fet
        if isempty(trPc2_spk), trPc2_spk = zeros(dimm_fet, 'single'); end
        trPc2_spk(:,:,viSpk1) = S1.trPc2_spk;
    end
    ccviSpk_site2_load{iLoad} = get_(S1, 'cviSpk2_site');
    
end
fclose_(fid_fet);
fclose_(fid_fet2);
S0 = makeStruct_(viSite_spk, viTime_spk, vrAmp_spk, mrVp_spk, ...
         viSite2_spk, trPc_spk, trPc2_spk, type_fet, dimm_fet, ...
         mrPos_spk, vrPow_spk, ccviSpk_site_load, ccviSpk_site2_load, vnSpk_load); 
end %func


%--------------------------------------------------------------------------
% Save trPc_spk and trPc2_spk and remove from the struct
function S_detect = detect_paged_save_(S_detect, P, iLoad)
% usage
% -----
% [S_detect, fid_fet, fid_fet2] = detect_paged_save_(S_detect, P)
%    Open _fet.irc and _fet2.irc files and save
% S_detect = detect_paged_save_(S_detect, P, fid_fet, fid_fet2)
%    Save to files
% S_detect = detect_paged_save_(S_detect, P, iLoad)
%    save to _fet_{iLoad}.irc and _fet2_{iLoad}.irc
% S_detect = detect_paged_save_('merge', P, nLoads)
%    merge _fet_{1:nLoads}.irc and also _fet2_{1:nLoads}.irc

fSave_fet = get_set_(P, 'fSave_fet', 1);
if ~fSave_fet
    S_detect.type_fet = class(S_detect.trPc_spk);
    S_detect.dimm_fet = size(S_detect.trPc_spk);
    return; 
end

vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');
write_bin_([vcFile_prm_, sprintf('_fet_%d.irc', iLoad)], S_detect.trPc_spk);    
S_detect.type_fet = class(S_detect.trPc_spk);
S_detect.dimm_fet = size(S_detect.trPc_spk);
% if get_set_(P, 'fLargeRecording', false)
S_detect.cviSpk_site = save_paged_fet_site_(...
    [vcFile_prm_, sprintf('_fet_site_%d.irc', iLoad)], S_detect.trPc_spk, S_detect.viSite_spk);
% end
S_detect.trPc_spk = [];

% save fet2
if isempty(get_(S_detect, 'trPc2_spk'))
    S_detect.cviSpk2_site = cell(size(S_detect.cviSpk_site));
    return; 
end   
write_bin_([vcFile_prm_, sprintf('_fet2_%d.irc', iLoad)], S_detect.trPc2_spk);    
% if get_set_(P, 'fLargeRecording', false)
S_detect.cviSpk2_site = save_paged_fet_site_(...
    [vcFile_prm_, sprintf('_fet2_site_%d.irc', iLoad)], S_detect.trPc2_spk, S_detect.viSite2_spk);
% end
S_detect.trPc2_spk = [];
end %func


%--------------------------------------------------------------------------
function cviSpk_site = save_paged_fet_site_(vcFile_out, trFet_spk, viSite_spk)
% t1=tic;
[cviSpk_site, nSites] = vi2cell_(viSite_spk);
fid_w = fopen(vcFile_out, 'w');
for iSite = 1:nSites
    write_bin_(fid_w, trFet_spk(:,:,cviSpk_site{iSite}));
end
fclose(fid_w);
% fprintf('save_paged_fet_site_: took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function varargout = cellfun_(varargin)
if nargout == 0
    cellfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==1
    varargout{1} = cellfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==2
    [varargout{1}, varargout{2}] = cellfun(varargin{:}, 'UniformOutput', 0);    
elseif nargout==3
    [varargout{1}, varargout{2}, varargout{3}] = cellfun(varargin{:}, 'UniformOutput', 0);    
else
    error('cellfun_: nargout exceeds 3');
end   
end %func


%--------------------------------------------------------------------------
function varargout = arrayfun_(varargin)
if nargout == 0
    arrayfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==1
    varargout{1} = arrayfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==2
    [varargout{1}, varargout{2}] = arrayfun(varargin{:}, 'UniformOutput', 0);    
elseif nargout==3
    [varargout{1}, varargout{2}, varargout{3}] = arrayfun(varargin{:}, 'UniformOutput', 0);    
else
    error('arrayfun_: nargout exceeds 3');
end   
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
mrWav2 = filter_transpose_(mrWav_T, P);
S_detect = get_spikes_(mrWav2, P, S_cache);
end %func


%--------------------------------------------------------------------------
function S_detect = get_spikes_(mrWav2, P, S_cache)
[vrThresh_site, nlim_wav1, mrPv_global, vrD_global] = ...
    struct_get_(S_cache, 'vrThresh_site', 'nlim_wav1', 'mrPv_global', 'vrD_global');

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
mrVp_spk = [];

% extract spike feaures
if isempty(mrPv_global)
    [mrPv_global, vrD_global] = get_pv_(trWav_spk, P); 
end
trPc_spk = gather_(project_pc_(trWav_spk, mrPv_global, P));

if get_set_(P, 'sort_mode', 1) == 1 && size(trWav_spk,2) > 1
    viSite2_spk = find_site_spk23_(trWav_spk, viSite_spk, P);
    trWav_spk = []; %clear memory
    trWav2_spk = mn2tn_wav_spk2_(mrWav2, viSite2_spk, viTime_spk, P);
    trPc2_spk = gather_(project_pc_(trWav2_spk, mrPv_global, P));
else
    [viSite2_spk, trPc2_spk] = deal([]);
end

% return struct
if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end
[mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk, viSite_spk, P); 
if false % average the first and second peak locations
    if ~isempty(trPc2_spk) && ~isempty(viSite2_spk)
        [mrPos2_spk, vrPow2_spk] = calc_pos_spk_(trPc2_spk, viSite2_spk, P);
        mrPos_spk = (mrPos_spk+mrPos2_spk)/2;
        vrPow_spk = (vrPow_spk+vrPow2_spk)/2;
    end
end
S_detect = makeStruct_(trPc_spk, mrPv_global, viTime_spk, vrAmp_spk, viSite_spk, ...
    mrVp_spk, trPc2_spk, viSite2_spk, vrThresh_site, mrPos_spk, vrPow_spk, vrD_global);
end %func


%--------------------------------------------------------------------------
% 5/28/2019 JJJ: returns a cell of matrix
function trPc_spk = project_pc_(trWav_spk, mrPv, P)
% [mrFet1, mrFet2, mrFet3] = deal([]);
project_ = @(x,y)reshape(x' * reshape(y, size(y,1), []), size(x,2), size(y,2), size(y,3));
trPc_spk = project_(mrPv, trWav_spk); % project 0-shift

% interpolate and shift
if get_set_(P, 'fInterp_fet', 0) == 0, return; end   
switch 1 % 1
    case 0, return;
    case 1, viShift = -1:.5:1;
    case 2, viShift = -1.5:.5:1.5;     
    case 3, viShift = -2:.5:2;       
    case 4, viShift = -1:.25:1;
end
viShift = [0, setdiff(viShift,0)]; % put 0 first
nPc = size(mrPv,2);
trPv_shift = arrayfun(@(x)vr2mr_shift_(mrPv(:,x), viShift), 1:nPc, 'UniformOutput', 0);
trPv_shift = permute(cat(3, trPv_shift{:}),[1,3,2]);

% find best alignment using first channel
mrWav1 = squeeze_(trWav_spk(:,1,:), 2); % first channel
mrPv1_shift = squeeze_(trPv_shift(:,1,:),2); %first pv shift
[~, viMax_spk] = max(abs(mrPv1_shift' * mrWav1));  % project chan1 shift   

for iShift=2:numel(viShift) % start from 2 to exclude shifting 0
    viSpk1 = find(viMax_spk == iShift);
    if isempty(viSpk1), continue; end
    trPc_spk(:,:,viSpk1) = project_(trPv_shift(:,:,iShift), trWav_spk(:,:,viSpk1));
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
function [mrPv1, vrD1] = get_pv_(tr, P)
%tr: nSamples x nSpikes x nChans
MAX_SAMPLE = 10000;        

viSpk_sub = subsample_vr_(1:size(tr,3), MAX_SAMPLE);
switch 3
    case 1, mr1 = squeeze(tr(:, 1, viSpk_sub)); % use peak chan only
    case 2, mr1 = reshape(tr(:, 1:P.nSites_fet, viSpk_sub), size(tr,1), []); 
    case 3, mr1 = reshape(tr(:, :, viSpk_sub), size(tr,1), []); % use all chan
end
nPc_spk = get_set_(P, 'nPc_spk', 9); % # components to compress spike waveforms

switch 2
    case 1
        % mrSpkWav1 = meanSubt_(mrSpkWav1);
        [mrPv1, vrD1] = eig(mr1 * mr1');
        mrPv1 = zscore_(fliplr(mrPv1)); % sort largest first
        vrD1 = flipud(diag(vrD1));
    case 2
        [mrPv1, ~, vrD1] = pca(gather_(mr1)','Centered',false, 'NumComponents', nPc_spk);
end
end %func


%--------------------------------------------------------------------------
function [A,B,r] = canoncorr_(X,Y, fCentered)
if nargin < 2
    error(message('stats:canoncorr:TooFewInputs'));
end
if nargin<3, fCentered = 1; end

[n,p1] = size(X);
if size(Y,1) ~= n
    error(message('stats:canoncorr:InputSizeMismatch'));
elseif n == 1
    error(message('stats:canoncorr:NotEnoughData'));
end
p2 = size(Y,2);

% Center the variables
if fCentered
    X = X - repmat(mean(X,1), n, 1);
    Y = Y - repmat(mean(Y,1), n, 1);
end
% Factor the inputs, and find a full rank set of columns if necessary
[Q1,T11,perm1] = qr(X,0);
rankX = sum(abs(diag(T11)) > eps(abs(T11(1)))*max(n,p1));
if rankX == 0
    error(message('stats:canoncorr:BadData', 'X'));
elseif rankX < p1
    warning(message('stats:canoncorr:NotFullRank', 'X'));
    Q1 = Q1(:,1:rankX); T11 = T11(1:rankX,1:rankX);
end
[Q2,T22,perm2] = qr(Y,0);
rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n,p2));
if rankY == 0
    error(message('stats:canoncorr:BadData', 'Y'));
elseif rankY < p2
    warning(message('stats:canoncorr:NotFullRank', 'Y'));
    Q2 = Q2(:,1:rankY); T22 = T22(1:rankY,1:rankY);
end

% Compute canonical coefficients and canonical correlations.  For rankX >
% rankY, the economy-size version ignores the extra columns in L and rows
% in D. For rankX < rankY, need to ignore extra columns in M and D
% explicitly. Normalize A and B to give U and V unit variance.
d = min(rankX,rankY);
[L,D,M] = svd(Q1' * Q2,0);
A = T11 \ L(:,1:d) * sqrt(n-1);
B = T22 \ M(:,1:d) * sqrt(n-1);
r = min(max(diag(D(:,1:d))', 0), 1); % remove roundoff errs

% Put coefficients back to their full size and their correct order
A(perm1,:) = [A; zeros(p1-rankX,d)];
B(perm2,:) = [B; zeros(p2-rankY,d)];
end %func


%--------------------------------------------------------------------------
function trWav_spk = get_spkwav_(mrWav, viSite_spk, viTime_spk, P)
nSpks = numel(viTime_spk);
% nSites = numel(P.viSite2Chan);
nSites = size(P.miSites,2);
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
% 2019/8/15 JJJ: can handle arrays or cells file handle
function fid = fclose_(fid, fVerbose)
% Sets fid = [] if closed properly
if nargin<2, fVerbose = 0; end
if isempty(fid), return; end
if ischar(fid), return; end
if numel(fid)>1
    for iFile=1:numel(fid)
        if iscell(fid)
            fclose_(fid{iFile}, fVerbose);
        else
            fclose_(fid(iFile), fVerbose);
        end
    end
    return;
end

if ~isempty(fid)
    try
        fclose(fid); 
    catch
        ;
    end
end
fid = [];
if fVerbose, disp('File closed.'); end
end %func


%--------------------------------------------------------------------------
function [out, nlim, fDone] = readmda_paged_(P)
% manage overlap
% [out, nlim, fDone] = readmda_paged_(P) % read header
% [out, nlim, fDone] = readmda_paged_() % read next page
% [vcFile, vS_load] = readmda_paged_('close'): close a file

persistent S_mda fid iLoad nLoads nSamples_load nSamples_last nSamples_pad vcFile P_ viSite2Chan_

if nargin==1
    if ischar(P)
        if strcmpi(P, 'close')
            fclose_(fid);
            vS_load = plan_load_pad_(S_mda, P_);
            [fid, out, nlim] = deal([], vcFile, vS_load);
            return;
        end
    end
    vcFile = P.vcFile;
    switch recording_type_(vcFile)
        case 'mda'
            [S_mda, fid] = readmda_header_(vcFile);
            viSite2Chan_ = [];
        case 'spikeglx'
            [S_mda, fid] = spikeglx_header_(vcFile);      
            viSite2Chan_ = P.viSite2Chan;
        otherwise, error('unsupported file format');
    end
    S_mda.dimm = S_mda.dimm(:)'; % formatting
    assert(S_mda.nBytes_missing<1, 'readmda_paged_: mda file is incomplete');
    iLoad = 0;
    [nLoads, nSamples_load, nSamples_last] = plan_load_(S_mda.nBytes_data, P);
    nSamples_pad = get_set_(P, 'nPad_filt', 300);
    viOffset_load = (0:nLoads-1) * nSamples_load; % return offset
    S_page = makeStruct_(nLoads, nSamples_load, nSamples_last, nSamples_pad, viOffset_load);
    [out, nlim, fDone, P_] = deal(S_page, [], 0, P);    
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
        if ~isempty(viSite2Chan_), out = out(P_.viSite2Chan,:); end
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
% determine fseek, fread bytes, and trim amount
% vector of struct {nBytes_offset, dimm1}
function vS_load = plan_load_pad_(S_mda, P)

[nLoads, nSamples_load, nSamples_last] = plan_load_(S_mda.nBytes_data, P);
nSamples_pad = get_set_(P, 'nPad_filt', 300);
bytes_per_sample = bytesPerSample_(S_mda.vcDataType);
nChans = S_mda.dimm(1);
% nSamples_file = S_mda.nBytes_data / bytes_per_sample / nChans;
vS_load = cell(nLoads, 1);
for iLoad = 1:nLoads
    offset1 = max((nSamples_load * (iLoad - 1) - nSamples_pad) * (bytes_per_sample * nChans), 0) + S_mda.nBytes_header;
    if nLoads == 1 % no padding if single load
        dimm1 = S_mda.dimm;
        nlim1 = [1, dimm1(end)]; 
    elseif iLoad == 1 % first
        dimm1 = [S_mda.dimm(1:end-1), nSamples_load + nSamples_pad];
        nlim1 = [1, nSamples_load];
    elseif iLoad == nLoads % last
        dimm1 = [S_mda.dimm(1:end-1), nSamples_last + nSamples_pad];
        nlim1 = [1, nSamples_last] + nSamples_pad;
    else % middle (double padded)
        dimm1 = [S_mda.dimm(1:end-1), nSamples_load + nSamples_pad*2];
        nlim1 = [1, nSamples_load] + nSamples_pad;
    end
    vS_load{iLoad} = struct('dimm', dimm1, 'nlim', nlim1, 'vcDataType', S_mda.vcDataType, 'offset', offset1);
end %for
vS_load = cell2mat(vS_load);
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
    error('unsupported test mode');
end
if ispc()
    vcDir_in = strrep(vcDir_in, '/', '\');    
    vcDir_in = fullfile('D:\Globus', vcDir_in); % RAID5
%     vcDir_in = fullfile('c:\tmp', vcDir_in); % SSD
elseif isunix()
    vcDir_in = fullfile('~/ceph', vcDir_in);
elseif ismac()
    vcDir_in = fullfile('~/ceph', vcDir_in);
else
    error('get_test_data_: unsupported OS');
end
end %func


%--------------------------------------------------------------------------
function P = makeParam_(vcDir_in, vcDir_out, vcFile_arg)
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end

[vcDir_, vcFile_, vcExt_] = fileparts(vcDir_in);
if strcmpi(vcExt_, '.mda')
    vcFile_raw = vcDir_in;
    vcDir_in = vcDir_;
elseif strcmpi(vcExt_, '.bin')
    vcFile_raw = vcDir_in;
    vcDir_in = vcDir_;
elseif strcmpi(vcExt_, '.dat')
    vcFile_raw = vcDir_in;
    vcDir_in = vcDir_;
else
    vcFile_raw = fullfile(vcDir_in, 'raw.mda');
    vcExt_ = '.mda';
end
if ~exist_file_(vcFile_raw)
    error('%s does not exist\n', vcFile_raw);
end
vcDir_out = fill_dir_out_(vcDir_in, vcDir_out);

% assume there is raw.mda, geom.csv, params.json, firings_true.mda
P = file2struct_(ircpath_(read_cfg_('default_prm', 0)));
P2 = file2struct_(ircpath_(read_cfg_('default2_prm', 0)));
P = struct_merge_(P, P2);

% now only supporting .mda file
P.vcFile = vcFile_raw;
P.vcDir_out = vcDir_out;
switch recording_type_(vcFile_raw)
    case 'mda'
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
        P.viShank_site = deal([]);
    case 'spikeglx'
        S_meta = read_meta_file_(strrep(vcFile_raw, '.bin', '.meta'));
        P.probe_file = fullfile([S_meta.Smeta.vcProbe, '.prb']);
        [P.sRateHz, P.uV_per_bit, P.vcDataType] = deal(S_meta.sRateHz, S_meta.uV_per_bit, S_meta.vcDataType);
        if contains(vcFile_raw, '.imec.')
            P.nChans = (S_meta.nChans+1)/2;
        else
            P.nChans = S_meta.nChans;
        end
        S_prb = load_prb_(P.probe_file);
        [P.mrSiteXY, P.vrSiteHW, P.viShank_site, P.viSite2Chan] = ...
            struct_get_(S_prb, 'mrSiteXY', 'vrSiteHW', 'viShank_site', 'viSite2Chan');
    otherwise
        error('Unsupported recording format');
end

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
function vcType = recording_type_(vcFile_raw)
vcType = '';
[vcDir_, vcFile_, vcExt_] = fileparts(vcFile_raw);
if strcmpi(vcExt_, '.mda')
    vcType = 'mda';
elseif strcmpi(vcExt_, '.bin')
    vcFile_meta = fullfile(vcDir_, [vcFile_, '.meta']);
    if exist_file_(vcFile_meta)
        vcType = 'spikeglx';
    end
end
end %func

%--------------------------------------------------------------------------
function vcDir_out = fill_dir_out_(vcDir_in, vcDir_out)
if isempty(vcDir_out)
    if ~contains(vcDir_out, 'groundtruth')
        vcDir_out = fullfile(vcDir_in, 'irc2');
    else
        vcDir_out = strrep(vcDir_in, 'groundtruth', 'irc2'); 
    end
end
if ~exist_dir_(vcDir_out), mkdir(vcDir_out); end
end %func


%--------------------------------------------------------------------------
function [S_mda, fid_r] = readmda_header_(fname)
fid_r = fopen(fname,'rb');

try
    code=fread(fid_r,1,'int32');
catch
    error('Problem reading file: %s',fname);
end
if (code>0) 
    num_dims=code;
    code=-1;
    nBytes_sample = 4;
else
    nBytes_sample = fread(fid_r,1,'int32');
    num_dims=fread(fid_r,1,'int32');    
end
dim_type_str='int32';
if (num_dims<0)
    num_dims=-num_dims;
    dim_type_str='int64';
end

dimm=zeros(1,num_dims);
for j=1:num_dims
    dimm(j)=fread(fid_r,1,dim_type_str);
end
nBytes_header = ftell(fid_r);

if (code==-1)
    vcDataType = 'single';
elseif (code==-2)
    vcDataType = 'uchar';
elseif (code==-3)
    vcDataType = 'single';
elseif (code==-4)
    vcDataType = 'int16';
elseif (code==-5)
    vcDataType = 'int32';
elseif (code==-6)
    vcDataType = 'uint16';
elseif (code==-7)
    vcDataType = 'double';
elseif (code==-8)
    vcDataType = 'uint32';
else
    vcDataType = ''; % unknown
end %if

% integrity check
nBytes_data = prod(dimm) * bytesPerSample_(vcDataType);
nBytes_file = filesize_(fname);
nBytes_missing = nBytes_data + nBytes_header - nBytes_file;

% deal with overflow case when one dimension exceeds 2^31
if nBytes_missing < 0 && strcmpi(dim_type_str, 'int32') && dimm(end) >= 2^31-1
    nSamples = (nBytes_file-nBytes_header) / bytesPerSample_(vcDataType);
    if num_dims == 1
        dimm = nSamples;
    else
        dimm(end) = round(nSamples / prod(dimm(1:end-1)));
    end
    nBytes_data = prod(dimm) * bytesPerSample_(vcDataType);
    nBytes_missing = nBytes_data + nBytes_header - nBytes_file;
end

S_mda = struct('dimm', dimm, 'vcDataType', vcDataType, ...
    'nBytes_header', nBytes_header, 'nBytes_sample', nBytes_sample, ...
    'nBytes_missing', nBytes_missing, 'nBytes_data', nBytes_data);

if nargout<2, fclose(fid_r); end
end %func


%--------------------------------------------------------------------------
function [S_mda, fid_r] = spikeglx_header_(fname)
S_meta = read_meta_file_(strrep(fname, '.bin', '.meta'));
vcDataType = S_meta.vcDataType;
nBytes_sample = bytesPerSample_(S_meta.vcDataType);
nBytes_data = filesize_(fname);
nSamples = nBytes_data / nBytes_sample;
nChans = S_meta.nChans;
if contains(fname, '.imec.')
    nChans = (nChans+1)/2;
end
dimm = [nChans, floor(nSamples/nChans)];
S_mda = struct('dimm', dimm, 'vcDataType', vcDataType, ...
    'nBytes_header', 0, 'nBytes_sample', nBytes_sample, ...
    'nBytes_missing', 0, 'nBytes_data', nBytes_data);

if nargout>=2, fid_r = fopen(fname,'rb'); end
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
%         P.nTime_clu = min(P.nTime_clu, get_set_(P, 'nBatch_max_drift', 32));
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
function varargout = gather_(varargin)
for i=1:nargin
    varargout{i} = varargin{i};
    if isa(varargin{i}, 'gpuArray')
        try
            varargout{i} = gather(varargin{i});
        catch
            ;
        end
    end
end
end %func


%--------------------------------------------------------------------------
function join_mda_(vcFile_txt, vcFile_out)
if isTextFile_(vcFile_txt)
    csFiles = load_batch_(vcFile_txt);
else
    csFiles = list_files_(vcFile_txt, 1);
end
if isempty(vcFile_out)
    vcFile_out = fullfile(fileparts(csFiles{1}), 'raw.mda');
end
csFiles = setdiff(csFiles, vcFile_out);
fprintf('Joining %d files\n', numel(csFiles)); t1=tic;
fid_w = fopen(vcFile_out, 'w');
for iFile = 1:numel(csFiles)
    fid_w = writemda_fid(fid_w, readmda_(csFiles{iFile}));
    fprintf('\tAppended %d/%d: %s\n', iFile, numel(csFiles), csFiles{iFile});
end
writemda_fid(fid_w, 'close');
fprintf('\tWrote to %s, took %0.1fs\n', vcFile_out, toc(t1));
end %func


%--------------------------------------------------------------------------
function test_mcc_(vcDir_in)
if nargin<1, vcDir_in = ''; end
if isempty(vcDir_in), vcDir_in = 'static'; end
if ~exist_dir_(vcDir_in)
    vcDir_in = get_test_data_(vcDir_in);
end
system(['run_irc ', vcDir_in]);
end %func


%--------------------------------------------------------------------------
function [mr, fGpu] = gpuArray_(mr, fGpu)
if nargin<2, fGpu = 1; end
% fGpu = 0; %DEBUG disable GPU array
if ~fGpu, return; end
try
    if ~isa(mr, 'gpuArray'), mr = gpuArray(mr); end
    fGpu = 1;
catch        
    try % retry after resetting the GPU memory
        gpuDevice(1); 
        mr = gpuArray(mr);
        fGpu = 1;
    catch % no GPU device found            
        fGpu = 0;
    end
end
end


%--------------------------------------------------------------------------
% 2018/6/1: Added support for .prm format
function mnWav = load_bin_(vcFile, vcDataType, dimm, header)
% mnWav = load_bin_(vcFile, dimm, vcDataType)
% mnWav = load_bin_(fid, dimm, vcDataType)
% mnWav = load_bin_(vcFile_prm)
% [Input arguments]
% header: header bytes

if nargin<2, vcDataType = []; end
if nargin<3, dimm = []; end
if nargin<4, header = 0; end
if isempty(vcDataType), vcDataType = 'int16'; end

mnWav = [];
switch numel(dimm)
    case 0, nChans = 1;
    case 1, nChans = dimm;
end
fTranspose_bin = 1;

if ischar(vcFile)
%     if matchFileExt_(vcFile, '.prm')
%         vcFile_prm = vcFile;    
%         if ~exist_file_(vcFile_prm, 1), return; end        
%         P = loadParam_(vcFile_prm);
%         [vcFile, vcDataType, header, nChans, fTranspose_bin] = ...
%             struct_get_(P, 'vcFile', 'vcDataType', 'header_offset', 'nChans', 'fTranspose_bin');
%     end
%     if ~exist_file_(vcFile, 1), return; end
    fid = fopen_(vcFile, 'r');
    if header>0, fseek(fid, header, 'bof'); end
    if isempty(fid) || fid<0, return; end      
    if numel(dimm) < 2
        nBytes_file = filesize_(vcFile) - header;
        nSamples = floor(nBytes_file / bytesPerSample_(vcDataType) / nChans);
        if fTranspose_bin
            dimm = [nChans, nSamples]; %return column
        else
            dimm = [nSamples, nChans]; %return column
        end
    end
else % fid directly passed
    fid = vcFile;
    if isempty(dimm), dimm = inf; end
end
try
    t1 = tic;
    mnWav = fread_(fid, dimm, vcDataType);
    if ischar(vcFile)
        fclose(fid);
%         fprintf('Loading %s took %0.1fs\n', vcFile, toc(t1)); 
    end
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
function mnWav1 = fread_(fid_bin, dimm_wav, vcDataType)
% Get around fread bug (matlab) where built-in fread resize doesn't work

% defensive programming practice
if strcmpi(vcDataType, 'float'), vcDataType = 'single'; end
if strcmpi(vcDataType, 'float32'), vcDataType = 'single'; end
if strcmpi(vcDataType, 'float64'), vcDataType = 'double'; end

try
    if isempty(dimm_wav)
        mnWav1 = fread(fid_bin, inf, ['*', vcDataType]);
    else
        if numel(dimm_wav)==1, dimm_wav = [dimm_wav, 1]; end
        mnWav1 = fread(fid_bin, prod(dimm_wav), ['*', vcDataType]);
        if numel(mnWav1) == prod(dimm_wav)
            mnWav1 = reshape(mnWav1, dimm_wav);
        else
            dimm2 = floor(numel(mnWav1) / dimm_wav(1));
            if dimm2 >= 1
                nSamples1 = dimm_wav(1) * dimm2;
                mnWav1 = reshape(mnWav1(1:nSamples1), dimm_wav(1), dimm2);
            else
                mnWav1 = [];
            end
        end
    end
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
function [fid, nBytes, header_offset] = fopen_(vcFile, vcMode)
if nargin < 2, vcMode = 'r'; end
[~,~,vcExt] = fileparts(vcFile);
try
    switch lower(vcExt)
        case {'.ns5', '.ns2'}, [fid, nBytes, header_offset] = fopen_nsx_(vcFile);
        case '.mda', [fid, nBytes, header_offset] = fopen_mda_(vcFile, vcMode);
        otherwise
            fid = fopen(vcFile, vcMode); 
            if nargout>1
                nBytes = filesize_(vcFile);   
                header_offset = 0;
            end
    end %switch
    if fid==-1, fid = []; end
catch
    disperr_();
    fid = []; 
    nBytes = [];
end
end %func


%--------------------------------------------------------------------------
% 9/24/2019 JJJ: Map index
function S_drift = calc_drift_(S0, P)
if nargin<2, P = []; end
if isempty(P), P = S0.P; end

fprintf('Calculating drift similarity...'); t1 = tic;
nTime_clu = get_set_(P, 'nTime_clu', 4);
nTime_drift = get_set_(P, 'nTime_drift', nTime_clu * 4);
viDrift_spk = []; %ones(1, numel(S0.viSite_spk), 'int64');
if nTime_clu == 1 || nTime_drift == 1 % no drift correlation analysis
    viLim_drift = [1, numel(S0.viSite_spk)];
    [miSort_drift, nTime_drift] = deal(1, 1);    
else
    nAmp_drift = get_set_(P, 'nQuantile_drift', 10);
    %nPos_drift = get_set_(P, 'nPos_drift', 100);
    nPos_drift = numel(P.viSite2Chan); % use number of sites
%     viSite_unique = unique(S0.viSite_spk);
%     nSites = numel(viSite_unique);
    nSpikes = numel(S0.viSite_spk);
    viLim_drift = int64([0, ceil((1:nTime_drift)/nTime_drift*nSpikes)])+1;

    % collect stats
    if isfield(S0, 'mrPos_spk')
        mrPos_spk = S0.mrPos_spk;
        vrAmp_spk = single(S0.vrAmp_spk(:));
    else
        [mrPos_spk, vrAmp_spk] = spk_pos_(S0);
        vrAmp_spk = single(vrAmp_spk(:));                
    end
    vrAmp_quantile = quantile_vr_(vrAmp_spk, (0:nAmp_drift)/nAmp_drift);
    vrPos_spk = single(mrPos_spk(:,2));
    mrCount_drift = zeros(nAmp_drift*nPos_drift, nTime_drift, 'single');
    vrPos_quantile = quantile_vr_(vrPos_spk, (0:nPos_drift)/nPos_drift);
        
    for iDrift = 1:nTime_drift
        viSpk1 = viLim_drift(iDrift):viLim_drift(iDrift+1)-1;
        mn_ = histcounts2(vrAmp_spk(viSpk1), vrPos_spk(viSpk1), vrAmp_quantile, vrPos_quantile);
        mrCount_drift(:,iDrift) = mn_(:);
    end
    
    mrDist_drift = squareform(pdist(mrCount_drift'));
    [mrSort_drift, miSort_drift] = sort(mrDist_drift, 'ascend');
    nSort_drift = ceil(nTime_drift / P.nTime_clu);
    miSort_drift = miSort_drift(1:nSort_drift,:);

    if read_cfg_('fPlot_drift')
        figure; imagesc(mrDist_drift); set(gcf,'Name', P.vcFile_prm);
        figure; imagesc(mrSort_drift); set(gcf,'Name', P.vcFile_prm);
        hold on; plot([0, size(mrSort_drift,1)], repmat(nSort_drift,1,2), 'r-');
    end
end
mlDrift = mi2ml_drift_(miSort_drift); %gpuArray_(mi2ml_drift_(miSort_drift), P.fGpu);
S_drift = makeStruct_(miSort_drift, nTime_drift, viDrift_spk, mlDrift, viLim_drift);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function vr_out = quantile_vr_(vr, vr_q, MAX_SAMPLE)
if nargin<3, MAX_SAMPLE = []; end
if isempty(MAX_SAMPLE), MAX_SAMPLE = 100000; end
vr = subsample_vr_(vr, MAX_SAMPLE);
n = numel(vr);
idx = max(min(round(n * vr_q), n), 1);
vr_out = sort(vr, 'ascend');
vr_out = vr_out(idx);
end %func


%--------------------------------------------------------------------------
function viDrift_spk = get_viDrift_spk_(S_drift)
viDrift_spk = struct_get_(S_drift, 'viDrift_spk');
if ~isempty(viDrift_spk), return; end

viLim_drift = struct_get_(S_drift, 'viLim_drift');
nSpk = viLim_drift(end);
nDrift = numel(viLim_drift) - 1;
viDrift_spk = ones(1, nSpk, 'int64');
for iDrift = 1:nDrift
    vi1 = (viLim_drift(iDrift)+1):viLim_drift(iDrift+1);
    viDrift_spk(vi1) = iDrift;
end
end %func
    

%--------------------------------------------------------------------------
% change the encoding scheme
function mlDrift = mi2ml_drift_(miMember)

[nMembers, nDrift] = size(miMember);
mlDrift = false(nDrift);
[viCol, viRow] = deal(miMember, repmat(1:nDrift, nMembers, 1));
mlDrift(sub2ind(size(mlDrift), viCol(:), viRow(:))) = true;
end %func


%--------------------------------------------------------------------------
function S_clu = postCluster_(S_clu, P, viSite_spk)

if isempty(S_clu), return; end
if isfield(S_clu, 'viClu')
    S_clu = rmfield(S_clu, 'viClu');
end

S_clu.icl = find(S_clu.delta(:) > get_set_(P, 'delta_cut', 1));
% Update P
S_clu.P.min_count = P.min_count;
S_clu.P.delta1_cut = P.delta1_cut;
S_clu.P.rho_cut = P.rho_cut;
S_clu.viClu = [];

if isempty(P.min_count), P.min_count = 0; end
if ~isfield(S_clu, 'viClu'), S_clu.viClu = []; end

fprintf('assigning clusters, nClu:%d\n', numel(S_clu.icl)); t1=tic;

[S_clu.viClu, S_clu.icl] = assignCluster_(S_clu.viClu, S_clu.ordrho, S_clu.nneigh, S_clu.icl);
[S_clu.viClu, S_clu.icl] = dpclus_remove_count_(S_clu.viClu, S_clu.icl, P.min_count);
switch get_set_(P, 'knn_merge_mode', 1)
    case 1, [~, S_clu, nClu_pre] = S_clu_peak_merge1_(S_clu, P, viSite_spk);  % knn overlap merging
    case 2, [~, S_clu, nClu_pre] = S_clu_peak_merge2_(S_clu, P, viSite_spk);  % knn overlap merging
end

S_clu = S_clu_refresh_(S_clu);
fprintf('\n\ttook %0.1fs. Pre-merged %d clusters: %d->%d\n', ...
    toc(t1), nClu_pre - S_clu.nClu, nClu_pre, S_clu.nClu);
end %func


%--------------------------------------------------------------------------
% 9/3/2019 JJJ: run_mode can be an array (sequential merge)
% 9/17/2018 JJJ: merge peaks based on their waveforms
function [viMap, S_clu, nClu_post] = S_clu_peak_merge2_(S_clu, P, viSite_spk) 
t1=tic;
miKnn = get_(S_clu, 'miKnn');
[vrRho_spk, viClu] = get_(S_clu, 'rho', 'viClu');
CORE_QUANTILE = .99;
knn_merge_thresh = get_set_(P, 'knn_merge_thresh', 1);
[cviSpk_clu, nClu] = vi2cell_(viClu);
cviSpk_core_clu = cell(nClu,1);
for iClu = 1:nClu    
    viSpk1 = int64(cviSpk_clu{iClu});
    vrRho1 = vrRho_spk(viSpk1);
    cviSpk_core_clu{iClu} = viSpk1(vrRho1 >= quantile_vr_(vrRho1, CORE_QUANTILE));    
end
if isempty(miKnn)
    cmiKnn_core_clu = load_miKnn_spk_(P, viSite_spk, cviSpk_core_clu);
end
mnKnn_clu = zeros(nClu);
for iClu1 = 1:nClu
    if ~isempty(miKnn)
        viSpk1 = miKnn(:,cviSpk_core_clu{iClu1});
    else
        viSpk1 = cmiKnn_core_clu{iClu1};
    end
    viSpk1 = viSpk1(:);
    for iClu2 = 1:nClu
        if iClu1 == iClu2, continue; end
        if ~isempty(miKnn)
            viSpk2 = miKnn(:,cviSpk_core_clu{iClu2});
        else
            viSpk2 = cmiKnn_core_clu{iClu2};
        end
        mnKnn_clu(iClu2, iClu1) = sum(ismember(viSpk2(:), viSpk1));
    end
end   
[viMap, viUniq_] = ml2map_(mnKnn_clu >= knn_merge_thresh);
nClu_post = numel(viUniq_);
viMap = viMap(:);
S_clu.viClu = map_index_(viMap, S_clu.viClu, 0);
S_clu = S_clu_prune_icl_(S_clu);

fprintf('S_clu_peak_merge_: %d->%d cluster centers (knn_merge_thresh=%d, took %0.1fs)\n', ...
    nClu, nClu_post, knn_merge_thresh, toc(t1));
end %func


%--------------------------------------------------------------------------
% 9/3/2019 JJJ: run_mode can be an array (sequential merge)
% 9/17/2018 JJJ: merge peaks based on their waveforms
function [viMap, S_clu, nClu_post] = S_clu_peak_merge1_(S_clu, P, viSite_spk) 

NUM_KNN_MERGE = 4;
miKnn = get_(S_clu, 'miKnn'); % to load from disk
knn_merge_thresh = get_set_(P, 'knn_merge_thresh', 1);
nClu = numel(setdiff(unique(S_clu.viClu), 0));
if ~isempty(miKnn)
    miKnn_clu = miKnn(:,S_clu.icl);
else
    miKnn_clu = load_miKnn_spk_(P, viSite_spk, S_clu.icl); 
    cmiKnn_clu = arrayfun_(@(x)miKnn_clu(:,x), (1:size(miKnn_clu,2))');
    cmiKnn_core_clu = load_miKnn_spk_(P, viSite_spk, cmiKnn_clu);
end

% find exact knn of the peaks using feature matrix
[mnKnn_lower_clu, mnKnn_upper_clu] = deal(zeros(nClu));
for iClu1 = 1:nClu
    viSpk1 = miKnn_clu(:,iClu1);     
    if ~isempty(miKnn)
        mi_ = miKnn(:,viSpk1);
    else
        mi_ = cmiKnn_core_clu{iClu1};
    end
    viSpk1 = sort(mi_(1:NUM_KNN_MERGE,:));
    if iClu1 > 1
        mnKnn_lower_clu(1:iClu1-1,iClu1) = sum(ismember(miKnn_clu(:,1:iClu1-1), viSpk1))';
    end
    if iClu1 < nClu
        mnKnn_upper_clu(iClu1+1:end,iClu1) = sum(ismember(miKnn_clu(:,iClu1+1:end), viSpk1))';
    end
end   
mnKnn_lower_clu = mnKnn_lower_clu + mnKnn_lower_clu';
mnKnn_upper_clu = mnKnn_upper_clu + mnKnn_upper_clu';
mnKnn_clu = min(mnKnn_lower_clu, mnKnn_upper_clu);

[viMap, viUniq_] = ml2map_(mnKnn_clu >= knn_merge_thresh);
nClu_post = numel(viUniq_);
viMap = viMap(:);
S_clu.viClu = map_index_(viMap, S_clu.viClu, 0);
S_clu = S_clu_prune_icl_(S_clu);

fprintf('S_clu_peak_merge_: %d->%d cluster centers (knn_merge_thresh=%d)\n', ...
    nClu, nClu_post, knn_merge_thresh);
end %func


%--------------------------------------------------------------------------
% Call from irc.m
function compile_cuda_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
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
function edit_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end

function out1 = load_prb_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = read_meta_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = isTextFile_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_batch_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = list_files_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = meta2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_merge_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ircpath_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = read_cfg_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = loadjson_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = filesize_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2ref_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = car_reject_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_copy_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cast_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2tr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsample_vr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cell2mat_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_default_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_filter_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gt2mda_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_dir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = meanSubt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = templateMatch_post_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = memory_matlab_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = recording_duration_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = squeeze_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refresh_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = find_site_spk23_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mn2tn_wav_spk2_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = matchFileExt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
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
function out1 = S_clu_quality_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = map_index_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_prune_icl_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

function [out1, out2] = ml2map_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = assignCluster_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = dpclus_remove_count_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = readmda_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = mr2thresh_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = filt_car_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = findNearSites_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = shift_range_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = wav_car_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = get_fig_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end

% function [out1, out2, out3] = S_clu_peak_merge_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = fopen_mda_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = fopen_nsx_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = plan_load_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = detect_spikes_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = cuda_delta_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = cuda_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = unique_count_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = rho_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = delta_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end


