% works on mda format only for now
% no manual interface now
% analyzes single mda file format
% todo: multiple mda files
% todo: multiple formats

%--------------------------------------------------------------------------
function varargout = irc2(vcDir_in, vcDir_out, vcFile_arg, vcArg3)
% irc2(vcDir_in, vcDir_out, vcFile_arg)
% irc2(vcCmd, vcArg1, vcArg2)
fDebug = 0;

if nargin<1, vcDir_in = ''; end
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end
if nargin<4, vcArg3 = ''; end
 
persistent vcFile_prm_

% batch processing. it uses default param for now
fDetect = read_cfg_('fForceRerun');

% run in batch mode
if iscell(vcDir_in) && iscell(vcDir_out)    
    [csDir_in, csDir_out, csFile_arg, fParfor] = deal(vcDir_in, vcDir_out, vcFile_arg, vcArg3);
    if isempty(fParfor), fParfor = false; end
    parfor iFile = 1:numel(csDir_in)
        try
            if exist_file_(fullfile(csDir_out{iFile}, 'firings.mda')) && ~fDetect
                continue;
            end
            fprintf('irc2 batch-processing %s (%d/%d)\n', csDir_in{iFile}, iFile, numel(csDir_in));
            irc2(csDir_in{iFile}, csDir_out{iFile}, csFile_arg{iFile}, fParfor);
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
vcFile_prm = dir2prm_(vcArg1);
if isempty(vcFile_prm)
    vcFile_prm = vcFile_prm_;
else
    vcFile_prm_ = vcFile_prm;
end
switch lower(vcCmd)
    case {'run-hs', 'run-herdingspikes', 'run-herdingspikes2'}
        run_spikeforest2_('herdingspikes2', vcArg1, vcArg2, vcArg3); return;
    case {'run-sc', 'run-spykingcircus'}
        run_spikeforest2_('spykingcircus', vcArg1, vcArg2, vcArg3); return;
    case {'run-ms4', 'run-mountainsort4', 'run-mountainsort'}
        run_spikeforest2_('mountainsort4', vcArg1, vcArg2, vcArg3); return;
    case {'run-tdc', 'run-tridesclous'}, run_spikeforest2_('tridesclous', vcArg1, vcArg2, vcArg3); return;
    case {'run-klusta'}, run_spikeforest2_('klusta', vcArg1, vcArg2, vcArg3); return;
    case {'run-jrclust', 'run-jrc'}, run_spikeforest2_('jrclust', vcArg1, vcArg2, vcArg3); return;
    case {'run-ks', 'run-kilosort'}, run_spikeforest2_('kilosort', vcArg1, vcArg2, vcArg3); return;
        
    case {'run-ks2', 'run-kilosort2'}, run_ksort2(vcArg1, vcArg2, vcArg3); return;
        
    case 'export-mda', save_firings_mda_(vcFile_prm, vcArg2); return;
    case {'which', 'select'}, fprintf('%s\n', vcFile_prm); return;
    case {'export-sf2', 'export-spikeforest2', 'export2spikeforest2'}, export_sf2_(); return;
    case {'export-ws', 'export-workspace', 'export'}, export_workspace_(vcFile_prm); return;
    case 'compile-deploy', compile_cuda_(vcArg1, '1'); return
    case {'compile', 'install'}, compile_cuda_(vcArg1, '0'); return
    case 'readmda_header', varargout{1} = readmda_header_(vcArg1); return;
    case 'mcc', irc('mcc'); return; 
    case {'join-mda', 'join_mda', 'joinmda'}
        join_mda_(vcArg1, vcArg2); return;
    case {'readmda', 'read-mda', 'read_mda'}
        mda = readmda_(vcArg1); assignWorkspace_(mda); return;
    case 'import-clip'
        [S0, P] = import_clip_(vcArg1); 
    case 'edit', edit_(vcFile_prm); return;
    case 'juxta'
        convert_mda_ui('english'); return;
    case 'version'
        if nargout==0, version_(); 
        else, varargout{1} = version_(); 
        end
        return;
    case 'scoreboard', irc2_scoreboard(); return;
    case {'spikesort', 'detectsort', 'detectsort-verify', 'detect-sort', ...
            'sort', 'auto', '', 'describe', 'manual', ...
            'validate', 'verify', 'auto-verify', 'sort-verify', 'spikesort-verify'}

        fprintf('irc2 (%s) opening %s\n', version_(), vcFile_prm);
        if isempty(vcFile_prm), fprintf(2, 'provide .prm file.\n'); return; end
        if ~exist_file_(vcFile_prm)
            fprintf(2, 'File does not exist: %s\n', vcFile_prm); return;
        end
        P = file2struct_(vcFile_prm);        
        switch lower(vcCmd)
            case {'detect-sort', 'spikesort', 'spikesort-verify'}, clear_(); fDetect = 1; fSort = 1;
            case {'sort', 'sort-verify'}, clear_('sort'); fDetect = 0; fSort = 1;   
            case {'auto', 'auto-verify'}, clear_('sort'); fDetect = 0; fSort = 0;
            case 'describe', describe_(vcFile_prm); return;
            case {'verify', 'validate'}, validate_(P, fPlot_gt); return;
            case {'manual', 'ui'}, irc2_ui(P); return;
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
            'test-monotrode1', 'test-monotrode2', 'test-monotrode3', ...
            'test-boyden', 'test-boyden2'}
        vcDir_in = get_test_data_(strsplit_get_(vcCmd,'-',2));
        [fDetect, fSort, fValidate] = deal(1, 1, 1);
    case 'test-all', test_all_(); return;
%     case 'export', irc('export', vcArg1); return;
    case {'export-phy', 'phy'}, irc2phy(vcFile_prm, vcArg2); return;
    case {'export-klusters', 'klusters', 'neurosuite'}, irc2klusters_v2(vcArg1, vcArg2); return;
    case 'call'
        switch nargout
            case 0, call_(vcArg1, vcArg2);
            case 1, varargout{1} = call_(vcArg1, vcArg2);
            case 2, [varargout{1}, varargout{2}] = call_(vcArg1, vcArg2);
            case 3, [varargout{1}, varargout{2}, varargout{3}] = call_(vcArg1, vcArg2);
            case 4, [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = call_(vcArg1, vcArg2);
        end %switch
        return;
    otherwise % directory running mode
        vcCmd=''; clear_(); 
        fValidate = exist_file_(fullfile(vcDir_in, 'firings_true.mda'));
end

fprintf('Running irc2.m (%s)\n', version_());
if isempty(P)
    fParfor = vcArg3;
    P = makeParam_(vcDir_in, vcDir_out, vcFile_arg, fParfor);    
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
    save_clu_(S0.S_clu, P);
end
try
    S0.S_auto = auto_(S0, P);
    save_auto_(S0.S_auto, P);
catch E
    fprintf(2, 'auto-merging error: \n\t%s\n', P.vcFile_prm);
    disp(lasterr);
    return;
end

% output
describe_(S0);
vcFile_firings_mda = fullfile(P.vcDir_out, 'firings.mda');
save_firings_mda_(S0, vcFile_firings_mda);

% Validate
if fValidate, validate_(P, fPlot_gt); end
end %func


%--------------------------------------------------------------------------
% 11/6/18 JJJ: Displaying the version number of the program and what's used. #Tested
function [vcVer, vcDate, vcHash] = version_()
vcVer = 'v5.5.3';
vcDate = '01/14/2020';
vcHash = file2hash_();

if nargout==0
    fprintf('%s (%s) installed, MD5: %s\n', vcVer, vcDate, vcHash);
    return;
end
end %func


%--------------------------------------------------------------------------
function A = readmda_(fname)
% Author: Jeremy Magland, modified by JJJ
% Jan 2015; Last revision: 15-Feb-2106
if ~exist_file_(fname)
    error('File does not exist: %s', fname);
end
F=fopen(fname,'r','l');

% read the first header sample: data type
try
    code=fread(F,1,'int32');
catch
    error('Problem reading file: %s',fname);
end

% read the second header sample: number of dimensions
if (code>0) 
    num_dims=code;
    code=-1;
else
    fread(F,1,'int32');
    num_dims=fread(F,1,'int32');    
end

% read the length per dimension
dim_type_str='int32';
if (num_dims<0)
    num_dims=-num_dims;
    dim_type_str='int64';
end

% read length per dimension
S = fread(F, num_dims, dim_type_str)';
N = prod(S);

switch code
    case -1
        A = fread(F,N*2,'*float');
        A = A(1:2:end) + sqrt(-1) * A(2:2:end);
    case -2, A = fread(F,N,'*uchar');
    case -3, A = fread(F,N,'*float');
    case -4, A = fread(F,N,'*int16');
    case -5, A = fread(F,N,'*int32');
    case -6, A = fread(F,N,'*uint16');
    case -7, A = fread(F,N,'*double');
    case -8, A = fread(F,N,'*uint32');
    otherwise, error('Unsupported data type code: %d',code);
end
A = reshape(A, S);
fclose(F);
end %func


%--------------------------------------------------------------------------
function test_all_()
irc2('test-monotrode');
irc2('test-tetrode');
irc2('test-static');
irc2('test-drift');
end %func


%--------------------------------------------------------------------------
% convert directory to prm file if directory path is given
function vcFile_prm = dir2prm_(vcDir_in)
[vcDir1, vcFile1, vcExt1] = fileparts(vcDir_in);
if ~strcmpi(vcExt1, '.prm')
    vcDir_out = fullfile(vcDir_in, 'irc2');
    S_prm = dir(fullfile(vcDir_out, '*.prm'));
    if numel(S_prm) == 1
        vcFile_prm = fullfile(S_prm.folder, S_prm.name);
    else
        vcFile_prm = '';
    end
else
    vcFile_prm = vcDir_in;
end
end %func


%--------------------------------------------------------------------------
function save_clu_(S_clu, P)
if isempty(S_clu), return; end

% save separately
% rho: nSpk x 1: single
% delta: nSpk x 1: single
% nneigh: nSpk x 1: int64
% ordrho: nSpk x 1: double
vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');
csVar = {'rho', 'delta', 'nneigh', 'ordrho'};
[S_clu, S_clu.S_var] = struct_save_bin_(S_clu, [vcFile_prm_, '_clu.irc'], csVar);
struct_save_(S_clu, [vcFile_prm_, '_clu_irc.mat'], 1);
end %func


%--------------------------------------------------------------------------
function save_auto_(S_auto, P)
if isempty(S_auto), return; end

% save separately
% viClu: nSpk x 1: int32
% cviSpk_clu: nClu x 1: contains int64

csVar = {'viClu', 'cviSpk_clu'};

vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');
[S_auto, S_auto.S_var] = struct_save_bin_(S_auto, [vcFile_prm_, '_auto.irc'], csVar);
struct_save_(S_auto, [vcFile_prm_, '_auto_irc.mat'], 1);
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
function export_workspace_(vcFile_prm)
if ~exist_file_(vcFile_prm)
    error('%s does not exist\n', vcFile_prm);
end
S0 = load0_(vcFile_prm);
S0.trPc_spk = load_fet_(S0, S0.P, 1);
S0.trPc2_spk = load_fet_(S0, S0.P, 2);
assignWorkspace_(S0);
end %func


%--------------------------------------------------------------------------
% used externally, do not remove
function trPc_spk = load_fet_(S0, P, iFet)
if nargin<3, iFet = 1; end

% return empty if trPc_spk should be loaded partially
% if get_set_(P, 'fLargeRecording', 0)
%     trPc_spk=[]; 
%     return; 
% end

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
        trPc_spk = load_bin_merge_(csFiles_in, S0.type_fet, S0.dimm_fet, S0.ccviSpk_site_load);
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
        trPc_spk = load_bin_merge_(csFiles_in, S0.type_fet, S0.dimm_fet, S0.ccviSpk_site2_load);
    end
else
    error('invalid feature id');
end
fprintf(' took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function trPc_spk = load_bin_merge_(csFiles_in, type_fet, dimm_fet, ccviSpk_site_load)
if nargin<4, ccviSpk_site_load={}; end
if ~exist_file_(csFiles_in{1})
    trPc_spk = []; return; % memory out
end
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
    if isempty(tr_), trPc_spk = []; return; end
    if ~isempty(ccviSpk_site_load)
        cviSpk_site1 = ccviSpk_site_load{iFile};
        vi1 = cat(1,cviSpk_site1{:}) + iOffset;
    else
        vi1 = (1:dimm_fet1(end)) + iOffset;
    end
    trPc_spk(:,:,vi1) = tr_;
    iOffset = vi1(end);
end %for
end %func


%--------------------------------------------------------------------------
function flag = is_sorted_(P)
% return true if already detected. .spkwav file must exist

vcFile_clu_mat = strrep(P.vcFile_prm, '.prm', '_clu_irc.mat');
flag = exist_file_(vcFile_clu_mat);
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

% save separately
% viTime_spk: nSpk x 1: int64
% viSite_spk, viSite2_spk: nSpk x 1: int32
% vrAmp_spk: nSpk x 1: single
% vrPow_spk: nSpk x 1: single
% mrPos_spk: nSpk x 2: single

csVar_spk = {'viTime_spk', 'viSite_spk', 'viSite2_spk', 'vrAmp_spk', ...
    'vrPow_spk', 'mrPos_spk', 'ccviSpk_site_load', 'ccviSpk_site2_load'};
[S0, S0.S_var] = struct_save_bin_(S0, [vcFile_prm_, '_spk.irc'], csVar_spk);

trPc_spk = gather_(get_(S0, 'trPc_spk'));
S0.trPc_spk = [];
fSave_fet = get_set_(S0.P, 'fSave_fet', 0);
if ~isempty(trPc_spk) && fSave_fet
    S0.trPc_spk = [];
    [S0.dimm_fet, S0.type_fet] = write_bin_([vcFile_prm_, '_fet.irc'], trPc_spk);
end

trPc2_spk = gather_(get_(S0, 'trPc2_spk'));
S0.trPc2_spk = [];
if ~isempty(trPc2_spk) && fSave_fet
    S0.trPc2_spk = [];
    write_bin_([vcFile_prm_, '_fet2.irc'], trPc2_spk);
end

S_clu = get_(S0, 'S_clu');  S0.S_clu = [];
if ~isempty(S_clu), save_clu_(S_clu, S0.P); end

S_auto = get_(S0, 'S_auto');  S0.S_auto = [];
if ~isempty(S_auto), save_auto_(S_auto, S0.P); end

struct_save_(S0, [vcFile_prm_, '_irc.mat'], 1);
end


%--------------------------------------------------------------------------
% 2020/1/13: can save cell of arrays to a binary file
function [S0, S_var] = struct_save_bin_(S0, vcFile, csName_var)
% save variables and clear field
S_var = struct('vcFile', vcFile);
S_var.csName_var = csName_var;
S_var.csType_var = cell(size(csName_var));
S_var.cDimm_var = cell(size(csName_var));

fid_w = fopen(vcFile, 'w');
for iVar = 1:numel(csName_var)
    name_ = csName_var{iVar};
    val_ = S0.(name_);
    S0.(name_) = []; % clear variable
    if ~iscell(val_)
        [S_var.cDimm_var{iVar}, S_var.csType_var{iVar}] = write_bin_(fid_w, val_);
    else
        [S_var.cDimm_var{iVar}, S_var.csType_var{iVar}] = save_cell_(fid_w, val_);
    end
end
fclose(fid_w);
end %func


%--------------------------------------------------------------------------
function [cDimm_cell1, csType_cell1] = save_cell_(fid_w, cVal)
if ~iscell(cVal)
    [cDimm_cell1, csType_cell1] = write_bin_(fid_w, cVal);
else
    [cDimm_cell1, csType_cell1] = deal(cell(size(cVal)));        
    for iCell1 = 1:numel(cVal)  
        val_ = cVal{iCell1};
        if ~iscell(val_)
            [cDimm_cell1{iCell1}, csType_cell1{iCell1}] = write_bin_(fid_w, val_);
        else
            [cDimm_cell1{iCell1}, csType_cell1{iCell1}] = save_cell_(fid_w, val_);
        end
    end 
end
end %func


%--------------------------------------------------------------------------
function S = struct_load_bin_(S_var, S)
% save variables and clear field
if nargin<2, S = struct(); end
import_struct_(S_var); % import all fields in this struct
fid_r = fopen(vcFile, 'r');
for iVar = 1:numel(csName_var)
    [type_, dimm_] = deal(csType_var{iVar}, cDimm_var{iVar});
    if iscell(type_)
        S.(csName_var{iVar}) = load_cell_(fid_r, type_, dimm_);
    else
        S.(csName_var{iVar}) = load_bin_(fid_r, type_, dimm_);
    end
end
fclose(fid_r);
end %func


%--------------------------------------------------------------------------
function cVal = load_cell_(fid_r, cType, cDimm)
if ~iscell(cType)
    cVal = load_bin_(fid_r, cType, cDimm);
else
    cVal = cell(cType);
    for iCell1 = 1:numel(cVal)  
        [type_, dimm_] = deal(cType{iCell1}, cDimm{iCell1});
        if ~iscell(type_)
            cVal{iCell1} = load_bin_(fid_r, type_, dimm_);
        else
            cVal{iCell1} = load_cell_(fid_r, type_, dimm_);
        end
    end 
end
end %func


%--------------------------------------------------------------------------
function S0 = load0_(vcFile_prm, fLoad_bin)

if nargin<2, fLoad_bin = 1; end

fprintf('Loading %s... ', vcFile_prm); t_fun = tic;
% set(0, 'UserData', []);
vcFile_mat = strrep(vcFile_prm, '.prm', '_irc.mat');
vcFile_clu_mat = strrep(vcFile_prm, '.prm', '_clu_irc.mat');
vcFile_auto_mat = strrep(vcFile_prm, '.prm', '_auto_irc.mat');
vcFile_knn = strrep(vcFile_prm, '.prm', '_knn.irc');

S0 = load(vcFile_mat);
if isfield(S0, 'S_var') && fLoad_bin
    S0 = struct_load_bin_(S0.S_var, S0);
end

S0.S_clu = get_(S0, 'S_clu');
if isempty(S0.S_clu)
    if exist_file_(vcFile_clu_mat)
        S0.S_clu = load(vcFile_clu_mat);
    end
    if isfield(S0.S_clu, 'S_var') && fLoad_bin
        S0.S_clu = struct_load_bin_(S0.S_clu.S_var, S0.S_clu);
    end
end

S0.S_auto = get_(S0, 'S_auto');
if isempty(S0.S_auto)
    if exist_file_(vcFile_auto_mat)
        S0.S_auto = load(vcFile_auto_mat);
    end
    if isfield(S0.S_auto, 'S_var') && fLoad_bin
        S0.S_auto = struct_load_bin_(S0.S_auto.S_var, S0.S_auto);
    end
end

if isempty(get_(S0.S_clu, 'miKnn'))
    if exist_file_(vcFile_knn)
        S0.S_clu.miKnn = load_bin_(vcFile_knn, S0.S_clu.type_knn, S0.S_clu.dimm_knn);
    end
end
% set(0, 'UserData', S0);
fprintf('took %0.1fs\n', toc(t_fun));
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

csFields = {'runtime_detect', 'memory_detect', 'memory_init', 'dimm_fet', 'P'};
csFields_clu = {'runtime_sort', 'memory_sort', 'nFeatures'};
csFields_auto = {'runtime_automerge', 'memory_auto', 'nSpk_unique', 'nClu'};

if ischar(S0)
    vcFile_prm = S0;
    vcFile_mat = strrep(vcFile_prm, '.prm', '_irc.mat');
    vcFile_clu_mat = strrep(vcFile_prm, '.prm', '_clu_irc.mat');
    vcFile_auto_mat = strrep(vcFile_prm, '.prm', '_auto_irc.mat');
    if exist_file_(vcFile_mat)
        S0 = load(vcFile_mat, csFields{:});
    end
    S0.S_clu = load_(vcFile_clu_mat, csFields_clu);
    S0.S_auto = load_(vcFile_auto_mat, csFields_auto);    
end
P=S0.P;

try
    [runtime_detect, runtime_sort, runtime_automerge] = ...
        deal(S0.runtime_detect, S0.S_clu.runtime_sort, S0.S_auto.runtime_automerge);
    runtime_total = runtime_detect + runtime_sort + runtime_automerge;
catch
    runtime_total = nan;
end
try
    tDur = recording_duration_(S0.P, S0); 
catch
    tDur = nan;
end
memory_detect = S0.memory_detect - S0.memory_init;
memory_sort = S0.S_clu.memory_sort - S0.memory_init;
memory_auto = S0.S_auto.memory_auto - S0.memory_init;
nSites = numel(P.viSite2Chan);
viShank_site = get_(P, 'viShank_site');
if isempty(viShank_site)
    nShanks = 1;
else
    nShanks = numel(unique(viShank_site));
end
nSites_spk = size(P.miSites,1);
nSpk = S0.dimm_fet(3);
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
    csDesc{end+1} = sprintf('    Whiten:                 %d', get_set_(P, 'fWhiten', 0));
    csDesc{end+1} = sprintf('    FFT threshold:          %d', get_set_(P, 'fft_thresh', 0));
    csDesc{end+1} = sprintf('    blank threshold:        %d', get_set_(P, 'blank_thresh', 0));    
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
end
try
   
    S_auto = get_(S0, 'S_auto');
    csDesc{end+1} = sprintf('Cluster');       
    csDesc{end+1} = sprintf('    #Clusters:              %d', get_(S_auto, 'nClu'));
    csDesc{end+1} = sprintf('    #Unique events:         %d', get_(S_auto, 'nSpk_unique'));
    csDesc{end+1} = sprintf('    min. spk/clu:           %d', P.min_count);
    csDesc{end+1} = sprintf('    Cluster method:         %s', P.vcCluster);
    csDesc{end+1} = sprintf('    knn:                    %d', P.knn);
    csDesc{end+1} = sprintf('    step_sec_drift:         %0.1fs', P.step_sec_drift);
    csDesc{end+1} = sprintf('    batch_sec_drift:        %0.1fs', P.batch_sec_drift);
%     csDesc{end+1} = sprintf('    nTime_drift:            %d', P.nTime_drift);
%     csDesc{end+1} = sprintf('    nTime_batch:            %d', P.nTime_batch);
    csDesc{end+1} = sprintf('Auto-merge');   
    csDesc{end+1} = sprintf('    delta_cut:              %0.3f', get_set_(P, 'delta_cut', 1));
    csDesc{end+1} = sprintf('    maxWavCor:              %0.3f', P.maxWavCor);
catch
end
try
    csDesc{end+1} = sprintf('Runtime (s)');
    csDesc{end+1} = sprintf('    Detect + feature (s):   %0.1fs', runtime_detect);    
    csDesc{end+1} = sprintf('    Cluster runtime (s):    %0.1fs', runtime_sort);
    csDesc{end+1} = sprintf('    merge runtime (s):      %0.1fs', runtime_automerge);
    csDesc{end+1} = sprintf('    Total runtime (s):      %0.1fs', runtime_total);
    csDesc{end+1} = sprintf('    Runtime speed:          x%0.1f realtime', tDur / runtime_total);    
    csDesc{end+1} = sprintf('    Processing speed:       %0.1f spikes/s', nSpk / runtime_total);    

    csDesc{end+1} = sprintf('memory usage (GiB):         %0.3f', max(memory_detect, memory_sort, memory_auto)/2^30);
    csDesc{end+1} = sprintf('    detect:                 %0.3f', memory_detect/2^30);
    csDesc{end+1} = sprintf('    sort:                   %0.3f', memory_sort/2^30);
    csDesc{end+1} = sprintf('    auto-merge:             %0.3f', memory_auto/2^30);

    csDesc{end+1} = sprintf('Execution');
    csDesc{end+1} = sprintf('    irc2 version:           %s', get_(P, 'vcVersion'));
    csDesc{end+1} = sprintf('    fGpu (GPU use):         %d', P.fGpu);
    csDesc{end+1} = sprintf('    fParfor (parfor use):   %d', P.fParfor);
    csDesc{end+1} = sprintf('    fLargeRecording:        %d', isLargeRecording_(P)); 
    csDesc{end+1} = sprintf('    Parameter file:         %s', P.vcFile_prm);
catch
end
if nargout==0
    cellfun(@(x)disp(x), csDesc);
end
end %func


%--------------------------------------------------------------------------
function S = load_(vcFile, csFields)
if nargin<2, csFields={}; end
S = [];
if exist_file_(vcFile)
    try
        if ~isempty(csFields)
            S = load(vcFile, csFields{:});
        else
            S = load(vcFile);
        end
    catch
    end
end
end %func


%--------------------------------------------------------------------------
function flag = isLargeRecording_(P)
flag = get_(P, 'nTime_drift') > get_set_(P, 'nTime_max_drift', 2^8);
end %func


%--------------------------------------------------------------------------
function save_firings_mda_(S0, vcFile_firings_mda)
% save_firings_mda_(S0, vcFile_firings_mda)
% save_firings_mda_(vcFile_prm, vcFile_firings_mda)
if nargin<2, vcFile_firings_mda = ''; end

try
    if ischar(S0)
        vcFile_prm = S0;
        S0 = load(strrep(vcFile_prm, '.prm', '_irc.mat'), 'S_var'); 
        S0 = struct_load_bin_(S0.S_var);
        S0.S_clu = load(strrep(vcFile_prm, '.prm', '_clu_irc.mat'), 'viClu'); 
    else
        vcFile_prm = S0.P.vcFile_prm;
    end
    if isempty(vcFile_firings_mda)
        vcFile_firings_mda = fullfile(fileparts(vcFile_prm), 'firings.mda');
    end
    mr = [S0.viSite_spk(:), S0.viTime_spk(:), S0.S_auto.viClu(:)];
    S0 = []; %free memory
catch
    error('save_firings_mda_: invalid format');
end
writemda_(vcFile_firings_mda, double(mr'));
fprintf('Wrote to %s\n', vcFile_firings_mda);
end %func


%--------------------------------------------------------------------------
% 64-bit addressing compatible
function writemda_(vcFile, var)
writemda_fid(vcFile, var);
end %func
            

%--------------------------------------------------------------------------
% auto merge
function S_auto = auto_(S0, P)

fDebug = 0;

fprintf('\nauto-merging...\n'); runtime_automerge = tic;
if fDebug, fprintf(2, 'auto_() is running in Debugging mode\n'); end

% refresh clu, start with fundamentals
S_auto = postCluster_(S0.S_clu, P, S0.viSite_spk); % peak merging

maxWavCor = get_set_(P, 'maxWavCor', .99);
nClu = get_(S_auto, 'nClu');
if maxWavCor<1 && nClu > 1
    try
        fprintf('\tMerging templates...\n\t'); t_merge=tic;
        nClu_pre = S_auto.nClu;
        S0.S_auto = S_auto;
        [mlDist_clu, viClu_delete] = wave_similarity_(S0, P);
        [S_auto, nClu_post] = ml_merge_clu_(S_auto, mlDist_clu, viClu_delete);    
        fprintf('\tMerged waveforms (%d->%d->%d), took %0.1fs\n', ...
            nClu_pre, nClu_post+numel(viClu_delete), nClu_post, toc(t_merge));
    catch E
        fprintf(2, 'Merging templates failed\n');
        disp(E);
%         disp(E.stack);
%         disp(E.message);
        assignWorkspace_(E);
        save0_(S0);
        rethrow(E);
    end
end

S_auto = S_clu_refrac_(S_auto, P, [], S0.viTime_spk); % refractory violation removal
S_auto = S_auto_refresh_(S_auto, 1, S0.viSite_spk);
S_auto = S_clu_sort_(S_auto, 'viSite_clu');
S_auto.runtime_automerge = toc(runtime_automerge);
S_auto.memory_auto = memory_matlab_();
fprintf('\tauto-merging took %0.1fs (fGpu=%d, fParfor=%d)\n', ...
    S_auto.runtime_automerge, P.fGpu, P.fParfor);
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Output message is added
% 8/2/17 JJJ: Test and documentation
function vcMsg = assignWorkspace_(varargin)
% Assign variables to the Workspace
vcMsg = {};
for i=1:numel(varargin)
    if ~isempty(varargin{i})
        assignin('base', inputname(i), varargin{i});
        vcMsg{end+1} = sprintf('assigned ''%s'' to workspace\n', inputname(i));        
    end
end
vcMsg = cell2mat(vcMsg);
if nargout==0, fprintf(vcMsg); end
end %func


%--------------------------------------------------------------------------
% 4/12/2019 JJJ: Template merging cluster
% no membership reassignment, no core calculation
function [S_auto, nClu_post] = ml_merge_clu_(S_auto, mlWavCor_clu, viClu_delete)
if nargin<3, viClu_delete=[]; end
if ~isempty(viClu_delete)
    mlWavCor_clu(viClu_delete,:) = false;
    mlWavCor_clu(:,viClu_delete) = false;
    S_auto.viClu(ismember(S_auto.viClu, viClu_delete)) = 0;
    fprintf('\tdeleted %d clusters with SNR below the detection threshold\n', numel(viClu_delete));
end
% nClu_pre = size(mlWavCor_clu,1);
[S_auto.viClu, viMap_new] = map_using_ml_(S_auto.viClu, mlWavCor_clu);
nClu_post = sum(viMap_new>0);
end %func


%--------------------------------------------------------------------------
function [viClu, viMap_new] = map_using_ml_(viClu, ml)
[viMap_clu, viMap_new] = ml2map_(ml);
vlPos = viClu > 0;
viClu(vlPos) = viMap_clu(viClu(vlPos)); %translate cluster number
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
function [mlDist_clu, viClu_remove] = wave_similarity_(S0, P)
S_auto = get_(S0, 'S_auto');
S_clu = get_(S0, 'S_clu');
t_fun = tic;
% fprintf('\tAutomated merging based on waveform similarity...\n'); t_template=tic;
viClu = S_auto.viClu;
vrRho = S_clu.rho;
[ccviSpk_site_load, ccviSpk_site2_load, type_fet, dimm_fet, mrPv, vrThresh_site] = ...
    get_(S0, 'ccviSpk_site_load', 'ccviSpk_site2_load', 'type_fet', ...
        'dimm_fet', 'mrPv_global', 'vrThresh_site');
nShift_max = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
viShift = -nShift_max:nShift_max;
[knn, nSpk_min, vcFile_prm, maxWavCor] = get_(P, 'knn', 'knn', 'vcFile_prm', 'maxWavCor');
nClu = S_auto.nClu;

% try loading the entire miKnn to RAM
nSites = size(P.miSites,2);
[viLim_drift, mlDrift] = get_(S_clu.S_drift, 'viLim_drift', 'mlDrift');
% nDrift = size(mlDrift,1);
S_auto = makeStruct_(nClu, nSites, knn, vcFile_prm, nSpk_min, ...
    vrRho, viClu, viLim_drift, ccviSpk_site_load, ccviSpk_site2_load, ...
    type_fet, dimm_fet, mrPv, vrThresh_site, viShift, mlDrift, maxWavCor);

% compute pairwise distance in parallel by sites
[cmlDist_site, cvlExist_site] = deal(cell(nSites, 1));
fParfor = get_set_(P, 'fParfor', 1) && nSites > 1;
if fParfor %&& ~isLargeRecording_(P)
    try
        parfor iSite = 1:nSites
            try
                [cmlDist_site{iSite}, cvlExist_site{iSite}] = wave_similarity_site_(iSite, S_auto);
            catch
            end
        end
    catch
    end
end

% merge cluster pairwise distance
mlDist_clu = false(nClu);
mlDist_clu(sub2ind([nClu,nClu], 1:nClu, 1:nClu)) = true; % self join
vlExist_clu = false(1, nClu);
for iSite = 1:nSites
    if isempty(cmlDist_site{iSite})
        [cmlDist_site{iSite}, cvlExist_site{iSite}] = wave_similarity_site_(iSite, S_auto);
    end
    mlDist_clu = mlDist_clu | cmlDist_site{iSite};
    vlExist_clu = vlExist_clu | cvlExist_site{iSite};
    [cmlDist_site{iSite}, cvlExist_site{iSite}] = deal([]); % clear memory
end
viClu_remove = find(~vlExist_clu);
fprintf('\n\twave_similarity_: took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
function [mlDist_clu, vlExist_clu] = wave_similarity_site_(iSite1, S_auto)

NUM_KNN = 10;
fUseSecondSite = 1;

import_struct_(S_auto);
nDrift = size(mlDrift, 1);
% fprintf('\twave_similarity_site_pre_: Site%d... ', iSite1); t_fun=tic;
miKnn1 = load_miKnn_site_(S_auto, iSite1);
miKnn1 = miKnn1(1:min(NUM_KNN, size(miKnn1,1)), :);
% [viLim_drift, nDrift, viClu, nClu, nSpk_min, vrThresh_site, mrPv, mlDrift] = ...
%     get_(S_auto, 'viLim_drift', 'nDrift', 'viClu', 'nClu', 'nSpk_min', 'vrThresh_site', 'mrPv', 'mlDrift');
thresh1 = vrThresh_site(iSite1);
[trPc1, viSpk1] = load_fet_site_(S_auto, 1, iSite1);
% cc1_drift_clu = cell(nDrift, nClu);
cvii1_drift = vi2cell_(discretize(viSpk1, viLim_drift), nDrift);
[vrRho1, viClu1] = deal(S_auto.vrRho(viSpk1), viClu(viSpk1));
vrRho_1 = copy_mask_(S_auto.vrRho, viSpk1);
[~, miiKnn1] = ismember(miKnn1, viSpk1);

if fUseSecondSite
    [trPc2, viSpk2] = load_fet_site_(S_auto, 2, iSite1);
else
    trPc2 = [];
end
fSecondSite = ~isempty(trPc2);
if fSecondSite   
    vrRho_2 = copy_mask_(S_auto.vrRho, viSpk2);
    [~, miiKnn2] = ismember(miKnn1, viSpk2);
end
% cc2_drift_clu = cell(nDrift, nClu);
cviClu_drift = cell(nDrift,1);
ctrPc_drift = cell(nDrift,1);

for iDrift = 1:nDrift    
    vii1 = cvii1_drift{iDrift};
    if isempty(vii1), continue; end
    [vrRho11, viClu11, miKnn11, miiKnn11] = ...
        deal(vrRho1(vii1), viClu1(vii1), miKnn1(:,vii1), miiKnn1(:,vii1));
    [cviiSpk_clu_, ~, viClu_uniq] = vi2cell_(viClu11, nClu);
    if fSecondSite, miiKnn21 = miiKnn2(:,vii1); end
    [viClu_drift1, cmrPc_drift1] = deal([], {}); 
    for iClu = viClu_uniq
        vii_ = cviiSpk_clu_{iClu};
        [miKnn11_, miiKnn11_] = deal(miKnn11(:,vii_), miiKnn11(:,vii_));
        vrRho11_T = vrRho11(vii_)';
        vii1_ = miiKnn11_(vrRho_1(miKnn11_) >= vrRho11_T);  
        mrPc1 = mean_conditional_(trPc1, vii1_, nSpk_min, mrPv, thresh1);
        if ~isempty(mrPc1)
            viClu_drift1(end+1) = iClu;
            cmrPc_drift1{end+1} = mrPc1;
        end
        if fSecondSite
            miiKnn21_ = miiKnn21(:,vii_);
            vii2_ = miiKnn21_(vrRho_2(miKnn11_) >= vrRho11_T);
            mrPc2 = mean_conditional_(trPc2, vii2_, nSpk_min, mrPv, thresh1);
            if ~isempty(mrPc2)
                viClu_drift1(end+1) = iClu;
                cmrPc_drift1{end+1} = mrPc2;
            end
        end         
    end
    cviClu_drift{iDrift} = viClu_drift1;
    ctrPc_drift{iDrift} = cat(3, cmrPc_drift1{:});
end
[trPc1, trPc2, miKnn1, miiKnn1, miiKnn2] = deal([]); % clear memory
vlExist_clu = false(1, nClu);
vlExist_clu([cviClu_drift{:}]) = true;
mlDist_clu = false(nClu);

norm_mr_ = @(mr)mr ./ sqrt(sum(mr.^2,1)); 
tr2mr_pv_norm_ = @(tr,mr)norm_mr_(reshape(mr*reshape(tr,size(tr,1),[]),[],size(tr,3))); 

% distance calculation
for iDrift = 1:nDrift
    viDrift1 = find(mlDrift(:,iDrift));
    viClu1 = cviClu_drift{iDrift};
    if isempty(viClu1), continue; end
    trPc_clu1 = cat(3, ctrPc_drift{iDrift});
    if isempty(trPc_clu1), continue; end
    viClu2 = [cviClu_drift{viDrift1}];
    trPc_clu2  = cat(3, ctrPc_drift{viDrift1});
    if isempty(trPc_clu2), continue; end
    mrWav_clu2 = tr2mr_pv_norm_(trPc_clu2, mrPv);
    for iiClu1 = 1:numel(viClu1)
        iClu1 = viClu1(iiClu1);
        mrWav11 = pc2wav_shift_(trPc_clu1(:,:,iiClu1), mrPv, viShift);
        for iiClu2 = 1:numel(viClu2)
            iClu2 = viClu2(iiClu2);         
            if iClu2 > iClu1 % symmetric
                dist12 = max(mrWav_clu2(:,iiClu2)' * mrWav11);
                mlDist_clu(iClu2, iClu1) = mlDist_clu(iClu2, iClu1) | (dist12>=maxWavCor);
            end
        end
    end
end

fprintf('.');
end %func


%--------------------------------------------------------------------------
function mrWav = pc2wav_shift_(trPc, mrPv, viShift)

if isempty(trPc), mrWav=[]; return; end
[nPc, nSites, nSpk, nT] = deal(size(trPc,1), size(trPc,2), size(trPc,3), size(mrPv,1));
trWav = reshape(mrPv*reshape(trPc, nPc,[]), [nT, nSites, nSpk]);
mrWav = reshape(shift_trWav_(trWav, viShift), nT*nSites, []);
mrWav = mrWav ./ sqrt(sum(mrWav.^2,1)); 
end %func


%--------------------------------------------------------------------------
function trB = shift_trWav_(trA, viShift)
% ctr = cell(numel(viShift), 1);
dimm_tr = size(trA);
if ismatrix(trA), dimm_tr = [dimm_tr, 1]; end
n = dimm_tr(1);
vi0 = 1:n;
nShift = numel(viShift);
qr = zeros([dimm_tr, nShift], 'like', trA);
for iShift = 1:nShift
    iShift_ = viShift(iShift);
    vi_ = min(max(vi0 + iShift_, 1), n);
    qr(:,:,:,iShift) = trA(vi_,:,:);
end %for
% tr = cat(3, ctr{:});
trB = reshape(qr, [dimm_tr(1), dimm_tr(2), dimm_tr(3)*nShift]);
end % func


%--------------------------------------------------------------------------
function mrPc1 = mean_conditional_(trPc1, vii1, nSpk_min, mrPv, thresh1) 
mrPc1 = [];
if numel(vii1) >= nSpk_min
    mrPc1 = mean(trPc1(:,:,vii1),3);
    if abs(min(mrPv * mrPc1(:,1))) < thresh1
        mrPc1 = [];
    end
end
end %func


%--------------------------------------------------------------------------
function vrB = copy_mask_(vrA, viA)
vrB = zeros(size(vrA), 'like', vrA);
vrB(viA) = vrA(viA);
end %func


%--------------------------------------------------------------------------
function vii2 = find_sorted_(vi1, vi2)
% find 1 in 2 and return index in 2
% assume vi2 is sorted
[vl1, vii1] = ismember(sort(vi1), vi2);
vii2 = vii1(vl1);
end %func


%--------------------------------------------------------------------------
function vii2 = find_sorted__(vi1, vi2)
% find 1 in 2 and return index in 2
% assume vi2 is sorted
vi1 = sort(vi1);
lim2 = [find(vi2>=vi1(1), 1, 'first'), find(vi2<=vi1(end), 1, 'last')];
[vl1, vii1] = ismember(vi1, vi2(lim2(1):lim2(2)));
vii2 = vii1(vl1);
end %func


%--------------------------------------------------------------------------
function [viClu_remove, ctrPc_clu, cviSite_clu, cviDrift_clu] = ...
    find_low_snr_clu_(S0, P, ctrPc_clu, cviSite_clu, cviDrift_clu)
if nargin<5, cviDrift_clu = {}; end
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
    if ~isempty(cviDrift_clu)
        viDrift_clu1 = cviDrift_clu{iClu};
        cviDrift_clu{iClu} = viDrift_clu1(vlKeep1);
    end
end %for
viClu_remove = find(cellfun(@isempty, ctrPc_clu));
% fprintf('merge: removed %d/%d templates below SNR=%0.3f\n', nRemoved, nTotal, P.qqFactor);    
end %func


%--------------------------------------------------------------------------
% long recording version
function [ctrPc_clu, cviSite_clu, viClu_remove] = merge_clu_pre_load_(cS_pre, S0, P)
cviSite_clu = cell(size(cS_pre));

[cviFet_clu, ccviSpk_clu] = deal(cell(size(cS_pre)));
for iClu = 1:numel(cS_pre)
    S_pre = cS_pre{iClu};
    cviSite_clu{iClu} = [S_pre.viSite1(:); S_pre.viSite2(:)];    
    n1 = numel(S_pre.viSite1);
    n2 = numel(S_pre.viSite2);
    cviFet_clu{iClu} = [ones(n1,1); repmat(2,n2,1)];
    if n2>0
        ccviSpk_clu{iClu} = [cS_pre{iClu}.cviSpk1, cS_pre{iClu}.cviSpk2];
    else
        ccviSpk_clu{iClu} = cS_pre{iClu}.cviSpk1;
    end
end
ctrPc_clu = mean_fet_site_load_(cviSite_clu, cviFet_clu, ccviSpk_clu, S0, P);
viClu_remove = find_low_snr_clu_(S0, P, ctrPc_clu, cviSite_clu);
end %func


%--------------------------------------------------------------------------
function ctrPc_clu = mean_fet_site_load_(cviSite_clu, cviFet_clu, ccviSpk_clu, S0, P)
fprintf('mean_fet_load_... '); t1=tic;
viSite_all = cell2mat_(cviSite_clu);
viFet_all = cell2mat_(cviFet_clu);
cviSpk_all = [ccviSpk_clu{:}];
ctrPc_all = cell(size(cviSpk_all));

% read the whole site and trim
S0.vcFile_prm = P.vcFile_prm;
[cviSpk_site, cviSpk2_site] = calc_cviSpk_site_(S0, P);
nSites = size(P.miSites, 2);
for iSite = 1:nSites
    vi_all1 = find(viSite_all == iSite & viFet_all == 1);
    trPc_ = load_fet_site_(S0, 1, iSite);    
    [~,cviSpk_all1] = cellfun_(@(x)ismember(x,cviSpk_site{iSite}), cviSpk_all(vi_all1));
    ctrPc_all(vi_all1) = cellfun_(@(x)mean(trPc_(:,:,x),3), cviSpk_all1);

    trPc_ = load_fet_site_(S0, 2, iSite); 
    if ~isempty(trPc_)
        vi_all2 = find(viSite_all == iSite & viFet_all == 2);
        [~,cviSpk_all2] = cellfun_(@(x)ismember(x,cviSpk2_site{iSite}), cviSpk_all(vi_all2));
        ctrPc_all(vi_all2) = cellfun_(@(x)mean(trPc_(:,:,x),3), cviSpk_all2);
    end
    trPc_ = [];
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
function vrDist_clu1 = wav_similarity_clu_post_(iClu1, S_post)
[cviSite_clu, ctrPc_clu, nClu, viShift, mrPv] = ...
    struct_get_(S_post, 'cviSite_clu', 'ctrPc_clu', 'nClu', 'viShift', 'mrPv');
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
function mrWav1 = normalize_tr_(trPc, mrPv)
[nPc, ~, nSpk] = size(trPc);
mrWav1 = reshape(mrPv * reshape(trPc, nPc, []), [], nSpk);
mrWav1 = bsxfun(@rdivide, mrWav1, sqrt(sum(mrWav1.^2)));
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
function [vr, vi] = cvr2vr_vi_(cvr)
vr = cell2mat_(cvr);
vn1 = cellfun(@(x)size(x,1), cvr);
vi = cell2mat_(arrayfun(@(x)repmat(x, vn1(x),1), 1:numel(cvr), 'UniformOutput', 0)');
end %func


%--------------------------------------------------------------------------
function vr = cell2mat_(cvr)
% create a matrix that is #vectors x # cells
% remove empty
vi = find(cellfun(@(x)~isempty(x), cvr));
vr = cell2mat(cvr(vi));
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

S_drift = calc_drift_(S0, P);
if true % disable parfor for large recordings
    if get_set_(P, 'fGpu', 1)
        P.fParfor = get_set_(P, 'fParfor', 1) && ~isLargeRecording_(P); % disable parfor
    end
end
[vrRho, vrDelta, viNneigh, memory_sort, nFeatures] = sort_page_(S0, P, S_drift);
miKnn = [];
% end

% output
vrRho = vrRho / max(vrRho);     % divide by 10 to be compatible with previous version displays
[~, ordrho] = sort(vrRho, 'descend');
S_clu = struct('rho', vrRho, 'delta', vrDelta, 'ordrho', ordrho, 'nneigh', viNneigh, ...
    'P', P, 'miKnn', miKnn, 'S_drift', S_drift, 'nFeatures', nFeatures);
S_clu.runtime_sort = toc(runtime_sort);
S_clu.memory_sort = memory_sort;
end %func


%--------------------------------------------------------------------------
function [mlPc, nFeatures] = get_mlPc_(S0, P)
nPcPerChan = get_set_(P, 'nPcPerChan', 0);
mlPc = calc_mlPc_(nPcPerChan, S0.dimm_fet);
nC_max = P.nC_max;
nFeatures = sum(mlPc(:));
if nFeatures > nC_max
    fprintf('get_pc_sort_: feature trimmed %d->%d\n', nFeatures, nC_max);
    nFeatures = nC_max;
    vi_ = find(mlPc(:));
    mlPc(vi_(nC_max+1:end)) = false;
end
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
function vcFile_knn_site1 = save_miKnn_site_(vcFile_prm, iSite1, miKnn1, fAppend)
% [Usage]
% save_miKnn_site_(vcFile_prm, iSite1, miKnn1)
% save_miKnn_site_(fid, iSite1, miKnn1)
if nargin<4, fAppend = false; end
vcFile_knn_site1 = strrep(vcFile_prm, '.prm', sprintf('_knn_%d.irc', iSite1));

if fAppend
    fid1 = fopen(vcFile_knn_site1, 'a');
    write_bin_(fid1, miKnn1);    
    fclose(fid1);
else
    write_bin_(vcFile_knn_site1, miKnn1);
end
end %func


%--------------------------------------------------------------------------
function [miKnn1, vl] = load_miKnn_site_(P, iSite1, viSite_spk, viSpk1)
if nargin<3, viSite_spk = []; end
if nargin<4, viSpk1 = []; end

vcFile_knn_site1 = strrep(P.vcFile_prm, '.prm', sprintf('_knn_%d.irc', iSite1));
if iscell(viSite_spk)
    viSpk_site1 = viSite_spk{iSite1};
    dimm1 = [P.knn, numel(viSpk_site1)];
elseif ~isempty(viSite_spk)
    viSpk_site1 = find(viSite_spk==iSite1);
    dimm1 = [P.knn, numel(viSpk_site1)];
else
    viSpk_site1 = [];
    dimm1 = P.knn;
end
miKnn1 = load_bin_(vcFile_knn_site1, 'int64', dimm1);
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
%    load given spikes
% cmiKnn1 = load_miKnn_spk_(P, viSite_spk, cviSpk1)
%    load cell
% miKnn1 = load_miKnn_spk_(P, viSite_spk)
%    load all

if nargin<3, viSpk1 = []; end
fCell = iscell(viSpk1);
if fCell
    vn_cell = cellfun(@numel, viSpk1);
    viSpk1 = cell2mat(viSpk1);
end
if isempty(viSpk1)
    viSite1 = viSite_spk;
else
    viSite1 = viSite_spk(viSpk1);
end
nSites = max(viSite1);
try
    miKnn1 = zeros(get_set_(P, 'knn', 30), numel(viSpk1), 'int64');
catch
    fprintf('load_miKnn_spk_: memory out\n');
    miKnn1 = []; return;
end
% cviSpk_site = 
for iSite = 1:nSites
    viiSpk1 = find(viSite1 == iSite);
    if isempty(viiSpk1), continue; end
    if isempty(viSpk1)
        miKnn1(:,viiSpk1) = load_miKnn_site_(P, iSite);
    else
        miKnn1(:,viiSpk1) = load_miKnn_site_(P, iSite, viSite_spk, viSpk1(viiSpk1));
    end
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
function [vrRho, vrDelta, viNneigh, memory_sort, nFeatures] = sort_page_(S0, P, S_drift)

fDebug = 0;
fSkip_rho = 0;

[viSite_spk, viSite2_spk] = get_(S0, 'viSite_spk', 'viSite2_spk');
t_fun = tic;
[mlPc, nFeatures] = get_mlPc_(S0, P);    
[viLim_drift, mlDrift] = get_(S_drift, 'viLim_drift', 'mlDrift');
nSpk = numel(viSite_spk);
nSites = size(P.miSites,2);
% fParfor = get_set_(P, 'fParfor', 0) && nSites>1;
% fGpu = get_set_(P, 'fGpu', 1);
vcFile_prm = P.vcFile_prm;
% if fDebug, nSites = 2; end

S_page = makeStruct_(P, mlDrift, mlPc, viLim_drift, nSites, vcFile_prm);
S_page = struct_merge_(S_page, struct_copy_(S0, ...
    'type_fet', 'dimm_fet', 'ccviSpk_site_load', 'ccviSpk_site2_load'));

% return schedules
% lim2range_ = @(lim)lim(1):lim(2);
[miSpk_lim_out, miSpk_lim_in, miDrift_lim_out, miDrift_lim_in] = plan_sort_page_(S_drift, P);
nPages = size(miSpk_lim_out,1);
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
viNneigh = zeros(nSpk, 1, 'int64');
S_global = makeStruct_(S_drift, miSpk_lim_out, miSpk_lim_in, miDrift_lim_out, viSite_spk, viSite2_spk);

% clear _miKnn_site_#.irc and append to these files
if ~fSkip_rho
    vcFile_miKnn = [strrep(vcFile_prm, '.prm', ''), '_knn_*.irc'];
    delete_(vcFile_miKnn);
%     assert(isempty(dir(vcFile_miKnn)), sprintf('sort_long_: %s must be deleted', vcFile_miKnn));
    fprintf('sort_page_: calculating Rho...\n'); t_rho = tic;
    for iPage = 1:nPages  
        fprintf('Page %d/%d ', iPage, nPages); t_ = tic;
        [S_page1, viSpk_in1, viSpk_out1] = prepare_page_(S_page, S_global, iPage);
        vrRho(viSpk_in1) = rho_page_(S_page1);
        fprintf('\n\ttook %0.1fs\n', toc(t_));
    end %for
    fprintf('calculating Rho took %0.1fs\n', toc(t_rho));
    if fDebug
        S_page.vrRho = vrRho;    
        struct_save_(S_page, 'S_sort.mat');
        fprintf(2, 'DEBUG: Saved to S_sort.mat\n');
    end
else
    load S_sort.mat
end

fprintf('sort_page_: calculating Delta...\n'); t_delta = tic;
for iPage = 1:nPages    
    fprintf('Page %d/%d ', iPage, nPages); t_ = tic;    
    [S_page1, viSpk_in1, viSpk_out1] = prepare_page_(S_page, S_global, iPage);
    [vrDelta(viSpk_in1), viNneigh(viSpk_in1)] = delta_page_(S_page1, vrRho(viSpk_out1));
%     [vrDelta(viSpk_in1), viNneigh(viSpk_in1)] = delta_page_(S_page1, vrRho);
    fprintf('\n\ttook %0.1fs\n', toc(t_));
end %for
fprintf('calculating Delta took %0.1fs\n', toc(t_delta));

fprintf('sort_long_: took %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t_fun), P.fGpu, P.fParfor);
memory_sort = memory_matlab_();
end %func


%--------------------------------------------------------------------------
function [S_sort1, viSpk_in1, viSpk_out1] = prepare_page_(S_sort, S_global, iPage)

import_struct_(S_global);
lim2range_page_ = @(miLim)miLim(iPage,1):miLim(iPage,end);
[viSpk_out1, viSpk_in1, viDrift_out1] = ...
    deal(lim2range_page_(miSpk_lim_out), lim2range_page_(miSpk_lim_in), lim2range_page_(miDrift_lim_out));
viLim_drift1 = S_drift.viLim_drift(viDrift_out1(1):viDrift_out1(end)+1);
mlDrift1 = S_drift.mlDrift(viDrift_out1, viDrift_out1);
[viSite_in1, viSite_out1] = deal(viSite_spk(viSpk_in1), viSite_spk(viSpk_out1));
if ~isempty(viSite2_spk)
    [viSite2_in1, viSite2_out1] = deal(viSite2_spk(viSpk_in1), viSite2_spk(viSpk_out1));
else
    [viSite2_in1, viSite2_out1] = deal([]);
end
S_sort1 = struct_add_(S_sort, viSpk_out1, viSpk_in1, mlDrift1, ...
    viSite_in1, viSite2_in1, viSite_out1, viSite2_out1, viLim_drift1, iPage);
end %func


%--------------------------------------------------------------------------
function vrRho_in1 = rho_page_(S_page1)

[viSpk_in1, nSites, P, viSite_in1, viSite_out1, viSite2_in1, viSite2_out1] = ...
    get_(S_page1, 'viSpk_in1', 'nSites', 'P', 'viSite_in1', 'viSite_out1', 'viSite2_in1', 'viSite2_out1');
cvrRho_in1 = cell(nSites, 1);
nSpk1 = numel(viSpk_in1);
cviiSpk_in1_site = vi2cell_(viSite_in1, nSites);
cviiSpk_out1_site = vi2cell_(viSite_out1, nSites);
cviiSpk2_in1_site = vi2cell_(viSite2_in1, nSites);
cviiSpk2_out1_site = vi2cell_(viSite2_out1, nSites);
csName_site1 = {'viiSpk_in1', 'viiSpk_out1', 'viiSpk2_in1', 'viiSpk2_out1'};
fParfor = get_set_(P, 'fParfor') && nSites>1;

if fParfor
    try
        parfor iSite = 1:nSites
            try
                S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
                    cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);
                cvrRho_in1{iSite} = rho_paged_site_(S_page1, S_site1, iSite);
            catch
                fprintf('x');
            end
        end %for
    catch
    end
end

vrRho_in1 = zeros(nSpk1, 1, 'single');
for iSite = 1:nSites
    if isempty(cviiSpk_in1_site{iSite}), continue; end
    if isempty(cvrRho_in1{iSite})
%         fprintf('\tSite %d: ', iSite); t1=tic;
        S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
            cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);        
        cvrRho_in1{iSite} = rho_paged_site_(S_page1, S_site1, iSite);        
%         fprintf('took %0.1fs\n', toc(t1));
        fprintf('.');
    end
    vrRho_in1(cviiSpk_in1_site{iSite}) = cvrRho_in1{iSite};
end
end %func


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1] = delta_page_(S_page1, vrRho_page1)

[viSpk_in1, nSites, P, viSite_in1, viSite_out1, viSite2_in1, viSite2_out1] = ...
    get_(S_page1, 'viSpk_in1', 'nSites', 'P', 'viSite_in1', 'viSite_out1', 'viSite2_in1', 'viSite2_out1');
[cvrDelta_in1, cviNneigh_in1] = deal(cell(nSites, 1));
nSpk1 = numel(viSpk_in1);
cviiSpk_in1_site = vi2cell_(viSite_in1, nSites);
cviiSpk_out1_site = vi2cell_(viSite_out1, nSites);
cviiSpk2_in1_site = vi2cell_(viSite2_in1, nSites);
cviiSpk2_out1_site = vi2cell_(viSite2_out1, nSites);
csName_site1 = {'viiSpk_in1', 'viiSpk_out1', 'viiSpk2_in1', 'viiSpk2_out1'};

if get_set_(P, 'fParfor', 1)
    try
        parfor iSite = 1:nSites
            try
                S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
                    cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);        
                [cvrDelta_in1{iSite}, cviNneigh_in1{iSite}] = ...
                    delta_paged_site_(S_page1, S_site1, vrRho_page1, iSite);
            catch
                fprintf('x');
            end
        end %for
    catch
    end
end

vrDelta1 = zeros(nSpk1, 1, 'single');
viNneigh1 = zeros(nSpk1, 1, 'int64');
for iSite = 1:nSites
    if isempty(cviiSpk_in1_site{iSite}), continue; end
    if isempty(cvrDelta_in1{iSite})
%         fprintf('\tSite %d: ', iSite); t1=tic;
        S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
            cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);                
        [cvrDelta_in1{iSite}, cviNneigh_in1{iSite}] = ...
            delta_paged_site_(S_page1, S_site1, vrRho_page1, iSite);
%         fprintf('took %0.1fs\n', toc(t1));
        fprintf('.');
    end
    vrDelta1(cviiSpk_in1_site{iSite}) = cvrDelta_in1{iSite};
    viNneigh1(cviiSpk_in1_site{iSite}) = cviNneigh_in1{iSite};
end
end %func


%--------------------------------------------------------------------------
function [csName, cVal] = import_struct_(varargin)
% [Usage]
% import_struct_(S)
%    import all fields from the struct S
% import_struct_(S, 'var1', 'var2', ...)
%    import specified variables from the struct S

S = varargin{1};
if ~isstruct(S), error('import_struct_: first argument must be a struct'); end

% assignes args to the caller
if nargin<2
    csName = fieldnames(S); 
else
    csName = varargin(2:end);
end
cVal = struct2cell(S);
for i=1:numel(csName)
    assignin('caller', csName{i}, cVal{i});
end
end %func


%--------------------------------------------------------------------------
function vrRho_in = rho_paged_site_(S_page, S_site, iSite)

vrRho_in=[]; 
if isempty(S_site.viiSpk_in1), return; end

csVar = import_struct_(prepare_page_site_(S_page, S_site, iSite));

fGpu = get_set_(P, 'fGpu', 1);
if fGpu
    try
        [vrRho_in, miKnn_in] = search_knn_drift_(vl_in, mrFet_out, viDrift_out, mlDrift1, P);
    catch
        fGpu=0; 
        fprintf(2, 'CPU');
    end
end
if ~fGpu
    knn = get_set_(P, 'knn', 30);
    [n_in, nDrift, viDrift_in] = deal(sum(vl_in), size(mlDrift1,1), viDrift_out(vl_in));
    cviiSpk_in_drift = vi2cell_(viDrift_in, nDrift);
    cviiSpk_out_drift = vi2cell_(viDrift_out, nDrift);
    vrRho_in = zeros(n_in,1,'single');
    miKnn_in = zeros(knn, n_in, 'int64');
    mrFet_in = mrFet_out(:,vl_in);
    for iDrift = 1:nDrift
        vii_in1 = cviiSpk_in_drift{iDrift};      
        if isempty(vii_in1), continue; end
        vii_out1 = cell2mat_(cviiSpk_out_drift(mlDrift1(:, iDrift)));
        [vrRho_in(vii_in1), miKnn_] = knn_cpu_(mrFet_out(:,vii_out1), mrFet_in(:,vii_in1), knn);
        miKnn_in(:,vii_in1) = vii_out1(miKnn_);
    end
end
vrRho_in = 1 ./ vrRho_in;
miKnn_in = int64(viSpk_out(miKnn_in));    

% save miKnn
n2 = size(miKnn_in,1);
if n2 < get_set_(P, 'knn', 30)       
    miKnn_in = [miKnn_in; repmat(miKnn_in(end,:), [knn-n2, 1])];
end
iPage = get_(S_page, 'iPage');
fAppend = iPage > 1;
save_miKnn_site_(S_page.vcFile_prm, iSite, miKnn_in, fAppend);
end %func


%--------------------------------------------------------------------------
function [vrDelta_in, viNneigh_in] = delta_paged_site_(S_page1, S_site1, vrRho_page1, iSite)

[vrDelta_in, viNneigh_in] = deal([]);
if isempty(S_site1.viiSpk_in1), return; end
SINGLE_INF = 3.402E+38;

csVar = import_struct_(prepare_page_site_(S_page1, S_site1, iSite));
vrRho_out = vrRho_page1(viiSpk_out);  
vrRho_in = vrRho_out(vl_in);

fGpu = get_set_(P, 'fGpu', 1);
if fGpu    
    try
        [vrDelta_in, viNneigh_in] = search_delta_drift_(vl_in, mrFet_out, vrRho_out, viDrift_out, mlDrift1, P);
    catch
        fGpu = 0;
        fprintf(2, 'CPU');
    end
end
if ~fGpu
    [n_in, nDrift, viDrift_in, vrRho_in] = ...
        deal(sum(vl_in), size(mlDrift1,1), viDrift_out(vl_in), vrRho_out(vl_in));
    mrFet_in = mrFet_out(:,vl_in);
    cviiSpk_in_drift = vi2cell_(viDrift_in, nDrift);
    cviiSpk_out_drift = vi2cell_(viDrift_out, nDrift);
    vrDelta_in = zeros(n_in, 1,'single');
    viNneigh_in = zeros(n_in, 1, 'int64');
    for iDrift = 1:nDrift
        vii_in1 = cviiSpk_in_drift{iDrift};        
        if isempty(vii_in1), continue; end        
        vii_out1 = cell2mat_(cviiSpk_out_drift(mlDrift1(:, iDrift)));
        [vrDelta_in(vii_in1), viNneigh_] = ...
            delta_cpu_(mrFet_out(:,vii_out1), mrFet_in(:,vii_in1), vrRho_out(vii_out1), vrRho_in(vii_in1));
        viNneigh_in(vii_in1) = vii_out1(viNneigh_);
    end
end
vrDelta_in = vrDelta_in .* vrRho_in;
viNneigh_in= int64(viSpk_out(viNneigh_in));

viNan = find(isnan(vrDelta_in) | isinf(vrDelta_in));
viNneigh_in(viNan) = viNan;
vrDelta_in(viNan) = sqrt(SINGLE_INF);
end %func


%--------------------------------------------------------------------------
% uses search_min_drift.cu
function [vrKnn_in, miKnn_in] = search_knn_drift_(vl_in, mrFet_out, viDrift_out, mlDrift1, P)

knn = get_set_(P, 'knn', 30);
n_in = sum(vl_in);
[CHUNK, nC_max, nThreads] = deal(8, 45, 1024); % tied to search_min_drift.cu
[nC, nF] = deal(size(mrFet_out,1), size(mrFet_out,2));

viDrift_in = viDrift_out(vl_in);
vrKnn_in = zeros(n_in, 1, 'single');
miKnn_in = zeros(knn, n_in, 'int32');
gmrFet_out = gpuArray_(mrFet_out, P.fGpu);
viAA = min(viDrift_in):max(viDrift_in);
viiBB = [1, find(diff(viDrift_out(:)'))+1, numel(viDrift_out)+1];
viBB = viDrift_out(viiBB(1:end-1));
vnBB = diff(viiBB);
in_offset = find(vl_in, 1, 'first') - 1;

CK = parallel.gpu.CUDAKernel('search_min_drift.ptx','search_min_drift.cu'); % auto-compile if ptx doesn't exist
CK.ThreadBlockSize = [nThreads, 1];          
CK.SharedMemorySize = 4 * CHUNK * (nC_max + 1); % @TODO: update the size

                
for iiAA = 1:numel(viAA)
    % determine index
    iAA1 = viAA(iiAA);
    viiAA1 = find(viDrift_in == iAA1);
    limAA1 = viiAA1([1, end]) + in_offset;
    viiBB1 = find(ismember(viBB, find(mlDrift1(:,iAA1))));
    [viBB1, vnBB1] = deal(viiBB(viiBB1), vnBB(viiBB1));  
    [iA1, nA1, nBB1] = deal(limAA1(1), diff(limAA1)+1, numel(viBB1));
    
    % call cuda function
    CK.GridSize = [ceil(nA1 / CHUNK / CHUNK), CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]    
    gmrMin1 = zeros([nThreads, nA1], 'single', 'gpuArray');
    gmiMin1 = zeros([nThreads, nA1], 'int32', 'gpuArray'); 
    gviBB1 = gpuArray(int32(viBB1)-1); % c indexing starts with 0
    gvnBB1 = gpuArray(int32(vnBB1));        
    gvnConst = gpuArray(int32([nC, nF, nBB1, iA1-1, nA1])); % c indexing sarts with zero
    [gmrMin1, gmiMin1] = feval(CK, gmrMin1, gmiMin1, gmrFet_out, gviBB1, gvnBB1, gvnConst);

    % process output
    gmiMin1 = gather(gmiMin1);
    [gmrMin1, miKnn1] = sort(gather(gmrMin1));
    vrKnn_in(viiAA1) = gmrMin1(knn,:)';
    miKnn_in(:,viiAA1) = gmiMin1(miKnn1(1:knn,:) + (0:nA1-1)*nThreads);    
%         disp(iiAA);
%         fprintf('%s\n', sprintf('%d ', viBB1));
%         fprintf('%s\n', sprintf('%d ', vnBB1));
%         disp(limAA1');
%         fprintf('\n');
end
end %func


%--------------------------------------------------------------------------
% uses search_min_drift.cu
function [vrDelta_in, viNneigh_in] = search_delta_drift_(vl_in, mrFet_out, vrRho_out, viDrift_out, mlDrift1, P)

knn = get_set_(P, 'knn', 30);
n_in = sum(vl_in);
[CHUNK, nC_max, nThreads] = deal(16, 45, 512); % tied to cuda_knn_index.cu
[nC, nF] = deal(size(mrFet_out,1), size(mrFet_out,2));

viDrift_in = viDrift_out(vl_in);
vrDelta_in = zeros(n_in, 1, 'single');
viNneigh_in = zeros(n_in, 1, 'int32');
gmrFet_out = gpuArray(mrFet_out);
gvrRho_out = gpuArray(vrRho_out);
viAA = min(viDrift_in):max(viDrift_in);
viiBB = [1, find(diff(viDrift_out(:)'))+1, numel(viDrift_out)+1];
viBB = viDrift_out(viiBB(1:end-1));
vnBB = diff(viiBB);
in_offset = find(vl_in, 1, 'first') - 1;

CK = parallel.gpu.CUDAKernel('search_delta_drift.ptx','search_delta_drift.cu'); % auto-compile if ptx doesn't exist
CK.ThreadBlockSize = [nThreads, 1];          
CK.SharedMemorySize = 4 * CHUNK * (nC_max + 2); % @TODO: update the size

                
for iiAA = 1:numel(viAA)
    % determine index
    iAA1 = viAA(iiAA);
    viiAA1 = find(viDrift_in == iAA1);
    limAA1 = viiAA1([1, end]) + in_offset;
    viiBB1 = find(ismember(viBB, find(mlDrift1(:,iAA1))));
    [viBB1, vnBB1] = deal(viiBB(viiBB1), vnBB(viiBB1));  
    [iA1, nA1, nBB1] = deal(limAA1(1), diff(limAA1)+1, numel(viBB1));
    
    % call cuda function
    CK.GridSize = [ceil(nA1 / CHUNK / CHUNK), CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]    
    gmrMin1 = zeros([nThreads, nA1], 'single', 'gpuArray');
    gmiMin1 = zeros([nThreads, nA1], 'int32', 'gpuArray'); 
    gviBB1 = gpuArray(int32(viBB1)-1); % c indexing starts with 0
    gvnBB1 = gpuArray(int32(vnBB1));        
    gvnConst = gpuArray(int32([nC, nF, nBB1, iA1-1, nA1])); % c indexing sarts with zero
    [gmrMin1, gmiMin1] = feval(CK, gmrMin1, gmiMin1, gmrFet_out, gvrRho_out, gviBB1, gvnBB1, gvnConst);

    % process output
    if false
        gmrMin1=gather(gmrMin1); gmiMin1=gather(gmiMin1);
        [vrMin1, viMin1] = min(gmrMin1);
        vrDelta_in(viiAA1) = vrMin1';
        viNneigh_in(viiAA1) = gmiMin1(viMin1 + (0:nA1-1)*nThreads)';    
    else
        [vrMin1, viMin1] = min(gmrMin1);
        vrDelta_in(viiAA1) = gather(vrMin1');
        viNneigh_in(viiAA1) = gather(gmiMin1(viMin1 + (0:nA1-1)*nThreads)');    
    end
%         disp(iiAA);
%         fprintf('%s\n', sprintf('%d ', viBB1));
%         fprintf('%s\n', sprintf('%d ', vnBB1));
%         disp(limAA1');
%         fprintf('\n');
end
end %func


%--------------------------------------------------------------------------
function S_page_site = prepare_page_site_(S_page, S_site, iSite)
% S_sort:     
%    'P'    'mlDrift'    'mlPc'    'viLim_drift'    'nSites'    'fParfor'    'fGpu'    'vcFile_prm'    'type_fet'    'dimm_fet'    'ccviSpk_site_load'    'ccviSpk_site2_load'
%    'viSpk_out1'    'viSpk_in1'    'mlDrift1'    'viSite_in1'    'viSite2_in1'    'viSite_out1'    'viSite2_out1'    'viLim_drift1'
% S_site : 
%    'viiSpk_in1'    'viiSpk_out1'    'viiSpk2_in1'    'viiSpk2_out1'

fDebug = 0;
[viSpk_in1, viSpk_out1, viLim_drift1] = get_(S_page, 'viSpk_in1', 'viSpk_out1', 'viLim_drift1');
[viiSpk_in1, viiSpk_out1, viiSpk2_in1, viiSpk2_out1] = ...
    get_(S_site, 'viiSpk_in1', 'viiSpk_out1', 'viiSpk2_in1', 'viiSpk2_out1');
[mlDrift1, P] = get_(S_page, 'mlDrift1', 'P');

% load fet
[viSpk_in, viSpk_out] = deal(viSpk_in1(viiSpk_in1), viSpk_out1(viiSpk_out1));
[viSpk_in, viSpk_out] = deal(viSpk_in(:), viSpk_out(:));
mrFet_out = load_fet_site_(S_page, 1, iSite, viSpk_out);
vl_in = viSpk_out >= viSpk_in(1) & viSpk_out <= viSpk_in(end);
% mrFet_in = mrFet_out(:,vl_in);
viiSpk_in = viiSpk_out1(vl_in);
% if fDebug
%     assert(all(viiSpk_in(:)==viiSpk_in1(:)), 'prepare_page_site_: must be the same'); 
% end
viSpk2_out = viSpk_out1(viiSpk2_out1);
mrFet2_out = load_fet_site_(S_page, 2, iSite, viSpk2_out);
if ~isempty(mrFet2_out)
    mrFet_out = [mrFet_out, mrFet2_out];    
    viSpk_out = [viSpk_out(:); viSpk2_out(:)];
    viiSpk_out = [viiSpk_out1(:); viiSpk2_out1(:)];
else
    viiSpk_out = viiSpk_out1(:);
end
if fDebug
    equal_ = @(a,b)all(a(:)==b(:));
    assert(equal_(mrFet_in, load_fet_site_(S_page, 1, iSite, viSpk_in)), 'prepare_paged_site_: must be the same');
    assert(size(mlDrift1,1) == size(mlDrift1,2), 'prepare_paged_site_: mlDrift1 must be square matrix');
end
% viDrift_in = discretize(viSpk_in, viLim_drift1); % 0-base for GPU
viDrift_out = discretize(viSpk_out, viLim_drift1); % 0-base for GPU
S_page_site = makeStruct_(mrFet_out, viDrift_out, ...
    mlDrift1, viSpk_in, viSpk_out, viiSpk_in, viiSpk_out, vl_in, P);
end %func


%--------------------------------------------------------------------------
function [miSpk_lim_out, miSpk_lim_in, miDrift_lim_out, miDrift_lim_in] = plan_sort_page_(S_drift, P)
fDebug = 0;

n = P.nTime_drift;
width = get_set_(P, 'nTime_max_drift', 2^8);
viLim_drift = S_drift.viLim_drift;

[quarter_, half_] = deal(width/4, width/2);
nPages = floor(n / half_);
if nPages==0
    [miSpk_lim_out, miSpk_lim_in] = deal(int64([viLim_drift(1), viLim_drift(end)-1]));
    [miDrift_lim_out, miDrift_lim_in] = deal(int32([1, numel(viLim_drift)-1]));
else
    ib_prev = 0;
    [miSpk_lim_out, miSpk_lim_in] = deal(zeros(nPages, 2, 'int64'));
    [miDrift_lim_out, miDrift_lim_in] = deal(zeros(nPages, 2, 'int32'));
    for iPage = 1:nPages
        iB_ = min((iPage-1) * half_ + width, n);
        iA_ = max(1, iB_ - width+1);
        ia_ = ifeq_(iPage==1, 1, min(iA_+quarter_, n));
        ia_ = max(ib_prev+1, ia_);
        ib_ = ifeq_(iPage==nPages, n, max(iB_-quarter_, 1));
        if fDebug
            fprintf('%d: %d-%d-%d-%d\n', iPage, iA_, ia_, ib_, iB_);
        end
        ib_prev = ib_;
        miDrift_lim_out(iPage, :) = [iA_, iB_];
        miDrift_lim_in(iPage, :) = [ia_, ib_];
        miSpk_lim_in(iPage,:) = [viLim_drift(ia_), viLim_drift(ib_+1)-1];
        miSpk_lim_out(iPage,:) = [viLim_drift(iA_), viLim_drift(iB_+1)-1];
    end
end
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
function miKnn = miKnn_format_(miKnn, knn)
n = size(miKnn,1);
if n < knn       
    miKnn = [miKnn; repmat(miKnn(end,:), [knn-n, 1])];
elseif n > knn
    miKnn = miKnn(1:knn,:);
end
end %func


%--------------------------------------------------------------------------
function c_mmf = save_fet_drift_(vcFile_fet1, mrFet1, cvii1_drift)
fid_w = fopen(vcFile_fet1, 'w');
for iFile = 1:numel(cvii1_drift)
    write_bin_(fid_w, mrFet1(:, cvii1_drift{iFile}));
end
fclose(fid_w);
if nargout==0, return; end

dimm1 = size(mrFet1);
nDrift = numel(cvii1_drift);
vn_drift = cellfun(@numel, cvii1_drift);
viOffset_drift = cumsum([0; vn_drift(:)]);
c_mmf = arrayfun_(@(i)memmapfile(...
    vcFile_fet1, 'Format', class(mrFet1), 'Offset', viOffset_drift(i), ...
    'Repeat', vn_drift(i) * prod(dimm1(1:end-1))), 1:nDrift);
end %func


%--------------------------------------------------------------------------
% 9/20/2018 JJJ: Memory-optimized knn for computing delta (dpclus)
function [vrDelta1, viNneigh1] = delta_cpu_(mrFet2, mrFet1, vrRho2, vrRho1)
nStep_knn = 1000;
n1 = numel(vrRho1);
vrDelta1 = zeros([n1, 1], 'single');
viNneigh1 = zeros([n1, 1], 'int32');
fh_dist_ = @(y)bsxfun(@plus, sum(y.^2), bsxfun(@minus, sum(mrFet2.^2)', 2*mrFet2'*y));
for i1 = 1:nStep_knn:n1
    vi1_ = i1:min(i1+nStep_knn-1, n1);
    mrD_ = fh_dist_(mrFet1(:,vi1_));
    mrD_(bsxfun(@le, vrRho2, vrRho1(vi1_)')) = inf;
    [vrDelta1(vi1_), viNneigh1(vi1_)] = min(mrD_);
end
vrDelta1 = sqrt(abs(vrDelta1));
end %func


%--------------------------------------------------------------------------
% load specific site and drift
% todo: speed up by using cache
function [dimm_fet1, viSpk1, vcFile_fet1] = dimm_fet_site_(S_fet, iFet, iSite)
% usage
% ----
% [out1, viSpk] = load_fet2pc_site_(S_fet, iFet, iSite) % returns all spk
% [out1] = load_fet2pc_site_(S_fet, iFet, iSite, viSpk) % returns subset

[dimm_fet1, ccviSpk_site_load, mlPc, vcFile_prm] = struct_get_(S_fet, 'dimm_fet', 'ccviSpk_site_load', 'mlPc', 'vcFile_prm');
nLoads = numel(ccviSpk_site_load);
vcFile_prm_ = strrep(vcFile_prm, '.prm', '');
switch iFet
    case 1
        cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site_load);
        vcFile_fet1 = sprintf('%s_fet1_drift_site_%d.irc', vcFile_prm_, iSite);
    case 2
        cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site2_load);
        vcFile_fet1 = sprintf('%s_fet2_drift_site_%d.irc', vcFile_prm_, iSite);
end
viSpk1 = cell2mat_(arrayfun_(@(iLoad)cviSpk_load{iLoad}, (1:nLoads)'));
if ~isempty(mlPc)
    dimm_fet1 = [sum(mlPc(:)), numel(viSpk1)];
else
    dimm_fet1 = [dimm_fet1(1), dimm_fet1(2), numel(viSpk1)];
end
end %func


%--------------------------------------------------------------------------
% load specific site and drift
% todo: speed up by using cache
function [out1, viSpk] = load_fet_site_(S_fet, iFet, iSite, viSpk)
% usage
% ----
% [out1, viSpk] = load_fet2pc_site_(S_fet, iFet, iSite) % returns all spk
% [out1] = load_fet2pc_site_(S_fet, iFet, iSite, viSpk) % returns subset

if nargin<4, viSpk = []; end

[type_fet, dimm_fet, vcFile_prm, mlPc] = ...
    struct_get_(S_fet, 'type_fet', 'dimm_fet', 'vcFile_prm', 'mlPc');
vcFile_prm_ = strrep(vcFile_prm, '.prm', '');
bytes_per_spk = bytesPerSample_(type_fet) * dimm_fet(1) * dimm_fet(2);
fh_sum = @(x,y)sum(cellfun(@numel, x(1:y-1)));
nLoads = numel(S_fet.ccviSpk_site_load);
switch iFet
    case 1
        vnSpk_load = cellfun(@(x)numel(x{iSite}), S_fet.ccviSpk_site_load);
        vnBytes_offset_load = cellfun(@(x)fh_sum(x, iSite), S_fet.ccviSpk_site_load) * bytes_per_spk;
        csFiles_fet = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet_%d.irc',x)], 1:nLoads);
        cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site_load);
    case 2
        vnSpk_load = cellfun(@(x)numel(x{iSite}), S_fet.ccviSpk_site2_load);
        csFiles_fet = arrayfun_(@(x)[vcFile_prm_, sprintf('_fet2_%d.irc',x)], 1:nLoads);
        vnBytes_offset_load = cellfun(@(x)fh_sum(x, iSite), S_fet.ccviSpk_site2_load) * bytes_per_spk;
        cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site2_load);
end
fLoad_all = isempty(viSpk);
if fLoad_all
    viSpk = cell2mat_(arrayfun_(@(iLoad)cviSpk_load{iLoad}, (1:nLoads)'));
end
if isempty(viSpk), out1=[]; return; end
if ~isempty(mlPc)
    nFeatures = sum(mlPc(:));
    tr2mr_ = @(x)reshape(x, [], size(x,3));
    fet2pc_ = @(x)reshape(x(mlPc(:),:), [nFeatures, size(x,2)]); 
    mrFet1 = zeros(nFeatures, numel(viSpk), type_fet);
else
    trPc1 = zeros(dimm_fet(1), dimm_fet(2), numel(viSpk), type_fet);
end
[iOffset_spk, viSpk_lim] = deal(0, [min(viSpk), max(viSpk)]);
for iLoad = 1:numel(csFiles_fet)
    if ~fLoad_all
        viSpk1 = [];
        viSpk_load1 = cviSpk_load{iLoad};
        if isempty(viSpk_load1), continue; end
        if a_in_b_(viSpk_lim, viSpk_load1, 0)
            [vl_, vi_] = ismember(viSpk, cviSpk_load{iLoad});
            if any(vl_)
                viiSpk_ = vi_(vl_);
                viSpk1 = (1:numel(viiSpk_)) + iOffset_spk;        
            end
        elseif viSpk_load1(1) > viSpk_lim(2)
            break;
        end
    else
        viSpk1 = (1:vnSpk_load(iLoad)) + iOffset_spk;
    end
    if isempty(viSpk1), continue; end    
    dimm_fet1 = [dimm_fet(1), dimm_fet(2), vnSpk_load(iLoad)];
    fid1=fopen(csFiles_fet{iLoad},'r'); 
    fseek(fid1, vnBytes_offset_load(iLoad), 'bof'); 
    trPc_load1 = fread_(fid1, dimm_fet1, type_fet);
    fclose(fid1);
    
    if ~isempty(mlPc)
        if fLoad_all
            mrFet1(:,viSpk1) = fet2pc_(tr2mr_(trPc_load1));            
        else
            mrFet1(:,viSpk1) = fet2pc_(tr2mr_(trPc_load1(:,:,viiSpk_)));
        end
    else
        if fLoad_all
            trPc1(:,:,viSpk1) = trPc_load1;
        else
            trPc1(:,:,viSpk1) = trPc_load1(:,:,viiSpk_);
        end
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
function flag = a_in_b_(vrA, vrB, fSorted)
if nargin<3, fSorted = []; end
if isempty(vrA) || isempty(vrB), flag = false; return; end
% assume sorrted
if fSorted
    limA = [vrA(1), vrA(end)];
    limB = [vrB(1), vrB(end)];
else
    limA = [min(vrA), max(vrA)];
    limB = [min(vrB), max(vrB)];    
end
flag = (limA(1) >= limB(1) && limA(1) <= limB(2)) || ...
       (limA(2) >= limB(1) && limA(2) <= limB(2)) || ...
       (limB(1) >= limA(1) && limB(1) <= limA(2)) || ...
       (limB(2) >= limA(1) && limB(2) <= limA(2));
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
    [mrKnn1, miKnn_] = pdist2(mrFet2_T, mrFet1(:,vi1_)', 'euclidean', 'Smallest', knn);
    vrKnn(vi1_) = mrKnn1(end,:);
    knn_ = size(miKnn_,1);
    if knn_ == knn
        miKnn(:,vi1_) = miKnn_;
    else
        miKnn(1:knn_,vi1_) = miKnn_;
        miKnn(knn_+1:end,vi1_) = repmat(miKnn_(end,:), [knn-knn_,1]);
    end    
end %for
end %func


%--------------------------------------------------------------------------
% 10/16/2018 JJJ: Using native function in matlab (pdist2 function)
% 9/20/2018 JJJ: Memory-optimized knn
function [mrKnn, miKnn] = knn_cpu_full_(mrFet2, mrFet1, knn)

nStep_knn = 1000;
mrFet2_T = mrFet2';
n1 = size(mrFet1,2);
miKnn = int32(zeros([knn, n1], 'like', mrFet1));
mrKnn = zeros([knn, n1], 'like', mrFet1);
for i1 = 1:nStep_knn:n1
    vi1_ = i1:min(i1+nStep_knn-1, n1);
    [mrKnn1, miKnn(:,vi1_)] = pdist2(mrFet2_T, mrFet1(:,vi1_)', 'euclidean', 'Smallest', knn);
    mrKnn(:,vi1_) = mrKnn1;
end %for
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
cS_detect = cell(nLoads, 1);
disp_load_ = @(iPage, var_size, t_load, nSpikes)...
    fprintf('\tDetecting %d/%d: %d spikes found (%0.1f spikes/s, %0.1f MB/s, took %0.1f s)\n', ...
    iPage, nLoads, nSpikes, nSpikes/t_load, var_size/t_load/1e6, t_load);

% process the first load
t_load1 = tic;
[mrWav_T1, nlim_wav1, fDone] = readmda_paged_(); var_size1 = var_size_(mrWav_T1);
cS_detect{1} = detect_paged_(mrWav_T1, P, makeStruct_(nlim_wav1)); mrWav_T1 = [];
disp_load_(1, var_size1, toc(t_load1), numel(get_(cS_detect{1}, 'viTime_spk')));

% extract info from the first load
[vrThresh_site, mrPv_global] = struct_get_(cS_detect{1}, 'vrThresh_site', 'mrPv_global');
S_cache = makeStruct_(vrThresh_site, mrPv_global);
% fprintf('\tMemory use: %0.3f GiB\n', (memory_matlab_()-memory_init)/2^30);
delete_file_fet_(P); % clear fet
[vcFile, vS_load] = readmda_paged_('close'); % close the file
cS_detect{1} = detect_paged_save_(cS_detect{1}, P, 1);    
viSite2Chan = get_(P, 'viSite2Chan');

if ~fDone
    if fParfor
        parfor iLoad = 2:nLoads  % change to for loop for debugging
            t_load1 = tic;
            S_load1 = vS_load(iLoad);
            mrWav_T1 = load_file_part_(vcFile, S_load1, viSite2Chan); var_size1 = var_size_(mrWav_T1);
            S_cache1 = setfield(S_cache, 'nlim_wav1', S_load1.nlim); t_detect1 = tic;
            cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache1);  mrWav_T1 = [];
            cS_detect{iLoad} = detect_paged_save_(cS_detect{iLoad}, P, iLoad);        
            disp_load_(iLoad, var_size1, toc(t_load1), numel(get_(cS_detect{iLoad}, 'viTime_spk')));
        end
    else
        for iLoad = 2:nLoads  % change to for loop for debugging
            t_load1 = tic;
            S_load1 = vS_load(iLoad);
            mrWav_T1 = load_file_part_(vcFile, S_load1, viSite2Chan); var_size1 = var_size_(mrWav_T1);
            S_cache1 = setfield(S_cache, 'nlim_wav1', S_load1.nlim);
            cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache1);  mrWav_T1 = [];
            cS_detect{iLoad} = detect_paged_save_(cS_detect{iLoad}, P, iLoad);      
            disp_load_(iLoad, var_size1, toc(t_load1), numel(get_(cS_detect{iLoad}, 'viTime_spk')));
        end
    end
end
S0 = detect_merge_(cS_detect, viOffset_load, P);
% fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
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
function [cviSpk_site, nSites, vi_site] = vi2cell_(viSite_spk, nSites)
if nargin<2, nSites = []; end
if isempty(nSites), nSites = max(viSite_spk); end

% based on unique() function, which sorts. faster than arrayfun.
cviSpk_site = cell(nSites, 1);
[vr, vi] = sort(viSite_spk);
vi_change = [1; find(diff(vr(:))>0)+1; numel(viSite_spk)+1];
if isempty(viSite_spk), vi_site=[]; return; end

vi_site = vr(vi_change(1:end-1));
vl_remove = vi_site < 1;
if any(vl_remove)
    vi_site(vl_remove) = [];
    vi_change(find(vl_remove)) = [];
end
for iStep = 1:numel(vi_site)
    cviSpk_site{vi_site(iStep)} = vi(vi_change(iStep):vi_change(iStep+1)-1);
end
vi_site = vi_site(:)';

% check equal condition
if false
    cviSpk_site0 = arrayfun(@(x)find(viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
    assert(all(cellfun(@(x,y)all(x==y), cviSpk_site0, cviSpk_site)), 'vi2cell_: must equal');
end
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
vnSpk_load = cellfun(@(x)numel(x.viSite_spk), cS_detect);
miSpk_load = [0; cumsum(vnSpk_load)];
miSpk_load = [miSpk_load(1:end-1)+1, miSpk_load(2:end)];
viSpk_offset_load = cumsum([0; vnSpk_load]);
cell_offset_ = @(x, iLoad)cellfun_(@(y)int64(y) + int64(viSpk_offset_load(iLoad)), x);
nSpk = sum(vnSpk_load);
[viSite_spk, viTime_spk, mrPos_spk] = ...
    deal(zeros(nSpk, 1, 'int32'), zeros(nSpk, 1, 'int64'), zeros(nSpk, 2, 'single'));
[vrPow_spk, vrAmp_spk] = deal(zeros(nSpk, 1, 'single'));
viOffset_load = int64(viOffset_load);
[mrVp_spk, type_fet, dimm_fet, viSite2_spk, trPc_spk, trPc2_spk] = deal([]);
[ccviSpk_site_load, ccviSpk_site2_load] = deal(cell(size(cS_detect)));
for iLoad = 1:numel(cS_detect)
    S1 = cS_detect{iLoad};
    if isempty(type_fet) && isempty(dimm_fet)
        type_fet = S1.type_fet;
        dimm_fet = [S1.dimm_fet(1), S1.dimm_fet(2), nSpk];
    end
    viSpk1 = miSpk_load(iLoad,1):miSpk_load(iLoad,2);
    viSite_spk(viSpk1) = S1.viSite_spk;
    viTime_spk(viSpk1) = int64(S1.viTime_spk) + viOffset_load(iLoad);
    vrAmp_spk(viSpk1) = S1.vrAmp_spk;
    if isfield(S1, 'mrPos_spk') && isfield(S1, 'vrPow_spk')
        [mrPos_spk(viSpk1,:), vrPow_spk(viSpk1)] = deal(S1.mrPos_spk, S1.vrPow_spk);    
    end
    ccviSpk_site_load{iLoad} = cell_offset_(get_(S1, 'cviSpk_site'), iLoad);
    
    % secondary peak   
    if ~isempty(get_(S1, 'viSite2_spk'))
        if isempty(viSite2_spk)
            viSite2_spk = zeros(nSpk, 1, 'int32');
        end
        viSite2_spk(viSpk1) = S1.viSite2_spk;
    end
    ccviSpk_site2_load{iLoad} = cell_offset_(get_(S1, 'cviSpk2_site'), iLoad);
    
end
S0 = makeStruct_(viSite_spk, viTime_spk, vrAmp_spk, mrVp_spk, ...
         viSite2_spk, trPc_spk, trPc2_spk, type_fet, dimm_fet, ...
         mrPos_spk, vrPow_spk, ccviSpk_site_load, ccviSpk_site2_load, vnSpk_load); 
end %func


%--------------------------------------------------------------------------
% Save trPc_spk and trPc2_spk and remove from the struct
function S_detect = detect_paged_save_(S_detect, P, iLoad)

fSave_fet = get_set_(P, 'fSave_fet', 0);
nSites = size(P.miSites, 2);

vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');
if fSave_fet
    write_bin_([vcFile_prm_, sprintf('_fet_%d.irc', iLoad)], S_detect.trPc_spk);    
end
S_detect.type_fet = class(S_detect.trPc_spk);
S_detect.dimm_fet = size(S_detect.trPc_spk);
S_detect.cviSpk_site = save_paged_fet_site_(...
    [vcFile_prm_, sprintf('_fet_%d.irc', iLoad)], ...
        S_detect.trPc_spk, S_detect.viSite_spk, nSites);
S_detect.trPc_spk = [];
if isempty(get_(S_detect, 'trPc2_spk'))
    S_detect.cviSpk2_site = cell(size(S_detect.cviSpk_site));
    return; 
end   
if fSave_fet
    write_bin_([vcFile_prm_, sprintf('_fet2_%d.irc', iLoad)], S_detect.trPc2_spk);    
end
S_detect.cviSpk2_site = save_paged_fet_site_(...
    [vcFile_prm_, sprintf('_fet2_%d.irc', iLoad)], ...
        S_detect.trPc2_spk, S_detect.viSite2_spk, nSites);
S_detect.trPc2_spk = [];
end %func


%--------------------------------------------------------------------------
function cviSpk_site = save_paged_fet_site_(vcFile_out, trFet_spk, viSite_spk, nSites)
% t1=tic;
[cviSpk_site] = vi2cell_(viSite_spk, nSites);
fid_w = fopen(vcFile_out, 'w');
for iSite = 1:nSites
    write_bin_(fid_w, trFet_spk(:,:,cviSpk_site{iSite}));
end
fclose(fid_w);
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
[mrWav_filt, vrWav_mean_filt] = filter_transpose_(mrWav_T, P);
S_detect = get_spikes_(mrWav_filt, vrWav_mean_filt, P, S_cache);
end %func


%--------------------------------------------------------------------------
function size_byte = var_size_(var)
size_byte = numel(var) * bytesPerSample_(class_(var));
end %func


%--------------------------------------------------------------------------
function S_detect = get_spikes_(mrWav_filt, vrWav_mean_filt, P, S_cache)
[vrThresh_site, nlim_wav1, mrPv_global, vrD_global] = ...
    struct_get_(S_cache, 'vrThresh_site', 'nlim_wav1', 'mrPv_global', 'vrD_global');

if isempty(nlim_wav1)
    [nPad_pre, nPad_post] = deal(0);
else
    nPad_pre = nlim_wav1(1)-1;
    nPad_post = size(mrWav_filt,1) - nlim_wav1(2);
end

% common mode rejection
nSites = size(P.miSites,2);
blank_thresh = get_set_(P, 'blank_thresh', 0);
if blank_thresh > 0 %% && nSites >= get_(P, 'nChans_min_car')
    if isempty(vrWav_mean_filt)
        vrWav_mean_filt = mean_excl_(mrWav_filt, P);
    end
    vlKeep_ref = car_reject_(vrWav_mean_filt(:), P); vrWav_mean_filt = [];
    fprintf('\tRejecting %0.3f %% of time due to motion\n', (1-mean(vlKeep_ref))*100 );
else
    vlKeep_ref = [];
end

% detect spikes or use the one passed from the input (importing)
if isempty(vrThresh_site), [vrThresh_site, fGpu] = mr2thresh_(mrWav_filt, P); end
[viTime_spk, vrAmp_spk, viSite_spk] = detect_spikes_(mrWav_filt, vrThresh_site, vlKeep_ref, P);
[viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)gather_(x), viTime_spk, vrAmp_spk, viSite_spk);    

% reject spikes within the overlap region
if ~isempty(nlim_wav1)
    viKeep_spk = find(viTime_spk >= nlim_wav1(1) & viTime_spk <= nlim_wav1(2));
    [viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)x(viKeep_spk), viTime_spk, vrAmp_spk, viSite_spk);    
end%if

% extract spike waveforms
trWav_spk = get_spkwav_(mrWav_filt, viSite_spk, viTime_spk, P);
mrVp_spk = [];

% extract spike feaures
if isempty(mrPv_global)
    [mrPv_global, vrD_global] = get_prinvec_(trWav_spk, P); 
end
trPc_spk = gather_(project_pc_(trWav_spk, mrPv_global, P));

if get_set_(P, 'sort_mode', 1) == 1 && size(trWav_spk,2) > 1
    viSite2_spk = find_site_spk23_(trWav_spk, viSite_spk, P);
    trWav_spk = []; %clear memory
    trWav2_spk = mn2tn_wav_spk2_(mrWav_filt, viSite2_spk, viTime_spk, P);
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
function [mrPv1, vrD1] = get_prinvec_(tr, P)
%tr: nSamples x nSpikes x nChans
% MAX_SAMPLE = 10000;        

% viSpk_sub = subsample_vr_(1:size(tr,3), MAX_SAMPLE);
% switch 4
%     case 1, mr1 = squeeze(tr(:, 1, viSpk_sub)); % use peak chan only
%     case 2, mr1 = reshape(tr(:, 1:P.nSites_fet, viSpk_sub), size(tr,1), []); 
%     case 3, mr1 = reshape(tr(:, :, viSpk_sub), size(tr,1), []); % use all chan
%     case 4, mr1 = reshape(tr, size(tr,1), []); % use all chan
% end
t_fun = tic;
nPc_spk = get_set_(P, 'nPc_spk', 9); % # components to compress spike waveforms
mr1 = gather_(reshape(tr, size(tr,1), []));
switch 3
    case 1
        % mrSpkWav1 = meanSubt_(mrSpkWav1);
        [mrPv1, vrD1] = eig(mr1 * mr1');
        mrPv1 = zscore_(fliplr(mrPv1)); % sort largest first
        vrD1 = flipud(diag(vrD1));
    case 2
        [mrPv1, ~, vrD1] = pca(mr1','Centered',false, 'NumComponents', nPc_spk);
    case 3
        [mrPv1, ~, vrD1] = pca(double(mr1'),'Centered',false, 'NumComponents', nPc_spk);
        [mrPv1, vrD1] = deal(single(mrPv1), single(vrD1));
end
fprintf('\tget_prinvec_: took %0.1fs\n', toc(t_fun));
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
function [mrWav_filt, vrWav_filt_mean] = filter_transpose_(mnWav_T, P)
%-----
% Filter
% fprintf('\tFiltering spikes (%s)...', get_(P, 'vcCommonRef')); t_fun = tic;
% if get_set_(P, 'fSmooth_spatial', 0)
%     mnWav_T = spatial_smooth_(mnWav_T, P);
% end
mrWav_filt = fft_filter_transpose(single(mnWav_T), P);

nChans = size(mrWav_filt, 2);
if nChans < get_set_(P, 'nChans_min_car', 32)
    P.vcCommonRef = 'none';
end
switch get_(P, 'vcCommonRef')
    case 'mean'
        vrWav_filt_mean = mean_excl_(mrWav_filt, P);                
    case 'median'
        vrWav_filt_mean = median_excl_(mrWav_filt, P);
    case {'trimmean', 'tmean'}
        vrWav_filt_mean = trimmean_excl_(mrWav_filt, P);   
    otherwise
        vrWav_filt_mean = [];
end
if ~isempty(vrWav_filt_mean)
    mrWav_filt = mrWav_filt - vrWav_filt_mean(:);
end
if get_set_(P, 'fWhiten', 0)
    fprintf('\tWhitening...');
    mrWav_filt = spatial_whiten_(mrWav_filt, P);
end
% fprintf(' took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
function mrB = spatial_whiten_(mrA, P)

MAX_SAMPLE = [];

nSites = size(mrA,2);
nSites_whiten = min(get_set_(P, 'nSites_whiten', 32), nSites);
if nSites <= nSites_whiten || nSites_whiten == 0
    % whiten globally using all sites
    mrB = mrA * whiten_matrix_(mrA, MAX_SAMPLE);
else
    % whiten locally, slower than the global method
    miSites = findNearestSites_(P.mrSiteXY, nSites_whiten);
    mrB = zeros(size(mrA), 'like', mrA);
    for iSite = 1:nSites
        mrA1 = mrA(:,miSites(:,iSite));
        mrW1 = whiten_matrix_(mrA1, MAX_SAMPLE);
        mrB(:,iSite) = mrA1 * mrW1(:,1);
    end
end
end %func


%--------------------------------------------------------------------------
function mrW = whiten_matrix_(mrA, nSamples)
% https://www.mathworks.com/matlabcentral/fileexchange/34471-data-matrix-whitening?s_tid=answers_rc2-2_p5_MLT
if nargin<2, nSamples=[]; end
if ~isempty(nSamples)
    mrA = subsample_row_(mrA, nSamples);
end
% mrA = mrA - mean(mrA,1); % subtract the channel mean
[mrE1, vrD1] = svd((mrA' * mrA)/size(mrA,1));
mrW = mrE1 * diag(1./sqrt(diag(vrD1) + eps)) * mrE1';    
end %func


%--------------------------------------------------------------------------
function miSites = findNearestSites_(mrSiteXY, nSites_spk)

nSites = size(mrSiteXY,1);
miSites = zeros(nSites_spk, nSites);
for iSite=1:nSites
    vrSiteDist = pdist2_(mrSiteXY(iSite,:), mrSiteXY);
    [vrSiteDist, viSrt] = sort(vrSiteDist, 'ascend');
    miSites(:,iSite) = viSrt(1:nSites_spk);
end
end %func


%--------------------------------------------------------------------------
function mr12 = pdist2_(mr1, mr2)
% mr12 = pdist2_(mr1) % self distance
% mr12 = pdist2_(mr1, mr2)
% mr1: n1xd, mr2: n2xd, mr12: n1xn2

% mr12 = sqrt(eucl2_dist_(mr1', mr2'));
% 20% faster than pdist2 for 10000x10 x 100000x10 single
if nargin==2
    mr12 = sqrt(bsxfun(@plus, sum(mr2'.^2), bsxfun(@minus, sum(mr1'.^2)', 2*mr1*mr2')));
else
    vr1 = sum(mr1'.^2);
    mr12 = sqrt(bsxfun(@plus, vr1, bsxfun(@minus, vr1', 2*mr1*mr1')));
end
end %func


%--------------------------------------------------------------------------
function [mrA, viCol] = subsample_row_(mrA, nSamples)
nRows = size(mrA,1);
if nRows < nSamples, return; end
viCol = round(linspace(1, nRows, nSamples));
mrA = mrA(viCol,:);
end %fun


%--------------------------------------------------------------------------
function vnWav1_mean = mean_excl_(mnWav1, P)
% calculate mean after excluding viSiteZero
viSiteZero = get_(P, 'viSiteZero');
if isempty(viSiteZero)
    vnWav1_mean = (mean(mnWav1,2));
else
    nSites_all = size(mnWav1, 2);
    nSites_excl = numel(viSiteZero);
    nSites = nSites_all - nSites_excl;
    vnWav1_mean = ((sum(mnWav1,2) - sum(mnWav1(:,nSites_excl),2)) / nSites);
end
end %func


%--------------------------------------------------------------------------
% trimmed mean
function vnWav1_mean = trimmean_excl_(mnWav1, P)
% calculate mean after excluding viSiteZero
viSiteZero = get_(P, 'viSiteZero');
if ~isempty(viSiteZero)
    mnWav1(:,viSiteZero) = [];
end
trimmean_pct = get_set_(P, 'trimmean_pct', 25);
vnWav1_mean = trimmean(mnWav1,trimmean_pct,2);
end %func


%--------------------------------------------------------------------------
function vnWav1_med = median_excl_(mnWav1, P)
% calculate mean after excluding viSiteZero
viSiteZero = get_(P, 'viSiteZero');
fGpu = isGpu_(mnWav1);
if fGpu, mnWav1 = gather_(mnWav1); end
if isempty(viSiteZero)
    vnWav1_med = median(mnWav1,2);
else
    viSites = setdiff(1:size(mnWav1,2), viSiteZero);
    vnWav1_med = median(mnWav1(:,viSites),2);
end
vnWav1_med = gpuArray_(vnWav1_med, fGpu);
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
    case 'boyden', vcDir_in = 'groundtruth/paired_recordings/boyden32c/915_10_1';
    case 'boyden2', vcDir_in = 'groundtruth/paired_recordings/boyden32c/509_1_1';
    otherwise, error('unsupported test mode');
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
function P = makeParam_(vcDir_in, vcDir_out, vcFile_arg, fParfor)
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end
if nargin<4, fParfor = ''; end

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
if ~isempty(fParfor)
    P.fParfor = fParfor; % override fParfor
end

% derived fields
P = fill_param_(P);
P.vcFile_prm = fullfile(vcDir_out, 'raw_geom.prm');

vcFile_gt = fullfile(vcDir_in, 'firings_true.mda');
if exist_file_(vcFile_gt)
    P.vcFile_gt = vcFile_gt;
end
P.vcVersion = version_();

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
fid_r = fopen(fname,'r','l');

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

dimm = fread(fid_r, num_dims, dim_type_str)';
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
elseif strcmpi(mask_out_artifacts, 'false')
    P.blank_thresh = 0;
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
step_sec_drift = get_set_(P, 'step_sec_drift', 20);
batch_sec_drift = get_set_(P, 'batch_sec_drift', 300);
if isempty(get_(P, 'nTime_batch'))
    P.nTime_batch = max(round(batch_sec_drift / step_sec_drift), 1);
    fprintf('\tnTime_batch = %d (batch_sec_drift = %0.1f s)\n', P.nTime_batch, batch_sec_drift);
end
if isempty(get_(P, 'nTime_drift'))
    try
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
S_cfg = read_cfg_();
[P.nC_max, P.nTime_max_drift] = get_(S_cfg, 'nC_max', 'nTime_max_drift'); % override nC_max (gpu parameter)
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
fid_w = vcFile_out;
for iFile = 1:numel(csFiles)
    [~,~,vcExt1] = fileparts(csFiles{iFile});
    switch lower(vcExt1)
        case '.mda', mr_ = readmda_(csFiles{iFile});
%         case '.bin', mr_ = load_spikeglx_(csFiles{iFile});
        otherwise, error('join_mda_: unsupported format');
    end
    fid_w = writemda_fid(fid_w, mr_); mr_ = [];
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
function mnWav = load_bin_(vcFile, vcDataType, dimm, header, vnLoad)
% mnWav = load_bin_(vcFile, dimm, vcDataType)
% mnWav = load_bin_(fid, dimm, vcDataType)
% mnWav = load_bin_(vcFile_prm)
% [Input arguments]
% header: header bytes

if nargin<2, vcDataType = []; end
if nargin<3, dimm = []; end
if nargin<4, header = 0; end
if nargin<5, vnLoad = []; end
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
    if numel(header) == 1
        if header>0, fseek(fid, header, 'bof'); end
    end
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
    if isempty(vnLoad)
        mnWav = fread_(fid, dimm, vcDataType);
    else % combine multiple skips
        viOffset = header;
        assert(numel(vnLoad) == numel(viOffset), 'load_bin_: vnLoad size does not match viOffset');
        dimm_wav = dimm; dimm_wav(end) = sum(vnLoad);
        mnWav = zeros(dimm_wav, vcDataType);  
        byte_per_sample1 = prod(dimm_wav(1:end-1)) * bytesPerSample_(vcDataType);
        iT_offset1 = 0;
        for iLoad = 1:numel(vnLoad)
            vi1 = (1:vnLoad(iLoad)) + iT_offset1;
            if isempty(vi1), continue; end
            fseek(fid, byte_per_sample1 * viOffset(iLoad), 'bof');
            dimm_wav(end) = vnLoad(iLoad);
            switch numel(dimm_wav)
                case 1, mnWav(vi1) = fread(fid, dimm_wav, vcDataType);
                case 2, mnWav(:,vi1) = fread(fid, dimm_wav, vcDataType);
                case 3, mnWav(:,:,vi1) = fread_(fid, dimm_wav, vcDataType);
                case 4, mnWav(:,:,:,vi1) = fread_(fid, dimm_wav, vcDataType);
                otherwise, error('load_bin_: unsupported dimensions');
            end
            iT_offset1 = vi1(end);
        end
    end
    if ischar(vcFile), fclose(fid); end
catch
    disp(lasterr());
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

fParfor = get_set_(P, 'fParfor', 0);
fprintf('Calculating drift similarity...'); t1 = tic;
[nTime_batch, nTime_drift] = get_(P, 'nTime_batch', 'nTime_drift');
if isempty(nTime_batch)
    nTime_batch = round(get_(P, 'batch_sec_drift') / get_(P, 'step_sec_drift'));    
end
if isempty(get_(P, 'nTime_max_drift'))
    P.nTime_max_drift = read_cfg_('nTime_max_drift');
end
viDrift_spk = []; %ones(1, numel(S0.viSite_spk), 'int64');
if isempty(nTime_batch) || isempty(nTime_drift) || nTime_batch >= nTime_drift
    % nDrift = 1;
    viLim_drift = [1, numel(S0.viSite_spk)+1];
    [mlDrift, nTime_drift] = deal(true, 1);    
else
    nAmp_drift = get_set_(P, 'nQuantile_drift', 10);
    nPos_drift = numel(P.viSite2Chan); % use number of sites
    nSpikes = numel(S0.viSite_spk);
    viLim_drift = int64([0, ceil((1:nTime_drift)/nTime_drift*nSpikes)])+1;

    % collect stats
    mrPos_spk = S0.mrPos_spk;
    vrAmp_spk = single(S0.vrAmp_spk(:));
    vrAmp_quantile = quantile_vr_(vrAmp_spk, (0:nAmp_drift)/nAmp_drift);
    vrPos_spk = single(mrPos_spk(:,2));
    mrCount_drift = zeros(nAmp_drift*nPos_drift, nTime_drift, 'single');
    vrPos_quantile = quantile_vr_(vrPos_spk, (0:nPos_drift)/nPos_drift);
    
    if fParfor
        try
            cvrAmp_spk_drift = arrayfun_(@(a,b)vrAmp_spk(a:b), viLim_drift(1:end-1), viLim_drift(2:end)-1);
            cvrPos_spk_drift = arrayfun_(@(a,b)vrPos_spk(a:b), viLim_drift(1:end-1), viLim_drift(2:end)-1);
            parfor iDrift = 1:nTime_drift
                mn_ = histcounts2(cvrAmp_spk_drift{iDrift}, cvrPos_spk_drift{iDrift}, vrAmp_quantile, vrPos_quantile);
                mrCount_drift(:,iDrift) = mn_(:);
            end
        catch
            fParfor = 0;
        end
    end
    if ~fParfor
        for iDrift = 1:nTime_drift
            viSpk1 = viLim_drift(iDrift):viLim_drift(iDrift+1)-1;
            mn_ = histcounts2(vrAmp_spk(viSpk1), vrPos_spk(viSpk1), vrAmp_quantile, vrPos_quantile);
            mrCount_drift(:,iDrift) = mn_(:);
        end
    end
    
    if nTime_drift > P.nTime_max_drift
        % large recording 
        quarter_ = P.nTime_max_drift/4;
        mlMask = abs((1:nTime_drift) - (1:nTime_drift)') <= quarter_;
        mlMask(1:quarter_*2+1,1:quarter_) = true;
        mlMask(end-quarter_*2:end,end-quarter_+1:end) = true;
        mrDist_drift = dist_mask_(mrCount_drift, mlMask);
    else
        mrDist_drift = squareform(pdist(mrCount_drift'));
    end
    [mrSort_drift, miSort_drift] = sort(mrDist_drift, 'ascend');
    miSort_drift = miSort_drift(1:nTime_batch,:);
    
    if read_cfg_('fPlot_drift')
        figure; imagesc(mrDist_drift); set(gcf,'Name', P.vcFile_prm);
        figure; imagesc(mrSort_drift); set(gcf,'Name', P.vcFile_prm);
        hold on; plot([0, size(mrSort_drift,1)], repmat(nTime_batch,1,2), 'r-');
    end
    mlDrift = mi2ml_drift_(miSort_drift); %gpuArray_(mi2ml_drift_(miSort_drift), P.fGpu);
end
S_drift = makeStruct_(nTime_drift, viDrift_spk, mlDrift, viLim_drift);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function mrD = dist_mask_(mrF, mlMask)
% return euclidean distance for masked portions

mrD = nan(size(mrF,2), 'like', mrF);
for iCol=1:size(mrF,2)
    viRow1 = find(mlMask(:,iCol));
    mrD(viRow1,iCol) = sqrt(sum((mrF(:,iCol) -  mrF(:,viRow1)).^2, 1));
end
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
function S_auto = postCluster_(S_clu, P, viSite_spk)
t_fun = tic;
if isempty(S_clu), S_auto=[]; return; end

fprintf('\tpostCluster_...\n');
P.min_count = get_set_(P, 'min_count', 0);

[viClu, icl, nClu_pre] = assignCluster_knn_(S_clu, viSite_spk, P);
nClu_post = numel(icl);
[viClu, icl] = dpclus_remove_count_(viClu, icl, P.min_count);

S_auto = makeStruct_(viClu, icl, P);
S_auto = S_auto_refresh_(S_auto, 0, []);
fprintf('\tpostCluster_: Pre-merged %d->%d->%d clusters, took %0.1fs\n', ...
    nClu_pre, nClu_post, S_auto.nClu, toc(t_fun));
end %func


%--------------------------------------------------------------------------
function [S_auto, vlKeep_clu] = S_auto_refresh_(S_auto, fRemoveEmpty, viSite_spk)

nClu = double(max(S_auto.viClu));
S_auto.nClu = nClu;
S_auto.cviSpk_clu = vi2cell_(S_auto.viClu, nClu);
S_auto.vnSpk_clu = cellfun(@numel, S_auto.cviSpk_clu); 
S_auto.nSpk_unique = sum(cellfun(@numel, S_auto.cviSpk_clu));
if ~isempty(viSite_spk)
    S_auto.viSite_clu = double(arrayfun(@(iClu)mode(viSite_spk(S_auto.cviSpk_clu{iClu})), 1:nClu));
    S_auto.viSite_clu = S_auto.viSite_clu(:);
end
if fRemoveEmpty
    [S_auto, vlKeep_clu] = S_clu_remove_empty_(S_auto); 
end
end %func


%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu)
vlKeep_clu = S_clu.vnSpk_clu>0;
if all(vlKeep_clu), return; end

% waveform
S_clu = S_clu_select_(S_clu, vlKeep_clu);
if min(S_clu.viClu) < 1
    S_clu.viClu(S_clu.viClu<1) = 0;
    [~,~,S_clu.viClu] = unique(S_clu.viClu+1);        
    S_clu.viClu = S_clu.viClu-1;
else
    [~,~,S_clu.viClu] = unique(S_clu.viClu);        
end
S_clu.viClu = int32(S_clu.viClu);
S_clu.nClu = double(max(S_clu.viClu));
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_select_(S_clu, viKeep_clu)
% automatically trim clusters
% 7/20/17 JJJ: auto selecting vectors and matrics
% excl vnSpk_clu, viSite_clu, vrPosX_clu, vrPosY_clu

% Quality
csNames = fieldnames(S_clu);
if isempty(csNames), return; end
viMatch_v = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^v\w*_clu$'), csNames, 'UniformOutput', false));
S_clu = struct_select_(S_clu, csNames(viMatch_v), viKeep_clu);

viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), csNames, 'UniformOutput', false));
S_clu = struct_select_(S_clu, csNames(viMatch_t), viKeep_clu, 3); 

viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), csNames, 'UniformOutput', false));
S_clu = struct_select_(S_clu, csNames(viMatch_c), viKeep_clu);
end %func


%--------------------------------------------------------------------------
function S = struct_select_(S, csNames, viKeep, iDimm)
if isempty(csNames), return; end
if nargin<4, iDimm = 1; end

% function test
if ischar(csNames), csNames = {csNames}; end
for i=1:numel(csNames)
    vcName_ = csNames{i};
    if ~isfield(S, vcName_), continue; end
    try
        val = S.(vcName_);
        if isempty(val), continue; end
        ndims_ = ndims(val);
        if ndims_==2 %find a column or row vectors
            if size(val,1)==1 || size(val,2)==1, ndims_=1; end %iscol or isrow
        end
        switch ndims_
            case 1, val = val(viKeep);
            case 2
                switch iDimm
                    case 1, val = val(viKeep,:);
                    case 2, val = val(:,viKeep);
                    case 3
                        val = val(:,:,viKeep);
                    otherwise,disperr_('struct_select_: invalid iDimm');
                end
            case 3
                switch iDimm
                    case 1, val = val(viKeep,:,:);
                    case 2, val = val(:,viKeep,:);
                    case 3, val = val(:,:,viKeep);
                    otherwise, disperr_('struct_select_: invalid iDimm');
                end
            otherwise, disperr_('struct_select_: invalid # of dimensions (1-3 supported)');
        end %switch
        S.(vcName_) = val;
    catch
        disperr_(sprintf('struct_select_: %s field error', vcName_));
    end
end %for
end %func


%--------------------------------------------------------------------------
function [viClu_spk, viSpk_peak, nClu_pre] = assignCluster_knn_(S_clu, viSite_spk, P)

nRepeat = 10;
[vrRho_spk, vrDelta_spk, ordrho_spk, viSpk_nneigh] = ...
    get_(S_clu, 'rho', 'delta', 'ordrho', 'nneigh');
viSpk_peak = find(vrDelta_spk(:) >= get_set_(P, 'delta_cut', 1));
nClu_pre = numel(viSpk_peak);

fParfor = get_set_(P, 'fParfor', 1);
nSpk = numel(ordrho_spk);
viClu_spk = zeros([nSpk, 1], 'int32');
if false
    viSpk_peak(vrDelta(viSpk_peak) > 1e10) = []; % remove super high 
end
[vrRho_peak, ix_] = sort(vrRho_spk(viSpk_peak), 'descend'); viSpk_peak = viSpk_peak(ix_);
miKnn_peak = load_miKnn_spk_(P, viSite_spk, viSpk_peak);    

% check for knn overlap    
cvi_peak = cell(numel(viSpk_peak), 1);  
if fParfor
    try
        parfor iPeak = 1:numel(viSpk_peak)
            viPeak1 = find(any(ismember(miKnn_peak(:,1:iPeak-1), miKnn_peak(:,iPeak))));        
            cvi_peak{iPeak} = [viPeak1(:); iPeak];
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iPeak = 1:numel(viSpk_peak)
        viPeak1 = find(any(ismember(miKnn_peak(:,1:iPeak-1), miKnn_peak(:,iPeak))));        
        cvi_peak{iPeak} = [viPeak1(:); iPeak];
    end
end
[viClu_spk(viSpk_peak), viiPeak] = cell2map_(cvi_peak);
[viSpk_peak, vrRho_peak] = deal(viSpk_peak(viiPeak), vrRho_peak(viiPeak));

if numel(viSpk_peak) == 0 || numel(viSpk_peak) == 1
    viClu_spk = ones([nSpk, 1], 'int32');
    viSpk_peak = ordrho_spk(1);
else
    % remove icl that is too close    
    nneigh1 = viSpk_nneigh(ordrho_spk);
    for iRepeat=1:nRepeat
        vi = find(viClu_spk(ordrho_spk)<=0);
        if isempty(vi), break; end
        vi=vi(:)';        
        for ii = vi
            viClu_spk(ordrho_spk(ii)) = viClu_spk(nneigh1(ii));        
        end
        n1 = sum(viClu_spk<=0);
        fprintf('\tnRepeat:%d, n0=%d\n', iRepeat, n1);
        if n1==0, break; end
    end
end
end %func


%--------------------------------------------------------------------------
function [viClu_new, icl_new] = dpclus_remove_count_(viClu, icl, min_count)
nClu = numel(icl);
viMap = zeros(nClu,1);
vnSpk_clu = cellfun(@numel, vi2cell_(viClu, nClu));
% vnSpk_clu = arrayfun(@(x)sum(viClu==x), 1:nClu);
vlClu_keep = vnSpk_clu >= min_count;
if all(vlClu_keep)
    viClu_new = viClu;
    icl_new = icl;
else
    viMap(vlClu_keep) = 1:sum(vlClu_keep);
    viClu_new = map_index_(viMap, viClu);
    icl_new = icl(vlClu_keep);
end
end %func


%--------------------------------------------------------------------------
% 190816 JJJ: Faster implementation
function [viMapClu_new, viUniq_, viMapClu] = cell2map_(cvi_clu)
nRepeat = 10;

t_fun = tic;
for iRepeat = 1:nRepeat
    cvi_clu1 = cell(size(cvi_clu));
    viClu_update = find(cellfun(@(x)numel(x)>1, cvi_clu));
    viClu_update = viClu_update(:)';
    for iClu = viClu_update
        viClu1 = cvi_clu{iClu};
        viClu1 = cat(1, viClu1(:), cvi_clu{viClu1});
        cvi_clu1(viClu1) = {viClu1};
    end    
    cvi_clu1 = cellfun_(@(x)unique(x), cvi_clu1);
    if all(cellfun(@numel, cvi_clu1) == cellfun(@numel, cvi_clu))
        break;
    else
        cvi_clu = cvi_clu1;
    end
end

nClu = numel(cvi_clu);
viMapClu = 1:nClu;
for iClu = 1:nClu
    iClu1 = min(cvi_clu{iClu});
    if ~isempty(iClu1), viMapClu(iClu) = iClu1; end
end

% Compact the map so the index doesn't have a gap
viUniq_ = unique(viMapClu);
viMap_(viUniq_) = 1:numel(viUniq_);
viMapClu_new = viMap_(viMapClu);
fprintf('\tcell2map_: nRepeat=%d, took %0.1fs\n', iRepeat, toc(t_fun));
end %func


%--------------------------------------------------------------------------
% do not map self it not found
function [viMapClu_new, viUniq_, viMapClu] = ml2map_(mlClu)
nClu = size(mlClu,1);
% mlClu = set_diag_(mlClu | mlClu', true(nClu,1));
% set diagonal to true

viClu = find(any(mlClu,1) | any(mlClu,2)');
for iClu = viClu
    vi1_ = find(mlClu(:,iClu) | mlClu(iClu,:)');
    mlClu(vi1_,vi1_) = true;
end    

viMapClu = zeros(1, nClu);
for iClu = 1:nClu
    iClu1 = find(mlClu(:,iClu), 1, 'first');
    if ~isempty(iClu1), viMapClu(iClu) = iClu1; end
end
% viMapClu = viMapClu(viMapClu>0);

% Compact the map so the index doesn't have a gap
viUniq_ = setdiff(unique(viMapClu), 0);
viMap_(viUniq_) = 1:numel(viUniq_);
viMapClu_new = zeros(size(viMapClu));
vlPos = viMapClu>0;
viMapClu_new(vlPos) = viMap_(viMapClu(vlPos));
end %func


%--------------------------------------------------------------------------
% 190816 JJJ: Faster implementation
function [viMapClu_new, viUniq_, viMapClu] = ml2map__(mlClu)
nClu = size(mlClu,1);
% mlClu = set_diag_(mlClu | mlClu', true(nClu,1));
% set diagonal to true

viClu = find(any(mlClu,1) | any(mlClu,2)');
for iClu = viClu
    vi1_ = find(mlClu(:,iClu) | mlClu(iClu,:)');
    mlClu(vi1_,vi1_) = true;
end    

viMapClu = 1:nClu;
for iClu = 1:nClu
    iClu1 = find(mlClu(:,iClu), 1, 'first');
    if ~isempty(iClu1), viMapClu(iClu) = iClu1; end
end

% Compact the map so the index doesn't have a gap
viUniq_ = unique(viMapClu);
viMap_(viUniq_) = 1:numel(viUniq_);
viMapClu_new = viMap_(viMapClu);
end %func


%--------------------------------------------------------------------------
% 10/12/17 JJJ: Works for non-square matrix and constant. Tested
function mr = set_diag_(mr, vr)
n = min(min(size(mr)), numel(vr));
mr(sub2ind(size(mr), 1:n, 1:n)) = vr(1:n);
end %func


%--------------------------------------------------------------------------
function [dimm_mr, type_mr] = write_bin_(vcFile, mr, fVerbose)
if nargin<3, fVerbose = []; end

t_fun = tic;
if isa(mr, 'gpuArray'), mr = gather(mr); end
dimm_mr = size(mr);
type_mr = class(mr);
if isempty(mr), return; end
if isstruct(mr)
    save(vcFile, '-struct', 'mr', '-v7.3'); % save to matlab file
else
    if ischar(vcFile)
        fid_w = fopen(vcFile, 'W'); 
    else
        fid_w = vcFile;
    end
    if iscell(mr)
        cmr = mr;
        ndimm_end = 0;
        for iCell = 1:numel(cmr)
            mr1 = cmr{iCell};
            fwrite(fid_w, mr1, class(mr1));
            dimm1 = size(mr1);
            ndimm_end = ndimm_end + dimm1(end);
        end
        dimm_mr = size(cmr{1});
        dimm_mr(end) = ndimm_end;
    else
        fwrite(fid_w, mr, type_mr);
    end
    if ischar(vcFile)
        fclose(fid_w); 
    else
        fVerbose = 0;
    end
end
if fVerbose
    fprintf('Writing to %s took %0.1fs\n', vcFile, toc(t_fun));
end
end %func


%--------------------------------------------------------------------------
% 17/9/13 JJJ: Behavior changed, if S==[], S0 is loaded
function val = get_set_(S, vcName, def_val)
% set a value if field does not exist (empty)

if isempty(S), S = get(0, 'UserData'); end
if isempty(S), val = def_val; return; end
if ~isstruct(S)
    val = []; 
    fprintf(2, 'get_set_: %s must be a struct\n', inputname(1));
    return;
end
val = get_(S, vcName);
if isempty(val), val = def_val; end
end %func


%--------------------------------------------------------------------------
% 11/19/2018 JJJ: improved matlab version check
% 7/13/17 JJJ: Version check routine
function struct_save_(S, vcFile, fVerbose)
nRetry = 3;
if nargin<3, fVerbose = 0; end
if fVerbose
    fprintf('Saving a struct to %s...\n', vcFile); t1=tic;
end
if version_matlab_() >= 2017
    for iRetry=1:nRetry
        try
            save(vcFile, '-struct', 'S', '-v7.3', '-nocompression'); %faster    
            break;
        catch
            pause(.5);
        end
        fprintf(2, 'Saving failed: %s\n', vcFile);
    end
else    
    for iRetry=1:nRetry
        try
            save(vcFile, '-struct', 'S', '-v7.3');   
            break;
        catch
            pause(.5);
        end
        fprintf(2, 'Saving failed: %s\n', vcFile);
    end    
end
if fVerbose
    fprintf('\ttook %0.1fs.\n', toc(t1));
end
end %func


%--------------------------------------------------------------------------
% 11/19/2018 JJJ: returns 2017 if R2017a and 2017.5 if R2017b
function ver_year = version_matlab_()

vcVer_matlab = version('-release');
ver_year = str2double(vcVer_matlab(1:end-1));
if lower(vcVer_matlab(end)) == 'b', ver_year = ver_year + .5; end
end %func


%--------------------------------------------------------------------------
% 01/08/2020 JJJ: compile, mcc, and copy to the container locatino
function export_sf2_()
irc2('compile');
irc2('mcc');
copyfile('run_irc', read_cfg_('spikeforest2_irc_path'));
end %func


%--------------------------------------------------------------------------
% 4/17/18 JJJ: documentation and testing
function varargout = call_(vcFunc, cell_Input)
% S_out = call_(vcFunc, cell_Input, nOutput)
% varargout = call_(vcFunc, cell_Input)

if vcFunc(end) ~= '_', vcFunc = [vcFunc, '_']; end

nOutput = nargout();
switch nOutput
    case 0, feval(vcFunc, cell_Input{:});
    case 1, varargout{1} = feval(vcFunc, cell_Input{:});
    case 2, [varargout{1}, varargout{2}] = feval(vcFunc, cell_Input{:});
    case 3, [varargout{1}, varargout{2}, varargout{3}] = feval(vcFunc, cell_Input{:});
    case 4, [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = feval(vcFunc, cell_Input{:});
end %switch
end %func


%--------------------------------------------------------------------------
% 17/12/1 JJJ: Load size is not limited by FFT cleanup process (fft_thresh>0)
function [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file, P)
% plan load file size according to the available memory and file size (nBytes_file1)
% LOAD_FACTOR = 8; % Memory usage overhead
vcDataType = get_set_(P, 'vcDataType_filter', P.vcDataType);
nSamples1 = floor(nBytes_file / bytesPerSample_(P.vcDataType) / P.nChans);
MAX_BYTES_LOAD = get_set_(P, 'MAX_BYTES_LOAD', .5e9);
nSamples_max = floor(MAX_BYTES_LOAD / P.nChans / bytesPerSample_(vcDataType));

fTranspose_bin = get_set_(P, 'fTranspose_bin', 1); %load all in one, Catalin's format
if ~fTranspose_bin 
    [nLoad1, nSamples_load1, nSamples_last1] = deal(1, nSamples1, nSamples1);
else
    [nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples1, nSamples_max);
end
end %func


%--------------------------------------------------------------------------
function [nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples1, nSamples_max)
nLoad1 = setlim_(ceil(nSamples1 / nSamples_max), [1, inf]); 
nSamples_load1 = min(nSamples1, nSamples_max);
if nLoad1 == 1
    nSamples_load1 = nSamples1;
    nSamples_last1 = nSamples1;
else
    nSamples_last1 = mod(nSamples1, nSamples_load1);
    if nSamples_last1==0
        nSamples_last1 = nSamples_load1;
    elseif nSamples_last1 < nSamples_load1/2
        % if last part is too small increase the size
        nLoad1 = nLoad1 - 1;
        if nLoad1==1
            nSamples_load1 = nSamples1;
            nSamples_last1 = nSamples1;
        else
            nSamples_last1 = nSamples_last1 + nSamples_load1;
        end
    end
    
end
end %func


%--------------------------------------------------------------------------
% 17/9/13 JJJ: Created and tested
function vr = setlim_(vr, lim_)
% Set low and high limits
vr = min(max(vr, lim_(1)), lim_(2));
end %func


%--------------------------------------------------------------------------
function [viTime_spk, vrAmp_spk, viSite_spk] = detect_spikes_(mnWav3, vnThresh_site, vlKeep_ref, P)
% fMerge_spk = 1;
fMerge_spk = get_set_(P, 'fMerge_spk', 1);

[n1, nSites, ~] = size(mnWav3);
[cviSpk_site, cvrSpk_site] = deal(cell(nSites,1));
if isempty(vnThresh_site)    
    vnThresh_site = mr2thresh_(mnWav3, P);
end

for iSite = 1:nSites
    [viSpk11, vrSpk11] = spikeDetectSingle_fast_(mnWav3(:,iSite), P, vnThresh_site(iSite));
    
    % Reject global mean
    if isempty(vlKeep_ref)
        cviSpk_site{iSite} = viSpk11;
        cvrSpk_site{iSite} = vrSpk11;        
    else
        [cviSpk_site{iSite}, cvrSpk_site{iSite}] = select_vr_(viSpk11, vrSpk11, find(vlKeep_ref(viSpk11)));
    end
end

% Group spiking events using vrWav_mean1. already sorted by time
if fMerge_spk
    [viTime_spk, vrAmp_spk, viSite_spk] = spikeMerge_(cviSpk_site, cvrSpk_site, P);
else
    viTime_spk = cell2mat_(cviSpk_site);
    vrAmp_spk = cell2mat_(cvrSpk_site);
    viSite_spk = cell2vi_(cviSpk_site);
    %sort by time
    [viTime_spk, viSrt] = sort(viTime_spk, 'ascend');
    [vrAmp_spk, viSite_spk] = multifun_(@(x)x(viSrt), vrAmp_spk, viSite_spk);
end
vrAmp_spk = gather_(vrAmp_spk);

end %func


%--------------------------------------------------------------------------
function varargout = select_vr_(varargin)
% [var1, var2, ...] = select_vr(var1, var2, ..., index)

% sort ascend
viKeep = varargin{end};
if islogical(viKeep), viKeep = find(viKeep); end
for i=1:(nargin-1)
    if isvector(varargin{i})
        varargout{i} = varargin{i}(viKeep);
    else
        varargout{i} = varargin{i}(viKeep, :);
    end
end
end %func


%--------------------------------------------------------------------------
function vi = cell2vi_(cvi)
% convert cell index to array of index
vn_site = cellfun(@(x)numel(x), cvi); %create uniform output
vi = cell(numel(cvi), 1);
for iSite=1:numel(cvi)
    vi{iSite} = iSite * ones(vn_site(iSite), 1);
end
vi = cell2mat_(vi);
end %func


%--------------------------------------------------------------------------
function [viSpk, vrSpk, viSite] = spikeMerge_(cviSpk, cvrSpk, P)
% provide spike index (cviSpk) and amplitudes (cvrSPk) per sites

nSites = numel(cviSpk);
viSpk = cell2mat_(cviSpk);      vrSpk = cell2mat_(cvrSpk);
viSite = cell2mat_(cellfun(@(vi,i)repmat(i,size(vi)), cviSpk, num2cell((1:nSites)'), 'UniformOutput', false));
[viSpk, viSrt] = sort(viSpk);   vrSpk = vrSpk(viSrt);   viSite = viSite(viSrt);
viSite = int32(viSite); 
viSpk = int64(viSpk);    % deal with longer recording

[cviSpkA, cvrSpkA, cviSiteA] = deal(cell(nSites,1));
%fParfor = get_set_(P, 'fParfor', 1);
fParfor = 0;
if fParfor
    try
        parfor iSite = 1:nSites %parfor speedup: 2x %parfor
            try
                [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                    spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P);            
            catch
                fprintf(2, 'Parfor error, retrying using a single thread.\n');
            end
        end
    catch
        fParfor = 0; 
    end
end %if
if ~fParfor
    for iSite = 1:nSites
        try
            [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P);            
        catch
            disperr_();
        end
    end
end

% merge parfor output and sort
viSpk = cell2mat_(cviSpkA);
vrSpk = cell2mat_(cvrSpkA);
viSite = cell2mat_(cviSiteA);
[viSpk, viSrt] = sort(viSpk); %sort by time
vrSpk = gather_(vrSpk(viSrt));
viSite = viSite(viSrt);
end %func


%--------------------------------------------------------------------------
function [viTime_spk2, vnAmp_spk2, viSite_spk2] = spikeMerge_single_(viTime_spk, vnAmp_spk, viSite_spk, iSite1, P)
% maxDist_site_um = get_set_(P, 'maxDist_site_spk_um', 75);
maxDist_site_um = get_set_(P, 'maxDist_site_um', 50);
nlimit = int64(abs(P.spkRefrac));
% spkLim = [-nlimit, nlimit];

% Find spikes from site 1
viSpk1 = find(viSite_spk == iSite1); % pre-cache
[viTime_spk1, vnAmp_spk1] = deal(viTime_spk(viSpk1), vnAmp_spk(viSpk1));

% Find neighboring spikes
viSite1 = findNearSite_(P.mrSiteXY, iSite1, maxDist_site_um);
viSpk12 = int32(find(ismember(viSite_spk, viSite1)));

% Coarse selection
[viTime_spk12, vnAmp_spk12, viSite_spk12] = deal(viTime_spk(viSpk12), vnAmp_spk(viSpk12), viSite_spk(viSpk12));

% Fine selection
vlKeep_spk1 = true(size(viSpk1));
for iDelay = -nlimit:nlimit    
    [vi12_, vi1_] = ismember(viTime_spk12, viTime_spk1 + iDelay);    
    vi12_ = find(vi12_);
    if iDelay == 0 % remove self if zero delay
        vi12_(viSpk12(vi12_) == viSpk1(vi1_(vi12_))) = [];
    end
    vi12_(vnAmp_spk12(vi12_) > vnAmp_spk1(vi1_(vi12_))) = []; % keep more negative spikes
    vlAmpEq = vnAmp_spk12(vi12_) == vnAmp_spk1(vi1_(vi12_));
    if any(vlAmpEq)
        if iDelay > 0 % spk1 occurs before spk12, thus keep 
            vi12_(vlAmpEq) = [];
        elseif iDelay == 0 % keep only if site is lower
            vlAmpEq(iSite1 > viSite_spk12(vi12_(vlAmpEq))) = 0;
            vi12_(vlAmpEq) = []; %same site same time same ampl is not possible
        end
    end
    vlKeep_spk1(vi1_(vi12_)) = 0;
end %for

% Keep the peak spikes only
viiSpk1 = find(vlKeep_spk1); %speed up since used multiple times
[viTime_spk2, vnAmp_spk2] = deal(viTime_spk1(viiSpk1), vnAmp_spk1(viiSpk1));
viSite_spk2 = repmat(int32(iSite1), size(viiSpk1));
end %func


%--------------------------------------------------------------------------
% 11/15/17 JJJ: Cast the threshold like the vrWav1
function [viSpk1, vrSpk1, thresh1] = spikeDetectSingle_fast_(vrWav1, P, thresh1)
% P: spkThresh, qqSample, qqFactor, fGpu, uV_per_bit
% vrWav1 can be either single or int16
% 6/27/17 JJJ: bugfix: hard set threshold is applied

% Determine threshold
MAX_SAMPLE_QQ = 300000; 
% fSpikeRefrac_site = 0;
if nargin < 3, thresh1 = []; end
if nargin < 2, P = struct('spkThresh', [], 'qqFactor', 5); end
if ~isempty(get_(P, 'spkThresh')), thresh1 = P.spkThresh; end

if thresh1==0, [viSpk1, vrSpk1] = deal([]); return; end % bad site
if isempty(thresh1)    
    thresh1 = median(abs(subsample_vr_(vrWav1, MAX_SAMPLE_QQ)));
    thresh1 = single(thresh1)* P.qqFactor / 0.6745;
end
thresh1 = cast(thresh1, 'like', vrWav1); % JJJ 11/5/17

% detect valley turning point. cannot detect bipolar
% pick spikes crossing at least three samples
nneigh_min = get_set_(P, 'nneigh_min_detect', 0);
viSpk1 = find_peak_(vrWav1, thresh1, nneigh_min);
if get_set_(P, 'fDetectBipolar', 0)
   viSpk1 = [viSpk1; find_peak_(-vrWav1, thresh1, nneigh_min)]; 
   viSpk1 = sort(viSpk1);
end
if isempty(viSpk1)
    viSpk1 = double([]);
    vrSpk1 = int16([]);
else
    vrSpk1 = vrWav1(viSpk1);
    % Remove spikes too large
    spkThresh_max_uV = get_set_(P, 'spkThresh_max_uV', []);
    if ~isempty(spkThresh_max_uV)
        thresh_max1 = abs(spkThresh_max_uV) / get_set_(P, 'uV_per_bit', 1);
        thresh_max1 = cast(thresh_max1, 'like', vrSpk1);
        viA1 = find(abs(vrSpk1) < abs(thresh_max1));
        viSpk1 = viSpk1(viA1);
        vrSpk1 = vrSpk1(viA1); 
    end        
end

if isGpu_(viSpk1)
    [viSpk1, vrSpk1, thresh1] = multifun_(@gather, viSpk1, vrSpk1, thresh1);
end
end %func


%--------------------------------------------------------------------------
% 9/13/17 JJJ: Edge case error fixed (vi2 indexing) for small block many loads
% 8/17/17 JJJ: fix for missing spikes
function viSpk1 = find_peak_(vrWav1, thresh1, nneigh_min)
% nneigh_min: number of neighbors around the spike below the threshold
%  0,1,2. # neighbors of minimum point below negative threshold 
% thresh1: absolute value. searching for negative peaks only

if nargin<3, nneigh_min = []; end
if isempty(nneigh_min), nneigh_min = 1; end

viSpk1 = [];
if isempty(vrWav1), return; end
vl1 = vrWav1 < -abs(thresh1);
vi2 = find(vl1);
%vi2 = find(vrWav1 < -thresh1);
if isempty(vi2), return; end

if vi2(1)<=1
    if numel(vi2) == 1, return; end
    vi2(1) = []; 
end    
if vi2(end)>=numel(vrWav1)
    if numel(vi2) == 1, return; end
    vi2(end) = []; 
end
vrWav12 = vrWav1(vi2);
viSpk1 = vi2(vrWav12 <= vrWav1(vi2+1) & vrWav12 <= vrWav1(vi2-1));
if isempty(viSpk1), return; end

switch nneigh_min
    case 1
        viSpk1 = viSpk1(vl1(viSpk1-1) | vl1(viSpk1+1));
    case 2
        viSpk1 = viSpk1(vl1(viSpk1-1) & vl1(viSpk1+1));
end
end %func


%--------------------------------------------------------------------------
function viSiteNear = findNearSite_(mrSiteXY, iSite, maxDist_site_um)
vrDist = pdist2_(mrSiteXY(iSite,:), mrSiteXY);
viSiteNear = find(vrDist <= maxDist_site_um);
end %func


%--------------------------------------------------------------------------
% Run spikeforest2 through python docker interface
function run_spikeforest2_(vcSorter, vcDir_in, vcDir_out, vcArg)

csSorters_sf2 = {'mountainsort4', 'ironclust', 'kilosort2', 'kilosort', 'spykingcircus', 'herdingspikes2', 'tridesclous', 'klusta', 'waveclus', 'jrclust'};

t_fun=tic;
if isempty(vcDir_out)
    vcDir_out = fullfile(vcDir_in, vcSorter);
    mkdir(vcDir_out);
end
if contains(lower(vcSorter), csSorters_sf2)
    % check python path
    spikeforest2_python_path = read_cfg_('spikeforest2_python_path');
    [~, loaded_python_path, ~] = pyversion();    
    if ~strcmpi(spikeforest2_python_path, loaded_python_path)
        try
            pyversion(spikeforest2_python_path)
        catch
            fprintf(2, 'Using %s\n', loaded_python_path);
        end
    end
    
    fprintf('Running %s through spikeforest2...\n', vcSorter);
    vcFile_firings = fullfile(vcDir_out, 'firings.mda');
    py.spikeforest2.experimental.sort(lower(vcSorter), vcDir_in, vcFile_firings);
    fprintf('%s: wrote to %s, took %0.1fs\n', vcSorter, vcFile_firings, toc(t_fun));
else
    fprintf(2, '%s is not currently supported\n');
    return;
end

% validate
vcFile_true = fullfile(vcDir_in, 'firings_true.mda');
if exist_file_(vcFile_true)
    fPlot_gt = read_cfg_('fPlot_gt');
    vcFile_raw = fullfile(vcDir_in, 'raw.mda');
    S_score = irc('validate-mda', vcFile_true, vcFile_firings, vcFile_raw, fPlot_gt); % assume that groundtruth file exists
    struct_save_(S_score, fullfile(vcDir_out, 'raw_geom_score.mat'), 1);
end
end %func


%--------------------------------------------------------------------------
% Call from irc.m
function compile_cuda_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function frewind_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function disperr_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function edit_prm_file_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function delete_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
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
function out1 = filesize_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = car_reject_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_copy_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cast_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2tr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsample_vr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_default_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_filter_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
% function out1 = gt2mda_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_dir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = meanSubt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = memory_matlab_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = recording_duration_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = squeeze_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = find_site_spk23_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mn2tn_wav_spk2_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = matchFileExt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
% function out1 = shift_trWav_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_copyas_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = set_bool_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ifeq_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2hash_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = dir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_sort_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refrac_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = map_index_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

% function [out1, out2] = readmda_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = mr2thresh_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = findNearSites_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = shift_range_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end

function [out1, out2, out3] = fopen_mda_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = fopen_nsx_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = plan_load_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = detect_spikes_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end