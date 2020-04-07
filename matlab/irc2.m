%--------------------------------------------------------------------------
% IronClust: terabyte-scale drift-resistent spike sorter 
%    based on KNN-graph density clustering and anatomical similarity
%
% Author: James Jun (2018-2020)
%   Center for Computational Mathematics
%   Flatiron Institute, a division of Simons Foundation


function varargout = irc2(vcDir_in, vcDir_out, vcFile_arg, vcArg3, vcArg4)
% irc2(vcDir_in, vcDir_out, vcFile_arg)
% irc2(vcCmd, vcArg1, vcArg2)

if nargin<1, vcDir_in = ''; end
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end
if nargin<4, vcArg3 = ''; end
if nargin<5, vcArg4 = ''; end

persistent vcFile_prm_

% fast response to call
if strcmpi(vcDir_in, 'call')
    [vcArg1, vcArg2] = deal(vcDir_out, vcFile_arg);
    switch nargout
        case 0, call_(vcArg1, vcArg2);
        case 1, varargout{1} = call_(vcArg1, vcArg2);
        case 2, [varargout{1}, varargout{2}] = call_(vcArg1, vcArg2);
        case 3, [varargout{1}, varargout{2}, varargout{3}] = call_(vcArg1, vcArg2);
        case 4, [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = call_(vcArg1, vcArg2);
        case 5, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}] = call_(vcArg1, vcArg2);
        case 6, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}] = call_(vcArg1, vcArg2);
        case 7, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}] = call_(vcArg1, vcArg2);
        case 8, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}, varargout{8}] = call_(vcArg1, vcArg2);
        otherwise, error('call_: too many output');
    end %switch
    return;
end

% batch processing. it uses default param for now
S_cfg = read_cfg_();
[fDetect, fSort, fAuto] = deal(get_(S_cfg, 'fForceRerun', 0));

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
% [fDetect, fSort, fAuto] = deal(exist_file_(vcDir_in) || exist_dir_(vcDir_in)); % cmd mode
[P, S0, fPlot_gt, fValidate] = deal([]); 
vcFile_prm = dir2prm_(vcArg1);
if isempty(vcFile_prm)
    vcFile_prm = vcFile_prm_;
elseif exist_file_(vcFile_prm)
    vcFile_prm_ = vcFile_prm;
end
switch lower(vcCmd)     
    % probe reader
    case {'resort', 'reclust', 'redo', 'rerun', 'rerun-validate'}
        S0 = rerun_(vcArg1, vcArg2, vcArg3, vcArg4, vcFile_prm_); 
        if strcmpi(vcCmd, 'rerun-validate')            
            vcFile_firings_mda = get_(S0.S_auto, 'vcFile_firings_mda');
            validate_(S0.P, vcFile_firings_mda);
        end
        return;
    case 'probe-image', probe_image_(vcArg1, vcArg2, vcArg3); return;
    
    % ui
    case {'remove-lock', 'unlock'}, remove_lock_(vcFile_prm); return;
    case 'export-prb', export_prb_json_(vcArg1, vcArg2); return;
    case 'probe', irc('probe', vcArg1); return;
    case 'traces', irc2_traces(vcArg1); return;
            
    % SpikeGLX functions
    case {'import-recordings', 'import-spikeglx', 'import-rhd', 'import-intan', 'import-neurosuite'}
        import_recordings_(vcArg1, vcArg2, vcArg3); return;
        
    % Git functions
    case 'push-readme', push_readme_(); return;    
    case {'push', 'git-push'}, git_push_(vcArg1); return;
    case {'pull', 'update', 'git-pull'}, git_pull_(vcArg1); return;
        
    % optimize
    case 'optimize-clear', optimize_clear_(vcArg1, vcArg2); return;
    case 'optimize-status', optimize_status_(vcArg1, vcArg2); return;
    case {'optimize', 'optimize-param', 'optimize-prmset', 'optimize-prm'}
        optimize_prmset_(vcArg1, vcArg2, vcArg3); return;
    case 'optimize-preview', optimize_prmset_(vcArg1, vcArg2, 1); return;
    case {'readme', 'edit-readme'}, edit_readme_(); return;
    % spikeforest2 interface
    case 'clear-jobs', clear_jobs_(vcArg1); return;
        
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
    case {'run-irc', 'run-ironclust'}, run_spikeforest2_('ironclust', vcArg1, vcArg2, vcArg3); return;
    case {'run-ks2', 'run-kilosort2'}, run_spikeforest2_('kilosort2', vcArg1, vcArg2, vcArg3); return;
        
    % direcly run kilosort2 from the source
    case {'ks2', 'kilosort2'}, run_ksort2(vcArg1, vcArg2, vcArg3); return;
        
    case 'describe-mda', describe_mda_(vcArg1, vcArg2); return;
    case 'compare-mda', compare_mda_(vcArg1, vcArg2); return;
    case 'export-mda', save_firings_mda_(vcFile_prm); return;
    case {'which', 'select', 'set'}
        if exist_file_(vcFile_prm)
            fprintf('Selected %s\n', vcFile_prm); 
        else
            fprintf(2, 'None selected\n');
        end
        return;
    case {'export-sf2', 'export-spikeforest2', 'export2spikeforest2'}, export_sf2_(); return;
    case {'export-ws', 'export-workspace', 'export'}, export_workspace_(vcArg1); return;
    case 'compile-deploy', compile_cuda_(vcArg1, '0'); return
    case 'compile-system', compile_cuda_(vcArg1, '1'); return
    case {'compile', 'install'}, compile_cuda_(vcArg1, vcArg2); return
    case {'readmda_header', 'mda-info'}, disp(readmda_header_(vcArg1)); return;
    case 'trim-mda', trim_mda_(vcArg1, vcArg2, vcArg3, vcArg4); return;
    case 'mcc', irc('mcc'); return; 
    case {'join-mda', 'join_mda', 'joinmda'}, join_mda_(vcArg1, vcArg2); return;
    case {'readmda', 'read-mda', 'read_mda'}
        mda = readmda_(vcArg1); assignWorkspace_(mda); return;
    case 'import-clip'
        [S0, P] = import_clip_(vcArg1); 
    case 'edit', edit_(vcFile_prm); return;
    case 'juxta', convert_mda_ui('english'); return;
    case 'emouse2mda', convert_mda('emouse', vcArg1, vcArg2); return;
    case {'extract-mda', 'crop-mda'}, convert_mda('extract-mda', vcArg1, vcArg2, vcArg3, vcArg4); return;
    case 'version'
        if nargout==0, version_(); 
        else, varargout{1} = version_(); 
        end
        return;
    case 'scoreboard', irc2_scoreboard(); return;
    case {'all', 'spikesort', 'detectsort', 'detect-sort', 'sort', 'auto', ...
            'describe', 'manual', 'validate', 'verify', 'spikesort-cache', ...
            'auto-verify', 'sort-verify', 'spikesort-verify', 'detectsort-verify', ...
            'auto-validate', 'sort-validate', 'spikesort-validate', 'detectsort-validate'}

        fprintf('irc2 (%s) opening %s\n', version_(), vcFile_prm);
        if isempty(vcFile_prm), fprintf(2, 'provide .prm file.\n'); return; end
        if ~exist_file_(vcFile_prm)
            fprintf(2, 'File does not exist: %s\n', vcFile_prm); return;
        end
        P = file2struct_(vcFile_prm);        
        switch lower(vcCmd)
            case {'spikesort-cache'}, clear_(); fDetect = 0; fSort = 0; fAuto=0;
            case {'detect-sort', 'spikesort', 'spikesort-verify', 'all'}, clear_(); fDetect = 1; fSort = 1; fAuto=1;
            case {'sort', 'sort-verify'}, fDetect = 0; fSort = 1; fAuto=1;   
            case {'auto', 'auto-verify'}, fDetect = 0; fSort = 0; fAuto=1;
            case 'describe', describe_(vcFile_prm); return;
            case {'verify', 'validate'}, validate_(P); return;
            case {'manual', 'ui'}, irc2_manual(P); return;
        end
        fValidate = contains(vcCmd, {'auto-verify', 'sort-verify', 'spikesort-verify', 'all'});
        vcDir_out = '';
    case 'benchmark'
        if nargout==0, benchmark_(vcArg1, vcArg2, vcArg3); 
        else, varargout{1} = benchmark_(vcArg1, vcArg2, vcArg3); 
        end
        return;
    case 'plot-drift', plot0_('drift', vcFile_prm, vcArg2); return;
    case 'plot-driftcorr', plot0_('driftcorr', vcFile_prm, vcArg2); return;
    case 'clear', clear_(vcArg1); vcFile_prm_=[]; return;
    case {'test-mcc', 'test_mcc', 'testmcc'}, test_mcc_(vcArg1); return;
    case {'test', 'test-clear'}
        vcDir_in = get_test_data_(vcArg1);
        if strcmpi(vcCmd, 'test-clear'), clear_(vcDir_in); end
        remove_lock_(vcDir_in);
        vcDir_out = '';
        fValidate = 1;
    case 'test-all', test_all_(); return;
%     case 'export', irc('export', vcArg1); return;
    case {'export-phy', 'phy'}, open_phy_(irc2phy(vcFile_prm, vcArg2)); return;
    case {'export-klusters', 'klusters', 'neurosuite'}
        open_klusters_(irc2klusters_v2(vcArg1, vcArg2));
        return;
    case {'export-jrclust', 'jrclust'}
        addpath(S_cfg.path_jrclust);
        cellfun(@(x,y)copyfile_(ircpath_(x), fullfile(S_cfg.path_jrclust, y)), ...
            S_cfg.patch_from_jrclust, S_cfg.patch_to_jrclust);
        geom2mat_(vcFile_prm);
        vcFile_prm_new = strrep(vcFile_prm, '.prm', '_jrclust.prm');
        copyfile(vcFile_prm, vcFile_prm_new, 'f');
        try
            jrc('import-irc', vcFile_prm_new);
            jrc('manual', vcFile_prm_new);
        catch
            disp(lasterr());
        end
%         open_jrc_(export_jrclust(vcFile_prm, vcArg2)); 
        return;
    otherwise % directory running mode
        if ~exist_file_(vcDir_in) && ~exist_dir_(vcDir_in)
            fprintf(2, 'invalid command: %s\n', vcCmd);
            return;
        end
        if matchFileExt_(vcDir_in, '.prm')
            P = file2struct_(vcDir_in);
            P.vcFile_prm = vcDir_in;
        end
        remove_lock_(vcDir_in);
        vcCmd=''; 
        fValidate = exist_file_(fullfile(vcDir_in, 'firings_true.mda')) && nargout==0;
end

fprintf('Running irc2.m (%s)\n', version_());
if isempty(P)
    fParfor = vcArg3;
    [P, vcDir_out] = makeParam_(vcDir_in, vcDir_out, vcFile_arg, fParfor);    
end
vcFile_prm_ = P.vcFile_prm;

% detect
[S0, fCached_detect] = detect_cache_(P, logical_(fDetect));

% sort
try
    [S0.S_clu, fCached_sort] = sort_cache_(S0, P, logical_(fSort));
catch ME % remove cache
    if fCached_detect
        fprintf(2, 'sort failed, retrying detect\n');
        pause(rand());
        S0 = detect_cache_(P, 1);
        [S0.S_clu, fCached_sort] = sort_cache_(S0, P, 1);
    else
        rethrow(ME);
    end
end

% auto
try
    [S0.S_auto, fCached_auto] = auto_cache_(S0, P, logical_(fAuto));
catch ME % remove cache
    if fCached_sort
        fprintf(2, 'auto failed, retrying sort\n');
        pause(rand());
        S0.S_clu = sort_cache_(S0, P, 1);
        S0.S_auto = auto_cache_(S0, P, 1);
    else
        rethrow(ME);
    end
end

vcFile_firings_mda = get_(S0.S_auto, 'vcFile_firings_mda');
if exist_file_(vcFile_firings_mda)
    if nargout==0
        copyfile(vcFile_firings_mda, fullfile(vcDir_out, 'firings.mda'));
    else
        varargout{1} = vcFile_firings_mda;
    end
end

% describe
S0.P=P;
describe_(S0);

% Validate
if fValidate, validate_(P, vcFile_firings_mda, S0.S_auto); end
end %func


%--------------------------------------------------------------------------
function [vcHash, csParam] = get_hash_(P, vcMode)
% usage
% -----
% vcHash = get_hash_(P_detect, vcMode) : vcMode={'detect','sort','auto'}

csParam_detect = {'version', 'vcFilter', 'freqLim', 'maxDist_site_spk_um', 'maxDist_site_um', ...
    'vcCommonRef', 'trimmean_pct', 'fWhiten', 'nSites_whiten', 'nChans_min_car', ...
    'qqFactor', 'nPc_spk', 'spkLim_ms', 'fInterp_fet', 'spkRefrac_ms', 'fft_thresh', ...
    'blank_thresh', 'blank_period_ms', 'fMatchedFilter_detect', 'prinvec_mode'};
csParam_sort = {'version', 'nPcPerChan', 'step_sec_drift', 'batch_sec_drift', ...
    'knn', 'nTime_max_drift', 'fMode_mlPc'};
csParam_auto = {'version', 'maxWavCor', 'merge_thresh_cc', 'spkJitter_ms_cc', ...
    't_burst_ms', 'min_snr_clu', 'spkRefrac_merge_ms', 'fUseSecondSite_merge', ...
    'merge_dist_thresh', 'merge_overlap_thresh'};

if ischar(P)
    P = file2struct_(P);
end
P.version = version_();
switch lower(vcMode)
    case 'detect', csParam = csParam_detect;
    case 'sort', csParam = csParam_sort;
    case 'auto', csParam = csParam_auto;
    otherwise, error('get_hash_: invalid mode'); 
end
vcHash = [vcMode, '_', struct2hash_(struct_copy_(P, csParam))];
end %func


%--------------------------------------------------------------------------
function [S0, fCached] = detect_cache_(P, fForce_detect)
if nargin<2, fForce_detect=0; end

csVar_bin = {'viTime_spk', 'viSite_spk', 'viSite2_spk', 'vrAmp_spk', ...
    'vrPow_spk', 'mrPos_spk', 'ccviSpk_site_load', 'ccviSpk_site2_load'};

vcDir_detect = fullfile(fileparts(P.vcFile_prm), get_hash_(P, 'detect'));
mkdir_(vcDir_detect);
vcFile_detect = fullfile(vcDir_detect, 'detect_irc.mat');
waitfor_lock_(vcDir_detect);
[S0, fCached] = deal([], 0);

if exist_file_(vcFile_detect) && ~fForce_detect
    try
        S0 = load(vcFile_detect);
        S0 = struct_load_bin_(S0.S_var, S0);
        fprintf('Loaded from cache: %s\n', vcFile_detect);
        fCached = 1;
    catch
    end
end
if isempty(S0)
    try
        lock_dir_(vcDir_detect);
        P.vcFile_prm = fullfile(vcDir_detect, 'detect.prm');
        S0 = detect_(P); 
        S0.P.vcFile_prm = P.vcFile_prm;
        S0.vcFile_detect = vcFile_detect;
        [S_, S_.S_var] = struct_save_bin_(S0, strrep(vcFile_detect,'_irc.mat','.irc'), csVar_bin);
        struct_save_(S_, vcFile_detect, 1);
        unlock_dir_(vcDir_detect);
    catch ME
        disp(ME.message());
        unlock_dir_(vcDir_detect);
        rethrow(ME);
    end
end
end %func


%--------------------------------------------------------------------------
function [S_clu, fCached] = sort_cache_(S0, P, fForce_sort)
if nargin<3, fForce_sort = 0; end

csVar_bin = {'rho', 'delta', 'nneigh', 'ordrho'};

vcDir_sort = fullfile(fileparts(S0.vcFile_detect), get_hash_(P, 'sort'));
mkdir_(vcDir_sort);
vcFile_sort = fullfile(vcDir_sort, 'sort_irc.mat');
waitfor_lock_(vcDir_sort);
[S_clu, fCached] = deal([], 0);

if exist_file_(vcFile_sort) && ~fForce_sort
    try
        S_clu = load(vcFile_sort);
        S_clu = struct_load_bin_(S_clu.S_var, S_clu);
        fprintf('Loaded from cache: %s\n', vcFile_sort);
        fCached = 1;
    catch
    end
end
if isempty(S_clu)
    try
        lock_dir_(vcDir_sort);
        P.vcFile_prm = strrep(vcFile_sort, '_irc.mat', '.prm');
        S_clu = sort_(S0, P);   
        S_clu.vcFile_sort = vcFile_sort;
        [S_, S_.S_var] = struct_save_bin_(S_clu, strrep(vcFile_sort,'_irc.mat','.irc'), csVar_bin);
        struct_save_(S_, vcFile_sort, 1);
        unlock_dir_(vcDir_sort);
    catch ME
        unlock_dir_(vcDir_sort);
%         disp(ME.message());
        rethrow(ME);
    end        
end
end %func


%--------------------------------------------------------------------------
% no need to use lock on auto since it's the final step
function [S_auto, fCached] = auto_cache_(S0, P, fForce_auto)

csVar_bin = {'viClu', 'cviSpk_clu'};

if nargin<3, fForce_auto = 0; end

vcDir_auto = fullfile(fileparts(S0.S_clu.vcFile_sort), get_hash_(P, 'auto'));
mkdir_(vcDir_auto);
vcFile_auto = fullfile(vcDir_auto, 'auto_irc.mat');
vcFile_firings_mda = fullfile(vcDir_auto, 'firings.mda');
% waitfor_lock_(vcDir_auto);

[S_auto, fCached] = deal([], 0);
if exist_file_(vcFile_auto) && ~fForce_auto
    try
        S_auto = load(vcFile_auto);
        S_auto = struct_load_bin_(S_auto.S_var, S_auto);
        fprintf('Loaded from cache: %s\n', vcFile_auto);
        if ~exist_file_(vcFile_firings_mda)
            S0.S_auto = S_auto;
            save_firings_mda_(S0, vcFile_firings_mda);
            fCached = 1;
        end
    catch
    end
end
if isempty(S_auto)
    try
%         lock_dir_(vcDir_auto);
        P.vcFile_prm = fullfile(vcDir_auto, 'auto.prm');
        S_auto = auto_(S0, P);

        % export to mda format
        S_auto.vcFile_firings_mda = vcFile_firings_mda;
        S0.S_auto = S_auto;
        save_firings_mda_(S0, vcFile_firings_mda);        

        [S_, S_.S_var] = struct_save_bin_(S_auto, strrep(vcFile_auto,'_irc.mat','.irc'), csVar_bin);
        struct_save_(S_, vcFile_auto, 1);  
%         unlock_dir_(vcDir_auto);
    catch ME
%         disp(ME.message());
%         unlock_dir_(vcDir_auto);
        rethrow(ME);
    end
end
end %func


%--------------------------------------------------------------------------
function waitfor_lock_(vcPath, lock_timeout, check_period)

if nargin<2, lock_timeout=3600; end
if nargin<3, check_period=1; end

vcFile_lock = lock_file_(vcPath);
t1=tic;
i=0;
while exist_file_(vcFile_lock)
    if i==0, fprintf('Waiting for lock: %s\n', vcFile_lock); end
    if toc(t1) > lock_timeout
        fprintf(2, 'Timeout (%0.1f s) for lock: %s\n', lock_timeout, vcFile_lock);
        break; 
    else        
        pause(check_period);
    end
    i=i+1;
end
end %func


%--------------------------------------------------------------------------
% create a name for a lock file
function vcFile_lock = lock_file_(vcPath)
[vcDir, vcFile, vcExt] = fileparts(vcPath);
vcFile_lock = fullfile(vcDir, ['.', vcFile, vcExt, '.lock']);
end %func


%--------------------------------------------------------------------------
function vcFile_lock = lock_dir_(vcPath)
vcFile_lock = lock_file_(vcPath);
if ~exist_file_(vcFile_lock)
    fid = fopen(vcFile_lock, 'w');
    fclose(fid);
end
end %func


%--------------------------------------------------------------------------
function vcFile_lock = unlock_dir_(vcPath)
vcFile_lock = lock_file_(vcPath);
if exist_file_(vcFile_lock)
    try
        delete(vcFile_lock, 'f');
    catch
    end
end
end %func


%--------------------------------------------------------------------------
function remove_lock_(csDir_rec, fParfor)
% recursively remove locks
% delete files associated with locks (they are incomplete)
if nargin<2, fParfor=0; end

if ischar(csDir_rec)
    if exist_dir_(csDir_rec)
        csDir_rec={csDir_rec}; 
    elseif exist_file_(csDir_rec)
        csDir_rec={fileparts(csDir_rec)};
    else
        fprintf(2, 'Does not exist: %s\n', csDir_rec);
        return;
    end
end
if fParfor
    parfor iDir=1:numel(csDir_rec)
        vS_dir = dir(fullfile(csDir_rec{iDir}, '**', '.*.lock'));
        csFiles_lock1 = arrayfun_(@(x)fullfile(fullfile(x.folder, x.name)), vS_dir);
        csDir_locked1 = cellfun_(@(x)strrep(x(2:end), '.lock', ''), csFiles_lock1);
        vnDelete(iDir) = numel(csDir_locked1);
        delete_(csFiles_lock1);
        rmdir_(csDir_locked1);
    end
else
    for iDir=1:numel(csDir_rec)
        vS_dir = dir(fullfile(csDir_rec{iDir}, '**', '.*.lock'));
        csFiles_lock1 = arrayfun_(@(x)fullfile(fullfile(x.folder, x.name)), vS_dir);
        csDir_locked1 = cellfun_(@(x)strrep(x(2:end), '.lock', ''), csFiles_lock1);
        vnDelete(iDir) = numel(csDir_locked1);
        delete_(csFiles_lock1);
        rmdir_(csDir_locked1);
    end
end
fprintf('Removed %d lock(s).\n', sum(vnDelete));
end %func


%--------------------------------------------------------------------------
function copyfile_(csPath_from, csPath_to)

if ischar(csPath_from), csPath_from={csPath_from}; end
if ischar(csPath_to), csPath_to={csPath_to}; end
assert(numel(csPath_from) == numel(csPath_to), 'number of elements must match');
for iFile = 1:numel(csPath_from)
    try
        copyfile(csPath_from{iFile}, csPath_to{iFile}, 'f');
    catch
        fprintf('Copy failed: %s to %s\n', csPath_from{iFile}, csPath_to{iFile});
    end
end
end %func


%--------------------------------------------------------------------------
function vcFile_mat = geom2mat_(vcFile_prm)
P = file2struct_(vcFile_prm);
vcFile_geom = P.probe_file;
geometry = csvread(vcFile_geom);
channels = 1:size(geometry,1);
probePad = [12, 12];
vcFile_mat = fullfile(fileparts(vcFile_prm), 'geom.mat');
struct_save_(makeStruct_(geometry, channels, probePad), vcFile_mat, vcFile_mat);
end %func


%--------------------------------------------------------------------------
function open_klusters_(csFile_par)
system_('%s %s', read_cfg_path_('path_klusters'), csFile_par{1});
end %func


%--------------------------------------------------------------------------
function open_phy_(vcFile_py)
system_('%s template-gui %s', read_cfg_path_('path_phy'), vcFile_py);
end %func


%--------------------------------------------------------------------------
function open_jrc_(obj)
addpath(read_cfg_path_('path_jrclust'));
obj.isCurate = 1;
obj.run();
end %func


%--------------------------------------------------------------------------
function [code, msg] = system_(varargin)
vcCmd = sprintf(varargin{:});
try
    disp(vcCmd);
    [code, msg] = system(vcCmd);
catch E
    code = -1;
    msg = E.message;    
end
if nargout==0
    if code==0
        fprintf('%s\n', msg);
    else
        fprintf(2, '%s\n', msg);
    end
end
end %func


%--------------------------------------------------------------------------
function vcPath = read_cfg_path_(vcKey)
if ispc()
    vcPath = read_cfg_([vcKey, '_pc']);
    if any(vcPath==' '), vcPath = ['"', vcPath, '"']; end
elseif ismac()
    vcPath = read_cfg_([vcKey, '_mac']);
elseif isunix()
    vcPath = read_cfg_([vcKey, '_lin']);
else
    fprintf(2, 'read_cfg_path_: unsupported OS: %s\n', vcKey);
    vcPath = []; % unsupported os
end
end


%--------------------------------------------------------------------------
% 11/6/18 JJJ: Displaying the version number of the program and what's used. #Tested
function [vcVer, vcDate, vcHash] = version_()

S_version = meta2struct_(ircpath_('version.txt'));
[vcVer, vcDate] = get_(S_version, 'version', 'date');
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
[~,~,vcExt] = fileparts(fname);
if ~strcmpi(vcExt, '.mda')
    error('File does not have .mda extension: %s', fname);
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
assert(numel(A) == prod(S), 'readmda_: dimension mismatch, incomplete data?');
A = reshape(A, S);
fclose(F);
end %func


%--------------------------------------------------------------------------
function test_all_()
irc2('test', 'monotrode');
irc2('test', 'tetrode');
irc2('test', 'static');
irc2('test', 'drift');
end %func


%--------------------------------------------------------------------------
% convert directory to prm file if directory path is given
function vcFile_prm = dir2prm_(vcDir_in)
if isempty(vcDir_in), vcDir_in=''; end
[vcDir1, vcFile1, vcExt1] = fileparts(vcDir_in);
switch lower(vcExt1)
    case '.prm', vcFile_prm = vcDir_in;        
    case '.mat', vcFile_prm = fullfile(vcDir1, strrep(vcFile1, '_irc', ''), '.prm');    
    case {'.bin', '.dat', '.mda'}
        vcFile_prm = dir_(fullfile(vcDir1, 'irc2', '*.prm'));
        if numel(vcFile_prm)==1
            vcFile_prm = vcFile_prm{1};
        else
            vcFile_prm = [];
        end
    case ''
        vcDir_out = fullfile(vcDir_in, 'irc2');
        S_prm = dir(fullfile(vcDir_out, '*.prm'));
        if numel(S_prm) == 1
            vcFile_prm = fullfile(S_prm.folder, S_prm.name);
        else
            vcFile_prm = fullfile(vcDir_out, 'raw_geom.prm');
            if ~exist_file_(vcFile_prm)
                vcFile_prm = '';
            end
        end
    otherwise
        vcFile_prm = vcDir_in;
end
end %func


%--------------------------------------------------------------------------
function save_clu_(S_clu, P)
% usage
% -----
% save_clu_(S_clu, P)

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

[P, S_cfg] = load_default_prm_();
S_mat = load(vcFile_mat);
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
P.nC_max = S_cfg.nC_max; % override nC_max (gpu parameter)
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
[~,~,vcExt] = fileparts(vcFile_prm);
switch lower(vcExt)
    case '.prm'
        S0 = load0_(vcFile_prm);
        S0.trPc_spk = load_fet_(S0, S0.P, 1);
        S0.trPc2_spk = load_fet_(S0, S0.P, 2);
    case '.mat'
        S0 = load(vcFile_prm);
        if isfield(S0, 'S_var')
            S0 = struct_load_bin_(S0.S_var, S0);
        end
end %switch
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
    elseif isnumeric(val1)
        varargout{iArg} = val1;
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
function [vcFile_score, S_score] = validate_(P, vcFile_firings_mda, S_auto)

t_fun=tic;
if nargin<2, vcFile_firings_mda=[]; end
if nargin<3, S_auto=[]; end

S_cfg = read_cfg_();
P.jitter = round(get_set_(S_cfg, 'spkJitter_ms_gt', 1) * P.sRateHz / 1000);
P.fCompute_snr_mda = get_set_(S_cfg, 'fCompute_snr_mda', 0);

vcFile_gt_mda = get_(P, 'vcFile_gt');
if ~exist_file_(vcFile_gt_mda), return; end
if isempty(vcFile_firings_mda)
    vcFile_firings_mda = fullfile(P.vcDir_out, 'firings.mda');
    vcFile_score = fullfile(P.vcDir_out, 'raw_geom_score.mat');
    [P.freqLim, P.freqLim_width] = deal(S_cfg.freqLim_gt, S_cfg.freqLim_width_gt);
    S_score = compare_mda_(vcFile_gt_mda, vcFile_firings_mda, P);
    struct_save_(S_score, vcFile_score, 1);
else
    vcFile_score = fullfile(fileparts(vcFile_firings_mda), 'score_irc.mat');
%     waitfor_lock_(vcFile_score);
%     if exist_file_(vcFile_score)
%         S_score = load(vcFile_score);
%     else
        S_score = compare_mda_(vcFile_gt_mda, vcFile_firings_mda, P);
        struct_save_(S_score, vcFile_score, 1);
%     end
end
if false
    S0 = load0_(P.vcFile_prm);
    S_score.vrSnr_clu = (S0.S_auto.vrScore_clu) * 20;
    S_score.vrSnr_gt = S_score.vrSnr_clu(S_score.viClu_gt);
end
hFig = S_score_plot_(S_score, S_cfg);
set(hFig, 'UserData', makeStruct_(P, S_cfg, vcFile_firings_mda, S_score));

fprintf('validate_: took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
function hFig = S_score_plot_(S_score, S_cfg)

if nargin<2, S_cfg=[]; end

if isempty(S_cfg)
    S_cfg = get_(S_score, 'S_cfg');
end
if isempty(S_cfg), S_cfg = read_cfg_(); end
fShow_text = get_set_(S_cfg, 'show_text_plot_gt', 1);
snr_max_plot = get_set_(S_cfg, 'snr_max_plot_gt', 20);
vrSnr_gt = get_(S_score, 'vrSnr_gt');
if ~isempty(vrSnr_gt)
    vcSnr_gt = sprintf(' >= SNR%0.1f', S_cfg.snr_thresh_gt);
    vlGt = vrSnr_gt>=S_cfg.snr_thresh_gt;   
    vcX_gt = 'vrSnr_gt';
else
    vrSnr_gt = rankorder_(S_score.vrAccuracy_gt, 'ascend');
    vcSnr_gt = '';
    vlGt = true(S_score.nGt, 1);
    vcX_gt = 'viGt_ordered';
    snr_max_plot = 0;
end

vrSnr_clu = get_(S_score, 'vrSnr_clu');
if ~isempty(vrSnr_clu)    
    vcSnr_clu = sprintf(' >= SNR%0.1f', S_cfg.snr_thresh_clu);
    vlClu = vrSnr_clu>=S_cfg.snr_thresh_clu;
    vcX_clu = 'vrSnr_clu';
else
    vrSnr_clu = rankorder_(S_score.vrAccuracy_clu, 'ascend');
    vcSnr_clu = '';
    vlClu = true(S_score.nClu,1);
    vcX_clu = 'viClu_ordered';
    snr_max_plot = 0;
end
snr_max = max([snr_max_plot, max(vrSnr_gt), max(vrSnr_clu)]);

disp_stats_(); % show caption        
disp_stats_(S_score.vrAccuracy_gt(vlGt), ['vrAccuracy_gt', vcSnr_gt]);
disp_stats_(S_score.vrF1_gt(vlGt), ['vrF1_gt', vcSnr_gt]);
disp_stats_(S_score.vrPrecision_gt(vlGt), ['vrPrecision_gt', vcSnr_gt]);
disp_stats_(S_score.vrRecall_gt(vlGt), ['vrRecall_gt', vcSnr_gt]);

fprintf('\n');      
disp_stats_(S_score.vrAccuracy_clu(vlClu), ['vrAccuracy_clu', vcSnr_clu]);
disp_stats_(S_score.vrF1_clu(vlClu), ['vrF1_clu', vcSnr_clu]);
disp_stats_(S_score.vrPrecision_clu(vlClu), ['vrPrecision_clu', vcSnr_clu]);
disp_stats_(S_score.vrRecall_clu(vlClu), ['vrRecall_clu', vcSnr_clu]);
fprintf('\n');

if fShow_text
    csText_gt = arrayfun_(@(x)sprintf('%d',x), 1:numel(vrSnr_gt))';  
    csText_clu = arrayfun_(@(x)sprintf('%d',x), 1:numel(vrSnr_clu))';
    text_gt_ = @(x)text_(vrSnr_gt, S_score.(x), csText_gt, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text_clu_ = @(x)text_(vrSnr_clu, S_score.(x), csText_clu, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
else
    text_gt_ = @(x)'';
    text_clu_ = @(x)'';
end
title_str_ = @(x)sprintf('n=%d, %0.1f/%0.1f, [%0.1f, %0.1f, *%0.1f, %0.1f, %0.1f]', ...
    numel(x), nanmean(x), nanstd(x), quantile(x, [.1,.25,.5,.75,.9]));
plot_gt_ = @(x){...
    plot(vrSnr_gt, S_score.(x), 'k.', vrSnr_gt(vlGt), S_score.(x)(vlGt), 'b.'), ...
    xylabel_([], vcX_gt, x, title_str_(S_score.(x)(vlGt)),1), text_gt_(x)};
plot_clu_ = @(x){...
    plot(vrSnr_clu, S_score.(x), 'k.', vrSnr_clu(vlClu), S_score.(x)(vlClu), 'b.'), ...
    xylabel_([], vcX_clu, x, title_str_(S_score.(x)(vlClu)),1), text_clu_(x)};

hFig = figure('Color','w', 'Name', ...
    sprintf('Groundtruth: %s, Sorted: %s', ...
        S_score.vcFile_gt_mda, S_score.vcFile_clu_mda));     
AX=[];
AX(end+1)=subplot(421); plot_gt_('vrAccuracy_gt'); 
AX(end+1)=subplot(423); plot_gt_('vrF1_gt');
AX(end+1)=subplot(425); plot_gt_('vrPrecision_F1_gt');
AX(end+1)=subplot(427); plot_gt_('vrRecall_F1_gt');
AX(end+1)=subplot(422); plot_clu_('vrAccuracy_clu');
AX(end+1)=subplot(424); plot_clu_('vrF1_clu');
AX(end+1)=subplot(426); plot_clu_('vrPrecision_F1_clu');
AX(end+1)=subplot(428); plot_clu_('vrRecall_F1_clu');
linkaxes(AX, 'xy');
axis(AX(1), [0 snr_max 0 100]);
% title_(AX(1), S_score.vcFile_gt_mda);
% title_(AX(5), S_score.vcFile_clu_mda);
end %func


%--------------------------------------------------------------------------
function [vi, viSort] = rankorder_(vr, vcOrder)
% warning: 32 bit addressing
if nargin<2, vcOrder = 'ascend'; end
n=numel(vr);
[~,viSort] = sort(vr, vcOrder);
if isGpu_(vr)
    vi = zeros(n,1,'int32', 'gpuArray');
    vi(viSort) = 1:n;
else
    vi=zeros(n,1,'int32');
    vi(viSort) = 1:n;
end
end %func


%--------------------------------------------------------------------------
function hAx = xylabel_(hAx, vcXLabel, vcYLabel, vcTitle, fGrid)

if isempty(hAx), hAx = gca; end
if nargin<4, vcTitle=''; end
if nargin<5, fGrid=[]; end

xlabel(hAx, vcXLabel, 'Interpreter', 'none');
ylabel(hAx, vcYLabel, 'Interpreter', 'none');
title_(hAx, vcTitle);
if fGrid, grid(hAx,'on'); end
end %func


%--------------------------------------------------------------------------
function hTitle = title_(hAx, vc)
% title_(vc)
% title_(hAx, vc)

if nargin==1, vc=hAx; hAx=[]; end
% Set figure title

if isempty(hAx), hAx = gca; end
hTitle = get_(hAx, 'Title');
if isempty(hTitle)
    hTitle = title(hAx, vc, 'Interpreter', 'none', 'FontWeight', 'normal');
else
    set_(hTitle, 'String', vc, 'Interpreter', 'none', 'FontWeight', 'normal');
end
end %func


%--------------------------------------------------------------------------
function vc = set_(vc, varargin)
% Set handle to certain values
% set_(S, name1, val1, name2, val2)

if isempty(vc), return; end
if isstruct(vc)
    for i=1:2:numel(varargin)        
        vc.(varargin{i}) = varargin{i+1};
    end
    return;
end
if iscell(vc)
    for i=1:numel(vc)
        try
            set(vc{i}, varargin{:});
        catch
        end
    end
elseif numel(vc)>1
    for i=1:numel(vc)
        try
            set(vc(i), varargin{:});
        catch
        end
    end
else
    try
        set(vc, varargin{:});
    catch
    end 
end
end %func


%--------------------------------------------------------------------------
function vc = disp_stats_(vr, vcCaption, nBlanks)

if nargin<3, nBlanks=[]; end
if isempty(nBlanks)
    nBlanks = 32; 
elseif nBlanks==0
    nBlanks = numel(vcCaption) + 2;
end
vcCaption_pad = repmat(' ', [1,nBlanks]);
if nargin<2, vcCaption = ''; end

if nargin==0
    vc = sprintf('%sn, mu/sd, (10,25,*50,75,90%%), [min-max]\t: ', vcCaption_pad);
else
    nCaption = min(nBlanks, numel(vcCaption));
    vcCaption_pad(1:nCaption) = vcCaption(1:nCaption);
    vr = vr(~isnan(vr));
    vr = vr(:);
    vc = sprintf('%s%d, %0.1f/%0.1f, (%0.1f, %0.1f, *%0.1f, %0.1f, %0.1f), [%0.1f-%0.1f]', ...
        vcCaption_pad, numel(vr), mean(vr), std(vr), quantile(vr, [.1,.25,.5,.75,.9]), min(vr), max(vr));
end
if nargout==0
    fprintf('%s\n', vc);
end
end %func


%--------------------------------------------------------------------------
function describe_mda_(vcFile_clu_mda, vcFile_raw_mda)
% usage
% -----
% describe_mda_(vcFile_clu_mda): describe in adc samples
% describe_mda_(vcFile_clu_mda, vcFile_raw_mda)

if nargin<2, vcFile_raw_mda=[]; end
if isempty(vcFile_raw_mda)
    % find out about the duration
    vcDir = fileparts(fileparts(vcFile_clu_mda));
    vcFile_raw_mda = fullfile(vcDir, 'raw.mda');
    if ~exist_file_(vcFile_raw_mda)
        vcFile_raw_mda = fullfile(fileparts(vcFile_clu_mda), 'raw.mda');
    end
    if ~exist_file_(vcFile_raw_mda)
        vcFile_raw_mda = ''; 
    end
end
if ~isempty(vcFile_raw_mda)
    vcFile_json = fullfile(fileparts(vcFile_raw_mda), 'params.json');
    S_json = loadjson_(vcFile_json);
    sRateHz = get_(S_json, 'samplerate');
else
    sRateHz = [];
end

if ~exist_file_(vcFile_clu_mda)
    error('%s does not exist', vcFile_clu_mda);
end

mr_clu = readmda_(vcFile_clu_mda)';
[viSite_spk, viTime_spk, viClu_spk] = ...
    deal(int32(mr_clu(:,1)), int64(mr_clu(:,2)), int32(mr_clu(:,3))); 
mr_clu = [];

cviSpk_clu = vi2cell_(viClu_spk);
vnSpk_clu = cellfun(@numel, cviSpk_clu);

% exclude zero-member clusters
viClu_keep = find(vnSpk_clu>0);
[nClu, cviSpk_clu, vnSpk_clu] = ...
    deal(numel(viClu_keep), cviSpk_clu(viClu_keep), vnSpk_clu(viClu_keep));


[cviSpk_site, nSites] = vi2cell_(viSite_spk);
vnSpk_site = cellfun(@numel, cviSpk_site);
vnDur_clu = cellfun(@(x)diff(viTime_spk([min(x), max(x)])), cviSpk_clu);
nDuration = max(viTime_spk) - min(viTime_spk);
nSpk = numel(viTime_spk);

csDesc = {};
csDesc{end+1} = sprintf('');    
csDesc{end+1} = sprintf('------------------------------');    
csDesc{end+1} = sprintf('%s', vcFile_clu_mda);
csDesc{end+1} = sprintf('------------------------------');    
csDesc{end+1} = sprintf('  # Spikes:               %d', nSpk);
csDesc{end+1} = sprintf('  # Clusters:             %d', nClu);
csDesc{end+1} = sprintf('  # Sites:                %d', nSites);
csDesc{end+1} = sprintf('  Duration (# samples):   %d', nDuration);   

if ~isempty(sRateHz)
    vrDur_clu = double(vnDur_clu) ./ sRateHz;
    vrRate_clu = vnSpk_clu ./ vrDur_clu;
    csDesc{end+1} = disp_stats_();
    csDesc{end+1} = disp_stats_(vrDur_clu, 'vrDur_clu');
    csDesc{end+1} = disp_stats_(vrRate_clu, 'vrRate_clu');
    csDesc{end+1} = disp_stats_(vnSpk_clu, 'vnSpk_clu');
    csDesc{end+1} = disp_stats_(vnSpk_site, 'vnSpk_site');
end

csDesc{end+1} = sprintf('------------------------------');
disp_cs_(csDesc);
end
    
    
%--------------------------------------------------------------------------
function S_score = compare_mda_(vcFile_gt_mda, vcFile_clu_mda, P)
% usage
% -----
% S_score = compare_mda_(vcFile_gt_mda, vcFile_clu_mda)
% S_score = compare_mda_(vcFile_gt_mda, vcFile_clu_mda, P)
nSamples_max = 2^10;

if nargin<3, P=[]; end
assert(exist_file_(vcFile_gt_mda) && exist_file_(vcFile_clu_mda), 'compare_mda_: Files must exist');

S_cfg = read_cfg_();

if isempty(P)
    P.vcFile = fullfile(fileparts(vcFile_gt_mda), 'raw.mda'); 
    S_json = loadjson_(fullfile(fileparts(vcFile_gt_mda), 'params.json'));
    P.sRateHz = get_(S_json, 'samplerate'); 
    [P.freqLim, P.freqLim_width, P.spkLim] = get_(S_cfg, 'freqLim_gt', 'freqLim_width_gt', 'spkLim_ms_gt');
    P.jitter = round(get_(S_cfg, 'spkJitter_ms_gt') * P.sRateHz / 1000); %1 ms jitter    
    P.fParfor=0;
end

fprintf('Validating cluster...\n'); t_fun = tic;
fCompute_snr_mda = get_set_(P, 'fCompute_snr_mda', get_(S_cfg, 'fCompute_snr_mda'));
fParfor = get_set_(P, 'fParfor', 1);

% read from firings.mda files
mr_gt = readmda_(vcFile_gt_mda); [viTimeGt, viGt] = deal(mr_gt(2,:)', mr_gt(3,:)'); mr_gt = [];
mr_clu = readmda_(vcFile_clu_mda); [viTimeClu, viClu] = deal(mr_clu(2,:)', mr_clu(3,:)'); mr_clu = [];

% compare cluster
[S_score, cviSpk_gt, cviSpk_clu] = compare_clu_(viGt, viTimeGt, viClu, viTimeClu, P.jitter, fParfor);

% compute the SNR
if fCompute_snr_mda
    t1=tic;
    vcFile_filt = strrep(P.vcFile, '.mda', '_filt.mda');
    mrWav_filt = [];
    if exist_file_(vcFile_filt)
        fprintf('\tLoading from %s...', vcFile_filt); 
        try mrWav_filt = readmda_(vcFile_filt); catch; end
        fprintf(' took %0.1fs\n', toc(t1));
    end
    if isempty(mrWav_filt)
        fprintf('\tFiltering...');
        mrWav_filt = fft_filter_transpose(readmda_(P.vcFile), P);
        writemda_(vcFile_filt , mrWav_filt);
        fprintf(' took %0.1fs, saved to %s\n', toc(t1), vcFile_filt);
    end
    tr2vr_peak_mean_ = @(tr)min(min(median(tr,2)));
    mr2tr_subsample_ = @(x,y)mr2tr_(mrWav_filt, P.spkLim, x(subsample_vr_(y, nSamples_max)));
    [vrPeak_gt, viSite_gt] = cellfun(@(x)tr2vr_peak_mean_(mr2tr_subsample_(viTimeGt, x)), cviSpk_gt);
    [vrPeak_clu, viSite_clu] = cellfun(@(x)tr2vr_peak_mean_(mr2tr_subsample_(viTimeClu, x)), cviSpk_clu);
    vrNoise_site = mr2rms_(mrWav_filt, 1e5)';
    vrSnr_gt = abs(vrPeak_gt(:)) ./ vrNoise_site(viSite_gt);
    vrSnr_clu = abs(vrPeak_clu(:)) ./ vrNoise_site(viSite_clu);
    mrWav_filt = [];
else
    [vrSnr_gt, vrSnr_clu, viSite_gt, viSite_clu, vrNoise_site] = deal([]);
end

S_score = struct_add_(S_score, vcFile_gt_mda, vcFile_clu_mda, P, S_cfg, ...
    vrSnr_gt, vrSnr_clu, viSite_gt, viSite_clu, vrNoise_site);

if nargout==0
    S_score_plot_(S_score, S_cfg);
    vcFile_score = fullfile(fileparts(vcFile_clu_mda), 'raw_geom_score.mat');
    struct_save_(S_score, vcFile_score, 1);
end
fprintf('\tValidation took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
function [S_compare, cviSpk_gt, cviSpk_clu] = compare_clu_(viGt, viTimeGt, viClu, viTimeClu, jitter, fParfor)
% usage
% -----
% [S_compare, cviSpk_gt, cviSpk_clu] = compare_clu_(viGt, viTimeGt, [], [], jitter, fParfor)
%    self-comparison
% [S_compare, cviSpk_gt, cviSpk_clu] = compare_clu_(viGt, viTimeGt, viClu, viTimeClu, jitter, fParfor)
%    cross-comparison

if nargin<4, fParfor = 0; end

[cviSpk_gt, nGt] = vi2cell_(int32(viGt));
cviTime_gt = cellfun(@(vi_)viTimeGt(vi_), cviSpk_gt, 'UniformOutput', 0);

if isempty(viClu) && isempty(viTimeClu)
    [cviSpk_clu, ~, cviTime_clu] = deal(cviSpk_gt, nGt, cviTime_gt);
else
    [cviSpk_clu, ~] = vi2cell_(int32(viClu));
    cviTime_clu = cellfun(@(vi_)viTimeClu(vi_), cviSpk_clu, 'UniformOutput', 0);
end
cviTime_clu = cviTime_clu(~cellfun(@isempty, cviTime_clu));
nClu = numel(cviTime_clu);

% Compute intersection
mnIntersect = zeros(nClu, nGt);
count_overlap_ = @(x,y)numel(find_overlap_(cviTime_gt{x}, cviTime_clu{y}, jitter));
for iGt=1:nGt
    for iClu = 1:nClu
        mnIntersect(iClu, iGt) = count_overlap_(iGt, iClu);
    end
end

% compute accuracy, precision, recall
vnSpk_clu = cellfun(@numel, cviTime_clu); vnSpk_clu = vnSpk_clu(:);
vnSpk_gt = cellfun(@numel, cviTime_gt); vnSpk_gt = vnSpk_gt(:)';
mrPrecision = mnIntersect ./ vnSpk_clu * 100;
mrRecall = mnIntersect ./ vnSpk_gt * 100;
mrSum_clu_gt = vnSpk_clu+vnSpk_gt;
mrAccuracy = mnIntersect ./ (mrSum_clu_gt-mnIntersect) * 100;
mrF1 = 2*mnIntersect ./ mrSum_clu_gt * 100;

[vrAccuracy_gt, viClu_gt, vrPrecision_gt, vrRecall_gt] = ...
    find_best_score_(mrAccuracy, mrPrecision, mrRecall);
[vrAccuracy_clu, viClu_clu, vrPrecision_clu, vrRecall_clu] = ...
    find_best_score_(mrAccuracy', mrPrecision', mrRecall');
[vrF1_gt, viClu_F1_gt, vrPrecision_F1_gt, vrRecall_F1_gt] = ...
    find_best_score_(mrF1, mrPrecision, mrRecall);
[vrF1_clu, viClu_F1_clu, vrPrecision_F1_clu, vrRecall_F1_clu] = ...
    find_best_score_(mrF1', mrPrecision', mrRecall');

S_compare = makeStruct_(vnSpk_clu, vnSpk_gt, jitter, nGt, nClu, ...
    mrPrecision, mrRecall, mrAccuracy, mrF1, ...
    vrAccuracy_gt, viClu_gt, vrPrecision_gt, vrRecall_gt, ...
    vrAccuracy_clu, viClu_clu, vrPrecision_clu, vrRecall_clu, ...
    vrF1_gt, viClu_F1_gt, vrPrecision_F1_gt, vrRecall_F1_gt, ...
    vrF1_clu, viClu_F1_clu, vrPrecision_F1_clu, vrRecall_F1_clu);
end %func


%--------------------------------------------------------------------------
function [viA_overlap, viB_overlap] = find_overlap_(A, B, tol)
% assume A,B are unique and sorted
% tol: maximum difference of valuse allowed to be considered as overlap

[na, nb] = deal(numel(A), numel(B));
L = false(na+nb,1); L(end-nb+1:end)=true; % membership: false is A, true is B
I = [1:na, 1:nb]';

% merge A and B and sort
[AB_srt, vi_srt] = sort([A(:); B(:)]);
L_srt = L(vi_srt);
I_srt = I(vi_srt);

vlOverlap = diff(AB_srt)<=tol;
viDiff = diff(L_srt);
viAB_srt = find(vlOverlap & viDiff>0); % find when no change is observed
viBA_srt = find(vlOverlap & viDiff<0); % find when no change is observed

viA_overlap = unique([I_srt(viAB_srt); I_srt(viBA_srt+1)]);
if nargout==1, return; end
viB_overlap = unique([I_srt(viAB_srt+1); I_srt(viBA_srt)]);
end %func


%--------------------------------------------------------------------------
function [S_compare, cviSpk_gt, cviSpk_clu] = compare_clu1_(viGt, viTimeGt, viClu, viTimeClu, jitter, fParfor)
% usage
% -----
% [S_compare, cviSpk_gt, cviSpk_clu] = compare_clu_(viGt, viTimeGt, [], [], jitter, fParfor)
%    self-comparison
% [S_compare, cviSpk_gt, cviSpk_clu] = compare_clu_(viGt, viTimeGt, viClu, viTimeClu, jitter, fParfor)
%    cross-comparison

if nargin<4, fParfor = 0; end

[cviSpk_gt, nGt] = vi2cell_(int32(viGt));
cviTime_gt = cellfun(@(vi_)unique(int32(viTimeGt(vi_)/jitter)), cviSpk_gt, 'UniformOutput', 0);

if isempty(viClu) && isempty(viTimeClu)
    [cviSpk_clu, ~, cviTime_clu] = deal(cviSpk_gt, nGt, cviTime_gt);
else
    [cviSpk_clu, ~] = vi2cell_(int32(viClu));
    cviTime_clu = cellfun(@(vi_)unique(int32(viTimeClu(vi_)/jitter)), cviSpk_clu, 'UniformOutput', 0);
end
cviTime_clu = cviTime_clu(~cellfun(@isempty, cviTime_clu));
nClu = numel(cviTime_clu);
viTime_min_clu = cellfun(@min, cviTime_clu);
viTime_max_clu = cellfun(@max, cviTime_clu);
miLim_clu = [viTime_min_clu(:), viTime_max_clu(:)];

% Compute intersection
mnIntersect = zeros(nClu, nGt);
% fParfor = get_set_(P, 'fParfor', 1);
if fParfor
    try
        parfor iGt=1:nGt
            mnIntersect(:, iGt) = compare_mda_gt_(cviTime_gt{iGt}, cviTime_clu, miLim_clu);
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iGt=1:nGt
        mnIntersect(:, iGt) = compare_mda_gt_(cviTime_gt{iGt}, cviTime_clu, miLim_clu);
    end
end

% compute accuracy, precision, recall
vnSpk_clu = cellfun(@numel, cviTime_clu); vnSpk_clu = vnSpk_clu(:);
vnSpk_gt = cellfun(@numel, cviTime_gt); vnSpk_gt = vnSpk_gt(:)';
mrPrecision = mnIntersect ./ vnSpk_clu * 100;
mrRecall = mnIntersect ./ vnSpk_gt * 100;
mrSum_clu_gt = vnSpk_clu+vnSpk_gt;
mrAccuracy = mnIntersect ./ (mrSum_clu_gt-mnIntersect) * 100;
mrF1 = 2*mnIntersect ./ mrSum_clu_gt * 100;

[vrAccuracy_gt, viClu_gt, vrPrecision_gt, vrRecall_gt] = ...
    find_best_score_(mrAccuracy, mrPrecision, mrRecall);
[vrAccuracy_clu, viClu_clu, vrPrecision_clu, vrRecall_clu] = ...
    find_best_score_(mrAccuracy', mrPrecision', mrRecall');
[vrF1_gt, viClu_F1_gt, vrPrecision_F1_gt, vrRecall_F1_gt] = ...
    find_best_score_(mrF1, mrPrecision, mrRecall);
[vrF1_clu, viClu_F1_clu, vrPrecision_F1_clu, vrRecall_F1_clu] = ...
    find_best_score_(mrF1', mrPrecision', mrRecall');

S_compare = makeStruct_(vnSpk_clu, vnSpk_gt, jitter, nGt, nClu, ...
    mrPrecision, mrRecall, mrAccuracy, mrF1, ...
    vrAccuracy_gt, viClu_gt, vrPrecision_gt, vrRecall_gt, ...
    vrAccuracy_clu, viClu_clu, vrPrecision_clu, vrRecall_clu, ...
    vrF1_gt, viClu_F1_gt, vrPrecision_F1_gt, vrRecall_F1_gt, ...
    vrF1_clu, viClu_F1_clu, vrPrecision_F1_clu, vrRecall_F1_clu);
end %func


%--------------------------------------------------------------------------
function [vrVrms_site, vrVsd_site] = mr2rms_(mr, max_sample)
% uses median to estimate RMS
if nargin<2, max_sample = []; end
if ~isempty(max_sample), mr = subsample_mr_(mr, max_sample, 1); end
vrVrms_site = median(abs(mr));
vrVrms_site = single(vrVrms_site) / 0.6745;
if nargout>=2, vrVsd_site = std(single(mr)); end
end % func


%--------------------------------------------------------------------------
function vi = subsample_vr_(vi, nMax)
if numel(vi)>nMax
    nSkip = floor(numel(vi)/nMax);
    if nSkip>1, vi = vi(1:nSkip:end); end
    if numel(vi)>nMax
        try
            nRemove = numel(vi) - nMax;
            viRemove = round(linspace(1, numel(vi), nRemove));
            viRemove = min(max(viRemove, 1), numel(vi));
            vi(viRemove) = [];
        catch
            vi = vi(1:nMax);
        end
    end
end
end %func


%--------------------------------------------------------------------------
function [mr, vi] = subsample_mr_(mr, nMax, dimm)
%[mr, vi] = subsample_mr_(mr, nMax, dimm)
% subsample the column
if nargin<3, dimm = 2; end
if isempty(nMax), return; end

n = size(mr,dimm);
nSkip = max(floor(n / nMax), 1);
vi = 1:nSkip:n;
if nSkip==1, return; end
vi = vi(1:nMax);

switch dimm
    case 2
        mr = mr(:,vi);
    case 1
        mr = mr(vi,:);
end

if nargout>=2
    if n > nMax
        vi = 1:nSkip:n;
        vi = vi(1:nMax);
    else
        vi = 1:n;
    end
end
end %func


%--------------------------------------------------------------------------
function [tr, miRange] = mr2tr_(mr, spkLim, viTime, viSite, fMeanSubt)
% tr: nSamples x nSpikes x nChans

if nargin<4, viSite=[]; end %faster indexing
if nargin<5, fMeanSubt=0; end

% JJJ 2015 Dec 24
% vr2mr2: quick version and doesn't kill index out of range
% assumes vi is within range and tolerates spkLim part of being outside
% works for any datatype
if isempty(viTime), tr=[]; return; end
[N, M] = size(mr);
if ~isempty(viSite), M = numel(viSite); end
if iscolumn(viTime), viTime = viTime'; end
if isGpu_(mr)
    [viTime, viSite] = deal(gpuArray_(viTime), gpuArray_(viSite)); 
else
    [viTime, viSite] = deal(gather_(viTime), gather_(viSite)); 
end
viTime0 = [spkLim(1):spkLim(end)]'; %column
miRange = bsxfun(@plus, int32(viTime0), int32(viTime));
miRange = min(max(miRange, 1), N);
miRange = miRange(:);
if isempty(viSite)
    tr = mr(miRange,:);
else
    tr = mr(miRange, viSite);
end
tr = reshape(tr, [numel(viTime0), numel(viTime), M]);

if fMeanSubt
%     trWav1 = single(permute(trWav1, [1,3,2])); 
    tr = single(tr);
    dimm1 = size(tr);
    tr = reshape(tr, size(tr,1), []);
    tr = bsxfun(@minus, tr, mean(tr)); %mean subtract
    tr = reshape(tr, dimm1);    
end
end %func


%--------------------------------------------------------------------------
function [vrScore, viRow, vr1, vr2] = find_best_score_(mrScore, mr1, mr2)
[vrScore, viRow] = max(mrScore, [], 1);
vi_mr = sub2ind(size(mrScore), viRow, 1:size(mrScore,2));
[vr1, vr2] = deal(mr1(vi_mr), mr2(vi_mr));
end %func


%--------------------------------------------------------------------------
function vnIntersect = compare_mda_gt_(viTimeA, cviTimeB, miLimB)
% usage
% ----
% vnIntersect = compare_mda_gt_(viTime1, cviTime_clu)
% mnIntersect = compare_mda_gt_(cviTime_gt, cviTime_clu)
% viTimeA and cviTimeB must be sorted (satisfied when unique() used)

nB = numel(cviTimeB);

if iscell(viTimeA)
    cviTimeA = viTimeA;
    vnIntersect = zeros(nB, numel(cviTimeA));
    for iGt=1:numel(cviTimeA)
        vnIntersect(:, iGt) = compare_mda_gt_(cviTimeA{iGt}, cviTimeB, miLimB);
    end
    return;
end

vnIntersect = zeros(nB,1);
if isempty(viTimeA), return; end
% viTimeA = sort(viTimeA);
[minA, maxA] = deal(viTimeA(1), viTimeA(end)); 
overlap_ = @(x)sum(ismembc(viTimeA,x) | ismembc(viTimeA+1,x) | ismembc(viTimeA-1,x));
[minB, maxB] = deal(miLimB(:,1), miLimB(:,2));
viB_update = find(...
    (minA >= minB & minA <= maxB) | (maxA >= minB & maxA <= maxB) | ...
    (minB >= minA & minB <= maxA) | (maxB >= minA & maxB <= maxA));
vnIntersect(viB_update) = arrayfun(@(x)overlap_(cviTimeB{x}), viB_update);
end %func


% %--------------------------------------------------------------------------
% function save0_(S0)
% if nargin<1, S0 = []; end
% if isempty(S0), S0 = get(0, 'UserData'); end
% 
% vcFile_prm_ = strrep(S0.P.vcFile_prm, '.prm', '');
% 
% % save separately
% % viTime_spk: nSpk x 1: int64
% % viSite_spk, viSite2_spk: nSpk x 1: int32
% % vrAmp_spk: nSpk x 1: single
% % vrPow_spk: nSpk x 1: single
% % mrPos_spk: nSpk x 2: single
% 
% csVar_spk = {'viTime_spk', 'viSite_spk', 'viSite2_spk', 'vrAmp_spk', ...
%     'vrPow_spk', 'mrPos_spk', 'ccviSpk_site_load', 'ccviSpk_site2_load'};
% [S0, S0.S_var] = struct_save_bin_(S0, [vcFile_prm_, '_spk.irc'], csVar_spk);
% 
% trPc_spk = gather_(get_(S0, 'trPc_spk'));
% S0.trPc_spk = [];
% fSave_fet = get_set_(S0.P, 'fSave_fet', 0);
% if ~isempty(trPc_spk) && fSave_fet
%     S0.trPc_spk = [];
%     [S0.dimm_fet, S0.type_fet] = write_bin_([vcFile_prm_, '_fet.irc'], trPc_spk);
% end
% 
% trPc2_spk = gather_(get_(S0, 'trPc2_spk'));
% S0.trPc2_spk = [];
% if ~isempty(trPc2_spk) && fSave_fet
%     S0.trPc2_spk = [];
%     write_bin_([vcFile_prm_, '_fet2.irc'], trPc2_spk);
% end
% 
% S_clu = get_(S0, 'S_clu');  S0.S_clu = [];
% if ~isempty(S_clu), save_clu_(S_clu, S0.P); end
% 
% S_auto = get_(S0, 'S_auto');  S0.S_auto = [];
% if ~isempty(S_auto), save_auto_(S_auto, S0.P); end
% 
% struct_save_(S0, [vcFile_prm_, '_irc.mat'], 1);
% end


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
% 2020/1/17: `csVar_load` added to select variables to load
function S = struct_load_bin_(S_var, S, csVar_load)
% Usage
% -----
% S = struct_load_bin_(S_var, S)
% S = struct_load_bin_(S_var, S, csVar_load)
% S = struct_load_bin_(S_var, [], csVar_load)

% save variables and clear field
if nargin<2, S = struct(); end
if nargin<3, csVar_load={}; end

import_struct_(S_var); % import all fields in this struct
fid_r = fopen(vcFile, 'r');
for iVar = 1:numel(csName_var)    
    [type_, dimm_, name_] = deal(csType_var{iVar}, cDimm_var{iVar}, csName_var{iVar});    
    if ~isempty(csVar_load)
        fSkip = ~contains(name_, csVar_load);
    else
        fSkip = 0;
    end
    if iscell(type_)
        S.(csName_var{iVar}) = load_cell_(fid_r, type_, dimm_, fSkip);
    else
        if fSkip
            fseek(fid_r, prod(dimm_) * bytesPerSample_(type_), 'cof');
        else
            S.(csName_var{iVar}) = load_bin_(fid_r, type_, dimm_);            
        end
    end
end
fclose(fid_r);
end %func


%--------------------------------------------------------------------------
function cVal = load_cell_(fid_r, cType, cDimm, fSkip)
if nargin<4, fSkip = 0; end

if ~iscell(cType)
    cVal = load_bin_(fid_r, cType, cDimm);
else
    cVal = cell(cType);
    for iCell1 = 1:numel(cVal)  
        [type_, dimm_] = deal(cType{iCell1}, cDimm{iCell1});
        if iscell(type_)
            cVal{iCell1} = load_cell_(fid_r, type_, dimm_);
        else
            if fSkip
                fseek(fid_r, prod(dimm_) * bytesPerSample_(type_), 'cof');
            else
                cVal{iCell1} = load_bin_(fid_r, type_, dimm_);
            end
        end
    end 
end
end %func


%--------------------------------------------------------------------------
function S0 = load0_(vcFile_prm)
t_fun=tic;
P = file2struct_(vcFile_prm);
S0 = detect_cache_(P);
S0.S_clu = sort_cache_(S0, P);
S0.S_auto = auto_cache_(S0, P);
fprintf('took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
function S0 = rerun_(vcArg1, vcArg2, vcArg3, vcArg4, vcFile_prm)
% usage
% rerun_(vcFile_prm, param1, param2, param3
% rerun_(vcArg1, vcArg2, vcArg3, vcArg4, vcFile_prm)

if nargin<5, vcFile_prm=''; end
S0=[];
    
% determine vcFile_prm
vcFile_prm1 = dir2prm_(vcArg1);
if exist_file_(vcArg1)
    vcFile_prm = vcFile_prm1;
    csParam = {vcArg2, vcArg3, vcArg4};
else
    fprintf('Using %s\n', vcFile_prm); 
    csParam = {vcArg1, vcArg2, vcArg3, vcArg4};
end
if ~exist_file_(vcFile_prm)
    fprintf(2, 'Specify a parameter file\n');
    return;
end

% update param
remove_lock_(fileparts(vcFile_prm));
P = file2struct_(vcFile_prm);
P = struct_merge_(P, cs2struct_(csParam));
P.vcFile_prm = vcFile_prm;
edit_prm_file_(P, P.vcFile_prm);
S0 = detect_cache_(P);
S0.S_clu = sort_cache_(S0, P);
S0.S_auto = auto_cache_(S0, P);
describe_(vcFile_prm);
end %func


%--------------------------------------------------------------------------
% negative index means from the end, 0 index means end index
function vc1 = strsplit_get_(vc,delim,idx)
cs = strsplit(vc, delim);
idx = mod(idx-1,numel(cs))+1;
vc1 = cs{idx};
end %func


%--------------------------------------------------------------------------
% clear ironclust output from specified directories
function clear_(vcFile_prm)
% usage
% vcDir_in = clear_(vcDir)
% vcDir_in = clear_(vcFile_prm)
% vcDir_in = clear_(vcFile_list_txt)

% clear specified recording
dir_prm_ = @(x)dir(fullfile(x, '*.prm'));
try
    if exist_dir_(vcFile_prm) 
        % directory is passed
        vcDir = vcFile_prm;
        vcDir_irc2 = fullfile(vcDir, 'irc2');
        if exist_dir_(vcDir_irc2)
            % output dir passed
            rmdir_(vcDir_irc2);
            fprintf('Removed %s\n', vcDir_irc2);
        end
    elseif exist_file_(vcFile_prm)
        [vcDir, vcFile, vcExt] = fileparts(vcFile_prm);
        switch lower(vcExt)
            case '.txt'
                csDir_rec = load_batch_(vcFile_prm);
                rmdir_(cellfun_(@(x)fullfile(x, 'irc2'), csDir_rec));
            case '.prm'
                % preserves .prm file
                vcFile_prm_ = fullfile(vcDir, vcFile);
                delete([vcFile_prm_, '*.irc']);
                delete([vcFile_prm_, '*_irc.mat']);
                delete_(fullfile(vcDir, '*_score.mat'));
                rmdir_(fullfile(vcDir, 'detect_*'));
                fprintf('Cleared %s\n', vcFile_prm);
        end
    else
        fprintf(2, 'Does not exist: %s\n', vcFile_prm);
    end
catch
    fprintf(2, 'Nothing is cleared.\n');
end
end %func


%--------------------------------------------------------------------------
function rmdir_(csDir)
if ischar(csDir)
    if exist_dir_(csDir)
        try
            rmdir(csDir, 's');
        catch        
        end
    end
else
    try
        parfor iDir = 1:numel(csDir)
            if exist_dir_(csDir{iDir})
                try
                    rmdir(csDir{iDir}, 's');
                catch        
                end
            end
        end
    catch
        for iDir = 1:numel(csDir)
            if exist_dir_(csDir{iDir})
                try
                    rmdir(csDir{iDir}, 's');
                catch        
                end
            end
        end
    end
end
end %func

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
    if exist_file_(vcFile_prm)
        S0 = load0_(vcFile_prm);
    else
        fprintf(2, 'File does not exist: %s\n\tSort the file first.\n', vcFile_mat);
        return;
    end
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
nPcPerChan = S0.P.nPc_spk;

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
    csDesc{end+1} = sprintf('    Matched Filter:         %d', get_set_(P, 'fMatchedFilter_detect', 0));    
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
%     csDesc{end+1} = sprintf('    #Features/event:        %d', nFeatures);    
    csDesc{end+1} = sprintf('    #PC/chan:               %d', nPcPerChan);
catch
end
try
   
    S_auto = get_(S0, 'S_auto');
    csDesc{end+1} = sprintf('Cluster');       
    csDesc{end+1} = sprintf('    #Clusters:              %d', get_(S_auto, 'nClu'));
    csDesc{end+1} = sprintf('    #Unique events:         %d', get_(S_auto, 'nSpk_unique'));
    csDesc{end+1} = sprintf('    min. spk/clu:           %d', S_auto.P.min_count);
    csDesc{end+1} = sprintf('    Cluster method:         %s', S_auto.P.vcCluster);
    csDesc{end+1} = sprintf('    knn:                    %d', S_auto.P.knn);
    csDesc{end+1} = sprintf('    step_sec_drift:         %0.1fs', S_auto.P.step_sec_drift);
    csDesc{end+1} = sprintf('    batch_sec_drift:        %0.1fs', S_auto.P.batch_sec_drift);
    csDesc{end+1} = sprintf('Auto-merge');   
    csDesc{end+1} = sprintf('    merge_overlap_thresh:   %0.3f', get_(S_auto.P, 'merge_overlap_thresh'));    
    csDesc{end+1} = sprintf('    delta_cut:              %0.3f', get_set_(S_auto.P, 'delta_cut', 1));
    csDesc{end+1} = sprintf('    merge_thresh_cc:        %0.3f', get_(S_auto.P, 'merge_thresh_cc'));
    csDesc{end+1} = sprintf('    maxWavCor:              %0.3f', get_(S_auto.P, 'maxWavCor'));
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

    csDesc{end+1} = sprintf('memory usage (GiB):         %0.3f', max([memory_detect, memory_sort, memory_auto])/2^30);
    csDesc{end+1} = sprintf('    detect:                 %0.3f', memory_detect/2^30);
    csDesc{end+1} = sprintf('    sort:                   %0.3f', memory_sort/2^30);
    csDesc{end+1} = sprintf('    auto-merge:             %0.3f', memory_auto/2^30);

    csDesc{end+1} = sprintf('Execution');
    csDesc{end+1} = sprintf('    irc2 version:           %s', get_(P, 'vcVersion'));
    csDesc{end+1} = sprintf('    fGpu (GPU use):         %d', P.fGpu);
    csDesc{end+1} = sprintf('    fParfor (parfor use):   %d', P.fParfor);
    csDesc{end+1} = sprintf('    fLargeRecording:        %d', isLargeRecording_(P)); 
    csDesc{end+1} = sprintf('    Parameter file:         %s', P.vcFile_prm);
    csDesc{end+1} = sprintf('------------------------------');        
catch
end
if nargout==0
    disp_cs_(csDesc); %cellfun(@(x)disp(x), csDesc);
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
t_fun = tic;
try
    if ischar(S0)
        vcFile_prm = S0;
        S_ = load(strrep(vcFile_prm, '.prm', '_irc.mat'), 'S_var'); 
        S0 = struct_load_bin_(S_.S_var, [], {'viSite_spk', 'viTime_spk'});
        S_ = load(strrep(vcFile_prm, '.prm', '_auto_irc.mat'), 'S_var', 'viClu'); 
        S0.S_auto = struct_load_bin_(S_.S_var, [], {}); 
    else
        vcFile_prm = S0.P.vcFile_prm;
    end
    if isempty(vcFile_firings_mda)
        vcFile_firings_mda = fullfile(fileparts(vcFile_prm), 'firings.mda');
    end
    mr = [double(S0.viSite_spk(:)), double(S0.viTime_spk(:)), double(S0.S_auto.viClu(:))];
    S0 = []; %free memory
catch
    error('save_firings_mda_: invalid format');
end
writemda_(vcFile_firings_mda, mr');
fprintf('Wrote to %s, took %0.1fs\n', vcFile_firings_mda, toc(t_fun));
end %func


%--------------------------------------------------------------------------
% 64-bit addressing compatible
function writemda_(vcFile, var)
writemda_fid(vcFile, var);
end %func
            

%--------------------------------------------------------------------------
% merge mutual neighbors
function S_auto = knn_overlap_merge_(S_auto, miKnn_spk, thresh)

nClu_pre = S_auto.nClu;
[vrScore_clu, viClu_nn_clu] = isolation_score_(S_auto, miKnn_spk);
cviClu_clu = cell(nClu_pre,1);
for iClu=1:nClu_pre
    if vrScore_clu(iClu) < thresh
        cviClu_clu{iClu} = [iClu; viClu_nn_clu(iClu)];
    else
        cviClu_clu{iClu} = iClu;
    end
end
viMapClu_new = cell2map_(cviClu_clu);
vlUpdate = S_auto.viClu > 0;
S_auto.viClu(vlUpdate) = viMapClu_new(S_auto.viClu(vlUpdate));
nClu = max(viMapClu_new);
fprintf('\tknn_overlap_merge_: %d->%d units\n', nClu_pre, nClu);
end %func


%--------------------------------------------------------------------------
% auto merge
function S_auto = auto_(S0, P)

fprintf('\nauto-merging...\n'); runtime_automerge = tic;

% Merge based on KNN-graph
S_auto = postCluster_(S0.S_clu, P, S0.viSite_spk); % peak merging

% Merge based on knn overlap
merge_overlap_thresh = get_set_(P, 'merge_overlap_thresh', 1);
miKnn_spk=[];
try
    if merge_overlap_thresh>0 && merge_overlap_thresh<1
        miKnn_spk = load_miKnn_spk_(P, S0.viSite_spk);
        S_auto = knn_overlap_merge_(S_auto, miKnn_spk, merge_overlap_thresh);
        S_auto = S_auto_refrac_(S_auto, P, S0.viTime_spk); % refractory violation removal
        S_auto = S_auto_refresh_(S_auto, 1, S0.viSite_spk);
    end
catch    
end

S_auto = wave_ccm_merge_(S0, S_auto, P);
S_auto = wave_similarity_merge_(S0, S_auto, P);

S_auto = S_auto_refrac_(S_auto, P, S0.viTime_spk); % refractory violation removal
S_auto = S_auto_refresh_(S_auto, 1, S0.viSite_spk);
S_auto = S_auto_sort_(S_auto, 'viSite_clu', S0.viSite_spk);
S_auto.memory_auto = memory_matlab_();
S_auto.runtime_automerge = toc(runtime_automerge);
fprintf('\tauto-merging took %0.1fs (fGpu=%d, fParfor=%d)\n', ...
    S_auto.runtime_automerge, P.fGpu, P.fParfor);

% compute quality metrics
if ~isempty(miKnn_spk)
    S_auto.vrScore_clu = isolation_score_(S_auto, miKnn_spk);
end
end %func


%--------------------------------------------------------------------------
function S_auto = wave_ccm_merge_(S0, S_auto, P)
% calculate cross-correlogram

P.merge_thresh_cc = get_set_(P, 'merge_thresh_cc', .9);
if P.merge_thresh_cc > 0 && P.merge_thresh_cc <1 && S_auto.nClu > 1
    S_auto0 = S_auto;
    nClu_pre = S_auto.nClu;
    fprintf('\tMerging based on cross-correlogram...\n\t'); t_cc=tic;
    switch 2
        case 1, cviClu_clu = calc_ccm_(S_auto.viClu, S0.viTime_spk, P);
        case 2, cviClu_clu = calc_ccm_space_(S_auto.viClu, S0.viTime_spk, S0.viSite_spk, P);
        case 3, cviClu_clu = calc_ccm_drift_(S0.S_clu.S_drift, ...
                    S_auto.viClu, S0.viTime_spk, S0.viSite_spk, P);
    end
    % remap clusters
    fParfor = get_set_(P, 'fParfor', 1);
    viMapClu_new = cell2map_(cviClu_clu, [], fParfor);
    vlUpdate = S_auto.viClu > 0;
    S_auto.viClu(vlUpdate) = viMapClu_new(S_auto.viClu(vlUpdate));
    S_auto.nClu = sum(unique(viMapClu_new)>0);
    if S_auto.nClu>0
        fprintf('\tMerged clusters cc>=%0.2f (%d->%d), took %0.1fs\n', ...
            P.merge_thresh_cc, nClu_pre, S_auto.nClu, toc(t_cc));    
    else
        S_auto = S_auto0;
        fprintf('\tNo waveforms were merged\n');
    end
end
end %func


%--------------------------------------------------------------------------
% calculate cross-correlogram
function cviClu_clu = calc_ccm_space_(viClu_spk, viTime_spk, viSite_spk, P)

merge_thresh_cc = get_set_(P, 'merge_thresh_cc', .5);
jitter = round(get_set_(P,'spkJitter_ms_cc',1) * P.sRateHz / 1000); %1 ms jitter

[cviSpk_clu, nClu] = vi2cell_(int32(viClu_spk));
cviTime_clu = cellfun(@(vi_)viTime_spk(vi_), cviSpk_clu, 'UniformOutput', 0);
cviClu_clu = cell(nClu, 1);
vnSpk_clu = cellfun(@numel, cviSpk_clu); vnSpk_clu = vnSpk_clu(:);
viSite_clu = cellfun(@(x)mode(viSite_spk(x)), cviSpk_clu); viSite_clu = viSite_clu(:);
miSites = P.miSites;
S_fun = makeStruct_(cviTime_clu, viSite_clu, miSites, vnSpk_clu, jitter, merge_thresh_cc, nClu);
% Compute intersection
if false
    [mrPrecision, mrRecall] = precision_recall_(S_fun, 1);
end
fParfor = get_set_(P, 'fParfor', 1) && nClu>1;
if fParfor
    try
        parfor iClu = 1:nClu
            cviClu_clu{iClu} = calc_ccm_space_for_(S_fun, iClu);
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iClu = 1:nClu
        cviClu_clu{iClu} = calc_ccm_space_for_(S_fun, iClu);
    end
end
end


%--------------------------------------------------------------------------
function [mrPrecision, mrRecall] = precision_recall_(S_fun, fSpatialMask)
if nargin<2, fSpatialMask=1; end
% treat iClu as groundtruth
t_fun=tic;
csVar_imported = import_struct_(S_fun);
[mrPrecision, mrRecall] = deal(nan(nClu,'single'));
[cviTime_clu, miSites, jitter] = get_(S_fun, 'cviTime_clu', 'miSites', 'jitter');
jitter = jitter * 1; 
parfor iClu = 1:nClu
    count_overlap_ = @(j)numel(find_overlap_(cviTime_clu{iClu}, cviTime_clu{j}, jitter));
    if fSpatialMask
        vnIntersect = zeros(nClu,1,'single');
        viClu2 = find(ismember(viSite_clu, miSites(:, viSite_clu(iClu))));
        vnIntersect(viClu2) = arrayfun(@(j)count_overlap_(j), viClu2);
    else
        vnIntersect = arrayfun(@(j)count_overlap_(j), 1:nClu)';
    end
    mrPrecision(:,iClu) = vnIntersect ./ vnSpk_clu;
    mrRecall(:,iClu) = vnIntersect ./ vnSpk_clu(iClu);
end  
mrPrecision1 = set_diag_(mrPrecision, zeros(nClu,1));
mrRecall1 = set_diag_(mrRecall, nan(nClu,1));
fprintf('precision_recall_: took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
function viClu = calc_ccm_space_for_(S_fun, iClu)
% treat iClu as groundtruth

csVar_imported = import_struct_(S_fun);
count_overlap_ = @(j)numel(find_overlap_(cviTime_clu{iClu}, cviTime_clu{j}, jitter));
viClu2 = find(ismember(viSite_clu, miSites(:, viSite_clu(iClu))));
vnIntersect2 = arrayfun(@(j)count_overlap_(j), viClu2);
% vrPrecision2 = vnIntersect2 ./ vnSpk_clu(viClu2);
vrRecall2 = vnIntersect2 ./ vnSpk_clu(iClu);
vrCC2 = vrRecall2;
% vrCC2 = max(vrPrecision2, vrRecall2);
viClu = viClu2(vrCC2 >= merge_thresh_cc);
end %func


%--------------------------------------------------------------------------
% calculate cross-correlogram
function cviClu_clu = calc_ccm_drift_(S_drift, viClu_spk, viTime_spk, viSite_spk, P)

miSites = P.miSites;
merge_thresh_cc = get_set_(P, 'merge_thresh_cc', .9);
jitter = round(get_set_(P,'spkJitter_ms_cc',1) * P.sRateHz / 1000); %1 ms jitter

t_fun=tic;
viLim_drift = S_drift.viLim_drift;
nDrift = S_drift.nTime_drift;
nClu = max(viClu_spk);
cviClu_clu = cell(nClu, 1);
for iDrift = 1:nDrift
    viDrift1 = find(S_drift.mlDrift(:,iDrift));
    % find spikes in this drift zone
    viSpk1 = cell2mat_(arrayfun_(@(x)(viLim_drift(x):(viLim_drift(x+1)-1))', viDrift1));
    [viTime_spk1, viClu_spk1, viSite_spk1] = multifun_(@(x)x(viSpk1), viTime_spk, viClu_spk, viSite_spk);
    cviSpk_clu1 = vi2cell_(viClu_spk1, nClu);
    vnSpk_clu1 = cellfun(@numel, cviSpk_clu1);
    viClu1 = int32(find(vnSpk_clu1>=P.min_count));
    [vnSpk_clu1, cviSpk_clu1] = deal(vnSpk_clu1(viClu1), cviSpk_clu1(viClu1));
    cviTime_clu1 = cellfun_(@(vi_)viTime_spk1(vi_), cviSpk_clu1, 'UniformOutput', 0);
    count_overlap_ = @(i,j)numel(find_overlap_(cviTime_clu1{i}, cviTime_clu1{j}, jitter));
    nClu1 = numel(viClu1);
    cviClu_clu1 = cell(nClu1,1);
    viSite_clu1 = cellfun(@(x)mode(viSite_spk1(x)), cviSpk_clu1);
    parfor iiClu1 = 1:nClu1
        viiClu2 = find(ismember(viSite_clu1, miSites(:, viSite_clu1(iiClu1))));
        vnIntersect2 = arrayfun(@(j)count_overlap_(iiClu1,j), viiClu2);
        vrCC2 = max(vnIntersect2 ./ vnSpk_clu1(viiClu2), vnIntersect2 / vnSpk_clu1(iiClu1));
        cviClu_clu1{iiClu1} = viClu1(viiClu2(vrCC2 >= merge_thresh_cc));
    end
    cviClu_clu(viClu1) = cellfun_(@(x,y)[x;y], cviClu_clu(viClu1), cviClu_clu1);
    fprintf('.');
end
cviClu_clu = cellfun_(@unique, cviClu_clu);
fprintf('\n\tcalc_ccm_drift_: took %0.1fs\n', toc(t_fun));
end


%--------------------------------------------------------------------------
% calculate cross-correlogram
function cviClu_clu = calc_ccm_(viClu_spk, viTime_spk, P)

merge_thresh_cc = get_set_(P, 'merge_thresh_cc', .5);
jitter = round(get_set_(P,'spkJitter_ms_gt',1) * P.sRateHz / 1000); %1 ms jitter

[cviSpk_clu, nClu] = vi2cell_(int32(viClu_spk));
cviTime_clu = cellfun(@(vi_)viTime_spk(vi_), cviSpk_clu, 'UniformOutput', 0);
cviClu_clu = cell(nClu, 1);
vnSpk_clu = cellfun(@numel, cviSpk_clu); vnSpk_clu = vnSpk_clu(:);

% Compute intersection
fParfor = get_set_(P, 'fParfor', 1);
count_overlap_ = @(i,j)numel(find_overlap_(cviTime_clu{i}, cviTime_clu{j}, jitter));
if fParfor
    try
        parfor iClu = 1:nClu
            vnIntersect1 = arrayfun(@(j)count_overlap_(iClu,j), 1:nClu)';
            vrCC1 = max(vnIntersect1 ./ vnSpk_clu, vnIntersect1 / vnSpk_clu(iClu));
            cviClu_clu{iClu} = find(vrCC1 >= merge_thresh_cc);
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iClu = 1:nClu
        vnIntersect1 = arrayfun(@(j)count_overlap_(iClu,j), 1:nClu)';
        vrCC1 = max(vnIntersect1 ./ vnSpk_clu, vnIntersect1 / vnSpk_clu(iClu));
        cviClu_clu{iClu} = find(vrCC1 >= merge_thresh_cc);
    end
end
end


%--------------------------------------------------------------------------
% calculate cross-correlogram
function cviClu_clu = calc_ccm1_(viClu_spk, viTime_spk, P)

merge_thresh_cc = get_set_(P, 'merge_thresh_cc', .5);
jitter = round(get_set_(P,'spkJitter_ms_cc',1) * P.sRateHz / 1000); %1 ms jitter

[cviSpk_clu, nClu] = vi2cell_(int32(viClu_spk));
cviTime_clu = cellfun(@(vi_)unique(int32(viTime_spk(vi_)/jitter)), cviSpk_clu, 'UniformOutput', 0);
viTime_min_clu = cellfun(@min, cviTime_clu);
viTime_max_clu = cellfun(@max, cviTime_clu);
miLim_clu = [viTime_min_clu(:), viTime_max_clu(:)];
cviClu_clu = cell(nClu, 1);
vnSpk_clu = cellfun(@numel, cviSpk_clu); vnSpk_clu = vnSpk_clu(:);

% Compute intersection
fParfor = get_set_(P, 'fParfor', 1);
if fParfor
    try
        parfor iClu = 1:nClu
            vnIntersect1 = compare_mda_gt_(cviTime_clu{iClu}, cviTime_clu, miLim_clu);
            vrCC1 = max(vnIntersect1 ./ vnSpk_clu, vnIntersect1 / vnSpk_clu(iClu));
            cviClu_clu{iClu} = find(vrCC1 >= merge_thresh_cc);
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iClu = 1:nClu
        vnIntersect1 = compare_mda_gt_(cviTime_clu{iClu}, cviTime_clu, miLim_clu);
        vrCC1 = max(vnIntersect1 ./ vnSpk_clu, vnIntersect1 / vnSpk_clu(iClu));
        cviClu_clu{iClu} = find(vrCC1 >= merge_thresh_cc);
    end
end
end


%--------------------------------------------------------------------------
function [S_auto, nRemoved] = S_auto_refrac_(S_auto, P, viTime_spk)

t_fun = tic;

% remove refractory spikes
nSkip = get_set_(P, 'nSkip_refrac', 4);
nRefrac = round(get_set_(P, 'spkRefrac_merge_ms', .25) * P.sRateHz / 1000);
nClu = max(S_auto.viClu);
nTotal = sum(S_auto.viClu>0);
fParfor = get_set_(P, 'fParfor', 1);
cviSpk_remove_clu = cell(nClu,1);
cviSpk_clu = vi2cell_(S_auto.viClu, S_auto.nClu);
cviTime_clu = cellfun_(@(x)viTime_spk(x), cviSpk_clu);
if fParfor
    try
        parfor iClu = 1:nClu
            cviSpk_remove_clu{iClu} = S_auto_refrac_clu_(...
                cviTime_clu{iClu}, cviSpk_clu{iClu}, nRefrac, nSkip);
        end
    catch
    end
end
if ~fParfor
    for iClu = 1:nClu
        cviSpk_remove_clu{iClu} = S_auto_refrac_clu_(...
            cviTime_clu{iClu}, cviSpk_clu{iClu}, nRefrac, nSkip);
    end    
end
viSpk_remove = cat(1, cviSpk_remove_clu{:});
S_auto.viClu(viSpk_remove) = 0;
nRemoved = numel(viSpk_remove);
% assert(nRemoved == numel(unique(viSpk_remove)))

fprintf('Removed %d/%d (%0.1f%%) duplicate spikes, took %0.1fs\n', ...
    nRemoved, nTotal, nRemoved/nTotal*100, toc(t_fun));
end %func


%--------------------------------------------------------------------------
function viSpk_remove = S_auto_refrac_clu_(viTime1, viSpk, nRefrac, nSkip)
% removal loop
vlKeep1 = true(size(viTime1));
while true
    viKeep1 = find(vlKeep1);
    viRefrac_ = find(diff(viTime1(viKeep1)) < nRefrac) + 1;
    if isempty(viRefrac_), break; end
    vlKeep1(viKeep1(viRefrac_(1:nSkip:end))) = false;
end
viSpk_remove = viSpk(~vlKeep1);
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
function S_auto = wave_similarity_merge_(S0, S_auto, P)


P.maxWavCor = get_set_(P, 'maxWavCor', .98);
if P.maxWavCor>=1 || P.maxWavCor<=0 || S_auto.nClu == 1
    return; 
end

S_auto0 = S_auto;
fprintf('\tMerging templates...\n\t'); t_fun=tic;
nClu_pre = S_auto.nClu;
S_clu = get_(S0, 'S_clu');
vrRho = S_clu.rho;
[viLim_drift, mlDrift] = get_(S_clu.S_drift, 'viLim_drift', 'mlDrift');

% fprintf('\tAutomated merging based on waveform similarity...\n'); t_template=tic;
viClu = S_auto.viClu;
[ccviSpk_site_load, ccviSpk_site2_load, type_fet, dimm_fet, mrPv, vrThresh_site, viTime_spk] = ...
    get_(S0, 'ccviSpk_site_load', 'ccviSpk_site2_load', 'type_fet', ...
        'dimm_fet', 'mrPv_global', 'vrThresh_site', 'viTime_spk');
nShift_max = round(P.spkRefrac_merge_ms/2 * P.sRateHz / 1000);
viShift = -nShift_max:nShift_max;
[knn, nSpk_min, vcFile_prm, maxWavCor] = get_(P, 'knn', 'knn', 'vcFile_prm', 'maxWavCor');
nClu = S_auto.nClu;

nSites = size(P.miSites,2);
S_param = makeStruct_(nClu, nSites, knn, vcFile_prm, nSpk_min, ...
    vrRho, viClu, viTime_spk, viLim_drift, ccviSpk_site_load, ccviSpk_site2_load, ...
    type_fet, dimm_fet, mrPv, vrThresh_site, viShift, mlDrift, maxWavCor, P);

% drift compensation
S_param.S_pos_clu = [];
if get_(P, 'merge_dist_thresh') > 0
    try
        S_param.S_pos_clu = static_position_clu_(setfield(S0, 'S_auto', S_auto));
        S_param.merge_dist_thresh = get_(P, 'merge_dist_thresh');
    catch
    end
end

% compute pairwise distance in parallel by sites
[cviClu_clu_site, cvlExist_site] = deal(cell(nSites, 1));
fParfor = get_set_(P, 'fParfor', 1) && nSites > 1;
if fParfor && nSites>1
    try
        parfor iSite = 1:nSites
            [cviClu_clu_site{iSite}, cvlExist_site{iSite}] = wave_similarity_site_(iSite, S_param);
        end
    catch
    end
end

% merge cluster pairwise distance
cviClu_clu = cell(nClu, 1);
vlExist_clu = false(1, nClu);
for iSite = 1:nSites
    if isempty(cviClu_clu_site{iSite})
        [cviClu_clu_site{iSite}, cvlExist_site{iSite}] = wave_similarity_site_(iSite, S_param);
    end    
    cviClu_clu = cellfun_(@(x,y)[x(:);y(:)], cviClu_clu, cviClu_clu_site{iSite});
    vlExist_clu = vlExist_clu | cvlExist_site{iSite};
end
cviClu_clu = cellfun_(@unique, cviClu_clu);
viClu_remove = find(~vlExist_clu);
fprintf('\n')

% merge and remap clusters
viMapClu_new = cell2map_(cviClu_clu, viClu_remove, fParfor);
vlUpdate = S_auto.viClu > 0;
S_auto.viClu(vlUpdate) = viMapClu_new(S_auto.viClu(vlUpdate));
S_auto.nClu = sum(unique(viMapClu_new)>0);

if S_auto.nClu>0
    fprintf('\tMerged waveforms (%d->%d->%d), took %0.1fs\n', ...
        nClu_pre, S_auto.nClu+numel(viClu_remove), S_auto.nClu, toc(t_fun));    
else
    S_auto = S_auto0;
    fprintf('\tNo waveforms were merged\n');
end
S_auto.P=P;
end %func


%--------------------------------------------------------------------------
function viSpk1 = get_viSpk_site_(S_fet, iFet, iSite)
switch iFet
    case 1, cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site_load);
    case 2, cviSpk_load = cellfun_(@(x)x{iSite}, S_fet.ccviSpk_site2_load);
end
viSpk1 = cat(1, cviSpk_load{:});
end %func


%--------------------------------------------------------------------------
function vi = unique_(vi)
if isempty(vi), return; end
vi = sort(vi);
vi(diff(vi)==0)=[];
end %func


%--------------------------------------------------------------------------
function cviSpk = separate_burst_(viSpk, viTime, nSamples_burst)
viiPre = find(diff(sort(viTime(viSpk))) >= nSamples_burst);
cviSpk = {viSpk(viiPre), viSpk(viiPre+1)};
end %func


%--------------------------------------------------------------------------
function [cviClu_clu, vlExist_clu] = wave_similarity_site_(iSite1, S_auto)
% free of KNN

% Load KNN and identify neighbors per cluster
csVar_imported = import_struct_(S_auto);

fUseSecondSite = get_set_(P, 'fUseSecondSite_merge', 1);
nDrift = size(mlDrift, 1);
viSpk1 = get_viSpk_site_(S_auto, 1, iSite1);
cvii1_drift = vi2cell_(discretize(viSpk1, viLim_drift), nDrift);
[vrRho1, viClu1, viTime1] = deal(S_auto.vrRho(viSpk1), viClu(viSpk1), viTime_spk(viSpk1));

thresh1 = -abs(vrThresh_site(iSite1)); % negative detection
min_snr_clu = get_(P, 'min_snr_clu');
if min_snr_clu>0
    thresh1 = thresh1 * min_snr_clu / P.qqFactor;
end

iT_peak = 1 - P.spkLim(1);
mrPv_peak = mrPv(iT_peak,:);
nSamples_burst = round(get_set_(P, 't_burst_ms', 20) * P.sRateHz / 1000);
if fUseSecondSite
    viSpk2 = get_viSpk_site_(S_auto, 2, iSite1);    
else
    viSpk2 = [];
end
fSecondSite = ~isempty(viSpk2);
if fSecondSite   
    cvii2_drift = vi2cell_(discretize(viSpk2, viLim_drift), nDrift);
    [vrRho2, viClu2, viTime2] = deal(S_auto.vrRho(viSpk2), viClu(viSpk2), viTime_spk(viSpk2));
    cvrRho_fet = {vrRho1, vrRho2};
    cviClu_fet = {viClu1, viClu2};
    cviTime_fet = {viTime1, viTime2};
    ccvii_drift_fet = {cvii1_drift, cvii2_drift};
else
    cvrRho_fet = {vrRho1};
    cviClu_fet = {viClu1};
    cviTime_fet = {viTime1};
    ccvii_drift_fet = {cvii1_drift};
end

[cviClu1_drift, cviClu2_drift, ccviSpk1_drift, ccviSpk2_drift] = deal(cell(nDrift,1));
for iFet=1:(fSecondSite+1)
    [vrRho1, viClu1, cvii1_drift] = deal(...
        cvrRho_fet{iFet}, cviClu_fet{iFet}, ccvii_drift_fet{iFet});    
    for iDrift = 1:nDrift      
        vii1 = cat(1, cvii1_drift{mlDrift(:,iDrift)});
        if numel(vii1) < nSpk_min, continue; end
        [vrRho11, viClu11] = deal(vrRho1(vii1), viClu1(vii1));
        [cviiSpk_clu_, ~, viClu_uniq] = vi2cell_(viClu11, nClu);
        nClu1 = numel(viClu_uniq);
        [viClu_drift1, cviSpk_drift1] = deal(nan(1,nClu1), cell(1,nClu1));
        for iiClu1 = 1:nClu1
            iClu1 = viClu_uniq(iiClu1);
            vii_ = cviiSpk_clu_{iClu1};   
            vrRho11_ = vrRho11(vii_);
            vii1_ = vii1(vii_(vrRho11_>=median(vrRho11_)));
            if numel(vii1_) >= nSpk_min
                cviSpk_drift1{iiClu1} = vii1_;
                viClu_drift1(iiClu1) = iClu1;
            end    
        end
        vi1_keep = find(~isnan(viClu_drift1));   
        switch iFet
            case 1
                cviClu1_drift{iDrift} = viClu_drift1(vi1_keep);
                ccviSpk1_drift{iDrift} = cviSpk_drift1(vi1_keep);    
            case 2
                cviClu2_drift{iDrift} = viClu_drift1(vi1_keep);
                ccviSpk2_drift{iDrift} = cviSpk_drift1(vi1_keep);    
        end
    end
end

% load trPc1 and trPc2
[cviClu_drift, ctrPc_drift] = deal(cell(nDrift,1));
for iFet = 1:numel(cviTime_fet)
    trPc1 = load_fet_site_(S_auto, iFet, iSite1);
    viTime1 = cviTime_fet{iFet};
    for iDrift = 1:nDrift
        [viClu_drift2, cmrPc_drift2] = deal([], {});
        if iFet==1
            viClu_drift1 = cviClu1_drift{iDrift};
            cviSpk_drift1 = ccviSpk1_drift{iDrift};
        else
            viClu_drift1 = cviClu2_drift{iDrift};
            cviSpk_drift1 = ccviSpk2_drift{iDrift};
        end 
        for ic = 1:numel(cviSpk_drift1)
            [iClu1, vii1] = deal(viClu_drift1(ic), cviSpk_drift1{ic}); 
            cvii1 = separate_burst_(vii1, viTime1, nSamples_burst);
            for iBurst=1:numel(cvii1)
                if numel(cvii1{iBurst}) < nSpk_min, continue; end
                tr_ = trPc1(:,:,cvii1{iBurst});
                switch 1
                    case 1, mrPc1 = mean(tr_,3);
                    case 2, mrPc1 = median(tr_,3);
                end
                if mrPv_peak*mrPc1(:,1) <= thresh1 % && peak_amp1 <= min(min(mrPv * mrPc1))/2
                    viClu_drift2(end+1) = iClu1;
                    cmrPc_drift2{end+1} = mrPc1;
                end
            end
        end
        if ~isempty(viClu_drift2)
            if iFet==1
                cviClu_drift{iDrift} = viClu_drift2;
                ctrPc_drift{iDrift} = cat(3, cmrPc_drift2{:});
            else
                cviClu_drift{iDrift} = [cviClu_drift{iDrift}, viClu_drift2];
                ctrPc_drift{iDrift} = cat(3, ctrPc_drift{iDrift}, cmrPc_drift2{:});
            end
        end
    end
    trPc1 = [];
end

% distance calculation
vlExist_clu = false(1, nClu);
vlExist_clu([cviClu_drift{:}]) = true;
cviClu_clu = num2cell((1:nClu)');
norm_mr_ = @(mr)mr ./ sqrt(sum(mr.^2,1)); 
tr2mr_pv_norm_ = @(tr,mr)norm_mr_(reshape(mr*reshape(tr,size(tr,1),[]),[],size(tr,3))); 

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
    if ~isempty(S_pos_clu)
%         mrPos_clu_drift1 = S_pos_clu.mrPos_clu;
        mrPos_clu_drift1 = [nanmedian(S_pos_clu.mrX_clu_drift(:,viDrift1),2), ...
            nanmedian(S_pos_clu.mrY0_clu_drift(:,viDrift1),2)]; 
    end
    for iiClu1 = 1:numel(viClu1)
        iClu1 = viClu1(iiClu1);
        mrWav11 = pc2wav_shift_(trPc_clu1(:,:,iiClu1), mrPv, viShift);
        viClu2_ = viClu2(max(mrWav11' * mrWav_clu2, [], 1) >= maxWavCor);
        if isempty(viClu2_), continue; end
        viClu2_ = setdiff(viClu2_, iClu1);
        if isempty(viClu2_), continue; end
        if ~isempty(S_pos_clu)
            vrDist2_ = pdist2_(mrPos_clu_drift1(iClu1,:), mrPos_clu_drift1(viClu2_,:));
            viClu2_(vrDist2_ > merge_dist_thresh) = [];
        end        
        if isempty(viClu2_), continue; end
        cviClu_clu{iClu1} = [cviClu_clu{iClu1}; viClu2_(:)];
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
% waveform correlation measures
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
end %func


%--------------------------------------------------------------------------
function [vr, vi] = cvr2vr_vi_(cvr)
vr = cell2mat_(cvr);
vn1 = cellfun(@(x)size(x,1), cvr);
vi = cell2mat_(arrayfun(@(x)repmat(x, vn1(x),1), 1:numel(cvr), 'UniformOutput', 0)');
end %func


%--------------------------------------------------------------------------
function vr = cell2mat_(cvr, empty_val)
if nargin<2, empty_val = []; end
% create a matrix that is #vectors x # cells
% remove empty
vi = find(cellfun(@(x)~isempty(x), cvr));
vr = cell2mat(cvr(vi));
if ~isempty(empty_val) && numel(vi) < numel(cvr)
    vr1 = repmat(empty_val, size(cvr));
    vr1(vi) = vr;
    vr = vr1;
end
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
if false % disable parfor for large recordings
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
S_clu.memory_sort = memory_sort;
S_clu.runtime_sort = toc(runtime_sort);
end %func


%--------------------------------------------------------------------------
function [mlPc, nFeatures] = get_mlPc_(S0, P)
nPcPerChan = get_set_(P, 'nPcPerChan', 0);
nPc_max = P.nC_max;
[nPc_spk, nSites_spk] = deal(S0.dimm_fet(1), S0.dimm_fet(2));

if isempty(nPcPerChan), nPcPerChan = 0; end
if numel(nPcPerChan) == 1
    if nPcPerChan == 0   
        switch get_set_(P, 'fMode_mlPc', 2)
            case 1, vnPc_site = nPc_spk:-1:1;
            case 2, vnPc_site = nPc_spk * ones(1, nSites_spk);
            case 3, vnPc_site = nPc_spk * ones(1, P.nSites_fet);
        end
        for iSite = numel(vnPc_site):-1:1
            if sum(vnPc_site) <= nPc_max, break; end
            vnPc_site(iSite:end) = vnPc_site(iSite:end) - 1;
            vnPc_site(vnPc_site<0) = 0; % non-negative
        end
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

nFeatures = sum(mlPc(:));
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

[~, vcFile_, ~] = fileparts(P.vcFile_prm);
if strcmp(vcFile_, 'auto')
    vcFile_prm = fullfile(dir_up_(P.vcFile_prm, 1), 'sort.prm');
else
    vcFile_prm = P.vcFile_prm;
end

vcFile_knn_site1 = strrep(vcFile_prm, '.prm', sprintf('_knn_%d.irc', iSite1));
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


[viSite_spk, viSite2_spk] = get_(S0, 'viSite_spk', 'viSite2_spk');
t_fun = tic;
[mlPc, nFeatures] = get_mlPc_(S0, P);    
[viLim_drift, mlDrift] = get_(S_drift, 'viLim_drift', 'mlDrift');
nSpk = numel(viSite_spk);
nSites = size(P.miSites,2);
vcFile_prm = P.vcFile_prm;

S_page = makeStruct_(P, mlDrift, mlPc, viLim_drift, nSites, vcFile_prm);
S_page = struct_merge_(S_page, struct_copy_(S0, ...
    'type_fet', 'dimm_fet', 'ccviSpk_site_load', 'ccviSpk_site2_load'));

% return schedules
[miSpk_lim_out, miSpk_lim_in, miDrift_lim_out, miDrift_lim_in] = plan_sort_page_(S_drift, P);
nPages = size(miSpk_lim_out,1);
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
viNneigh = zeros(nSpk, 1, 'int64');
S_global = makeStruct_(S_drift, miSpk_lim_out, miSpk_lim_in, miDrift_lim_out, viSite_spk, viSite2_spk);

vcFile_miKnn = [strrep(vcFile_prm, '.prm', ''), '_knn_*.irc'];
delete_(vcFile_miKnn);
fprintf('sort_page_: calculating Rho...\n'); t_rho = tic;
for iPage = 1:nPages  
    fprintf('Page %d/%d ', iPage, nPages); t_ = tic;
    [S_page1, viSpk_in1, viSpk_out1] = prepare_page_(S_page, S_global, iPage);
    vrRho(viSpk_in1) = rho_page_(S_page1);
    fprintf(' took %0.1fs\n', toc(t_));
end %for
fprintf('calculating Rho took %0.1fs\n', toc(t_rho));

fprintf('sort_page_: calculating Delta...\n'); t_delta = tic;
for iPage = 1:nPages    
    fprintf('Page %d/%d ', iPage, nPages); t_ = tic;    
    [S_page1, viSpk_in1, viSpk_out1] = prepare_page_(S_page, S_global, iPage);
    [vrDelta(viSpk_in1), viNneigh(viSpk_in1)] = delta_page_(S_page1, vrRho(viSpk_out1));
    fprintf(' took %0.1fs\n', toc(t_));
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
fParfor = get_set_(P, 'fParfor',1) && nSites>1;

if fParfor
    try
        parfor iSite = 1:nSites
            S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
                cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);
            cvrRho_in1{iSite} = rho_paged_site_(S_page1, S_site1, iSite);
        end %for
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iSite = 1:nSites
        S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
            cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);
        cvrRho_in1{iSite} = rho_paged_site_(S_page1, S_site1, iSite);
        fprintf('.');
    end
end

vrRho_in1 = zeros(nSpk1, 1, 'single');
for iSite = 1:nSites
    if isempty(cviiSpk_in1_site{iSite}), continue; end
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
fParfor = get_set_(P, 'fParfor', 1);
if fParfor
    try
        parfor iSite = 1:nSites
            S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
                cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);        
            [cvrDelta_in1{iSite}, cviNneigh_in1{iSite}] = ...
                delta_paged_site_(S_page1, S_site1, vrRho_page1, iSite);
        end %for
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iSite = 1:nSites
        S_site1 = cell2struct({cviiSpk_in1_site{iSite}, cviiSpk_out1_site{iSite}, ...
            cviiSpk2_in1_site{iSite}, cviiSpk2_out1_site{iSite}}, csName_site1, 2);        
        [cvrDelta_in1{iSite}, cviNneigh_in1{iSite}] = ...
            delta_paged_site_(S_page1, S_site1, vrRho_page1, iSite);
        fprintf('.');
    end %for
end

vrDelta1 = zeros(nSpk1, 1, 'single');
viNneigh1 = zeros(nSpk1, 1, 'int64');
for iSite = 1:nSites
    if isempty(cviiSpk_in1_site{iSite}), continue; end
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
        assert(numel(vrRho_in)==sum(vl_in), 'rho_paged_site_: vrRho_in size mismatch');
    catch
        fGpu=0;             
        fprintf(2, 'C');
%         disp(lasterr)
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
        fGpu=0; 
        fprintf(2, 'C');
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
persistent CK

knn = get_set_(P, 'knn', 30);
n_in = sum(vl_in);
nThreads = get_set_(P, 'nThreads_rho', 1024);
[CHUNK, nC_max] = deal(8, 45); % tied to search_min_drift.cu
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

if isempty(CK)
    CK = parallel.gpu.CUDAKernel('search_min_drift.ptx','search_min_drift.cu');
end
for iRetry = 1:2
    try
        nThreads = min(nThreads, CK.MaxThreadsPerBlock);
        CK.ThreadBlockSize = [nThreads, 1];          
        CK.SharedMemorySize = 4 * CHUNK * (nC_max + 1); % @TODO: update the size
        break;
    catch
        CK = parallel.gpu.CUDAKernel('search_min_drift.ptx','search_min_drift.cu');
    end
end
                
for iiAA = 1:numel(viAA)
    % determine index
    iAA1 = viAA(iiAA);
    viiAA1 = find(viDrift_in == iAA1);
    if isempty(viiAA1), continue; end
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
end
end %func


%--------------------------------------------------------------------------
% uses search_min_drift.cu
function [vrDelta_in, viNneigh_in] = search_delta_drift_(vl_in, mrFet_out, vrRho_out, viDrift_out, mlDrift1, P)

persistent CK
knn = get_set_(P, 'knn', 30);
n_in = sum(vl_in);
nThreads = get_set_(P, 'nThreads_delta', 512);
[CHUNK, nC_max] = deal(16, 45); % tied to cuda_knn_index.cu
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

if isempty(CK)
    CK = parallel.gpu.CUDAKernel('search_delta_drift.ptx','search_delta_drift.cu');
end
for iRetry = 1:2
    try
        nThreads = min(nThreads, CK.MaxThreadsPerBlock);
        CK.ThreadBlockSize = [nThreads, 1];          
        CK.SharedMemorySize = 4 * CHUNK * (nC_max + 2); % @TODO: update the size
        break;
    catch
        CK = parallel.gpu.CUDAKernel('search_delta_drift.ptx','search_delta_drift.cu');
    end
end

                
for iiAA = 1:numel(viAA)
    % determine index
    iAA1 = viAA(iiAA);
    viiAA1 = find(viDrift_in == iAA1);
    if isempty(viiAA1), continue; end
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
function vcDir = dir_up_(vcPath, nUp)
if nargin<2, nUp=0; end
vcDir = fileparts(vcPath);
for iUp=1:nUp
    vcDir = fileparts(vcDir);
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
if isempty(vcFile_prm), vcFile_prm = S_fet.P.vcFile_prm; end

[~,vcFile_,~] = fileparts(vcFile_prm);
if strcmp(vcFile_, 'sort')
    vcFile_prm_ = fullfile(dir_up_(vcFile_prm, 1), 'detect');
elseif strcmp(vcFile_, 'auto')
    vcFile_prm_ = fullfile(dir_up_(vcFile_prm, 2), 'detect');
elseif strcmpi(vcFile_, 'detect')
    vcFile_prm_ = strrep(vcFile_prm, '.prm', '');
else
    vcFile_prm_ = fullfile(fileparts(vcFile_prm), get_hash_(vcFile_prm, 'detect'), 'detect');
end
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
for iLoad = 1:nLoads
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
    
    fid1 = fopen(csFiles_fet{iLoad},'r'); 
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
% cached access to feature files to save fopen time for single-threads
function [fid, fCached] = fid_fet_cache_(vcCmd, arg1, arg2)
% usage
% -----
% fid_fet_cache_('clear')
% fid = fid_fet_cache_('set', v_fid_fet, v_fid_fet2)
% fid = fid_fet_cache_(vcFile_prm, iLoad, iFet)

persistent c_fid_load c_fid2_load
fCached = 0;

switch lower(vcCmd)
    case 'set'
        c_fid_load = arg1;
        c_fid2_load = arg2;
    case 'clear'
        fclose_(c_fid_load);
        fclose_(c_fid2_load);
        c_fid_load = [];
        c_fid2_load = [];
    otherwise
        [vcFile_prm, iLoad, iFet] = deal(vcCmd, arg1, arg2);
        vcFile_ = strrep(vcFile_prm, '.prm', '');
        switch iFet
            case 1
                if isempty(c_fid_load)
                    fid = fopen(sprintf('%s_fet_%d.irc', vcFile_, iLoad), 'r');
                else
                    fid = c_fid_load{iLoad};
                    fCached = 1;
                end
            case 2
                if isempty(c_fid2_load)
                    fid = fopen(sprintf('%s_fet2_%d.irc', vcFile_, iLoad), 'r');
                else
                    fid = c_fid2_load{iLoad};
                    fCached = 1;
                end
        end % switch
end % switch
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
% fCache_fid = false;

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
viSite2Chan = get_(P, 'viSite2Chan');
% if fCache_fid, fid_fet_cache_('clear'); end
if fParfor
    cS_detect{1} = detect_paged_save_(cS_detect{1}, P, 1);    
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
    cS_detect{1} = detect_paged_save_(cS_detect{1}, P, 1); 
    for iLoad = 2:nLoads  % change to for loop for debugging
        t_load1 = tic;
        S_load1 = vS_load(iLoad);
        mrWav_T1 = load_file_part_(vcFile, S_load1, viSite2Chan); var_size1 = var_size_(mrWav_T1);
        S_cache1 = setfield(S_cache, 'nlim_wav1', S_load1.nlim);
        cS_detect{iLoad} = detect_paged_(mrWav_T1, P, S_cache1);  mrWav_T1 = [];
        cS_detect{iLoad} = detect_paged_save_(cS_detect{iLoad}, P, iLoad);
        disp_load_(iLoad, var_size1, toc(t_load1), numel(get_(cS_detect{iLoad}, 'viTime_spk')));
    end
%     if fCache_fid, fid_fet_cache_('set', c_fid_fet, c_fid_fet2); end
end
S0 = detect_merge_(cS_detect, viOffset_load, P);
memory_detect = memory_matlab_();
runtime_detect = toc(runtime_detect);
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
function [cviSpk_site, nSites, vi_site] = vi2cell_(viSite_spk, nSites, fRow)
if nargin<2, nSites = []; end
if nargin<3, fRow=0; end 

if isempty(nSites), nSites = max(viSite_spk); end

% based on unique() function, which sorts. faster than arrayfun.
if fRow
    cviSpk_site = cell(1, nSites);
else
    cviSpk_site = cell(nSites, 1);
end
[vr, vi] = sort(viSite_spk);
vi_change = [1; find(diff(vr(:))>0)+1; numel(viSite_spk)+1];
if isempty(viSite_spk), vi_site=[]; return; end

vi_site = vr(vi_change(1:end-1));
vl_remove = vi_site < 1;
if any(vl_remove)
    vi_site(vl_remove) = [];
    vi_change(vl_remove) = [];
end
for iStep = 1:numel(vi_site)
    vr_ = vi(vi_change(iStep):vi_change(iStep+1)-1);
    if fRow
        cviSpk_site{vi_site(iStep)} = vr_(:)';
    else
        cviSpk_site{vi_site(iStep)} = vr_(:);
    end
end
vi_site = vi_site(:)';
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
function [S_detect, fid_fet, fid_fet2] = detect_paged_save_(S_detect, P, iLoad)

fKeep_fid = nargout>1;

nSites = size(P.miSites, 2);
vcFile_prm_ = strrep(P.vcFile_prm, '.prm', '');

fid_fet = fopen([vcFile_prm_, sprintf('_fet_%d.irc', iLoad)], 'w+');
S_detect.type_fet = class(S_detect.trPc_spk);
S_detect.dimm_fet = size(S_detect.trPc_spk);
S_detect.cviSpk_site = save_paged_fet_site_(...
    fid_fet, S_detect.trPc_spk, S_detect.viSite_spk, nSites);
S_detect.trPc_spk = [];
if ~fKeep_fid, fclose(fid_fet); end

if isempty(get_(S_detect, 'trPc2_spk'))
    S_detect.cviSpk2_site = cell(size(S_detect.cviSpk_site));
    fid_fet2 = [];
    S_detect.trPc2_spk = [];
else
    fid_fet2 = fopen([vcFile_prm_, sprintf('_fet2_%d.irc', iLoad)], 'w+');
    S_detect.cviSpk2_site = save_paged_fet_site_(...
        fid_fet2, S_detect.trPc2_spk, S_detect.viSite2_spk, nSites);
    S_detect.trPc2_spk = [];
    if ~fKeep_fid, fclose(fid_fet2); end
end
end %func


%--------------------------------------------------------------------------
function cviSpk_site = save_paged_fet_site_(fid_w, trFet_spk, viSite_spk, nSites)
cviSpk_site = vi2cell_(viSite_spk, nSites);
for iSite = 1:nSites
    write_bin_(fid_w, trFet_spk(:,:,cviSpk_site{iSite}));
end
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

if get_set_(P, 'fMatchedFilter_detect', 0)
    if isempty(mrPv_global)
        [viTime_spk, vrAmp_spk, viSite_spk, vrThresh_site] = detect_spikes_(mrWav_filt, vrThresh_site, vlKeep_ref, nlim_wav1, P);
        trWav_spk = get_spkwav_(mrWav_filt, viSite_spk, viTime_spk, P);
        [mrPv_global, vrD_global] = get_prinvec_(trWav_spk, P); 
    end
    S_detect = denoise_detect_(mrWav_filt, mrPv_global, nlim_wav1, vlKeep_ref, P);
    S_detect = struct_add_(S_detect, vrThresh_site, mrPv_global, vrD_global);
else
    % detect spikes or use the one passed from the input (importing)
    [viTime_spk, vrAmp_spk, viSite_spk, vrThresh_site] = detect_spikes_(mrWav_filt, vrThresh_site, vlKeep_ref, nlim_wav1, P);
    trWav_spk = get_spkwav_(mrWav_filt, viSite_spk, viTime_spk, P);
    % extract spike feaures
    if isempty(mrPv_global)
        [mrPv_global, vrD_global] = get_prinvec_(trWav_spk, P); 
    end    
    trPc_spk = project_pc_(trWav_spk, mrPv_global, P); trWav_spk = [];
    if get_set_(P, 'sort_mode', 1) == 1 && size(trPc_spk,2) > 1
        [viSite2_spk, viTime2_spk] = ...
            find_site2_spk_(trPc_spk, mrPv_global, viSite_spk, viTime_spk, P); 
        trWav2_spk = mr2tr_wav_spk2_(mrWav_filt, viSite2_spk, viTime2_spk, P);
        trPc2_spk = project_pc_(trWav2_spk, mrPv_global, P);
    else
        [viSite2_spk, trPc2_spk] = deal([]);
    end
    % return struct
    if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end
    [mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk, viSite_spk, P); 
    mrVp_spk = [];
    S_detect = makeStruct_(trPc_spk, mrPv_global, viTime_spk, vrAmp_spk, viSite_spk, ...
        mrVp_spk, trPc2_spk, viSite2_spk, vrThresh_site, mrPos_spk, vrPow_spk, vrD_global);
end
end %func


%--------------------------------------------------------------------------
function trWav_spk1 = mr2tr_wav_spk2_(mrWav_filt, viSite_spk, viTime_spk, P)

if nargin<5, mrWav_detect=[]; end

nSpks = numel(viSite_spk);
nSites = numel(P.viSite2Chan);
spkLim_wav = P.spkLim;
nSites_spk = size(P.miSites,1);
trWav_spk1 = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'like', mrWav_filt);
% nTime_search = 2;
for iSite = 1:nSites
    viiSpk1 = find(viSite_spk == iSite);
    if isempty(viiSpk1), continue; end
    viTime_spk1 = viTime_spk(viiSpk1); %already sorted by time
    viSite1 = P.miSites(:,iSite);        
    trWav_spk1(:,:,viiSpk1) = permute(mr2tr_(mrWav_filt, spkLim_wav, viTime_spk1, viSite1), [1,3,2]); %raw    
end
end %func


%--------------------------------------------------------------------------
% 6/20/2018 JJJ
% Search the min time from the range started
function [viSpk, viUpdated, viShifted] = search_local_min_(vrWav, viSpk, nSearch)
% center the wave at the min

nSearch = cast(nSearch, 'like', viSpk);
mi_ = viSpk(:)' + (-nSearch:nSearch)';
mi_(mi_<1) = 1;
mi_(mi_>numel(vrWav)) = numel(vrWav);
[~, viShift_spk1] = min(vrWav(mi_));
viShift_spk1 = cast(viShift_spk1, 'like', viSpk);
iTime0 = nSearch + 1;
viUpdated = find(viShift_spk1 ~= iTime0);
viShifted = viShift_spk1(viUpdated) - iTime0;
viSpk(viUpdated) = viSpk(viUpdated) + viShifted(:);
end %func


%--------------------------------------------------------------------------
function S_detect = denoise_detect_(mrWav_filt, mrPv_global, nlim_wav, vlKeep_wav, P)
if nargin<3, nlim_wav=[]; end
if nargin<4, vlKeep_wav=[]; end

% error('denoise_detect_: not implemented');
[nT, nSite] = size(mrWav_filt);
nPc = size(mrPv_global,2);
mad_ = @(x)x./median(abs(x),1);

% convolve with the first PC component and detect
nShift = round(mean(P.spkLim));
vrConv = mrPv_global(end:-1:1,1); % time-inverse
conv_detect_ = @(i)spikeDetectSingle_fast_(conv(mrWav_filt(:,i), vrConv,'same'), P);
[cviSpk_site, cvrSpk_site] = arrayfun_(@(i)conv_detect_(i), (1:nSite)');
[viTime_spk, vrAmp_spk, viSite_spk] = spikeMerge_(cviSpk_site, cvrSpk_site, P);
[cviSpk_site, cvrSpk_site] = deal([]); % free memory
viTime_spk = viTime_spk - nShift+1; % apply shift for PC projection

% exclude edge cases adn motion detection
if ~isempty(vlKeep_wav)
    [viTime_spk, vrAmp_spk] = select_vr_(viTime_spk, vrAmp_spk, find(vlKeep_wav(viTime_spk)));
end
if ~isempty(nlim_wav)
    vi_ = find(viTime_spk >= nlim_wav(1) & viTime_spk <= nlim_wav(2));
    [viTime_spk, vrAmp_spk, viSite_spk] = deal(viTime_spk(vi_), vrAmp_spk(vi_), viSite_spk(vi_));    
end%if

% extract waveforms and PC
trWav_spk = get_spkwav_(mrWav_filt, viSite_spk, viTime_spk, P);
trPc_spk = project_pc_(trWav_spk, mrPv_global, P);
trWav_spk = [];
% extract secondary PC
if get_set_(P, 'sort_mode', 1) == 1 && size(trPc_spk,2) > 1
    [viSite2_spk, viTime2_spk] = ...
        find_site2_spk_(trPc_spk, mrPv_global, viSite_spk, viTime_spk, P); 
    trWav2_spk = mr2tr_wav_spk2_(mrWav_filt, viSite2_spk, viTime2_spk, P);
    trPc2_spk = project_pc_(trWav2_spk, mrPv_global, P);
else
    [viSite2_spk, trPc2_spk] = deal([]);
end

% output
if ~isempty(nlim_wav), viTime_spk = viTime_spk - (nlim_wav(1)-1); end
[mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk, viSite_spk, P);
mrVp_spk = [];
S_detect = makeStruct_(viTime_spk, viSite_spk, vrAmp_spk, trPc_spk, ...
    viSite2_spk, trPc2_spk, mrPos_spk, vrPow_spk, mrVp_spk);   
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
switch get_set_(P, 'prinvec_mode', 1)
    case 1, mr1 = gather_(reshape(tr, size(tr,1), []));
    case 2, mr1 = gather_(reshape(tr(:,1,:), size(tr,1), []));
    otherwise, error('get_prinvec_: unsupported');
end
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
        catch
            disperr_('mr2tr_: failed'); 
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
function [mrWav_filt, vrWav_filt_mean] = filter_car_(mnWav, P)
%-----
% Filter
% fprintf('\tFiltering spikes (%s)...', get_(P, 'vcCommonRef')); t_fun = tic;
% if get_set_(P, 'fSmooth_spatial', 0)
%     mnWav_T = spatial_smooth_(mnWav_T, P);
% end
mrWav_filt = fft_filter(single(mnWav), P);

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
    viSite2Chan_ = get_(P, 'viSite2Chan');
    switch recording_type_(vcFile)        
        case 'mda'
            [S_mda, fid] = readmda_header_(vcFile);
            viSite2Chan_ = [];
        case 'spikeglx', [S_mda, fid] = spikeglx_header_(vcFile);      
        case 'rhd', [S_mda, fid] = rhd_header_(vcFile); 
        case 'neuroscope', [S_mda, fid] = neuroscope_header_(vcFile);
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
    case 'drift', vcDir_in = 'hybrid_synth/drift_siprobe/rec_64c_1200s_11'; 
    case 'static', vcDir_in = 'hybrid_synth/static_siprobe/rec_64c_1200s_11'; 
    case 'tetrode', vcDir_in = 'hybrid_synth/static_tetrode/rec_4c_1200s_11'; 
    case 'tetrode1', vcDir_in = 'hybrid_synth/static_tetrode/rec_4c_1200s_21'; 
    case 'tetrode2', vcDir_in = 'hybrid_synth/static_tetrode/rec_4c_1200s_31'; 
    case 'bionet', vcDir_in = 'bionet/bionet_static/static_8x_A_2A';
    case 'bionet1', vcDir_in = 'bionet/bionet_drift/drift_8x_A_2A';     
    case 'monotrode', vcDir_in = 'waveclus_synth/quiroga_difficult1/C_Difficult1_noise005';
    case 'monotrode1', vcDir_in = 'waveclus_synth/quiroga_difficult1/C_Difficult1_noise01';
    case 'monotrode2', vcDir_in = 'waveclus_synth/quiroga_difficult1/C_Difficult1_noise02';
    case 'boyden', vcDir_in = 'paired_recordings/boyden32c/419_1_7';
    case 'boyden1', vcDir_in = 'paired_recordings/boyden32c/419_1_8';
    case 'boyden2', vcDir_in = 'paired_recordings/boyden32c/915_8_1';
    case 'boyden3', vcDir_in = 'paired_recordings/boyden32c/915_10_1';
    case 'boyden4', vcDir_in = 'paired_recordings/boyden32c/915_18_1';
    case 'neuropix', vcDir_in = read_cfg_('path_neuropix_sample');
    case 'kampff', vcDir_in = 'paired_recordings/kampff/2015_09_03_Pair_9_0A';
    case 'kampff1', vcDir_in = 'paired_recordings/kampff/2015_09_03_Pair_9_0B';
    case 'kampff2', vcDir_in = 'paired_recordings/kampff/c28';
    case 'kampff3', vcDir_in = 'paired_recordings/kampff/c45';
    case 'kampff4', vcDir_in = 'paired_recordings/kampff/c46';
    case 'english', vcDir_in = '../DanEnglish/juxta_cell_curated/m139_200114_230220';
    case 'english1', vcDir_in = '../DanEnglish/juxta_cell_curated/m139_200114_222743';
    case 'english2', vcDir_in = '../DanEnglish/juxta_cell_curated/m113_191125_213423';
    case 'english3', vcDir_in = '../DanEnglish/juxta_cell_curated/m14_190326_160710_cell1';
    case 'english4', vcDir_in = '../DanEnglish/juxta_cell_curated/m15_190315_152315_cell1';        
    otherwise, error('unsupported test mode');
end
vcDir_in = fullfile(read_cfg_('path_groundtruth'), vcDir_in);
if ispc()
    vcDir_in = strrep(vcDir_in, '/', '\');    
end
end %func


%--------------------------------------------------------------------------
function [P, S_cfg] = load_default_prm_(S_cfg)
if nargin<1, S_cfg=[]; end
if isempty(S_cfg), S_cfg = read_cfg_(); end

read_struct_ = @(x)file2struct_(ircpath_(S_cfg.(x)));
P = read_struct_('default_prm');
P = struct_merge_(P, read_struct_('default2_prm'));
P_user = read_struct_('user_prm');
if ~isempty(P_user)
    P = struct_merge_(P, P_user);
    csName_user = fieldnames(P_user);
    fprintf('Default parameters were overrode from %s: %s\n', ...
        S_cfg.user_prm, sprintf('%s, ', csName_user{:}));
end
end %func


%--------------------------------------------------------------------------
function [P, vcDir_out] = makeParam_(vcDir_in, vcDir_out, vcFile_arg, fParfor)
% usage
% -----
% makeParam_(vcDir_in, vcDir_out, vcFile_arg, fParfor)
%    Provide the output directory, optional argumentts, and optional parfor
%    flag
% makeParam_(vcDir_in, myprobe.prb)
%    provide myprobe.prb file. output to `myprobe` folder under the current
%    directory
% makeParam_(vcDir_in, myprobe.prb, vcDir_out)
% makeParam_(vcDir_in, vcFile_prm, vcDir_out)

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

% parse probe file
[vcFile_prb, vcFile_prm] = deal('', '');
[~,~,vcExt_out] = fileparts(vcDir_out);
switch lower(vcExt_out)
    case '.prb'
        vcFile_prb = vcDir_out;
        if exist_dir_(vcFile_arg)
            vcDir_out = vcFile_arg;
        else
            [~, vcProbe, ~] = fileparts(vcFile_prb);
            vcDir_out = fullfile(vcDir_in, vcProbe);
        end
    case '.prm'
        vcFile_prm = vcDir_out;
        vcDir_out = fileparts(vcDir_out);
end
vcDir_out = fill_dir_out_(vcDir_in, vcDir_out);

% assume there is raw.mda, geom.csv, params.json, firings_true.mda
% group.csv
P = load_default_prm_();

% now only supporting .mda file
P.vcFile = vcFile_raw;
P.vcDir_out = vcDir_out;
vcRecordingType = recording_type_(vcFile_raw);
switch vcRecordingType
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
        if exist_file_(fullfile(vcDir_in, 'group.csv'))
            P.viShank_site = csvread(fullfile(vcDir_in, 'group.csv'));
            assert(numel(P.viShank_site)==numel(P.viSite2Chan), 'number of sites should match');
        end

        % load json file
        S_json = loadjson_(fullfile(vcDir_in, 'params.json'));
        P.sRateHz = get_set_(S_json, 'samplerate', P.sRateHz);
        P.fInverse_file = get_set_(S_json, 'spike_sign', -1) > 0;
    case {'spikeglx', 'rhd', 'neuroscope'}
        P.probe_file = vcFile_prb;
        switch vcRecordingType
            case 'spikeglx' % load neuropixels probe if meta has it
                S_meta = read_meta_spikeglx_(strrep(vcFile_raw, '.bin', '.meta')); 
                if isempty(vcFile_prb) && ~isempty(get_(S_meta, 'vcProbe'))
                    P.probe_file = fullfile([S_meta.vcProbe, '.prb']);
                end  
            case 'rhd', S_meta = rhd_header_(vcFile_raw);
            case 'neuroscope', S_meta = neuroscope_header_(vcFile_raw);
        end                
        S_prb = load_prb_(P.probe_file);
        assert(~isempty(S_prb), 'probe file must be provided');
        [P.sRateHz, P.uV_per_bit, P.vcDataType, P.nChans] = ...
            get_(S_meta, 'sRateHz', 'uV_per_bit', 'vcDataType', 'nChans');
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
if ~isempty(fParfor), P.fParfor = fParfor; end % override parfor

% derived fields
P = fill_param_(P);
if isempty(vcFile_prm)
    P.vcFile_prm = fullfile(vcDir_out, 'raw_geom.prm');
else
    P.vcFile_prm = vcFile_prm;
end

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
elseif strcmpi(vcExt_, '.dat')
    if exist_file_(fullfile(vcDir_, [vcFile_, '.xml']))
        vcType = 'neuroscope';
    elseif exist_file_(fullfile(vcDir_, 'info.rhd'))
        vcType = 'rhd';
    end
end
end %func

%--------------------------------------------------------------------------
function vcDir_out = fill_dir_out_(vcDir_in, vcDir_out)
if isempty(vcDir_out)
    vcDir_out = fullfile(vcDir_in, 'irc2');
end
mkdir_(vcDir_out);
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
S_meta = read_meta_spikeglx_(strrep(fname, '.bin', '.meta'));
vcDataType = S_meta.vcDataType;
nBytes_sample = bytesPerSample_(S_meta.vcDataType);
nBytes_data = filesize_(fname);
nSamples = nBytes_data / nBytes_sample;
nChans = S_meta.nChans;
dimm = [nChans, floor(nSamples/nChans)];
S_mda = struct('dimm', dimm, 'vcDataType', vcDataType, ...
    'nBytes_header', 0, 'nBytes_sample', nBytes_sample, ...
    'nBytes_missing', 0, 'nBytes_data', nBytes_data);

if nargout>=2, fid_r = fopen(fname,'rb'); end
end %func


%--------------------------------------------------------------------------
% intan format
function [S_mda, fid_r] = rhd_header_(vcFile_dat)
% usage
% S_mda = read_info_rhd_(fname)
% [S_mda, fid_r] = read_info_rhd_(fname)

vcFile_info = fullfile(fileparts(vcFile_dat), 'info.rhd');
assert(exist_file_(vcFile_info), 'read_info_rhd_: info.rhd does not exist');
S_rhd = import_rhd(vcFile_info);

sRateHz = S_rhd.frequency_parameters.amplifier_sample_rate;
vcDataType = 'int16';
nBytes_sample = bytesPerSample_(vcDataType);
nBytes_data = filesize_(vcFile_dat);
nSamples = nBytes_data / nBytes_sample;
nChans = numel(S_rhd.amplifier_channels);
dimm = [nChans, floor(nSamples/nChans)];
S_mda = struct('dimm', dimm, 'vcDataType', vcDataType, ...
    'nBytes_header', 0, 'nBytes_sample', nBytes_sample, ...
    'nBytes_missing', 0, 'nBytes_data', nBytes_data, ...
    'nChans', nChans, 'sRateHz', sRateHz, 'uV_per_bit', .195);

if nargout>=2, fid_r = fopen(vcFile_dat,'rb'); end
end %func


%--------------------------------------------------------------------------
% intan format
function [S_mda, fid_r] = neuroscope_header_(vcFile_dat)
% usage
% S_mda = neuroscope_header_(fname): fname is the binary file
% [S_mda, fid_r] = neuroscope_header_(fname)

[vcDir, vcFile, vcExt] = fileparts(vcFile_dat);
if strcmpi(vcExt, '.xml') % XML is provided
    vcFile_dat = fullfile(vcDir, [vcFile, '.dat']);
end
vcFile_xml = fullfile(vcDir, [vcFile, '.xml']);
assert(exist_file_(vcFile_xml), 'neuroscope_header_: .xml does not exist');
assert(exist_file_(vcFile_dat), 'neuroscope_header_: .dat does not exist');

S_xml = xml2struct(vcFile_xml); % depends on an external file
sRateHz = str2num_(S_xml.parameters.acquisitionSystem.samplingRate.Text);
nChans = str2num_(S_xml.parameters.acquisitionSystem.nChannels.Text);

nBits = str2num_(S_xml.parameters.acquisitionSystem.nBits.Text);
voltageRange = str2num_(S_xml.parameters.acquisitionSystem.voltageRange.Text);
amplification = str2num_(S_xml.parameters.acquisitionSystem.amplification.Text);
uV_per_bit = (voltageRange / amplification * 1e6) / 2 ^ nBits;
vcDataType = sprintf('int%d', nBits);

nBytes_sample = bytesPerSample_(vcDataType);
nBytes_data = filesize_(vcFile_dat);
nSamples = nBytes_data / nBytes_sample;
dimm = [nChans, floor(nSamples/nChans)];
S_mda = struct('dimm', dimm, 'vcDataType', vcDataType, ...
    'nBytes_header', 0, 'nBytes_sample', nBytes_sample, ...
    'nBytes_missing', 0, 'nBytes_data', nBytes_data, ...
    'nChans', nChans, 'sRateHz', sRateHz, 'uV_per_bit', uV_per_bit);

if nargout>=2, fid_r = fopen(vcFile_dat,'rb'); end
end %func


%--------------------------------------------------------------------------
function S_meta = read_meta_spikeglx_(vcFile_meta)
% Import SpikeGLX meta file format
% https://github.com/billkarsh/SpikeGLX/blob/fc28efc2ff464be45339180424fad378c83eda92/Markdown/Metadata_20.md                    

S_meta = [];
assert(exist_file_(vcFile_meta), sprintf('%s does not exist\n', vcFile_meta));
try
    %Read Meta
    S_sglx = meta2struct_(vcFile_meta);
    S_meta.vcDataType = get_set_(S_sglx, 'vcDataType', 'int16');
    S_meta.ADC_bits = bytesPerSample_(S_meta.vcDataType) * 8;    
    
    %convert new fields to old fields   
    if isfield(S_sglx, 'niSampRate')        
        % SpikeGLX
        S_meta.nChans = S_sglx.nSavedChans;
        S_meta.sRateHz = S_sglx.niSampRate;
        S_meta.rangeMax = S_sglx.niAiRangeMax;
        S_meta.rangeMin = S_sglx.niAiRangeMin;
        S_meta.auxGain = S_sglx.niMNGain;
        try
            S_meta.outputFile = S_sglx.fileName;
            S_meta.sha1 = S_sglx.fileSHA1;      
%             S_meta.vcProbe = 'imec2';
        catch
            S_meta.outputFile = '';
            S_meta.sha1 = [];      
        end
    elseif isfield(S_sglx, 'imSampRate')        
        % IMEC Phase 3
        S_meta.sRateHz = S_sglx.imSampRate;
        S_meta.rangeMax = S_sglx.imAiRangeMax;
        S_meta.rangeMin = S_sglx.imAiRangeMin;
        S_meta.nChans = S_sglx.nSavedChans;
        [mr_imro, vrHeader_imro] = parse_imroTbl_(get_(S_sglx, 'imroTbl'));  
        S_meta.auxGain = median(mr_imro(:,4));
        S_meta.auxGain_lfp = median(mr_imro(:,5));          
        if isfield(S_sglx, 'imProbeOpt')
            % Neuropix 3A            
            S_meta.vcProbe = sprintf('imec3_opt%d', S_sglx.imProbeOpt);            
            S_meta.ADC_bits = 10;
        else            
            probe_type = vrHeader_imro(1);
            switch probe_type
                case 0 % Neuropix 3B2
                    S_meta.vcProbe = 'imec3B2';  
                    S_meta.ADC_bits = 10;
                case [21, 24] % Neuropix Phase 2
                    S_meta.auxGain = 80; %fixed
                    S_meta.ADC_bits = 14;
                    error('Neuropixels 2.0 is not supported');
                otherwise
                    error('read_meta_spikeglx_: unsupported probe: %d', probe_type);
            end            
        end
    else
        error('read_meta_spikeglx_: unsupported format');
    end
    
     %number of bits of ADC [was 16 in Chongxi original]
    try
        S_meta.uV_per_bit = ((S_meta.rangeMax-S_meta.rangeMin)/(2^S_meta.ADC_bits))/S_meta.auxGain * 1e6;  %uVolts
    catch
        S_meta.uV_per_bit = 1;
    end    
catch
    error('Parsing error: %s\n\t%s\n', vcFile_meta, lasterr());
end
end %func


%--------------------------------------------------------------------------
function [mrBody, vrHeader] = parse_imroTbl_(vc_imroTbl)
try
    vnIMRO = strrep(strrep(vc_imroTbl(2:end-1),')(',';'),' ',',');
    iBody_ = find(vnIMRO==';',1,'first');
    mrBody = str2num(sprintf('[%s]',vnIMRO(iBody_+1:end)));
    vrHeader = str2num(sprintf('[%s]',vnIMRO(1:iBody_-1)));
catch
    [mrBody, vrHeader] = deal([]);
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function S = meta2struct_(vcFile)
% Convert text file to struct
S = struct();
if ~exist_file_(vcFile, 1), return; end

fid = fopen(vcFile, 'r');
mcFileMeta = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', false);
fclose(fid);
csName = mcFileMeta{1};
csValue = mcFileMeta{2};
for i=1:numel(csName)
    vcName1 = strtrim(csName{i});
    if vcName1(1) == '~', vcName1(1) = []; end
    vcValue1 = csValue{i};
    try   
        if isempty(vcValue1)
            eval(sprintf('%s = [];', vcName1));
        elseif all(vcValue1([1,end])=='''')
            eval(sprintf('%s = %s;', vcName1, vcValue1));
        else
            eval(sprintf('%s = ''%s'';', vcName1, vcValue1));
        end
        eval(sprintf('num = str2double(%s);', vcName1));
        if ~isnan(num)
            eval(sprintf('%s = num;', vcName1));
        end
        eval(sprintf('S = setfield(S, ''%s'', %s);', vcName1, vcName1));
    catch
        fprintf(2, 'error: %s = %s\n', csName{i}, csValue{i});
    end
end
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
if P.maxDist_site_spk_um<=0
    P.maxDist_site_spk_um = inf; 
elseif P.maxDist_site_spk_um < P.maxDist_site_um
    P.maxDist_site_spk_um = P.maxDist_site_um;
end

% set frequency
freq_min = get_set_(S_txt, 'freq_min', []);
freq_max = get_set_(S_txt, 'freq_max', []);
if ~isempty(freq_min) && ~isempty(freq_max)
    P.freqLim = [freq_min, freq_max]; 
elseif isempty(freq_min) && ~isempty(freq_max)
    P.freqLim = [nan, freq_max]; 
elseif ~isempty(freq_min) && isempty(freq_max)
    P.freqLim = [freq_min, nan]; 
elseif isempty(freq_min) && isempty(freq_max)
    P.freqLim = [nan, nan];
end

% spkLim_ms
clip_pre = get_set_(S_txt, 'clip_pre', []);
clip_post = get_set_(S_txt, 'clip_post', []);
if ~isempty(clip_pre) && ~isempty(clip_post)
    P.spkLim_ms = [-abs(clip_pre), abs(clip_post)];
end

% disable filter
fFilter = logical_(get_set_(S_txt, 'filter', 1));
if ~fFilter, P.freqLim = [nan, nan]; end

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
    {'qqFactor', 'nPc_spk', 'maxWavCor', 'uV_per_bit'}, 1);

% String parameters
P = struct_copyas_(P, S_txt, {'filter_type', 'feature_type'}, {'vcFilter', 'vcFet'});

% same name
P = struct_copyas_(P, S_txt, ...
    {'knn', 'batch_sec_drift', 'step_sec_drift', 'min_count', 'nSites_whiten', ...
    'fft_thresh', 'delta_cut', 'fft_thresh_low', 'post_merge_mode', 'sort_mode', ...
    'merge_thresh_cc', 'merge_overlap_thresh'});

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
if step_sec_drift>0 && batch_sec_drift>0
    P.nTime_batch = max(round(batch_sec_drift / step_sec_drift), 1);
    P.nTime_drift = max(round(recording_duration_(P) / step_sec_drift), 1);
else
    P.nTime_batch = 1;
    P.nTime_drift = 1;
end
% if isempty(get_(P, 'nTime_batch'))
%     P.nTime_batch = max(round(batch_sec_drift / step_sec_drift), 1);
%     fprintf('\tnTime_batch = %d (batch_sec_drift = %0.1f s)\n', P.nTime_batch, batch_sec_drift);
% end
% if isempty(get_(P, 'nTime_drift'))
%     try
%         P.nTime_drift = max(round(recording_duration_(P) / step_sec_drift), 1);
%         fprintf('\tnTime_drift = %d (step_sec_drift = %0.1f s)\n', P.nTime_drift, step_sec_drift);
%     catch
%         P.nTime_drift = 64;
%     end
% end

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
    if contains(vcDataType, {'uint16', 'uint32'})
        mnWav = uint2int_(mnWav);
    end
    if ischar(vcFile), fclose(fid); end
catch
    disp(lasterr());
end
end %func


%--------------------------------------------------------------------------
function mn = uint2int_(mn)
if isa(mn, 'uint16')
    mn = int16(single(mn)-2^15);
elseif isa(mn, 'uint32')
    mn = int32(double(mn)-2^31);
end
end


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
        mlDrift = dist_mask_ml_(mrCount_drift, mlMask, nTime_batch);
    else
        mrDist_drift = squareform(pdist(mrCount_drift'));
        [~, miSort_drift] = sort(mrDist_drift, 'ascend');
        miSort_drift = miSort_drift(1:nTime_batch,:);
        mlDrift = mi2ml_drift_(miSort_drift); %gpuArray_(mi2ml_drift_(miSort_drift), P.fGpu);
    end
    if false % widen the diagonal
        n_ = nTime_drift;
        mlDrift(2:n_+1:n_^2) = true;
        mlDrift(n_+2:n_+1:n_^2) = true;
    end
    if read_cfg_('fPlot_drift')
        figure; imagesc(mrDist_drift); set(gcf,'Name', P.vcFile_prm);
%         figure; imagesc(mrSort_drift); set(gcf,'Name', P.vcFile_prm);
        hold on; plot([0, size(mrSort_drift,1)], repmat(nTime_batch,1,2), 'r-');
    end    
end
S_drift = makeStruct_(nTime_drift, viDrift_spk, mlDrift, viLim_drift);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function ml = dist_mask_ml_(mrF, mlMask, nneigh)
% return euclidean distance for masked portions

ml = false(size(mrF,2));
for iCol=1:size(mrF,2)
    viMask1 = find(mlMask(:,iCol));
    vrD1 = sqrt(sum((mrF(:,iCol) -  mrF(:,viMask1)).^2, 1));
    [~, viSrt1] = sort(vrD1, 'ascend');
    ml(viMask1(viSrt1(1:nneigh)), iCol) = true;
end
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
            cvi_peak{iPeak} = find(any(ismember(miKnn_peak, miKnn_peak(:,iPeak))))';      
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iPeak = 1:numel(viSpk_peak)        
        cvi_peak{iPeak} = find(any(ismember(miKnn_peak, miKnn_peak(:,iPeak))))';
    end
end
[viClu_spk(viSpk_peak), viiPeak] = cell2map_(cvi_peak, [], fParfor);
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
% compute cluster isolation distance
function [vrScore_clu, viClu_nn_clu] = isolation_score_(S_auto, miKnn_spk)
% usage
% [vrScore_clu, viClu_nn_clu] = isolation_score_(S_auto, viSite_spk, miKnn_spk)
% [vrScore_clu, viClu_nn_clu] = isolation_score_(S_auto, viSite_spk, miKnn_spk)

t_fun=tic;

[viClu_spk, P] = get_(S_auto, 'viClu', 'P');
try
   [cviSpk_clu, nClu] = deal(S_auto.cviSpk_clu, S_auto.nClu);
catch
    [cviSpk_clu, nClu] = vi2cell_(viClu_spk);
end
[vrScore_clu, viClu_nn_clu] = deal(nan(nClu,1));
try
    get_miClu_ = @(x)viClu_spk(miKnn_spk(:,cviSpk_clu{x}));
    for iClu1 = 1:nClu
        miClu1 = get_miClu_(iClu1);
        n11 = sum(miClu1(:) == iClu1);    
        iClu2 = mode(miClu1(miClu1~=iClu1 & miClu1>0));
        n12 = sum(miClu1(:) == iClu2);
        if isempty(iClu2) % no matching cluster found
            vrScore_clu(iClu1) = 1;
        else
            vrScore_clu(iClu1) = (n11) / (n11+n12);
            viClu_nn_clu(iClu1) = iClu2;
        end
    end
    fprintf('\tisolation_score_: took %0.1fs\n', toc(t_fun));
catch
    fprintf(2, 'isolation_score_: failed\n');
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
function [viMapClu_new, viUniq_, viMapClu] = cell2map_(cvi_clu, viClu_remove, fParfor)
if nargin<2, viClu_remove = []; end
if nargin<3, fParfor=0; end

nRepeat = 10;

if ~isempty(viClu_remove)
    cvi_clu = cellfun_(@(x)setdiff(x, viClu_remove), cvi_clu);
    cvi_clu(viClu_remove) = {[]};
end
nClu = numel(cvi_clu);
fprintf('\tcell2map_: '); t_fun=tic;
for iRepeat = 1:nRepeat
    cvi_clu1 = cell(size(cvi_clu));
    viClu_update = find(cellfun(@(x)numel(x)>1, cvi_clu));
    viClu_update = viClu_update(:)';
    for iClu = viClu_update
        viClu1 = cvi_clu{iClu};
        viClu1 = cat(1, viClu1(:), cvi_clu{viClu1});
        viClu1 = int32(unique(viClu1));
        cvi_clu1(viClu1) = {viClu1};
    end    
    if all(cellfun(@numel, cvi_clu1) == cellfun(@numel, cvi_clu))
        break;
    else
        cvi_clu = cvi_clu1;
    end
    fprintf('.');
end

nClu = numel(cvi_clu);
viMapClu = 1:nClu;
for iClu = 1:nClu
    iClu1 = min(cvi_clu{iClu});    
    if ~isempty(iClu1)
        viMapClu(iClu) = iClu1; 
    end
end
if ~isempty(viClu_remove)
    viMapClu(viClu_remove) = 0;
end

% Compact the map so the index doesn't have a gap
viUniq_ = setdiff(unique(viMapClu), 0);
viMap_(viUniq_) = 1:numel(viUniq_);
vlUpdate = viMapClu>0;
viMapClu_new = zeros(size(viMapClu));
viMapClu_new(vlUpdate) = viMap_(viMapClu(vlUpdate));
fprintf(' nRepeat=%d, took %0.1fs\n', iRepeat, toc(t_fun));
end %func


%--------------------------------------------------------------------------
function cvi_clu1 = func_parfor_(fh, cvi_clu1, fParfor)
if nargin<2, fParfor = 1; end
nClu = numel(cvi_clu1);
if fParfor
    try
        parfor iClu=1:nClu
            cvi_clu1{iClu} = fh(cvi_clu1{iClu});
        end
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iClu=1:nClu
        cvi_clu1{iClu} = fh(cvi_clu1{iClu});
    end
end
end %func


%--------------------------------------------------------------------------
% do not map self it not found
function [viMapClu_new, viUniq_, viMapClu] = ml2map_(mlClu)
nClu = size(mlClu,1);
% mlClu = set_diag_(mlClu | mlClu', true(nClu,1));
% set diagonal to true

mlClu = mlClu | mlClu'; % make symmetric
viClu = find(any(mlClu,1));
for iClu = viClu % very slow, how to speed it up?
    vl_ = mlClu(:,iClu);
    mlClu(vl_,vl_) = true;
end

viMapClu = zeros(1, nClu);
for iClu = 1:nClu 
    iClu1 = find(mlClu(:,iClu), 1, 'first');
    if ~isempty(iClu1), viMapClu(iClu) = iClu1; end
end

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
    fprintf('Saving a struct to %s: ', vcFile); t1=tic;
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
    fprintf('took %0.1fs.\n', toc(t1));
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
try
    irc2('compile');
catch
    fprintf(2, 'CUDA compilation failed\n');
end
irc2('mcc');
S_cfg = read_cfg_();
[sf2_path, sf2_test_path, sf2_python_path, sf2_test_script, sf2_docker_name] = ...
    get_(S_cfg, 'spikeforest2_irc_path', 'spikeforest2_test_path', ...
        'spikeforest2_python_path', 'spikeforest2_test_script', 'spikeforest2_docker_name');

% build successful, now copy and edit docker files
vc_irc_ver = sprintf('%s:%s', sf2_docker_name, version_());
save_cs_(fullfile(sf2_path, 'build_docker.sh'), {'#!/bin/bash','',['docker build -t ',vc_irc_ver,' .']});
save_cs_(fullfile(sf2_path, 'push_docker.sh'), {'#!/bin/bash','',['docker push ', vc_irc_ver]});
copyfile('run_irc', sf2_path);
system(sprintf('(cd %s && ./build_docker.sh && ./push_docker.sh)', sf2_path), '-echo');
system(sprintf('(cd %s && %s %s)', sf2_test_path, sf2_python_path, sf2_test_script), '-echo');
end %func


%--------------------------------------------------------------------------
% 4/17/18 JJJ: documentation and testing
function varargout = call_(vcFunc, cell_Input)
% S_out = call_(vcFunc, cell_Input, nOutput)
% varargout = call_(vcFunc, cell_Input)

if vcFunc(end) ~= '_', vcFunc = [vcFunc, '_']; end

nOutput = nargout();
try
    switch nOutput
        case 0, feval(vcFunc, cell_Input{:});
        case 1, varargout{1} = feval(vcFunc, cell_Input{:});
        case 2, [varargout{1}, varargout{2}] = feval(vcFunc, cell_Input{:});
        case 3, [varargout{1}, varargout{2}, varargout{3}] = feval(vcFunc, cell_Input{:});
        case 4, [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = feval(vcFunc, cell_Input{:});
        case 5, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}] = feval(vcFunc, cell_Input{:});
        case 6, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}] = feval(vcFunc, cell_Input{:});
        case 7, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}] = feval(vcFunc, cell_Input{:});
        case 8, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}, varargout{8}] = feval(vcFunc, cell_Input{:});
        otherwise, error('call_: too many output');
    end %switch
catch E
    % try irc 
    switch nOutput
        case 0, irc('call', vcFunc, cell_Input);
        case 1, varargout{1} = irc('call', vcFunc, cell_Input);
        case 2, [varargout{1}, varargout{2}] = irc('call', vcFunc, cell_Input);
        case 3, [varargout{1}, varargout{2}, varargout{3}] = irc('call', vcFunc, cell_Input);
        case 4, [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = irc('call', vcFunc, cell_Input);
        case 5, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}] = irc('call', vcFunc, cell_Input);
        case 6, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}] = irc('call', vcFunc, cell_Input);
        case 7, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}] = irc('call', vcFunc, cell_Input);
        case 8, [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}, varargout{8}] = irc('call', vcFunc, cell_Input);
        otherwise, error('call_: too many output');
    end %switch
end
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
function [viTime_spk, vrAmp_spk, viSite_spk, vrThresh_site] = ...
    detect_spikes_(mrWav_filt, vrThresh_site, vlKeep_wav, nlim_wav, P)

% fMerge_spk = 1;
fMerge_spk = get_set_(P, 'fMerge_spk', 1);

[n1, nSites, ~] = size(mrWav_filt);
[cviSpk_site, cvrSpk_site] = deal(cell(nSites,1));
if isempty(vrThresh_site)    
    vrThresh_site = mr2thresh_(mrWav_filt, P);
end

for iSite = 1:nSites
    [viSpk1, vrSpk1] = spikeDetectSingle_fast_(mrWav_filt(:,iSite), P, vrThresh_site(iSite));   
    
    % Reject global mean
    if isempty(vlKeep_wav)
        cviSpk_site{iSite} = viSpk1;
        cvrSpk_site{iSite} = vrSpk1;        
    else
        [cviSpk_site{iSite}, cvrSpk_site{iSite}] = select_vr_(viSpk1, vrSpk1, find(vlKeep_wav(viSpk1)));
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

% reject spikes within the overlap region
if ~isempty(nlim_wav)
    viiKeep_spk = find(viTime_spk >= nlim_wav(1) & viTime_spk <= nlim_wav(2));
    [viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)x(viiKeep_spk), viTime_spk, vrAmp_spk, viSite_spk);    
end%if
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
MAX_SAMPLE_QQ = 2^16; %300000; 
% fSpikeRefrac_site = 0;
if nargin < 3, thresh1 = []; end
if nargin < 2, P = struct('spkThresh', [], 'qqFactor', 5); end
if ~isempty(get_(P, 'spkThresh')), thresh1 = P.spkThresh; end

if thresh1==0, [viSpk1, vrSpk1] = deal([]); return; end % bad site
if isempty(thresh1)  
    vr_ = subsample_vr_(vrWav1, MAX_SAMPLE_QQ);
    thresh1 = median(abs(vr_ - median(vr_))); vr_=[];
%     thresh1 = median(abs(subsample_vr_(vrWav1, MAX_SAMPLE_QQ)));
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
function run_spikeforest2_(vcSorter, vcDir_in, vcDir_out, vcArg, fValidate_sf2)
if nargin<4, vcArg = ''; end
if nargin<5, fValidate_sf2 = []; end

S_cfg = read_cfg_();

if isempty(fValidate_sf2)
    fValidate_sf2 = get_(S_cfg, 'spikeforest2_validate'); 
end
csSorters_sf2 = {'mountainsort4', 'ironclust', 'kilosort2', 'kilosort', 'spykingcircus', 'herdingspikes2', 'tridesclous', 'klusta', 'waveclus', 'jrclust'};

t_fun=tic;
if isempty(vcDir_out)
    vcDir_out = fullfile(vcDir_in, vcSorter);    
end
mkdir_(vcDir_out);

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
    job_timeout = get_(S_cfg, 'spikeforest2_timeout');
    sf2_params = {'algorithm', lower(vcSorter), 'recording_path', vcDir_in, ...
        'sorting_out', vcFile_firings, 'job_timeout', job_timeout};
    fContainer = get_(S_cfg, 'spikeforest2_use_container');   
    if fContainer == 1
        sf2_params = {sf2_params{:}, 'container', 'default'};
        fprintf('irc2: spikeforest2_use_container=True\n');
    elseif fContainer == 0
        sf2_params = {sf2_params{:}, 'container', py.NoneType};
        fprintf('irc2: spikeforest2_use_container=False\n');
    end
    if ~isempty(vcArg)
        sf2_params = {sf2_params{:}, 'params', sf2_params_(vcArg, vcSorter)};        
    end    
    py.spikeforest2.experimental.sort(pyargs(sf2_params{:}));
    fprintf('%s: wrote to %s, took %0.1fs\n', vcSorter, vcFile_firings, toc(t_fun));
else
    fprintf(2, '%s is not currently supported\n', vcSorter);
    return;
end

% validate
vcFile_true = fullfile(vcDir_in, 'firings_true.mda');
if exist_file_(vcFile_true) && fValidate_sf2
    fPlot_gt = read_cfg_('fPlot_gt');
    vcFile_raw = fullfile(vcDir_in, 'raw.mda');
%     S_score = irc('validate-mda', vcFile_true, vcFile_firings, vcFile_raw, fPlot_gt); % assume that groundtruth file exists
    S_score = compare_mda_(vcFile_true, vcFile_firings);
    S_score_plot_(S_score);
    struct_save_(S_score, fullfile(vcDir_out, 'raw_geom_score.mat'), 1);
end
end %func


%--------------------------------------------------------------------------
function fCreatedDir = mkdir_(vcDir)
% make only if it doesn't exist. provide full path for dir
fCreatedDir = exist_dir_(vcDir);
if ~fCreatedDir
    try
        mkdir(vcDir); 
    catch
        fCreatedDir = 0;
    end
end
end %func


%--------------------------------------------------------------------------
function params = sf2_params_(vcArg, vcSorter)

params = py.dict;
if isempty(vcArg), return; end

if isstruct(vcArg)
    S_arg = vcArg;
elseif ischar(vcArg)
    if ~exist_file_(vcArg), error('%s does not exist', vcArg); end
    S_arg = meta2struct_(vcArg); % name = value type
end

% global parameters
cs_common_double = {'detect_threshold', 'detect_sign'};
params = pydict_add_(params, S_arg, cs_common_double, 'double');

% sorter-specific parameters
[cs_sorter_double, cs_sorter_char, cs_sorter_logical, cs_sorter_int, cs_sorter_pylist] = deal({});
switch lower(vcSorter)
    case 'ironclust'
        cs_sorter_double = {'freq_min', 'adjacency_radius_out', 'merge_thresh', ...
            'fft_thresh', 'knn', 'min_count', 'delta_cut', 'pc_per_chan', ...
            'batch_sec_drift', 'step_sec_drift', 'freq_max', 'clip_pre', 'clip_post', 'adjacency_radius'};
        cs_sorter_char = {'common_ref_type'};
        cs_sorter_logical = {'whiten', 'fGpu', 'filter'};
    case 'mountainsort4'
        cs_sorter_double = {'freq_min', 'detect_interval', 'freq_max', 'adjacency_radius'};
        cs_sorter_logical = {'curation', 'whiten', 'filter'};
        cs_sorter_int = {'clip_size'};
    case 'kilosort2'
        cs_sorter_double = {'freq_min', 'sigmaMask', 'nPCs', 'minFR', 'Nt', 'preclust_threshold'};
        cs_sorter_logical = {'car'};
        cs_sorter_pylist = {'projection_threshold'};
    case 'spykingcircus'
        cs_sorter_double = {'auto_merge', 'template_width_ms', 'whitening_max_elts', 'clustering_max_elts', 'adjacency_radius'};
        cs_sorter_logical = {'filter', 'merge_spikes'};                
    case 'klusta'
        cs_sorter_double = {'threshold_strong_std_factor', 'threshold_weak_std_factor'};
        cs_sorter_int = {'n_features_per_channel', 'num_starting_clusters', 'extract_s_before', 'extract_s_after'};
    case 'jrclust'        
    case 'waveclus'        
    case 'kilosort'        
    case 'tridesclous'
    otherwise, error('unsupported sorter');
end
params = pydict_add_(params, S_arg, cs_sorter_double, 'double');
params = pydict_add_(params, S_arg, cs_sorter_char, 'char');
params = pydict_add_(params, S_arg, cs_sorter_logical, @(x)logical_(x));
params = pydict_add_(params, S_arg, cs_sorter_int, 'int32');
params = pydict_add_(params, S_arg, cs_sorter_pylist, @(x)py.list(x));
end %func


%--------------------------------------------------------------------------
function x = logical_(x)
if isempty(x)
    x = false;
elseif ischar(x)
    if strcmpi(x, 'true') || strcmpi(x, '1')
        x = true;
    elseif strcmpi(x, 'false') || strcmpi(x, '0')
        x = false;
    else
        error('logical_: undefined: %s', x);
    end
else
    x = logical(x);
end
end %func


%--------------------------------------------------------------------------
function pydict = pydict_add_(pydict, S, csNames, fh_cast)

if nargin<4, fh_cast = []; end
if isempty(csNames), return; end

if isempty(fh_cast)
    fh_cast = @(x)x;
elseif ischar(fh_cast)
    fh_cast = @(x)cast(x, fh_cast);
end

if ischar(csNames), csNames = {csNames}; end
for iCell = 1:numel(csNames)
    vcName_ = csNames{iCell};
    if isfield(S, vcName_)
        pydict{vcName_} = fh_cast(S.(vcName_));
    end
end
end %func


%--------------------------------------------------------------------------
% Compile CUDA codes for IronClust
% 10/5/17 JJJ: Error messages converted to warning
% 7/26/17 JJJ: Code cleanup and testing
function fSuccess = compile_cuda_(csFiles_cu, vcSystem)
% usage
% -----
% compile_cuda_(csFiles_cu): compile cuda
% compile_cuda_(csFiles_cu, '1'): compile cuda for local system

if nargin<2, vcSystem=''; end

S_cfg = read_cfg_();
if nargin<1 || isempty(csFiles_cu)     
    csFiles_cu = S_cfg.csFiles_cuda; %version 3 cuda
elseif ischar(csFiles_cu)
    csFiles_cu = {csFiles_cu};
end
fSystem = ifeq_(isempty(vcSystem), 0, str2num_(vcSystem));

t1 = tic;
disp('Compiling CUDA codes...');
fSuccess = 1;
S_gpu = gpuDevice(1);
vcPath_nvcc = get_set_(S_cfg, 'vcPath_nvcc', get_nvcc_path_(S_gpu));
if isempty(vcPath_nvcc), fSuccess = 0; return; end
if fSuccess
    try % determine compute capability version
        vc_cuda_arch_target = get_set_(S_cfg, 'cuda_arch', 'sm_35');
        sm_ver_target = str2double(strrep(vc_cuda_arch_target, 'sm_', '')) / 10;     
        sm_ver_system = str2double(S_gpu.ComputeCapability);
        vc_cuda_arch_system = sprintf('sm_%d', round(sm_ver_system*10));
        if fSystem
            vc_cuda_arch = vc_cuda_arch_system;
        else
            % support less than 4GB GPU memory
            if sm_ver_system < sm_ver_target
                vc_cuda_arch = 'sm_30';
            else
                vc_cuda_arch = vc_cuda_arch_target;
            end
        end
    catch
        fSuccess = 0;
    end
    for i=1:numel(csFiles_cu)
        vcFile_ = ircpath_(csFiles_cu{i});
        for iRetry = 1:2
            vcCmd1 = sprintf('%s -ptx -m 64 -arch %s "%s"', vcPath_nvcc, vc_cuda_arch, vcFile_);
            fprintf('\t%s\n\t', vcCmd1);            
            try          
                status = system(vcCmd1);
                fSuccess1 = (status==0);        
            catch
                fprintf('\tWarning: CUDA could not be compiled: %s\n', vcFile_); 
                fSuccess1 = 0;
            end
            if fSuccess1, break; end
            fSuccess = fSuccess && fSuccess1;
            vc_cuda_arch = vc_cuda_arch_target;
        end
    end
end
if ~fSuccess
    disp_cs_({'    Warning: CUDA could not be compiled but IronClust may work fine.', 
    sprintf('    If not, install CUDA toolkit v%0.1f and run "irc install".', S_gpu.ToolkitVersion),
    '    If you encounter error:  nvcc fatal   : Microsoft Visual Studio configuration file ''vcvars64.bat'' could not be found...', 
    '      1. Copy $VS12/VC/bin/x86_amd64 to $VS12/VC/bin/amd64', 
    '      2. Rename $VS12/VC/bin/amd64/vcvarsx86_amd64.bat to vcvars64.bat'});
end

% clear cuda code persistent cache
clear cuda_knn_ search_knn_drift_ search_delta_drift_ 

fprintf('\tFinished compiling, took %0.1fs\n', toc(t1));
end %func


%---------------------------------------------------------------------------
function out = ifeq_(if_, true_, false_)
if (if_)
    out = true_;
else
    out = false_;
end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile_full = ircpath_(vcFile, fConditional)
% make it a irc path
% Add a full path if the file doesn't exist in the current folder
% 

if nargin<1, vcFile = ''; end
if nargin<2, fConditional = 0; end

ircpath = fileparts(mfilename('fullpath'));
if fConditional
    if exist_file_(vcFile)
        vcFile_full = vcFile;
        return;
    end
end
vcFile_full = [ircpath, filesep(), vcFile];
end %func


%--------------------------------------------------------------------------
function disp_cs_(cs)
% display cell string
cellfun(@(s)fprintf('%s\n',s), cs);
end %func


%--------------------------------------------------------------------------
function vcPath_nvcc = get_nvcc_path_(S_gpu)
if ispc()    
    csDir = sub_dir_('C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\');
    if isempty(csDir), vcPath_nvcc=[]; return; end
    vrVer = cellfun(@(x)str2num(x(2:end)), csDir);
    [~,imin] = min(abs(vrVer-S_gpu.ToolkitVersion));    
    vcPath_nvcc = sprintf('"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v%0.1f\\bin\\nvcc"', vrVer(imin));    
else
    vcPath_nvcc = read_cfg_('nvcc_path');
    if ~exist_file_(vcPath_nvcc)
        vcPath_nvcc = '/usr/local/cuda/bin/nvcc'; % try default installation location
    end
end
end %func


%--------------------------------------------------------------------------
function csDir = sub_dir_(vc, fFullPath)
if nargin<2, fFullPath = 0; end
S_dir = dir(vc);
[csDir, csDir_up] = deal({S_dir.name}, {S_dir.folder});
vlKeep = ~ismember(csDir, {'.', '..'}) & [S_dir.isdir];
csDir = csDir(vlKeep);
if fFullPath
    csDir = cellfun_(@(x,y)fullfile(x,y), csDir_up(vlKeep), csDir);
end
end %func


%--------------------------------------------------------------------------
% 8/7/2018 JJJ
function flag = exist_dir_(vcDir)
if isempty(vcDir)
    flag = 0;
elseif ischar(vcDir)
    S_dir = dir(vcDir);
    if isempty(S_dir)
        flag = 0;
    else
        flag = sum([S_dir.isdir]) > 0;
    end
else
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
% loads mda format
function S_firings = load_firings_(vcFile, sRateHz)
if nargin<2, sRateHz = []; end
if isempty(sRateHz), sRateHz = 30000; end

% load firings.mda
mr_ = readmda_(vcFile);
[viSite_spk, viTime_spk, viClu_spk] = deal(int32(mr_(1,:))', int64(mr_(2,:)'), int32(mr_(3,:)'));
mr_ = [];

cviSpk_clu = vi2cell_(viClu_spk);
vnSpk_clu = cellfun(@(x)numel(x), cviSpk_clu);
cviTime_clu = cellfun_(@(x)viTime_spk(x), cviSpk_clu);
cviSite_clu = cellfun_(@(x)viSite_spk(x), cviSpk_clu);
nClu_exist = sum(vnSpk_clu>0);

nSites = numel(unique(viSite_spk));
nClu = numel(cviSpk_clu);
vrDur_clu = cell2mat_(cellfun_(@(x)double(range(x))/sRateHz, cviTime_clu), nan);
vrRate_clu = arrayfun(@(x)vnSpk_clu(x) / vrDur_clu(x), 1:nClu); % .17 hz rate
duration_sec = double(range(viTime_spk)) / sRateHz; % .916 hr median, 35 hr max

S_firings = makeStruct_(sRateHz, nClu, nClu_exist, nSites, vcFile, duration_sec, ...
    viSite_spk, viTime_spk, viClu_spk, ...
    cviSpk_clu, cviTime_clu, cviSite_clu, vrDur_clu, vnSpk_clu, vrRate_clu);
end %func


%--------------------------------------------------------------------------
% 2020/jan/23, run parameter optimizer for ironclust
function optimize_prmset_(vcDir_rec, vcFile_prmset, fPreview)
% usage
% -----
% optimize_prmset_(vcFile_list)
% optimize_prmset_(vcDir_rec)
% optimize_prmset_(..., vcFile_prmset)
% optimize_prmset_(..., vcFile_prmset, vcFile_out)

fParfor_rec = 0; % disable parfor on individual recording

if nargin<3, fPreview=0; end
fPreview = logical_(fPreview);

% if fDebug, fParfor = 0; end
S_cfg = read_cfg_();
fUse_cache = get_set_(S_cfg, 'fUse_cache_optimize', 0);
csSorter_gpu = {'irc2', 'kilosort2', 'kilosort', 'jrclust'};
P_prmset = makeStruct_(fUse_cache, csSorter_gpu);
P_prmset.fParfor = fParfor_rec;

if nargin<2, vcFile_prmset=''; end
vcFile_out = ''; 
if isempty(vcFile_prmset), vcFile_prmset = S_cfg.ironclust_prmset; end
assert(exist_file_(vcFile_prmset), '.prmset file must exist');

[~,vcPostfix_] = fileparts(vcFile_prmset); 
vcSorter = infer_sorter_(vcPostfix_);
if exist_file_(vcDir_rec)   
    vcFile_out = strrep(vcDir_rec, '.txt', sprintf('_scores_prmset_%s.mat', vcPostfix_));
    csDir_rec = load_batch_(vcDir_rec);
    vcDir_rec = fileparts(vcDir_rec);
elseif exist_dir_(vcDir_rec)
    csDir_rec = {vcDir_rec};
end
if isempty(vcFile_out)
    vcFile_out = fullfile(vcDir_rec, sprintf('scores_prmset_%s.mat', vcPostfix_));
end


if ~fUse_cache, optimize_clear_(vcDir_rec,vcFile_prmset); end
t_fun = tic;
if isempty(csDir_rec)
    csDir_rec = sub_dir_(vcDir_rec, 1);
end
nRec = numel(csDir_rec);
S_prmset = file2struct_ordered_(vcFile_prmset);
[csName_prm, cVal_prm] = deal(fieldnames(S_prmset), struct2cell(S_prmset));
nPrmset = prod(cellfun(@numel, cVal_prm));
ccScore_prmset_rec = cell(nRec*nPrmset, 1);
% remove lock files

nRuns = numel(ccScore_prmset_rec);
parfor iRun = 1:nRuns
    [iRec, iPrmset] = ind2sub([nRec,nPrmset], iRun);
    cVal_prm1 = permute_prm_(cVal_prm, iPrmset);
    ccScore_prmset_rec{iRun} = load_score_prmset_(...
        vcSorter, csDir_rec{iRec}, csName_prm, cVal_prm1, P_prmset, iPrmset);
end
% rebalance the parfor loop            
viRun1 = find(cellfun(@(x)isempty(x), ccScore_prmset_rec));
ccScore_prmset_rec1 = cell(size(viRun1));
nRuns_loaded = nRuns - numel(viRun1);
fprintf(2, 'Loaded %d/%d (%0.1f%%) from cache\n', nRuns_loaded, nRuns, nRuns_loaded/nRuns*100);
if ~fPreview
    remove_lock_(csDir_rec, 1);
    parfor iRun1 = 1:numel(viRun1)
        [iRec, iPrmset] = ind2sub([nRec,nPrmset], viRun1(iRun1));
        cVal_prm1 = permute_prm_(cVal_prm, iPrmset);
        ccScore_prmset_rec1{iRun1} = score_prmset_(...
            vcSorter, csDir_rec{iRec}, csName_prm, cVal_prm1, P_prmset, iPrmset);
    end      
    ccScore_prmset_rec(viRun1) = ccScore_prmset_rec1;
    ccScore_prmset_rec1 = {}; % clear
end
ccScore_prmset_rec = reshape(ccScore_prmset_rec, [nRec,nPrmset]);
t_fun = toc(t_fun);
fprintf('took %0.1fs\n', t_fun);

% display
S_prmset_rec = makeStruct_(S_prmset, ccScore_prmset_rec, vcFile_prmset, ...
    csDir_rec, cVal_prm, csName_prm, t_fun, vcDir_rec, vcFile_out);

[csDesc, S_best_score] = optimize_param_show_(S_prmset_rec);
assignWorkspace_(S_prmset_rec, S_best_score);
cellstr2file_(strrep(vcFile_out, '.mat', '.txt'), csDesc, 1);
edit(strrep(vcFile_out, '.mat', '.txt'));
end %func


%--------------------------------------------------------------------------
% load from existing cache
function S_score1 = load_score_prmset_(...
    vcSorter, vcDir_in, csName_prm, cVal_prm, P_prmset, iPrmset)

vcDir_out = fullfile(vcDir_in, vcSorter);
fUse_cache = get_set_(P_prmset, 'fUse_cache', 1);
vcFile_true_mda = fullfile(vcDir_in, 'firings_true.mda');
try    
    vcFile_out1 = fullfile(vcDir_out, sprintf('firings_p%d.mda', iPrmset));
    vcFile_score1 = strrep(vcFile_out1, '.mda', '_score.mat');
    if exist_file_(vcFile_score1) && fUse_cache   
        S_score1 = load(vcFile_score1);
        fprintf('Loaded from cache: %s\n', vcFile_score1);
    elseif exist_file_(vcFile_out1) && fUse_cache
        S_score1 = compare_mda_(vcFile_true_mda, vcFile_out1);
        struct_save_(S_score1, vcFile_score1, 1);        
    else
        S_score1 = [];
    end
catch ME
    S_score1 = [];
    fprintf(2, '%s: paramset#%d failed:\n\t%s\n', vcDir_in, iPrmset, ME.message());
end
end %func


%--------------------------------------------------------------------------
function S_score1 = score_prmset_(...
    vcSorter, vcDir_in, csName_prm, cVal_prm, P_prmset, iPrmset)

t_fun = tic;
vcDir_out = fullfile(vcDir_in, vcSorter);
fUse_cache = get_set_(P_prmset, 'fUse_cache', 1);
fParfor = get_set_(P_prmset, 'fParfor', 0);
vcFile_true_mda = fullfile(vcDir_in, 'firings_true.mda');
vcDir_tmp = [];
pause(rand()); % create phase delay
S_score1 = [];
try    
    vcFile_out1 = fullfile(vcDir_out, sprintf('firings_p%d.mda', iPrmset));
    vcFile_score1 = strrep(vcFile_out1, '.mda', '_score.mat');
    if exist_file_(vcFile_score1) && fUse_cache   
        S_score1 = load(vcFile_score1);
        fprintf('Loaded from cache: %s\n', vcFile_score1);
        return;
    end    
    if exist_file_(vcFile_out1) && fUse_cache
        fprintf('Loaded from cache: %s\n', vcFile_out1);
    else
        S_prm = cell2struct(cVal_prm, csName_prm, 1);
        switch lower(vcSorter)
            case 'irc2' % reuse last shared parameter output
                vcFile_prm = fullfile(vcDir_out, sprintf('raw_geom_p%d.prm', iPrmset));
                makeParam_(vcDir_in, vcFile_prm, S_prm, fParfor);
                vcFile_firings = irc2('spikesort-cache', vcFile_prm); 
                vcFile_score = strrep(vcFile_firings, '.mda', '_score.mat');
                if exist_file_(vcFile_score)
                    S_score1 = load(vcFile_score);
                else
                    S_score1 = compare_mda_(vcFile_true_mda, vcFile_firings);
                    struct_save_(S_score1, vcFile_score);
                end
            case 'kilosort2'
                vcDir_tmp = fullfile(vcDir_out, sprintf('p%d', iPrmset));
                kilosort2_(vcDir_in, vcDir_tmp, S_prm);   
                vcFile_firings = fullfile(vcDir_tmp, 'firings.mda');
            case {'mountainsort4', 'spykingcircus', 'tridesclous', 'herdingspikes2', 'klusta', 'waveclus', 'jrclust', 'kilosort'}
                vcDir_tmp = fullfile(vcDir_out, sprintf('p%d', iPrmset));
                vcFile_firings = fullfile(vcDir_tmp, 'firings.mda');
                run_spikeforest2_(vcSorter, vcDir_in, vcDir_tmp, S_prm, 0);
            otherwise
                error('unsupported sorters: %s', vcSorter);
        end
        copyfile(vcFile_firings, vcFile_out1, 'f');
        delete_(vcDir_tmp);
        fprintf('Wrote to %s (took %0.1fs)\n', vcFile_out1, toc(t_fun));
    end
    if isempty(S_score1)
        S_score1 = compare_mda_(vcFile_true_mda, vcFile_out1);
    end
    struct_save_(S_score1, vcFile_score1, 1);
catch ME
    S_score1 = [];
    fprintf(2, '%s: paramset#%d failed:\n\t%s\n', vcDir_in, iPrmset, ME.message());
end
end %func


%--------------------------------------------------------------------------
function vcSorter = infer_sorter_(vcName)
vcName = lower(vcName);
if contains(vcName, {'ironclust','irc2'}), vcSorter = 'irc2';
elseif contains(vcName, {'jrclust'}), vcSorter = 'jrclust';
elseif contains(vcName, {'kilosort2', 'ksort2'}), vcSorter = 'kilosort2';
elseif contains(vcName, {'kilosort', 'ksort'}), vcSorter = 'kilosort';
elseif contains(vcName, {'mountainsort4', 'mountainsort', 'ms4'}), vcSorter = 'mountainsort4';
elseif contains(vcName, {'spykingcircus'}), vcSorter = 'spykingcircus';    
elseif contains(vcName, {'herdingspikes2'}), vcSorter = 'herdingspikes2';
elseif contains(vcName, {'klusta', 'klustakwik'}), vcSorter = 'klusta';
elseif contains(vcName, {'tridesclous'}), vcSorter = 'tridesclous';
elseif contains(vcName, {'waveclus'}), vcSorter = 'waveclus';
else, error('infer_sorter_: unsupported sorter %s', vcName);
end %if
end %func


%--------------------------------------------------------------------------
function mr = mr_set_col_(mr, vi, val)
mr(:,vi) = val;
end %func


%--------------------------------------------------------------------------
function mr = cell_struct_fun_(cS, fh, vcName)
% cc: cell, fh: function handle, vcName: name of the field
mr = nan(size(cS));
for ic=1:numel(cS)
    S1 = cS{ic};
    if isfield(S1, vcName)
        mr(ic) = fh(S1.(vcName));
    end
end
end %func


%--------------------------------------------------------------------------
function cvr_out = cell_struct_join_(cc, vcName, iDimm)
% cc: cell, fh: function handle, vcName: name of the field
cc_out = cell(size(cc));
for ic=1:numel(cc)
    S1 = cc{ic};
    if isfield(S1, vcName)
        cc_out{ic} = S1.(vcName);
    end
end
if iDimm==2, cc_out=cc_out'; end
cvr_out = cell(1, size(cc_out,2));
for iCol=1:size(cc_out,2)
    cc1 = cc_out(:,iCol);
    cvr_out{iCol} = cvr2vr_(cc1);
end
if iDimm==2, cvr_out = cvr_out'; end
end %func


%--------------------------------------------------------------------------
% join cell, return column vector
function vr_column = cvr2vr_(cvr1)
vr_column=[]; 
if isempty(cvr1), return; end
vnRow1 = cellfun(@(x)size(x,1), cvr1);
vnCol1 = cellfun(@(x)size(x,2), cvr1);
if max(vnRow1)==1
    vr_column = cat(2, cvr1{:})';
elseif max(vnCol1)==1
    vr_column = cat(1, cvr1{:});
end
end %func


%--------------------------------------------------------------------------
function plot_prmset_analysis_(vr_prmset, csName_prm, cVal_prm)
% nVal_prm = cellfun(@numel,cVal_prm)';
% viPrm1 = find(nVal_prm>1);
% csName_prm1 = csName_prm(viPrm1);
% nPrm1 = numel(viPrm1);
% [nVal_prm1, csName_prm1] = deal(nVal_prm(viPrm1), csName_prm(viPrm1));

vnType_prm = cellfun(@numel, cVal_prm);
[nPrmset, nPrm] = deal(numel(vr_prmset), numel(cVal_prm));
miPrm_prmset = zeros(nPrm, nPrmset);
for iPrmset = 1:nPrmset
    [~, miPrm_prmset(:,iPrmset)] = permute_prm_(cVal_prm, iPrmset);
end

cScores_prm = cell(nPrm,1);
for iPrm = 1:nPrm
    nType1 = vnType_prm(iPrm);
    if nType1>1
        vr1 = zeros(nType1, 1);
        for iType1 = 1:vnType_prm(iPrm)
            vr_ = vr_prmset(miPrm_prmset(iPrm, :) == iType1);
            vr1(iType1) = nanmean(vr_);
        end
        cScores_prm{iPrm} = vr1;
    end
end

for iPrm = 1:nPrm
    vr1 = cScores_prm{iPrm};
    if isempty(vr1), continue; end
    disp(csName_prm{iPrm});
    cVal_prm1 = cVal_prm{iPrm};
    for iType1 = 1:vnType_prm(iPrm)
        fprintf('\t%s: %f\n', num2str(cVal_prm1{iType1}), vr1(iType1));
    end
end
end %func    
    
    
%--------------------------------------------------------------------------
function [csDesc, S_best_score] = optimize_param_show_(S_prmset)

S_cfg = read_cfg_();
[MAX_PRMSET, SNR_THRESH, vcSnr_mode] = deal(inf, S_cfg.snr_thresh_gt, 'vrSnr_min_clu'); 
SNR_THRESH = 5;
THRESH_SCORE = 80;

% vrVpp_clu, vrVmin_clu, vrSnr_min_clu
csVar_imported = import_struct_(S_prmset);
[~,vcSorter] = fileparts(vcFile_out); vcSorter = strrep(vcSorter, 'scores_prmset_','');
switch vcSorter
    case {'ironclust', 'irc2'}, cColor='r';
    case 'mountainsort4', cColor='b';
    case 'kilosort2', cColor='g';
    case 'klusta', cColor='k';
end

% get unit amplitudes
try
    % paired ground truth
    cS_gt = cellfun_(@(x)load_(fullfile(x, 'raw_geom_gt1.mat'), {vcSnr_mode}), csDir_rec);
    vrSnr_gt = cellfun(@(S)S.(vcSnr_mode), cS_gt);
    [vrSnr_gt_srt, viSnr_gt_srt] = sort(vrSnr_gt);
catch  % SNR not saved
    csScore = {'vrAccuracy_gt', 'vrF1_gt', 'vrPrecision_gt', 'vrRecall_gt'};
    switch lower(S_cfg.vcMode_optimize)
        case 'mean'
            mr_func_ = @(x)cell_struct_fun_(ccScore_prmset_rec, @(y)nanmean(y), x);
            vcScore = 'mean';
        case 'median'
            mr_func_ = @(x)cell_struct_fun_(ccScore_prmset_rec, @(y)nanmedian(y), x);
            vcScore = 'median';            
        case 'count'
            mr_func_ = @(x)cell_struct_fun_(ccScore_prmset_rec, @(y)sum(y>=THRESH_SCORE), x);
            vcScore = sprintf('count|x>%0.1f', THRESH_SCORE);
        case {'mean_pooled'}
            mr_func_ = @(x)cellfun(@(y)nanmean(y), cell_struct_join_(ccScore_prmset_rec, x, 1));
            vcScore = 'mean_pooled';
        case {'mean_pooled_clu'}
            mr_func_ = @(x)cellfun(@(y)nanmean(y), cell_struct_join_(ccScore_prmset_rec, strrep(x,'_gt','_clu'), 1));
            vcScore = 'mean_pooled';            
        case {'median_pooled'}
            mr_func_ = @(x)cellfun(@(y)nanmedian(y), cell_struct_join_(ccScore_prmset_rec, x, 1));
            vcScore = 'median_pooled';            
        case {'count_pooled'}
            mr_func_ = @(x)cellfun(@(y)sum(y>=THRESH_SCORE), cell_struct_join_(ccScore_prmset_rec, x, 1));
            vcScore = 'count_pooled';            
        case {'count_pooled_clu'}
            mr_func_ = @(x)cellfun(@(y)sum(y>=THRESH_SCORE), cell_struct_join_(ccScore_prmset_rec, strrep(x,'_gt','_clu'), 1));
            vcScore = 'count_pooled_clu';                        
        otherwise, error('optimize_param_show_: unsupported `vcMode_optimize`');
    end
    cmrScore_prmset_gt = cellfun_(@(x)mr_func_(x), csScore);
    [vrAccuracy_prmset, vrF1_prmset, vrPrecision_prmset, vrRecall_prmset] = ...
        multifun_(@(x)nanmean(cmrScore_prmset_gt{x},1), 1, 2, 3, 4);    
    plot_prmset_analysis_(vrAccuracy_prmset, csName_prm, cVal_prm);
    vrAccuracy_prmset(isnan(vrAccuracy_prmset)) = 0;
    [~, viPrmset_srt] = sort(vrAccuracy_prmset, 'descend');
    csDesc = {};
    for iiPrmset = 1:numel(viPrmset_srt)
        iPrmset = viPrmset_srt(iiPrmset);
        cVal_prm1 = permute_prm_(cVal_prm, iPrmset);
        csDesc{end+1} = '----------';
        csDesc{end+1} = sprintf('Prmset-rank %d (p#%d):', iiPrmset, iPrmset);        
        csDesc{end+1} = sprintf('  (%s) Accuracy:%0.1f, F1:%0.1f, Precision:%0.1f, Recall:%0.1f', ...
            vcScore, vrAccuracy_prmset(iPrmset), vrF1_prmset(iPrmset), vrPrecision_prmset(iPrmset), vrRecall_prmset(iPrmset));
        for iPrm = 1:numel(cVal_prm)
            csDesc{end+1} = sprintf('  %s: %s', csName_prm{iPrm}, numstr_(cVal_prm1{iPrm}));
        end
    end
%     disp_cs_(csDesc);
    S_best_score = [];
    return;
end

vlPlot_gt = vrSnr_gt>=SNR_THRESH;
get_score_ = @(vc)cellfun(@(x)nanmean(get_(x, vc)), ccScore_prmset_rec);
mean_mr_snr_ = @(mr)nanmean(mr_set_col_(mr,~vlPlot_gt,nan),2);

% display the optimized param
[nPrm, nRec, nPrmset] = deal(numel(cVal_prm), size(ccScore_prmset_rec,2), size(ccScore_prmset_rec,1));
mrF1_prmset_rec = get_score_('vrF1_gt'); vrF1_mean_prmset = mean_mr_snr_(mrF1_prmset_rec); 
[~, iPrmset_best] = max(vrF1_mean_prmset);
[~, iPrmset_worst] = min(vrF1_mean_prmset);
mrAccuracy_prmset_rec = get_score_('vrAccuracy_gt'); vrAccuracy_mean_prmset = mean_mr_snr_(mrAccuracy_prmset_rec);
mrPrecision_prmset_rec = get_score_('vrPrecision_F1_gt'); vrPrecision_mean_prmset = mean_mr_snr_(mrPrecision_prmset_rec);
mrRecall_prmset_rec = get_score_('vrRecall_F1_gt'); vrRecall_mean_prmset = mean_mr_snr_(mrRecall_prmset_rec);


% do the SNR analysis by recording?
[vrF1_prmset_srt, viPrmset_srt] = sort(vrF1_mean_prmset, 'descend');
cell_ = cellfun_(@(x)x(viPrmset_srt), {vrAccuracy_mean_prmset, vrPrecision_mean_prmset, vrRecall_mean_prmset});
[vrAccuracy_srt, vrPrecision_srt, vrRecall_srt] = deal(cell_{:});
cmrScore_prmset_rec = {mrF1_prmset_rec, mrPrecision_prmset_rec, mrRecall_prmset_rec, mrAccuracy_prmset_rec};
get_scores_prmset_ = @(i)cellfun_(@(x)x(i,:), cmrScore_prmset_rec);

csName_score = {'F1', 'Precision', 'Recall', 'Accuracy'};
cVal_best_score = get_scores_prmset_(iPrmset_best);
cVal_worst_score = get_scores_prmset_(iPrmset_worst);
S_best_score = cell2struct(cVal_best_score, csName_score, 2);
S_worst_score = cell2struct(cVal_worst_score, csName_score, 2);
S_best_score = struct_add_(S_best_score, vrSnr_gt, vcSnr_mode);
S_worst_score = struct_add_(S_worst_score, vrSnr_gt, vcSnr_mode);

nPrmset_show = min(MAX_PRMSET, nPrmset);
viX = 1:nPrmset_show;

fig_ = @()figure('Color','w','Name', vcSorter, 'NumberTitle', 'off');
csLegend = cellfun_(@(x)sprintf('%s.%s', vcSorter, x), {'best','worst', 'diff'});
ax1=axes(fig_()); hold(ax1,'on'); 
% vcShapes = 'ox+*';
calc_mean_snr_ = @(vr)arrayfun(@(x)nanmean(vr(vrSnr_gt>=x)), vrSnr_gt_srt);
yyaxis(ax1, 'left');
[vrY1_L1, vrY1_L2] = deal(calc_mean_snr_(S_best_score.F1), calc_mean_snr_(S_worst_score.F1));
plot(ax1, vrSnr_gt_srt, vrY1_L1, [cColor, '-'], vrSnr_gt_srt, vrY1_L2, [cColor, '--']); 
xylabel_(ax1, 'SNR threshold (=x)', '<F1-score> | SNR>=x'); ylim(ax1, [0 100]);
yyaxis(ax1, 'right'); 
plot(ax1, vrSnr_gt_srt, abs(vrY1_L1 - vrY1_L2), [cColor,':']); ylabel(ax1, '|F1.best - F1.worst|');
grid(ax1,'on'); ylim(ax1, [0 100/4]);
legend(ax1, csLegend, 'Location', 'SE');

ax2=axes(fig_()); hold(ax2,'on');
excl_nan_ = @(x)x(~isnan(x));
[vrY2, vrY3] = deal(excl_nan_(S_best_score.F1), excl_nan_(S_worst_score.F1));
[vrX2_L1, vrX2_L2, vrY2_L1, vrY2_L2] = deal(sort(vrY2), sort(vrY3), numel(vrY2):-1:1, numel(vrY3):-1:1);
vrX2_R = linspace(max(vrX2_L1(1),vrX2_L2(1)), min(vrX2_L1(end),vrX2_L2(end)), 100);
yyaxis(ax2, 'left');
plot(ax2, vrX2_L1, vrY2_L1, [cColor,'-'], vrX2_L2, vrY2_L2, [cColor,'--']);
xylabel_(ax2, 'F1 threshold (=x)', '# GT-units | F1>=x'); ylim(ax2, [0 nRec]);   
yyaxis(ax2, 'right'); 
vrY2_R = abs(interp1(vrX2_L1, vrY2_L1, vrX2_R, 'linear') - interp1(vrX2_L2, vrY2_L2, vrX2_R, 'linear'));
plot(ax2, vrX2_R, vrY2_R, [cColor,':']); ylabel(ax2, '|count.best - count.worst|');
grid(ax2,'on'); ylim(ax2, [0 nRec/4]);
legend(ax2, csLegend, 'Location', 'NE');

% create a bar plot
ax3 = axes(fig_());
stairs(ax3, linspace(0,1,numel(viX)), vrF1_prmset_srt(viX), cColor);
% bar(ax3, [vrF1_prmset_srt(viX), vrPrecision_srt(viX), vrRecall_srt(viX), vrAccuracy_srt(viX)]); 
xylabel_(ax3, 'prmset # (sorted by F1-score)', sprintf('F1-score|SNR>=%0.1f',SNR_THRESH), ...
    disp_stats_(S_best_score.F1(vlPlot_gt), sprintf('%s|SNR>=%0.1f', csName_score{1}, SNR_THRESH), 0), 1);
axis(ax3, [0 1 0 100]);
% legend(ax3, {'F1','Precision','Recall','Accuracy'}, 'Location', 'SE');

fig_();
nScores = numel(csName_score);
for iScore=1:nScores
    subplot(nScores, 1, iScore);
    vrScore_ = S_best_score.(csName_score{iScore});
    plot(vrSnr_gt, vrScore_, [cColor,'.']);
    xylabel_(gca, vcSnr_mode, csName_score{iScore}, disp_stats_(vrScore_(vlPlot_gt), ...
        sprintf('%s|SNR>=%0.1f', csName_score{iScore}, SNR_THRESH), 0), 1);
    hold on; plot(SNR_THRESH*[1,1], [0 100], 'k-');
end

% output summary text
csDesc = {};
csDesc{end+1} = sprintf('');    
csDesc{end+1} = sprintf('------------------------------');    
csDesc{end+1} = sprintf('  vcDir_rec:              %s', vcDir_rec);
csDesc{end+1} = sprintf('  vcFile_prmset:          %s', vcFile_prmset);
csDesc{end+1} = sprintf('  # Recordings:           %d', nRec);
csDesc{end+1} = sprintf('  # Recordings|SNR>=%0.1f:  %d', SNR_THRESH, sum(vlPlot_gt));
csDesc{end+1} = sprintf('  # parameters:           %d', nPrm);
csDesc{end+1} = sprintf('  # parameter sets:       %d', nPrmset);
csDesc{end+1} = sprintf('  run-time (s):           %0.1f', t_fun);
csDesc{end+1} = sprintf('------------------------------');
csDesc{end+1} = sprintf('  Best scores:');
csDesc{end+1} = disp_stats_();
for iScore = 1:numel(csName_score)
    csDesc{end+1} = disp_stats_(...
        cVal_best_score{iScore}(vlPlot_gt), csName_score{iScore});
end
csDesc{end+1} = sprintf('  Worst scores:');
csDesc{end+1} = disp_stats_();
for iScore = 1:numel(csName_score)
    csDesc{end+1} = disp_stats_(...
        cVal_worst_score{iScore}(vlPlot_gt), csName_score{iScore});
end

% show prmset
for iiPrmset = 1:nPrmset_show
    iPrmset1 = viPrmset_srt(iiPrmset);
    csDesc{end+1} = sprintf('------------------------------');
    csDesc{end+1} = sprintf('  prmset rank #%d (p#%d):', iiPrmset, iPrmset1);
    csDesc{end+1} = sprintf('    F1:%0.1f Precision:%0.1f Reccall:%0.1f Accuracy:%0.1f', ...
        cellfun(@(x)nanmean(x(vlPlot_gt)), get_scores_prmset_(iPrmset1)));
    cVal_prm1 = permute_prm_(cVal_prm, iPrmset1);
    for iPrm = 1:nPrm
        csDesc{end+1} = sprintf('  %s: %s', cName_prm{iPrm}, numstr_(cVal_prm1{iPrm}));
    end
end

csDesc{end+1} = sprintf('');    
% disp_cs_(csDesc);

% fprintf(2, '@TODO: show # units above accuracy 80\n');
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function cellstr2file_(vcFile, csLines, fVerbose)
% Write a cellstring to a text file
if nargin<3, fVerbose = 0; end
vcDir = fileparts(vcFile);
mkdir_(vcDir);
fid = fopen(vcFile, 'w');
for i=1:numel(csLines)
    fprintf(fid, '%s\n', csLines{i});
end
fclose(fid);
if fVerbose
    fprintf('Wrote to %s\n', vcFile);
end
end %func


%--------------------------------------------------------------------------
function vc = numstr_(vr)
if ischar(vr)
    vc = vr;
elseif isnumeric(vr)
    vc = num2str(vr);
    if numel(vr)>1, vc = sprintf('[%s]', vc); end
else
    vc = vr;
end
end %func


%--------------------------------------------------------------------------
function P = file2struct_ordered_(vcFile_file2struct)
% James Jun 2017 May 23
% Run a text file as .m script and result saved to a struct P
% _prm and _prb can now be called .prm and .prb files
P = []; 
if ~exist_file_(vcFile_file2struct), return; end
    
% load text file. trim and line break. remove comments.  replace 
csLines = file2lines_(vcFile_file2struct);
P = cs2struct_(csLines);
end %func


%--------------------------------------------------------------------------
function P = cs2struct_(csLines)
% cell stringt o struct
% csLines = {'name1=value1', 'name2=value2', ...}

P = [];    
csLines = strip_comments_(csLines);
if isempty(csLines), return; end

P = struct();
for iLine = 1:numel(csLines)
    try
        vcLine1 = strtrim(csLines{iLine});
        if ~isempty(vcLine1)
            eval(sprintf('P.%s;', vcLine1));
        end
    catch
        fprintf(2, '%s\n', lasterr());
    end
end
end %func


%--------------------------------------------------------------------------
% Strip comments from cell string
% 7/24/17 JJJ: Code cleanup
function csLines = strip_comments_(csLines)
csLines = csLines(cellfun(@(x)~isempty(x), csLines));
csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
csLines = csLines(cellfun(@(x)x(1)~='%', csLines));

% remove comments in the middle
for i=1:numel(csLines)
    vcLine1 = csLines{i};
    iComment = find(vcLine1=='%', 1, 'first');
    if ~isempty(iComment)
        vcLine1 = vcLine1(1:iComment-1);
    end
    vcLine1 = strrep(vcLine1, '...', '');
    if ismember(strsplit(vcLine1), {'for', 'end', 'if'})
        csLines{i} = [strtrim(vcLine1), ', ']; %add blank at the end
    else
        csLines{i} = [strtrim(vcLine1), ' ']; %add blank at the end
    end
end
% csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
csLines = csLines(cellfun(@(x)~isempty(x), csLines));
end %func


%--------------------------------------------------------------------------
% Read a text file and output cell strings separated by new lines
% 7/24/17 JJJ: Code cleanup
function csLines = file2lines_(vcFile_file2struct)
csLines = {};
if ~exist_file_(vcFile_file2struct, 1), return; end

fid = fopen(vcFile_file2struct, 'r');
csLines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

csLines = csLines{1};
end %func


%--------------------------------------------------------------------------
function kilosort2_(vcDir_in, vcDir_out, S_prm)

fprintf('Running kilosort2: %s\n', vcDir_in); t_fun = tic;

% add path
S_cfg = read_cfg_();
kilosort_src = S_cfg.path_ksort2;
ironclust_src = fileparts(mfilename('fullpath'));
addpath(genpath(ironclust_src));
addpath(genpath(kilosort_src));

% convert input path
vcFile_raw = fullfile(vcDir_in, 'raw.mda');
vcFile_firings = fullfile(vcDir_out, 'firings.mda');
ops = import_ksort2_(vcFile_raw, vcFile_firings, S_prm);

%-----
% call kilosort2 functions

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
% save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];
close all; % close all open figures
fprintf('\n\tfound %d good units \n', sum(rez.good>0))
%-----

% Export kilosort
export_ksort2_(rez, vcFile_firings);
fprintf('kilosort2 finished (took %0.1fs): %s\n', toc(t_fun), vcFile_firings);

end %func


%--------------------------------------------------------------------------
function ops = import_ksort2_(vcFile_raw, vcFile_firings, S_prm)

vcDir_in = fileparts(vcFile_raw);
vcDir_out = fileparts(vcFile_firings);
mkdir_(vcDir_out); % create dir if doesn't exist

% convert .bin file
S_json = loadjson_(fullfile(vcDir_in, 'params.json'));
[sRateHz, detect_sign] = get_(S_json, 'samplerate', 'spike_sign');
vcFile_bin = strrep(vcFile_raw, '.mda', '.bin');
[nChans, ~] = mda2bin_(vcFile_raw, vcFile_bin, detect_sign);

% create a probe file
mrXY_site = csvread(fullfile(vcDir_in, 'geom.csv'));
vcFile_chanMap = fullfile(vcDir_out, 'chanMap.mat');
createChannelMapFile_(vcFile_chanMap, nChans, mrXY_site(:,1), mrXY_site(:,2));

% copy from S_prm struct
detect_threshold = get_set_(S_prm, 'detect_threshold', 6);
freq_min = get_set_(S_prm, 'freq_min', 150);
nPCs = get_set_(S_prm, 'nPCs', 3);
sigmaMask = get_set_(S_prm, 'sigmaMask', 30);
minFR = get_set_(S_prm, 'minFR', 1/50);
preclust_threshold = get_set_(S_prm, 'preclust_threshold', 8);
projection_threshold = get_set_(S_prm, 'projection_threshold', [10, 4]);
car = get_set_(S_prm, 'car', 1);
Nt = get_set_(S_prm, 'Nt', 128 * 1024 * 5 + 64);

%-----
% assign kilosort2 parameters
ops.fpath = vcDir_out;
ops.fproc = fullfile(vcDir_out, 'temp_wh.dat'); % proc file on a fast SSD  ;
ops.trange = [0 Inf]; % time range to sort
ops.NchanTOT = nChans; % total number of channels in your recording
ops.fbinary = vcFile_bin;   % the binary file is in this folder
ops.chanMap = vcFile_chanMap;
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if no chanMap file
ops.fs = sRateHz;           % sample rate
ops.fshigh = freq_min;      % frequency for high pass filtering (150)
ops.minfr_goodchannels = 0.1; % minimum firing rate on a "good" channel (0 to skip)d
ops.Th = projection_threshold;           % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.lam = 10;               % how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.AUCsplit = 0.9;         % splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.minFR = minFR;           % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.momentum = [20 400];    % number of samples to average over (annealed from first to second value) 
ops.sigmaMask = sigmaMask;         % spatial constant in um for computing residual variance of spike
ops.ThPre = preclust_threshold;              % threshold crossings for pre-clustering (in PCA projection space)

% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh = -abs(detect_threshold);      % spike threshold in standard deviations (-6)
ops.reorder = 1;            % whether to reorder batches for drift correction. 
ops.nskip = 25;             % how many batches to skip for determining spike PCs
ops.GPU = 1;                % has to be 1, no CPU version yet, sorry
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor = 4;       % max number of clusters per good channel (even temporary ones)
ops.ntbuff = 64;            % samples of symmetrical buffer for whitening and spike detection
ops.NT = Nt;                % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
ops.whiteningRange = 32;    % number of channels to use for whitening each channel
ops.nSkipCov = 25;          % compute whitening matrix from every N-th batch
ops.scaleproc = 200;        % int16 scaling of whitened data
ops.nPCs = nPCs;            % how many PCs to project the spikes into
ops.useRAM = 0;             % not yet available
ops.CAR = car;              % perform CAR
end %func


%--------------------------------------------------------------------------
function S_chanMap = createChannelMapFile_(vcFile_channelMap, Nchannels, xcoords, ycoords, shankInd)
if nargin<6, shankInd = []; end

connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = xcoords(:);
ycoords   = ycoords(:);

if isempty(shankInd)
    shankInd   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
end
[~, name, ~] = fileparts(vcFile_channelMap);
S_chanMap = makeStruct_(chanMap, connected, xcoords, ycoords, shankInd, chanMap0ind, name);
save(vcFile_channelMap, '-struct', 'S_chanMap')
end %func


%--------------------------------------------------------------------------
% convert mda to int16 binary format, flip polarity if detect sign is
% positive
function [nChans, nSamples] = mda2bin_(vcFile_raw, vcFile_bin, detect_sign)
if nargin<3, detect_sign=-1; end

% if exist_file_(vcFile_bin)
%     S_mda = readmda_header_(vcFile_raw);
%     [nChans, nSamples] = deal(S_mda.dimm(1), S_mda.dimm(2));
%     return;
% else
    mr = readmda_(vcFile_raw);
    % adjust scale to fit int16 range with a margin
    if isa(mr,'single') || isa(mr,'double')
        uV_per_bit = 2^14 / max(abs(mr(:)));
        mr = int16(mr * uV_per_bit);
    end
    [nChans, nSamples] = size(mr);
    if detect_sign > 0, mr = -mr; end % force negative detection
    write_bin_(vcFile_bin, mr, 1);
% end
end %func


%--------------------------------------------------------------------------
function export_ksort2_(rez, firings_out_fname)

mr_out = zeros(size(rez.st3,1), 3, 'double'); 
mr_out(:,2) = rez.st3(:,1); %time
mr_out(:,3) = rez.st3(:,2); %cluster
writemda_(firings_out_fname, mr_out');
end %func


%--------------------------------------------------------------------------
function [out1, viPrm1] = permute_prm_(cVal_prm, iPrmset)
% usage
% -----
% permute_prm_(cVal_prm)
% permute_prm_

if nargin<2, iPrmset=[]; end
viPrm1 = [];
dimm_prmset = cellfun(@numel, cVal_prm);
if isempty(iPrmset)
    out1 = prod(dimm_prmset);
else    
    viPrm1 = flipud(ind2sub_(flipud(dimm_prmset), iPrmset));
    out1 = cellfun_(@(x,y)x{y}, cVal_prm, arrayfun_(@(x)x, viPrm1));
end
end %func


%--------------------------------------------------------------------------
function vi_out = ind2sub_(siz,ndx)

vi_out = zeros(size(siz));
k = cumprod(siz);
for i = numel(siz):-1:2
    vi = rem(ndx-1, k(i-1)) + 1;
    vj = (ndx - vi)/k(i-1) + 1;
    vi_out(i) = double(vj);
    ndx = vi;
end
vi_out(1) = ndx;
end% func


%--------------------------------------------------------------------------
% clear slurm jobs
function clear_jobs_(vcArg1)

if nargin<1, vcArg1=''; end
if isempty(vcArg1), vcArg1 = 'local'; end
myCluster = parcluster(vcArg1);
Jobs = myCluster.Jobs;
nJobs = numel(Jobs);
if nJobs>0
    delete(Jobs); 
    fprintf('Cleared %d jobs\n', nJobs);
else
    fprintf('No jobs to clear\n');
end
end


%--------------------------------------------------------------------------
function [mrX_clu_drift, mrY_clu_drift] = position_clu_drift_(S0)

S_auto = S0.S_auto;
S_drift = S0.S_clu.S_drift;
min_spk = S0.P.knn/2;

switch 2
    case 1, fh = @(x)trimmean(x,20);
    case 2, fh = @(x)median(x);
end

[viLim_drift, nDrift] = get_(S_drift, 'viLim_drift', 'nTime_drift');
[mrX_clu_drift, mrY_clu_drift] = deal(nan(S_auto.nClu, nDrift));
get_lim_ = @(x,i)find(x>=viLim_drift(i) & x<viLim_drift(i+1));
for iClu=1:S_auto.nClu
    viSpk1 = S_auto.cviSpk_clu{iClu};
    [vrX1,vrY1] = deal(S0.mrPos_spk(viSpk1,1), S0.mrPos_spk(viSpk1,2));
    for iDrift = 1:nDrift
        vi_ = get_lim_(viSpk1,iDrift);
        if numel(vi_) >= min_spk
            mrX_clu_drift(iClu, iDrift) = fh(vrX1(vi_));
            mrY_clu_drift(iClu, iDrift) = fh(vrY1(vi_));
        end
    end
end
end %func


%--------------------------------------------------------------------------
function S_pos_clu = static_position_clu_(S0, fPlot)
if nargin<2, fPlot=[]; end
[mrX_clu_drift, mrY_clu_drift] = position_clu_drift_(S0);

% select most smooth transitions
mrDY_clu_drift = mrY_clu_drift(:,2:end)-mrY_clu_drift(:,1:end-1);
vrDY_drift = trimmean(mrDY_clu_drift, 20);
vrGain_clu = nanmedian(abs(mrDY_clu_drift ./ vrDY_drift),2);
vrY_drift = cumsum([0, vrDY_drift]);
[nClu, nDrift] = size(mrY_clu_drift);
mrY0_clu_drift = mrY_clu_drift;
mrY_clu_drift = mrY_clu_drift-vrY_drift.*vrGain_clu;
mrPos_clu = [nanmedian(mrX_clu_drift,2), nanmedian(mrY_clu_drift,2)];
S_pos_clu = makeStruct_(mrPos_clu, mrX_clu_drift, mrY_clu_drift, mrY0_clu_drift, vrY_drift, vrGain_clu);
if fPlot
    viClu_plot = find(~any(isnan(mrY0_clu_drift),2));
    viClu_plot = viClu_plot(:)';
    plot_ = @(x,y,c)plot(x,y,'.','color',c,'MarkerSize',5);
    figure('Color','w'); ax=[]; 
    rand('seed',0); mrColor_clu = rand(nClu,3);
    for iMode=1:3
        ax(end+1) = subplot(1,3,iMode);
        hold on; grid on;
        for iClu=viClu_plot
            vrColor1 = mrColor_clu(iClu,:);
            switch iMode
                case 3
                    plot_(mrX_clu_drift(iClu,:), mrY_clu_drift(iClu,:), vrColor1); 
                    plot(mrPos_clu(iClu,1), mrPos_clu(iClu,2), 'k*');
%                     text(mrPos_clu(iClu,1), mrPos_clu(iClu,2), sprintf('%d', iClu));
                    title('xy (corrected)');                    
                case 1, plot_(1:nDrift, mrY0_clu_drift(iClu,:), vrColor1); title('original');
                case 2, plot_(1:nDrift, mrY_clu_drift(iClu,:), vrColor1); title('corrected');
            end
        end
    end
    linkaxes(ax,'y');
end
end %func


%--------------------------------------------------------------------------
function plot0_(vcMode, vcFile_prm, vcArg1)

S0 = load0_(vcFile_prm);
P = S0.P;
S_auto = S0.S_auto;
S_drift = S0.S_clu.S_drift;

% [mrX_clu_drift, mrY_clu_drift] = position_clu_drift_(S0);

% viClu_use = find(~isnan(vrJitter_clu));
% figure; plot(mrY_clu_drift(viClu_use,:)');

nClu = max(S_auto.viClu);
cviSpk_clu = arrayfun_(@(x)find(S_auto.viClu==x), (1:nClu)');
mrXY_clu = cell2mat_(cellfun_(@(x)median(S0.mrPos_spk(x,:)), cviSpk_clu));
if isempty(vcArg1)
    viClu_plot = [];
else
    viClu_plot = str2num(vcArg1);
    vl_spk = ismember(S_auto.viClu, viClu_plot);
end

switch lower(vcMode)
    case 'scatter'
        hFig = create_figure_([], [0 0 .5 1], P.vcFile, 1, 1);
        [vrX1, vrY1] = deal(S0.mrPos_spk(:,1), S0.mrPos_spk(:,2));
        if isempty(viClu_plot)
            scatter(vrX1, vrY1, 1, S_auto.viClu);
            text_(mrXY_clu(:,1), mrXY_clu(:,2), arrayfun_(@num2str, 1:nClu), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        else
            plot(vrX1(~vl_spk), vrY1(~vl_spk), '.', 'MarkerSize', 1, 'Color', ones(1,3)*.75); hold on;
            scatter(vrX1(vl_spk), vrY1(vl_spk), 1, S_auto.viClu(vl_spk));
            text_(mrXY_clu(viClu_plot,1), mrXY_clu(viClu_plot,2), arrayfun_(@num2str, viClu_plot), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end
        xylabel_([], 'x pos (um)', 'y pos (um)', P.vcFile); grid on;
        
    case 'driftcorr', S_pos_clu = static_position_clu_(S0, 1);

    case 'drift'
        S_pos_clu = static_position_clu_(S0, 1);
        % plot high-occupancy clusters only
        viClu_plot = find(sum(isnan(S_pos_clu.mrY0_clu_drift),2) == 0);
        vl_spk = ismember(S_auto.viClu, viClu_plot);
        viRand_clu = randperm(S_auto.nClu);
        
        hFig = create_figure_([], [0 0 .5 1], P.vcFile, 1, 1);
        [vrX1, vrY1] = deal(double(S0.viTime_spk)/P.sRateHz, S0.mrPos_spk(:,2));
        if isempty(viClu_plot)
            vl1 = S_auto.viClu>0;            
            rand_map = randperm(max(S_auto.viClu));
            viClu1 = rand_map(S_auto.viClu(vl1));
            scatter(vrX1(vl1), vrY1(vl1), 1, viClu1);
            text_(mrXY_clu(:,1), mrXY_clu(:,2), arrayfun_(@num2str, 1:nClu), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        else
%             plot(vrX1(~vl_spk), vrY1(~vl_spk), '.', 'MarkerSize', 1, 'Color', ones(1,3)*.75); hold on;
            scatter(vrX1(vl_spk), vrY1(vl_spk), 1, viRand_clu(S_auto.viClu(vl_spk)));
%             text_(mrXY_clu(viClu_plot,1), mrXY_clu(viClu_plot,2), arrayfun_(@num2str, viClu_plot), ...
%                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');            
        end
        xylabel_([], 'Time (s)', 'y pos (um)', P.vcFile); grid on;
        colormap jet;
        
    case 'drift-clip'
        ylim1 = str2num(vcArg1);
        hFig = create_figure_([], [0 0 .5 1], P.vcFile, 1, 1);
        [vrX1, vrY1] = deal(double(S0.viTime_spk)/P.sRateHz, S0.mrPos_spk(:,2));
        vl1 = S_auto.viClu>0 & vrY1>=ylim1(1) & vrY1<=ylim1(2); 
        viClu1 = S_auto.viClu(vl1);
        viClu_uniq = unique(viClu1);
        linMap = 1:max(viClu1); linMap(viClu_uniq) = 1:numel(viClu_uniq);
        rand_map = randperm(max(linMap));
        viClu1 = rand_map(linMap(viClu1));
        scatter(vrX1(vl1), vrY1(vl1), 1, viClu1);
        set(gca, 'XTickLabel',{},'YTickLabel',{},'XTick',[],'YTick',[]);
        colormap jet;
        
    case 'drift-clu'
        hFig = create_figure_([], [0 0 .5 1], P.vcFile, 1, 1);
        find_spk_ = @(x,y)intersect(S_auto.cviSpk_clu{x}, S_drift.cviSpk_drift{y});
        [xx,yy]=meshgrid(1:S_auto.nClu, 1:S_auto.S_drift.nTime_drift);
        vrY1 = arrayfun(@(x,y)median(S0.mrPos_spk(find_spk_(x,y),2)), xx, yy);
        vrX1 = (1:S_auto.S_drift.nTime_drift) * P.step_sec_drift;
        if isempty(viClu_plot)
            plot(vrX1, vrY1, '.-');
            text_(zeros(S_auto.nClu,1), mrXY_clu(:,2), arrayfun_(@num2str, 1:nClu), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        else
            plot(vrX1, vrY1(:,viClu_plot), '.-');
            text_(zeros(numel(viClu_plot),1), mrXY_clu(viClu_plot,2), arrayfun_(@num2str, viClu_plot), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end
        xylabel_([], 'Time (s)', 'y pos (um)', P.vcFile); grid on;        
        
    case 'drift3-clu'
        hFig = create_figure_([], [0 0 .5 1], P.vcFile, 1, 1);
        find_spk_ = @(x,y)intersect(S_auto.cviSpk_clu{x}, S_drift.cviSpk_drift{y});
        [xx,yy]=meshgrid(1:S_auto.nClu, 1:S_auto.S_drift.nTime_drift);
        vrZ1 = arrayfun(@(x,y)median(S0.mrPos_spk(find_spk_(x,y),1)), xx, yy);
        vrY1 = arrayfun(@(x,y)median(S0.mrPos_spk(find_spk_(x,y),2)), xx, yy);
        vrX1 = (1:S_auto.S_drift.nTime_drift) * P.step_sec_drift;
        if isempty(viClu_plot)
            plot3(vrX1, vrY1, vrZ1, '.-');
            text3_(zeros(S_auto.nClu,1), mrXY_clu(:,2), mrXY_clu(:,1), arrayfun_(@num2str, 1:nClu), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        else
            plot3(vrX1, vrY1(:,viClu_plot), vrZ1(:,viClu_plot), '.-');
            text3_(zeros(numel(viClu_plot),1), mrXY_clu(viClu_plot,2), mrXY_clu(viClu_plot,1), arrayfun_(@num2str, viClu_plot), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end
        xylabel_([], 'Time (s)', 'y pos (um)', P.vcFile); grid on;       
        zlabel('x pos (um)');
        
    case 'scatter3'
        hFig = create_figure_([], [0 0 .5 1], P.vcFile, 1, 1);
        vrLogP_clu = cellfun(@(x)median(log(S0.vrPow_spk(x))), cviSpk_clu);
        [vrX1, vrY1, vrZ1] = deal(S0.mrPos_spk(:,1), S0.mrPos_spk(:,2), log(S0.vrPow_spk));
        if isempty(viClu_plot)
            scatter3(vrX1, vrY1, vrZ1, 1, S_auto.viClu);
            text3_(vrX1, vrY1, vrZ1, arrayfun_(@num2str, 1:nClu), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        else
            plot3(vrX1(~vl_spk), vrY1(~vl_spk), vrZ1(~vl_spk), '.', 'MarkerSize', 1, 'Color', ones(1,3)*.75); hold on;
            scatter3(vrX1(vl_spk), vrY1(vl_spk), vrZ1(vl_spk), 1, S_auto.viClu(vl_spk));
            text3_(mrXY_clu(viClu_plot,1), mrXY_clu(viClu_plot,2), vrZ1(viClu_plot), arrayfun_(@num2str, viClu_plot), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');        
        end
        xylabel_([], 'x pos (um)', 'y pos (um)', P.vcFile); grid on;
        zlabel('log power');
    case 'probe'
        hFig = create_figure_([], [0 0 .5 1], P.vcFile, 1, 1);
        plot_probe_(P.mrSiteXY, [12 12], P.viSite2Chan);
        axis equal;
    case {'verify', 'validate'}
        validate_(P);
    otherwise
        error('plot0_: invalid option');
end %switch
end %func


%--------------------------------------------------------------------------
function edit_readme_()
path_code = read_cfg_path_('path_vscode');
system_('%s "%s"', path_code, ircpath_(fullfile('..', 'README.md')));
% S_cfg = read_cfg_();
% if ispc()
%     path_code = ['"', get_(S_cfg, 'path_vscode_win'), '"'];
% elseif ismac()
%     path_code = get_(S_cfg, 'path_vscode_mac');
% elseif isunix()
%     path_code = get_(S_cfg, 'path_vscode_lin');
% else
%     path_code = 'code';
% end
% path_irc = fileparts(fileparts(mfilename('fullpath')));
% path_readme = fullfile(path_irc, 'README.md');
% system(sprintf('%s "%s"', path_code, path_readme));
end %func


%--------------------------------------------------------------------------
function push_readme_()
system('git add ../README.md');
system('git commit -m "README.md updated"');
system('git push');
end %func


%--------------------------------------------------------------------------
function git_push_(vcArg1)
system('git add -u .');
system('git add -u ../');
system('git status');
system(sprintf('git commit -m "%s"', vcArg1));
system('git push');
end %func


%--------------------------------------------------------------------------
function optimize_status_(vcDir_rec, vcFile_prmset, hAx)
if nargin<3, hAx=[]; end

% usage
% optimize_status_(vcFile_list, vcFile_prmset)
% optimize_status_(vcDir_rec, vcFile_prmset)
if nargin<2, vcFile_prmset=''; end
if isempty(vcFile_prmset), vcFile_prmset=read_cfg_('ironclust_prmset'); end
try
    hMsg = [];
    if ~isempty(hAx), hMsg = msgbox('updating...'); end
    vcArg1 = vcDir_rec;
    if exist_file_(vcDir_rec)
        csDir_rec = load_batch_(vcDir_rec);
        vcDir_rec=fileparts(vcDir_rec); 
    elseif exist_dir_(vcDir_rec)
        csDir_rec = {vcDir_rec};
    else
        error('not found: %s', vcDir_rec);
    end
    
    assert(exist_file_(vcFile_prmset) && exist_dir_(vcDir_rec), 'file or dir does not exist');

    S_prmset = file2struct_ordered_(vcFile_prmset);
    [cName_prm, cVal_prm] = deal(fieldnames(S_prmset), struct2cell(S_prmset));
    nPrmset = prod(cellfun(@numel, cVal_prm));        
    vcSorter = lower(strrep(vcFile_prmset, '.prmset', ''));
    vS_dir = cellfun_(@(x)dir(fullfile(x, vcSorter, 'firings_p*.mda')), csDir_rec);
    vS_dir = cat(1, vS_dir{:}); 
    csName_file = arrayfun_(@(x)x.name, vS_dir);
    csPath_file = arrayfun_(@(x)fullfile(x.folder, x.name), vS_dir);    
    viPrmset = cellfun(@(x)str2num(strrep(strrep(x,'firings_p',''), '.mda','')), csName_file);
    vlSelect = viPrmset >= 1 & viPrmset <= nPrmset; 
    nOutput = sum(vlSelect);
    vrDatenum_file = [vS_dir.datenum];
    vrDatenum_file = sort(vrDatenum_file(vlSelect));
    
    t_passed = range(vrDatenum_file) * 24 * 60;
    nRec = numel(csDir_rec);
    nOutput_total = nRec * nPrmset;   
    t_left = t_passed * nOutput_total / nOutput - t_passed;
    
%     cellfun(@(x)fprintf('\t%s\n',x), csPath_file(vlSelect));

    fprintf('%s on %s:\n\t%d recordings, %d parameters\n\t%d/%d (%0.1f%%) completed (%0.1f min passed, %0.1f min remaining)\n', ...
        vcDir_rec, vcFile_prmset, nRec, nPrmset, nOutput, nOutput_total, nOutput/nOutput_total*100, t_passed, t_left);
    
    
    % plot
    if isempty(vrDatenum_file), return; end
    if isempty(hAx)
        hFig = figure('Color','w','Name', [vcDir_rec, ' - ', vcFile_prmset, ' (click to update)']);
        hAx = axes(hFig);
    end
    cla(hAx);
    plot(hAx, (vrDatenum_file - vrDatenum_file(end))*24*60, nOutput_total - (1:nOutput), '.-'); 
    vcTitle = sprintf('%d/%d (%0.1f%%) completed (%0.1f min passed, %0.1f min remaining)', ...
        nOutput, nOutput_total, nOutput/nOutput_total*100, t_passed, t_left);
    xylabel_(hAx, 'Time (minutes)', '# param left', vcTitle); 
    grid(hAx, 'on'); axis(hAx, 'tight');
    ylim(hAx, [0, nOutput_total]);
    set(hAx, 'ButtonDownFcn', @(h,e)optimize_status_(vcArg1, vcFile_prmset, hAx));
    close_(hMsg);
catch
    fprintf(2, '%s\n', lasterr());
end
end %func


%--------------------------------------------------------------------------
function close_(h)
try
    close(h);
catch
end
end %func


%--------------------------------------------------------------------------
function optimize_clear_(vcDir_rec, vcFile_prmset)
% usage
% -----
% optimize_clear_(vcFile_list, vcFile_prmset)
% optimize_clear_(vcDir_rec, vcFile_prmset)

if nargin<2, vcFile_prmset=[]; end
if isempty(vcFile_prmset), vcFile_prmset = read_cfg_('ironclust_prmset'); end

vcSorter = lower(strrep(vcFile_prmset, '.prmset', ''));
if exist_file_(vcDir_rec)
    vcFile_mat = strrep(vcDir_rec, '.txt', sprintf('_scores_prmset_%s.mat', vcSorter));
    csDir_rec = load_batch_(vcDir_rec);
    vcDir_rec=fileparts(vcDir_rec); 
elseif exist_dir_(vcDir_rec)
    csDir_rec = {vcDir_rec};
    vcFile_mat = fullfile(vcDir_rec, sprintf('scores_prmset_%s.mat', vcSorter));
end
assert(exist_file_(vcFile_prmset) && exist_dir_(vcDir_rec), 'file or dir does not exist');

vS_dir = cellfun_(@(x)dir(fullfile(x, vcSorter, 'firings_p*')), csDir_rec);
vS_dir = cat(1, vS_dir{:}); 
csFiles_remove = arrayfun_(@(x)fullfile(x.folder, x.name), vS_dir);

try
    parfor ic=1:numel(csDir_rec)
        delete_(fullfile(csDir_rec{ic}, vcSorter, 'firings_p*'));
    end
catch
    for ic=1:numel(csDir_rec)
        delete_(fullfile(csDir_rec{ic}, vcSorter, 'firings_p*'));
    end    
end
  
delete_(vcFile_mat);
disp(csFiles_remove(:));
fprintf(2, 'Deleted %d previous outputs\n', numel(csFiles_remove));
end %func


%--------------------------------------------------------------------------
% 7/21/2018 JJJ: rejecting directories, strictly search for flies
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile, fVerbose)
if nargin<2, fVerbose = 0; end
if isempty(vcFile)
    flag = false; 
elseif iscell(vcFile)
    flag = cellfun(@(x)exist_file_(x, fVerbose), vcFile);
    return;
else
    S_dir = dir(vcFile);
    if numel(S_dir) == 1
        flag = ~S_dir.isdir;
    else
        flag = false;
    end
end
if fVerbose && ~flag
    fprintf(2, 'File does not exist: %s\n', vcFile);
end
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
function import_recordings_(vcFile_bin, vcFile_prb, vcDir_out)

if nargin<3, vcDir_out=''; end
if nargin<2, vcFile_prb=''; end
if nargin<1, vcFile_bin = []; end

if isempty(vcFile_bin)
    csHelp = {'Usages', '-----', ...
    '    irc2 import-recordings [path_to_bin_file] [path_to_prb_file] (output_path)`', ...    
    '    irc2 import-recordings [path_to_*bin_files] [path_to_prb_file] (output_path)`', ...
    '    irc2 import-recordings [path_to_bin_list_txt] [path_to_prb_file] (output_path)`'};
    disp_cs_(csHelp);
    return; 
end

% batch mode
if isTextFile_(vcFile_bin)
    csFiles_bin = load_batch_(vcFile_bin);
%     if isempty(vcDir_out), vcDir_out = fileparts(vcFile_bin); end
else
    csFiles_bin = list_files_(vcFile_bin, 1);
end
if numel(csFiles_bin)==0, fprintf(2, 'No files were found.\n'); return; end

vlExist = cellfun(@(x)exist_file_(x), csFiles_bin);
if ~all(vlExist)
    fprintf(2, 'Following files are not found:');
    arrayfun_(@(x)fprintf(2, '\t%s\n', x), csFiles_bin(~vlExist));
    return;
end

% write params.json
vcRecordingType = recording_type_(csFiles_bin{1});
switch vcRecordingType
    case 'spikeglx'
        vcFile_meta = strrep(csFiles_bin{1}, '.bin', '.meta');
        assert(exist_file_(vcFile_meta), sprintf('%s does not exist', vcFile_meta));
        S_meta = read_meta_spikeglx_(vcFile_meta);
    case 'rhd', S_meta = rhd_header_(csFiles_bin{1});   
    case 'neuroscope', S_meta = neuroscope_header_(csFiles_bin{1});           
end
S_json = struct('spike_sign', -1, 'samplerate', S_meta.sRateHz, 'scale', S_meta.uV_per_bit);
if isempty(vcFile_prb) && ~isempty(get_(S_meta, 'vcProbe'))
    vcFile_prb = [S_meta.vcProbe, '.prb']; 
end
[~, vcProbe, ~] = fileparts(vcFile_prb);
if isempty(vcDir_out), vcDir_out = fullfile(fileparts(vcFile_bin), vcProbe); end

mkdir_(vcDir_out);
struct2json_(S_json, fullfile(vcDir_out, 'params.json'));

% write geom.csv (probe file)
S_prb = load_prb_(vcFile_prb);    
assert(~isempty(S_prb), sprintf('Probe file is not found: %s', vcFile_prb));
csvwrite(fullfile(vcDir_out, 'geom.csv'), S_prb.mrSiteXY);
csvwrite(fullfile(vcDir_out, 'group.csv'), S_prb.viShank_site(:)');
fprintf('Wrote to %s\n', fullfile(vcDir_out, 'geom.csv'));

% write raw.mda
switch vcRecordingType
    case {'spkleglx', 'neuroscope'}
        append_bin2mda_(csFiles_bin, fullfile(vcDir_out, 'raw.mda'),...
            S_meta.nChans, S_meta.vcDataType, S_prb.viSite2Chan);
    case 'rhd' %  join channels
        mr = []; 
        assert(numel(csFiles_bin)==S_meta.nChans, 'Number of files must be equal to the number of channels');
        for iFile = 1:numel(csFiles_bin)
            vr_ = load_bin_(csFiles_bin{iFile}, S_meta.vcDataType);
            if isempty(mr)
                mr = zeros(numel(vr_), numel(csFiles_bin), 'like', vr_);
            end
            mr(:,iFile) = vr_;
        end
        writemda_(fullfile(vcDir_out, 'raw.mda'), mr(:,S_prb.viSite2Chan)');
end
fprintf('\nWrote to %s\n', fullfile(vcDir_out, 'raw.mda'));
if numel(csFiles_bin)>1
    save_cs_(fullfile(vcDir_out, 'bin_files.txt'), csFiles_bin);
    fprintf('\tCombined %d files: %s\n', numel(csFiles_bin), fullfile(vcDir_out, 'bin_files.txt')); 
end
end %func


%--------------------------------------------------------------------------
function append_bin2mda_(csFiles_bin, vcFile_mda, nChans, vcDataType, viSite2Chan)
% viSite2Chan: channel reordering
CACHE_SIZE = 1e9; % half a gig

if nargin<3, nChans = 1; end
if nargin<4, vcDataType = 'uint8'; end
if nargin<5, viSite2Chan=[]; end

bytes_per_sample = bytesPerSample_(vcDataType) * nChans;
nSamples_copy = floor(CACHE_SIZE/bytes_per_sample); % number of sammples to copy at a time

fid_mda = writemda_fid(vcFile_mda);
for iFile = 1:numel(csFiles_bin)    
    t1=tic;
    assert(exist_file_(csFiles_bin{iFile}), 'file must exist');
    fid_r = fopen(csFiles_bin{iFile}, 'r');
    b1 = ftell(fid_mda);
    while ~feof(fid_r) % read one sample at a time        
        mr_ = fread(fid_r, nChans*nSamples_copy, ['*', vcDataType]);
        if mod(numel(mr_), nChans) > 0
            nSamples_ = floor(numel(mr_) / nChans);
            if nSamples_==0, break; end
            mr_ = reshape(mr_(1:nChans * nSamples_), nChans, nSamples_);
        else
            mr_ = reshape(mr_, nChans, []);
        end
        if ~isempty(viSite2Chan)
            writemda_fid(fid_mda, mr_(viSite2Chan,:));    
        else
            writemda_fid(fid_mda, mr_);
        end
    end
    [t1, b1] = deal(toc(t1), ftell(fid_mda)-b1);
    fprintf('\t%d/%d: took %0.1fs (%0.1f MB/s)\n', ...
        iFile, numel(csFiles_bin), t1, b1/1e6/t1);
    fclose(fid_r);
end

writemda_fid(fid_mda, 'close');
end %func


%--------------------------------------------------------------------------
function save_cs_(vcFile, cs)
% display cell string
try
    fid = fopen(vcFile, 'w');
    if numel(cs)>1
        cellfun(@(s)fprintf(fid, '%s\n',s), cs(1:end-1));
    end
    if numel(cs)>=1
        fprintf(fid, '%s', cs{end});
    end
    fclose(fid);
catch ME
    fprintf(2, 'File write error: %s\n', vcFile);
end
end %func


%--------------------------------------------------------------------------
% 2018/6/26 JJJ
function vcStr = struct2json_(S, vcFile_json)
if nargin<2, vcFile_json = ''; end

csName = fieldnames(S);
vcStr = '{';
for iField = 1:numel(csName)
    vcStr_ = sprintf('\t"%s": %s', csName{iField}, field2str_(S.(csName{iField}),true));
    if iField < numel(csName)
        vcStr = sprintf('%s\n%s,', vcStr, vcStr_);
    else
        vcStr = sprintf('%s\n%s\n}', vcStr, vcStr_);
    end
end

if ~isempty(vcFile_json)
    fid = fopen(vcFile_json, 'w');
    fprintf(fid, '%s', vcStr);
    fclose(fid);
    fprintf('Wrote to %s\n', vcFile_json);
end
end


%---------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function vcStr = field2str_(val, fDoubleQuote)
% convert a value to a strong
if nargin<2, fDoubleQuote = false; end

switch class(val)
    case {'int', 'int16', 'int32', 'uint16', 'uint32'}
        vcFormat = '%d';
    case {'double', 'single'}
        vcFormat = '%g';
        if numel(val)==1
            if mod(val(1),1)==0, vcFormat = '%d'; end
        end
    case 'char'
        if fDoubleQuote
            vcStr = sprintf('"%s"', val);
        else
            vcStr = sprintf('''%s''', val);
        end
        return;
    case 'cell'
        vcStr = '{';
        for i=1:numel(val)
            vcStr = [vcStr, field2str_(val{i})];
            if i<numel(val), vcStr = [vcStr, ', ']; end
        end
        vcStr = [vcStr, '}'];
        return;
    case 'logical'
        vcFormat = '%s';
        if val
            vcStr = '1';
        else
            vcStr = '0';
        end
        return;
    otherwise
        vcStr = '';
        fprintf(2, 'field2str_: unsupported format: %s\n', class(val));
        return;
end

if numel(val) == 1
    vcStr = sprintf(vcFormat, val);
else % Handle a matrix or array
    vcStr = '[';
    for iRow=1:size(val,1)
        for iCol=1:size(val,2)
            vcStr = [vcStr, field2str_(val(iRow, iCol))];
            if iCol < size(val,2), vcStr = [vcStr, ', ']; end
        end
        if iRow<size(val,1), vcStr = [vcStr, '; ']; end
    end
    vcStr = [vcStr, ']'];
end
end %func


%--------------------------------------------------------------------------
function val = read_cfg_(vcName, fVerbose)
% read configuration file that stores path to folder
% load from default.cfg but override with user.cfg if it exists
if nargin<2, fVerbose = 0; end

S_cfg = file2struct_(ircpath_('default.cfg'));
vcFile_user = ircpath_(S_cfg.user_cfg);
if exist_file_(vcFile_user)
    S_cfg = struct_merge_(S_cfg, file2struct_(vcFile_user)); %, {'path_dropbox', 'path_backup', 'default_prm'});
    if fVerbose, fprintf('Configuration overrode from user.cfg\n'); end
else
    if fVerbose, fprintf('Configuration loaded from default.cfg\n'); end
end

if nargin==0
    val = S_cfg;
else
    try
        val = S_cfg.(vcName);
    catch
        disperr_(['read_cfg_: error reading ', ircpath_('default.cfg')]);
        switch lower(vcName)
            case 'default_prm'
                val = 'default.prm';
            otherwise
                val = [];
        end
    end
end
end %func


%--------------------------------------------------------------------------
% used by outside
function [trWav_clu, viSite_clu, mrWavCor] = clu_wav_(S0, S_auto, trPc_spk, trPc2_spk)
if nargin<2, S_auto = []; end
if nargin<3, trPc_spk=[]; end
if nargin<4, trPc2_spk=[]; end

if isempty(S_auto), S_auto = S0.S_auto; end
if isempty(trPc_spk), trPc_spk = load_fet_(S0, S0.P, 1); end
[P, viSite_spk, viSite2_spk, mrPv_global] = ...
    get_(S0, 'P', 'viSite_spk', 'viSite2_spk', 'mrPv_global');
if isempty(viSite2_spk)
    cviSite_spk_fet = {viSite_spk};
    ctrPc_spk_fet = {trPc_spk};
else
    cviSite_spk_fet = {viSite_spk, viSite2_spk};
    if isempty(trPc2_spk), trPc2_spk = load_fet_(S0, S0.P, 2); end
    ctrPc_spk_fet = {trPc_spk, trPc2_spk};
end
[cviSite_clu_fet, ctrWav_clu_fet] = deal(cell(size(cviSite_spk_fet)));
cviSpk_clu = vi2cell_(S_auto.viClu, S_auto.nClu);
for iFet = 1:numel(cviSite_spk_fet)
    viSite_spk = cviSite_spk_fet{iFet};
    cviSite_clu = cellfun_(@(x)mode(viSite_spk(x)), cviSpk_clu);    
    cviSpk1_clu = cellfun_(@(x,y)x(viSite_spk(x)==y), cviSpk_clu, cviSite_clu);
    cviSite_clu_fet{iFet} = cell2mat(cviSite_clu);
    trPc_spk = ctrPc_spk_fet{iFet};
    cmrWav_clu = cellfun_(@(x)mrPv_global * mean(trPc_spk(:,:,x),3), cviSpk1_clu);
    ctrWav_clu_fet{iFet} = cat(3, cmrWav_clu{:});
end
[viSite_clu, trWav_clu] = deal(cviSite_clu_fet{1}, ctrWav_clu_fet{1});
if isempty(viSite2_spk)
    [viSite2_clu, trWav2_clu] = deal([]);
else
    [viSite2_clu, trWav2_clu] = deal(cviSite_clu_fet{2}, ctrWav_clu_fet{2});
end

mrWavCor = zeros(S_auto.nClu, 'single');
nShift = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
for iClu=1:S_auto.nClu
    [iSite1, mrWav_clu1] = deal(viSite_clu(iClu), trWav_clu(:,:,iClu));
    viClu2 = find(viSite_clu==iSite1);
    if ~isempty(viClu2)        
        vrWavCor2 = wavcor_(mrWav_clu1, trWav_clu(:,:,viClu2), nShift);
        mrWavCor(viClu2,iClu) = max(mrWavCor(viClu2,iClu), vrWavCor2);
    end
    if ~isempty(viSite2_spk)
        [iSite2, mrWav2_clu1] = deal(viSite2_clu(iClu), trWav2_clu(:,:,iClu));
        viClu2A = find(viSite_clu==iSite2);
        if ~isempty(viClu2A)
            vrWavCor2A = wavcor_(mrWav2_clu1, trWav_clu(:,:,viClu2A), nShift);
            mrWavCor(viClu2A,iClu) = max(mrWavCor(viClu2A,iClu), vrWavCor2A);
        end
        viClu2B = find(viSite2_clu==iSite1);
        if ~isempty(viClu2B)
            vrWavCor2B = wavcor_(mrWav_clu1, trWav2_clu(:,:,viClu2B), nShift);
            mrWavCor(viClu2B,iClu) = max(mrWavCor(viClu2B,iClu), vrWavCor2B);
        end
    end
end
mrWavCor = max(mrWavCor, mrWavCor');
end %func


%--------------------------------------------------------------------------
function [S0, S_auto, trPc_spk, trPc2_spk, P] = load_irc2_(vcFile_prm)

S0 = load0_(vcFile_prm);
if isempty(S0), error('%s: output is not found', vcFile_prm); end
S_auto = get_(S0, 'S_auto');
assert(~isempty(S_auto), 'S_auto does not exist');
P = S0.P;
trPc_spk = load_fet_(S0, P, 1);
trPc2_spk = load_fet_(S0, P, 2);

% export valid clusters only
vlKeep = S_auto.viClu>0;
[S0.vrAmp_spk, S0.viTime_spk, S0.viSite_spk, S_auto.viClu] = ...
    deal(S0.vrAmp_spk(vlKeep), S0.viTime_spk(vlKeep), S0.viSite_spk(vlKeep), S_auto.viClu(vlKeep));
if ~isempty(get_(S0, 'viSite2_spk')), S0.viSite2_spk = S0.viSite2_spk(vlKeep); end
trPc_spk = trPc_spk(:,:,vlKeep);
if ~isempty(trPc2_spk), trPc2_spk = trPc2_spk(:,:,vlKeep); end
end %fnc


%--------------------------------------------------------------------------
% 10/18/2018 JJJ: Cluster order fixed
function S_auto = S_auto_sort_(S_auto, vcField_sort, viSite_spk)
% sort clusters by the centroid position
% vcField_sort: {'', 'vrPosY_clu + vrPosX_clu'}

if nargin<2, vcField_sort = ''; end

% Sort clusters by its sites
if isempty(vcField_sort), vcField_sort = 'viSite_clu'; end

switch vcField_sort
    case 'vrPosY_clu + vrPosX_clu'
        [~, viCluSort] = sort(S_auto.vrPosY_clu + S_auto.vrPosX_clu, 'ascend');
    otherwise
        [~, viCluSort] = sort(S_auto.(vcField_sort), 'ascend');
end
S_auto.viClu = mapIndex_(S_auto.viClu, viCluSort); % fixed
S_auto = struct_reorder_(S_auto, viCluSort, ...
    'cviSpk_clu', 'vrPosX_clu', 'vrPosY_clu', 'vnSpk_clu', 'viSite_clu', 'cviTime_clu', 'csNote_clu');
S_auto = S_clu_refresh_(S_auto, 1, viSite_spk);
end %func


%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_refresh_(S_clu, fRemoveEmpty, viSite_spk)

if nargin<2, fRemoveEmpty=1; end
nClu = double(max(S_clu.viClu));
S_clu.nClu = nClu;
S_clu.cviSpk_clu = vi2cell_(S_clu.viClu, nClu);
S_clu.vnSpk_clu = cellfun(@numel, S_clu.cviSpk_clu); 
if ~isempty(viSite_spk)
    S_clu.viSite_clu = double(arrayfun(@(iClu)mode(viSite_spk(S_clu.cviSpk_clu{iClu})), 1:nClu));
end
if fRemoveEmpty, [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu); end
end %func


%--------------------------------------------------------------------------
function [vi, nClu, viA, viAB] = mapIndex_(vi, viA, viB)
% change the index of vi according to the map (viA)

if nargin<2, viA = setdiff(unique(vi), 0); end %excl zero
if nargin<3, viB = 1:numel(viA); end
if isempty(viA), viA = 1:max(vi); end
nClu = viB(end);
viAB(viA) = viB; %create a translation table A->B
vl = vi>0;
vi(vl) = viAB(vi(vl)); %do not map zeros
end %func


%--------------------------------------------------------------------------
function S = struct_reorder_(S, viKeep, varargin)
for i=1:numel(varargin)
    try
        vcVar = varargin{i};
        if ~isfield(S, vcVar), continue; end %ignore if not
        vr1 = S.(vcVar);
        if isvector(vr1)
            vr1 = vr1(viKeep);
        elseif ismatrix(vr1)
            vr1 = vr1(viKeep, :);
        else
            vr1 = vr1(viKeep, :, :);
        end
        S.(vcVar) = vr1;
    catch
        ;
    end
end
end %func


%--------------------------------------------------------------------------
function [viSite2_spk, viTime2_spk] = find_site2_spk_(trPc_spk, mrPv, viSite_spk, viTime_spk, P)
% find second min, excl local ref sites

nSites_fet = P.nSites_fet;

% vrDist_site = pdist(P.mrSiteXY);
% vrDist_site_uniq = unique(pdist(P.mrSiteXY));
% mrDist_site = squareform(vrDist_site);
% median(sum(mrDist_site <= vrDist_site_uniq(1),1))
% median(sum(mrDist_site <= vrDist_site_uniq(2),1))
% median(sum(mrDist_site <= vrDist_site_uniq(3),1))
% median(sum(mrDist_site <= vrDist_site_uniq(8),1))

% denoised waveform
nT_wav = abs(P.spkLim(1))*2+1;
trWav_spk = pc2wav_(mrPv(1:nT_wav,:), trPc_spk);

[mrMin_spk, miMin_spk] = min(trWav_spk(:,2:nSites_fet,:),[],1);
miMin_spk = permute(miMin_spk, [2,3,1]);
[~, viMin_spk] = min(mrMin_spk,[],2);

miSites2 = P.miSites(2:nSites_fet, viSite_spk);
viSite2_spk = int32(mr2vr_sub2ind_(miSites2, viMin_spk(:), []));    

iT_peak = 1 - P.spkLim(1);
viTime2_offset_spk = mr2vr_sub2ind_(miMin_spk, viMin_spk(:), []) - iT_peak;
viTime2_spk = viTime_spk + int64(viTime2_offset_spk(:));
end %func


%--------------------------------------------------------------------------
% Remove leading singular dimension
% 12/15/17 JJJ: squeeze out specific dimension
% 7/26/17 JJJ: code cleanup and testing
function val = squeeze_(val, idimm)
% val = squeeze_(val) : when squeezeing matrix, transpose if leading dimm is 1
% val = squeeze_(val, idimm): permute specified dimension out
size_ = size(val);
if nargin>=2
    dimm_ = [setdiff(1:ndims(val), idimm), idimm];
    val = permute(val, dimm_);
elseif numel(size_)==2 && size_(1) == 1
    val = val';
else
    val = squeeze(val);
end
end


%--------------------------------------------------------------------------
function vr = mr2vr_sub2ind_(mr, vi1, vi2)
if isempty(mr), vr = []; return; end
if nargin<3, vi2=[]; end
if isempty(vi1), vi1 = 1:size(mr,1); end
if isempty(vi2), vi2 = 1:size(mr,2); end
vr = mr(sub2ind(size(mr), vi1(:), vi2(:)));
end %func


%--------------------------------------------------------------------------
function [viSite_spk2] = find_site2_spk__(tnWav_spk, viSite_spk, P)
% find second min, excl local ref sites

fUse_min = 1;
imin0 = 1 - P.spkLim(1);
nSites_spk = size(P.miSites, 1);
viSites2 = 2:P.nSites_fet;
miSites2 = P.miSites(viSites2, viSite_spk);
tnWav_spk2 = tnWav_spk(:,viSites2,:);
if fUse_min
    mnMin_spk = squeeze_(min(tnWav_spk2));
else
%     mnMin_spk = -squeeze_(std(single(tnWav_spk(:,viSites2,:))));
    mnMin_spk = squeeze_(min(tnWav_spk2) - max(tnWav_spk2)); % use Vpp to determine second peak site
end
[~, viSite_spk] = min(mnMin_spk);
viSite_spk2 = int32(mr2vr_sub2ind_(miSites2, viSite_spk, []));
    
end %func


%--------------------------------------------------------------------------
% Call from irc.m
function cout = call_irc_(dbstack1, cell_input, nargout)
vcFunc = dbstack1(1).name;
try
    switch nargout
        case 0, cout{1} = []; irc('call', vcFunc, cell_input);
        case 1, cout{1} = irc('call', vcFunc, cell_input);
        case 2, [cout{1}, cout{2}] = irc('call', vcFunc, cell_input);
        case 3, [cout{1}, cout{2}, cout{3}] = irc('call', vcFunc, cell_input);
        case 4, [cout{1}, cout{2}, cout{3}, cout{4}] = irc('call', vcFunc, cell_input);
        otherwise, error('call_irc_: undefined func: %s', vcFunc);
    end
catch ME
    fprintf(2, 'call_irc_: %s\n', ME.message);
    rethrow(ME);
end
end %func


%--------------------------------------------------------------------------
function [miSites, nSites_fet] = findNearSites_(mrSiteXY, P)
% find nearest sites
% if nargin<2, maxSite = []; end
% if nargin<3, viSiteZero = []; end
% if nargin<4, viShank_site = []; end
% if numel(unique(viShank_site)) <= 1, viShank_site = []; end
if nargin<2, P=[]; end
if ~isempty(P)
    [maxDist_site_um, maxDist_site_spk_um, viSiteZero, viShank_site] = ...
        struct_get_(P, 'maxDist_site_um', 'maxDist_site_spk_um', 'viSiteZero', 'viShank_site');
else
    [maxDist_site_um, maxDist_site_spk_um, viSiteZero, viShank_site] = deal([]);
end
if ~isempty(viSiteZero), mrSiteXY(viSiteZero,:) = inf; end
% max_dist = max(pdist(mrSiteXY));
nSites = size(mrSiteXY,1);
if ~isempty(viShank_site)
    [vi_uniq, vn_uniq] = unique_count_(viShank_site);    
    nSites_shank = min(vn_uniq);
else
    nSites_shank = nSites;
end
nSites_spk = nSites_within_radius_(mrSiteXY, maxDist_site_spk_um);
nSites_spk = min(nSites_spk, nSites_shank);

miSites = zeros(nSites_spk, nSites);
for iSite=1:nSites
    vrSiteDist = pdist2_(mrSiteXY(iSite,:), mrSiteXY);
    if ~isempty(viShank_site)
        vrSiteDist(viShank_site~=viShank_site(iSite)) = inf;
    end
    [vrSiteDist, viSrt] = sort(vrSiteDist, 'ascend');
    miSites(:,iSite) = viSrt(1:nSites_spk);
end

if ~isempty(maxDist_site_um)
    nSites_fet = nSites_within_radius_(mrSiteXY, maxDist_site_um);
    nSites_fet = min(nSites_fet, nSites_spk);
else
    nSites_fet = nSites_spk;
end
end %func


%--------------------------------------------------------------------------
function nSites = nSites_within_radius_(mrSiteXY, radius)
mrDist_site = squareform(pdist(mrSiteXY));
nSites = floor(median(sum(mrDist_site <= radius,1)));
end %func


%--------------------------------------------------------------------------
% export .prb files to .json (flatiron-0.1) format
function S_json = export_prb_json_(vcFile_prb, vcDir_out)

S_cfg = read_cfg_();
format_version = S_cfg.flatiron_probe_version;
probe_type = S_cfg.flatiron_probe_default;

% batch export
if any(vcFile_prb=='*')
    csFile_prb = list_files_(vcFile_prb);
    S_json = cellfun_(@(x)export_prb_json_(x, vcDir_out), csFile_prb);
    return;
end

[~, probe_name, vcExt_] = fileparts(vcFile_prb);
assert(strcmpi(vcExt_, '.prb'), 'must be a valid probe format');

S_prb = load_prb_(vcFile_prb);
assert(~isempty(S_prb), 'probe file must exist');
nSites = size(S_prb.mrSiteXY,1);
[x, y] = deal(S_prb.mrSiteXY(:,1)', S_prb.mrSiteXY(:,2)');
z = zeros(size(x));

if ~isempty(get_(S_prb, 'viSite2Chan'))
    channel = S_prb.viSite2Chan(:)' - 1; %zero-base
else
    channel = 0:nSites-1;
end
if ~isempty(get_(S_prb, 'vrSiteHW'))
    [site_height, site_width] = deal(S_prb.vrSiteHW(1), S_prb.vrSiteHW(2));
else
    [site_height, site_width] = deal([]);
end
group = get_(S_prb, 'viShank_site');
if isempty(group)
    group = zeros(nSites,1);
elseif min(group) > 0
    group = group - 1;
end
group = group(:)';

% write to a file
S_json = makeStruct_(format_version, probe_name, probe_type, channel, x, y, z, group, site_width, site_height);
mkdir_(vcDir_out);
vcFile_json = fullfile(vcDir_out, [probe_name, '.json']);
struct2json_(S_json, vcFile_json);
if nargout==0, edit_(vcFile_json); end
end %func


%--------------------------------------------------------------------------
function probe_image_(vcFile_img, vcShankPitch, vcSitePitch)
% vcFile_img: file name of probe image
% vcShankPitch: nominal distance between shank tips in microns
% vcSitePitch: nominal vertical distance between adjacent sites in microns

if nargin<2, vcShankPitch=''; end
if nargin<2, vcSitePitch=''; end

assert(exist_file_(vcFile_img), 'probe image file not found');
imtool(vcFile_img);
S_cfg = read_cfg_();
color_prbimg = get_(S_cfg, 'color_prbimg', [0 0 0]);
img = imread(vcFile_img);
BW = img(:,:,1)==color_prbimg(1) & img(:,:,2)==color_prbimg(2) & img(:,:,3)==color_prbimg(3);
s=regionprops(bwlabel(BW,8), 'centroid'); 
mrXY = cat(1, s.Centroid);
figure('Color','w','Name',vcFile_img);
imshow(~BW);
hold on
plot(mrXY(:,1), mrXY(:,2), 'go-');

% shank pitch and site pitch scaling
mrXY1 = mrXY;
mrXY1(:,2) = max(mrXY1(:,2))-mrXY1(:,2);
viSite_tip = find(mrXY1(:,2)<.5); % sites at the bottom tips
[~, imin] = min(mrXY1(viSite_tip,1)); iSite0 = viSite_tip(imin(1));
mrXY1(:,1) = mrXY1(:,1) - mrXY1(iSite0,1);
plot(mrXY(iSite0,1), mrXY(iSite0,2), 'g*-');

if ~isempty(vcShankPitch) % scale x
    shankPitch = mean(diff(mrXY1(viSite_tip,1)));
    mrXY1(:,1) = mrXY1(:,1) / shankPitch * str2num(vcShankPitch);
end

if ~isempty(vcSitePitch)
    sitePitch = mean(diff(unique(round(mrXY1(:,2)))));
    mrXY1(:,2) = mrXY1(:,2) / sitePitch * str2num(vcSitePitch);
end

mrXY1 = round(mrXY1*2)/2; 
csText_xy = arrayfun_(@(x,y)sprintf('%0.1f, %0.1f',x,y), mrXY1(:,1), mrXY1(:,2));
text(mrXY(:,1), mrXY(:,2), csText_xy, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','r');

end %func


%--------------------------------------------------------------------------
function trim_mda_(vcFile_in, vcFile_out, nSamples_copy, iOffset)
if nargin<4, iOffset=''; end

if ischar(nSamples_copy), nSamples_copy = str2num(nSamples_copy); end
if isempty(iOffset), iOffset=0; end

% vcTrim: str2num samples
t_fun = tic;
[S_mda, fid_r] = readmda_header_(vcFile_in);
fid_mda = writemda_fid(vcFile_out);
[vcDataType, dimm] = get_(S_mda, 'vcDataType', 'dimm');
nChans = dimm(1);

b1 = ftell(fid_r);
mr_ = fread_(fid_r, [nChans, nSamples_copy], vcDataType);
writemda_fid(fid_mda, mr_);
b1 = ftell(fid_r) - b1;
t_fun = toc(t_fun);

fclose(fid_r);
writemda_fid(fid_mda, 'close');
fprintf('\ttook %0.1fs (%0.1f MB/s)\n', t_fun, b1/1e6/t_fun);
end %func


%--------------------------------------------------------------------------
% 11/6/2018 JJJ: compute hash
% uses DataHash.m
% https://www.mathworks.com/matlabcentral/fileexchange/31272-datahash
function [vcHash, csHash] = file2hash_(csFiles)
% Usage
% ----
% file2hash_(vcFile)
% file2hash_(csFiles)
if nargin<1, csFiles = []; end
if isempty(csFiles)
    csFiles = [mfilename('fullpath'), '.m'];
end
if ischar(csFiles), csFiles = {csFiles}; end
csHash = cell(size(csFiles));
for iFile = 1:numel(csFiles)
    csHash{iFile} = DataHash(csFiles{iFile}, struct('Method', 'MD5', 'Input', 'file', 'Format', 'hex'));
    if iFile == 1
        vcHash = csHash{iFile};
    else
        vcHash = bitxor(vcHash, csHash{iFile});
    end
end
end %func


%--------------------------------------------------------------------------
% 11/6/2018 JJJ: compute hash
% uses DataHash.m
% https://www.mathworks.com/matlabcentral/fileexchange/31272-datahash
function vcHash = struct2hash_(S)
% Usage
% ----
% file2hash_(vcFile)
% file2hash_(csFiles)
vcHash = DataHash(S, struct('Method', 'MD5', 'Input', 'array', 'Format', 'hex'));
end %func


%--------------------------------------------------------------------------
% 8/22/18 JJJ: changed from the cell output to varargout
% 9/26/17 JJJ: Created and tested
function S_copy = struct_copy_(varargin)
% Obtain a member of struct
% cvr = cell(size(varargin));
S = varargin{1};
S_copy = struct();
for iArg=2:nargin
    csName1 = varargin{iArg};
    if ischar(csName1)
        csName1 = {csName1};
    elseif ~iscell(csName1)
        continue; %skip
    end
    for iCell=1:numel(csName1)
        vcName1 = csName1{iCell};
        if isfield(S, vcName1)
            S_copy.(vcName1) = S.(vcName1);
        else
            S_copy.(vcName1) = [];
        end
    end
end %for
end %func


%--------------------------------------------------------------------------
function S_txt = loadjson_(vcArg_txt)
S_txt=[]; 
if ~exist_file_(vcArg_txt), return; end
addpath_('jsonlab-1.5/');
S_txt = loadjson(vcArg_txt);
end %func


function varargout = git_pull_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = frewind_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = disperr_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = edit_prm_file_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = edit_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = text3_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = text_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = create_figure_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_prb_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = isTextFile_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_batch_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = list_files_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_merge_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = read_cfg_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = file2struct_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = loadjson_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = filesize_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = car_reject_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = struct_copy_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = cast_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_default_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_filter_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = meanSubt_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = memory_matlab_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = recording_duration_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_set_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = find_site_spk23_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = matchFileExt_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_copyas_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = set_bool_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = file2hash_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = dir_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = S_clu_refresh_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = map_index_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = mr2thresh_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = findNearSites_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = shift_range_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = fopen_mda_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = fopen_nsx_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_probe_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct2file_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = addpath_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
