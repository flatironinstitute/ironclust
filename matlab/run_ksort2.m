function run_ksort2(vcDir_in, vcDir_out, arg_fname)
% usage
% -----
% run_ksort2(command)
% run_ksort2(vcDir_in, vcDir_out, vcFile_template)
% run_ksort2('mcc', vcDir_ksort2, vcDir_out)
%    compile run_ksort2
%
% arguments
% -----
% vcDir_in: input directory
% vcDir_out: output directory
% vcFile_template: template file (optional)

if nargin<1, error('provide vcDir_in'); end
if nargin<2, vcDir_out = ''; end
if nargin<3, arg_fname = []; end

if strcmpi(vcDir_in, 'mcc')
    [vcDir_ksort2, vcDir_mcc_out] = deal(vcDir_out, arg_fname);
    mcc_ksort2_(vcDir_ksort2, vcDir_mcc_out);
    return;
end

if isempty(arg_fname), arg_fname = 'ksort2.txt'; end
if isempty(vcDir_out)
    vcDir_out = strrep(vcDir_in, 'groundtruth', 'ksort2');
end

source_path = fileparts(mfilename('fullpath'));
ironclust_src = fileparts(source_path);
if ~exist_file_(arg_fname)
    arg_fname = fullfile(source_path, arg_fname);
end

if ~isdeployed()
    addpath(genpath(fullfile(source_path))); 
    S_cfg = file2struct(fullfile(source_path, 'default.cfg'));
    kilosort_src = S_cfg.path_ksort2;
else
    kilosort_src = '';
    fprintf('Running run_ksort2 in deployed mode.\n');
end

    
% inferred from the path
firings_out_fname = fullfile(vcDir_out, 'firings_out.mda');
raw_fname = fullfile(vcDir_in, 'raw.mda');

if ~exist_file_(firings_out_fname)
    irc('call', 'mkdir', {vcDir_out}); % create temp output directory
    geom_fname = fullfile(vcDir_in, 'geom.csv');
    prm_fname = fullfile(vcDir_in, 'params.json');
    p_kilosort2_(kilosort_src, ironclust_src, vcDir_out, raw_fname, geom_fname, firings_out_fname, arg_fname);
    fprintf('Clustering result wrote to %s\n', firings_out_fname);
end
if ~exist_file_(firings_out_fname), fprintf(2, 'No output file found\n'); end


% validate
vcFile_gt_mda = fullfile(vcDir_in, 'firings_true.mda');
if exist_file_(vcFile_gt_mda)
    try
        irc('validate-mda', vcFile_gt_mda, firings_out_fname, raw_fname); % assume that groundtruth file exists
    catch
        fprintf(2, 'Validation failed\n');
    end
end

% Exit
exit_deployed_();
end %func


%--------------------------------------------------------------------------
% 11/2/2018 JJJ: matlab compiler, generates run_irc
% @TODO: get the dependency list from sync_list
function mcc_ksort2_(vcDir_ksort2, vcDir_out)
if nargin<1, vcDir_ksort2 = ''; end
if nargin<2, vcDir_out = ''; end
assert(~isempty(vcDir_ksort2), 'vcDir_ksort2 must be specified');
if ~exist_func_('mcc')
   fprintf(2, 'Matlab Compiler Toolbox is not installed.\n');
   return; 
end

compile_ksort2_(vcDir_ksort2);

fprintf('Compiling run_ksort2.m\n'); t1=tic;
vcEval1 = ['mcc -m -v -R ''-nodesktop, -nosplash -nojvm'' -a ' vcDir_ksort2];
vcEval2 = ' -a ./mdaio/* -a ./jsonlab-1.5/* -a ./npy-matlab/* -a default.cfg -a irc.m run_ksort2.m';
if ~isempty(vcDir_out)
    mkdir_(vcDir_out);
    vcEval = [vcEval1, vcEval2, ' -d ', vcDir_out];
else
    vcEval = [vcEval1, vcEval2];
end
disp(vcEval);
eval(vcEval);
fprintf('\n\trun_ksort2.m is compiled by mcc, took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
% assume mexGPUall.m is already in the path
function compile_ksort2_(vcDir_ksort2)
vcPath_pwd = pwd();
try
    cd(fullfile(vcDir_ksort2, 'CUDA'));
    [~,path_nvcc_] = system('which nvcc');
    path_nvcc_ = strrep(path_nvcc_, 'nvcc', '');
    disp(['path_nvcc_: ', path_nvcc_]);
    setenv('MW_NVCC_PATH', path_nvcc_);
    run('mexGPUall.m');
catch
    disp(lasterr());
    fprintf(2, 'Problem running mexGPUall.\n');
end
cd(vcPath_pwd);
end %func


%--------------------------------------------------------------------------
% 12/13/2018 JJJ: checks if function exists
function flag = exist_func_(vcFunc)
flag = ~isempty(which(vcFunc));
end %func


%--------------------------------------------------------------------------
function exit_deployed_()
try
    if isdeployed() || ismcc(), exit(); end
catch
    ;
end
end %func


%--------------------------------------------------------------------------
% 7/21/2018 JJJ: rejecting directories, strictly search for flies
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile, fVerbose)
if nargin<2, fVerbose = 0; end
if isempty(vcFile)
    flag = 0; 
else
    S_dir = dir(vcFile);
    if numel(S_dir) == 1
        flag = ~S_dir.isdir;
    else
        flag = 0;
    end
end
if fVerbose && ~flag
    fprintf(2, 'File does not exist: %s\n', vcFile);
end
end %func


%--------------------------------------------------------------------------
% 2019/4/12 JJJ: copied from spikeforest
function mr_out = p_kilosort2_(kilosort_src, ironclust_src, temp_path, raw_fname, geom_fname, firings_out_fname, arg_fname)
% cmdstr2 = sprintf("p_ironclust('$(tempdir)','$timeseries$','$geom$','$firings_out$','$(argfile)');");

if exist(temp_path, 'dir') ~= 7
    mkdir(temp_path);
end

% prepare for kilosort execution
if ~isdeployed()
    addpath(genpath(kilosort_src));
    addpath(fullfile(ironclust_src, 'matlab'), fullfile(ironclust_src, 'matlab/mdaio'), fullfile(ironclust_src, 'matlab/npy-matlab'));    
end
ops = import_ksort_(raw_fname, geom_fname, arg_fname, temp_path);

% Run kilosort
t1=tic;
fprintf('Running kilosort on %s\n', raw_fname);
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

fprintf('\n\tfound %d good units \n', sum(rez.good>0))

fprintf('\n\ttook %0.1fs\n', toc(t1));

% Export kilosort
mr_out = export_ksort_(rez, firings_out_fname);

fprintf('Clustering result wrote to %s\n', firings_out_fname);

end %func


%--------------------------------------------------------------------------
function mr_out = export_ksort_(rez, firings_out_fname)

mr_out = zeros(size(rez.st3,1), 3, 'double'); 
mr_out(:,2) = rez.st3(:,1); %time
mr_out(:,3) = rez.st3(:,2); %cluster
writemda(mr_out', firings_out_fname, 'float64');
end %func


%--------------------------------------------------------------------------
function ops = import_ksort_(raw_fname, geom_fname, arg_fname, fpath)
% fpath: output path
fprintf('import_ksort_...\n'); t1=tic;
S_txt = irc('call', 'meta2struct', {arg_fname});
useGPU = 1;
[freq_min, pc_per_chan, sRateHz, spkTh] = struct_get_(S_txt, 'freq_min', 'pc_per_chan', 'samplerate', 'detect_threshold');
spkTh = -abs(spkTh);

% convert to binary file (int16)
fbinary = strrep(raw_fname, '.mda', '.bin');
[Nchannels, ~] = mda2bin_(raw_fname, fbinary, S_txt.detect_sign);

% create a probe file
mrXY_site = csvread(geom_fname);
vcFile_chanMap = fullfile(fpath, 'chanMap.mat');
createChannelMapFile_(vcFile_chanMap, Nchannels, mrXY_site(:,1), mrXY_site(:,2));
nChans = size(mrXY_site,1);

S_ops = makeStruct_(fpath, fbinary, nChans, vcFile_chanMap, spkTh, useGPU, sRateHz, pc_per_chan, freq_min);
ops = config_kilosort2_(S_ops); %obtain ops
fprintf('import_ksort_ took %0.1fs\n', toc(t1));
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
function [nChans, nSamples] = mda2bin_(raw_fname, fbinary, detect_sign)

mr = readmda(raw_fname);
% adjust scale to fit int16 range with a margin
if isa(mr,'single') || isa(mr,'double')
    uV_per_bit = 2^14 / max(abs(mr(:)));
    mr = int16(mr * uV_per_bit);
end
[nChans, nSamples] = size(mr);
if detect_sign > 0, mr = -mr; end % force negative detection
fid = fopen(fbinary, 'w');
fwrite(fid, mr, 'int16');
fclose(fid);
end %func


%--------------------------------------------------------------------------
function ops = config_kilosort2_(S_opt)
% S_opt: fpath, fbinary, nChans, vcFile_chanMap, spkTh, useGPU, sRateHz,
%       pc_per_chan, freq_min

% rootH = '~/kilosort';
ops.fproc       = fullfile(S_opt.fpath, 'temp_wh.dat'); % proc file on a fast SSD  ;
ops.trange = [0 Inf]; % time range to sort
ops.NchanTOT    = S_opt.nChans; % total number of channels in your recording

% the binary file is in this folder
ops.fbinary = S_opt.fbinary;

ops.chanMap = S_opt.vcFile_chanMap;
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if no chanMap file

% sample rate
ops.fs = S_opt.sRateHz;  

% frequency for high pass filtering (150)
ops.fshigh = S_opt.freq_min;   

% minimum firing rate on a "good" channel (0 to skip)d
ops.minfr_goodchannels = 0.1; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];  

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 10;  

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; 

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/50; 

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8;

% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = S_opt.spkTh;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = S_opt.useGPU; % has to be 1, no CPU version yet, sorry
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = S_opt.pc_per_chan; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

end %func


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