function ret = p_ironclust(params) % params already passing string
% question: exception handling. if file doesn't exist what to do?

% internal parameters

if (ischar(params))&&(strcmp(params,'spec'))
    ret=get_spec();
    return;
else
    ret = [];
end

P = load_params_(params); %@TODO

t1 = tic;
% Open file : @TODO: input support omre than MDA format
[S_mda_r, fid_r] = readmda_header(P.vcFile); % 
nChans = S_mda_r.dimm(1); % specify number of channels
P.nChans = nChans;

if ~P.fPaged_load
    nLoad = 1;
else
    [nLoad, nSamples_load, nSamples_last] = plan_load_(S_mda_r, P);
end

if nLoad == 1
    S0 = ironclust_nopage_(P);
else
    for iLoad = 1:nLoad
        % Read        
        fprintf('Processing %d/%d...\n', iLoad, nLoad);
        nSamples_ = ifeq_(iLoad == nLoad, nSamples_last, nSamples_load);
        if iLoad == 1
            viT_ = 1:nSamples_;
            nSamples_load_ = nSamples_ + P.nPad_filt;   
            S0 = struct();
            ctrFet = cell(nLoad, 1);
        elseif iLoad == nLoad
            viT_ = (1+P.nPad_filt):(nSamples_+P.nPad_filt);
            nSamples_load_ = nSamples_ + P.nPad_filt;
        else
            viT_ = (1+P.nPad_filt):(nSamples_+P.nPad_filt);
            nSamples_load_ = nSamples_ + 2 * P.nPad_filt;
        end
%         fprintf('\tLoading from file...'); t_load_ = tic();
        mnWav_raw = fread(fid_r, [nChans, nSamples_load_], S_mda_r.vcDataType);
%         fprintf('took %0.1fs\n', toc(t_load_));
        fseek(fid_r, -2 *nChans*P.nPad_filt*S_mda_r.nBytes_sample, 'cof');

        % run processing
        [ctrFet{iLoad}, mnWav_filt, S0] = wav2fet_(mnWav_raw, S0, P);        

        % Write
        if P.fSave_filt
            if iLoad == 1
                S_mda_w = setfield(S_mda_r, 'vcDataType', class(mnWav_filt));
                fid_w = writemda_header(params.timeseries_out, S_mda_w);    
            end
            fwrite(fid_w, mnWav_filt(:,viT_), S_mda_w.vcDataType);    
            if iLoad == nLoad, fclose(fid_w); end
        end
        
        % clear
        [mnWav_raw, mnWav_filt] = deal([]);
    end %for    
    S0 = cluster_(cell2mat(ctrFet), S0, P);
end %if 
fclose(fid_r); % reading done

% save to firings_out (mda file)
write_S0_(S0, P);

% Export to 
fprintf('Done with slopefilter, took %0.1fs (nLoad=%d)\n', ...
    toc(t1), nLoad);
end % func


%--------------------------------------------------------------------------
function [trFet, mnWav_filt, S0] = wav2fet_(mnWav_raw, S0, P)

disp('Filtering...'); % @TODO: provide more filtering options
mnWav_filt = slopefilter(mnWav_raw, P.nDiff_filt);

disp('Spike detection...'); % provide more options @TODO

disp('Event grouping...'); 

disp('Feature extraction...'); % provide more options @TODO

end %func


%--------------------------------------------------------------------------
function S0 = cluster_(trFet, S0, P)
% cluster the data

end %func


%--------------------------------------------------------------------------
function P = load_params_(params)
% load params from a params file, fill default
% @TODO: read from a text file. params.paramfile = vcFile_prm (jrclust
% format)

P = call('file2struct', get_(params, 'paramfile'));

P = struct('fPaged_load', 1, 'MAX_BYTES_LOAD', 1e8, 'nPad_filt', 10, 'fSave_filt', 0, ...
    'vcFile', params.timeseries, 'nDiff_filt', str2double(params.nfilt), ...
    'params', params);
end %func


%--------------------------------------------------------------------------
function val = get_(S, vcField)
if isfield(S, vcField)
    val = S.(vcField);
else
    val = [];
end
end %func


%--------------------------------------------------------------------------
function S0 = ironclust_nopage_(params, P)
disp('Reading...');
mrWav_ = readmda(params.timeseries); % Later we will do things in chunks

[trFet, S0] = wav2fet_(mrWav_, P);
viClu = fet2clu_(trFet, S0, P);

disp('Writing...');
if P.fSave_filt, writemda32(Y, params.timeseries_out); end

end %func


%--------------------------------------------------------------------------
function spec = get_spec()

name='ironclust.ironclust';
version='0.1';
% 0.1: initial. mar 1, 2018

csInputList = {'timeseries', 'geom', 'script'};
csOutputList = {'timeseries_out', 'firings_out', 'metrics_out', 'templates_out'};
csParamList = {'paramfile', ...
    'samplerate', 'freq_min', 'frq_max', 'freq_wid', ...
    'detect_sign', 'detect_threshold', 'detect_interval', 'clip_size', ...
    'adjacency_radius', 'mask_out_artifacts', 'whiten', 'curate', ...
    'duration', 'waveform_upsamplefrac', ...
    'quantization_unit', 'consolidate_clusters', 'consolidation_factor', ...
    'subsample_factor', 'fit_stage', ...
    'waveforms_true', 'firings_true', ...
    '_params'...
    };

inputs = cellfun(@(vc)struct('name', vc), csInputList, 'UniformOutput', 0);
outputs = cellfun(@(vc)struct('name', vc), csOutputList, 'UniformOutput', 0);
parameters = cellfun(@(vc)struct('name', vc), csParamList, 'UniformOutput', 0);

% inputs={};
% inputs{end+1}=struct('name','timeseries');

% outputs={};
% outputs{end+1}=struct('name','firings');

% parameters={};
% parameters{end+1}=struct('name','paramfile');
% parameters{end+1}=struct('name','detect_sign'); %default: 0
% parameters{end+1}=struct('name','detect_threshold'); %default: 3
% parameters{end+1}=struct('name','freq_min'); %default: 300
% parameters{end+1}=struct('name','freq_max'); %default: 6000
% parameters{end+1}=struct('name','adjacency_radius'); %default: -1
% parameters{end+1}=struct('name','mask_out_artifacts'); %default: false
% parameters{end+1}=struct('name','whiten'); %default: true
% parameters{end+1}=struct('name','curate'); %default: false
% parameters{end+1}=struct('name','geom'); %default: false

opts=struct();

spec=struct( ...
    'name',name, ...
    'version',version, ...
    'inputs',{inputs}, ...
    'outputs',{outputs}, ...
    'parameters',{parameters}, ...
    'opts',opts ...
); 

end %func


%--------------------------------------------------------------------------
function [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(S_mda, P)
% plan load file size according to the available memory and file size (nBytes_file1)
% LOAD_FACTOR = 8; % Memory usage overhead

nSamples1 = S_mda.dimm(2);
if ~isfield(P, 'MAX_BYTES_LOAD'), P.MAX_BYTES_LOAD = []; end
if isempty(P.MAX_BYTES_LOAD), P.MAX_BYTES_LOAD = floor(mem_max_(P)); end

nSamples_max = floor(P.MAX_BYTES_LOAD / S_mda.dimm(1) / S_mda.nBytes_sample);

[nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples1, nSamples_max);
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
        nSamples_last1 = nSamples_last1 + nSamples_load1;
    end
end
end %func


%--------------------------------------------------------------------------
% 17/9/13 JJJ: Created and tested
function vr = setlim_(vr, lim_)
% Set low and high limits
vr = min(max(vr, lim_(1)), lim_(2));
end %func


%---------------------------------------------------------------------------
function out = ifeq_(if_, true_, false_)
if (if_)
    out = true_;
else
    out = false_;
end
end %func