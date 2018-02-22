function ret = p_slopefilter(params) % params already passing string
% question: exception handling. if file doesn't exist what to do?

% internal parameters
P = struct('fPaged_load', 1, 'MAX_BYTES_LOAD', 1e8, 'nPad_filt', 10);

if (ischar(params))&&(strcmp(params,'spec'))
    ret=get_spec();
    return;
else
    ret = [];
end

P.nfilt = str2double(params.nfilt);
t1 = tic;

% Open file 
[S_mda_r, fid_r] = readmda_header(params.timeseries);
nChans = S_mda_r.dimm(1); % specify number of channels
P.nChans = nChans;

if ~P.fPaged_load
    nLoad = 1;
else
    [nLoad, nSamples_load, nSamples_last] = plan_load_(S_mda_r, P);
end

if nLoad == 1
    slopefilter_nopage_(params, P);
else
    for iLoad = 1:nLoad
        % Read        
        fprintf('Processing %d/%d...\n', iLoad, nLoad);
        nSamples_ = ifeq_(iLoad == nLoad, nSamples_last, nSamples_load);
        if iLoad == 1
            viT_ = 1:nSamples_;
            nSamples_load_ = nSamples_ + P.nPad_filt;            
        elseif iLoad == nLoad
            viT_ = (1+P.nPad_filt):(nSamples_+P.nPad_filt);
            nSamples_load_ = nSamples_ + P.nPad_filt;
        else
            viT_ = (1+P.nPad_filt):(nSamples_+P.nPad_filt);
            nSamples_load_ = nSamples_ + 2 * P.nPad_filt;
        end

%         fprintf('\tLoading from file...'); t_load_ = tic();
        mnWav_ = fread(fid_r, [nChans, nSamples_load_], S_mda_r.vcDataType);
%         fprintf('took %0.1fs\n', toc(t_load_));
        fseek(fid_r, -2*nChans*P.nPad_filt*S_mda_r.nBytes_sample, 'cof');

        % Filter
        mnWav_ = slopefilter(mnWav_, P.nfilt);

        % Write
        if iLoad == 1
            S_mda_w = setfield(S_mda_r, 'vcDataType', class(mnWav_));
            fid_w = writemda_header(params.timeseries_out, S_mda_w);    
        end
        fwrite(fid_w, mnWav_(:,viT_), S_mda_w.vcDataType);        
        mnWav_ = []; % clear memory    
        if iLoad == nLoad, fclose(fid_w); end
    end %for
end %if 
fclose(fid_r);

fprintf('Done with slopefilter, took %0.1fs (nLoad=%d)\n', ...
    toc(t1), nLoad);

end % func


%--------------------------------------------------------------------------
function slopefilter_nopage_(params, P)
disp('Reading...');
X = readmda(params.timeseries); % Later we will do things in chunks

disp('Filtering...');
Y = slopefilter(X, P.nfilt);

disp('Writing...');
writemda32(Y, params.timeseries_out);
end %func


%--------------------------------------------------------------------------
function spec = get_spec()

name='ironclust.slopefilter';
version='0.11';
% 0.1: single memory loading
% 0.11: paged memory loading

inputs={};
inputs{end+1}=struct('name','timeseries');

outputs={};
outputs{end+1}=struct('name','timeseries_out');

parameters={};
parameters{end+1}=struct('name','nfilt');

opts=struct;

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