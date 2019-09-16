%--------------------------------------------------------------------------
% load from the file and construct waveforms
% phase 1: identify spike ID to load
% phase 2: paged loading and accumulation
% phase 3: waveform simiarity comparison

function mrDist_clu = clu_wave_similarity_paged(S_clu, P)

[KNN, MAX_SAMPLE, NUM_PC, fUse_raw] = deal(16, 4000, 3, 0);
nAve_knn = min(KNN, get_set_(P, 'knn', 30));
MIN_COUNT = get_set_(P, 'min_count');

fprintf('Automated merging (post-hoc)\n'); t_func = tic;
[viClu_spk, miKnn, vrRho_spk] = deal(S_clu.viClu, S_clu.miKnn, S_clu.rho);
miKnn = miKnn(1:nAve_knn,:);
[viSite_spk, viTime_spk] = get0_('viSite_spk', 'viTime_spk');
nShift_max = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
viShift = -nShift_max:nShift_max;
nClu = S_clu.nClu;
MAX_BYTES_LOAD = get_set_(P, 'MAX_BYTES_LOAD', 1e9);


% phase 1. identify spikes to load per drift per cluster
fprintf('\tFinding spike indices '); t1=tic;
cviSite_drift_clu = cell(nClu, 1);
ccviSpk_drift_clu = cell(nClu, 1);
nDrift = get_set_(P, 'nTime_drift', 64);
nSpk_min = get_set_(P, 'knn', 30);
for iClu = 1:nClu
    viSpk1 = find(viClu_spk == iClu);
    viSpk1 = viSpk1(:);
    viiSpk1 = round(linspace(1, numel(viSpk1), nDrift+1));
    [vlKeep_clu1, viSite_clu1] = deal(true(nDrift, 1), zeros(nDrift,1));
    cviSpk_clu1 = cell(nDrift, 1);
    [miKnn1, vrRho1] = deal(miKnn(:,viSpk1), vrRho_spk(viSpk1)');
    for iDrift = 1:nDrift
        vii1 = viiSpk1(iDrift):viiSpk1(iDrift+1);
        [vrRho11, miKnn11] = deal(vrRho1(vii1), miKnn1(:,vii1));
        
        viSpk11 = miKnn11(vrRho_spk(miKnn11) >= vrRho11);
        viSpk11 = unique(viSpk11);
        iSite11 = mode(viSite_spk(viSpk11));
        vl_ = viSite_spk(viSpk11) == iSite11;
        viSpk11 = viSpk11(vl_);   
        viSite_clu1(iDrift) = iSite11;
        if numel(viSpk11) < nSpk_min
            vlKeep_clu1(iDrift) = false;
            continue;
        end        
        cviSpk_clu1{iDrift} = subsample_vr_(viSpk11, MAX_SAMPLE);
    end
    ccviSpk_drift_clu{iClu} = cviSpk_clu1(vlKeep_clu1);
    cviSite_drift_clu{iClu} = viSite_clu1(vlKeep_clu1);
    fprintf('.');
end
% turn ccviSpk_clu into matrix: nSpk x 3 dimension
miSpk_drift_clu = cell2mat_vi_(cellfun(@(x)cell2mat_vi_(x), ccviSpk_drift_clu, 'UniformOutput', 0));
miSpk_drift_clu = sort_col_(miSpk_drift_clu, 1); % sort by spike index
vnDrift_clu = cellfun(@numel, ccviSpk_drift_clu);
cvnSpk_drift_clu = cellfun_(@(x)cellfun(@numel, x), ccviSpk_drift_clu);
fprintf('took %0.1fs\n', toc(t1));



% phase 2. load spikes from file (memory paging)
vcFile_spkwav = strrep(P.vcFile_prm, '.prm', '_spkwav.jrc');
[dimm_spk, type_spk] = get0_('dimm_spk', 'type_spk');
load_spkwav_('open', vcFile_spkwav, dimm_spk, type_spk);
ctrWav_clu = cell(nClu, 1);
iiSpk_end = 0;
[nLoads, nSpk_load, nSpk_last] = load_spkwav_('plan', MAX_BYTES_LOAD); % give buffer
for iLoad = 1:nLoads
    if iLoad == nLoads
        nSpk_load1 = nSpk_last;
    else
        nSpk_load1 = nSpk_load;
    end
    [trWav_spk1, viSpk1] = load_spkwav_('load', nSpk_load1);
    
    iiSpk_start = iiSpk_end + 1;
    iiSpk_end = find(miSpk_drift_clu(:,1) <= viSpk1(end), 1, 'last');
    miSpk_drift_clu1 = miSpk_drift_clu(iiSpk_start:iiSpk_end, :);
    
    % Accumulate waveforms per cluster per drift
    for iClu = 1:nClu
        trWav11 = ctrWav_clu{iClu};
        nDrift11 = vnDrift_clu(iClu);
        if isempty(trWav11)
            trWav11 = zeros(dimm_spk(1), dimm_spk(2), nDrift11); 
        end        
        vii11 = find(miSpk_drift_clu1(:,3) == iClu);
        viiSpk11 = miSpk_drift_clu1(vii11,1) - viSpk1(1) + 1;
        viDrift11 = miSpk_drift_clu1(vii11,2);
        cviiSpk11 = arrayfun(@(x)viiSpk11(viDrift11==x), 1:nDrift11, 'UniformOutput', 0);        
        cmr_ = arrayfun(@(i)trWav11(:,:,i) + sum(trWav_spk1(:,:,cviiSpk11{i}), 3), ...
                1:nDrift11, 'UniformOutput', 0);
        ctrWav_clu{iClu} = cat(3, cmr_{:});
    end 
end %for
load_spkwav_('close');


% phase 3. compute waveform similarity
fh_norm = @(x)bsxfun(@rdivide, x, std(x,1)*sqrt(size(x,1)));
fh_norm_tr = @(x)fh_norm(reshape(x, [], size(x,3)));
mrDist_clu = nan(nClu, 'single');
for iClu1 = 1:nClu
    viSite_clu1 = cviSite_drift_clu{iClu1};
    if isempty(viSite_clu1), continue; end
    mr1_ = fh_norm_tr(shift_trWav_(ctrWav_clu{iClu1}, viShift));
    viSite_clu1 = repmat(viSite_clu1(:), numel(viShift), 1);
    vrDist_clu1 = zeros(nClu, 1, 'single');
    for iClu2 = iClu1+1:nClu
        viSite2 = cviSite_drift_clu{iClu2};
        viSite12 = intersect(viSite_clu1, viSite2);
        if isempty(viSite12), continue; end
        mr2_ = fh_norm_tr(ctrWav_clu{iClu2});   
        for iSite12_ = 1:numel(viSite12)
            iSite12 = viSite12(iSite12_);
            mrDist12 = mr2_(:, viSite2==iSite12)' * mr1_(:, viSite_clu1==iSite12);
            vrDist_clu1(iClu2) = max(vrDist_clu1(iClu2), max(mrDist12(:)));
        end
    end
    mrDist_clu(:, iClu1) = vrDist_clu1;
%     fprintf('.');
end %for
end %func


%--------------------------------------------------------------------------
% cell2mat and add a column of cell index at the end
% concatenate in the row direction (dimension 1)
function mr1 = cell2mat_vi_(cmr1)
% input
% ----
% cmr1: cell (c x 1) of 2D matrix (n_k x m, k=1..c)
% mr1: 2D matrix ( n_sum x (m+1)), n_sum = sum(n_k, k=1..c)

vn1 = cellfun(@(x)size(x,1), cmr1);
arrayfun_ = @(fh,A)cell2mat(arrayfun(fh, A(:), 'UniformOutput', 0));
vi1 = arrayfun_(@(x)repmat(x, vn1(x), 1), 1:numel(cmr1));
mr1 = [cell2mat(cmr1(:)), vi1];
end %func


%--------------------------------------------------------------------------
function B = sort_col_(A, iCol, vcDirection)
if nargin<3, vcDirection = 'ascend'; end

[~, ix] = sort(A(:, iCol), vcDirection);
B = A(ix, :);
end %func


%--------------------------------------------------------------------------
% load binary file
function varargout = load_spkwav_(varargin)
% syntax
% fid = load_spkwav_('open', vcFile, dimm, type)
% nLoad = load_spkwav_('plan', nBytes_buffer)
% [mr_loaded, viSpk] = load_spkwav_('load', nSamples_requested)
% fSuccess = load_spkwav_('close')

persistent fid vcFile dimm vcDataType nBytes_file iSpk_start

vcMode = varargin{1};
switch lower(vcMode)
    case 'open'
        [vcFile, dimm, vcDataType] = deal(varargin{2}, varargin{3}, varargin{4});
        try
            [fid, nBytes_file] = fopen_(vcFile, 'r');
        catch
            fprintf(2, 'File open error: %s\n', lasterr()); 
            fid = [];
            nBytes_file = 0;
        end
        varargout{1} = fid;
        varargout{2} = nBytes_file;
        iSpk_start = 1;
        
    case 'load'
        dimm1 = dimm;
        if nargin >=2
            dimm1(end) = varargin{2}; 
        end
        try
            assert(~isempty(fid), 'file is invalid');
            mr1 = fread(fid, prod(dimm1), ['*', vcDataType]);
            mr1 = reshape(mr1, dimm1);
        catch
            fprintf(2, 'File load: %s\n', lasterr()); 
            mr1 = [];
        end
        varargout{1} = mr1;
        varargout{2} = [iSpk_start:(iSpk_start+dimm1(end)-1)]; % tell number of samples loaded
        iSpk_start = iSpk_start + size(mr1, 1);
        
    case 'close'
        try
            fclose(fid);
            fSuccess = 1;
        catch
            fSuccess = 0;
            fprintf(2, 'File close error: %s\n', lasterr()); 
        end
        [fid, vcFile] = deal([]);
        varargout{1} = fSuccess;
        
    case 'plan'
        if nargin>=2
            nBytes_buffer = varargin{2};
        else
            nBytes_buffer = 1e9;
        end
        if numel(dimm)  == 1
            nBytes_per_spk = bytesPerSample_(vcDataType);
        else
            nBytes_per_spk = prod(dimm(1:end-1)) * bytesPerSample_(vcDataType);
        end
        nSamples_file = floor(nBytes_file / nBytes_per_spk);
        nSamples_load1 = floor(nBytes_buffer / nBytes_per_spk);
        nByte_load1 = nSamples_load1 * nBytes_per_spk;
        nLoads = ceil(nBytes_file / nByte_load1);
        if nLoads > 1
            nSamples_last1 = nSamples_file - (nLoads-1) * nSamples_load1;
        else
            nSamples_last1 = nSamples_file;
            nSamples_load1 = nSamples_last1;
        end
        [varargout{1}, varargout{2}, varargout{3}] = deal(nLoads, nSamples_load1, nSamples_last1);
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function varargout = get0_(varargin)
% returns get(0, 'UserData') to the workspace
% [S0, P] = get0_();
% [S0, P, S_clu] = get0_();
% [var1, var2] = get0_('var1', 'var2'); %sets [] if doesn't exist
S0 = get(0, 'UserData'); 
if ~isfield(S0, 'S_clu'), S0.S_clu = []; end
if nargin==0
    varargout{1} = S0; 
    if nargout==0, assignWorkspace_(S0); return; end
    if nargout>=1, varargout{1} = S0; end
    if nargout>=2, varargout{2} = S0.P; end
    if nargout>=3, varargout{3} = S0.S_clu; end
    return;
end
for i=1:nargin
    try                
        eval(sprintf('%s = S0.%s;', varargin{i}, varargin{i}));
        varargout{i} = S0.(varargin{i});
    catch
        varargout{i} = [];
    end
end
end %func


%==========================================================================
% call irc.m
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = bytesPerSample_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cellfun_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsample_vr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = shift_trWav_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function [out1, out2] = fopen_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
