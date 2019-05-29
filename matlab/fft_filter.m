%--------------------------------------------------------------------------
% 4/10/2019 JJJ: GPU memory use improved
% 3/11/2019 JJJ: GPU performance improved
% 10/15/2018 JJJ: Modified from ms_bandpass_filter (MountainLab) for
% memory-access efficiency
function mrWav_filt = fft_filter(mrWav, P, vcMode)
% [Usages]
% fft_filter()
%   clear cache
% mrWav_filt = fft_filter(mrWav, P)
% mrWav_filt = fft_filter(mrWav, P, vcMode)

if nargin<3, vcMode = ''; end
if isempty(vcMode), vcMode = 'bandpass'; end
fDebug = 0;

if ndims(mrWav) == 3
    for iSpk = 1:size(mrWav,3)
        mrWav(:,:,iSpk) = fft_filter(mrWav(:,:,iSpk), P, vcMode);
    end
    return;
end

if get_set_(P, 'fGpu', 1) && fDebug
    mrWav = gather(mrWav);
    P.fGpu = 0;
    fGpu_out = 1;
else
    fGpu_out = 0;
end

NSKIP_MAX = 2^19; % fft work length
scale_filter = get_set_(P, 'scale_filter', 1);
nPad = 300;
[nT, nC] = size(mrWav);
nSkip = min(nT, NSKIP_MAX);
[sRateHz, freqLim, freqLim_width, fGpu] = ...
    struct_get_(P, 'sRateHz', 'freqLim', 'freqLim_width', 'fGpu');
if isempty(freqLim), mrWav_filt = mrWav; return; end

try
    mrWav_filt = zeros(size(mrWav), 'like', mrWav); 
catch
    mrWav_filt = zeros(size(mrWav), class_(mrWav)); 
    fGpu = 0;
end
switch lower(vcMode)
    case 'bandpass', fh_make_filter = @fft_bandpass_;
    case 'highpass', fh_make_filter = @fft_highpass_;
    case 'lowpass', fh_make_filter = @fft_lowpass_;
    case {'fftdiff', 'banddiff'}, fh_make_filter = @fft_diff_;
    case 'wiener'
        vrFilt_ = fft_wiener_(mrWav, P); % build cache if not built yet
        fh_make_filter = @fft_wiener_;
    otherwise, error(['fftfilt_: invalid option: ', vcMode]);
end
if ~fGpu, mrWav = gather_(mrWav); end
fh_filter = @(x,f)real(ifft(bsxfun(@times, fft(single(x)), f)));
n_prev = nan;
fprintf('Running fft filter (%s)\n\t', vcMode); t1=tic;
if nT <= nSkip, nPad = 0; end
for iStart = 1:nSkip:nT
    iEnd = min(iStart+nSkip-1, nT);
    iStart1 = iStart - nPad;
    iEnd1 = iEnd + nPad;
    vi1 = iStart1:iEnd1;
    if iStart1 < 1 % wrap the filter (reflect boundary)
        vl_ = vi1 < 1;
        vi1(vl_) = 2 - vi1(vl_);
    end
    if iEnd1 > nT % wrap the filter (reflect boundary)
        vl_ = vi1 > nT;
        vi1(vl_) = 2*nT - vi1(vl_);
    end
    [mrWav1, fGpu] = gpuArray_(mrWav(vi1,:), fGpu);
    n1 = size(mrWav1,1);
    if n1 ~= n_prev
        vrFilt1 = fh_make_filter(n1, freqLim, freqLim_width, sRateHz);
        vrFilt1 = gpuArray_(vrFilt1, fGpu);
        if scale_filter ~= 1
            vrFilt1 = vrFilt1 * scale_filter;
        end
        n_prev = n1;
    end    
    if ~isempty(vrFilt1)        
        try
            mrWav1 = fh_filter(mrWav1, vrFilt1);  
        catch
            mrWav1 = fh_filter(gather_(mrWav1), vrFilt1);  
        end
    end
    if ~isGpu_(mrWav_filt)
        mrWav_filt(iStart:iEnd,:) = gather_(mrWav1(nPad+1:end-nPad,:));
    else
        mrWav_filt(iStart:iEnd,:) = mrWav1(nPad+1:end-nPad,:);
    end
    mrWav1 = []; % clear memory
    fprintf('.');
end
if ~isGpu_(mrWav), mrWav_filt = gather_(mrWav_filt); end
if fGpu_out
    mrWav_filt = gpuArray(mrWav_filt);
end
fprintf('\n\ttook %0.1fs (fGpu=%d)\n', toc(t1), fGpu);
end %func


%--------------------------------------------------------------------------
% 10/15/2018 JJJ: Modified from ms_bandpass_filter (MountainLab) 
function filt = fft_bandpass_(N, freqLim, freqLim_width, sRateHz)
% Usage
% [Y, filt] = bandpass_fft_(X, freqLim, freqLim_width, sRateHz)
% [filt] = 
% sRateHz: sampling rate
% freqLim: frequency limit, [f_lo, f_hi]
% freqLim_width: frequency transition width, [f_width_lo, f_width_hi]
if isempty(freqLim), freqLim = [nan, nan]; end
[flo, fhi] = deal(freqLim(1), freqLim(2));
[fwid_lo, fwid_hi] = deal(freqLim_width(1), freqLim_width(2));

[n1, n2, f] = get_freq_(N, sRateHz);

if ~isnan(flo) && ~isnan(fhi)
    filt = sqrt((1+erf((abs(f)-flo)/fwid_lo)) .* (1-erf((abs(f)-fhi)/fwid_hi)))/2;
elseif ~isnan(flo)
    filt = sqrt((1+erf((abs(f)-flo)/fwid_lo))/2);
elseif ~isnan(fhi)
    filt = sqrt((1-erf((abs(f)-fhi)/fwid_hi))/2);
else
    filt = [];
end
end %func


%--------------------------------------------------------------------------
function filt = fft_lowpass_(N, freqLim, freqLim_width, sRateHz)
[flo, fhi] = deal(freqLim(1), freqLim(2));
[fwid_lo, fwid_hi] = deal(freqLim_width(1), freqLim_width(2));
[n1, n2, f] = get_freq_(N, sRateHz);
filt = sqrt((1-erf((abs(f)-fhi)/fwid_hi))/2);
end %func


%--------------------------------------------------------------------------
function filt = fft_highpass_(N, freqLim, freqLim_width, sRateHz)
[flo, fhi] = deal(freqLim(1), freqLim(2));
[fwid_lo, fwid_hi] = deal(freqLim_width(1), freqLim_width(2));
[n1, n2, f] = get_freq_(N, sRateHz);
filt = sqrt((1+erf((abs(f)-flo)/fwid_lo))/2);
end %func


%--------------------------------------------------------------------------
% 4/12/2019 JJJ: band limited differentiator
function vrFilter_diff = fft_diff_(N, freqLim, freqLim_width, sRateHz)
% Usage
% [Y, filt] = bandpass_fft_(X, freqLim, freqLim_width, sRateHz)
% sRateHz: sampling rate
% freqLim: frequency limit, [f_lo, f_hi]
% freqLim_width: frequency transition width, [f_width_lo, f_width_hi]
if mod(N-1,2)==0
    n0 = (N-1)/2;
    vrFilter_diff = [0, linspace(0,2,n0), linspace(2,0,n0)]';
else
    n0 = (N-2)/2;
    vrFilter_diff = [0, linspace(0,2,n0), 2, linspace(2,0,n0)]';
end
if nargin==1, return; end

[flo, fhi] = deal(freqLim(1), freqLim(2));
[fwid_lo, fwid_hi] = deal(freqLim_width(1), freqLim_width(2));
if ~isnan(fhi)
    [n1, n2, f] = get_freq_(N, sRateHz);
    vrFilter_diff = vrFilter_diff .* sqrt(1-erf((abs(f)-fhi)/fwid_hi));
    vrFilter_diff = vrFilter_diff / max(vrFilter_diff) * 2;
end
end %funcmr1


%--------------------------------------------------------------------------
% 4/12/2019 JJJ: band limited differentiator
function vrFilter_wiener = fft_wiener_(N, freqLim, freqLim_width, sRateHz)
% Usage
% filt = fft_wiener_(mrWav, P)
%   set the persistent variable
% sRateHz: sampling rate
% freqLim: frequency limit, [f_lo, f_hi]
% freqLim_width: frequency transition width, [f_width_lo, f_width_hi]

persistent P_ filt0_
% computes the cache only if the cache is abscent
if nargin==2
    if ~isempty(filt0_)
        vrFilter_wiener = filt0_;
        return;
    end
    [mr_, P_] = deal(N, freqLim);
    mr1 = single(mr_(1:end/8, :)); % make a smaller copy
    mr1 = fft_filter(mr1, P_, 'bandpass');
    % reshape to the spike size and perform fft
    nSamples = (diff(P_.spkLim_raw) + 1)*2; %*2;
    mr2 = reshape_(mr1(:), nSamples);
    mrFft = fft(mr2-mean(mr2));
    switch 2
        case 4, vr2 = std(mr2(3:end,:) - mr2(1:end-2,:));
        case 3, vr2 = std(diff(mr2));
        case 2, vr2 = std(mr2);    
        case 1, vr2 = range(diff(mr2)); % find window containing spikes
    end
    threshLim = quantile(vr2, [.1, .99]);
    switch 2
        case 2, fh1 = @(x)mean(abs(mrFft(:,x)),2);
        case 1, fh1 = @(x)median(abs(mrFft(:,x)),2);
    end
    vrFft_noise = fh1(vr2<=threshLim(1));
    vrFft_noise_signal = fh1(vr2>=threshLim(2));
    vrFilt = rectify_(vrFft_noise_signal.^2 - vrFft_noise.^2) ./ vrFft_noise_signal.^2;
    vrFilt(vrFilt<0) = 0;
%     switch 3
%         case 3, scale = 1;
%         case 2, scale = quantile(vrFft_signal ./ vrFft_noise, 1/4);
%         case 1, scale = min(vrFft_signal(2:end)) / min(vrFft_noise(2:end));
%     end
% 
%     vrFft_signal = vrFft_signal(:) / scale - vrFft_noise(:);
%     vrFft_signal(vrFft_signal<0) = 0;
%     switch 1
%         case 2, vrFilt = vrFft_signal ./ (vrFft_signal + vrFft_noise);
%         case 1, vrFilt = vrFft_signal.^2 ./ (vrFft_signal.^2 + vrFft_noise.^2);
%     end
    vrFilt(1) = vrFilt(end);
    filt0_ = vrFilt;
    fprintf(2, 'fft_filter: wiener filter recomputed `filt0_`\n');
    vrFilter_wiener = vrFilt;
    return;
end


nSamples = numel(filt0_);
vrFilter_wiener = filt_symmetrize_(interp1(1:nSamples, filt0_, linspace(1,nSamples,N), 'linear'))';
vrFilter_wiener(1) = 0; % set zero dc gain

% apply high pass filter
switch 2
    case 3 %combine with diff
        vcFilter_diff = fft_diff_(N);
        vrFilter = vcFilter_diff .* vrFilter_wiener;
        vrFilter = vrFilter / max(vrFilter) * 2;
    case 2 %combine with High-pass
        vrFilter_bp = fft_bandpass_(N, P_.freqLim, P_.freqLim_width, P_.sRateHz);
        vrFilter = vrFilter_bp(:) .* vrFilter_wiener;
    case 1
        vrFilter = vrFilter_wiener;
end
filt_ = vrFilter;
end %func


%--------------------------------------------------------------------------
function x = rectify_(x)
x(x<0) = 0;
end %func


%--------------------------------------------------------------------------
% copied from irc.m on Apr 15, 2019

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
function flag = isGpu_(vr)
try
    flag = isa(vr, 'gpuArray');
catch
    flag = 0;
end
end


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
function mr = reshape_(vr, n1)
% n1: leading dimension
n = numel(vr);
n2 = floor(n / n1);
if n == (n1*n2)
    mr = reshape(vr, n1, n2);
else
    mr = reshape(vr(1:(n1*n2)), n1, n2);
end
end %func


%--------------------------------------------------------------------------
function [n1, n2, f] = get_freq_(N, sRateHz)
% [n1, n2, f] = get_freq_(N, sRateHz)
if mod(N,2)==0
    if nargin>=2
        df = sRateHz/N;
        f = df * [0:N/2, -N/2+1:-1]';
    end
    n1 = N/2+1;
else
    if nargin>=2
        df = sRateHz/N;
        f = df * [0:(N-1)/2, -(N-1)/2:-1]'; 
    end
    n1 = (N-1)/2+1;
end
n2 = N-n1;
end %func


%--------------------------------------------------------------------------
function vrFilt1 = filt_symmetrize_(vrFilt)
vrFilt1 = vrFilt;
% vrFilt1(1) = 0; % no DC gain
N = numel(vrFilt);
if mod((N-1),2) == 0 %even
    n1 = (N-1)/2;
else
    n1 = (N-2)/2;
    vrFilt1(1+n1+1) = (vrFilt(1+n1)+vrFilt(1+n1+2))/2; %set the middle to zero
end
vi1 = 2:(n1+1);
vi2 = N:-1:(N-n1+1);
vrFilt12 = (vrFilt(vi1) + vrFilt(vi2))/2;
vrFilt1(vi1) = vrFilt12;
vrFilt1(vi2) = vrFilt12;
% vrFilt1(1) = 0; % zero DC gain
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
