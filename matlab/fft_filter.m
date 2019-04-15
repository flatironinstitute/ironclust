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

% clear the cache
if nargout == 0
    clear fft_wiener_
    return;
end

NSKIP_MAX = 2^19; % fft work length
nPad = 300;
[nT, nC] = size(mrWav);
nSkip = min(nT, NSKIP_MAX);
[sRateHz, freqLim, freqLim_width, fGpu] = ...
    struct_get_(P, 'sRateHz', 'freqLim', 'freqLim_width', 'fGpu');
if isempty(freqLim), mrWav_filt = mrWav; return; end

try
    mrWav_filt = zeros(size(mrWav), 'like', mrWav); catch
    mrWav_filt = zeros(size(mrWav), class_(mrWav)); mrWav_filt
    fGpu = 0;
end
switch lower(vcMode)
    case 'bandpass', fh_filter = @fft_bandpass_;
    case {'fftdiff', 'banddiff'}, fh_filter = @fft_diff_;
    case 'wiener'
        fft_wiener_(mrWav, P); % set persistent variable
        fh_filter = @fft_wiener_;
    otherwise, error(['fftfilt_: invalid option: ', vcMode]);
end
if ~fGpu, mrWav = gather_(mrWav); end
fh_filt = @(x,f)real(ifft(bsxfun(@times, fft(single(x)), f)));
n_prev = nan;
fprintf('Running fft filter (%s)\n\t', vcMode); t1=tic;
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
        vrFilt1 = fh_filter(n1, freqLim, freqLim_width, sRateHz);
        vrFilt1 = gpuArray_(vrFilt1, fGpu);
        n_prev = n1;
    end    
    try
        mrWav1 = fh_filt(mrWav1, vrFilt1);  
    catch
        mrWav1 = fh_filt(gather_(mrWav1), vrFilt1);  
    end
    if ~isGpu_(mrWav_filt)
        mrWav_filt(iStart:iEnd,:) = gather_(mrWav1(nPad+1:end-nPad,:));
    else
        mrWav_filt(iStart:iEnd,:) = mrWav1(nPad+1:end-nPad,:);
    end
    mrWav1 = []; % clear memory
    fprintf('.');
end
clear fft_wiener_
if ~isGpu_(mrWav), mrWav_filt = gather_(mrWav_filt); end
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
[flo, fhi] = deal(freqLim(1), freqLim(2));
[fwid_lo, fwid_hi] = deal(freqLim_width(1), freqLim_width(2));

df = sRateHz/N;
if mod(N,2)==0
    f = df * [0:N/2, -N/2+1:-1]';
else
    f = df * [0:(N-1)/2, -(N-1)/2:-1]'; 
end
% if isa_(X, 'single'), f = single(f); end
% if isGpu_(X), f = gpuArray_(f); end

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
% 4/12/2019 JJJ: band limited differentiator
function filt = fft_diff_(N, freqLim, freqLim_width, sRateHz)
% Usage
% [Y, filt] = bandpass_fft_(X, freqLim, freqLim_width, sRateHz)
% sRateHz: sampling rate
% freqLim: frequency limit, [f_lo, f_hi]
% freqLim_width: frequency transition width, [f_width_lo, f_width_hi]
[flo, fhi] = deal(freqLim(1), freqLim(2));
[fwid_lo, fwid_hi] = deal(freqLim_width(1), freqLim_width(2));

df = sRateHz/N;
if mod(N,2)==0
    f = df * [0:N/2, -N/2+1:-1]';
    n1 = N/2+1;
else
    f = df * [0:(N-1)/2, -(N-1)/2:-1]'; 
    n1 = (N-1)/2+1;
end
n2 = numel(f)-n1;

% if isa_(X, 'single'), f = single(f); end
% if isGpu_(X), f = gpuArray_(f); end

filt = [linspace(0, 2, n1), linspace(2,0,n2)]';% [filt] = 

if ~isnan(fhi)
    filt = filt .* sqrt(1-erf((abs(f)-fhi)/fwid_hi));
    filt = filt / max(filt) * 2;
end
end %func


%--------------------------------------------------------------------------
% 4/12/2019 JJJ: band limited differentiator
function vrFilter = fft_wiener_(N, freqLim, freqLim_width, sRateHz)
% Usage
% filt = fft_wiener_(mrWav, P)
%   set the persistent variable
% sRateHz: sampling rate
% freqLim: frequency limit, [f_lo, f_hi]
% freqLim_width: frequency transition width, [f_width_lo, f_width_hi]

persistent mr_ P_ filt_
if nargin==2, [mr_, P_, v] = deal(N, freqLim, []); return; end

% use cache
if ~isempty(filt_), 
    if N ~= numel(filt_)
        filt_ = [];
    else
        vrFilter = filt_;
        return;
    end
end

% use a cache to 
% to do: do it channel by channel
mr1 = single(mr_(1:end/8, :)); % make a smaller copy
% reshape to the spike size and perform fft
nSamples = (diff(P_.spkLim_raw) + 1)*2;
mr2 = reshape_(mr1(:), nSamples);
mrFft = fft(mr2-mean(mr2));
% vrFft = abs(mean(mrFft,2));
% vrFft = mean(abs(mrFft),2);
% fit noise power
% 
switch 2
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
vrFft_signal = fh1(vr2>=threshLim(2));
scale = min(vrFft_signal(2:end)) / min(vrFft_noise(2:end));
vrFft_signal = abs(vrFft_signal(:) / scale - vrFft_noise(:));
switch 1
    case 2, vrFilt = vrFft_signal ./ (vrFft_signal + vrFft_noise);
    case 1, vrFilt = vrFft_signal.^2 ./ (vrFft_signal.^2 + vrFft_noise.^2);
end
vrFilt(1) = 0; % no DC gain
vrFilter = interp1(1:nSamples, vrFilt, linspace(1,nSamples,N), 'linear')';

% apply high pass filter
if 1    
    P_.freqLim(2) = nan;
    vrFilter_hp = fft_bandpass_(N, P_.freqLim, P_.freqLim_width, P_.sRateHz);
    vrFilter = vrFilter_hp(:) .* vrFilter;
end
if isempty(filt_), filt_ = vrFilter; end
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