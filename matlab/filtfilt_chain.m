function mr = filtfilt_chain(mr, varargin)

P = funcDefStr_(funcInStr_(varargin{:}), ...
    'filtOrder', 3, 'sRateHz', 25000, 'fEllip', 1, ...
    'freqLimStop', [], 'freqLim', [], 'freqLimNotch', [], ...
    'nPad_filt', 100, 'fGpu_filt', 0, 'gain_boost', 1, 'nDiff_filt', 0);
if isfield(P, 'fGpu')
    P.fGpu_filt = P.fGpu;
end


[cvrA, cvrB] = deal({});
if ~isempty(P.freqLim)
   [cvrB{end+1}, cvrA{end+1}] = makeFilt_(P.freqLim, 'bandpass', P);
end
if ~isempty(P.freqLimStop)
   [cvrB{end+1}, cvrA{end+1}] = makeFilt_(P.freqLimStop, 'stop', P);
end  
if ~isempty(P.freqLimNotch)
    if ~iscell(P.freqLimNotch)
        [cvrB{end+1}, cvrA{end+1}] = makeFilt_(P.freqLimNotch, 'notch', P);
    else
        for iCell=1:numel(P.freqLimNotch)
            [cvrB{end+1}, cvrA{end+1}] = makeFilt_(P.freqLimNotch{iCell}, 'notch', P);
        end
    end
end    
    
%----------------
% Run the filter chain
fInt16 = isa(mr, 'int16');
if P.fGpu_filt, mr = gpuArray(mr); end
mr = filt_pad_('add', mr, P.nPad_filt); %slow
if fInt16, mr = single(mr); end   %double for long data?

% first pass
for iFilt=1:numel(cvrA)
%     mr = filtfilt_single(cvrB{iFilt}, cvrA{iFilt}, mr);
%     mr = filter(cvrB{iFilt}, cvrA{iFilt}, mr);
%     if P.fFir_filt, continue; end %digital filter
%     mr = flipud(filter(cvrB{iFilt}, cvrA{iFilt}, flipud(mr)));
    mr = flipud(filter(cvrB{iFilt}, cvrA{iFilt}, ...
            flipud(filter(cvrB{iFilt}, cvrA{iFilt}, mr))));
end
% mr = flipud(mr);

% second pass, backward
% for iFilt=numel(cvrA):-1:1
%     mr = filter(cvrB{iFilt}, cvrA{iFilt}, mr);
% end
% mr = flipud(mr);
if P.gain_boost ~= 1
    mr = mr * P.gain_boost;
end
if fInt16, mr = int16(mr); end
mr = filt_pad_('remove', mr, P.nPad_filt); %slow    
%if P.fGpu_filt, mr = gather(mr); end
end %func


%--------------------------------------------------------------------------
function [vrFiltB, vrFiltA] = makeFilt_(freqLim, vcType, P)
if nargin<2, vcType = 'bandpass'; end
freqLim = freqLim / P.sRateHz * 2;
if ~strcmpi(vcType, 'notch')
    if P.fEllip  %copied from wave_clus
        if isinf(freqLim(1)) || freqLim(1) <= 0
            [vrFiltB, vrFiltA]=ellip(P.filtOrder,0.1,40, freqLim(2), 'low');
        elseif isinf(freqLim(2))
            [vrFiltB, vrFiltA]=ellip(P.filtOrder,0.1,40, freqLim(1), 'high');
        else
            [vrFiltB, vrFiltA]=ellip(P.filtOrder,0.1,40, freqLim, vcType);
        end    
    else
        if isinf(freqLim(1)) || freqLim(1) <= 0
            [vrFiltB, vrFiltA] = butter(P.filtOrder, freqLim(2),'low');        
        elseif isinf(freqLim(2))
            [vrFiltB, vrFiltA] = butter(P.filtOrder, freqLim(1),'high');
        else
            [vrFiltB, vrFiltA] = butter(P.filtOrder, freqLim, vcType);    
        end
    end
else
    [vrFiltB, vrFiltA] = iirnotch_(mean(freqLim), diff(freqLim));
end
% vrFiltA = single(vrFiltA);
% vrFiltB = single(vrFiltB);  
if P.fGpu_filt  
    vrFiltB = gpuArray(vrFiltB);
    vrFiltA = gpuArray(vrFiltA);
end
end %func


%--------------------------------------------------------------------------
function mrWav = filt_pad_(vcMode, mrWav, nPad)
% add padding by reflection. removes artifact at the end

if isempty(nPad), return; end
if nPad==0, return; end
nPad = min(nPad, size(mrWav,1));

switch lower(vcMode)
    case 'add'
        mrWav = [flipud(mrWav(1:nPad,:)); mrWav; flipud(mrWav(end-nPad+1:end,:))];
    case 'remove'
        mrWav = mrWav(nPad+1:end-nPad,:);
end %switch
end %func
    

%--------------------------------------------------------------------------
function [num,den] = iirnotch_(Wo,BW)
    % Define default values.
Ab = abs(10*log10(.5)); % 3-dB width
% Design a 2nd-order notch digital filter.

% Inputs are normalized by pi.
BW = BW*pi;
Wo = Wo*pi;

Gb   = 10^(-Ab/20);
beta = (sqrt(1-Gb.^2)/Gb)*tan(BW/2);
gain = 1/(1+beta);

num  = gain*[1 -2*cos(Wo) 1];
den  = [1 -2*gain*cos(Wo) (2*gain-1)];
end %func


%--------------------------------------------------------------------------
function P = funcDefStr_(P, varargin)
csNames = varargin(1:2:end);
csValues = varargin(2:2:end);

for iField=1:numel(csNames)
    if ~isfield(P, csNames{iField})
        P = setfield(P, csNames{iField}, csValues{iField});
    end
end
end %func


%--------------------------------------------------------------------------
function P = funcInStr_(varargin)

if isempty(varargin), P=struct(); return; end
if isstruct(varargin{1}), P = varargin{1}; return; end

csNames = varargin(1:2:end);
csValues = varargin(2:2:end);
P = struct();
for iField=1:numel(csNames)
    if ~isfield(P, csNames{iField})
%         v1 = csValues{iField};
%         eval(sprintf('P.%s = v1;', csNames{iField}));
        P = setfield(P, csNames{iField}, csValues{iField});
    end
end
end %func



