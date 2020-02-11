%--------------------------------------------------------------------------
function irc2_traces(vcFile_prm)
% show raw traces
% If file format is nChans x nSamples, load subset of file (fTranspose=1)
% If file format is nSamples x nChans, load all and save to global (fTranspose=0)
% 2017/6/22 James Jun: Added multiview (nTime_traces )
global fDebug_ui mnWav mnWav1 % only use if P.fTranspose=0


S0 = load0_(vcFile_prm);
P = S0.P;
set(0, 'UserData', S0);
set0_(fDebug_ui);

% get file to show
iFile_show = 1; %files to display for clustered together
if ~isempty(get_(P, 'csFile_merge'))
    csFiles_bin = filter_files_(P.csFile_merge);
    if numel(csFiles_bin)==1
        vcFile_bin = csFiles_bin{1}; 
    else %show multiple files        
        if isempty(vcFileId)
            arrayfun(@(i)fprintf('%d: %s\n', i, csFiles_bin{i}), 1:numel(csFiles_bin), 'UniformOutput', 0);
            fprintf('---------------------------------------------\n');
            vcFileId = input('Please specify Fild ID from the list above:', 's');
        end
        if isempty(vcFileId), return; end
        iFile_show = str2num(vcFileId);        
        try
            vcFile_bin = csFiles_bin{iFile_show};
        catch
            return;
        end
    end
else
    vcFile_bin = P.vcFile; % if multiple files exist, load first
end
set0_(iFile_show);
tlim_bin = P.tlim;

% Open file
fprintf('Opening %s\n', vcFile_bin);
[fid_bin, nBytes_bin, header_offset] = fopen_(vcFile_bin, 'r');
if P.header_offset ~= header_offset && header_offset > 0
    P.header_offset = header_offset;
    set0_(P); % update header_offset for .ns5 format
end
if isempty(fid_bin), fprintf(2, '.bin file does not exist: %s\n', vcFile_bin); return; end
nSamples_bin = floor(nBytes_bin / bytesPerSample_(P.vcDataType) / P.nChans);
nLoad_bin = min(round(diff(tlim_bin) * P.sRateHz), nSamples_bin);
if tlim_bin(1)>0
    iSample_bin = ceil(tlim_bin(1) * P.sRateHz) + 1; %offset sample number    
else
    iSample_bin = 1; %sample start location
end    
nlim_bin = [0,nLoad_bin-1] + iSample_bin;
if nlim_bin(1) < 1, nlim_bin = [1, nLoad_bin]; end
if nlim_bin(2) > nSamples_bin, nlim_bin = [-nLoad_bin+1, 0] + nSamples_bin; end

nTime_traces = get_(P, 'nTime_traces');
[cvn_lim_bin, viRange_bin] = sample_skip_(nlim_bin, nSamples_bin, nTime_traces);    
if P.fTranspose_bin   
    mnWav = [];
    fseek_(fid_bin, iSample_bin, P);
    if nTime_traces > 1            
        mnWav1 = load_bin_multi_(fid_bin, cvn_lim_bin, P)';        
    else        
        mnWav1 = load_bin_(fid_bin, P.vcDataType, [P.nChans, nLoad_bin])'; %next keypress: update tlim_show    
    end
%     @TODO: load from cvn_lim_bin specifiers. check for end or beginning when keyboard command
else %load whole thing
    mnWav = load_bin_(fid_bin, P.vcDataType, [nSamples_bin, P.nChans]); %next keypress: update tlim_show
    fclose(fid_bin);
    fid_bin = []; 
    %mnWav1 = mnWav((nlim_bin(1):nlim_bin(2)), :);
    mnWav1 = mnWav(viRange_bin, :);
    disp('Entire raw traces are cached to RAM since fTranspose=0.');
end %if
mnWav1 = uint2int_(mnWav1);

% full screen width
hFig_traces = create_figure_('Fig_traces', [0 0 .5 1], vcFile_bin, 0, 1); %remove all other figure traces
hAx = axes_new_(hFig_traces); % create axis
hPlot = line(hAx, nan, nan, 'Color', [1 1 1]*.5, 'Parent', hAx, 'LineWidth', .5);
hPlot_edges = plot_(nan, nan, 'Color', [1 0 0]*.5, 'Parent', hAx, 'LineWidth', 1);
set(hAx, 'Position',[.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual');
S_fig = makeStruct_(hAx, hPlot, nlim_bin, fid_bin, nSamples_bin, nLoad_bin, hPlot_edges);
S_fig.maxAmp = P.maxAmp;
S_fig.vcTitle = '[H]elp; (Sft)[Up/Down]:Scale(%0.1f uV); (Sft)[Left/Right]:Time; [F]ilter; [J]ump T; [C]han. query; [R]eset view; [P]SD; [S]pike; [A]ux chan; [E]xport; [T]race; [G]rid';
S_fig.csHelp = { ...    
    'Left/Right: change time (Shift: x4)', ...
    '[J]ump T', ...
    '[Home/End]: go to beginning/end of file', ...
    '---------', ...
    'Up/Down: change scale (Shift: x4)', ...
    'Zoom: Mouse wheel', ...
    '[x/y/ESC]: zoom direction', ...
    'Pan: hold down the wheel and drag', ...
    '[R]eset view', ...
    '---------', ...
    '[F]ilter toggle', ...
    '[S]pike toggle', ...
    'Gri[D] toggle', ...
    '[T]races toggle', ... 
    '---------', ...    
    '[C]hannel query', ...
    '[A]ux channel display', ...
    '[P]ower spectrum', ...    
    '[E]xport to workspace', ...
    };
S_fig = struct_merge_(S_fig, ...
    struct('vcGrid', 'on', 'vcFilter', 'off', 'vcSpikes', 'on', 'vcTraces', 'on'));
set(hFig_traces, 'UserData', S_fig);
set(hFig_traces, 'color', 'w', 'KeyPressFcn', @keyPressFcn_Fig_traces_, 'BusyAction', 'cancel', 'CloseRequestFcn', @close_hFig_traces_);
mouse_figure(hFig_traces);

Fig_traces_plot_(1); % Plot spikes and color clusters 
end %func


%--------------------------------------------------------------------------
% plot data
function Fig_traces_plot_(fAxis_reset)
% fAxis_reset: reset the axis limit
% [usage]
% Fig_traces_plot_()
% Fig_traces_plot_(fAxis_reset)
% 2017/06/22 James Jun
% 6/22 JJJ: added seperator lines, fixed the reset view and spike view

global mnWav1 mrWav1 % current timeslice to plot
if nargin<1, fAxis_reset = 0; end
fWait = msgbox_('Plotting...',0,1);
fShuttleOrder = 1; %shuffle cluster color
S0 = get(0, 'UserData'); 
[P, S_auto] = get_(S0, 'P', 'S_auto');
[hFig, S_fig] = get_fig_cache_('Fig_traces'); 
figure_wait_(1, hFig);
sRateHz = P.sRateHz / P.nSkip_show;
viSamples1 = 1:P.nSkip_show:size(mnWav1,1);
spkLim = round(P.spkLim / P.nSkip_show); %show 2x of range
P.vcFilter = get_filter_(P);
if strcmpi(S_fig.vcFilter, 'on')
    P1=P; P1.sRateHz = sRateHz; P1.fGpu = 0;
    P1.vcFilter = get_set_(P, 'vcFilter_show', P.vcFilter);
%     if P.fft_thresh>0, mnWav1 = fft_clean_(mnWav1, P); end
    mrWav1 = bit2uV_(filter_car_(mnWav1(viSamples1, P.viSite2Chan), P1), P1);
    vcFilter_show = P1.vcFilter;
else
    mrWav1 = meanSubt_(single(mnWav1(viSamples1, P.viSite2Chan))) * P.uV_per_bit;    
    vcFilter_show = 'off';
end
viSites = 1:numel(P.viSite2Chan);
% mrWav1 = meanSubt_(single(mnWav1(:, P.viSite2Chan))) * P.uV_per_bit;
% hide bad channels
nTime_traces = get_(P, 'nTime_traces');
if isempty(nTime_traces) || nTime_traces==1
    vrTime_bin = ((S_fig.nlim_bin(1):P.nSkip_show:S_fig.nlim_bin(end))-1) / P.sRateHz;    
    vcXLabel = 'Time (s)';
else    
    vrTime_bin = (0:(size(mrWav1,1)-1)) / (P.sRateHz / P.nSkip_show) + (S_fig.nlim_bin(1)-1) / P.sRateHz;
    [cvn_lim_bin, viRange_bin, viEdges] = sample_skip_(S_fig.nlim_bin, S_fig.nSamples_bin, nTime_traces);
    tlim_show = (cellfun(@(x)x(1), cvn_lim_bin([1,end]))) / P.sRateHz;
    vcXLabel = sprintf('Time (s), %d segments merged (%0.1f ~ %0.1f s, %0.2f s each)', nTime_traces, tlim_show, diff(P.tlim));
    mrX_edges = vrTime_bin(repmat(viEdges(:)', [3,1]));
    mrY_edges = repmat([0;numel(P.viSite2Chan)+1;nan],1,numel(viEdges));
    set(S_fig.hPlot_edges, 'XData', mrX_edges(:), 'YData', mrY_edges(:));
    csTime_bin = cellfun(@(x)sprintf('%0.1f', x(1)/P.sRateHz), cvn_lim_bin, 'UniformOutput', 0);
    set(S_fig.hAx, {'XTick', 'XTickLabel'}, {vrTime_bin(viEdges), csTime_bin});
end
multiplot(S_fig.hPlot, S_fig.maxAmp, vrTime_bin, mrWav1, viSites);
% axis(S_fig.hAx, [vrTime_bin(1), vrTime_bin(end), viSites(1)-1, viSites(end)+1]);
grid(S_fig.hAx, S_fig.vcGrid);
set(S_fig.hAx, 'YTick', viSites);
title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp));
xlabel(S_fig.hAx, vcXLabel); 
ylabel(S_fig.hAx, 'Site #'); 
set(S_fig.hPlot, 'Visible', S_fig.vcTraces);

% Delete spikes from other threads
S_fig_ = get(hFig, 'UserData');
if isfield(S_fig_, 'chSpk'), delete_multi_(S_fig_.chSpk); end
if isfield(S_fig, 'chSpk'), delete_multi_(S_fig.chSpk); end

% plot spikes
if strcmpi(S_fig.vcSpikes, 'on') && isfield(S0, 'viTime_spk')
%     viTime_spk = int32(S0.viTime_spk) - int32(S0.viT_offset_file(S0.iFile_show));
    viTime_spk = (S0.viTime_spk);
    if nTime_traces > 1
        viSpk1 = find(in_range_(viTime_spk, cvn_lim_bin));
        [viSite_spk1, viTime_spk1] = multifun_(@(vr)vr(viSpk1), S0.viSite_spk, viTime_spk);
        viTime_spk1 = round(reverse_lookup_(viTime_spk1, viRange_bin) / P.nSkip_show);
    else
        viSpk1 = find(viTime_spk >= S_fig.nlim_bin(1) & viTime_spk < S_fig.nlim_bin(end));
        [viSite_spk1, viTime_spk1] = multifun_(@(vr)vr(viSpk1), S0.viSite_spk, viTime_spk);
        viTime_spk1 = round((viTime_spk1 - S_fig.nlim_bin(1) + 1) / P.nSkip_show); %time offset
    end        
    t_start1 = single(S_fig.nlim_bin(1) - 1) / P.sRateHz;
    viSite_spk1 = single(viSite_spk1);
    % check if clustered
    if isempty(S_auto)
        nSites = size(mrWav1,2);
        chSpk = cell(nSites, 1);
        for iSite=1:nSites %deal with subsample factor
            viSpk11 = find(viSite_spk1 == iSite);
            if isempty(viSpk11), continue; end
            viTime_spk11 = viTime_spk1(viSpk11);
            [mrY11, mrX11] = vr2mr_(mrWav1(:,iSite), viTime_spk11, spkLim); %display purpose x2
%             vr2mr_spk_(mrWav1(:,iSite), viTime_spk11, P);
            mrT11 = single(mrX11-1) / sRateHz + t_start1;
            chSpk{iSite} = line(nan, nan, 'Color', [1 0 0], 'LineWidth', 1.5, 'Parent', S_fig.hAx);
            multiplot(chSpk{iSite}, S_fig.maxAmp, mrT11, mrY11, iSite);
        end        
    else % different color for each clu
        viClu_spk1 = S_auto.viClu(viSpk1);        
        mrColor_clu = [jet(S_auto.nClu); 0 0 0];        
        vrLineWidth_clu = (mod((1:S_auto.nClu)-1, 3)+1)'/2 + .5;  %(randi(3, S_clu.nClu, 1)+1)/2;
        if fShuttleOrder
            mrColor_clu = shuffle_static_(mrColor_clu, 1);
            vrLineWidth_clu = shuffle_static_(vrLineWidth_clu, 1);
        end
        nSpk1 = numel(viTime_spk1);
        chSpk = cell(nSpk1, 1);
        for iSpk1 = 1:nSpk1
            iTime_spk11 = viTime_spk1(iSpk1);
            iSite11 = viSite_spk1(iSpk1);
            [mrY11, mrX11] = vr2mr_(mrWav1(:,iSite11), iTime_spk11, spkLim); %display purpose x2
            mrT11 = double(mrX11-1) / sRateHz + t_start1;
            iClu11 = viClu_spk1(iSpk1);
            if iClu11<=0, continue; end
%                 vrColor1 = [0 0 0]; 
%                 lineWidth1 = .5;
%             else
            vrColor1 = mrColor_clu(iClu11,:);
            lineWidth1 = vrLineWidth_clu(iClu11);
%             end
            chSpk{iSpk1} = line(nan, nan, 'Color', vrColor1, 'LineWidth', lineWidth1, 'Parent', S_fig.hAx);
            multiplot(chSpk{iSpk1}, S_fig.maxAmp, mrT11, mrY11, iSite11);
        end
    end
    S_fig.chSpk = chSpk;
else
    % delete spikes    
    S_fig.chSpk = [];    
end
if fAxis_reset, fig_traces_reset_(S_fig); end
set(hFig, 'UserData', S_fig, 'Name', sprintf('%s: filter: %s', P.vcFile_prm, (vcFilter_show)));
figure_wait_(0, hFig);
close_(fWait);
end %func


%--------------------------------------------------------------------------
function [mr, miRange] = vr2mr_(vr, vi, spkLim)
% JJJ 2015 Dec 24
% vr2mr2: quick version and doesn't kill index out of range
% assumes vi is within range and tolerates spkLim part of being outside
% works for any datatype

% prepare indices
if size(vi,2)==1, vi=vi'; end %row
viSpk = int32(spkLim(1):spkLim(end))';
miRange = viSpk + int32(vi);
miRange(miRange<1) = 1; 
miRange(miRange > numel(vr)) = numel(vr); %keep # sites consistent
% miRange = int32(miRange);

% build spike table
nSpks = numel(vi);
mr = reshape(vr(miRange(:)), [], nSpks);
end %func


%--------------------------------------------------------------------------
function fig_traces_reset_(S_fig)
global mnWav1

if nargin<1, [hFig, S_fig] = get_fig_cache_('Fig_traces');  end
% axis(S_fig.hAx, [S_fig.nlim_bin / P.sRateHz, 0, nSites+1]);
S0=get(0,'UserData'); P=S0.P;
nTime_traces = get_(P, 'nTime_traces');
if nTime_traces > 1
    tlim1 = ([0, size(mnWav1,1)] + S_fig.nlim_bin(1) - 1) / P.sRateHz;
    tlim1 = round(tlim1*1000)/1000;
    axis_(S_fig.hAx, [tlim1, 0, numel(P.viSite2Chan)+1]);
else
    axis_(S_fig.hAx, [S_fig.nlim_bin / P.sRateHz, 0, numel(P.viSite2Chan)+1]);
end 
end %func


%--------------------------------------------------------------------------
function keyPressFcn_Fig_traces_(hFig, event)
% 2017/6/22 James Jun: Added nTime_traces multiview

global mnWav1 mrWav1 mnWav
S0 = get(0, 'UserData'); 
P = S0.P;
S_fig = get(hFig, 'UserData');
factor = 1 + 3 * key_modifier_(event, 'shift');
nSites = numel(P.viSite2Chan);

switch lower(event.Key)
    case 'h', msgbox_(S_fig.csHelp, 1);
        
    case {'uparrow', 'downarrow'}
        if isfield(S_fig, 'chSpk')
            S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, S_fig.hPlot, S_fig.chSpk);
        else
            S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, S_fig.hPlot);
        end
        title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp));        
        set(hFig, 'UserData', S_fig);
        
    case {'leftarrow', 'rightarrow', 'j', 'home', 'end'}
        switch lower(event.Key)
            case 'leftarrow'
                nlim_bin = S_fig.nlim_bin - (S_fig.nLoad_bin) * factor; %no overlap
                if nlim_bin(1)<1
                    msgbox_('Beginning of file', 1); 
                    nlim_bin = [1, S_fig.nLoad_bin]; 
                end
            case 'rightarrow'
                nlim_bin = S_fig.nlim_bin + (S_fig.nLoad_bin + 1) * factor; %no overlap
                if nlim_bin(2) > S_fig.nSamples_bin
                    msgbox_('End of file', 1); 
                    nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
                end
            case 'home' %beginning of file
                nlim_bin = [1, S_fig.nLoad_bin];
            case 'end' %end of file
                nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
            case 'j'
                vcAns = inputdlg_('Go to time (s)', 'Jump to time', 1, {'0'});
                if isempty(vcAns), return; end
                try
                    nlim_bin = round(str2double(vcAns)*P.sRateHz) + [1, S_fig.nLoad_bin];
                catch
                    return;
                end
        end %switch
        nTime_traces = get_(P, 'nTime_traces');
        [cvn_lim_bin, viRange_bin] = sample_skip_(nlim_bin, S_fig.nSamples_bin, nTime_traces); 
        if P.fTranspose_bin
            fseek_(S_fig.fid_bin, nlim_bin(1), P);
            if nTime_traces > 1
                mnWav1 = load_bin_multi_(S_fig.fid_bin, cvn_lim_bin, P)';
            else
                mnWav1 = load_bin_(S_fig.fid_bin, P.vcDataType, [P.nChans, S_fig.nLoad_bin])';            
            end
        else
            mnWav1 = mnWav(viRange_bin, :);
        end
        mnWav1 = uint2int_(mnWav1);
        S_fig.nlim_bin = nlim_bin;
        set_fig_(hFig, S_fig);
        Fig_traces_plot_(1); %redraw
        
    case 'f' %apply filter
        S_fig.vcFilter = str_toggle_(S_fig.vcFilter, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();        
        
    case 'g' %grid toggle on/off
        S_fig.vcGrid = str_toggle_(S_fig.vcGrid, 'on', 'off');
        grid(S_fig.hAx, S_fig.vcGrid);
        set(hFig, 'UserData', S_fig);
        
    case 'r' %reset view
        fig_traces_reset_(S_fig);        
        
    case 'e' %export current view        
        assignWorkspace_(mnWav1, mrWav1);
        disp('mnWav1: raw traces, mrWav1: filtered traces');
        
    case 'p' %power spectrum
        iSite_show = inputdlg_num_(sprintf('Site# to show (1-%d, 0 for all)', nSites), 'Site#', 0);
        if isnan(iSite_show), return; end        
        hFig = create_figure_('FigPsd', [.5 0 .5 1], P.vcFile_prm, 1, 1); %show to the right
        % ask user which channels to plot
        if iSite_show>0
            mrWav2 = mrWav1(:, iSite_show);
        else
            mrWav2 = mrWav1;
        end
        plotMedPower_(mrWav2, 'sRateHz', P.sRateHz/P.nSkip_show, 'viChanExcl', P.viSiteZero);
        
    case 's' %show/hide spikes
        S_fig.vcSpikes = str_toggle_(S_fig.vcSpikes, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();
        
    case 't' %show/hide traces
        S_fig.vcTraces = str_toggle_(S_fig.vcTraces, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();         
        
    case 'c' %channel query
        msgbox_('Draw a rectangle', 1);
        hRect = imrect_();
        if isempty(hRect), return; end
        vrPos_rect = getPosition(hRect);            
        S_plot = get(S_fig.hPlot, 'UserData');
        vrX = get(S_fig.hPlot, 'XData');
        vrY = get(S_fig.hPlot, 'YData');
        viIndex = find(vrX >= vrPos_rect(1) & vrX <= sum(vrPos_rect([1,3])) & vrY >= vrPos_rect(2) & vrY <= sum(vrPos_rect([2,4])));
        if isempty(viIndex), delete_multi_(hRect); return; end
        index_plot = round(median(viIndex));
        [time1, iSite] = ind2sub(size(mrWav1), index_plot);        
        mrX = reshape(vrX, S_plot.dimm);
        mrY = reshape(vrY, S_plot.dimm);
        hold(S_fig.hAx, 'on');
        hPoint = plot_(vrX(index_plot), vrY(index_plot), 'r*');
        hLine = plot_(S_fig.hAx, mrX(:,iSite), mrY(:,iSite), 'r-');
        hold(S_fig.hAx, 'off');
        iChan = P.viSite2Chan(iSite);
        msgbox_(sprintf('Site: %d/ Chan: %d', iSite, iChan), 1);
        delete_multi_(hRect, hLine, hPoint);
end %return if S_fig didn't change
end %func


%--------------------------------------------------------------------------
function close_hFig_traces_(hFig, event)
try
    if ~ishandle(hFig), return; end
    if ~isvalid(hFig), return; end
    S_fig = get(hFig, 'UserData');    
    fclose_(S_fig.fid_bin);
    try delete(hFig); catch; end %close one more time
catch
    disperr_();
    close(hFig);
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function S0 = set0_(varargin)
% Set(0, 'UserData')

S0 = get(0, 'UserData'); 
% set(0, 'UserData', []); %prevent memory copy operation
for i=1:nargin
    try
        S0.(inputname(i)) = varargin{i};
    catch
        disperr_();
    end
end
set(0, 'UserData', S0);
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
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
end %func


%--------------------------------------------------------------------------
% Call from irc2.m
function cout = call_irc2_(dbstack1, cin, nargout)
vcFunc = dbstack1(1).name;
try
    switch nargout
        case 0, cout{1} = []; irc2('call', vcFunc, cin);
        case 1, cout{1} = irc2('call', vcFunc, cin);
        case 2, [cout{1}, cout{2}] = irc2('call', vcFunc, cin);
        case 3, [cout{1}, cout{2}, cout{3}] = irc2('call', vcFunc, cin);
        case 4, [cout{1}, cout{2}, cout{3}, cout{4}] = irc2('call', vcFunc, cin);
        case 5, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}] = irc2('call', vcFunc, cin);
        case 6, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}, cout{6}] = irc2('call', vcFunc, cin);
        case 7, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}, cout{6}, cout{7}] = irc2('call', vcFunc, cin);
        case 8, [cout{1}, cout{2}, cout{3}, cout{4}, cout{5}, cout{6}, cout{7}, cout{8}] = irc2('call', vcFunc, cin);
        otherwise, error('call_irc2_: too many output');
    end
catch ME
    fprintf(2, 'call_irc2_: %s\n', ME.message);
    rethrow ME;
end
end %func


function varargout = load_fet_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = axes_new_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = bytesPerSample_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = create_figure_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = disperr_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = filter_files_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = fopen_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = fseek_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_bin_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_bin_multi_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load0_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = sample_skip_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_merge_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = uint2int_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = bit2uV_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = close_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = delete_multi_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = figure_wait_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = filter_car_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_filter_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = in_range_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = meanSubt_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = msgbox_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = multifun_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = reverse_lookup_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = shuffle_static_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = title_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = vr2mr3_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = axis_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = get0_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = assignWorkspace_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = change_amp_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = imrect_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = inputdlg_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = inputdlg_num_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = key_modifier_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plotMedPower_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = set_fig_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = str_toggle_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = fclose_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_fig_cache_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
