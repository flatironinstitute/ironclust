% current working condition
%  1. penetrating shank
%  2. one shank
% @TODO
%  1. multi-shank selection
%  2. 2D array display 
function show_drift_view(S0, hFig)
% [P, S_clu] = get_(S0, 'P', 'S_clu');
[iClu1, iClu2] = get_(S0, 'iCluCopy', 'iCluPaste');
[S_clu, mrPos_spk, viTime_spk, P] = ...
    get_(S0, 'S_clu', 'mrPos_spk', 'viTime_spk', 'P');
S_fig = get(hFig, 'UserData');
switch 2
    case 2
        hAx_y = get_(S_fig, 'hAx_y');
        if isempty(hAx_y)
            hAx_y = axes_new_(hFig);
            set(hAx_y, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');            
            S_fig.hAx_y = hAx_y;
            set(hFig, 'UserData', S_fig);
            ylabel(hAx_y, 'Y-position (um)'); xlabel(hAx_y, 'Time (s)'); grid(hAx_y, 'on');
        end
    case 1
        [hAx_x, hAx_y] = get_(S_fig, 'hAx_x', 'hAx_y');
        % create a new axes if it doesn't exist
        if isempty(hAx_y)
            hTabGroup = uitabgroup(hFig);
            htab_y = uitab(hTabGroup, 'Title', 'Y-position', 'BackgroundColor', 'w'); 
            htab_x = uitab(hTabGroup, 'Title', 'X-position', 'BackgroundColor', 'w'); 
            hAx_y = axes('Parent', htab_y, 'OuterPosition', [0 0 1 1]); 
            hAx_x = axes('Parent', htab_x, 'OuterPosition', [0 0 1 1]);  
            S_fig.hAx_y = hAx_y;
            S_fig.hAx_x = hAx_x;
            set(hFig, 'UserData', S_fig);
        %     set(hAx_y, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');
            ylabel(hAx_y, 'Y-position (um)'); xlabel(hAx_y, 'Time (s)'); grid(hAx_y, 'on');
            ylabel(hAx_x, 'X-position (um)'); xlabel(hAx_x, 'Time (s)'); grid(hAx_x, 'on');
        end
end %switch

%-----
% collect info

% clu1
[iShank1, iShank2] = deal(1);
viSpk1 = S_clu.cviSpk_clu{iClu1};
viTime1 = viTime_spk(viSpk1);
vrPos1 = mrPos_spk(viSpk1,2);
% try iShank1 = P.viShank_site(S_clu.viSite_clu(iClu1)); catch, end

% clu2
if ~isempty(iClu2)
    viSpk2 = S_clu.cviSpk_clu{iClu2};
    viTime2 = viTime_spk(S_clu.cviSpk_clu{iClu2});
    vrPos2 = mrPos_spk(viSpk2,2);
    try iShank2 = P.viShank_site(S_clu.viSite_clu(iClu2)); catch, end
else
    [viSpk2, viTime2, vrPos2, iShank2] = deal([]);
end
if iShank1 ~= iShank2, iClu2 = []; end

% background spikes:
um_per_pix = get_set_(P, 'um_per_pix', 25);
ypos_lim1 = [min(vrPos1), max(vrPos1)] + um_per_pix * [-1,1];
viSpk0 = find(mrPos_spk(:,2) >= ypos_lim1(1) & mrPos_spk(:,2) <= ypos_lim1(end));
vrPos0 = mrPos_spk(viSpk0, 2);
viTime0 = viTime_spk(viSpk0);


%-----
% Plot spk0, spk1, spk2
S_ax = get(hAx_y, 'UserData');
[hPlot0, hPlot1, hPlot2] = deal(get_(S_ax, 'hPlot0'), get_(S_ax, 'hPlot1'), get_(S_ax, 'hPlot2'));
% create line objects
if isempty(hPlot0) || isempty(hPlot1) || isempty(hPlot2)
    mrC_ = get_set_(P, 'mrColor_proj', [.75 .75 .75; 0 0 0; 1 0 0]);
    hPlot0 = line(hAx_y, nan, nan, 'Marker', '.', 'Color', mrC_(1,:), 'MarkerSize', 5, 'LineStyle', 'none');
    hPlot1 = line(hAx_y, nan, nan, 'Marker', '.', 'Color', mrC_(2,:), 'MarkerSize', 5, 'LineStyle', 'none');
    hPlot2 = line(hAx_y, nan, nan, 'Marker', '.', 'Color', mrC_(3,:), 'MarkerSize', 5, 'LineStyle', 'none');   %place holder  
    S_ax.hPlot0 = hPlot0;
    S_ax.hPlot1 = hPlot1;
    S_ax.hPlot2 = hPlot2;
    set(hAx_y, 'UserData', S_ax);
end

set(hPlot0, 'XData', single(viTime0)/P.sRateHz, 'YData', vrPos0);
set(hPlot1, 'XData', single(viTime1)/P.sRateHz, 'YData', vrPos1);
if isempty(iClu2)
    set(hPlot2, 'XData', nan, 'YData', nan);
else
    set(hPlot2, 'XData', single(viTime2)/P.sRateHz, 'YData', vrPos2);
end
axis(hAx_y, [0, viTime_spk(end)/P.sRateHz, ypos_lim1(1), ypos_lim1(2)]);
drawnow_();
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
function val = get_userdata_(h, vcName, fDelete)
% fDelete: delete the field after retrieving
val = [];
if nargin<3, fDelete = 0; end
try
    S_userdata = get(h, 'UserData');
    if ~isfield(S_userdata, vcName), return; end
    val = S_userdata.(vcName);
    if fDelete
        set(h, 'UserData', rmfield(S_userdata, vcName));
    end
catch
    val = [];
end
end %func


%--------------------------------------------------------------------------
function S_userdata = set_userdata_(varargin)
% Usage
% -----
% S_userdata = set_userdata_(h, val)
% S_userdata = set_userdata_(h, name1, val1, name2, val2, ...)
h = varargin{1};
S_userdata = get(h, 'UserData');
if nargin==2
    [val, vcName] = deal(varargin{2}, inputname(2));
    try
        S_userdata.(vcName) = val;
    catch
        disp(lasterr());
    end
elseif nargin>2
    
    for iArg_in = 2:2:nargin
        [vcName1, val1] = deal(varargin{iArg_in}, varargin{iArg_in+1});
        S_userdata.(vcName1) = val1;
    end %for
else
    return;
end
set(h, 'UserData', S_userdata);
end %func


%==========================================================================
% call irc.m
function xylabel_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function drawnow_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function out1 = title_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = axes_new_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
