%--------------------------------------------------------------------------
% current working condition
%  1. penetrating shank
%  2. one shank
% @TODO
%  1. multi-shank selection
%  2. 2D array display 
function show_drift_view(S0, hFig)
% [P, S_clu] = get_(S0, 'P', 'S_clu');
[iClu1, iClu2] = get_(S0, 'iCluCopy', 'iCluPaste');
[S_clu, mrPos_spk, viTime_spk, viSite_spk, P] = ...
    get_(S0, 'S_clu', 'mrPos_spk', 'viTime_spk', 'viSite_spk', 'P');
switch 2
    case 1
        [hAx_x, hAx_y] = get_userdata_(hFig, 'hAx_x', 'hAx_y');
        if isempty(hAx_y)
            hAx_y = axes_new_(hFig);
            set(hAx_y, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');  
            set_userdata_(hFig, hAx_y);
            ylabel(hAx_y, 'Y-position (um)'); xlabel(hAx_y, 'Time (s)'); grid(hAx_y, 'on');
        end
    case 2
        [hAx_x, hAx_y] = get_userdata_(hFig, 'hAx_x', 'hAx_y');
        if isempty(hAx_y) || isempty(hAx_x)
            hTabGroup = uitabgroup(hFig);
            htab_y = uitab(hTabGroup, 'Title', 'Y-position', 'BackgroundColor', 'w'); 
            htab_x = uitab(hTabGroup, 'Title', 'X-position', 'BackgroundColor', 'w'); 
            hAx_y = axes('Parent', htab_y, 'Position', [.05 .2 .9 .7], 'Tag', 'ax_drift_y');
            hAx_x = axes('Parent', htab_x, 'Position', [.05 .2 .9 .7], 'Tag', 'ax_drift_x');
            set_userdata_(hFig, hAx_x, hAx_y);
            ylabel(hAx_y, 'Y-position (um)'); xlabel(hAx_y, 'Time (s)'); grid(hAx_y, 'on');
            ylabel(hAx_x, 'X-position (um)'); xlabel(hAx_x, 'Time (s)'); grid(hAx_x, 'on');
        end
end %switch
if isempty(get(hFig, 'KeyPressFcn'))
    hFig.KeyPressFcn = @(h,e)keyPressFcn_(h,e);
end

%-----
% collect info

% clu1
[iShank1, iShank2] = deal(1);
viSpk1 = S_clu.cviSpk_clu{iClu1};
viTime1 = viTime_spk(viSpk1);
mrPos1 = mrPos_spk(viSpk1,:);
try iShank1 = P.viShank_site(S_clu.viSite_clu(iClu1)); catch, end

% clu2
if ~isempty(iClu2)
    viSpk2 = S_clu.cviSpk_clu{iClu2};
    viTime2 = viTime_spk(S_clu.cviSpk_clu{iClu2});
    mrPos2 = mrPos_spk(viSpk2,:);
    try iShank2 = P.viShank_site(S_clu.viSite_clu(iClu2)); catch, end
else
    [viSpk2, viTime2, mrPos2, iShank2] = deal([]);
end
if iShank1 ~= iShank2, iClu2 = []; end

% background spikes:
um_per_pix = get_set_(P, 'um_per_pix', 25);
ypos_lim1 = [min(mrPos1(:,2)), max(mrPos1(:,2))];
xpos_lim1 = [min(mrPos1(:,1)), max(mrPos1(:,1))];
vcTitle = '[H]elp; [M]erge; [S]plit; [B]ackground; [G]rid';
if ~isempty(iClu2)
    ypos_lim2 = [min(mrPos2(:,2)), max(mrPos2(:,2))];
    ypos_lim1 = [min(ypos_lim1(1), ypos_lim2(1)), max(ypos_lim1(2), ypos_lim2(2))];    
    xpos_lim2 = [min(mrPos2(:,1)), max(mrPos2(:,1))];
    xpos_lim1 = [min(xpos_lim1(1), xpos_lim2(1)), max(xpos_lim1(2), xpos_lim2(2))];    
    vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', iClu1, iClu2, vcTitle);
else
    vcTitle = sprintf('Clu%d (black); %s', iClu1, vcTitle);
end
ypos_lim1 = ypos_lim1 + um_per_pix * [-.5,.5];
viSpk0 = find(mrPos_spk(:,2) >= ypos_lim1(1) & mrPos_spk(:,2) <= ypos_lim1(end));
try
    viSpk0 = viSpk0(P.viShank_site(viSite_spk(viSpk0)) == iShank1);
catch
   	;
end
mrPos0 = mrPos_spk(viSpk0,:);
viTime0 = viTime_spk(viSpk0);


%-----
% Plot spk0, spk1, spk2
vhAx = [hAx_x, hAx_y];
cell_pos_lim1 = {xpos_lim1, ypos_lim1}; 
for iTab = 1:2
    hAx_ = vhAx(iTab);        
    [hPlot0, hPlot1, hPlot2] = get_userdata_(hAx_, 'hPlot0', 'hPlot1', 'hPlot2');
    % create line objects
    if isempty(hPlot0) || isempty(hPlot1) || isempty(hPlot2)
        mrC_ = get_set_(P, 'mrColor_proj', [.75 .75 .75; 0 0 0; 1 0 0]);
        hPlot0 = line(hAx_, nan, nan, 'Marker', '.', 'Color', mrC_(1,:), 'MarkerSize', 3, 'LineStyle', 'none');
        hPlot1 = line(hAx_, nan, nan, 'Marker', '.', 'Color', mrC_(2,:), 'MarkerSize', 5, 'LineStyle', 'none');
        hPlot2 = line(hAx_, nan, nan, 'Marker', '.', 'Color', mrC_(3,:), 'MarkerSize', 5, 'LineStyle', 'none');   %place holder  
        set_userdata_(hAx_, hPlot0, hPlot1, hPlot2);
    end
    set(hPlot0, 'XData', single(viTime0)/P.sRateHz, 'YData', mrPos0(:,iTab));
    set(hPlot1, 'XData', single(viTime1)/P.sRateHz, 'YData', mrPos1(:,iTab));
    if isempty(iClu2)
        set(hPlot2, 'XData', nan, 'YData', nan);
    else
        set(hPlot2, 'XData', single(viTime2)/P.sRateHz, 'YData', mrPos2(:,iTab));
    end
    pos_lim1 = cell_pos_lim1{iTab};
    axes_lim0 = [0, viTime_spk(end)/P.sRateHz, pos_lim1(1), pos_lim1(2)];
    axis(hAx_, axes_lim0);
    set_userdata_(hAx_, axes_lim0);
    title_(hAx_, vcTitle);
end
drawnow_();

% add callback functions
csHelp = {...
    '[R]eset scale', ...
    'Zoom: mouse wheel', ...
    'Pan: Drag while pressing wheel', ...
    '[G]rid on/off'};
set_userdata_(hFig, csHelp);
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
function varargout = get_userdata_(varargin)
% fDelete: delete the field after retrieving
if nargin<1, return; end
try
    S = get(varargin{1}, 'UserData');
    for iArg = 2:nargin        
        vcName1 = varargin{iArg};
        if isfield(S, vcName1)
            varargout{iArg-1} = S.(vcName1);
        else
            varargout{iArg-1} = [];
        end
    end
catch
    return;
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
    for iArg_in = 2:nargin
        [vcName1, val1] = deal(inputname(iArg_in), varargin{iArg_in});
        S_userdata.(vcName1) = val1;
    end %for
else
    return;
end
set(h, 'UserData', S_userdata);
end %func


%--------------------------------------------------------------------------
function keyPressFcn_(hFig, event, S0)
if nargin<3, S0 = get(0, 'UserData'); end
[P, S_clu] = deal(S0.P, S0.S_clu);
uiTab_ = hFig.Children(1).SelectedTab;
hAx = uiTab_.Children(1);

switch lower(event.Key)
    case 'r' %reset view
        axis(hAx, get_userdata_(hAx, 'axes_lim0'));

    case 'm' %merge
        ui_merge_(S0);
        
    case 'h' %help
        msgbox_(get_userdata_(hFig, 'csHelp'), 1);
        
    case 'b' %background spike toggle
        toggleVisible_(get_userdata_(hAx, 'hPlot0'));
        
    case 'g' %grid toggle
        fGrid = get_userdata_(hAx, 'fGrid');
        if isempty(fGrid), fGrid = true; end
        fGrid = ~fGrid;
        grid_(hAx, fGrid);
        set_userdata_(hAx, fGrid);
        
    case 's' %split. draw a polygon
        if ~isempty(S0.iCluPaste)
            msgbox_('Select one cluster'); return;
        end
        try
            hPoly = impoly_();
            if isempty(hPoly); return; end
            mrPolyPos = getPosition(hPoly);
            hPlot1 = get_userdata_(hAx, 'hPlot1');
            vrX1 = double(get(hPlot1, 'XData'));
            vrY1 = double(get(hPlot1, 'YData'));
            vlIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
            hSplit = line(vrX1(vlIn), vrY1(vlIn), 'Color', [1 0 0], 'Marker', '.', 'LineStyle', 'none');
            if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes'), 'yes')
                split_clu_(S0.iCluCopy, vlIn);
            end
            delete_multi_(hPoly, hSplit);
        catch
            disp(lasterror());
        end
end %switch
end %func


%==========================================================================
% call irc.m
function xylabel_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function drawnow_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function ui_merge_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function toggleVisible_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function grid_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function split_clu_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function delete_multi_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end

function out1 = axes_new_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = title_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = msgbox_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = impoly_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = questdlg_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
