function vh = fastplot(varargin)
% fastplot(hAx1, vrX1, vrY1, vcLineStyle1, vrX2, vrY2, vcLineStyle2);
% fastplot(hAx1, vrX1, vrY1, vrX2, vrY2)
% fastplot(hAx1, vrX1, vrY1, vcLineStyle1)
% fastplot(hAx1, vrX1, vrY1)
% fastplot(hAx1, vrX1, vrY1)

% fastplot(hAx1, 'draw')
% fastplot(hAx1)

% fastplot(hAx1, 'draw', xlim1)
% fastplot(hAx1, 'draw', xlim1, ylim1)

nPoints_max = 10000;

[hAx, vcLineStyle1, vcLineStyle2, vrX1, vrY1, vrX2, vrY2] = deal([]);

if nargin<1, return; end
hAx = varargin{1};
if ~isa(hAx, 'matlab.graphics.axis.Axes')
    error('fastplot: first argument should be an axes handle');
end

% command parsing
fDraw = 0;
if nargin==1
    fDraw = 1;
else
    if ischar(varargin{2})
        if strcmp(varargin{2}, 'draw')
            fDraw = 1; 
        else
            error('fastplot: invalid command.');
        end
    end
end
if ~fDraw
    switch nargin
        case 2, vrY1 = varargin{2}; vrX1 = (1:numel(vrY1))';
        case 3, vrX1 = varargin{2}; vrY1 = varargin{3};
        case 4, vrX1 = varargin{2}; vrY1 = varargin{3}; vcLineStyle1 = varargin{4};
        case 5
            vrX1 = varargin{2}; vrY1 = varargin{3}; 
            vrX2 = varargin{4}; vrY2 = varargin{5}; 
        case 7
            vrX1 = varargin{2}; vrY1 = varargin{3}; vcLineStyle1 = varargin{4};
            vrX2 = varargin{5}; vrY2 = varargin{6}; vcLineStyle2 = varargin{7};
        otherwise, error('fastplot: unsupported number of arguments');
    end %switch    
    
    % draw first time
    [vrX1_, vrY1_, vrX2_, vrY2_] = subsample_(vrX1, vrY1, vrX2, vrY2, nPoints_max);    
    if isempty(vrY2)
        if isempty(vcLineStyle1)
            hLine1 = plot(hAx, vrX1_, vrY1_);
        else
            hLine1 = plot(hAx, vrX1_, vrY1_, vcLineStyle1);
        end
        hLine2 = [];
    else
        if ~isempty(vcLineStyle1)
            vhPlot = plot(hAx, vrX1_, vrY1_, vcLineStyle1, vrX2_, vrY2_, vcLineStyle2);            
        else
            vhPlot = plot(hAx, vrX1_, vrY1_, vrX2_, vrY2_);
        end
        [hLine1, hLine2] = deal(vhPlot(1), vhPlot(2));
    end
    if isempty(vrX2)
        xlim0 = [min(vrX1)-eps(), max(vrX1)+eps()];
        ylim0 = [min(vrY1)-eps(), max(vrY1)+eps()];
    else
        xlim0 = [min(min(vrX1), min(vrX2))-eps(), max(max(vrX1), max(vrX2))+eps()];
        ylim0 = [min(min(vrY1), min(vrY2))-eps(), max(max(vrY1), max(vrY2))+eps()];
    end
    hLine1.UserData = struct('vrX', vrX1, 'vrY', vrY1, 'fastplot', 1);
    if ~isempty(hLine2)        
        hLine2.UserData = struct('vrX', vrX2, 'vrY', vrY2, 'fastplot', 1);
    end    
    fastplot = 1;
    S_ax = makeStruct_(hLine1, hLine2, xlim0, ylim0, fastplot);
    set(hAx, 'UserData', S_ax, 'XLim', xlim0, 'YLim', ylim0);
else
    % draw
    S_ax = get(hAx, 'UserData');    
    [hLine1, hLine2, xlim0, ylim0] = get_(S_ax, 'hLine1', 'hLine2', 'xlim0', 'ylim0');
    
    % filter axes
    switch nargin
        case {1, 2}, xlim1 = xlim0; ylim1 = ylim0;
        case 3, xlim1 = varargin{3}; ylim1 = [];
        case 4, xlim1 = varargin{3}; ylim1 = varargin{4};
        otherwise, error('fastplot: unsupported number of arguments');
    end
    [vrX1, vrY1] = get_(hLine1.UserData, 'vrX', 'vrY');
    if ~isempty(hLine2)
        [vrX2, vrY2] = get_(hLine2.UserData, 'vrX', 'vrY');
    else
        [vrX2, vrY2] = deal([]);
    end
    [vrX1_, vrY1_, vrX2_, vrY2_] = filter_lim_(vrX1, vrY1, vrX2, vrY2, xlim1);    
    [vrX1_, vrY1_, vrX2_, vrY2_] = subsample_(vrX1_, vrY1_, vrX2_, vrY2_, nPoints_max);    
    if ~isempty(hLine1), set(hLine1, 'XData', vrX1_, 'YData', vrY1_); end
    if ~isempty(hLine2), set(hLine2, 'XData', vrX2_, 'YData', vrY2_); end
    set(hAx, 'XLim', xlim1);
    if ~isempty(ylim1), set(hAx, 'YLim', ylim1); end
end

if isempty(hLine2)
    vh = hLine1;
else
    vh = [hLine1, hLine2];
end    
end %func


%--------------------------------------------------------------------------
function [vrX1, vrY1, vrX2, vrY2, vi] = subsample_(vrX1, vrY1, vrX2, vrY2, nPoints_max)
n = numel(vrX1); % assert all points have same size
if n <= nPoints_max
    return;
else
    vi = round(linspace(1, n, nPoints_max));
    if ~isempty(vrY2)     
        vi = sort([vi(:); vrX2(vrX2 <= n)]);
    end
    vrX1 = vrX1(vi);
    vrY1 = vrY1(vi);
end
end %func


%--------------------------------------------------------------------------
function [vrX1, vrY1, vrX2, vrY2] = filter_lim_(vrX1, vrY1, vrX2, vrY2, xlim)
if max(vrX1) > xlim(end) || min(vrX1) < xlim(1)
    vi1 = find(vrX1>=xlim(1) & vrX1 <= xlim(end));
    vrX1 = vrX1(vi1);
    vrY1 = vrY1(vi1);
end
if ~isempty(vrY2)
    if max(vrX2) > xlim(end) || min(vrX2) < xlim(1)
        vi2 = find(vrX2>=xlim(1) & vrX2 <= xlim(end));
        vrX2 = vrX2(vi2);
        vrY2 = vrY2(vi2);
    end
end
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