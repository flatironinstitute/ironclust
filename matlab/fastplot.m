function vh = fastplot(varargin)
% [Usages]
% fastplot(hAx1, vrX1, vrY1, vcLineStyle1, vrX2, vrY2, vcLineStyle2);
% fastplot(hAx1, vrX1, vrY1, vrX2, vrY2)
% fastplot(hAx1, vrX1, vrY1, vcLineStyle1)
% fastplot(hAx1, vrX1, vrY1)
% fastplot(hAx1, vrX1, vrY1)
% fastplot(hAx1, 'draw')
% fastplot(hAx1)    
%   same as fastplot(hAx1, 'draw')
% fastplot(hAx1, 'draw', xlim1)
% fastplot(hAx1, 'draw', xlim1, ylim1)
% fastplot(hAx1, 'draw', xlim1, ylim1)
% fastplot(var_name, var_value)
%   set `nPoints_max` (default is 10000)


persistent nPoints_max
if isempty(nPoints_max), nPoints_max = 10000; end

[hAx, vcLineStyle1, vcLineStyle2, vrX1, vrY1, vrX2, vrY2] = deal([]);

if nargin<1, return; end
if isa(varargin{1}, 'matlab.graphics.axis.Axes')
    hAx = varargin{1};
elseif ischar(varargin{1})
    val_ = varargin{2};
    switch varargin{1}
        case 'nPoints_max', nPoints_max = val_;
        otherwise, error('fastplot: invalid var_name `%s`', val_);
    end% switch
else
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
    [vrX1_, vrY1_, vrX2_, vrY2_, ~, scale, dimm] = subsample_(vrX1, vrY1, vrX2, vrY2, nPoints_max); 
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
        xlim0 = [min(vrX1_)-eps(), max(vrX1_)+eps()];
        ylim0 = [min(vrY1_)-eps(), max(vrY1_)+eps()];
    else
        xlim0 = [min(min(vrX1_), min(vrX2_))-eps(), max(max(vrX1_), max(vrX2_))+eps()];
        ylim0 = [min(min(vrY1_), min(vrY2_))-eps(), max(max(vrY1_), max(vrY2_))+eps()];
    end    
    S_line1 = struct('vrX', vrX1, 'vrY', vrY1, 'fastplot', 1, 'scale', scale, 'dimm', dimm);
    hLine1.UserData = S_line1;
    if ~isempty(hLine2)        
        hLine2.UserData = struct('vrX', vrX2, 'vrY', vrY2, 'fastplot', 1);
    end    
    fastplot = 1;
    S_ax = makeStruct_(hLine1, hLine2, xlim0, ylim0, fastplot, scale);
    set(hAx, 'UserData', S_ax, 'XLim', xlim0, 'YLim', ylim0);
else
    % draw
    S_ax = get(hAx, 'UserData');    
    [hLine1, hLine2, xlim0, ylim0, scale] = get_(S_ax, 'hLine1', 'hLine2', 'xlim0', 'ylim0', 'scale');
    
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
    [vrX1_, vrY1_, vrX2_, vrY2_, ~, scale, dimm] = subsample_(vrX1_, vrY1_, vrX2_, vrY2_, nPoints_max);    
    if ~isempty(hLine1)
        try
            set(hLine1, 'XData', vrX1_, 'YData', vrY1_);
            set_userdata_(hLine1, 'scale', scale, 'dimm', dimm);
        catch
            disp(lasterr());
        end
    end
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


%--------------------------------------------------------------------------
function [vrX, vrY, scale, dimm] = multiplot_(vrX, mrY, scale)
if numel(vrX) == numel(mrY), return; end
if nargin<3, scale = []; end
if isempty(scale)
    scale = single(max(abs(max(mrY(:))), abs(min(mrY(:)))));
end

nChans = size(mrY,2);
viY = 1:nChans;
mrX = repmat(vrX(:), 1, nChans);
mrX(end,:) = nan;
vrX = mrX(:);

mrY = bsxfun(@plus, single(mrY)/scale, viY(:)');
vrY = mrY(:);
dimm = size(mrY);
end %func


%--------------------------------------------------------------------------
function [viX1, vrY1, viX2, vrY2, vi, scale, dimm] = subsample_(viX1, vrY1, viX2, vrY2, nPoints_max, scale)
% vrY1: can be a matrix
if nargin<6, scale = []; end
dimm = [];

n = numel(viX1); % assert all points have same size
if n <= nPoints_max
    vi = 1:numel(viX1);
    if numel(viX1) ~= numel(vrY1)
        [viX1, vrY1, scale, dimm] = multiplot_(viX1, vrY1, scale); 
    end
else
    vi = round(linspace(1, n, nPoints_max));
    if ~isempty(vrY2)          
        [~,vi12] = intersect(viX1, viX2);
        vi = sort([vi(:); vi12(:)]);        
    end    
    if numel(viX1) == numel(vrY1)
        [viX1, vrY1] = deal(viX1(vi), vrY1(vi));
    else
        [viX1, vrY1, scale, dimm] = multiplot_(viX1(vi), vrY1(vi,:), scale); 
    end
end
end %func


%--------------------------------------------------------------------------
function [vrX1, vrY1, vrX2, vrY2] = filter_lim_(vrX1, vrY1, vrX2, vrY2, xlim)
if max(vrX1) > xlim(end) || min(vrX1) < xlim(1)
    vi1 = find(vrX1>=xlim(1) & vrX1 <= xlim(end));
    if numel(vrX1) == numel(vrY1)
        [vrX1, vrY1] = deal(vrX1(vi1), vrY1(vi1));
    else
        vrX1 = vrX1(vi1);
        vrY1 = vrY1(vi1,:);
    end
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