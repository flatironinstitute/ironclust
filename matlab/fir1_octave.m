% Copied from octave "signal" package
% Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
% Code converted using oct2ml.m

% Examples:
% freqz (fir1 (40, 0.3));
% freqz (fir1 (15, [0.2, 0.5], "stop"));  # note the zero-crossing at 0.1
% freqz (fir1 (15, [0.2, 0.5], "stop", "noscale"));

function b = fir1_octave(n, w, varargin)

% if nargin < 2 || nargin > 5
%     print_usage;
% end

%# Assign default window, filter type and scale.
%# If single band edge, the first band defaults to a pass band to
%# create a lowpass filter.  If multiple band edges, the first band
%# defaults to a stop band so that the two band case defaults to a
%# band pass filter.  Ick.
window  = [];
scale   = 1;
ftype   = (length(w)==1);

%# sort arglist, normalize any string
for i=1:length(varargin)
    arg = varargin{i};
    if ischar(arg), arg=lower(arg);end
    if isempty(arg) continue; end  % octave bug---can't switch on []
    switch arg
        case {'low','stop','dc-1'},             ftype  = 1;
        case {'high','pass','bandpass','dc-0'}, ftype  = 0;
        case {'scale'},                         scale  = 1;
        case {'noscale'},                       scale  = 0;
        otherwise                               window = arg;
    end
end

%# build response function according to fir2 requirements
bands = length(w)+1;
f = zeros(1,2*bands);
f(1) = 0; f(2*bands)=1;
f(2:2:2*bands-1) = w;
f(3:2:2*bands-1) = w;
m = zeros(1,2*bands);
m(1:2:2*bands) = rem([1:bands]-(1-ftype),2);
m(2:2:2*bands) = m(1:2:2*bands);

%# Increment the order if the final band is a pass band.  Something
%# about having a nyquist frequency of zero causing problems.
if rem(n,2)==1 && m(2*bands)==1,
    warning('n must be even for highpass and bandstop filters. Incrementing.');
    n = n+1;
    if isvector(window) && isreal(window) && ~ischar(window)
        %# Extend the window using interpolation
        M = length(window);
        if M == 1,
            window = [window; window];
        elseif M < 4
            window = interp1_(linspace(0,1,M),window,linspace(0,1,M+1),'linear');
        else
            window = interp1_(linspace(0,1,M),window,linspace(0,1,M+1),'spline');
        end
    end
end

%# compute the filter
b = fir2_(n, f, m, [], 2, window);

%# normalize filter magnitude
if scale == 1
    %# find the middle of the first band edge
    %# find the frequency of the normalizing gain
    if m(1) == 1
        %# if the first band is a passband, use DC gain
        w_o = 0;
    elseif f(4) == 1
        %# for a highpass filter,
        %# use the gain at half the sample frequency
        w_o = 1;
    else
        %# otherwise, use the gain at the center
        %# frequency of the first passband
        w_o = f(3) + (f(4)-f(3))/2;
    end
    
    %# compute |h(w_o)|^-1
    renorm = 1/abs(polyval_(b, exp(-1i*pi*w_o)));
    
    %# normalize the filter
    b = renorm*b;
end

end %func


%--------------------------------------------------------------------------
function [y, dy] = polyval_(p, x, s, mu)

if nargin<3, s=[]; end
% if (nargin < 2 || nargin > 4 || (nargout == 2 && nargin < 3))
%     print_usage ();
% end

if (isempty (x))
    y = [];
    return;
elseif (isempty (p))
    y = zeros (size (x));
    return;
elseif (~ isvector (p))
    error ('polyval: first argument must be a vector');
end

if (nargin > 3)
    x = (x - mu(1)) / mu(2);
end

n = length (p) - 1;
y = p(1) * ones (size (x));
for i = 2:n+1
    y = y .* x + p(i);
end

if (nargout == 2)
    %# Note: the F-Distribution is generally considered to be single-sided.
    %# http://www.itl.nist.gov/div898/handbook/eda/section3/eda3673.htm
    %#   t = finv (1-alpha, s.df, s.df);
    %#   dy = t * sqrt (1 + sumsq (A/s.R, 2)) * s.normr / sqrt (s.df)
    %# If my inference is correct, then t must equal 1 for polyval.
    %# This is because finv (0.5, n, n) = 1.0 for any n.
    sumsq_ = @(x)sum(x.^2);
    try
        k = numel (x);
        A = (x(:) * ones (1, n+1)) .^ (ones (k, 1) * (n:-1:0));
        dy = sqrt (1 + sumsq_(A/s.R, 2)) * s.normr / sqrt (s.df);
        dy = reshape (dy, size (x));
    catch
        if (isempty (s))
            error ('polyval: third input is required.');
        elseif (isstruct (s) ...
                && all (ismember ({'R', 'normr', 'df'}, fieldnames (s))))
            error (lasterr ());
        elseif (isstruct (s))
            error ('polyval: third input is missing the required fields.');
        else
            error ('polyval: third input is not a structure.');
        end
    end
end %if
end %func


%--------------------------------------------------------------------------
function b = fir2_(n, f, m, grid_n, ramp_n, window)

%   if nargin < 3 || nargin > 6
%     print_usage;
%   end

  %# verify frequency and magnitude vectors are reasonable
  t = length(f);
  if t<2 || f(1)~=0 || f(t)~=1 || any(diff(f)<0)
    error ('fir2: frequency must be nondecreasing starting from 0 and ending at 1');
  elseif t ~= length(m)
    error ('fir2: frequency and magnitude vectors must be the same length');
  %# find the grid spacing and ramp width
  elseif (nargin>4 && length(grid_n)>1) || ...
         (nargin>5 && (length(grid_n)>1 || length(ramp_n)>1))
    error ('fir2: grid_n and ramp_n must be integers');
  end
  if nargin < 4, grid_n=[]; end
  if nargin < 5, ramp_n=[]; end

  %# find the window parameter, or default to hamming
  w=[];
  if length(grid_n)>1, w=grid_n; grid_n=[]; end
  if length(ramp_n)>1, w=ramp_n; ramp_n=[]; end
  if nargin < 6, window=w; end
  if isempty(window), window=hamming_(n+1); end
  if ~isreal(window) || ischar(window), window=feval(window, n+1); end
  if length(window) ~= n+1, error ('fir2: window must be of length n+1'); end

  %# Default grid size is 512... unless n+1 >= 1024
  if isempty (grid_n)
    if n+1 < 1024
      grid_n = 512;
    else
      grid_n = n+1;
    end
  end

  %# ML behavior appears to always round the grid size up to a power of 2
  grid_n = 2 ^ nextpow2 (grid_n);

  %# Error out if the grid size is not big enough for the window
  if 2*grid_n < n+1
    error ('fir2: grid size must be greater than half the filter order');
  end

  if isempty (ramp_n), ramp_n = fix (grid_n / 25); end

  %# Apply ramps to discontinuities
  if (ramp_n > 0)
    %# remember original frequency points prior to applying ramps
    basef = f(:); basem = m(:);

    %# separate identical frequencies, but keep the midpoint
    idx = find (diff(f) == 0);
    f(idx) = f(idx) - ramp_n/grid_n/2;
    f(idx+1) = f(idx+1) + ramp_n/grid_n/2;
    f = [f(:);basef(idx)]';

    %# make sure the grid points stay monotonic in [0,1]
    f(f<0) = 0;
    f(f>1) = 1;
    f = unique([f(:);basef(idx)]'); ...

    %# preserve window shape even though f may have changed
    m = interp1_(basef, basem, f);

    %# axis([-.1 1.1 -.1 1.1])
    %# plot(f,m,'-xb;ramped;',basef,basem,'-or;original;'); pause;
  end

  %# interpolate between grid points
  grid = interp1_(f,m,linspace(0,1,grid_n+1)'); ...
  %# hold on; plot(linspace(0,1,grid_n+1),grid,'-+g;grid;'); hold off; pause;

  %# Transform frequency response into time response and
  %# center the response about n/2, truncating the excess
  if (rem(n,2) == 0)
    b = ifft([grid ; grid(grid_n:-1:2)]);
    mid = (n+1)/2;
    b = real ([ b([end-floor(mid)+1:end]) ; b(1:ceil(mid)) ]);
  else
    %# Add zeros to interpolate by 2, then pick the odd values below.
    b = ifft([grid ; zeros(grid_n*2,1) ;grid(grid_n:-1:2)]);
    b = 2 * real([ b([end-n+1:2:end]) ; b(2:2:(n+1))]);
  end

  %# Multiplication in the time domain is convolution in frequency,
  %# so multiply by our window now to smooth the frequency response.
  %# Also, for matlab compatibility, we return return values in 1 row
  b = b(:)' .* window(:)';

end


%--------------------------------------------------------------------------
function c = hamming_(m, opt)

% if (nargin < 1 || nargin > 2)
%     print_usage ();
% end

if (~ (isscalar (m) && (m == fix (m)) && (m > 0)))
    error ('hamming: M must be a positive integer');
end

N = m - 1;
if (nargin == 2)
    switch (opt)
        case 'periodic'
            N = m;
        case 'symmetric'
            %# Default option, same as no option specified.
        otherwise
            error ('hamming: window type must be either ''periodic'' or ''symmetric''');
    end
end

if (m == 1)
    c = 1;
else
    m = m - 1;
    c = 0.54 - 0.46 * cos (2 * pi * (0 : m)' / N); ...
end

end %func


%--------------------------------------------------------------------------
function yi = interp1_(x, y, varargin)

%   if (nargin < 2 || nargin > 6)
%     print_usage ();
%   end

method = 'linear';
extrap = [];
xi = [];
ispp = false;
firstnumeric = true;
rightcontinuous = NaN;

if (nargin > 2)
    for i = 1:length (varargin)
        arg = varargin{i};
        if (ischar (arg))
            arg = lower(arg);
            switch (arg)
                case 'extrap'
                    extrap = 'extrap';
                case 'pp'
                    ispp = true;
                case {'right', '-right'}
                    rightcontinuous = true;
                case {'left', '-left'}
                    rightcontinuous = false;
                otherwise
                    method = arg;
            end
        else
            if (firstnumeric)
                xi = arg;
                firstnumeric = false;
            else
                extrap = arg;
            end
        end
    end
end

if (isempty (xi) && firstnumeric && ~ ispp)
    xi = y;
    y = x;
    if (isvector (y))
        x = 1:numel (y);
    else
        x = 1:size(y,1);
    end
end

%# reshape matrices for convenience
x = x(:);
nx = size(x,1);
szx = size (xi);
if (isvector (y))
    y = y(:);
end

szy = size (y);
y = y(:,:);
[ny, nc] = size (y);
xi = xi(:);

%# determine sizes
if (nx < 2 || ny < 2)
    error ('interp1: minimum of 2 points required in each dimension');
end

%# check whether x is sorted; sort if not.
if (~ issorted (x, 'strictmonotonic'))
    [x, p] = sort (x);
    y = y(p,:);
end

if (any (strcmp (method, {'previous', '*previous', 'next', '*next'})))
    rightcontinuous = NaN; % needed for these methods to work
end

if (isnan (rightcontinuous))
    %# If not specified, set the continuity condition
    if (x(end) < x(1))
        rightcontinuous = false;
    else
        rightcontinuous = true;
    end
elseif ((rightcontinuous && (x(end) < x(1))) ...
        || (~ rightcontinuous && (x(end) > x(1))))
    %# Switch between left-continuous and right-continuous
    x = flipud (x);
    y = flipud (y);
end

%# Because of the way mkpp works, it's easiest to implement 'next'
%# by running 'previous' with vectors flipped.
if (strcmp (method, 'next'))
    x = flipud (x);
    y = flipud (y);
    method = 'previous';
elseif (strcmp (method, '*next'))
    x = flipud (x);
    y = flipud (y);
    method = '*previous';
end

starmethod = method(1) == '*';

if (starmethod)
    dx = x(2) - x(1);
else
    jumps = x(1:end-1) == x(2:end);
    have_jumps = any (jumps);
    if (have_jumps)
        if (strcmp (method, 'linear') || strcmp (method, ('nearest')))
            if (any (jumps(1:nx-2) & jumps(2:nx-1)))
                warning ('interp1: multiple discontinuities at the same X value');
            end
        else
            error ('interp1: discontinuities not supported for method ''%s''', method);
        end
    end
end

%# Proceed with interpolating by all methods.
switch (method)
    
    case 'nearest'
        pp = mkpp ([x(1); (x(1:nx-1)+x(2:nx))/2; x(nx)], ...
            shiftdim (y, 1), szy(2:end));
        pp.orient = 'first';
        
        if (ispp)
            yi = pp;
        else
            yi = ppval (pp, reshape (xi, szx));
        end
        
    case '*nearest'
        pp = mkpp ([x(1), x(1)+[0.5:(nx-1)]*dx, x(nx)], ...
            shiftdim (y, 1), szy(2:end));
        pp.orient = 'first';
        
        if (ispp)
            yi = pp;
        else
            yi = ppval (pp, reshape (xi, szx));
        end
        
    case 'previous'
        pp = mkpp ([x(1:nx); 2*x(nx)-x(nx-1)], ...
            shiftdim (y, 1), szy(2:end));
        pp.orient = 'first';
        
        if (ispp)
            yi = pp;
        else
            yi = ppval (pp, reshape (xi, szx));
        end
        
    case '*previous'
        pp = mkpp (x(1)+[0:nx]*dx, ...
            shiftdim (y, 1), szy(2:end));
        pp.orient = 'first';
        
        if (ispp)
            yi = pp;
        else
            yi = ppval (pp, reshape (xi, szx));
        end
        
    case 'linear'
        
        xx = x;
        nxx = nx;
        yy = y;
        dy = diff (yy);
        if (have_jumps)
            %# Omit zero-size intervals.
            xx(jumps) = [];
            nxx = size(xx,1);
            yy(jumps, :) = [];
            dy(jumps, :) = [];
        end
        
        dx = diff (xx);
        size_dy_ = size(dy);
        dx = repmat (dx, [1 size_dy_(2:end)]);
        
        coefs = [(dy./dx).', yy(1:nxx-1, :).'];
        
        pp = mkpp (xx, coefs, szy(2:end));
        pp.orient = 'first';
        
        if (ispp)
            yi = pp;
        else
            yi = ppval (pp, reshape (xi, szx));
        end
        
    case '*linear'
        dy = diff (y);
        toVec_ = @(x)x(:);
        coefs = [toVec_(dy/dx), toVec_(y(1:nx-1, :))]; % may need to transpose
        pp = mkpp (x, coefs, szy(2:end));
        pp.orient = 'first';
        
        if (ispp)
            yi = pp;
        else
            yi = ppval (pp, reshape (xi, szx));
        end
        
    case {'pchip', '*pchip', 'cubic', '*cubic'}
        if (nx == 2 || starmethod)
            x = linspace (x(1), x(nx), ny);
        end
        
        if (ispp)
            y = shiftdim (reshape (y, szy), 1);
            yi = pchip (x, y);
            yi.orient = 'first';
        else
            y = shiftdim (y, 1);
            yi = pchip (x, y, reshape (xi, szx));
            if (~ isvector (y))
                yi = shiftdim (yi, 1);
            end
        end
        
    case {'spline', '*spline'}
        if (nx == 2 || starmethod)
            x = linspace (x(1), x(nx), ny);
        end
        
        if (ispp)
            y = shiftdim (reshape (y, szy), 1);
            yi = spline (x, y);
            yi.orient = 'first';
        else
            y = shiftdim (y, 1);
            yi = spline (x, y, reshape (xi, szx));
            if (~ isvector (y))
                yi = shiftdim (yi, 1);
            end
        end
        
    otherwise
        error ('interp1: invalid method ''%s''', method);
        
end

if (~ ispp && isnumeric (extrap))
    %# determine which values are out of range and set them to extrap,
    %# unless extrap == 'extrap'.
    minx = min (x(1), x(nx));
    maxx = max (x(1), x(nx));
    
    xi = reshape (xi, szx);
    outliers = xi < minx | ~ (xi <= maxx); % this even catches NaNs
    if numel(outliers) == numel(yi)
        yi(outliers) = extrap;
        yi = reshape (yi, szx);
    elseif (~ isvector (yi))
        yi(outliers, :) = extrap;
    else
        yi(outliers.') = extrap; ...
    end

end

end %func

