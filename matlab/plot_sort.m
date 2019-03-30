function plot_sort(varargin)
hold on;
for iArg = 1:nargin
    vy = varargin{iArg};
    plot((1:numel(vy))/numel(vy), sort(vy));
end