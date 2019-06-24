function jrclust(varargin)
% calls irc.m

warning off;
fprintf('Running ''%s%sirc.m''\n', fileparts(mfilename('fullpath')), filesep());
irc(varargin{:});