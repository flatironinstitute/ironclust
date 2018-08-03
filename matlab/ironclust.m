function ironclust(varargin)
% calls jrclust
warning off;
%jrclust(varargin{:});
fprintf('Running ''%s%sirc.m''\n', fileparts(mfilename('fullpath')), filesep());
irc(varargin{:});