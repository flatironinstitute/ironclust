function [S_mda, fid_r] = readmda_header(fname)
%READMDADIMS - read only dimensions of a .mda file. MDA stands for multi-dimensional array.
%
% See http://magland.github.io//articles/mda-format/
%
% Syntax: dims=readmdadims(fname)
%
% Inputs:
%    fname - path to the .mda file
%
% Outputs:
%    dims - row vector of dimension sizes of multi-dimensional array
%
% Other m-files required: none
%
% See also: writemda

% Author: Alex Barnett 7/22/16

fid_r = fopen(fname,'rb');

try
    code=fread(fid_r,1,'int32');
catch
    error('Problem reading file: %s',fname);
end
if (code>0) 
    num_dims=code;
    code=-1;
    nBytes_sample = 4;
else
    nBytes_sample = fread(fid_r,1,'int32');
    num_dims=fread(fid_r,1,'int32');    
end
dim_type_str='int32';
if (num_dims<0)
    num_dims=-num_dims;
    dim_type_str='int64';
end

dimm=zeros(1,num_dims);
for j=1:num_dims
    dimm(j)=fread(fid_r,1,dim_type_str);
end

S_mda = struct('dimm', dimm, 'vcDataType', getDataType_(code), ...
    'nBytes_header', ftell(fid_r), 'nBytes_sample', nBytes_sample);

if nargout<2, fclose(fid_r); end
end %func


%--------------------------------------------------------------------------
function vcDataType = getDataType_(code)
if (code==-1)
    vcDataType = 'float';
elseif (code==-2)
    vcDataType = 'uchar';
elseif (code==-3)
    vcDataType = 'float';
elseif (code==-4)
    vcDataType = 'int16';
elseif (code==-5)
    vcDataType = 'int32';
elseif (code==-6)
    vcDataType = 'uint16';
elseif (code==-7)
    vcDataType = 'double';
elseif (code==-8)
    vcDataType = 'uint32';
else
    vcDataType = ''; % unknown
end %if
end %func