function FF = writemda_fid(FF, X)
%WRITEMDA_FID - write to a .mda file. MDA stands for
%multi-dimensional array.
% outputs to 64-bit dimension format
%
% Usage
% -----
% FF = writemda_fid(FF, X)
%    write the data to the file
% FF = writemda_fid(file_name, X)
%    write the data to the file
% FF = writemda_fid(file_name)
%    open file with 'w+' mode
% writemda_fid(FF, 'close')
%    close the file and update the header
% writemda_fid(file_name, X)
%    Write X to 'file_name' and close
%
% Inputs:
%    X - the multi-dimensional array
%    FF - fid (file pointer, write permission)
%
% Other m-files required: none
%
% See also: readmda

% Author: James Jun
% 11/25/2019: Modified from writemda.m
% 12/10/2019: Persistent variables removed, using the header instead

dim_type_str = 'int64'; %either int64 or int32
fClose = 0;
if nargin<2, X=[]; end
switch dim_type_str
    case 'int32', DIMM_BYTE = 4;
    case 'int64', DIMM_BYTE = 8;
    otherwise, error('dim_type_str must be either ''int32'' or ''int64'''); 
end
if ischar(FF)
    FF = fopen(FF, 'w+', 'l');
    if nargout==0, fClose = 1; end
    if nargin==1, return; end
end
nbytes_now = ftell(FF); % update header if 

% close file, update the header (last dimension)
if ischar(X)
    if strcmpi(X, 'close')
        % update the last index    
        frewind(FF);
        vec_ = fread(FF,3,'int32');
        [type_, bytesPerSample_, ndims_] = deal(vec_(1), vec_(2), abs(vec_(3)));
        dimm_ = fread(FF, ndims_, dim_type_str);
        
        % write the last sample
        fseek(FF, -DIMM_BYTE, 'cof');
        nBytes_header_ = 3 * 4 + ndims_ * DIMM_BYTE;
        nSamples = (nbytes_now - nBytes_header_) / bytesPerSample_;
        if numel(dimm_) > 1
            last_dimm = nSamples / prod(dimm_(1:end-1));
        else
            last_dimm = nSamples;
        end
        fwrite(FF, last_dimm, dim_type_str); % update last dimension             
        fclose(FF);        
        return;
    end
end

% write the file header
if nbytes_now == 0
    if isa(X,'single')
        fwrite(FF,[-3,4],'int32');
    elseif isa(X,'double')
        fwrite(FF,[-7,8],'int32');
    elseif isa(X,'int32')
        fwrite(FF,[-5,4],'int32');
    elseif isa(X,'int16')
        fwrite(FF,[-4,2],'int32');
    elseif isa(X,'uint16')
        fwrite(FF,[-6,2],'int32');
    elseif isa(X,'uint32')
        fwrite(FF,[-8,4],'int32');
    else
        error('Unknown dtype %s', class(X));
    end
    factor_ = strcmpi(dim_type_str, 'int32')*2-1;
    fwrite(FF, factor_ * ndims(X), 'int32');
    fwrite(FF, size(X), dim_type_str);
end

fwrite(FF, X, class(X));
if fClose, writemda_fid(FF, 'close'); end
end

