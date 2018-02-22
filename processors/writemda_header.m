function fid_w = writemda_header(vcFile, S_mda)
fid_w = fopen(vcFile, 'w');

fwrite(fid_w, getCode_(S_mda.vcDataType), 'int32');
fwrite(fid_w, numel(S_mda.dimm), 'int32');
fwrite(fid_w, S_mda.dimm, 'int32');
if nargout==0, fclose(fid_w); end
end %func



%--------------------------------------------------------------------------
function vnCode = getCode_(dtype)
switch(lower(dtype))
    case 'complex32', vnCode = [-1,8];
    case {'float32', 'single'}, vnCode = [-3 4];
    case {'float64', 'double'}, vnCode = [-7,8];
    case 'int32', vnCode = [-5,4];
    case 'int16', vnCode = [-4 2];
    case 'uint16', vnCode = [-6,2];
    case 'uint32', vnCode = [-8,4];
    otherwise, vnCode = [];
end %switch
end %func