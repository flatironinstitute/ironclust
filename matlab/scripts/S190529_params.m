

csDir = irc('call', 'dir_', {'/mnt/home/jjun/ceph/groundtruth/paired_recordings/boyden32c/*/params.json'});

for iDir = 1:numel(csDir)
    vcDir1 = csDir{iDir};
    vcDir2 = strrep(vcDir1, 'boyden32c', 'boyden');
    S_json2 = loadjson(vcDir2);
    S_json2.scale_factor = .195;
    savejson(vcDir1, S_json2);
end

%% all json files

cs_json = irc('call', 'dir', {'/mnt/ceph/users/jjun/groundtruth/**/params.json'});

vl_scale = false(size(cs_json));
for i=1:numel(cs_json)
    S1 = loadjson(cs_json{i});
    
    vl_scale(i) = isfield(S1, 'scale_factor');    
end %for




%% save kampff dataset as integer




