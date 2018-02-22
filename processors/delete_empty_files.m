function delete_empty_files(vcDir)
if nargin<1, vcDir=[]; end
delete_files_(find_empty_files_(vcDir));
end %func


function csFiles = find_empty_files_(vcDir)
% find files with 0 bytes
if nargin==0, vcDir = []; end
if isempty(vcDir), vcDir = pwd(); end
vS_dir = dir(vcDir);
viFile = find([vS_dir.bytes] == 0 & ~[vS_dir.isdir]);
csFiles = {vS_dir(viFile).name};
csFiles = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end %func


function delete_files_(csFiles)
for iFile = 1:numel(csFiles)
    try
        if exist(csFiles{iFile}, 'file')
            delete(csFiles{iFile});
            fprintf('\tdeleted %s.\n', csFiles{iFile});
        end
    catch
        disperr_();
    end
end
end %func