% copy english files

path_from = '/mnt/ceph/users/jjun/DanEnglish/juxta_cell_curated/';
path_to = '/mnt/ceph/users/jjun/groundtruth/paired_english/';

vS_dir = dir(path_from);
vl_dir = [vS_dir.isdir];
csName = {vS_dir.name};
csName = csName(vl_dir);
csName = setdiff(csName, {'.','..'});
csPath_from = cellfun(@(x)fullfile(path_from, x), csName, 'UniformOutput',0);
csPath_to = cellfun(@(x)fullfile(path_to, x), csName, 'UniformOutput',0);

csFiles_copy = {'raw.mda', 'params.json', 'firings_true.mda', 'geom.csv'};
for iRec = 1:numel(csPath_from)
    vcDir_from1 = csPath_from{iRec};
    vcDir_to1 = csPath_to{iRec};
    if ~exist(vcDir_to1, 'dir'), mkdir(vcDir_to1); end
    try
        cellfun(@(x)copyfile(fullfile(vcDir_from1,x), fullfile(vcDir_to1,x), 'f'), csFiles_copy);
    catch
        fprintf(2, 'error copying %s\n', vcDir_from1);
    end
    cellfun(@(x)fprintf('copyfile %s %s\n', fullfile(vcDir_from1,x), fullfile(vcDir_to1,x)), csFiles_copy);
    fprintf('\n');
end