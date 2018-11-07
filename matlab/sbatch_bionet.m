% Non-blocking call
csStudies = {'static', 'drift', 'shuffle'};
% csStudies = {'static'};
vcDir_in = '/mnt/ceph/users/jjun/groundtruth/bionet/bionet_%s';
vcDir_out = '/mnt/ceph/users/jjun/groundtruth_irc/bionet/bionet_%s';
vcFile_template = '/mnt/home/jjun/src/ironclust/matlab/drift_study.prm'; 
fRun_remote = 1;

[csFile_batch, csFile_end] = deal({});
for iStudy = 1:numel(csStudies)
    vcDir_in1 = sprintf(vcDir_in, csStudies{iStudy});
    vcDir_out1 = sprintf(vcDir_out, csStudies{iStudy});
    if fRun_remote
        [csFile_batch{end+1}, csFile_end{end+1}] = irc('sbatch-mda', vcDir_in1, vcDir_out1, vcFile_template, 0); % non-blocking call
    else
        csFile_batch{end+1} = irc('batch-mda', vcDir_in1, vcDir_out1, vcFile_template);
    end %switch
end

%% Wait
if fRun_remote
    disp('-----------------------------------------------------------------');
    irc('waitfor', csFile_end, 3600); % wait for 1 hour max
end

%% Plot
disp('-----------------------------------------------------------------');
for iFile=1:numel(csFile_end)
    irc('batch-verify', csFile_batch{iFile}, 'skip');
end

% fRun_remote: 383.5s, 92.8%
