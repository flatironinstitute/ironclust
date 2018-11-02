% Non-blocking call
[csFile_batch, csFile_end] = deal({});
for noise = [10,20]
    for K = [10, 20]
        for C = [4 8]
            vcDir_in = sprintf('~/ceph/groundtruth/magland_synth/datasets_noise%d_K%d_C%d', noise, K, C);
            vcDir_out = sprintf('~/ceph/groundtruth_irc/magland_synth/datasets_noise%d_K%d_C%d', noise, K, C);
            [csFile_batch{end+1}, csFile_end{end+1}] = irc('sbatch-mda', vcDir_in, vcDir_out, 'tetrode_template.prm', 0); % non-blocking call
        end
    end
end


% Wait
disp('-----------------------------------------------------------------');
irc('waitfor', csFile_end, 3600); % wait for 1 hour max


% Plot
disp('-----------------------------------------------------------------');
for iFile=1:numel(csFile_end)
    irc('batch-verify', csFile_batch{iFile}, 'skip');
end