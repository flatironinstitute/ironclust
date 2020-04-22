%% 1. Generate linear drift
% copy and run scripts below
% /mnt/ceph/users/jjun/groundtruth/hybrid_synth/linear_drift/make_eMouseData_drift.m
% /mnt/ceph/users/jjun/groundtruth/hybrid_synth/linear_drift/master_eMouse_drift.m


%% 2. convert the output of eMouse simulation to mda format
% /mnt/ceph/users/jjun/groundtruth/hybrid_synth/linear_drift
vcDir_ksort = '/mnt/ceph/users/jjun/groundtruth/hybrid_synth/linear_static/';
irc2('emouse2mda', vcDir_ksort); 


%% 3. cut and generate 10 recordings each
nRecordings = 10;
nChans = 64;
vnChans_out = [8, 16];
vcDir_out = '/mnt/ceph/users/jjun/groundtruth/long_static/';
vnT_seg = 2.^(0:4);
t_dur_sec = 80*60;
fTest = 0;

for nChans_out = vnChans_out
    viChanA = round(linspace(1, nChans-nChans_out+1, nRecordings));
    for iT_seg = 1:numel(vnT_seg)
        nT_seg = vnT_seg(iT_seg);
        vcDir0 = sprintf('%ds_%dc', round(t_dur_sec/nT_seg), nChans_out);
        for iRec = 1:nRecordings
            chanLim1 = [0, nChans_out-1] + viChanA(iRec);
            vcDir_out1 = fullfile(vcDir_out, vcDir0, sprintf('%03d_synth', iRec));
            vcArg1 = sprintf('%s %s %d:%d [0,1/%d]', vcDir_ksort, vcDir_out1, chanLim1(1), chanLim1(2), nT_seg);
            vcCmd1 = sprintf('irc2 extract-mda %s', vcArg1);
            fprintf(vcCmd1);
            fprintf('\n');
            if ~fTest
                t1=tic;
                eval(vcCmd1);
                fprintf('Took %0.1fs\n', toc(t1));
            end
        end
        fprintf('\n');
    end
    fprintf('\n');
end
