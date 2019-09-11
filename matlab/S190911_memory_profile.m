% memory profile loop
% paramteric test
%   fSave_spkwav {0, 1}
%   nChans, timeDuration
%
%   ToDo: extend time and channels, add noise
%       plan: use compiled matlab (./run_irc myfileloc)

% loop over the files, exract values

% change paramter

% loop over the files
vcDir = '~/raid/groundtruth/hybrid_synth/static_siprobe';
csFiles = {'rec_16c_1200s_11', 'rec_16c_600s_11', 'rec_32c_1200s_11', 'rec_32c_600s_11', 'rec_64c_1200s_11', 'rec_64c_600s_11'};
nParams = 1;

for iParam = 1:nParams
    for iFile = 1:numel(csFiles)
        vcFile1 = fullfile(vcDir, csFiles{iFile});
        [~, out1] = system(sprintf('./run_irc %s', vcFile1));
        cs1 = strsplit(out1, '\n');
        cvi1 = strfind(cs1, 'memory usage (GiB):');
        iLine1 = find(~cellfun(@isempty, cvi1), 1, 'first');
        cs2 = strsplit(cs1{iLine1}, ' ');
        memory_gb1 = str2double(cs2{end});
        mrPeakMem(iParam, iFile) = memory_gb1;
    end
end

disp(csFiles)
disp(mrPeakMem)