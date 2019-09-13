
% load data
% todo: ues distributed file loading from ceph


% detect
cS_detect = cell(nBatch,1);
plan_detect_(P); % if multifile, parfor file. if long recording, parfor ...
cviLim_batch = cell(nBatch, 1);

parfor iBatch = 1:nBatch
    cS_detect{iBatch} = detect_batch_(P, cviLim_batch{iBatch});
end
% sort



% merge