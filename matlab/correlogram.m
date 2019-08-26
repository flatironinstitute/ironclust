function trCC = correlogram(S_clu, viTime_spk, P)
burst_interval_ms = 20;
binsize_ms = 1;

binsize = round(binsize_ms * P.sRateHz/1000); % 1 ms bin
nbins_shift = ceil(burst_interval_ms / binsize_ms);
cviTime_clu = cellfun_(@(x)viTime_spk(x), S_clu.cviSpk_clu);
max_time = double(max(cellfun(@max, cviTime_clu)));
viShift = int32([-nbins_shift:nbins_shift]); % 20 msec window
nShift = numel(viShift);
nClu = numel(cviTime_clu);
trCC = zeros(nShift, nClu, nClu, 'single');
cviTime_b_clu = cellfun_(@(x)int32(double(x(:))/binsize), cviTime_clu);
fprintf('Computing correlogram\n\t'); t1=tic;
for iClu1 = 1:nClu
%     vlTime1 = sparse(max_time,1);
%     vlTime1(cviTime_b_clu{iClu1}) = 1;
    vlTime1 = false(max_time,1);
    vlTime1(cviTime_b_clu{iClu1}) = true;
    mr_ = zeros(nShift, nClu, 'single');
    for iClu2 = 1:nClu
        if iClu1 == iClu2, continue; end
%         vr_ = arrayfun(@(x)sum(ismembc(cviTime_b_clu{iClu2} + x, viTime1)), viShift);
        miTime2 = cviTime_b_clu{iClu2} + viShift;
        miTime2 = min(max(miTime2,1),max_time);
        mr_(:,iClu2) = sum(vlTime1(miTime2)) / numel(cviTime_b_clu{iClu2}); 
    end
    trCC(:,:,iClu1) = conv2(mr_, [.25;.5;.25], 'same');
    fprintf('.');
end
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function a=cellfun_(varargin)
a=cellfun(varargin{:}, 'UniformOutput',0);
end %func