function trCC = burst_merge(cviTime_clu_in)
persistent cviTime_clu
if isempty(cviTime_clu_in)
    if isempty(cviTime_clu)
        load cviTime_clu
    end
else
    cviTime_clu = cviTime_clu_in;
end

max_time = max(cellfun(@max, cviTime_clu));
binsize = 32 * 1; %2 msec bin
viShift = int32([-20:20]);
nShift = numel(viShift);
nClu = numel(cviTime_clu);
trCC = zeros(nShift, nClu, nClu);
cviTime_b_clu = cellfun(@(x)ceil(x(:)/binsize), cviTime_clu, 'UniformOutput', 0);
tic
mr_ = zeros(nShift, nClu);
for iClu1 = 1:nClu
    vlTime1 = false(max_time,1);
    vlTime1(cviTime_b_clu{iClu1}) = true;
    n1 = numel(cviTime_b_clu{iClu1});    
    for iClu2 = 1:nClu
        if iClu1 == iClu2, continue; end
        miTime2 = cviTime_b_clu{iClu2} + viShift;
        miTime2= max(min(miTime2,max_time),1);
        n2 = numel(cviTime_b_clu{iClu2});        
        mr_(:,iClu2) = sum(vlTime1(miTime2)) / n2;
%         vr_ = sum(vlTime1(miTime2));
%         vr_ = conv(vr_,[1 2 1]/4/n2,'same');
%         trCC(:,iClu2,iClu1) = vr_;
    end
    trCC(:,:,iClu1) = conv2(mr_, [.25;.5;.25], 'same');
end

end %func