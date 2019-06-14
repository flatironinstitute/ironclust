function mrDist_clu = post_merge_drift(S0, tnWav_spk)

if 0
    load S_clu_drift      % contains S0 P tnWav_spk
end

fprintf('post_merge_drift\n'); 
S_clu = S0.S_clu;
P = S0.P;

viSite_spk = S0.viSite_spk;
S_drift = S_clu.S_drift;
nBatch_drift = S_drift.nTime_drift;
miBatch_drift = S_drift.miSort_drift;
nClu = S_clu.nClu;
nSpk_min = P.min_count;

% do a better merging scheme
% compute average waveforms for each cluster for each drift state
fprintf('\tComputing average waveforms per unit per batch\n\t\t'); t1=tic;
cm_vrWav_batch_clu = cell(nBatch_drift, nClu);
miSite_batch_clu = zeros(nBatch_drift, nClu);
for iClu = 1:nClu    
    viSpk_clu1 = S_clu.cviSpk_clu{iClu};
    viDrift_spk1 = S_drift.viDrift_spk(viSpk_clu1);
    for iBatch = 1:nBatch_drift
        viSpk1 = viSpk_clu1(viDrift_spk1 == iBatch);
        if numel(viSpk1) < nSpk_min, continue; end
        viSite1 = viSite_spk(viSpk1);
        iSite1 = mode(viSite1);
        viSpk2 = viSpk1(viSite1==iSite1);
        if numel(viSpk2) < nSpk_min, continue; end
        miSite_batch_clu(iBatch, iClu) = iSite1;
        mr_ = mean(single(tnWav_spk(:,:,viSpk2)), 3);
        vr_ = mr_(:);
        cm_vrWav_batch_clu{iBatch, iClu} = (vr_ - mean(vr_)) / std(vr_,1); % / sqrt(numel(vr_));
    end %for
    fprintf('.');
end %for
fprintf('\n\t\ttook %0.1fs\n', toc(t1));

% cm_vrWav_batch_clu = cellfun(@(x)zscore(x(:)), cm_mrWav_batch_clu, 'UniformOutput', 0);

% compute waveform similarities
fprintf('\tComputing waveform similarity between unit pairs\n\t\t'); t2=tic;
mrDist_clu = zeros(nClu);
for iClu1 = 1:nClu
    vrDist_clu1 = zeros(nClu,1);
    viBatch1 = find(miSite_batch_clu(:,iClu1) > 0);
    for iiBatch1 = 1:numel(viBatch1)
        iBatch11 = viBatch1(iiBatch1);
        iSite11 = miSite_batch_clu(iBatch11, iClu1);        
        vrWav_clu1_T = cm_vrWav_batch_clu{iBatch11, iClu1}';
        viBatch11 = miBatch_drift(:,iBatch11); % accept clusters from this        
        [viiBatch12, viClu12] = find(miSite_batch_clu(viBatch11,:) == iSite11);
        vl_ = viClu12 ~= iClu1;
        viBatch12 = viBatch11(viiBatch12(vl_));
        viClu12 = viClu12(vl_);
        for ii12 = 1:numel(viClu12)
            iClu2 = viClu12(ii12);
%             if iClu2 == iClu1, continue; end
            iBatch2 = viBatch12(ii12);
            vrWav_clu2 = cm_vrWav_batch_clu{iBatch2, iClu2};
            vrDist_clu1(iClu2) = max(vrDist_clu1(iClu2), vrWav_clu1_T * vrWav_clu2);
        end
    end
    mrDist_clu(:,iClu1) = vrDist_clu1;
    fprintf('.');
end
mrDist_clu = mrDist_clu / (size(tnWav_spk,1)*size(tnWav_spk,2));
fprintf('\n\t\ttook %0.1fs\n', toc(t2));
end %func