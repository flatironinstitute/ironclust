%% drift merge
[S_clu, P] = deal(S0.S_clu, S0.P);
% cluster 7 and 9
cviSpk12 = S0.S_clu.cviSpk_clu([7,9]);
[viSpk1, viSpk2] = deal(cviSpk12{1}, cviSpk12{2});
miKnn1 = S0.S_clu.miKnn(:, viSpk1);
miKnn2 = S0.S_clu.miKnn(:, viSpk2);
fh_col = @(x)x(:);
cviSpk_drift = S0.S_clu.S_drift.cviSpk_drift;

% if drift correction is not done
if numel(cviSpk_drift) == 1
    nSpk = numel(S0.viTime_spk);
    t_dur = single(diff(S0.viTime_spk([1,end]))) / P.sRateHz;
    t_step = P.step_sec_drift;
    nDrift = ceil(t_dur / t_step);
    viSpk_drift = round(linspace(0, nSpk, nDrift+1));
    cviSpk_drift = arrayfun(@(x,y) [x+1:y], viSpk_drift(1:end-1), ...
        viSpk_drift(2:end), 'UniformOutput', 0);
end


%% two cluster examples
knn = min(30, P.knn);
[vpDrift12, vpDrift21, vpDrift11, vpDrift22] = deal(nan(size(cviSpk_drift)));
for iDrift = 1:numel(cviSpk_drift)
    viiSpk1 = find(ismember(viSpk1, cviSpk_drift{iDrift}));
    viiSpk2 = find(ismember(viSpk2, cviSpk_drift{iDrift}));
    
    if numel(viiSpk1) >= P.min_count
        viSpk_out1 = fh_col(miKnn1(1:knn,viiSpk1));
        vpDrift12(iDrift) = mean(ismember(viSpk_out1, viSpk2));
        vpDrift11(iDrift) = mean(ismember(viSpk_out1, viSpk1));        
    end
    if numel(viiSpk2) >= P.min_count
        viSpk_out2 = fh_col(miKnn2(1:knn,viiSpk2));
        vpDrift21(iDrift) = mean(ismember(viSpk_out2, viSpk1));
        vpDrift22(iDrift) = mean(ismember(viSpk_out2, viSpk2));    
    end
end
figure; bar([vpDrift12; vpDrift21; vpDrift11; vpDrift22]'); legend({'12','21','11','22'});
figure; bar([vpDrift12+vpDrift11; vpDrift21+vpDrift22]'); legend({'12+11','21+22'});

%% non-drift merge
mean(ismember(miKnn1(:), viSpk2))
mean(ismember(miKnn2(:), viSpk1))


%% global matrix
tic
[nDrift, nClu] = deal(numel(cviSpk_drift), S_clu.nClu);
[trDist_out, trDist_in] = deal(nan(nClu, nDrift, nClu));
mpSelf_drift_clu = nan(nDrift, nClu);
fh_member = @(x,y)mean(ismember(x,y));
for iClu = 1:nClu
    viSpk1 = S_clu.cviSpk_clu{iClu};
    miKnn1 = S0.S_clu.miKnn(:, viSpk1);
    for iDrift = 1:nDrift        
        viiSpk1 = find(ismember(viSpk1, cviSpk_drift{iDrift}));
        if numel(viiSpk1) >= P.min_count
            viSpk_out1 = fh_col(miKnn1(1:knn,viiSpk1));
            vr_ = cellfun(@(y)fh_member(viSpk_out1, y), S_clu.cviSpk_clu);
            trDist_out(:,iDrift,iClu) = vr_;
%             assert(sum(trDist_out(:,iDrift,iClu)) <= 1, 'should be');
%             viSpk_in1 = viSpk1(viiSpk1);
%             trDist_in(:,iDrift,iClu) = cellfun(@(y)fh_member(viSpk_in1, y), S_clu.cviSpk_clu);
        end
    end %for    
    fprintf('.');
end %for
toc

%% find the highest pair
[vrPre_clu, vrPos_clu, vpRatio_clu, viMerge_clu] = deal(zeros(nClu,1));
for iClu = 1:nClu
    mr1 = trDist_out(:,:,iClu);    
    vr1 = mr1(iClu,:);
    mr1(iClu,:) = nan;
    [v1,i1] = max(mr1(:), [], 'omitnan');
    [iClu1, iDrift1] = ind2sub(size(mr1), i1);
    v2 = vr1(iDrift1) + v1; 
    
    [vrPre_clu(iClu), vrPost_clu(iClu), vpRatio_clu(iClu), viMerge_clu(iClu)] = ...
        deal(v1, v2, v1 / vr1(iDrift1), iClu1);
end
figure; plot(vrPre_clu, vrPre_clu(viMerge_clu), '.');
axis([0 1 0 1]);

figure; bar([vrPre_clu(:), vrPost_clu(:)]);
figure; bar(vpRatio_clu); ylabel('ratio out/in');

viClu_merge = find(vpRatio_clu>1/8);
disp(viClu_merge(:)');
disp(viMerge_clu(viClu_merge)');


%%
mrDist_clu = squeeze(max(trDist_out,[],2));
fh_diag = @(x)x(sub2ind(size(x),1:size(x,1), 1:size(x,1)));
fh_diag_idx = @(x)sub2ind(size(x),1:size(x,1), 1:size(x,1));
mrDist_clu1 = mrDist_clu;
mrDist_clu1(fh_diag_idx(mrDist_clu)) = 0;
figure; bar(fh_diag(mrDist_clu),1);
% figure; plot(fh_diag(mrDist_clu) + max(mrDist_clu1));