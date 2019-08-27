%--------------------------------------------------------------------------
% 7/2/2019 JJJ: Merge using knn relationship, not using waveforms (fast)
function S_clu = post_merge_knn1(S_clu, P)

KNN = 16;
merge_thresh = get_set_(P, 'out_in_ratio_merge', 1/8);
knn = min(KNN, size(S_clu.miKnn,1));
miKnn = S_clu.miKnn(1:knn,:);
vrRho_spk = S_clu.rho(:)';

fprintf('Automated merging using KNN (post-hoc)\n\t'); t_merge=tic;

%-----
% compute distances
try
    cviSpk_drift = S_clu.S_drift.cviSpk_drift;
catch
    cviSpk_drift = {};
end

% if drift correction is not done
if numel(cviSpk_drift) <= 1
    nSpk = numel(S_clu.viClu);
    t_dur = irc('call', 'recording_duration_', {P});
    t_step = P.step_sec_drift;
    nDrift = ceil(t_dur / t_step);
    viSpk_drift = round(linspace(0, nSpk, nDrift+1));
    cviSpk_drift = arrayfun(@(x,y) [x+1:y], viSpk_drift(1:end-1), ...
        viSpk_drift(2:end), 'UniformOutput', 0);
end


%-----
% global matrix
[nDrift, nClu] = deal(numel(cviSpk_drift), S_clu.nClu);
trDist_clu = nan(nClu, nDrift, nClu);
fh_member = @(x,y)mean(ismember(x,y));    
miSites_clu = P.miSites(1:P.nSites_fet, S_clu.viSite_clu);
for iClu = 1:nClu
    viSpk1 = S_clu.cviSpk_clu{iClu};    
    miKnn1 = miKnn(:,viSpk1);
    viClu2 = find(ismember(S_clu.viSite_clu, miSites_clu(:,iClu)));
    vrRho1 = vrRho_spk(viSpk1);
    for iDrift = 1:nDrift        
        viiSpk1 = find(ismember(viSpk1, cviSpk_drift{iDrift}));
        if numel(viiSpk1) >= P.min_count            
            switch 2
                case 2
                    miSpk_out11 = miKnn1(:,viiSpk1);
                    viSpk_out11 = miSpk_out11(vrRho_spk(miSpk_out11) >=  vrRho1(viiSpk1));
                case 1
                    viSpk_out11 = miKnn1(:,viiSpk1);
            end
            vr_ = cellfun(@(y)fh_member(viSpk_out11(:), y), S_clu.cviSpk_clu(viClu2));
            trDist_clu(viClu2,iDrift,iClu) = vr_;
        end
    end %for    
    fprintf('.');
end %for

%-----
% find the highest pair
mrDist_clu = zeros(nClu);
for iClu = 1:nClu
    mr1 = trDist_clu(:,:,iClu);    
    vr1 = mr1(iClu,:);        
    switch 2
        case 2
            mr2 = mr1 ./ vr1;
            mrDist_clu(:,iClu) = max(mr2,[],2, 'omitnan');
        case 1
            mr1(iClu,:) = nan;
            [p_out1, i1] = max(mr1(:), [], 'omitnan');
            if isnan(p_out1), continue; end
            [iClu1, iDrift1] = ind2sub(size(mr1), i1);
            p_in1 = vr1(iDrift1);
            if p_in1==0, continue; end
            mrDist_clu(iClu1,iClu) = p_out1 / p_in1;
    end %switch
end


%-----
% Merge
mlWavCor_clu = mrDist_clu >= merge_thresh;
viMap_clu = int32(ml2map_(mlWavCor_clu));
vlPos = S_clu.viClu > 0;
S_clu.viClu(vlPos) = viMap_clu(S_clu.viClu(vlPos)); %translate cluster number
S_clu = S_clu_refresh_(S_clu);
nClu_post = S_clu.nClu;
nClu_pre = nClu;
fprintf('\nMerged %d waveforms (%d->%d), took %0.1fs\n', nClu-nClu_post, nClu, nClu_post, toc(t_merge));
end %func


%--------------------------------------------------------------------------
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ml2map_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refresh_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
