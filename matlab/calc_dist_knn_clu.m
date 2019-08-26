function mrDist_clu = calc_dist_knn_clu(S_clu, P, run_mode)
if nargin<3, run_mode = []; end
if isempty(run_mode), run_mode = 14; end
    
DEG_SEPARATION = 2;
NUM_KNN = 16;
MIN_COUNT = P.min_count;

nClu = numel(setdiff(unique(S_clu.viClu), 0));
mrDist_clu = zeros(nClu);
nknn = min(size(S_clu.miKnn,1), NUM_KNN);
miKnn = int32(S_clu.miKnn);
miKnn1 = miKnn(1:nknn,:);
viSpk_nn_spk = int32(S_clu.nneigh);
viClu_spk = S_clu.viClu;
S_drift = S_clu.S_drift;

% deprecated
miKnn_clu = miKnn(:,S_clu.icl);
viSpk_clu = int32(S_clu.icl);
% miKnn_nn_clu = [miKnn_clu; viSpk_nn_spk(miKnn_clu)];
fprintf('calc_dist_knn_clu (run_mode=%d)\n\t', run_mode); 
t1 = tic;

switch run_mode % 11, 5, run_mode
    case 15        
        nDrift = numel(S_drift.cviSpk_drift);        
        for iDrift = 1:nDrift            
            viDrift1 = find(S_drift.mlDrift(:,iDrift));
            viSpk1 = cell2mat(S_drift.cviSpk_drift(viDrift1));    
            viClu1 = viClu_spk(viSpk1);
            [vnUniq_, viUniq_] = unique_count1_(viClu1);
            viClu_uniq1 = viUniq_(vnUniq_ >= MIN_COUNT & viUniq_ > 0);
            if isempty(viClu_uniq1), continue; end
            
            nClu1 = numel(viClu_uniq1);                        
            miKnn1 = miKnn(:,viSpk1);
            mnKnn_clu = zeros(nClu);
            for iiClu1 = 1:nClu1
                iClu11 = viClu_uniq1(iiClu1);
                viSpk11 = miKnn1(:,viClu1==iClu11);  
                viSpk11 = unique(viSpk11(:));
                [vn_,vi_] = unique_count1_(viClu_spk(viSpk11));
                
                mnKnn_clu(:, iiClu1) = histcounts(viClu_spk(viSpk_clu1), 1:nClu+1);
            end %for            
            [~,viMax] = max(mnKnn_clu);
            mrDist_clu(sub2ind([nClu,nClu], viMax, 1:nClu)) = 1;
            fprintf('.');
        end %for
        
    case 14
        mnKnn_clu = zeros(nClu);
        vnSpk_clu = zeros(nClu, 1);
        for iClu1 = 1:nClu
            vl_ = viClu_spk==iClu1;
            vnSpk_clu(iClu1) = sum(vl_);
            viSpk_clu1 = miKnn(:,vl_);
            viSpk_clu1 = unique(viSpk_clu1(:));
            mnKnn_clu(:, iClu1) = histcounts(viClu_spk(viSpk_clu1), 1:nClu+1);
        end
        [~,viMax] = max(mnKnn_clu);
        mrDist_clu = zeros(nClu);
        mrDist_clu(sub2ind([nClu,nClu], viMax, 1:nClu)) = 1;
        
    case 13
        mrDist_clu = zeros(nClu);
        nDrift = numel(S_drift.cviSpk_drift);        
        for iDrift = 1:nDrift
            viDrift1 = find(S_drift.mlDrift(:,iDrift));
            viSpk1 = cell2mat(S_drift.cviSpk_drift(viDrift1));    
            viClu1 = viClu_spk(viSpk1);
            [vnUniq_, viUniq_] = unique_count1_(viClu1);
            viClu_uniq1 = viUniq_(vnUniq_ >= MIN_COUNT & viUniq_ > 0);
            if isempty(viClu_uniq1), continue; end
            nClu1 = numel(viClu_uniq1);            
            miKnn1 = miKnn(:,viSpk1);
            for iiClu1 = 1:nClu1
                iClu11 = viClu_uniq1(iiClu1);
                viSpk11 = miKnn1(:,viClu1==iClu11);  
                viSpk11 = unique(viSpk11(:));
                [vn_,vi_] = unique_count1_(viClu_spk(viSpk11));
                vi = vi_(vn_ > vn_(vi_==iClu11) & vi_>0);    
                mrDist_clu(vi,iClu11) = 1;
            end %for
            fprintf('.');
        end %for
        
    case 12
        for iClu1 = 1:nClu
            viKnn1 = miKnn_clu(:,iClu1);
            mrDist_clu(:,iClu1) = sum(ismember(miKnn_clu, viKnn1));
        end 
    case 11
        mrDist_clu = zeros(nClu);
        for iClu1 = 1:nClu
            viSpk1 = miKnn(:,viClu_spk==iClu1); viSpk1 = unique(viSpk1(:));
%             viSpk1 = miKnn(:,viSpk1(:));  %viSpk1 = unique(viSpk1(:));
%             viSpk1 = miKnn1(:,viSpk1);  viSpk1 = (viSpk1(:));
            [vn_,vi_] = unique_count1_(viClu_spk(viSpk1));
            vi = vi_(vn_ > vn_(vi_==iClu1) & vi_>0);    
            mrDist_clu(vi,iClu1) = 1;
        end
%         mnKnn_clu = min(mnKnn_clu, mnKnn_clu');
    case 10 % diffusing population ratio
        mrDist_clu = zeros(nClu);
        for iClu1 = int32(1:nClu)
            viKnn1 = miKnn_clu(:,iClu1);
            for iDeg = 1:DEG_SEPARATION
                viKnn1 = miKnn(:,viKnn1);
                viKnn1 = (viKnn1(:));
            end
            switch 1
                case 2
                    viKnn1 = unique(viKnn1);
                    iClu2 = mode(viClu_spk(viKnn1));
                    if iClu2 ~= iClu1 && iClu2 > 0
                        mrDist_clu(iClu2,iClu1) = 1;
                    end
                case 1
                    [vn_,vi_] = unique_count1_(viClu_spk(viKnn1));
                    vi = vi_(vn_ > vn_(vi_==iClu1)/2);
                    vi = setdiff(vi,0);
                    mrDist_clu(vi,iClu1) = 1;
            end                        
        end
    case 9
        [vn, mi] = ismember(miKnn, viSpk_clu);
        vi = find(sum(vn) == 2);
        mi = sort(mi(:,vi)); mi = mi(end-1:end,:)';
        mrDist_clu = histcounts2(mi(:,1), mi(:,2), 1:nClu+1, 1:nClu+1);
    case 8
        for iClu1 = 1:nClu
            viKnn1 = miKnn_clu(:,iClu1);
            for iDeg = 1:DEG_SEPARATION   
                viKnn1 = miKnn1(:,viKnn1);
                viKnn1 = unique(viKnn1(:));
            end
            mrDist_clu(:,iClu1) = sum(ismember(miKnn_clu, viKnn1));
        end 
    case 7
        miKnn_clu = int32(S_clu.miKnn(:,S_clu.icl));
        for iClu1 = 1:nClu
            viKnn1 = miKnn_clu(:,iClu1);
            for iDeg = 2:DEG_SEPARATION   
                viKnn1 = miKnn1(:,viKnn1);
                viKnn1 = unique(viKnn1(:));
            end
            mrDist_clu(:,iClu1) = ismember(S_clu.icl, viKnn1);
        end   
    case 6
        miKnn_clu = int32(S_clu.miKnn(:,S_clu.icl));
        for iClu1 = 1:nClu
            viKnn1 = miKnn_clu(:,iClu1);   
            viKnn1 = miKnn1(:,viKnn1);
            mrDist_clu(:,iClu1) = sum(ismember(miKnn_clu, viKnn1(:)));
        end   
    case 5
        % find exact knn of the peaks using feature matrix
        [mnKnn_lower_clu, mnKnn_upper_clu] = deal(zeros(nClu));
        for iClu1 = 1:nClu
            vi_ = miKnn_clu(:,iClu1);            
            vi_ = miKnn(1:nknn,vi_);
            vi_ = sort(vi_(:));
            if iClu1 > 1
                mnKnn_lower_clu(1:iClu1-1,iClu1) = sum(ismember(miKnn_clu(:,1:iClu1-1), vi_))';
            end
            if iClu1 < nClu
                mnKnn_upper_clu(iClu1+1:end,iClu1) = sum(ismember(miKnn_clu(:,iClu1+1:end), vi_))';
            end
        end   
        mnKnn_lower_clu = mnKnn_lower_clu + mnKnn_lower_clu';
        mnKnn_upper_clu = mnKnn_upper_clu + mnKnn_upper_clu';
        mrDist_clu = min(mnKnn_lower_clu, mnKnn_upper_clu);
    case 4
        for iClu1 = 1:nClu
            vi_ = miKnn_clu(:,iClu1);            
            for iDegree = 1:(DEG_SEPARATION-1)
                vi_ = miKnn(1:nknn,vi_);
                vi_ = vi_(:);
            end
            mrDist_clu(:,iClu1) = mrDist_clu(:,iClu1) + sum(ismember(miKnn_clu, sort(vi_)))';
        end 
    case 3
        for iClu1 = 2:nClu
%             vi_ = miKnn_clu(1:nknn,iClu1);
            vi_ = miKnn_clu(:,iClu1);
            for iDegree = 1:(DEG_SEPARATION-1)
                vi_ = miKnn(1:nknn,vi_);
                vi_ = vi_(:);
            end
            mrDist_clu(1:iClu1-1,iClu1) = sum(ismember(miKnn_clu(:,1:iClu1-1), sort(vi_)));
        end 
    case 2
        for iClu1 = 2:nClu
            mrDist_clu(1:iClu1-1,iClu1) = sum(ismember(miKnn_clu(:,1:iClu1-1), miKnn_clu(:,iClu1)));
        end        
    case 1
        for iClu1 = 1:nClu
            viKnn1 = miKnn_clu(:,iClu1);
            for iClu2 = 1:iClu1-1        
                mrDist_clu(iClu2,iClu1) = sum(ismember(viKnn1, miKnn_clu(:,iClu2)));
            end
        end
end %switch
fprintf('\n\ttook %0.1fs\n', toc(t1));

if nargout>=2
    knn_merge_thresh = 1;
    [viMap, viUniq_] = ml2map_(mrDist_clu >= knn_merge_thresh);
    viMap = viMap(:);
    fprintf('S_clu_peak_merge_: %d->%d cluster centers (knn_merge_thresh=%d)\n', ...
        nClu, numel(viUniq_), knn_merge_thresh);
end
end %func


%--------------------------------------------------------------------------
function [vnUnique, viUnique] = unique_count1_(vn)
% vn = [1 1 1 1 4 4 4 8];
% vn = [1 1 1 1];
vn = sort(vn(:)');
vi_ = find(diff(sort(vn)) > 0);
if isempty(vi_)
    vnUnique = 1;
    viUnique = vn(1);
else
    vnUnique = diff([0, vi_, numel(vn)]);
    viUnique = [vn(1), vn(vi_+1)];
end
end %func


%--------------------------------------------------------------------------
function [out1, out2] = ml2map_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
