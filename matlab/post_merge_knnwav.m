%--------------------------------------------------------------------------
% old templateMatch code
function S_clu = post_merge_knnwav(S_clu, viSite_spk, P)

KNN = 8;
fUse_raw = 0;
nAve_knn = min(KNN, P.knn);

fprintf('Automated merging (post-hoc)\n'); t1=tic;
[viClu_spk, miKnn, vrRho_spk] = deal(S_clu.viClu, S_clu.miKnn, S_clu.rho);
miKnn = miKnn(1:nAve_knn,:);
% viSite_spk = get0_('viSite_spk');
% frac_thresh = P.thresh_core_knn;
nShift_max = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
viShift = -nShift_max:nShift_max;

tnWav_spk = get_spkwav_(P, fUse_raw); % use raw waveform

% create template (nTemplate per cluster)
[cviSpk_in_clu, cviSpk_out_clu, ctrWav_in_clu, cviSite_in_clu] = deal(cell(S_clu.nClu, 1));
fh_car = @(tr)tr - repmat(mean(tr,2), [1,size(tr,2),1]);
fh_wav = @(vi)single(tnWav_spk(:,:,vi));
switch 1
    case 3 % knn smoothed waveform returned
        fh_trimmean = @(vi)tr2mr_trimmean_(fh_wav(vi));
        fh_mr = @(vi)tr2mr_mean_knn_(tnWav_spk, miKnn, viSite_spk, vi); 
    case 2 % spatial ref
        fh_trimmean = @(vi)tr2mr_trimmean_(fh_car(fh_wav(vi)));
        fh_mr = @(vi)reshape(fh_car(fh_wav(vi)), [], numel(vi));        
    case 1 
        fh_trimmean = @(vi)tr2mr_trimmean_(fh_wav(vi));
        fh_mr = @(vi)single(reshape(tnWav_spk(:,:,vi), [], numel(vi)));
end
fh_med = @(vi)single(median(tnWav_spk(:,:,vi),3));
fh_mean = @(vi)single(mean(tnWav_spk(:,:,vi),3));
fh_pv1 = @(vi)tr_pv1_(single(tnWav_spk(:,:,vi)));
fh_meanalign = @(vi)tr_mean_align_(single(tnWav_spk(:,:,vi)));
fh_denoise = @(vi)tr2mr_denoise_(single(tnWav_spk(:,:,vi)));
nSites = max(viSite_spk);
switch 1
    case 1, nDrift = get_set_(P, 'nTime_drift', 64);
    case 2, nDrift = get_set_(P, 'nTime_clu', 1);
end
nSpk_min = get_set_(P, 'knn', 30);
fprintf('\tComputing template\n\t'); t_template = tic;


for iClu = 1:S_clu.nClu
    viSpk1 = find(S_clu.viClu == iClu);
    viSpk1 = viSpk1(:);
    viiSpk1 = round(linspace(1, numel(viSpk1), nDrift+1));
    [vlKeep_clu1, viSite_clu1] = deal(true(nDrift, 1), zeros(nDrift,1));
    trWav_clu1 = zeros(size(tnWav_spk,1), size(tnWav_spk,2), nDrift, 'single');
    [miKnn1, vrRho1] = deal(miKnn(:,viSpk1), vrRho_spk(viSpk1)');
    for iDrift = 1:nDrift
        vii1 = viiSpk1(iDrift):viiSpk1(iDrift+1);
        [viSpk11, vrRho11, miKnn11] = deal(viSpk1(vii1), vrRho1(vii1), miKnn1(:,vii1));
        
        switch 3 %3 % expand selection using miKnn
            case 5 % only connect to neighbors with higher density
                viSpk11 = miKnn11(vrRho_spk(miKnn11) >= vrRho11);
                iSite11 = mode(viSite_spk(viSpk11));
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11);   
            case 4
                viSpk11 = viSpk11(vrRho11 > median(vrRho11));
                viSpk11 = miKnn(:, viSpk11);
                viSpk11 = viSpk11(:); % density based selection                
                iSite11 = mode(viSite_spk(viSpk11));
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11);   
            case 3     
                viSpk11 = miKnn(:, viSpk11);
                viSpk11 = viSpk11(:); % density based selection                
                iSite11 = mode(viSite_spk(viSpk11));
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11);   
            case 2                
                viSpk11 = miKnn(:, viSpk11);
                viSpk11 = viSpk11(:);
                viSpk11 = viSpk11(vrRho11 > quantile(vrRho11, frac_thresh));
                viSpk11 = viSpk11(viClu_spk(viSpk11) == iClu);
                iSite11 = mode(viSite_spk(viSpk11));
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11);
            case 1 % find local peak
                vrRho1 = S_clu.rho(viSpk1);
                [~,ii11] = max(vrRho1(vii1)); % local peak
                iSpk11 = viSpk1(vii1(ii11)); % center spike index
                iSite11 = viSite_spk(iSpk11);
                % find neighbors around the local peak
                viSpk11 = miKnn(:, miKnn(:,iSpk11));
                viSpk11 = viSpk11(:);
                viSpk11 = viSpk11(viSite_spk(viSpk11) == iSite11 & viClu_spk(viSpk11) == iClu);     
        end
        viSite_clu1(iDrift) = iSite11;
        if numel(viSpk11) < nSpk_min
            vlKeep_clu1(iDrift) = false;
            continue;
        end             
        switch 1
            case 15
                tr_ = single(tnWav_spk(:,:,viSpk11));
                mr_ = mean(tr_,3);
                mr_ = mr_(:,1) * (mr_(:,1)' * mr_); %  projection to the largest chan
            case 14
                tr_ = single(tnWav_spk(:,:,viSpk11));
                [~,mr_] = pca(reshape(tr_,[],size(tr_,3)), 'NumComponents',1);
                mr_ = reshape(mr_, size(tr_,1),[]);
            case 13
                tr_ = single(tnWav_spk(:,:,viSpk11));
                vr_ = sum(sum((tr_ - median(tr_,3)).^2),2);
                [~, vi_] = sort(vr_(:), 'ascend');
                mr_ = mean(tr_(:,:,vi_(1:end/2)),3);
            case 12
                tr_ = single(tnWav_spk(:,:,viSpk11));
                mr_ = trimmean(tr_,50,3);
            case 11
                tr_ = single(tnWav_spk(:,:,viSpk11));
                mr_ = median(tr_,3);
            case 10
                tr_ = single(tnWav_spk(:,:,viSpk11));
                mr_ = median(tr_,3) ./ std(tr_,1,3);     
            case 9
                tr_ = single(tnWav_spk(:,:,viSpk11));
                mr1_ = median(tr_,3);
                mr2_ = median(abs(tr_ - mr1_), 3);
                mr_ = mr1_ ./ mr2_;
            case 8 % weigh means by sd
                tr_ = single(tnWav_spk(:,:,viSpk11));
                mr_ = mean(tr_,3) ./ std(tr_,1,3);
            case 7
                tr_ = single(tnWav_spk(:,:,viSpk11));
                mr_ = mean(tr_,3) ./ std(tr_,1,3);
            case 6, mr_ = meanSubt_(single(median(tnWav_spk(:,:,viSpk11),3)));
            case 5, mr_ = meanSubt_(mean(single(tnWav_spk(:,:,viSpk11)),3));
            case 4, mr_ = single(median(tnWav_spk(:,:,viSpk11),3));
            case 3, mr_ = fh_pv1(viSpk11);
            case 2, mr_ = fh_trimmean(viSpk11);
            case 1, mr_ = mean(single(tnWav_spk(:,:,viSpk11)),3);
        end
        trWav_clu1(:,:,iDrift) = mr_;
    end
    ctrWav_in_clu{iClu} = trWav_clu1(:,:,vlKeep_clu1);
    cviSite_in_clu{iClu} = viSite_clu1(vlKeep_clu1);
    fprintf('.');
end
fprintf('\n\ttook %0.1fs\n', toc(t_template));


% merge the templates: todo, faster code
nClu = S_clu.nClu;
fprintf('Merging templates\n\t'); t_merge=tic;
fh_norm = @(x)bsxfun(@rdivide, x, std(x,1)*sqrt(size(x,1)));
switch 2
    case 2, fh_norm_tr = @(x)fh_norm(reshape(x, [], size(x,3)));
    case 1, fh_norm_tr = @(x)fh_norm(reshape(meanSubt_(x), [], size(x,3)));
end
mrDist_clu = nan(S_clu.nClu, 'single');
for iClu1 = 1:S_clu.nClu
    viSite_clu1 = cviSite_in_clu{iClu1};
    if isempty(viSite_clu1), continue; end
    mr1_ = fh_norm_tr(ctrWav_in_clu{iClu1});
    if 1 %shift template in time
        mr1_ = fh_norm_tr(shift_trWav_(ctrWav_in_clu{iClu1}, viShift));
        viSite_clu1 = repmat(viSite_clu1(:), numel(viShift), 1);
    end
    vrDist_clu1 = zeros(S_clu.nClu, 1, 'single');
    for iClu2 = iClu1+1:S_clu.nClu
        viSite2 = cviSite_in_clu{iClu2};
        viSite12 = intersect(viSite_clu1, viSite2);
        if isempty(viSite12), continue; end
        mr2_ = fh_norm_tr(ctrWav_in_clu{iClu2});        
        for iSite12_ = 1:numel(viSite12)
            iSite12 = viSite12(iSite12_);
            mrDist12 = mr2_(:, viSite2==iSite12)' * mr1_(:, viSite_clu1==iSite12);
            vrDist_clu1(iClu2) = max(vrDist_clu1(iClu2), max(mrDist12(:)));
        end
    end
    mrDist_clu(:, iClu1) = vrDist_clu1;
%     fprintf('.');
end %for
mlWavCor_clu = mrDist_clu >= P.maxWavCor;
viMap_clu = int32(ml2map_(mlWavCor_clu));
vlPos = S_clu.viClu > 0;
S_clu.viClu(vlPos) = viMap_clu(S_clu.viClu(vlPos)); %translate cluster number
S_clu = S_clu_refresh_(S_clu);
nClu_post = S_clu.nClu;
nClu_pre = nClu;
fprintf('\nMerged %d waveforms (%d->%d), took %0.1fs\n', nClu-nClu_post, nClu, nClu_post, toc(t_merge));
end %func


%==========================================================================
% call irc.m

function out1 = get_spkwav_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = shift_trWav_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ml2map_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refresh_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end