function mrDist_clu = calc_dist_drift(S0, S_clu, P, tnWav_spk)

fprintf('calc_dist_drift\n\t'); t1=tic;
MIN_COUNT = P.min_count;
NUM_PC = 3;
NUM_KNN = 15;

mrDist_clu = zeros(S_clu.nClu);
viSite_spk = S0.viSite_spk;
[viClu_spk, nClu, S_drift] = deal(S_clu.viClu, S_clu.nClu, S_clu.S_drift);
% S_drift: {miSort_drift, cviSpk_drift, nTime_drift, viDrift_spk, mlDrift}
nSites = max(viSite_spk);
nDrift = numel(S_drift.cviSpk_drift);
dimm_wav = size(tnWav_spk);
[miKnn_spk, vrRho_spk] = deal(S_clu.miKnn(1:NUM_KNN,:), S_clu.rho);

% find site range for each peak channels
switch 2
    case 2, nSites_fet = dimm_wav(2);
    case 1, nSites_fet = S0.dimm_fet(1) / P.nPcPerChan;
end
fh_norm = @(x)(x-mean(x)) ./ std(x,1);
fh_corr = @(x,y)(x'*y)/size(x,1);

if true
    [cviSite_nn_site, cmiSite_nn_site] = deal(cell(nSites, 1));
    for iSite = 1:nSites
        viSite1 = P.miSites(1:nSites_fet,iSite);
        [a_, b_] = ismember(P.miSites, viSite1);
        viSite_nn1 = find(sum(a_>0) == nSites_fet);
        cviSite_nn_site{iSite} = viSite_nn1;

        b_ = b_(:,viSite_nn1);
        [~, c_] = sort(b_);
        c_ = c_(end-nSites_fet+1:end,:);  
        cmiSite_nn_site{iSite} = c_;
    end %for
end

nSamples11 = (dimm_wav(1) - 2) * nSites_fet; % single sample shift
% nSamples11 = (dimm_wav(1) - 4) * nSites_fet; % two-sample shift
for iDrift = 1:nDrift
    viDrift1 = find(S_drift.mlDrift(:,iDrift));
    viSpk1 = cell2mat(S_drift.cviSpk_drift(viDrift1));    
    viSpk1 = viSpk1(viClu_spk(viSpk1)>0);    
    [viClu1, viSite1] = deal(viClu_spk(viSpk1), viSite_spk(viSpk1));

    switch 1
        case 2
            % find clu in this snapshot
            [vnClu_unique1, viClu_unique1] = unique_count_(viClu1);
            viClu_unique1 = viClu_unique1(vnClu_unique1 >= MIN_COUNT);
            nClu1 = numel(viClu_unique1);
            trWav_clu1 = zeros(dimm_wav(1), dimm_wav(2), nClu1, 'single');
            viSite_clu1 = zeros(nClu1, 1);
            for iiClu1 = 1:nClu1
                iClu11 = viClu_unique1(iiClu1);
                viiSpk11 = find(viClu1 == iClu11);
                [viSpk11, viSite11] = deal(viSpk1(viiSpk11), viSite1(viiSpk11));
                iSite11 = mode(viSite11);                
                viSpk12 = viSpk11(viSite11 == iSite11);
                if numel(viSpk12) < MIN_COUNT, continue; end
                viSite_clu1(iiClu1) = iSite11;
                trWav_clu1(:,:,iiClu1) = svd_mean_(tnWav_spk(:,:,viSpk12), NUM_PC); 
            end %for
            % remove empty clu
            vii_ = find(viSite_clu1>0);
            [viSite_clu1, trWav_clu1, viClu_clu1] = deal(viSite_clu1(vii_), trWav_clu1(:,:,vii_), viClu_unique1(vii_));
            mrDist1 = calc_corr_(viSite_clu1, trWav_clu1, viClu_clu1, cviSite_nn_site, cmiSite_nn_site);
            mrDist_clu(viClu_clu1,viClu_clu1) = max(mrDist_clu(viClu_clu1,viClu_clu1), mrDist1);
            
        case 1
            for iSite = 1:nSites        
                viiSpk11 = find(viSite1 == iSite);
                [viSpk11, viClu11] = deal(viSpk1(viiSpk11), viClu1(viiSpk11));
                [vnUniq_, viUniq_] = unique_count_(viClu11);
                viClu_uniq11 = viUniq_(vnUniq_ >= MIN_COUNT);
                nClu11 = numel(viClu_uniq11);
                if nClu11 < 2, continue; end
                % collect average waveforms (svd cleaned)
%                 mr11 = zeros(nSamples11*nSites_fet, nClu11, 'single'); % stores waveforms to correlate
                [a12, b12, c12] = deal(zeros(nSamples11, nClu11, 'single'));
                for iiClu11 = 1:nClu11
                    iClu = viClu_uniq11(iiClu11);
                    viSpk12 = viSpk11(viClu11==iClu);
                    switch 1
                        case 2
                            viSpk12 = miKnn_spk(:,viSpk12);
                            viSpk12 = unique(viSpk12(:));
                            viSpk12 = viSpk12(viSite_spk(viSpk12) == iSite);
                        case 1 % expand selection
                            miKnn12 = miKnn_spk(:,viSpk12);
                            viSpk12 = miKnn12(vrRho_spk(miKnn12) >= vrRho_spk(viSpk12)');
                            viSpk12 = unique(viSpk12(:));
                            viSpk12 = viSpk12(viSite_spk(viSpk12) == iSite);
                    end
                    mr_ = svd_mean_(tnWav_spk(:,:,viSpk12), NUM_PC);
                    mr_ = mr_(:, 1:nSites_fet);
                    [a_,b_,c_] = deal(mr_(1:end-2,:), mr_(2:end-1,:), mr_(3:end,:)); % single sample shift
%                     [a_,b_,c_] = deal(mr_(1:end-4,:), mr_(3:end-2,:), mr_(5:end,:)); % two-sample shift
                    [a12(:,iiClu11), b12(:,iiClu11), c12(:,iiClu11)] = deal(a_(:), b_(:), c_(:));
                end
                % compute waveform correlation
                [a12, b12, c12] = deal(fh_norm(a12), fh_norm(b12), fh_norm(c12));
                mrC1 = max(max(fh_corr(b12,b12), fh_corr(b12,a12)), fh_corr(b12,c12));
                mrDist_clu(viClu_uniq11,viClu_uniq11) = max(mrDist_clu(viClu_uniq11,viClu_uniq11), mrC1);
            end
    end
    fprintf('.');
end %for
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function mrDist1 = calc_corr_(viSite_clu1, trWav_clu1, viClu_clu1, cviSite_nn_site, cmiSite_nn_site)

dimm1 = size(trWav_clu1);
nClu1 = numel(viSite_clu1);
mrDist1 = zeros(nClu1);
viSite_uniq = unique(viSite_clu1);
nSites1 = numel(viSite_uniq);
fh_norm = @(x)(x-mean(x)) ./ std(x,1);
fh_corr = @(x,y)(x'*y)/size(x,1);
nSamples12 = (dimm1(1)-2) * size(cmiSite_nn_site{1},1);
for iiSite = 1:nSites1
    iSite1 = viSite_clu1(iiSite);
    viSite_nn1 = cviSite_nn_site{iSite1};
    miSite_nn1 = cmiSite_nn_site{iSite1}; 
    
    [viClu12, viiSite12] = ismember(viSite_clu1, viSite_nn1);
    [viClu12, viiSite12] = deal(find(viClu12), viiSite12(viClu12));
    nClu12 = numel(viClu12);
    if nClu12 < 2, continue; end
    [a12, b12, c12] = deal(zeros(nSamples12, nClu12, 'single'));
    for iiClu12 = 1:nClu12
        viSite12 = miSite_nn1(:, viiSite12(iiClu12));
        mr_ = trWav_clu1(:, viSite12, viClu12(iiClu12));
        [a_,b_,c_] = deal(mr_(1:end-2,:), mr_(2:end-1,:), mr_(3:end,:));
        [a12(:,iiClu12), b12(:,iiClu12), c12(:,iiClu12)] = deal(a_(:), b_(:), c_(:));
    end
    [a12, b12, c12] = deal(fh_norm(a12), fh_norm(b12), fh_norm(c12));
    
    % calc correlation
    mrC1 = max(max(fh_corr(b12,b12), fh_corr(b12,a12)), fh_corr(b12,c12));
    mrDist1(viClu12,viClu12) = max(mrDist1(viClu12,viClu12), mrC1);
end

end %func

    
%     switch 1
%         case 2
%             [viSite_nn11, miSite_nn11] = deal(cviSite_nn_site{iSite}, cmiSite_nn_site{iSite});
%             nSites11 = numel(viSite_nn11);
%             cmrWav_site11 = cell(nClu, nSites11);
%             for iiSite11 = 1:nSites11
%                 iSite11 = viSite_nn11(iiSite11);
%                 viiSpk11 = find(viSite1 == iSite11);
%                 [viSpk11, viClu11] = deal(viSpk1(viiSpk11), viClu1(viiSpk11));                    
%                 [cvrWav_clu11, viClu_uniq11] = get_wav_clu_(tnWav_spk, viSpk11, viClu11, miSite_nn11(:,iiSite11), MIN_COUNT);
%                 cmrWav_site11(viClu_uniq11, iiSite11) = cvrWav_clu11;
%             end
%             mrC_ = corr_max_(cmrWav_site11); % apply correlation and take max
%             mrDist_clu(viClu11,viClu11) = max(mrDist_clu(viClu11,viClu11), mrC_);
% 
%         case 1


%--------------------------------------------------------------------------
function [cvrWav_clu11, viClu_uniq11] = get_wav_clu_(tnWav_spk, viSpk11, viClu11, viiSite_nn11, MIN_COUNT)

dimm_wav = size(tnWav_spk);
[vnUniq_, viUniq_] = unique_count_(viClu11);
viClu_uniq11 = viUniq_(vnUniq_ >= MIN_COUNT);
nClu11 = numel(viClu_uniq11);

% collect average waveforms (svd cleaned)
nSites_fet = numel(viiSite_nn11);
% cmrWav_clu11 = zeros(dimm_wav(1)*nSites_fet, nClu11, 'single'); % stores waveforms to correlate
cvrWav_clu11 = cell(nClu11, 1);
for iiClu11 = 1:nClu11
    iClu = viClu_uniq11(iiClu11);
    mr_ = svd_mean_(tnWav_spk(:,:,viSpk11(viClu11==iClu)));
    mr_ = mr_(:,viiSite_nn11);
    vr_ = mr_(:);
    cvrWav_clu11{iiClu11} = (vr_ - mean(vr_)) / std(vr_);
    %cmrWav_clu11(:,iiClu11) = mr_(:);
end

% normalize waveform
% cmrWav_clu11 = (cmrWav_clu11 - mean(cmrWav_clu11)) ./ std(cmrWav_clu11,1);
end %func


% %--------------------------------------------------------------------------
% function mrC_ = corr_max_(cmrWav_clu11)
% 
% end %func


%--------------------------------------------------------------------------
function [vnUnique, viUnique] = unique_count_(vn)
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
% 8/21/2019 JJJ: Avegage waveforms using svd
function [mr2, mrPv] = svd_mean_(tr1, NUM_PC)
if nargin<2, NUM_PC = 3; end

dimm1 = size(tr1);
switch 1
    case 2
        tr2 = reshape(meanSubt_(reshape(tr1, [], size(tr1,3))), dimm1);
        mr1 = single(reshape(tr2,dimm1(1),[]));
    case 1
        mr1 = single(reshape(tr1,dimm1(1),[]));
end
mrC = mr1 * mr1';
[U,S,V] = svd(mrC);
mrPv = U(:,1:NUM_PC);
mrPc = reshape(mrPv' * mr1, NUM_PC, dimm1(2), dimm1(3));
mr2 = mrPv * mean(mrPc,3);
end %func