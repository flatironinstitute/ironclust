function mrPos_spk = localize_moment(P, viSite_spk, tnWav_raw, nSites_fit, vrExp)
if nargin<4
    vrExp = [.5, 1, 2]; % exponents weighing factor
end
ndims_site = 2;
% maxDist_site_um = 50; %P.maxDist_site_um;

% cviSpk_site = S0.cviSpk_site;
miSites = P.miSites;
mrSiteXY = single(P.mrSiteXY);

nSpk = size(tnWav_raw,3);
nSites = max(viSite_spk);
nFet = ndims_site * numel(vrExp);
if numel(vrExp) == 1
    mrPos_spk = zeros(3, nSpk, 'single');
else
    mrPos_spk = zeros(nFet, nSpk, 'single');
end

% create a grid
% mrDist_site = pdist2(mrSiteXY, mrSiteXY);
% nSites_fit = median(sum(mrDist_site < maxDist_site_um));
nSites_fit = min(nSites_fit, size(tnWav_raw,2));

viSites_fit = 1:nSites_fit;
fh_vec = @(x)x(:);
for iSite = 1:nSites
    viSpk1 = find(viSite_spk==iSite);
    if isempty(viSpk1), continue; end
    viSite1 = miSites(viSites_fit,iSite);
    nSpk1 = numel(viSpk1);
    mrPos_spk1 = zeros(nFet, nSpk1, 'single');
    trWav_spk1 = single(tnWav_raw(:,viSites_fit,viSpk1));
    mrFet1 = tr2slope_(trWav_spk1, 3, [-3,3], 3);
    mrFet1 = abs(mrFet1);
    mrPos_site1 = mrSiteXY(viSite1,:);  
    if numel(vrExp) == 1
        iCase = 3;
    else
        iCase = 1;
    end
    switch iCase
        case 3
            % localize z using multiple expansion and site voltage difference
            mrPos_ = single(com_mr_(mrPos_site1, mrFet1, vrExp))';
            vrZ1 = zeros(nSpk1,1,'single');
            parfor iSpk1 = 1:nSpk1
                xy0_ = mrPos_(:,iSpk1)';
                vk_ = mrFet1(:,iSpk1);
                rk_ = pdist2(xy0_, mrPos_site1);
                rk2_ = rk_(:).^2; 
                vx2_ = fh_vec(rk2_-rk2_');
                vy_ = fh_vec(vk_(:) - vk_(:)');
%                 figure; plot(vx2_, vy_, '.');
                vrZ1(iSpk1) = abs((vy_) \ (vx2_));
            end
            mrPos_spk1 = [mrPos_; vrZ1(:)'];
        case 2
            for iExp = 1:numel(vrExp)
                exp1 = vrExp(iExp);
                viFet = (1:ndims_site) + (iExp-1) * ndims_site;
                mrPos_ = com_mr_(mrPos_site1.^exp1, mrFet1, 1) - mrPos_spk1(1:ndims_site,:)'.^exp1;
                mrPos_spk1(viFet,:) = (mrPos_.^(1/exp1))';
            end
        case 1
            for iExp = 1:numel(vrExp)
                exp1 = vrExp(iExp);
                viFet = (1:ndims_site) + (iExp-1) * ndims_site;
                mrPos_spk1(viFet,:) = com_mr_(mrPos_site1, mrFet1, exp1)';
            end
    end % switch
    mrPos_spk(:,viSpk1) = mrPos_spk1;
end %for
                
                
end %func


%--------------------------------------------------------------------------
function mrFet = tr2slope_(trWav, nDiff, nlim, nPc)
dimm = size(trWav);
[nSamples, nSites, nSpk] = deal(dimm(1), dimm(2), dimm(3));
% trWav1 = ndiff_(trWav, nDiff); % diff filter
% trWav1 = sgfilt_(ndiff_(trWav, 5), nDiff); % diff filter
% trWav1 = sgfilt_(ndiff_(trWav, 5), nDiff); % diff filter

switch 3
    case 3
%         no differentiation
        [~, imin1] = min(mean(mean(trWav,3),2)); 
        vi1 = max(1, imin1+nlim(1)):min(nSamples, imin1+nlim(end));
        trWav2 = pca_clean_(trWav(vi1,:,:), nPc);
    case 2
        trWav1 = ndiff_(trWav, nDiff);
        [~, imin1] = min(mean(mean(trWav1,3),2)); 
        vi1 = max(1, imin1+nlim(1)):min(nSamples, imin1+nlim(end));
        trWav2 = pca_clean_(trWav1(vi1,:,:), nPc);
    case 1        
        [~, imin1] = min(mean(mean(trWav,3),2)); 
        vi1 = max(1, imin1+nlim(1)):min(nSamples, imin1+nlim(end));
        trWav1 = trWav(vi1,:,:);
        trWav2 = ndiff_(pca_clean_(trWav1, nPc), nDiff); % diff filter
end %switch

% restrict search domain
nSamples1 = numel(vi1);

% denoise and slope

% trWav2 = pca_clean_(trWav2, 5);

% extract peak slope (summed)
[~, vimin_spk] = min(mean(trWav2,2));                          
vimin_spk = squeeze(vimin_spk);

mrFet = reshape(permute(trWav2, [2,3,1]), nSites, []);
mrFet = mrFet(:, sub2ind([nSpk, nSamples1], 1:nSpk, vimin_spk'));
end %func


%--------------------------------------------------------------------------
function trWav1 = pca_clean_(trWav, nPc)
nChans = size(trWav,2);
switch 2
    case 2
         % global channel pca
         [a,b,c] = pca(reshape(trWav,size(trWav,1),[])', 'Centered', 0, 'NumComponents', nPc);   
         trWav1 = reshape(a*b', size(trWav));         
    case 1        
        % individual channel pca
        trWav1 = permute(trWav, [1,3,2]);
        try
            parfor iChan = 1:nChans
                [a,b] = pca(trWav1(:,:,iChan)', 'Centered', 0, 'NumComponents', nPc);   
                trWav1(:,:,iChan) = a*b'; % denoised version
            end
        catch
            for iChan = 1:nChans
                [a,b] = pca(trWav1(:,:,iChan)', 'Centered', 0, 'NumComponents', nPc);   
                trWav1(:,:,iChan) = a*b'; % denoised version
            end
        end
        trWav1 = permute(trWav1, [1,3,2]);
end
end %func


%--------------------------------------------------------------------------
function mrPos_spk = com_mr_(mrPos, mrV, exponent)
if nargin<3, exponent = 1; end
if exponent == 1
    mrV2 = abs(mrV)';
else
    mrV2 = abs(mrV.^exponent)';
end
vrV2 = sum(mrV2,2);

if size(mrPos,2) > 1
    vrX = mrPos(:,1)';
    vrY = mrPos(:,2)';
    mrPos_spk = [sum(mrV2 .* vrX,2)./vrV2, sum(mrV2 .* vrY,2)./vrV2];
else
    vrX = mrPos(:,1)';
    mrPos_spk = sum(mrV2 .* vrX,2)./vrV2;
end
end %func


%--------------------------------------------------------------------------
% 9/3/2018 JJJ: nDiff filter generalized to any order (previously up to 4)
%   and uses less memory
function mn1 = ndiff_(mn, nDiff_filt)

if nDiff_filt==0, return; end
mn1 = zeros(size(mn), 'like', mn);
switch ndims(mn)
    case 2
        mn1(1:end-1,:) = diff(mn);
        for i = 2:nDiff_filt
            mn1(i:end-i,:) = 2*mn1(i:end-i,:) + (mn(i*2:end,:)-mn(1:end-i*2+1,:));
        end
    case 3
        mn1(1:end-1,:,:) = diff(mn);
        for i = 2:nDiff_filt
            mn1(i:end-i,:,:) = 2*mn1(i:end-i,:,:) + (mn(i*2:end,:,:)-mn(1:end-i*2+1,:,:));
        end
    otherwise
        error('ndiff_; unsupported dims');
end
end %func