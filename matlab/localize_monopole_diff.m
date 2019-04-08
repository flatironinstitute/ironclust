function mrPos_spk = localize_monopole_diff(P, viSite_spk, tnWav_raw, nSites_fit, exponent)
if nargin<4
    exponent = .5;
end
ndims_site = 2;
pad_dimm = 12;
% maxDist_site_um = 50; %P.maxDist_site_um;

% cviSpk_site = S0.cviSpk_site;
miSites = P.miSites;
mrSiteXY = single(P.mrSiteXY);

nSpk = size(tnWav_raw,3);
nSites = max(viSite_spk);
mrPos_spk = zeros(3, nSpk, 'single');

% create a grid
% mrDist_site = pdist2(mrSiteXY, mrSiteXY);
% nSites_fit = median(sum(mrDist_site < maxDist_site_um));
nSites_fit = min(nSites_fit, size(tnWav_raw,2));

viSites_fit = 1:nSites_fit;
fh_vec = @(x)x(:);
fh_normalize = @(x)x./sqrt(sum(x.^2));

x = -20:.5:20;
y = -20:.5:20;
z = 1;
[xx,yy,zz] = meshgrid(x,y,z);
mrPos_fit0 = [xx(:), yy(:), zz(:)];

for iSite = 1:nSites
    viSpk1 = find(viSite_spk==iSite);
    if isempty(viSpk1), continue; end
    viSite1 = miSites(viSites_fit,iSite);
    nSpk1 = numel(viSpk1);
    trWav_spk1 = single(tnWav_raw(:,viSites_fit,viSpk1));
    mrFet1 = tr2slope_diff_(trWav_spk1, 3, [-3,3], 3);
%     mrFet1 = abs(mrFet1);
    mrPos_site1 = [mrSiteXY(viSite1,:), zeros(numel(viSite1),1)];  
    vrX0_spk1 = com_mr_(mrPos_site1(:,1), mrFet1, exponent); % initial guess
    vrY0_spk1 = com_mr_(mrPos_site1(:,2), mrFet1, exponent); % initial guess
    mrPos0_spk1 = single([vrX0_spk1(:), vrY0_spk1(:), zeros(nSpk1,1)]');
    mrPos_spk1 = zeros(3, nSpk1, 'single');
    vrFit_spk1 = zeros(nSpk1, 1, 'single');
    if 0 % use finite site geometry. not helping and slow
        [mrPos_site2, mrSum2] = expand_site_(mrPos_site1, pad_dimm, 6);
    else
        [mrPos_site2, mrSum2] = deal(mrPos_site1, []);
    end

    % Fitting loop   
    parfor iSpk1 = 1:nSpk1
%     parfor iSpk1 = 1:10000
        xyz0_ = mrPos0_spk1(:,iSpk1)';
        mrVd_fit1 = create_fit_(mrPos_site2, mrSum2, mrPos_fit0 + xyz0_); % return normalized fit
        vrV_ = mrFet1(:,iSpk1);
        vrVd_ = vrV_ - vrV_';
        [vrFit_spk1(iSpk1), imax_] = max(abs(fh_normalize(vrVd_(:))' * fh_normalize(mrVd_fit1))); % max dot product        
        mrPos_spk1(:,iSpk1) = mrPos_fit0(imax_,:) + xyz0_;
%         fprintf('fit:%f, z:%f\n', vrFit_spk1(iSpk1), mrPos_spk1(3,iSpk1));
    end
%     median(vrFit_spk1(1:100))
%     mr_ = mrPos0_spk1(:,1:100) - mrPos_spk1(:,1:100);
    mrPos_spk(:,viSpk1) = mrPos_spk1;
    fprintf('.');
end %for
                
                
end %func


%--------------------------------------------------------------------------
function mrVd_fit2 = create_fit_(mrPos_site1, mrSum1, mrPos_fit1)
% create a position grid and 
% [mrVd_fit1, mrPos_fit1] = deal([]);

mrVd_fit1 = 1 ./ pdist2(mrPos_site1, mrPos_fit1);
if ~isempty(mrSum1)
    mrVd_fit1 = mrSum1 * mrVd_fit1;
end
nSites = size(mrVd_fit1,1);
[viB, viA] = meshgrid(1:nSites); 
mrVd_fit2 = mrVd_fit1(viA(:),:) - mrVd_fit1(viB(:),:);
end %func


%--------------------------------------------------------------------------
% use channel difference
function mrFet = tr2slope_diff_(trWav, nDiff, nlim, nPc)
dimm = size(trWav);
[nSamples, nSites, nSpk] = deal(dimm(1), dimm(2), dimm(3));
% trWav1 = ndiff_(trWav, nDiff); % diff filter
% trWav1 = sgfilt_(ndiff_(trWav, 5), nDiff); % diff filter
% trWav1 = sgfilt_(ndiff_(trWav, 5), nDiff); % diff filter

switch 4
    case 4
%         no differentiation
        [~, imin1] = min(mean(mean(trWav,3),2)); 
        vi1 = max(1, imin1+nlim(1)):min(nSamples, imin1+nlim(end));
        trWav2 = pca_clean_(trWav(vi1,:,:), nPc);
        
    case 3
        nlim=nlim*1;
        dimm = size(trWav);
        trWav1 = zeros(dimm + [0,-1,0],'like', trWav);
        for iSpk = 1:dimm(3)
            mr_ = trWav(:,:,iSpk);
            trWav1(:,:,iSpk) = mr_(:,1) - mr_(:,2:end);
        end        
%         trWav1 = ndiff_(trWav1, nDiff);
        [~, imin1] = min(mean(mean(trWav1,3),2)); 
        vi1 = max(1, imin1+nlim(1)):min(nSamples, imin1+nlim(end));
        [trWav2, vrL] = pca_clean_(trWav1(vi1,:,:), nPc);
        
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

mrFet = reshape(permute(trWav2, [2,3,1]), size(trWav2,2), []);
mrFet = mrFet(:, sub2ind([nSpk, nSamples1], 1:nSpk, vimin_spk'));
end %func


%--------------------------------------------------------------------------
function [trWav1, vrLambda] = pca_clean_(trWav, nPc)
nChans = size(trWav,2);
vrLambda= [];
switch 2
    case 2
         % global channel pca
         [a,b,vrLambda] = pca(reshape(trWav,size(trWav,1),[])', 'Centered', 0, 'NumComponents', nPc);   
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



%--------------------------------------------------------------------------
% expand each electrode positions by four positions
function [mrPos1, mrSum1] = expand_site_(mrPos, pad_dimm, nPolygon)
% pad_dimm: [width, height] or [diameter]

if numel(pad_dimm) == 1
    diameter = pad_dimm;
else
    diameter = sqrt(prod(pad_dimm) / pi) * 2;
end
% build polygon
r = diameter/3;
vrTh = linspace(0, 2*pi, nPolygon+1)';
vrTh(end) = [];

switch 1
    case 2
        [vrX0, vrY0, vrZ0] = deal(r * cos(vrTh)+mrPos(1,1), r * sin(vrTh)+mrPos(1,2), zeros(nPolygon,1)+mrPos(1,3));
        mrPos1 = [[vrX0(:), vrY0(:), vrZ0(:)]; mrPos];
        mrSum1 = blkdiag(repmat(1/(nPolygon+1),1, nPolygon+1), eye(size(mrPos,1)-1));
    case 1
        [vrX0, vrY0, vrZ0] = deal([r * cos(vrTh); 0], [r * sin(vrTh); 0], zeros(nPolygon+1,1));
        [vrX1, vrY1, vrZ1] = deal(mrPos(:,1)' + vrX0, mrPos(:,2)' + vrY0, mrPos(:,3)' + vrZ0);
        mrPos1 = [vrX1(:), vrY1(:), vrZ1(:)];
        mrSum1 = kron(eye(size(mrPos,1)), repmat(1/(nPolygon+1),1, nPolygon+1));
end
%  disp(mrSum1*mrPos1 - mrPos)
end %func