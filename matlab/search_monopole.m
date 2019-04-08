% goal: minimize cosine distance 
% apply fminbnd from the max position
% return strength
function [mrPos_spk, vrSource_spk, vrErr_spk] = search_monopole(mrPos_site, mrObs, xyz0, MaxIter)
% try excluding the nearest electrode
% if nargin<4, iFet = 1; end
% if nargin<5, MaxIter = 30; end

% thresh_error = .05;
thresh_error = [];

% problem based optimization


% initial position
% fh_cos = @(x,y)1-sum(x.*y); %assume normalized
fh_normalize = @(x)x./sqrt(sum(x.^2));
fh_norm = @(x)sqrt(sum(x.^2));
% fh_forward = @(xyz)fh_normalize(1./pdist2(xyz, mrPos_site));
% switch 3
%     case 3, fh_err = @(xyz, obs)fh_cos(fh_forward(xyz), obs);
%     case 2, fh_err = @(xyz, obs)pdist2(fh_forward(xyz), obs, 'correlation');
%     case 1, fh_err = @(xyz, obs)pdist2(fh_forward(xyz), obs, 'cosine');
% end
mrObs = double(abs(mrObs));
mrObs = mrObs ./ sqrt(sum(mrObs.^2)); % normalize output
mrObs_inv_norm = fh_normalize(1./mrObs);
mrObs_sq_inv_norm = fh_normalize(1./mrObs.^2);

% xyz0 = [mrPos_site(1,1), mrPos_site(1,2), 1];
nSpk = size(mrObs,2);

%%%%%
% use matrix fitting using lsqlin...


%%%%%
% use matrix inversion
% xyz0 = mrPos_site(1,:);
% dist_min = pdist(mrPos_site); 
% dist_min = min(dist_min(dist_min>0));
% vrX_field = xyz0(1) + linspace(dist_min*-2, dist_min*2, 64);
% vrY_field = xyz0(2) + linspace(dist_min*-2, dist_min*2, 64);
% vrZ_field = xyz0(3) + linspace(.01, dist_min*2, 32);
% [xx,yy,zz] = meshgrid(vrX_field, vrY_field, vrZ_field);
% mrPos_field = [xx(:), yy(:), zz(:)];
% mrB = pinv(1./pdist2(mrPos_site, mrPos_field));
% [~, viMax] = max(mrB * mrObs);
% [vx,vy,vz] = deal(mrPos_field(:,1), mrPos_field(:,2), mrPos_field(:,3));
% mrPos_spk = [vx(viMax), vy(viMax), vz(viMax)];
% vrErr_spk = arrayfun(@(i)fh_err(mrPos_spk(i,:), mrObs(:,i)'), 1:nSpk);

%%%%%
% use fitting 
% tic
mrPos_spk = zeros(3, nSpk);
[vrErr_spk, vrErr0_spk] = deal(zeros(nSpk,1));
if nargin<4
    S_opt = []; 
else
    S_opt = struct('Display', 'off', 'MaxIter', MaxIter, 'TolFun', .01);
end
if 1 % finite electrode size
    [mrPos1_site, mrSum1_site] = expand_site_(mrPos_site, [12, 12], 6);
    [mrPos1_site, mrSum1_site] = deal(single(mrPos1_site), single(mrSum1_site));
else
    [mrPos1_site, mrSum1_site] = deal([]);
end
switch 11 % 
    case 12 % project to 3D space
        mrPos_spk = zeros(3, nSpk, 'single');
        vrSource_spk = zeros(nSpk,1);
        [nx,ny,nz, nAve] = deal(16, 32, 64, 128);
        xlim = [-32,32] + xyz0(1);
        ylim = [-64,64] + xyz0(2);
        x = single(linspace(xlim(1),xlim(2),nx));
        y = single(linspace(ylim(1),ylim(2),ny));
        z = single(logspace(-1,2,nz));r
        [xx,yy,zz] = meshgrid(x, y, z);
        mrXYZ = [xx(:), yy(:), zz(:)]; 
        mrPos_site = single(mrPos_site);
        mrObs = single(mrObs);
        mrD = 1./pdist2(mrPos_site, mrXYZ);
        parfor iSpk = 1:nSpk
            vrMag = sum(mrD .* mrObs(:,iSpk));
            [vrSrt_, viSrt_] = sort(vrMag, 'descend');
            [vrSrt_, viSrt_] = deal(vrSrt_(1:nAve), viSrt_(1:nAve));            
            [mrPos_spk(:,iSpk), vrSource_spk(iSpk)] = com_mr_(mrXYZ(viSrt_,:), vrSrt_);
        end %for
        mrPos_spk=mrPos_spk';
        return;
        
    case 11 % recursive search
        % find min electrode distance and use that to determine the search
        nAve = 8;
        nDepth = 3;
        vcDist = 'euclidean'; % 'squaredeuclidean' for dipole
        vrD = pdist(mrPos_site);
        site_pitch = min(vrD(vrD>0));
        [nx, ny, nz] = deal(8, 16, 16);
        [xlim, ylim, zlim] = deal([-2,2]*site_pitch, [-2,2]*site_pitch, [-2,2]*site_pitch+eps());
        mrPos_spk = single(repmat(xyz0(:), 1, nSpk));
        mrObs_n = fh_normalize(mrObs);
        mrPos_site = single(mrPos_site);
        for iDepth = 1:nDepth
            factor = (1/2)^(iDepth-1);          
            [xlim1,ylim1,zlim1] = deal(xlim*factor, ylim*factor, zlim*factor);
            x = single(linspace(xlim1(1),xlim1(2),nx));
            y = single(linspace(ylim1(1),ylim1(2),ny));
            z = single(linspace(zlim1(1),zlim1(2),nz));
            if iDepth==1
                z = single(logspace(-1, log10(zlim(2)), nz));
            else
                z = single(linspace(zlim1(1),zlim1(2),nz));
            end
            [xx,yy,zz] = meshgrid(x, y, z);
            mrXYZ1 = [xx(:), yy(:), zz(:)];
%             tic
%             [mrPos_site, mrXYZ1, mrObs_n, mrPos_spk] = deal(gpuArray(mrPos_site), gpuArray(mrXYZ1), gpuArray(mrObs_n), gpuArray(mrPos_spk));
            parfor iSpk=1:nSpk
                vrPos_spk = mrPos_spk(:,iSpk)';                
                if isempty(mrPos1_site)                    
                    mrV_ = 1 ./ pdist2(mrPos_site, mrXYZ1+vrPos_spk, vcDist);
                else                    
                    mrV_ = mrSum1_site * (1 ./ pdist2(mrPos1_site, mrXYZ1+vrPos_spk, vcDist));
                end
%                 vrFit_ = fh_normalize(mrV_)' * mrObs_n(:,iSpk);
                vrFit_ = mrV_' * mrObs_n(:,iSpk) ./ sqrt(sum(mrV_.^2)');
                if iDepth < nDepth
                    [~, imax_] = max(vrFit_);
                    vrPos_ = mrXYZ1(imax_,:);
                else
                    if nAve>1
                        [vrSrt_, viSrt_] = sort(vrFit_, 'descend');
                        [vrSrt_, viSrt_] = deal(vrSrt_(1:nAve), viSrt_(1:nAve));
                        vrPos_ = mean(mrXYZ1(viSrt_,:));
                    else
                        [~,imax_] = max(vrFit_);
                        vrPos_ = mrXYZ1(imax_,:);
                    end                    
                end
                mrPos_spk(:,iSpk) = vrPos_ + vrPos_spk;
            end
%             toc
            mrPos_spk(3,:) = abs(mrPos_spk(3,:)); % make z positive
        end
        mrPos_spk = (mrPos_spk);
        mrD = pdist2(mrPos_site, mrPos_spk');
        vrSource_spk = mean(mrObs .* mrD)';        
        mrPos_spk = mrPos_spk';
        return;    
    case 10 % full pairwise correlation
        nx = 32;
        [xx,yy,zz] = meshgrid(xyz0(1) + linspace(-32,32,nx), xyz0(2) + linspace(-32,32,nx), xyz0(3) + linspace(0.1,32,nx));
        mrXYZ = [xx(:), yy(:), zz(:)];
        mrD = pdist2(mrXYZ, mrPos_site); mrObs_n = fh_normalize(mrObs);
        mrV_fit = mr2mr_mtimes_(1./mrD, [], 2);
        [vrError, viError] = min(1 - fh_normalize(mrV_fit')'*mrObs_n);
        mrPos_spk = mrXYZ(viError,:);
        vrSource_spk = mean(mrObs .* mrV_fit(viError,:)')';
        return;
    case 9  % center pairwise correlation
        nx = 32;
        [xx,yy,zz] = meshgrid(xyz0(1) + linspace(-32,32,nx), xyz0(2) + linspace(-32,32,nx), xyz0(3) + logspace(-1,2,nx));
        mrXYZ = [xx(:), yy(:), zz(:)];
        mrD = pdist2(mrXYZ, mrPos_site); mrObs_n = fh_normalize(mrObs);
        mrV_fit = 1./(mrD .* mrD(:,1));
        [vrError, viError] = min(1 - fh_normalize(mrV_fit')'*mrObs_n);
        mrPos_spk = mrXYZ(viError,:);
        vrSource_spk = mean(mrObs .* mrD(viError,:)')';
        return;
    case 8
        [xx,yy,zz] = meshgrid(xyz0(1) + linspace(-32,32,nx), xyz0(2) + linspace(-32,32,nx), xyz0(3) + logspace(-1,2,nx));
        mrXYZ = [xx(:), yy(:), zz(:)];
        mrD = pdist2(mrXYZ, mrPos_site); mrObs_n = fh_normalize(mrObs);
        [~, viError] = min(1 - fh_normalize(1./mrD')'*mrObs_n);
        mrPos_spk = mrXYZ(viError,:);
        vrSource_spk = mean(mrObs .* mrD(viError,:)')';
        return;
    case 7
        % first try forward model and find optimal
%         mrPos_spk = zeros(4, nSpk);
        nx = 20;
        [xx,yy,zz] = meshgrid(xyz0(1) + linspace(-20,20,nx), xyz0(2) + linspace(-20,20,nx), xyz0(3) + linspace(0.1,40,nx));
        mrXYZ = [xx(:), yy(:), zz(:)];
%         mrD = pdist2(mrXYZ, mrPos_site, 'squaredeuclidean'); mrObs_n = fh_normalize(mrObs.^2);
        mrD = pdist2(mrXYZ, mrPos_site); mrObs_n = fh_normalize(mrObs);
        [vrError, viError] = min(1 - fh_normalize(1./mrD')'*mrObs_n);
        vrAmp = mean(mrObs_n .* mrD(viError,:)');
        % initial guess
        mrXYZA0 = [mrXYZ(viError,:), vrAmp(:)]';
        
        parfor iSpk = 1:nSpk
            xyza0 = mrXYZA0(:,iSpk)';
            obs1 = mrObs_n(:,iSpk)';
            fh_forward1 = @(xyza)xyza(4)./sqrt(sum((mrPos_site-xyza(1:3))'.^2));
            fh_error1 = @(xyza)sum(abs(fh_forward1(xyza) - obs1));
            vrErr0_spk(iSpk) = fh_error1(xyza0);
            [xyza_, vrErr_spk(iSpk)] = fminsearch(@(xyza)fh_error1(xyza), xyza0, S_opt);
            mrPos_spk(:,iSpk) = xyza_(1:3);
            vrSource_spk(iSpk) = xyza_(4);
        end        
        % return
        mrPos_spk = mrPos_spk';
        vrSource_spk = vrSource_spk .* fh_norm(mrObs);
        return;
    case 6
        fh_forward1 = @(xyz)fh_normalize(sum((mrPos_site-xyz)'.^2));
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fminsearch(@(xyz)1-fh_forward1(xyz) * mrObs_sq_inv_norm(:,iSpk), xyz0, S_opt);
        end
    case 5
        fh_forward1 = @(xyz)(1./sqrt(sum((mrPos_site-xyz)'.^2)));
        fh_penalty1 = @(xyz)sum(abs(xyz-xyz0)) * 1e-12;
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fminsearch(@(xyz)fh_penalty1(xyz) - fh_forward1(xyz) * mrObs(:,iSpk), xyz0, S_opt);
        end
    case 4
        fh_forward1 = @(xyz)fh_normalize(1./sqrt(sum((mrPos_site-xyz)'.^2)));
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fminsearch(@(xyz)1-fh_forward1(xyz) * mrObs(:,iSpk), xyz0, S_opt);
        end
    case 3
        fh_search = @(obs_)fmincon(@(xyz)fh_err(xyz,obs_'), xyz0, [1 1 1], 1000);
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fh_search(mrObs(:,iSpk));
        end
    case 2
        fh_search = @(obs_)fminunc(@(xyz)fh_err(xyz,obs_'), xyz0, S_opt);
        S_opt = struct('Display', 'off', 'MaxIter', MaxIter*2, 'TolFun', .01);
        tic
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fh_search(mrObs(:,iSpk));
        end
        toc
    case 1
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fminsearch(@(xyz)fh_err(xyz,mrObs(:,iSpk)'), xyz0, S_opt);    
        end
%     case 2
%         fh_search = @(obs_)fminsearch(@(xyz)fh_err(xyz,obs_'), gpuArray(double(xyz0)), S_opt);              
%         mrObs = gpuArray(double(mrObs));
%         [mrPos_spk, vrErr_spk] = arrayfun(@(i)fh_search(mrObs(:,i)), 1:size(mrObs,2));
%         toc
%         [mrPos_spk, vrErr_spk] = arrayfun(@fh_search, mrObs);
end %switch
mrPos_spk=mrPos_spk';
vrSource_spk = sqrt(abs(sum(mrObs.^2) ./ sum(1./pdist2(mrPos_site, mrPos_spk, 'squaredeuclidean'))))';

% redo poor fit and fit two sources and take stronger 
if ~isempty(thresh_error)
    viSelect1 = find(vrErr_spk>thresh_error);
    [mrPos_spk1, vrErr_spk1] = deal(mrPos_spk(viSelect1,:), vrErr_spk(viSelect1));
    [mrPos_spk2, vrErr_spk2] = fit_two_monopoles_(mrObs(:,viSelect1), mrPos_site, xyz0, S_opt);
    viiReplace = find(vrErr_spk2 < vrErr_spk1);
    viReplace1 = viSelect1(viiReplace);
    mrPos_spk(:,viReplace1) = mrPos_spk2(:,viiReplace);
    vrErr_spk(viReplace1) = vrErr_spk2(viiReplace);
end
mrPos_spk(:,3) = abs(mrPos_spk(:,3));
% vrErr_spk = 1./vrErr_spk;

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


%--------------------------------------------------------------------------
function [pos, mag] = com_mr_(mrPos, vrV)
vrV2 = vrV(:).^2;
v2_sum = mean(vrV2);
mag = sqrt(v2_sum);
pos = mean(mrPos .* vrV2) ./ v2_sum;
end %func


%--------------------------------------------------------------------------
% multiply spikes by spikes
function mr12 = mr2mr_mtimes_(mr1, mr2, idimm)
if nargin<2, mr2=[]; end
if isempty(mr2), mr2=mr1; end
if idimm == 2
    mr1=mr1';
    mr2=mr2';
end
dimm1 = size(mr1);
dimm2 = size(mr2);
mr12 = zeros(dimm1(1) * dimm2(1), dimm1(2), 'single');
for iSpk = 1:dimm1(2)
    vr_ = mr1(:,iSpk) * mr2(:,iSpk)';
    mr12(:,iSpk) = vr_(:);
end %for
if idimm==2
    mr12 = mr12';
end
end %func


%--------------------------------------------------------------------------
function [mrPos_spk1, vrErr_spk1] = fit_two_monopoles_(mrObs1, mrPos_site, xyz0, S_opt)
nSpk1 = size(mrObs1,2);
mrPos_spk1 = zeros(3, nSpk1);
vrErr_spk1 = zeros(nSpk1,1);

fh_forward1 = @(xyzxyzp)1./pdist2(xyzxyzp(1:3), mrPos_site) + xyzxyzp(7)./pdist2(xyzxyzp(4:6), mrPos_site);
fh_err1 = @(xyz12p, obs)pdist2(fh_forward1(xyz12p), obs, 'cosine'); % todo: speed up cosine distance calculation
fh_search1 = @(obs_)fminsearch(@(xyzxyzp)fh_err1(xyzxyzp,obs_'), [xyz0, xyz0, 1], S_opt);

parfor iSpk1 = 1:nSpk1
    [y_, vrErr_spk1(iSpk1)] = fh_search1(mrObs1(:,iSpk1));
    if y_(end) > 1
        mrPos_spk1(:,iSpk1) = y_(4:6);
    else
        mrPos_spk1(:,iSpk1) = y_(1:3);
    end
end
mrPos_spk1 = mrPos_spk1';
end %func



%         
%         [mrSrt_, miSrt_] = sort(mr_, 'descend');
%         [mrSrt_, miSrt_] = deal(mrSrt_(1:nAve_localize,:), miSrt_(1:nAve_localize,:));
    
%     end

% 
% 
% 
% mrPos12 = mrB1 * abs(double(mrFet12));
% [vrMax12, viMax12] = max(mrPos12);
% mrPos12 = mrPos_field(viMax12,:);



%--------------------------------------------------------------------------
function mr12 = pdist2_(mr1, mr2)
mr12 = sqrt(bsxfun(@plus, sum(mr2'.^2), bsxfun(@minus, sum(mr1'.^2)', 2*mr1*mr2')));
end %func