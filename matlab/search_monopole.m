% goal: minimize cosine distance 
% apply fminbnd from the max position
% return strength
function [mrPos_spk, vrSource_spk, vrErr_spk] = search_monopole(mrPos_site, mrObs, xyz0, MaxIter)
% try excluding the nearest electrode
% if nargin<4, iFet = 1; end
% if nargin<5, MaxIter = 30; end



% initial position
fh_cos = @(x,y)1-sum(x.*y); %assume normalized
fh_norm = @(x)x./sqrt(sum(x.^2));
fh_forward = @(xyz)fh_norm(1./pdist2(xyz, mrPos_site));
switch 3
    case 3, fh_err = @(xyz, obs)fh_cos(fh_forward(xyz), obs);
    case 2, fh_err = @(xyz, obs)pdist2(fh_forward(xyz), obs, 'correlation');
    case 1, fh_err = @(xyz, obs)pdist2(fh_forward(xyz), obs, 'cosine');
end
mrObs = double(abs(mrObs));
mrObs = mrObs ./ sqrt(sum(mrObs.^2)); % normalize output

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
vrErr_spk = zeros(nSpk,1);
S_opt = struct('Display', 'off', 'MaxIter', MaxIter, 'TolFun', .01);
switch 1
    case 2
        fh_search = @(obs_)fminunc(@(xyz)fh_err(xyz,obs_'), xyz0, S_opt);
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fh_search(mrObs(:,iSpk));
        end
    case 1        
        fh_search = @(obs_)fminsearch(@(xyz)fh_err(xyz,obs_'), xyz0, S_opt);
        parfor iSpk = 1:nSpk
            [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fh_search(mrObs(:,iSpk));
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
if 0
    thresh_error = .08;
    viSelect = find(vrErr_spk>thresh_error);
    [mrPos_spk(viSelect,:), vrErr_spk(viSelect)] = fit_two_monopoles_(mrObs(:,viSelect), mrPos_site, xyz0, S_opt);
end
mrPos_spk(:,3) = abs(mrPos_spk(:,3));
% vrErr_spk = 1./vrErr_spk;

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
