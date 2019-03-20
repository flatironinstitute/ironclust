% goal: minimize cosine distance 
% apply fminbnd from the max position
% return strength
function [mrPos_spk, vrSource_spk, vrErr_spk] = search_monopole(mrPos_site, trFet, xyz0, iFet, MaxIter)
% try excluding the nearest electrode
if nargin<4, iFet = 1; end
if nargin<5, MaxIter = 30; end



% initial position
fh_forward = @(xyz)1./pdist2(xyz, mrPos_site);
switch 1
    case 2, fh_err = @(xyz, obs)pdist2(fh_forward(xyz), obs, 'correlation');
    case 1, fh_err = @(xyz, obs)pdist2(fh_forward(xyz), obs, 'cosine');
end
mrObs = double(abs(trFet(:,:,iFet)));

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
% mrPos_spk=mrPos_spk';
vrErr_spk = zeros(nSpk,1);
parfor iSpk = 1:nSpk
% iSpk = 1;
    obs_ = mrObs(:,iSpk);
%     switch 1
%         case 1, [~,imax_] = max(obs_); xyz0_ = mrPos_site(imax_,:); xyz0_(3) = .1;
% %         case 2, xyz0_ = mrPos_spk(:,iSpk)';
%     end
%     
    [mrPos_spk(:,iSpk), vrErr_spk(iSpk)] = fminsearch(@(xyz)fh_err(xyz,obs_'), xyz0, struct('Display', 'off', 'MaxIter', MaxIter, 'TolFun', .01));
%     mrPos(:,iSpk) = lsqnonneg(, mrObs(:,iSpk))
%     fprintf('.');
end
% toc
mrPos_spk=mrPos_spk';

vrSource_spk = sqrt(abs(sum(mrObs.^2) ./ sum(1./pdist2(mrPos_site, mrPos_spk, 'squaredeuclidean'))))';

mrPos_spk(:,3) = abs(mrPos_spk(:,3));
% vrErr_spk = 1./vrErr_spk;

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
