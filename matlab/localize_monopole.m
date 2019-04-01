% goal: minimize cosine distance 
% apply fminbnd from the max position

function [mrPos_spk, vrSource_spk] = localize_monopole(mrPos_site, mrPos_field, trFet12, nAve)
if nargin<4, nAve = []; end
if isempty(nAve), nAve=8; end


% inverse matrix
% mrDist = pdist2(mrPos_site, mrPos_field);
% mrA = 1./mrDist;
mrB = pinv(mrA); % monopole assumption
[vx_field, vy_field, vz_field] = deal(mrPos_field(:,1), mrPos_field(:,2), mrPos_field(:,3)); 

nFet = size(trFet12,3);
trFet12 = double(trFet12);
mrSource = [];
for iFet = 1:nFet
    mrObs_ = abs(trFet12(:,:,iFet));
    mr_ = mrB * mrObs_;
    if isempty(mrSource)
        mrSource = mr_;
    else
        mrSource = mr_ + mrSource;
    end
    
    % compute error and pick least error points   
end      
[mrSrt1, miSrt1] = sort(single(mrSource), 'descend');
[mrSrt1, miSrt1] = deal(mrSrt1(1:nAve,:), miSrt1(1:nAve,:));
nSpk = size(trFet12,2);

fh_com = @(x,y)sum(x.*y.^2) ./ sum(y.^2);
mrPos_spk = [fh_com(vx_field(miSrt1), mrSrt1); fh_com(vy_field(miSrt1), mrSrt1); fh_com(vz_field(miSrt1), mrSrt1)]';
vrSource_spk = sqrt(sum(mrSrt1));



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