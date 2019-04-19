function [mrB, mrW, miSites] = spatial_whiten(mrA, miSites)
% spatial_whiten(mrA, miSites)
% mrB = spatial_whiten(mrA)

if nargin<2, miSites = []; end

persistent mrW_ miSites_

if isempty(mrW_)
    miSites_ = miSites;
    mrW_ = compute_whiten_(mrA, miSites);
    fprintf(2, 'spatial_whiten: computed whitening matrix\n');
end
if nargout==0, return; end % just compute cache and exit

[nSamples, nSites] = size(mrA);
fGpu = isGpu_(mrA);
try
    mrB = zeros(size(mrA), 'like', mrA);
catch % no gpu memory
    mrA = gather_(mrA);
    mrB = zeros(size(mrA), 'like', mrA);
    fGpu = 0;
end
[mrW, fGpu] = gpuArray_(mrW_, fGpu);
miSites = miSites_;

for iSite = 1:nSites
    viSites1 = miSites(:,iSite);
    mrB(:,iSite) = single(mrA(:,viSites1)) * mrW(:,iSite); % single to int16 may occur
end

end %func


%--------------------------------------------------------------------------
% assume data is mean-subtracted already
function [mrW, miSites] = compute_whiten_(mr, miSites)
eps = 1e-6;
scaleproc = 100; % int16 scaling of whitened data

[nSamples, nSites] = size(mr);
nSites1 = size(miSites,1);
mrW = zeros(nSites1, 'single');
for iSite = 1:nSites
    viSites1 = miSites(:,iSite);
    mr1 = single(mr(:,viSites1));    
    [mrE1, vrD1] = svd((mr1' * mr1) / nSamples);
    vrD1 = diag(vrD1);
    switch 2
        case 2, scale1 = scaleproc;
        case 1, scale1  = mean((vrD1+eps) ./ (vrD1+1e-6)).^.5 * scaleproc;
    end
    mrW1 = mrE1 * diag(scale1./(vrD1 + eps).^.5) * mrE1';    
    mrW(:,iSite)  = gather_(mrW1(:,1));
end %for
mrW = mrW; 
end %func


%--------------------------------------------------------------------------
% copyed from irc.m
% isGpu, gpuArray_
%--------------------------------------------------------------------------
function mr = gather_(mr)
try
    mr=gather(mr);
catch
    ;
end
end %func



%--------------------------------------------------------------------------
function flag = isGpu_(vr)
try
    flag = isa(vr, 'gpuArray');
catch
    flag = 0;
end
end


%--------------------------------------------------------------------------
function [mr, fGpu] = gpuArray_(mr, fGpu)
if nargin<2, fGpu = 1; end
% fGpu = 0; %DEBUG disable GPU array
if ~fGpu, return; end
try
    if ~isa(mr, 'gpuArray'), mr = gpuArray(mr); end
    fGpu = 1;
catch        
    try % retry after resetting the GPU memory
        gpuDevice(1); 
        mr = gpuArray(mr);
        fGpu = 1;
    catch % no GPU device found            
        fGpu = 0;
    end
end
end