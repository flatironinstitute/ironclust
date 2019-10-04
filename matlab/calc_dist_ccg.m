% optimized from https://github.com/MouseLand/Kilosort2/blob/c3573704a18bc3ca42f866a45145418b40a9cdb9/postProcess/ccg.m
% Original author: Marius Pachitariu, Feb 18 2019
% Modified by James Jun on Oct 4, 2019
function mrDist_clu = calc_dist_ccg(S0, P)
if nargin<2, P=[]; end
if isempty(P), P = S0.P; end

fprintf('Running calc_dist_ccg:\n\t'); t1=tic;
[Q_THRESH, R_THRESH, NBINS, DT] = deal(.2, .05, 500, .001);
[sRateHz, S_clu, miSites] = deal(P.sRateHz, S0.S_clu, P.miSites);

nClu = max(S_clu.viClu);
cviSpk_clu = arrayfun(@(x)find(S_clu.viClu==x), 1:nClu, 'UniformOutput', 0);
vnSpk_clu = cellfun(@(x)numel(x), cviSpk_clu);
cvrTime_clu = cellfun(@(x)double(S0.viTime_spk(x))/sRateHz, cviSpk_clu, 'UniformOutput', 0);
viSite_clu = cellfun(@(x)mode(S0.viSite_spk(x)), cviSpk_clu);
miSites_clu = miSites(1:P.nSites_fet,viSite_clu);
[mrQ_clu, mrR_clu] = deal(nan(nClu, 'single'));
for iClu1 = 1:nClu
    vrTime1 = cvrTime_clu{iClu1};
    
    % preselection of clusters
    viClu2 = find(vnSpk_clu > vnSpk_clu(iClu1) & ismember(viSite_clu, miSites_clu(:,iClu1)));
    if isempty(viClu2), continue; end
    
    [vrQ1, vrR1] = deal(zeros(size(viClu2), 'single'));
    for iiClu2 = 1:numel(viClu2)
        vrTime2 = cvrTime_clu{viClu2(iiClu2)};
        [~, Qi, Q00, Q01, rir] = ccg_(vrTime1, vrTime2, NBINS, DT);
        vrQ1(iiClu2) = min(Qi/(max(Q00, Q01)));
        vrR1(iiClu2) = min(rir);    
    end
    mrQ_clu(viClu2,iClu1) = vrQ1;
    mrR_clu(viClu2,iClu1) = vrR1;
    fprintf('.');
end %for
mrDist_clu = mrQ_clu < Q_THRESH & mrR_clu < R_THRESH;
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function [K, Qi, Q00, Q01, Ri] = ccg_(st1, st2, nbins, tbin)
% cross correlogram code
% 
% Inputs
% -----
% st1, st2: spike times in seconds unit
% nbins: number of bins in cross correlogram
% tbin: time bin resolution in sec
%
% Outputs
% -----
% K, Q00, Q01: 2*nbins+1 x 1
% Qi, Ri: 10 x 1

st1 = sort(st1(:));
st2 = sort(st2(:));

dt = nbins*tbin;


T = (max([st1; st2])-min([st1; st2]));

ilow = 1;
ihigh = 1;
j = 1;

K = zeros(2*nbins+1, 1); 

while j<=numel(st2)
%     disp(j)
    while (ihigh<=numel(st1)) && (st1(ihigh) < st2(j)+dt)
        ihigh = ihigh + 1;
    end
    while (ilow<=numel(st1)) && st1(ilow) <= st2(j)-dt
        ilow = ilow + 1;
    end
    if ilow>numel(st1)
        break;
    end
    if st1(ilow) > st2(j)+dt
        j = j+1;
        continue;
    end
    for k = ilow:(ihigh-1)
        ibin = round((st2(j) - st1(k))/tbin);
%         disp(ibin)
        K(ibin+ nbins+1) = K(ibin + nbins+1) + 1;
    end
    j = j+1;
end

irange1 = [2:nbins/2 (3/2*nbins+1):(2*nbins)];
Q00 = sum(K(irange1)) / (numel(irange1)*tbin* numel(st1) * numel(st2)/T);

irange2 = [nbins+1-50:nbins-10];
irange3 = [nbins+12:nbins+50];

R00 = max(mean(K(irange2)), mean(K(irange3)));
R00 = max(R00, mean(K(irange1)));

Q01 = sum(K(irange2)) / (numel(irange2)*tbin* numel(st1) * numel(st2)/T);
Q01 = max(Q01, sum(K(irange3)) / (numel(irange3)*tbin* numel(st1) * numel(st2)/T));

% disp(R00)

a = K(nbins+1);
K(nbins+1) = 0;
for i = 1:10
    irange = [nbins+1-i:nbins+1+i];
    Qi0 = sum(K(irange)) / (2*i*tbin* numel(st1) * numel(st2)/T);
    Qi(i) = Qi0;
    
    n = sum(K(irange))/2;
    lam = R00 * i;
    
%     logp = log(lam) * n - lam - gammaln(n+1);
    p = 1/2 * (1+ erf((n - lam)/sqrt(2*lam)));
    
    Ri(i) = p;
    
end

K(nbins+1) = a;
end %func