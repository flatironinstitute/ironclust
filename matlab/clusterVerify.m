function [mrMiss, mrFp, vnSpk_gt, miCluMatch, S_score_clu] = clusterVerify(viCluGt, viTimeGt, viClu, viTime, jitter)
% usage
% -----
% [mrMiss, mrFp, vnCluGt, miCluMatch, S_score_clu] = clusterVerify(viCluGt, viTimeGt, viClu, viTime, jitter)
% [S] = clusterVerify(viCluGt, viTimeGt, viClu, viTime, jitter)

if nargin<5, jitter = 12; end

[viTimeGt, viTime] = deal(double(viTimeGt), double(viTime));
[viCluGt, viClu] = deal(int32(viCluGt), int32(viClu));
[cviSpk_gt, nCluGt] = vi2cell_(viCluGt);
[cviSpk_clu, nClu] = vi2cell_(viClu);
cviTime_gt = cellfun(@(vi_)unique(int32(viTimeGt(vi_)/jitter)), cviSpk_gt, 'UniformOutput', 0);
cviTime_clu = cellfun(@(vi_)unique(int32(viTime(vi_)/jitter)), cviSpk_clu, 'UniformOutput', 0);
vnSpk_clu = cellfun(@numel, cviTime_clu);
vnSpk_gt = cellfun(@numel, cviTime_gt);
[mrMiss, mrFp, mrAccuracy] = deal(nan(nClu, nCluGt, 'single'));
vnDetected_gt = zeros(nCluGt,1,'int32');

[cviHit_gt, cviMiss_gt, cviHit_clu, cviMiss_clu, ...
    cviSpk_gt_hit, cviSpk_gt_miss, cviSpk_clu_hit, cviSpk_clu_miss] = ...
    deal(cell(1, nCluGt));
viTime0 = (int32(viTime/jitter));
fprintf('Validating cluster\n\t');
t1 = tic;
for iCluGt1=1:nCluGt
    n_gt1 = vnSpk_gt(iCluGt1);
    if n_gt1==0, continue; end
    rGT1 = cviTime_gt{iCluGt1};
    lim_gt1 = [min(rGT1), max(rGT1)];    
    vnDetected_gt(iCluGt1) = count_overlap_(rGT1, viTime0);
    [vrMiss_, vrFp_, vrAccuracy_] = deal(zeros(nClu, 1, 'single'));
    for iClu=1:nClu        
        rComp1 = cviTime_clu{iClu};
        n_clu1 = vnSpk_clu(iClu);
        if a_in_b_(lim_gt1, rComp1, 0)
            n_gt1_clu1 = count_overlap_(rGT1, rComp1);
        else
            n_gt1_clu1 = 0;
        end
        vrMiss_(iClu) = (n_gt1-n_gt1_clu1)/n_gt1;
        vrFp_(iClu) = (n_clu1-n_gt1_clu1)/n_clu1;  %it was /n2 before sep 1 2018 
        vrAccuracy_(iClu) = n_gt1_clu1 / (n_clu1+n_gt1-n_gt1_clu1);
    end
    [mrMiss(:,iCluGt1), mrFp(:,iCluGt1), mrAccuracy(:,iCluGt1)] = ...
        deal(vrMiss_, vrFp_, vrAccuracy_);
%     fprintf('.');
end
fprintf('\n\ttook %0.1fs.\n', toc(t1));
[mrMiss, mrFp, vnDetected_gt, vnSpk_gt] = multifun_(@gather, mrMiss, mrFp, vnDetected_gt, vnSpk_gt);
vrDetected = double(vnDetected_gt) ./ double(vnSpk_gt);
% mrScore = 1-mrFp-mrMiss;
mrScore = mrAccuracy;
[vrScore, viCluMatch] = max(mrScore, [], 1);
[~, miCluMatch] = sort(mrScore, 'descend');
vi_match = sub2ind(size(mrMiss), viCluMatch, 1:numel(viCluMatch));
[vrMiss, vrFp] = deal(mrMiss(vi_match), mrFp(vi_match));
% vrScore = 1-vrMiss-vrFp; vrScore(vrScore<0)=0;

for iCluGt=1:nCluGt
    viSpk_gt1 = cviSpk_gt{iCluGt};
    viTime_gt1 = viTimeGt(viSpk_gt1);    
    iClu1 = viCluMatch(iCluGt);
    viSpk_clu1 = cviSpk_clu{iClu1};
    viTime_clu1 = viTime(viSpk_clu1);
    [vlSpk_gt1, vlSpk_clu1, viiSpk_gt1, viiSpk_clu1] = time_match2_(viTime_gt1, viTime_clu1, jitter);    
    cviHit_gt{iCluGt} = viTime_gt1(viiSpk_gt1);
    cviMiss_gt{iCluGt} = viTime_gt1(~vlSpk_gt1);
    cviHit_clu{iCluGt} = viTime_clu1(viiSpk_clu1);
    cviMiss_clu{iCluGt} = viTime_clu1(~vlSpk_clu1);
    cviSpk_gt_hit{iCluGt} = viSpk_gt1(vlSpk_gt1);
    cviSpk_gt_miss{iCluGt} = viSpk_gt1(~vlSpk_gt1);
    cviSpk_clu_hit{iCluGt} = viSpk_clu1(vlSpk_clu1);
    cviSpk_clu_miss{iCluGt} = viSpk_clu1(~vlSpk_clu1);  
    cvlHit_gt{iCluGt} = vlSpk_gt1;
    cvlHit_clu{iCluGt} = vlSpk_clu1;
end
[vrAccuracy, viCluMatch_accuracy] = max(mrAccuracy);
S_score_clu = makeStruct_(vrScore, vrMiss, vrFp, viCluMatch, cviHit_gt, ...
    cviMiss_gt, cviHit_clu, cviMiss_clu, cviSpk_gt_hit, cviSpk_gt_miss, ...
    cviSpk_clu_hit, cviSpk_clu_miss, vrAccuracy, viCluMatch_accuracy, ...
    cvlHit_gt, cvlHit_clu);
% viCluMatch1 = viClu_unique(viCluMatch);
% vrScore = 1-vrMiss-vrFp;
func1=@(x)quantile(x, [.25,.5,.75])*100;
fprintf('Validation summary\n');
fprintf('\tEvents detected (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrDetected)*100, func1(vrDetected), sprintf('%0.1f ', vrDetected*100));
fprintf('\tFalse-positives (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrFp)*100, func1(vrFp), sprintf('%0.1f ', vrFp*100));
fprintf('\tFalse-negatives (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrMiss)*100, func1(vrMiss), sprintf('%0.1f ', vrMiss*100));
fprintf('\tAccuracy (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrAccuracy)*100, func1(vrAccuracy), sprintf('%0.1f ', vrAccuracy*100));
fprintf('\tScore (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrScore)*100, func1(vrScore), sprintf('%0.1f ', vrScore*100));
fprintf('\tCluster-size: %s\n', sprintf('%d, ', vnSpk_gt));
fprintf('\tMatching clu: %s\n', sprintf('%d, ', viCluMatch));

if nargout==1
    mrMiss = makeStruct_(mrMiss, mrFp, vnSpk_gt, miCluMatch, S_score_clu);
end
end %func


%--------------------------------------------------------------------------
function flag = a_in_b_(vrA, vrB, fSorted)
if nargin<3, fSorted = []; end
if isempty(vrA) || isempty(vrB), flag = false; return; end
% assume sorrted
if fSorted
    limA = [vrA(1), vrA(end)];
    limB = [vrB(1), vrB(end)];
else
    limA = [min(vrA), max(vrA)];
    limB = [min(vrB), max(vrB)];    
end
flag = (limA(1) >= limB(1) && limA(1) <= limB(2)) || ...
       (limA(2) >= limB(1) && limA(2) <= limB(2)) || ...
       (limB(1) >= limA(1) && limB(1) <= limA(2)) || ...
       (limB(2) >= limA(1) && limB(2) <= limA(2));
end %func


%--------------------------------------------------------------------------
function [cviSpk_site, nSites, vi_site] = vi2cell_(viSite_spk, nSites)
if nargin<2, nSites = []; end
if isempty(nSites), nSites = max(viSite_spk); end

% based on unique() function, which sorts. faster than arrayfun.
cviSpk_site = cell(nSites, 1);
[vr, vi] = sort(viSite_spk);
vi_change = [1; find(diff(vr(:))>0)+1; numel(viSite_spk)+1];
if isempty(viSite_spk), vi_site=[]; return; end

vi_site = vr(vi_change(1:end-1));
vl_remove = vi_site < 1;
if any(vl_remove)
    vi_site(vl_remove) = [];
    vi_change(vl_remove) = [];
end
for iStep = 1:numel(vi_site)
    cviSpk_site{vi_site(iStep)} = vi(vi_change(iStep):vi_change(iStep+1)-1);
end
vi_site = vi_site(:)';
end %func


%--------------------------------------------------------------------------
function nOverlap = count_overlap_(viGt, viTest)
nOverlap = sum(ismember(viGt, viTest) | ismember(viGt, viTest+1) | ismember(viGt, viTest-1));
% viGt_diff = setdiff(setdiff(setdiff(viGt, viTest), viTest+1), viTest-1);
% nOverlap = numel(viGt) - numel(viGt_diff);
end %func


%--------------------------------------------------------------------------
% function nOverlap = count_overlap_(viGt, viTest)
% viGt_diff = setdiff(setdiff(setdiff(viGt, viTest), viTest+1), viTest-1);
% nOverlap = numel(viGt) - numel(viGt_diff);
% end %func


%--------------------------------------------------------------------------
function [vlA, vlB, viA1, viB1] = time_match2_(viA, viB, jitter)
% A: ground truth, B: matching unit
if nargin<3, jitter=25; end

if jitter>0
    viA = (double(viA)/jitter);
    viB = (double(viB)/jitter);
end
viA = int32(viA);
viB = int32(viB);
vlA = false(size(viA));
vlB = false(size(viB));
for i1=-1:1
    for i2=-1:1
        vlA = vlA | ismember(viA+i1, viB+i2);
    end
end
if nargout==1, return; end

%viA_match = find(vlA);
viA1 = find(vlA);
viB1 = zeros(size(viA1));
viA11 = viA(viA1);
for iA=1:numel(viA1)
    [~, viB1(iA)] = min(abs(viB - viA11(iA)));
end
vlB(viB1) = 1;
end %func


%--------------------------------------------------------------------------
function [ S ] = makeStruct_( varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name

S=[];
for i=1:nargin
    S = setfield(S, inputname(i), varargin{i});
end
end %func


%--------------------------------------------------------------------------
function varargout = multifun_(hFun, varargin)
% apply same function to the input, unary function only

if nargout ~= numel(varargin), error('n arg mismatch'); end
for i=1:nargout
    try
        varargout{i} = hFun(varargin{i});
    catch
        varargout{i} = varargin{i};
    end
end
end %func