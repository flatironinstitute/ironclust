% works on mda format only for now
% no manual interface now
% analyzes single mda file format
% todo: multiple mda files
% todo: multiple formats

%--------------------------------------------------------------------------
function irc2(vcDir_in, vcDir_out, vcFile_arg)
if nargin<1, vcDir_in = ''; end
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end

% batch processing. it uses default param for now
if iscell(vcDir_in) && iscell(vcDir_out)
    [csDir_in, csDir_out] = deal(vcDir_in, vcDir_out);
    for iFile = 1:numel(csDir_in)
        set(0,'UserData',[]); % clear
        try
            fprintf('irc2 batch-processing %s (%d/%d)\n', csDir_in{iFile}, iFile, numel(csDir_in));
            irc2(csDir_in{iFile}, csDir_out{iFile});
        catch
            disp(lasterr());
        end
    end 
    return;
end


[vcCmd, vcArg1, vcArg2] = deal(vcDir_in, vcDir_out, vcFile_arg); 

if isempty(vcDir_in),  vcDir_in = get_test_data_('static'); end
if isempty(vcDir_out), vcDir_out = strrep(vcDir_in, 'groundtruth', 'irc'); end
if isempty(vcFile_arg)
    vcFile_arg = struct(...
        'nPcPerChan',1,'fGpu',1,'fParfor',1,'qqFactor',4,'nPc_spk',6,'MAX_BYTES_LOAD',.5e9,'spkLim_ms',[-.25,.75]*2,'maxWavCor',.9, 'step_sec_drift', 20, 'batch_sec_drift', 300); 
end
if ~exist_dir_(vcDir_out), mkdir(vcDir_out); end

% special commands
switch lower(vcCmd)
    case 'plot', irc('plot', vcArg1, vcArg2); return;
    case 'verify', irc('plot', 'verify'); return;
    case 'clear', clear_(); return;
    case 'clear-sort', clear_('sort'); return;
    case 'test-static', clear_(); irc2(get_test_data_('static'), [], []); return;
    case 'test-drift', clear_(); irc2(get_test_data_('drift'), [], []); return;
    case 'sort', irc2(clear_('sort')); return;
    case 'export', irc('export', vcArg1); return;   
end

P = makeParam_(vcDir_in, vcDir_out, vcFile_arg);
S0 = get(0, 'UserData');
if isempty(struct_get_(S0, 'trPc_spk'))
    S0 = detect_(P); 
    set(0, 'UserData', S0);
end
if isempty(struct_get_(S0, 'S_clu'))
    S0 = sort_(S0, P);
    set(0, 'UserData', S0);
end
S0 = auto_(S0, P);
set(0, 'UserData', S0);
describe_(S0);

% output
struct_save_(S0, fullfile(vcDir_out, 'out_irc.mat'), 1);
vcFile_firings_mda = fullfile(vcDir_out, 'firings.mda');
save_firings_mda_(S0, vcFile_firings_mda);

% Validate
vcFile_gt_mda = fullfile(vcDir_in, 'firings_true.mda');
if exist_file_(vcFile_gt_mda)
    raw_fname = fullfile(vcDir_in, 'raw.mda');
    S_score = irc('validate-mda', vcFile_gt_mda, vcFile_firings_mda, raw_fname); % assume that groundtruth file exists
    struct_save_(S_score, fullfile(vcDir_out, 'raw_geom_score.mat'), 1);
end
end %func


%--------------------------------------------------------------------------
function vcDir_in = clear_(vcMode)
if nargin<1, vcMode = 'all'; end

switch vcMode
    case 'all'
        set(0, 'UserData', []);
        disp('irc2: cleared all');
        vcDir_in = '';
    case 'sort'
        S0 = set(0, 'UserData');
        S0.S_clu = [];
        set(0,'UserData',S0);
        disp('irc2: cleared sort');
        try
            vcDir_in = fileparts(S0.P.vcFile);
        catch
            vcDir_in = '';
        end
end
end

%--------------------------------------------------------------------------
function csDesc = describe_(S0, vcFile_firings_mda)
runtime_total = S0.runtime_detect + S0.runtime_sort + S0.runtime_automerge;
tDur = recording_duration_(S0.P, S0); 
memory_sort = S0.memory_sort - S0.memory_init;
memory_detect = S0.memory_detect - S0.memory_init;

csDesc = {};

csDesc{end+1} = sprintf('Runtime (s):');
csDesc{end+1} = sprintf('    Detect + feature (s):   %0.1fs', S0.runtime_detect);    
csDesc{end+1} = sprintf('    Cluster (s):    %0.1fs', S0.runtime_sort);
csDesc{end+1} = sprintf('    Automerge (s):    %0.1fs', S0.runtime_automerge);
csDesc{end+1} = sprintf('    Total runtime (s):      %0.1fs', runtime_total);
csDesc{end+1} = sprintf('    Runtime speed           x%0.1f realtime', tDur / runtime_total);    

csDesc{end+1} = sprintf('memory usage (GiB):');
csDesc{end+1} = sprintf('    detect(GiB):     %0.3f', memory_detect/2^30);
csDesc{end+1} = sprintf('    sort(GiB):     %0.3f', memory_sort/2^30);

if nargout==0
    cellfun(@(x)disp(x), csDesc);
end
end %func


%--------------------------------------------------------------------------
function save_firings_mda_(S0, vcFile_firings_mda)
S_gt1 = struct('viTime', S0.viTime_spk, 'viSite', S0.viSite_spk, 'viClu', S0.S_clu.viClu);
gt2mda_(S_gt1, vcFile_firings_mda);
end %func
            

%--------------------------------------------------------------------------
% auto merge
function S0 = auto_(S0, P)
if nargin<2, P=[]; end
if nargin<1, S0=[]; end
if isempty(S0)
    S0 = get(0, 'UserData');
elseif isempty(P)
    P = S0.P;
end

fprintf('\nauto-merging...\n'); runtime_automerge=tic;
% auto merge
P.post_merge_mode0 = [12];
S0.S_clu = postCluster_(S0.S_clu, P); % peak merging

% fMode = 6;  P.maxWavCor = .97; %82.1, RMS
fMode = 3; %P.maxWavCor = .90;
% fMode = 3;  P.maxWavCor = .97; %82.1, RMS
% fMode = 0;  P.maxWavCor = .90; %wavcor
% fMode = 6;  P.maxWavCor = .97; % 81.6, correlation
% fMode = 7;  P.maxWavCor = .95;  % 81.6, both
switch fMode 
    case 0
        mrDist_clu = wave_similarity_clu_(S0, P, 1);
    otherwise
        mrDist_clu = calc_dist_clu_(S0, P, fMode);

end
S0.S_clu = templateMatch_post_(S0.S_clu, P, mrDist_clu);

% P.post_merge_mode0 = [12, 19];
% S0.S_clu = postCluster_(S0.S_clu, P); % peak merging

S0.S_clu = S_clu_refresh_(S0.S_clu);

S0.runtime_automerge = toc(runtime_automerge);
fprintf('\n\tauto-merging took %0.1fs (fGpu=%d, fParfor=%d)\n', ...
    S0.runtime_automerge, P.fGpu, P.fParfor);
end %func


%--------------------------------------------------------------------------
function mrDist_clu = calc_dist_clu_(S0, P, fMode)
% output
% -----
% mrDist_clu: 0..1, 1 being most similar, 0 being most dissimilar
S_clu = S0.S_clu;
MIN_COUNT = P.min_count;

% compute average waveforms by clusters
nClu = max(S_clu.viClu);
nSites = max(S0.viSite_spk);
if 0
    trPc_spk = S0.trPc_spk(:,1:P.nSites_fet,:);
else
    trPc_spk = S0.trPc_spk;
end
% compute average clu waveforms (sort by centers)
[viClu, vrPos_spk, mrPos_spk] = deal(S_clu.viClu, S0.vrPow_spk, S0.mrPos_spk);
[nPc_spk, nSites_spk, nSpk] = size(S0.trPc_spk);
[cviSpk_drift, mlDrift] = struct_get_(S0.S_clu.S_drift, 'cviSpk_drift', 'mlDrift');
unique_vec_ = @(x)unique(x(:));
if P.fParfor
    gcp_ = gcp();
else
    gcp_ = [];
end

fprintf('calc_dist_clu_ (fParfor=%d): \n\t', P.fParfor); t1=tic;
if ~isempty(cviSpk_drift)
    mrDist_clu = zeros(nClu, 'single');
    cmrDist_clu = cell(size(cviSpk_drift));
    for iDrift = 1:numel(cviSpk_drift)
        switch 1
            case 4
                viSpk1 = cell2mat_(cviSpk_drift(find(mlDrift(:, iDrift))));
                viSpk1 = unique_vec_(S0.S_clu.miKnn(:,viSpk1));
            case 3, viSpk1 = unique_vec_(S0.S_clu.miKnn(:,cviSpk_drift{iDrift}));
            case 2, viSpk1 = cviSpk_drift{iDrift};
            case 1, viSpk1 = cell2mat_(cviSpk_drift(find(mlDrift(:, iDrift))));
        end
        [viClu1, viSite1, trPc1] = deal(viClu(viSpk1), S0.viSite_spk(viSpk1), trPc_spk(:,:,viSpk1)); 
        if isempty(gcp_)
            cmrDist_clu{iDrift} = wav_dist_clu_(viClu1, viSite1, trPc1, fMode, nClu, P);
        else
            vS_out(iDrift) = parfeval(gcp_, ...
                @(x,y,z)wav_dist_clu_(x,y,z, fMode, nClu, P), 1, viClu1, viSite1, trPc1);
        end
        fprintf('.');
    end %for
    if ~isempty(gcp_)
        for iDrift1 = 1:numel(cviSpk_drift)
            [iDrift, mrDist_clu1] = fetchNext(vS_out);
            cmrDist_clu{iDrift} = mrDist_clu1;
        end
    end
    mrDist_clu = max(cat(3,cmrDist_clu{:}), [], 3);
else
    mrDist_clu = pos_dist_clu_(viClu, mrPos_spk, vrPos_spk, nClu, MIN_COUNT);
end
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function mrDist_clu = pos_dist_clu_(viClu, mrPos_spk, vrPow_spk, nClu, MIN_COUNT)

switch 2
    case 3, fh_pow_ = @(x)sqrt(x);
    case 2, fh_pow_ = @(x)log(x);
    case 1, fh_pow_ = @(x)x;
end

% compute average pc per cluster per site
cviSpk_clu = arrayfun(@(x)find(viClu==x), 1:nClu, 'UniformOutput', 0);
[mrPos_med_clu, mrPos_sd_clu] = deal(nan(size(mrPos_spk,2)+1, nClu));
for iClu = 1:nClu
    viSpk1 = cviSpk_clu{iClu};
    if numel(viSpk1) >= MIN_COUNT
        mrPos_spk1 = [mrPos_spk(viSpk1,:), fh_pow_(vrPow_spk(viSpk1))];
        [mrPos_sd_clu(:,iClu), mrPos_med_clu(:,iClu)] = sd_mad_(mrPos_spk1);
    end
end

% compare cluster by cluster
sd_global = nanmedian(mrPos_sd_clu,2);
mrDist_clu = zeros(nClu);
for iClu1 = 1:nClu    
    switch 2
        case 2, vrCorr12 = (mrPos_med_clu - mrPos_med_clu(:,iClu1)) ./ sd_global;  
        case 1
            [mu1, sd1] = deal(mrPos_med_clu(:,iClu1), mrPos_sd_clu(:,iClu1));
            vrCorr12 = (mrPos_med_clu - mu1) ./ sd1;  
    end
    mrDist_clu(:,iClu1) = sqrt(sum(vrCorr12.^2))';
end
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wav_dist_clu_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P)
switch 1
    case 2
        % this uses new method of site referencing
        mrDist_clu = wav_dist_clu2_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P);
    case 1
        mrDist_clu = wav_dist_clu1_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P);
end %switch
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wav_dist_clu2_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P)

% compute average waveforms by clusters
MIN_COUNT = P.min_count;
nSites = max(viSite_spk);
mrDist_clu = zeros(nClu, 'single');

% site loop
cviSpk_site = arrayfun(@(x)find(viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
for iSite = 1:nSites
    viSpk1 = find(viSite_spk==iSite);
    [mrFet12, viSpk12, ~, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, [], P, iSite);
    if isempty(viSpk12), continue; end
    viClu12 = viClu(viSpk12);
    [viClu_uniq1, vnClu_uniq1, cviClu_uniq1] = unique_count_(viClu12);    
    [viClu_uniq1, vnClu_uniq1, cviClu_uniq1] = ...
        multifun_(@(x)x(vnClu_uniq1 >= MIN_COUNT & viClu_uniq1 ~= 0), viClu_uniq1, vnClu_uniq1, cviClu_uniq1);
    if numel(viClu_uniq1)<2, continue; end
    mrPc12 = cell2matfun_(@(x)median2_(mrFet12(:,x),2), cviClu_uniq1');
%     mrPc12 = cell2matfun_(@(x)mean(mrFet12(:,x),2), cviClu_uniq1');
    switch 2
        case 2
            mrDist1 = zeros(size(mrPc12,2), 'single');
            for iClu1 = 1:size(mrPc12,2)
                mrDist1(:,iClu1) = sum((mrPc12 - mrPc12(:,iClu1)).^2) ./ sum(mrPc12(:,iClu1).^2);
            end
            mrDist1 = 1 - min(mrDist1,1);
        case 1, mrDist1 = corr(mrPc12);
    end
    mrDist_clu(viClu_uniq1,viClu_uniq1) = max(mrDist1, mrDist_clu(viClu_uniq1,viClu_uniq1));
end %for
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wav_dist_clu1_(viClu, viSite_spk, trPc_spk, mrPv_global, nClu, P)

% compute average waveforms by clusters
MIN_COUNT = P.min_count;
[ctrPc_mu_clu, ctrPc_sd_clu, cviSite_clu] = deal(cell(nClu,1));
dimm_pc = size(trPc_spk);

% compute average pc per cluster per site
% find centroid of the cluster
for iClu = 1:nClu
    viSpk1 = find(viClu==iClu);
    if isempty(viSpk1), continue; end
    viSite1 = viSite_spk(viSpk1);   
    
    [viSite1_uniq1, vn_uniq1, cviSpk1_uniq1] = unique_count_(viSite1);    
    [viSite1_uniq1, cviSpk1_uniq1] = ...
        multifun_(@(x)x(vn_uniq1>=MIN_COUNT), viSite1_uniq1, cviSpk1_uniq1);
    if isempty(cviSpk1_uniq1), continue; end    
    cviSite_clu{iClu} = viSite1_uniq1;
    [ctrPc_mu_clu{iClu}, ctrPc_sd_clu{iClu}] = med_rms_tr_(trPc_spk(:,:,viSpk1), cviSpk1_uniq1);
end
% construct arrays and sites
[viSite_all, viClu_all] = cvr2vr_vi_(cviSite_clu);
trPc_all = cat(3, ctrPc_mu_clu{:});

% compare cluster by cluster
mrDist_clu = zeros(nClu);
for iClu1 = 1:nClu
    [trPc_mu1, trPc_sd1, viSite1] = deal(ctrPc_mu_clu{iClu1}, ctrPc_sd_clu{iClu1}, cviSite_clu{iClu1});
    vrDist_clu1 = zeros(nClu,1);
    for iiSite1 = 1:numel(viSite1)        
        vi_all2 = find(viSite_all == viSite1(iiSite1) & viClu_all ~= iClu1);
        if isempty(vi_all2), continue; end
        mrPc_mu1 = trPc_mu1(:,:,iiSite1);
        mrPc_sd1 = trPc_sd1(:,:,iiSite1);
        viClu2 = viClu_all(vi_all2);
        trPc2 = trPc_all(:,:,vi_all2);
        vrCorr12 = wavcor_pc_(mrPc_mu1, trPc2, mrPv_global, mrPc_sd1);
        vrDist_clu1(viClu2) = max(vrDist_clu1(viClu2), vrCorr12(:));
    end
    mrDist_clu(:,iClu1) = vrDist_clu1;
end
end %func


%--------------------------------------------------------------------------
function [trPc_med_site, trPc_sd_site] = med_rms_tr_(trPc, cviSpk_site)
if numel(cviSpk_site) == 1
    trPc_med_site = median3_(trPc,3);
    if nargout>=2
        trPc_sd_site = std(trPc,[],3);
    end
else        
    trPc1 = permute(trPc, [3,1,2]);
    trPc_med_site = cell2matfun_(@(x)median3_(trPc1(x,:,:)), cviSpk_site);    
    trPc_med_site = permute(trPc_med_site, [2,3,1]);
    if nargout>=2
        trPc_sd_site = cell2matfun_(@(x)std(trPc1(x,:,:)), cviSpk_site);
        trPc_sd_site = permute(trPc_sd_site, [2,3,1]);
    end
end   
end %func


%--------------------------------------------------------------------------
function B = median3_(A,idimm)
if nargin<2, idimm=1; end

A = sort(A,idimm);
imid = ceil(size(A,idimm)/2);
if idimm==1
    B = A(imid,:,:);
elseif idimm==2
    B = A(:,imid,:);
else
    B = A(:,:,imid);
end
end %func


%--------------------------------------------------------------------------
function mat = cell2matfun_(varargin)
mat = cell2mat(cellfun(varargin{1}, varargin{2:end}, 'UniformOutput', 0));
end %func


%--------------------------------------------------------------------------
function mrDist_clu = wave_similarity_clu_(S0, P, mode_sim)
S_clu = S0.S_clu;
MIN_COUNT = P.min_count;
if nargin<3, mode_sim = []; end
if isempty(mode_sim), mode_sim = 4; end

% compute average waveforms by clusters
nClu = max(S_clu.viClu);
cviSpk_clu = arrayfun(@(x)find(S_clu.viClu==x), 1:nClu, 'UniformOutput', 0);
% compute average clu waveforms (sort by centers)
[ctrPc_clu, ctrPos_clu, cviSite_clu] = deal(cell(nClu,1));
dimm_pc = size(S0.trPc_spk); 
dimm_pc(2) = min(dimm_pc(2), P.nSites_fet);
% dimm_pc(1) = min(dimm_pc(1), P.nPcPerChan);
dimm_pos = size(S0.mrPos_spk);
switch 2
    case 2, fh_pow_ = @(x)log(x);
    case 1, fh_pow_ = @(x)x;
end
% compute average pc per cluster per site
% find centroid of the cluster
for iClu = 1:nClu
    viSpk1 = cviSpk_clu{iClu};
    viSite1 = S0.viSite_spk(viSpk1);
    [viSite1_uniq1, vn_uniq1, cviSpk1_uniq1] = unique_count_(viSite1);
    trPc_spk1 = S0.trPc_spk(1:dimm_pc(1), 1:dimm_pc(2), viSpk1);
    mrPos_spk1 = [S0.mrPos_spk(viSpk1,:), fh_pow_(S0.vrPow_spk(viSpk1))];
    [viSite1_uniq1, cviSpk1_uniq1] = ...
        multifun_(@(x)x(vn_uniq1>=MIN_COUNT), viSite1_uniq1, cviSpk1_uniq1);
    trPc1 = zeros(dimm_pc(1), dimm_pc(2), numel(cviSpk1_uniq1), 'single');
    trPos1 = zeros(dimm_pos(2)+1, 2, numel(cviSpk1_uniq1), 'single');    
    for iiSite1 = 1:numel(cviSpk1_uniq1)
        % compute average waveforms by peak sites
        viSpk11 = cviSpk1_uniq1{iiSite1};
        trPc1(:,:,iiSite1) = mean(trPc_spk1(:,:,viSpk11), 3);
        [sd_, med_] = sd_mad_(mrPos_spk1(viSpk11,:));
        trPos1(:,:,iiSite1) = [med_(:), sd_(:)];
    end
    [ctrPc_clu{iClu}, ctrPos_clu{iClu}, cviSite_clu{iClu}] = deal(trPc1, trPos1, viSite1_uniq1);
end
% construct arrays and sites
[viSite_all, viClu_all] = cvr2vr_vi_(cviSite_clu);
[trPc_all, trPos_all] = deal(cat(3, ctrPc_clu{:}), cat(3, ctrPos_clu{:}));

% compare cluster by cluster
mrPv_global = S0.mrPv_global(:,1:dimm_pc(1));
mr2vr_ = @(x)x(:);
tr2mr_ = @(x)reshape(x,[],size(x,3));
mrDist_clu = zeros(nClu);
norm_mean_ = @(x)sqrt(mean(x.^2));
norm_sum_ = @(x)sqrt(sum(x.^2));
normalize_mean_ = @(x)x ./ sqrt(mean(x.^2));
normalize_sum_ = @(x)x ./ sqrt(sum(x.^2));
for iClu1 = 1:nClu
    [trPc1, viSite1, trPos1] = deal(ctrPc_clu{iClu1}, cviSite_clu{iClu1}, ctrPos_clu{iClu1});
    vrDist_clu1 = zeros(nClu,1);
    for iiSite1 = 1:numel(viSite1)
        [mrPc1, mrPos1] = deal(trPc1(:,:,iiSite1), trPos1(:,:,iiSite1));
        vi_all2 = find(viSite_all == viSite1(iiSite1) & viClu_all ~= iClu1);
        if isempty(vi_all2), continue; end
        viClu2 = viClu_all(vi_all2);
        [trPc2, trPos2] = deal(trPc_all(:,:,vi_all2), trPos_all(:,:,vi_all2));
        switch mode_sim
            case 7 % centroid merging
                [mu1, sd1] = deal(mrPos1(:,1), mrPos1(:,2));
                mu2 = squeeze_(trPos2(:,1,:), 2);
                vrCorr12 = norm_mean_((mu1-mu2) ./ sd1);
                vrCorr12 = 1 - min(vrCorr12/10,1);
            case 6
                a = normalize_sum_(mr2vr_(pc2wav_(mrPv_global, mrPc1)));
                b = normalize_sum_(tr2mr_(pc2wav_(mrPv_global, trPc2)));
                c = a'*a;
                vrCorr12 = 1-abs(c-a'*b) ./ c;
            case 5 % rms after normalization
                a = normalize_sum_(mr2vr_(pc2wav_(mrPv_global, mrPc1)));
                b = normalize_sum_(tr2mr_(pc2wav_(mrPv_global, trPc2)));
                vrCorr12 = 1 - std(a-b, 1); % standardized rms
            case 4 % rms after normalization, .97 threshold
                a = normalize_sum_(mrPc1(:));
                b = normalize_sum_(reshape(trPc2, [], size(trPc2,3)));
                vrCorr12 = 1 - std(a-b, 1); % standardized rms
            case 3
                vrCorr12 = normalize_sum_(mrPc1(:))' * normalize_sum_(tr2mr_(trPc2));
            case 2
                vrCorr12 = 1-pdist2(mrPc1(:)', reshape(trPc2, [], size(trPc2,3))', 'cosine');
            case 1
                vrCorr12 = wavcor_(pc2wav_(mrPv_global, mrPc1), pc2wav_(mrPv_global, trPc2));
            otherwise, error('wave_similarity_clu_: invalid mode');
        end
        vrDist_clu1(viClu2) = max(vrDist_clu1(viClu2), vrCorr12(:));
    end
    mrDist_clu(:,iClu1) = vrDist_clu1;
end
end %func


%--------------------------------------------------------------------------
function [vr, vr_med] = sd_mad_(mr)
vr_med = median2_(mr);
vr = median2_(abs(mr - vr_med)) / .6745;
end % func


%--------------------------------------------------------------------------
function B = median2_(A,idimm)
if nargin<2, idimm=1; end

A = sort(A,idimm);
imid = ceil(size(A,idimm)/2);
if idimm==1
    B = A(imid,:);
else
    B = A(:,imid);
end
end %func


%--------------------------------------------------------------------------
function vrCorr12 = wavcor_pc_(mrPc1, trPc2, fMode, mrPc_sd1)
if numel(fMode)>1
    mrPv_global = fMode;
    fMode = 1;
end
if nargin<4, mrPc_sd1 = []; end

switch fMode 
    case 1 % correlation distance of waveform
        mrWav1 = mrPv_global * mrPc1;
        a_ = zscore_(mrWav1(:), 1);
        mrWav2 = mrPv_global * reshape(trPc2, size(trPc2,1), []);
        b_ = zscore_(reshape(mrWav2, [], size(trPc2,3)), 1);
        vrCorr12 = (a_'*b_) / size(mrWav1,1) / size(mrWav1,2);
        
    case 2 % RMS distance of PC
        vr1 = mrPc1(:);
        mrD = vr1 - reshape(trPc2,[],size(trPc2,3));
        vrCorr12 = sqrt(sum(mrD.^2) ./ sum(vr1.^2));
        vrCorr12 = 1 - min(vrCorr12/10,1);
        
    case 3 % RMS distance of PC
        vr1 = mrPc1(:);
        mrD = vr1 - reshape(trPc2,[],size(trPc2,3));
        vrCorr12 = 1 - min(sum(mrD.^2) ./ sum(vr1.^2), 1);
        
    case 4 % RMS distance of PC. same result as 2
        mrWav1 = mrPv_global * mrPc1;
        vr1 = mrWav1(:);
        mrWav2 = reshape(mrPv_global * reshape(trPc2, size(trPc2,1), []), [], size(trPc2,3));
        mrD = vr1 - mrWav2;
        vrCorr12 = sqrt(sum(mrD.^2) ./ sum(vr1.^2));
        vrCorr12 = 1 - min(vrCorr12,1);
        
    case 5
        mrWav1 = mrPv_global * mrPc1;
        a_ = zscore_(mrWav1(:), 1);
        mrWav2 = mrPv_global * reshape(trPc2, size(trPc2,1), []);
        b_ = zscore_(reshape(mrWav2, [], size(trPc2,3)), 1);
        vrCorr12 = (a_'*b_) / size(mrWav1,1) / size(mrWav1,2);
        
    case 6
        vr1 = mrPc1(:);
        mr1 = reshape(trPc2,[],size(trPc2,3));
        vr2 = vr1/sqrt(sum(vr1.^2));
        mr2 = mr1./sqrt(sum(mr1.^2));
        vrCorr12 = vr2' * mr2;
        
    case 7
        vr1 = mrPc1(:);
        mr1 = reshape(trPc2,[],size(trPc2,3));
        vrCorr1 = 1 - min(sum((vr1 - mr1).^2) ./ sum(vr1.^2), 1);
        vr2 = vr1/sqrt(sum(vr1.^2));
        mr2 = mr1./sqrt(sum(mr1.^2));
        vrCorr2 = vr2' * mr2;  
        vrCorr12 = (vrCorr1 + vrCorr2)/2;
        
    case 8
        vr1 = mrPc1(:);
        mrD = vr1 - reshape(trPc2,[],size(trPc2,3));
        vrCorr12 = 1 - min(sum(mrD.^2) ./ sum(mrPc_sd1(:).^2), 1);        
        
end %switch
end %func


%--------------------------------------------------------------------------
function vrCorr12 = wavcor_(mrWav1, trWav2, nShift)
% mrWav1 = meanSubt_(mrWav1);
% trWav2 = meanSubt_(trWav2);

if nargin<3, nShift = 2; end
% try time delay and match
switch 1
    case 1
        if nShift == 0
            a_ = zscore_(mrWav1(:));
            b_ = zscore_(reshape(trWav2, [], size(trWav2,3)),1);
            vrCorr12 = (a_'*b_) / size(mrWav1,1) / size(mrWav1,2);
        else
            vi0 = (nShift+1):(size(mrWav1,1)-nShift);
            b_ = zscore_(reshape(trWav2(vi0,:,:), [], size(trWav2,3)),1);
            mrCorr12 = zeros(nShift*2+1, size(trWav2,3));
            for iShift1 = 1:(nShift*2+1)
                a_ = mrWav1(vi0 + iShift1-1-nShift,:);
                mrCorr12(iShift1,:) = zscore_(a_(:))' * b_;
            end
            vrCorr12 = max(mrCorr12)' / numel(vi0) / size(mrWav1,2);
        end
    case 2
        nT = size(mrWav1,1);
        mrCorr12 = zeros(nShift*2+1, size(trWav2,3));
        [cvi1, cvi2] = shift_range_(nT, nShift);
        for iShift1 = 1:numel(cvi1)
            [vi1,vi2] = deal(cvi1{iShift1}, cvi2{iShift1});
            a_ = mrWav1(vi1,:);
            b_ = zscore_(reshape(trWav2(vi2,:,:), [], size(trWav2,3)),1);
            mrCorr12(iShift1,:) = zscore_(a_(:))' * b_ / numel(vi1);
        end
        vrCorr12 = max(mrCorr12)' / size(mrWav1,2);
end
end %func


%--------------------------------------------------------------------------
function out = pc2wav_(mrPv, mrPc)
switch ndims(mrPc)
    case 2, out = mrPv * mrPc;
    case 3        
        [nPc, nSites, nWav] = size(mrPc);
        out = reshape(mrPv * reshape(mrPc, nPc,[]), [], nSites, nWav);
end
end %func


%--------------------------------------------------------------------------
function [vr, vi] = cvr2vr_vi_(cvr)
vr = cell2mat_(cvr);
vn1 = cellfun(@(x)size(x,1), cvr);
vi = cell2mat_(arrayfun(@(x)repmat(x, vn1(x),1), 1:numel(cvr), 'UniformOutput', 0)');
end %func


%--------------------------------------------------------------------------
% todo: delayed execution, use parfeval
function S0 = sort_(S0, P)

% drift processing
fprintf('Clustering\n'); 
runtime_sort = tic;

[nPc_fet, nSites_fet, miSites, knn] = ...
    struct_get_(P, 'nPcPerChan', 'nSites_fet', 'miSites', 'knn');
nPc_fet = min(nPc_fet, size(S0.trPc_spk,1));
trPc_spk = S0.trPc_spk(1:nPc_fet, :, :);
nSites = max(S0.viSite_spk);
cviSpk_site = arrayfun(@(x)find(S0.viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
vnSpk_site = cellfun(@numel, cviSpk_site);
nSpk = numel(S0.viSite_spk);

% parfor loop
if get_set_(P, 'fParfor', 1)
    gcp_ = gcp();
else
    gcp_ = [];
end
if 0
    [S0.mrPos_spk, S0.vrPow_spk] = calc_pos_spk_(S0.trPc_spk, S0.viSite_spk, P);
end
S_drift = calc_drift_(S0, P);
[viDrift_spk, mlDrift] = struct_get_(S_drift, 'viDrift_spk', 'mlDrift');

% Calculate Rho
P_sort = struct_set_(P, 'mlDrift', mlDrift, 'fSort_drift', 1);
[vrRho, vrDelta] = deal(zeros(nSpk, 1, 'single'));
miKnn = zeros(knn, nSpk, 'int32');
fprintf('Calculating Rho\n\t'); t1=tic;
[cvrRho, cmiKnn, cvrDelta, cviNneigh] = deal(cell(nSites,1));
% send jobs
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, viDrift_spk, P_sort, iSite);
    if ~isempty(gcp_)
        vS_out(iSite) = parfeval(gcp_, ...
            @(x,y,z,a)rho_knn_(x,y,z,a, P_sort), 2, mrFet12, viSpk12, viDrift12, n1);
    else
        [cvrRho{iSite}, cmiKnn{iSite}, P_sort.fGpu] = ...
            rho_knn_(mrFet12, viSpk12, viDrift12, n1, P_sort);      
    end
end %for
% collect jobs
if ~isempty(gcp_)
    for iSite1 = 1:nSites
        [iSite, vrRho1, miKnn1] = fetchNext(vS_out);
        [cvrRho{iSite}, cmiKnn{iSite}] = deal(vrRho1, miKnn1);
    end
end
% assemble jobs
for iSite = 1:nSites     
    viSpk1 = cviSpk_site{iSite};
    if isempty(viSpk1), continue; end
    [vrRho(viSpk1), miKnn(:,viSpk1)] = deal(cvrRho{iSite}, cmiKnn{iSite});
end %for
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t1), P_sort.fGpu, P_sort.fParfor);


% Calculate Delta
fprintf('Calculating Delta\n\t'); t2=tic;
% send jobs
for iSite = 1:nSites
    [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, viDrift_spk, P_sort, iSite);
    vrRho12 = vrRho(viSpk12);
    if ~isempty(gcp_)
        vS_out(iSite) = parfeval(gcp_, ...
            @(x,y,z,a,b)delta_knn_(x,y,z,a,b, P_sort), 2, mrFet12, viSpk12, vrRho12, viDrift12, n1);
    else
        [cvrDelta{iSite}, cviNneigh{iSite}, P_sort.fGpu] = ...
            delta_knn_(mrFet12, viSpk12, vrRho12, viDrift12, n1, P_sort);
    end
end %for
% collect jobs
if ~isempty(gcp_)
    for iSite1 = 1:nSites
        [iSite, vrDelta1, viNneigh1] = fetchNext(vS_out);
        [cvrDelta{iSite}, cviNneigh{iSite}] = deal(vrDelta1, viNneigh1);
    end
end
% assemble jobs
for iSite = 1:nSites     
    viSpk1 = cviSpk_site{iSite};
    if isempty(viSpk1), continue; end
    [vrDelta(viSpk1), viNneigh(viSpk1)] = deal(cvrDelta{iSite}, cviNneigh{iSite});
end %for
vrRho = vrRho / max(vrRho) / 10;     % divide by 10 to be compatible with previous version displays
fprintf('\n\ttook %0.1fs (fGpu=%d, fParfor=%d)\n', toc(t2), P_sort.fGpu, P_sort.fParfor);

% output
runtime_sort = toc(runtime_sort);
[~, ordrho] = sort(vrRho, 'descend');
memory_sort = memory_matlab_();
S0.S_clu = struct('rho', vrRho, 'delta', vrDelta, 'ordrho', ordrho, 'nneigh', viNneigh, ...
    'P', P, 'miKnn', miKnn, 'S_drift', S_drift);
S0.runtime_sort = runtime_sort;
S0.memory_sort = memory_sort;
end %func


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1, fGpu] = delta_knn_(mrFet12, viSpk12, vrRho12, viDrift12, n1, P)
if isempty(mrFet12)
    [vrDelta1, viNneigh1] = deal([]);
    return;
end
switch get_set_(P, 'fSort_drift', 1)
    case 1
        [vrDelta1, viNneigh1, fGpu] = delta_drift_knn_(mrFet12, viDrift12, P.mlDrift, vrRho12, n1, P);
        viNneigh1 = viSpk12(viNneigh1);
    case 0
        SINGLE_INF = 3.402E+38;
        vrRho12 = vrRho12(:);
        [vi12, vi1] = deal(1:size(mrFet12,2), 1:n1);

        % do cuda
        [vrDelta1, viNneigh1, fGpu] = cuda_delta_knn_(mrFet12, vrRho12, vi12, vi1, P);

        vrDelta1 = gather_(vrDelta1) .* vrRho12(vi1);
        viNneigh1 = viSpk12(gather_(viNneigh1));
        viNan = find(isnan(vrDelta1) | isinf(vrDelta1));
        viNneigh1(viNan) = viNan;
        vrDelta1(viNan) = sqrt(SINGLE_INF);
end %switch
end %func


%--------------------------------------------------------------------------
function [vrRho1, miKnn1, fGpu] = rho_knn_(mrFet12, viSpk12, viDrift12, n1, P)
if isempty(mrFet12)
    [vrRho1, miKnn1] = deal([]);
    return;
end
switch get_set_(P, 'fSort_drift', 1)
    case 1
        [vrRho1, fGpu, miKnn1] = rho_drift_knn_(mrFet12, viDrift12, P.mlDrift, n1, P);  
        miKnn1 = viSpk12(miKnn1);
    case 0
        knn = get_set_(P, 'knn', 30);
        [vi12, vi1] = deal(1:size(mrFet12,2), 1:n1);
        [vrKnn1, fGpu, miKnn1] = cuda_knn_(mrFet12, vi12, vi1, P);
        vrRho1 = gather_(1./vrKnn1);

        n_ = size(miKnn1,1);
        if n_ == knn
            miKnn1(:,:) = gather_(miKnn1);
        else        
            miKnn1(1:n_,:) = gather_(miKnn1);
            miKnn1(n_+1:end,:) = repmat(gather_(miKnn1(end,:)), [knn-n_, 1]);
        end
        miKnn1 = viSpk12(miKnn1);    
end %switch
end %func


%--------------------------------------------------------------------------
function [mrFet12, viSpk12, viDrift12, n1, n2] = pc2fet_site_(trPc_spk, cviSpk_site, viDrift_spk, P, iSite)

if isempty(cviSpk_site{iSite})
    [mrFet12, viSpk12, viDrift12] = deal([]);
    [n1, n2] = deal(0); 
    return;
end
[nSites_fet, miSites] = struct_get_(P, 'nSites_fet', 'miSites');

mrFet1 = trPc_spk(:, 1:nSites_fet, cviSpk_site{iSite});
mrFet1 = reshape(mrFet1, [], size(mrFet1,3));
viSites1 = miSites(1:nSites_fet,iSite);

% find neighboring sites
[~, mi_] = ismember(miSites, viSites1);
mi_(:,iSite)=0; %exclude seif
[~,mi2] = sort(mi_);
viSites2 = find(sum(mi_>0) == nSites_fet); %neighbor sites
miSites2 = mi2(end-nSites_fet+1:end, viSites2); % sites to index

% extract fet and build distance matrix
switch 1
    case 2
        mrFet12 = mrFet1;
    case 1
        cmrFet2 = cell(1, numel(viSites2));
        for iiSite2 = 1:numel(viSites2)
            trFet2_ = trPc_spk(:, miSites2(:,iiSite2), cviSpk_site{viSites2(iiSite2)});
            cmrFet2{iiSite2} = reshape(trFet2_, [], size(trFet2_,3));
        end
        mrFet12 = [mrFet1, cell2mat_(cmrFet2)];
end
n1 = numel(cviSpk_site{iSite});
n2 = size(mrFet12,2) - n1;
viSpk12 = cell2mat(cviSpk_site([iSite, viSites2]));
if isempty(viDrift_spk)
    viDrift12 = [];
else
    viDrift12 = viDrift_spk(viSpk12);
end

% spatial mask
if 0
    if get_set_(P, 'fSpatialMask_clu', 1) && nSites_fet >= get_set_(P, 'nChans_min_car', 8)
        vrSpatialMask = spatialMask_(P, iSite, nSites_fet, P.maxDist_site_um);
        vrSpatialMask = repmat(vrSpatialMask(:)', [P.nPcPerChan, 1]);
        mrFet12 = bsxfun(@times, mrFet12, vrSpatialMask(:));
    end
end
end %func


%--------------------------------------------------------------------------
function S0 = detect_(P)
% keep all features in the memory, no disk storage

% parfor loop
if get_set_(P, 'fParfor', 1)
    gcp_ = gcp();
else
    gcp_ = [];
end

% P.fGpu = 1;
% else
%     gcp_ = [];
% end

% fParfor=0, fGpu=0: 90s
% fParfor=0, fGpu=1: 26s **GPU is fastest**
% fParfor=1, fGpu=0: 30s (local, 9 workers)
% fParfor=1, fGpu=1: 35s (local, 9 workers)
runtime_detect = tic; 
memory_init = memory_matlab_();

% load one
S_paged = readmda_paged_(P); % initialize
[nLoads, viOffset_load] = deal(S_paged.nLoads, S_paged.viOffset_load);
[mrWav1, nlim_wav1, fDone] = readmda_paged_(); % process first part
cS_detect = cell(nLoads, 1);
cS_detect{1} = detect_paged_(mrWav1, P, makeStruct_(nlim_wav1)); % process the first part
[vrThresh_site, mrPv_global] = struct_get_(cS_detect{1}, 'vrThresh_site', 'mrPv_global');
S_cache1 = makeStruct_(vrThresh_site, mrPv_global);
fprintf('Memory use: %0.3f GiB\n', memory_matlab_()/2^30);
for iLoad = 2:nLoads          
    [mrWav1, nlim_wav1, fDone] = readmda_paged_(); % process first part
    S_cache1.nlim_wav1 = nlim_wav1; % trim waveform
    if isempty(gcp_)            
        cS_detect{iLoad} = detect_paged_(mrWav1, P, S_cache1);            
    else
        vS_out(iLoad-1) = parfeval(gcp_, @(x,y)detect_paged_(x,P,y), 1, mrWav1, S_cache1);
    end    
    fprintf('\tMemory use: %0.3f GiB\n', memory_matlab_()/2^30);
end
if ~isempty(gcp_)
    for iLoad = 2:nLoads
        [completedIdx, S_] = fetchNext(vS_out);
        cS_detect{completedIdx+1} = S_;
    end
end
[viSite_spk, viTime_spk, vrAmp_spk, trPc_spk, mrVp_spk] = detect_merge_(cS_detect, viOffset_load);
switch 1
    case 2
        trPc_spk1 = reshape(mrVp_spk, 1, size(mrVp_spk,1), size(mrVp_spk,2));
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk1, viSite_spk, P);
    case 1
        [mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk, viSite_spk, P);
end %switch
% Save output
runtime_detect = toc(runtime_detect);
memory_detect = memory_matlab_();
S0 = makeStruct_(viSite_spk, viTime_spk, vrAmp_spk, trPc_spk, vrPow_spk, ...
    vrThresh_site, mrPv_global, runtime_detect, P, memory_detect, memory_init, mrPos_spk, mrVp_spk);
fprintf('detect_: took %0.1fs (fParfor=%d, fGpu=%d)\n', runtime_detect, P.fParfor, P.fGpu);
end %func


%--------------------------------------------------------------------------
function [mrPos_spk, vrPow_spk] = calc_pos_spk_(trPc_spk, viSite_spk, P)
[~, nSites_spk, nSpk] = size(trPc_spk);
switch 1
    case 4, mrW = squeeze_(abs(trPc_spk(1,:,:)),1);
    case 3, mrW = squeeze_(sum(trPc_spk(1:P.nPcPerChan,:,:).^2),1);
    case 2, mrW = squeeze_(sum(trPc_spk.^2),1);
    case 1, mrW = squeeze_(trPc_spk(1,:,:),1) .^2;
end
vrPow_spk = sum(mrW)';

fprintf('Calculating spike positions\n\t');
t1=tic;
toVec_ = @(x)x(:);
switch 1        
    case 1 % centroid       
        trSiteXY_spk = reshape(P.mrSiteXY(P.miSites(:, viSite_spk),:), nSites_spk, nSpk, 2);
        mrPos_spk = [sum(mrW .* trSiteXY_spk(:,:,1))', sum(mrW .* trSiteXY_spk(:,:,2))'] ./ vrPow_spk;
%         figure; plot(mrPos_spk(:,1), mrPos_spk(:,2), '.', 'MarkerSize',1); axis([10 60 10 210])
        % figure; plot3(mrPos_spk(:,1), mrPos_spk(:,2), log(mean(mrW)'), '.');
        
    case 2 % 2d peak (non-uniform kernel)
        mrPos_spk = zeros(nSpk, 2, 'single');
        nSites = max(viSite_spk);
        cviSpk_site = arrayfun(@(x)find(viSite_spk==x), 1:nSites, 'UniformOutput', 0)';
        for iSite = 1:nSites
            [viSpk1, viSite1] = deal(cviSpk_site{iSite}, P.miSites(:,iSite));
            [mrSiteXY1, mrW1] = deal(P.mrSiteXY(viSite1,:), mrW(:,viSpk1));
            mrPos_spk(viSpk1,:) = find_peak_2d_(mrSiteXY1, mrW1);
            fprintf('.');
        end
        
    case 3 % centroid by pc components
        [nPc, nSites_spk, nSpk] = size(trPc_spk);
        trPos_pc_spk = zeros(nSpk, 2, nPc, 'single');
        trSiteXY_spk = reshape(P.mrSiteXY(P.miSites(:, viSite_spk),:), nSites_spk, nSpk, 2);
        for iPc = 1:size(trPc_spk,1)
            mrW1 = squeeze_(trPc_spk(iPc,:,:),1).^2;
            mrPos1 = [sum(mrW1 .* trSiteXY_spk(:,:,1))', sum(mrW1 .* trSiteXY_spk(:,:,2))'] ./ sum(mrW1)';
            trPos_pc_spk(:,:,iPc) = mrPos1;
        end %for
        trPos_pc_spk = permute(trPos_pc_spk, [1,3,2]);        
end %switch
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function mrXY_spk1 = find_peak_2d_(mrXY1, mrW1)
% mrW1=double(mrW1);
% F = scatteredInterpolant(mrXY1(:,1), mrXY1(:,2), mrW1(:,1), 'natural');

% define grid
d0 = pdist(mrXY1); d0 = min(d0(d0>0));
dxy = 2; %4x interp
[xy_min, xy_max] = deal(min(mrXY1), max(mrXY1));
[xlim1, ylim1] = deal(...
        [xy_min(1)-d0, xy_max(1)+d0], ...
        [xy_min(2)-d0, xy_max(2)+d0]);
[xx,yy] = meshgrid(xlim1(1):dxy:xlim1(2), ylim1(1):dxy:ylim1(2));

% convolve and evaluate
mrXY2 = single([xx(:), yy(:)]);
A = exp(-(pdist2(mrXY1, mrXY2) / (d0)).^2); % 4 neighbors half
[~,viMax] = max(single(A)' * mrW1);
mrXY_spk1 = mrXY2(viMax,:);
mrXY_spk1 = mrXY_spk1 + (rand(size(mrXY_spk1))*2-1);

if 0
    mrW2 = reshape(mrW1(:,1)'*A, size(xx));
    figure; mesh(xx,yy,mrW2/10);
    hold on; plot3(mrXY1(:,1), mrXY1(:,2), mrW1(:,1), 'ro');
    hold on; plot3(mrXY_spk1(1,1), mrXY_spk1(1,2), mrW2(viMax(1)), 'r*');
    view(2)
end
end %func


%--------------------------------------------------------------------------
function [viSite_spk, viTime_spk, vrAmp_spk, trPc_spk, mrVp_spk] = detect_merge_(cS_detect, viOffset_load)

vnSpk_load = cellfun(@(x)numel(x.viSite_spk), cS_detect);
miSpk_load = [0; cumsum(vnSpk_load)];
miSpk_load = [miSpk_load(1:end-1)+1, miSpk_load(2:end)];

nSpk = sum(vnSpk_load);
[viSite_spk, viTime_spk, vrAmp_spk] = ...
    deal(zeros(nSpk, 1, 'int32'), zeros(nSpk, 1, 'int64'), zeros(nSpk, 1, 'single'));
viOffset_load = int64(viOffset_load);
trPc_spk = [];
mrVp_spk = [];

for iLoad = 1:numel(cS_detect)
    S1 = cS_detect{iLoad};
    viSpk1 = miSpk_load(iLoad,1):miSpk_load(iLoad,2);
    viSite_spk(viSpk1) = S1.viSite_spk;
    viTime_spk(viSpk1) = int64(S1.viTime_spk) + viOffset_load(iLoad);
    vrAmp_spk(viSpk1) = S1.vrAmp_spk;
    if isempty(trPc_spk)
        trPc_spk = zeros(size(S1.trPc_spk,1), size(S1.trPc_spk,2), nSpk, 'single');
    end
    trPc_spk(:,:,viSpk1) = S1.trPc_spk;
    if 0
        if isempty(mrVp_spk)
            mrVp_spk = zeros(size(S1.mrVp_spk,1), nSpk, 'single');
        end
        mrVp_spk(:,viSpk1) = S1.mrVp_spk;
    end
end
end %func


%--------------------------------------------------------------------------
function S_detect = detect_paged_(mrWav, P, S_cache)
if nargin<3, S_cache = []; end

% filter and trim 
% nlim_wav1 = struct_get_(S_cache, 'nlim_wav1');
mrWav2 = filter_(mrWav, P);
S_detect = get_spikes_(mrWav2, P, S_cache);

end %func


%--------------------------------------------------------------------------
function S_detect = get_spikes_(mrWav2, P, S_cache)
[vrThresh_site, nlim_wav1, mrPv_global] = ...
    struct_get_(S_cache, 'vrThresh_site', 'nlim_wav1', 'mrPv_global');

if isempty(nlim_wav1)
    [nPad_pre, nPad_post] = deal(0);
else
    nPad_pre = nlim_wav1(1)-1;
    nPad_post = size(mrWav2,1) - nlim_wav1(2);
end

% common mode rejection
if P.blank_thresh > 0
    if isempty(vnWav11)
        vnWav11 = mr2ref_(mrWav2, P.vcCommonRef, P.viSiteZero); %vrWav_mean1(:);    
    end
    vlKeep_ref = car_reject_(vnWav11(:), P);
    fprintf('Rejecting %0.3f %% of time due to motion\n', (1-mean(vlKeep_ref))*100 );
else
    vlKeep_ref = [];
end

% detect spikes or use the one passed from the input (importing)
if isempty(vrThresh_site), [vrThresh_site, fGpu] = mr2thresh_(mrWav2, P); end
[viTime_spk, vrAmp_spk, viSite_spk] = detect_spikes_(mrWav2, vrThresh_site, vlKeep_ref, P);
[viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)gather_(x), viTime_spk, vrAmp_spk, viSite_spk);    

% reject spikes within the overlap region
if ~isempty(nlim_wav1)
    viKeep_spk = find(viTime_spk >= nlim_wav1(1) & viTime_spk <= nlim_wav1(2));
    [viTime_spk, vrAmp_spk, viSite_spk] = multifun_(@(x)x(viKeep_spk), viTime_spk, vrAmp_spk, viSite_spk);    
end%if

% extract spike waveforms
trWav_spk = get_spkwav_(mrWav2, viSite_spk, viTime_spk, P);
if 0
    lim1 = round(size(trWav_spk,1) * [1/4,3/4]);
    mrVp_spk = min(trWav_spk(lim1(1):lim1(2),:,:));  
    mrVp_spk = gather_(squeeze_(mrVp_spk,1));
else
    mrVp_spk = [];
end
if 0
    dimm_spk = size(trWav_spk);
    trWav_spk = reshape(meanSubt_(reshape(trWav_spk, [], size(trWav_spk,3))), dimm_spk);
end

% extract spike feaures
if isempty(mrPv_global)
    [mrPv_global, vrD_global] = get_pv_(trWav_spk, P); 
end
trPc_spk = gather_(project_pc_(trWav_spk, mrPv_global));
%     trWav_spk1 = pc2spkwav_(trPc_spk, mrPv_global); % recover waveform

% return struct
if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end
S_detect = makeStruct_(trPc_spk, mrPv_global, viTime_spk, vrAmp_spk, viSite_spk, mrVp_spk);
end %func


%--------------------------------------------------------------------------
function trWav_spk1 = pc2spkwav_(trPc_spk, mrPv_global)
[nPc_spk, nSites_spk, nSpk] = size(trPc_spk);
nSamples_spk = size(mrPv_global,1);
dimm_spk = [nSamples_spk, nSites_spk, nSpk];
trWav_spk1 = reshape(mrPv_global * reshape(trPc_spk, size(trPc_spk,1), []), dimm_spk);
end %func


%--------------------------------------------------------------------------
function trPc_spk = project_pc_(trWav_spk, mrPv)
[nSamples, nSites, nSpk] = size(trWav_spk);
nPc = size(mrPv,2);
mr = reshape(trWav_spk, size(trWav_spk,1), []);
trPc_spk = reshape((mrPv' * mr), nPc, nSites, nSpk);
end %func


%--------------------------------------------------------------------------
function [mrPv1, vrD1, fh_proj] = get_pv_(tr, P)
%tr: nSamples x nSpikes x nChans
MAX_SAMPLE = 10000;        

viSpk_sub = subsample_vr_(1:size(tr,3), MAX_SAMPLE);
switch 2
    case 2, mr1 = reshape(tr(:, :, viSpk_sub), size(tr,1), []); 
    case 1, mr1 = squeeze(tr(:, 1, viSpk_sub)); 
end
nPc_spk = get_set_(P, 'nPc_spk', 6); % # components to compress spike waveforms

switch 2
    case 1
        % mrSpkWav1 = meanSubt_(mrSpkWav1);
        [mrPv1, vrD1] = eig(mr1 * mr1');
        mrPv1 = zscore_(fliplr(mrPv1)); % sort largest first
        vrD1 = flipud(diag(vrD1));
    case 2
        [mrPv1, ~, vrD1] = pca(gather_(mr1)','Centered',false, 'NumComponents', nPc_spk);
%         fh_proj = @(x)(mrPv' * x)';
%         mr2 = mrPv * mrPc'; % confirm the agreement
%         mr3 = mrPv * (mrPv' * mr1);
end
% spike center should be negative
% iMid = 1-P.spkLim(1);
% vrSign = (mrPv1(iMid,:) < 0) * 2 - 1; %1 or -1 depending on the sign
% mrPv = bsxfun(@times, mrPv1, vrSign);
% trFft = fft(tr);
% trFft = abs(trFft(2:end/2,:,:));
% figure; imagesc(std(trFft,1,3))
end %func


%--------------------------------------------------------------------------
function trWav_spk = get_spkwav_(mrWav, viSite_spk, viTime_spk, P)
nSpks = numel(viTime_spk);
nSites = numel(P.viSite2Chan);
spkLim_wav = P.spkLim;
% spkLim_raw = P.spkLim_raw;
nSites_spk = size(P.miSites,1); %(P.maxSite * 2) + 1;

% Realignment parameters
fRealign_spk = get_set_(P, 'fRealign_spk', 0); %0,1,2
if isempty(viSite_spk)
    trWav_spk = permute(mr2tr_(mrWav, spkLim_wav, viTime_spk), [1,3,2]);
else
    trWav_spk = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'like', mrWav);
    for iSite = 1:nSites
        viiSpk11 = find(viSite_spk == iSite);
        if isempty(viiSpk11), continue; end
        viTime_spk11 = viTime_spk(viiSpk11); %already sorted by time
        viSite11 = P.miSites(:,iSite);
        try
            tnWav_spk1 = mr2tr_(mrWav, spkLim_wav, viTime_spk11, viSite11);
            trWav_spk(:,:,viiSpk11) = permute(tnWav_spk1, [1,3,2]);
        catch % GPU failure
            disperr_('mn2tn_wav_: GPU failed'); 
        end
    end
end
% if 0 %10/19/2018 JJJ
%     trWav_spk = meanSubt_spk_(trWav_spk);
% end
end %func


%--------------------------------------------------------------------------
function mnWav2 = filter_(mnWav1, P)
%-----
% Filter
fprintf('\tFiltering spikes...'); t_filter = tic;
if get_set_(P, 'fSmooth_spatial', 0)
    mnWav1 = spatial_smooth_(mnWav1, P);
end
vcDataType_filter = get_set_(P, 'vcDataType_filter', 'single');
try    
    mnWav1_ = cast_(mnWav1, vcDataType_filter);
    [mnWav1_, P.fGpu] = gpuArray_(mnWav1_, P.fGpu);
    [mnWav2, vnWav11] = filt_car_(mnWav1_, P);    
catch % GPU failure
    P.fGpu = 0;
    mnWav1_ = cast_(mnWav1, vcDataType_filter);
    [mnWav2, vnWav11] = filt_car_(mnWav1_, P);    
end
fprintf(' took %0.1fs\n', toc(t_filter));
end %func


%--------------------------------------------------------------------------
function [out, nlim, fDone] = readmda_paged_(P)
% manage overlap
persistent S_mda fid iLoad nLoads nSamples_load nSamples_last nSamples_pad vcFile

if nargin==1
    vcFile = P.vcFile;
    [S_mda, fid] = readmda_header_(vcFile);
    S_mda.dimm = S_mda.dimm(:)'; % formatting
    assert(S_mda.nBytes_missing<1, 'readmda_paged_: mda file is incomplete');
    iLoad = 0;
    [nLoads, nSamples_load, nSamples_last] = plan_load_(S_mda.nBytes_data, P);
    nSamples_pad = get_set_(P, 'nPad_filt', 100);
    viOffset_load = (0:nLoads-1) * nSamples_load; % return offset
    S_page = makeStruct_(nLoads, nSamples_load, nSamples_last, nSamples_pad, viOffset_load);
    [out, nlim, fDone] = deal(S_page, [], 0);    
else
    if iLoad <= nLoads
        iLoad = iLoad + 1;
        if nLoads == 1 % no padding if single load
            dimm1 = S_mda.dimm;
            nlim = [1, dimm1(end)]; 
        elseif iLoad == 1 % first
            dimm1 = [S_mda.dimm(1:end-1), nSamples_load + nSamples_pad];
            nlim = [1, nSamples_load];
        elseif iLoad == nLoads % last
            dimm1 = [S_mda.dimm(1:end-1), nSamples_last + nSamples_pad];
            nlim = [1, nSamples_last] + nSamples_pad;
        else % middle (double padded)
            dimm1 = [S_mda.dimm(1:end-1), nSamples_load + nSamples_pad*2];
            nlim = [1, nSamples_load] + nSamples_pad;
        end
        t1=tic;
        out = fread_(fid, dimm1, S_mda.vcDataType);  % transpose MDA
        t_dur1 = toc(t1);
        out = out'; % transpose
        mb_loaded1 = prod(dimm1) * bytesPerSample_(S_mda.vcDataType) / 1e6;
        fprintf('Read %s (%d/%d), took %0.1fs (%0.1f MB/s, %0.1f MB)\n', ...
            vcFile, iLoad, nLoads, t_dur1, mb_loaded1/t_dur1, mb_loaded1);
        if iLoad == nLoads
            fclose(fid);
            fDone = 1;
        else
            fDone = 0;
            frewind_(fid, [S_mda.dimm(1:end-1), nSamples_pad*2], S_mda.vcDataType);
        end
    else
        [out, nlim] = deal([]);
        fDone = 1;
    end
end

end %func


%--------------------------------------------------------------------------
function vcDir_in = get_test_data_(vcMode)
if nargin<1, vcMode = 'static'; end

if ispc()
    vcDir_in = 'C:\tmp\groundtruth\hybrid_synth\%s_siprobe\rec_64c_1200s_11'; 
elseif isunix()
    vcDir_in = '~/ceph/groundtruth/hybrid_synth/%s_siprobe/rec_64c_1200s_11'; 
end
vcDir_in = strrep(vcDir_in, '%s', vcMode);
end %func


%--------------------------------------------------------------------------
function P = makeParam_(vcDir_in, vcDir_out, vcFile_arg)
if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_arg = ''; end

% assume there is raw.mda, geom.csv, params.json, firings_true.mda
P = file2struct_(ircpath_(read_cfg_('default_prm', 0)));  %P = defaultParam();
P.vcFile = fullfile(vcDir_in, 'raw.mda');
P.vcDir_out = vcDir_out;
S_mda = readmda_header_(P.vcFile);
P.nChans = S_mda.dimm(1);
P.vnSamples_file = S_mda.dimm(2); % scalar

% probe file
P.mrSiteXY = csvread(fullfile(vcDir_in, 'geom.csv'));
P.viSite2Chan = 1:size(P.mrSiteXY,1);

% load json file
S_json = loadjson_(fullfile(vcDir_in, 'params.json'));
P.sRateHz = get_set_(S_json, 'samplerate', P.sRateHz);
P.fInverse_file = get_set_(S_json, 'spike_sign', -1) == -1;

% read param
if isstruct(vcFile_arg)
    S_arg = vcFile_arg;
elseif exist_file_(vcFile_arg)
    S_arg = meta2struct_(vcFile_arg);
else
    S_arg = [];
end
P = struct_merge_(P, S_arg);

% derived fields
P = fill_param_(P);
end %func


%--------------------------------------------------------------------------
function P = fill_param_(P)
P.spkRefrac = round(P.spkRefrac_ms * P.sRateHz / 1000);
P.spkLim = round(P.spkLim_ms * P.sRateHz / 1000);
% P.spkLim_raw = calc_spkLim_raw_(P);
if ~isempty(get_(P, 'viChanZero')) && isempty(P.viSiteZero)
    [~, viSiteZero] = ismember(P.viChanZero, P.viSite2Chan);
    P.viSiteZero = viSiteZero(viSiteZero>0);
end
if ~isempty(get_(P, 'viSiteZero')), P.viSiteZero(P.viSiteZero > numel(P.viSite2Chan)) = []; end
if ~isfield(P, 'viShank_site'), P.viShank_site = []; end
[P.miSites, P.nSites_fet] = findNearSites_(P.mrSiteXY, P);
% LFP sampling rate
if ~isempty(get_(P, 'nSkip_lfp'))
    P.sRateHz_lfp = P.sRateHz / P.nSkip_lfp;
else
    P.sRateHz_lfp = get_set_(P, 'sRateHz_lfp', 2500);
    P.nSkip_lfp = round(P.sRateHz / P.sRateHz_lfp);
end
if isempty(get_(P, 'nTime_clu'))
    try
        batch_sec_drift = get_set_(P, 'batch_sec_drift', 300);
        P.nTime_clu = max(round(recording_duration_(P) / batch_sec_drift), 1);
        P.nTime_clu = min(P.nTime_clu, get_set_(P, 'nBatch_max_drift', 32));
        fprintf('\tnTime_clu = %d (batch_sec_drift = %0.1f s)\n', P.nTime_clu, batch_sec_drift);
    catch
        P.nTime_clu = 1;
    end
end
if isempty(get_(P, 'nTime_drift'))
    try
        step_sec_drift = get_set_(P, 'step_sec_drift', 10);
        P.nTime_drift = max(round(recording_duration_(P) / step_sec_drift), 1);
        fprintf('\tnTime_drift = %d (step_sec_drift = %0.1f s)\n', P.nTime_drift, step_sec_drift);
    catch
        P.nTime_drift = 64;
    end
end

P.bytesPerSample = bytesPerSample_(P.vcDataType);
P = struct_default_(P, 'vcFile_prm', subsFileExt_(P.vcFile, '.prm'));
% if ~isempty(get_(P, 'gain_boost')), P.uV_per_bit = P.uV_per_bit / P.gain_boost; end
P.spkThresh = P.spkThresh_uV / P.uV_per_bit;
P = struct_default_(P, 'cvrDepth_drift', {});
% P = struct_default_(P, {'maxSite_fet', 'maxSite_detect', 'maxSite_sort','maxSite_pix', 'maxSite_dip', 'maxSite_merge', 'maxSite_show'}, P.maxSite);
P = struct_default_(P, 'mrColor_proj', [.75 .75 .75; 0 0 0; 1 0 0]);
P.mrColor_proj = reshape(P.mrColor_proj(:), [], 3); %backward compatible
P = struct_default_(P, {'blank_thresh', 'thresh_corr_bad_site', 'tlim_load'}, []);
if numel(P.tlim)==1, P.tlim = [0, P.tlim]; end
if isfield(P, 'rejectSpk_mean_thresh'), P.blank_thresh = P.rejectSpk_mean_thresh; end
P.vcFilter = get_filter_(P);
if isempty(get_(P, 'vcFilter_show'))
    P.vcFilter_show = P.vcFilter;
end
P.nC_max = read_cfg_('nC_max'); % override nC_max (gpu parameter)
end %func


%--------------------------------------------------------------------------
% from irc.m (v4.9.11)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 8/22/18 JJJ: changed from the cell output to varargout
% 9/26/17 JJJ: Created and tested
function varargout = struct_get_(varargin)
% Obtain a member of struct
% cvr = cell(size(varargin));
% if varargin is given as cell output is also cell
S = varargin{1};
for iArg=1:nargout
    vcName = varargin{iArg+1};
    if iscell(vcName)
        csName_ = vcName;
        cell_ = cell(size(csName_));
        for iCell = 1:numel(csName_)
            vcName_ = csName_{iCell};
            if isfield(S, vcName_)
                cell_{iCell} = S.(vcName_);
            end
        end %for
        varargout{iArg} = cell_;
    elseif ischar(vcName)
        if isfield(S, vcName)
            varargout{iArg} = S.(vcName);
        else
            varargout{iArg} = [];
        end
    else
        varargout{iArg} = [];
    end
end %for
end %func


%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
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


%--------------------------------------------------------------------------
% 1/31/2019 JJJ: get the field(s) of a struct or index of an array or cell
function varargout = get_(varargin)
% same as struct_get_ function
% retrieve a field. if not exist then return empty
% [val1, val2] = get_(S, field1, field2, ...)
% [val] = get_(cell, index)

if nargin==0, varargout{1} = []; return; end
S = varargin{1};
if isempty(S), varargout{1} = []; return; end

if isstruct(S)
    for i=2:nargin
        vcField = varargin{i};
        try
            varargout{i-1} = S.(vcField);
        catch
            varargout{i-1} = [];
        end
    end
elseif iscell(S)
    try    
        varargout{1} = S{varargin{2:end}};
    catch
        varargout{1} = [];
    end
else
    try    
        varargout{1} = S(varargin{2:end});
    catch
        varargout{1} = [];
    end
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: added '.' if dir is empty
% 7/31/17 JJJ: Substitute file extension
function varargout = subsFileExt_(vcFile, varargin)
% Substitute the extension part of the file
% [out1, out2, ..] = subsFileExt_(filename, ext1, ext2, ...)
if isempty(vcFile), varargout{1}=[]; return; end
[vcDir_, vcFile_, ~] = fileparts(vcFile);
if isempty(vcDir_), vcDir_ = '.'; end
for i=1:numel(varargin)
    vcExt_ = varargin{i};    
    varargout{i} = [vcDir_, filesep(), vcFile_, vcExt_];
end
end %func


%--------------------------------------------------------------------------
function [vi_uniq, vn_uniq, cvi_uniq] = unique_count_(vi)
% count number of unique elements and sort by vn_uniq
vi_uniq = unique(vi);
cvi_uniq = arrayfun(@(x)find(vi==x), vi_uniq, 'UniformOutput', 0);
vn_uniq = cellfun(@numel, cvi_uniq);
[vn_uniq, ix] = sort(vn_uniq, 'descend');
vi_uniq = vi_uniq(ix);
cvi_uniq = cvi_uniq(ix);
end %func


%--------------------------------------------------------------------------
% viSite1_uniq1, vn_uniq1, cviSpk1_uniq1
function [viUnique, vnUnique, cvi_uniq] = unique_count__(vi)
% vn = [1 1 1 1 4 4 4 8];
% vn = [1 1 1 1];
[viUnique, vnUnique, cvi_uniq] = deal([], [], {});
if isempty(vi), return; end
vi1 = sort(vi(:)');
vi_ = find(diff(sort(vi1)) > 0);
if isempty(vi_)
    vnUnique = 1;
    viUnique = vi1(1);
else
    vnUnique = diff([0, vi_, numel(vi)]);
    viUnique = [vi1(1), vi1(vi_+1)];
end
if nargout >= 3
    cvi_uniq = arrayfun(@(x)find(vi==x), viUnique, 'UniformOutput', 0);
end
end %func


%--------------------------------------------------------------------------
function z = zscore_(x, flag, dim)
if isempty(x), z=[]; return; end
if nargin < 2, flag = 0; end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
mu = mean(x,dim);
sigma = std(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
end %func


%--------------------------------------------------------------------------
function frewind_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function disperr_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function struct_save_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end

function out1 = meta2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_merge_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ircpath_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = read_cfg_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = loadjson_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = filesize_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = bytesPerSample_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = fread_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gather_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2ref_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = car_reject_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_copy_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cast_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2tr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsample_vr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cell2mat_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_default_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_filter_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = set0_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = postCluster_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gt2mda_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_dir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = meanSubt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
% function out1 = zscore_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = templateMatch_post_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = memory_matlab_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = calc_drift_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = recording_duration_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = squeeze_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refresh_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end

function [out1, out2] = readmda_header_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = mr2thresh_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = gpuArray_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = filt_car_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = findNearSites_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = spatialMask_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = shift_range_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end

function [out1, out2, out3] = plan_load_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = detect_spikes_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = cuda_delta_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = cuda_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = unique_count_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = rho_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = delta_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end



