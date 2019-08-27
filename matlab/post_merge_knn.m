%--------------------------------------------------------------------------
% 7/2/2019 JJJ: Merge using knn relationship, not using waveforms (fast)
function S_clu = post_merge_knn(S_clu, P)

% waveform based merging. find clusters within maxSite
% also removes duplicate spikes
if nargin<3, fPostCluster=1; end

nRepeat_merge = get_set_(P, 'nRepeat_merge', 10);
% refresh clu, start with fundamentals
S_clu = struct_copy_(S_clu, 'rho', 'delta', 'ordrho', 'nneigh', 'P', ...
    't_runtime', 'halo', 'viiSpk', 'trFet_dim', 'vrDc2_site', 'miKnn', ...
    'viClu', 'icl', 'S_drift');

if fPostCluster, S_clu = postCluster_(S_clu, P); end

S_clu = S_clu_refresh_(S_clu);
S_clu.viClu_premerge = S_clu.viClu;

switch get_set_(P, 'post_merge_mode', 1) %1 previously 1
    case 16, S_clu = post_merge_local_(S_clu, P);
    case 15, S_clu = post_merge_similarity_cc_(S_clu, P);
    case 14, S_clu = post_merge_cc_(templateMatch_post_(S_clu, P), P);
    case 13, S_clu = templateMatch_post_burst_(S_clu, P);
    case 12, S_clu = post_merge_knnwav(S_clu, get0_('viSite_spk'), P);
    case 11, S_clu = post_merge_knn1(S_clu, P);
    case 10, S_clu = post_merge_knn1(templateMatch_post_(post_merge_knn1(S_clu, P), P), P);
    case 9, S_clu = post_merge_knn1(templateMatch_post_(S_clu, P), P);
    case 8, S_clu = templateMatch_post_(post_merge_knn1(S_clu, P), P);
    case 7, S_clu = post_merge_drift_(post_merge_drift_(S_clu, P), P);
    case 6, S_clu = post_merge_drift_(S_clu, P);
    case 5, S_clu = drift_merge_post_(S_clu, P);
    case 4, S_clu = graph_merge_(S_clu, P);
    case 3
        for iRepeat = 1:2
            S_clu = driftMatch_post_(S_clu, P); 
        end
    case 2, S_clu = featureMatch_post_(S_clu, P);
    case 1, S_clu = templateMatch_post_(S_clu, P);
    otherwise, fprintf('No post-merge is performed\n'); 
end %switch

S_clu = post_merge_wav_(S_clu, 0, P);
S_clu = S_clu_sort_(S_clu, 'viSite_clu');
S_clu = S_clu_refrac_(S_clu, P); % refractory violation removal
S_clu = S_clu_update_wav_(S_clu, P);
S_clu.mrCC = correlogram_(S_clu, get0('viTime_spk'), P);

% set diagonal element
[S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, S_clu_self_corr_(S_clu, [], S0));
S_clu.P = P;
S_clu = S_clu_position_(S_clu);
S_clu.csNote_clu = cell(S_clu.nClu, 1);  %reset note
S_clu = S_clu_quality_(S_clu, P);
[S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
end %func


%--------------------------------------------------------------------------
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ml2map_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refresh_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
