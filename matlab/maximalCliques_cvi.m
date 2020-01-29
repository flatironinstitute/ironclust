function [ cvi_out ] = maximalCliques_cvi( cvi_in )
% Modified by James Jun on 2020 Jan 28
%   Accept cell of vector of index (cvi) instead of binary matrix (to big)
%   Language updated for Matlab R2019b
%   jamesjun@gmail.com

%MAXIMALCLIQUES Find maximal cliques using the Bron-Kerbosch algorithm
%   Given a graph's boolean adjacency matrix, A, find all maximal cliques 
%   on A using the Bron-Kerbosch algorithm in a recursive manner.  The 
%   graph is required to be undirected and must contain no self-edges.
%
%   V_STR ='v2' used (faster)
%
%   MC is the output matrix that contains the maximal cliques in its 
%   columns.
%
%   Note: This function can be used to compute the maximal independent sets
%   of a graph A by providing the complement of A as the input graph.  
%
%   Note: This function can be used to compute the maximal matchings of a 
%   graph A by providing the complement of the line graph of A as the input
%   graph.
%
%   Ref: Bron, Coen and Kerbosch, Joep, "Algorithm 457: finding all cliques
%   of an undirected graph", Communications of the ACM, vol. 16, no. 9, 
%   pp: 575–577, September 1973.
%
%   Ref: Cazals, F. and Karande, C., "A note on the problem of reporting 
%   maximal cliques", Theoretical Computer Science (Elsevier), vol. 407,
%   no. 1-3, pp: 564-568, November 2008.
%
%   Jeffrey Wildman (c) 2011
%   jeffrey.wildman@gmail.com

% test mode
if nargin==0 && nargout==0
    test_(); return;
end


% make it a row shape and int32
cvi_in = cellfun(@(x)unique_(int32(x(:)')), cvi_in(:)', 'UniformOutput', 0); 
P = int32(1):int32(numel(cvi_in));      
% remove self to satisfy no self edges
vr2cell_ = @(vr)arrayfun(@(x)x, P, 'UniformOutput', 0);
cvi_in = cellfun(@(x,y)x(x~=y), cvi_in, vr2cell_(P), 'UniformOutput', 0);
% make sure cvi_in is diagonal?

% output: storage for maximal cliques
cvi_out = {};   

intersect_ = @(a,b)a(ismember_(a,b));
setdiff_ = @(a,b)a(~ismember_(a,b));
union_ = @(a,b)unique([a,b]);

BKv2_(int32([]), P, int32([]));


% version 2 of the Bron-Kerbosch algo
function BKv2_(R, P, X)
% P : prospective nodes connected to all nodes in R
% R : currently growing clique
% X : nodes already processed

if isempty(P) && isempty(X)           
    cvi_out{end+1} = unique_(R); % report R as a maximal clique
else % choose pivot
    ppivots = union_(P, X); % potential pivots
    [~, ind] = max(cellfun(@(x)sum(ismember_(x,P)), cvi_in(ppivots)));            
    vi_u = setdiff_(P, cvi_in{ppivots(ind)}); % select one of the ppivots with the largest count
    vi_u = vi_u(:)'; % make sure it's a row vector
    for u = vi_u
        vi_ = find(P==u);
        if isempty(vi_)
            P(vi_) = [];
        else
            P = [P, u];
        end
%         P = setxor_(P,u);
        Nu = cvi_in{u};
        BKv2_([R, u], intersect_(P,Nu), intersect_(X,Nu));
        X = [X, u];
    end
end % if
end % bkv2       
end % maximalCliques


%--------------------------------------------------------------------------
function vl = ismember_(a,b)
vl = builtin('_ismemberhelper',a,b);
end %func


%--------------------------------------------------------------------------
function c = setxor_(a,b)
c = unique([a,b]); % union
d = a(ismember_(a,b)); %intersect
c = c(~ismember_(c,d));
end %func


%--------------------------------------------------------------------------
function vr = unique_(vr)
vr = sort(vr);
vi = find(diff(vr)>0);
if numel(vi)+1 < numel(vr)
    vr = vr([1, vi]);
end
end %func


%--------------------------------------------------------------------------
function test_()
% ml_in = logical([0 1; 1 0]);
ml_in = logical([0 1 0; 1 0 0; 0 0 0]);
% ml_in = logical([0 1 1; 1 0 0; 1 0 0]);
cvi_in = arrayfun(@(x)find(ml_in(:,x)), 1:size(ml_in,2), 'UniformOutput', 0);
disp('maximalCliques() output:')
b = maximalCliques(ml_in); arrayfun(@(x)disp(find(b(:,x))'), 1:size(b,2));

disp('maximalCliques_cvi() output:')
b = maximalCliques_cvi(cvi_in); cellfun(@(x)disp(x), b);
end %func