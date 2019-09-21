%% measure read overhead

%vcDir_in = 'C:\tmp\hybrid_synth\static_siprobe\rec_64c_1200s_11';
vcDir_in = 'C:\tmp\hybrid_synth\drift_siprobe\rec_64c_1200s_11';
vcFile_raw = fullfile(vcDir_in, 'raw.mda');

%%
t1=tic;
readmda_ = @(x)irc('call','readmda', {x}, 1); % 4.22 s, 10x faster
% readmda_ = @(x)readmda(x); % 49.367929
a = readmda_(vcFile_raw);
toc(t1)

%% benchmark disk IO speed
t2=tic
run_irc(vcDir_in);
toc(t2)
% Runtime (s):
%     Detect + feature (s):   69.2s
%     Cluster + merge (s):    57.9s
%       Cluster runtime (s):  18.7s
%       merge runtime (s):    39.2s
%     Total runtime (s):      127.1s
%     Runtime speed           x9.4 realtime


%% goal: run four threads simultaneously 

%----
% 0. create a param file
P = loadParam_('default.prm');
P.vcFile = vcFile_raw;

%-----
% 1. detect loop: divide by time
% To request multiple evaluations, use a loop.
read_par_('init', vcFile_raw, P); % detect threshold and eigenvalue
mrWav1 = read_par_('next');
[S0, S_out0] = detect_init_(mrWav1, P); % determine threshold and basis function
if S0.nLoads == 1
    S0 = struct_add_(S0, S_out0); % add the initial output
else
    % Launch a distributed task
    p = gcp();
    for iLoad = 1:S0.nLoads      
        vS_out(iLoad) = parfeval(p, @detect_par_, mrWav1, S0); % Square size determined by idx    
        mrWav1 = read_par_('next');
    end

    % Collect the results as they become available.
    cell_detect = cell(1,nLoads);
    for iLoad = 1:nLoads
      % fetchNext blocks until next results are available.
      [completedIdx, S_out1] = fetchNext(vS_out);
      cell_detect{completedIdx} = S_out1;
      fprintf('Got result with index: %d.\n', completedIdx);
    end
end

%-----
% 2. sort loop: divide by channels
% build knn-graph
for iSite = 1:nSites
    
end %for


%% make a smart use of multithread (parfeval). now it's divided to 9 parts to process 4.5GB

%% parfeval loop is utilized 