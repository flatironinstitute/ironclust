ret = p_ironclust('spec')

%%

params = struct( ...
    'timeseries', 'D:\Dropbox\flatiron\ironclust\processors\raw.mda', ...
    'nfilt', '2');
if isunix()
    params.timeseries = '/home/jamesjun/Dropbox/flatiron/ironclust/processors/raw.mda';
end
params.timeseries_out = strrep(params.timeseries, 'raw.mda', 'filt.mda');

p_slopefilter(params);

return;

%%
P = call('file2struct', 'default.prm');

%%
params = struct('paramfile', 'raw_full.prm', 'timeseries', 'raw.mda', 'geom', 'geom_2A.csv');

p_ironclust(params);

%% jrclust to mountainsort converter
% geom.prb to geom.csv
call('prb2geom', 'geom_2A_.prb')

call('get_set', struct('a', 1), 'a', 2)
call('get_set', struct('a', 1), 'b', 2)
call('get_set', struct('a', []), 'b', 2)
call('get_set', struct('a', []), 'a', 2)