ret = p_slopefilter('spec');
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
mr_raw = readmda(params.timeseries);
mr_filt = readmda(params.timeseries_out);

figure; hold on; plot(mr_raw(1,:)); plot(mr_filt(1,:));

%%
S_mda = readmda_header(params.timeseries);

%% 
FF = writemda_header(params.timeseries_out, S_mda);
nBytes_header = ftell(FF)
fclose(FF);

S_mda_out = readmda_header(params.timeseries_out);

%%
mr_filt_6 = readmda(params.timeseries_out);
mr_filt_1 = readmda(params.timeseries_out);
max(max(abs(mr_filt_1-mr_filt_6)))