function p_jrclust(temp_path,raw_fname,prb_fname,param_fname,settings_fname,firings_out_fname)

% temp_path='/tmp/test_jrclust'

copyfile(raw_fname,[temp_path,'/raw.mda'],'f');
copyfile(prb_fname,[temp_path,'/raw.prb'],'f');
copyfile(param_fname,[temp_path,'/raw.meta'],'f');
copyfile(settings_fname,[temp_path,'/settings.prm'],'f');

jrc('makeprm',[temp_path,'/raw.mda'],[temp_path,'/raw.prb'],[temp_path,'/settings.prm']); %creates tmp_jrclust/raw_raw.prm
jrc3('call', 'edit_prm_file_', {struct('header_offset', 20, 'fSavePlot_RD', 0), [temp_path,'/raw_raw.prm']}, 0);

jrc('detectsort',[temp_path,'/raw_raw.prm']);
jrc('export-csv',[temp_path,'/raw_raw.prm']); % writes to tmp_jrclust/raw_raw.csv

disp('======================================================================');

firings=csvread([temp_path,'/raw_raw.csv']);

disp('======================================================================');
disp('Shape of firings:');
disp(size(firings));

firings=firings(:,[3,1,2])';
S_out = jrc3('call', 'file2struct_', {[temp_path, '/raw_raw.prm']}, 1);
P = S_out.out1;
firings(2,:)=firings(2,:) * P.sRateHz;

writemda(firings,firings_out_fname);
