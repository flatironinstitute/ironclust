function p_ironclust(temp_path, raw_fname, geom_fname, prm_fname, firings_out_fname, arg_fname)
% cmdstr2 = sprintf("p_ironclust('$(tempdir)','$timeseries$','$geom$','$firings_out$','$(argfile)');");

% temp_path='/tmp/test_jrclust'
if exist(temp_path, 'dir') ~= 7
    mkdir(temp_path);
end
    
subs_dir_(raw_fname, temp_path);

copyfile_(raw_fname, subsdir_(raw_fname)); [temp_path, '/raw.mda']);
copyfile_(prb_fname, [temp_path, filesep(), prb_fname]);
copyfile_(param_fname,[temp_path,'/raw.meta']);
copyfile_(settings_fname,[temp_path,'/settings.prm']);

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

end %func


%--------------------------------------------------------------------------
function copyfile_(vcFrom, vcDest)
try
    copyfile(vcFrom, vcDest, 'f');
catch
    ;
end
end %func


%--------------------------------------------------------------------------
function copyfile_(vcFrom, vcDest)

end %func