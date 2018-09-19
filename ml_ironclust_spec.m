function ml_ironclust_spec


curpath = fileparts(mfilename('fullpath'));
curpath = fullfile(curpath, filesep());
addpath([curpath,'jsonlab-1.5']);

S_json = loadjson('ironclust_spec.json');
  
% version
% S_json.processors{1}.version = '0.1';
      
% command string
% exe_command: '[needs to be wrapped in matlab call syntax] p_ironclust_sort('$timeseries$','$geom$','$prb$','$firings_out$','$(argfile)');'
cmdstr1 = sprintf("addpath('%s'); addpath('%smatlab'); addpath('%smdaio');", curpath,curpath,curpath); 
cmdstr2 = sprintf("p_ironclust('$(tempdir)','$timeseries$','$geom$','$prm$','$firings_true$','$firings_out$','$(argfile)');");
cmdstr = sprintf("%s %s quit;", cmdstr1, cmdstr2);
S_json.processors{1}.exe_command = sprintf('matlab -nosplash -nodisplay -r "%s"',cmdstr);

json=savejson('', S_json);
fprintf('%s\n',json);

