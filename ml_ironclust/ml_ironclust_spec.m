function ml_ironclust_spec


curpath = fileparts(mfilename('fullpath'));
curpath = fullfile(curpath, filesep());
addpath([curpath,'/jsonlab-1.5']);

S_json = loadjson('ironclust_spec.json');
  
% version
% S_json.processors{1}.version = '0.1';
      
% command string
% exe_command: '[needs to be wrapped in matlab call syntax] p_ironclust_sort('$timeseries$','$geom$','$prb$','$firings_out$','$(argfile)');'
cmdstr1 = sprintf("addpath('%s'); addpath('%sJRCLUST'); addpath('%smdaio');", curpath,curpath,curpath); 
cmdstr2 = sprintf("p_ironclust('$(tempdir)','$timeseries$','$geom$','$firings_out$','$(argfile)');");
cmdstr = sprintf("%s try %s catch, fprintf(2,lasterr()); exit(-1); end; quit;", cmdstr1, cmdstr2);
S_json.processors{1}.exe_command = sprintf('matlab -nodisplay -r "%s"',cmdstr);

% pspec=struct;
% pspec.name='test.jrclust';
% pspec.version='0.11';
% pspec.description='';
% inputs0={};
% inputs0{end+1}=struct('name','timeseries');
% inputs0{end+1}=struct('name','prb');
% inputs0{end+1}=struct('name','prb');
% pspec.inputs=inputs0;
% outputs0={};
% outputs0{end+1}=struct('name','firings_out');
% pspec.outputs=outputs0;
% parameters0={};
% pspec.parameters=parameters0;
% %cmdstr=sprintf("addpath('%s'); addpath('%s/JRCLUST'); addpath('%s/mdaio'); try p_jrclust('$(tempdir)','$timeseries$','$prb$','$meta$','$firings_out$'); catch, fprintf(2,lasterr()); exit(-1); end; quit;",curpath,curpath,curpath);
% cmdstr=sprintf("addpath('%s'); addpath('%s/JRCLUST'); addpath('%s/mdaio'); try p_jrclust('$(tempdir)','$timeseries$','$prb$','$(argfile)','$firings_out$'); catch, fprintf(2,lasterr()); exit(-1); end; quit;",curpath,curpath,curpath);
% exe_command=sprintf('matlab -nodisplay -r "%s"',cmdstr);
% pspec.exe_command=exe_command;
% 
% processors={};
% processors{end+1}=pspec;
% spec=struct('processors',{	});

json=savejson('', S_json);
fprintf('%s\n',json);

