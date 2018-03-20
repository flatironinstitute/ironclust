function jrclust_spec

curpath=fileparts(mfilename('fullpath'));
addpath([curpath,'/jsonlab-1.5']);

pspec=struct;
pspec.name='test.jrclust';
pspec.version='0.1.7';
pspec.description='';
inputs0={};
inputs0{end+1}=struct('name','timeseries');
inputs0{end+1}=struct('name','prb');
inputs0{end+1}=struct('name','meta');
inputs0{end+1}=struct('name','settings');
pspec.inputs=inputs0;
outputs0={};
outputs0{end+1}=struct('name','firings_out');
pspec.outputs=outputs0;
parameters0={};
pspec.parameters=parameters0;
cmdstr=sprintf("addpath('%s'); addpath('%s/JRCLUST'); addpath('%s/mdaio'); try p_jrclust('$(tempdir)','$timeseries$','$prb$','$meta$','$settings$','$firings_out$'); catch, fprintf(2,lasterr()); exit(-1); end; quit;",curpath,curpath,curpath);
exe_command=sprintf('matlab -nodisplay -r "%s"',cmdstr);
pspec.exe_command=exe_command;

processors={};
processors{end+1}=pspec;
spec=struct('processors',{processors});

json=savejson('',spec);
fprintf('%s\n',json);

