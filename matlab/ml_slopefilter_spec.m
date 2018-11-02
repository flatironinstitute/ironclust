function ml_slopefilter_spec

curpath=fileparts(mfilename('fullpath'));
addpath([curpath,'/processors']);
addpath([curpath,'/jsonlab-1.5']);
addpath([curpath,'/matlab']);

processors={};

pfunctions={...
    'p_slopefilter',...
};

for ii=1:length(pfunctions)
    func_name=pfunctions{ii};
    func=str2func(func_name);
    spec0=func('spec');
    spec0.exe_command=create_exe_command(func_name,spec0,curpath);
    processors{end+1}=spec0;
end

spec=struct('processors',{processors});

json=savejson('',spec);
disp(json);

function exe_command=create_exe_command(func_name,spec0,curpath)

pnames={};
for ii=1:length(spec0.inputs)
    pnames{end+1}=spec0.inputs{ii}.name;
end
for ii=1:length(spec0.outputs)
    pnames{end+1}=spec0.outputs{ii}.name;
end
for ii=1:length(spec0.parameters)
    pnames{end+1}=spec0.parameters{ii}.name;
end

struct_string='';
for ii=1:length(pnames)
    pname=pnames{ii};
    if (ii>1)
        struct_string=[struct_string,','];
    end
    struct_string=[struct_string,"'",pname,"','$",pname,"$'"];
end

str=sprintf("addpath('%s/processors'); %s(struct(%s));",curpath,func_name,struct_string);
exe_command=sprintf('octave --eval "%s"',str);
