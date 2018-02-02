function ret=p_slopefilter(params)

if (ischar(params))&&(strcmp(params,'spec'))
    ret=get_spec();
    return;
end;

X_path=params.timeseries;
Y_path=params.timeseries_out;
nfilt=params.nfilt;

disp('Reading...');
X=readmda(X_path); % Later we will do things in chunks

disp('Filtering...');
Y=slopefilter(X,nfilt);

disp('Writing...');
writemda32(Y,Y_path);

disp('Done with slopefilter');

end % func


%-------------------------------
function spec=get_spec

name='ironclust.slopefilter';
version='0.1';

inputs={};
inputs{end+1}=struct('name','timeseries');

outputs={};
outputs{end+1}=struct('name','timeseries_out');

parameters={};
parameters{end+1}=struct('name','nfilt');

opts=struct;

spec=struct( ...
    'name',name, ...
    'version',version, ...
    'inputs',{inputs}, ...
    'outputs',{outputs}, ...
    'parameters',{parameters}, ...
    'opts',opts ...
); 

end %func