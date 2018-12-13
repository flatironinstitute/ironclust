function run_irc(vcDir_in, vcDir_out, vcFile_template)
% usage
% -----
% run_irc(command)
% run_irc(vcDir_in, vcDir_out, vcFile_template)
%
% arguments
% -----
% vcDir_in: input directory
% vcDir_out: output directory
% vcFile_template: template file (optional)

if nargin<3, vcFile_template = ''; end

if ~isdeployed()
    source_path = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(source_path, '/matlab'))); 
end
if nargin==1
    vcCmd = vcDir_in;
    fprintf('%s\n', irc(vcCmd)); 
    return; 
end

irc('call', 'mkdir', {vcDir_out}); % create temp output directory
    
% inferred from the path
raw_fname = fullfile(vcDir_in, 'raw.mda');
geom_fname = fullfile(vcDir_in, 'geom.csv');
firings_out_fname = fullfile(vcDir_out, 'firings_out.mda');
prm_fname = fullfile(vcDir_in, 'params.json');
vcFile_gt_mda = fullfile(vcDir_in, 'firings_true.mda');

vcFile_prm = irc('makeprm-mda', raw_fname, geom_fname, prm_fname, vcDir_out, vcFile_template);
irc('clear', vcFile_prm);
irc('run', vcFile_prm);
irc('export-mda', vcFile_prm, firings_out_fname);

% create a ground truth
try
    irc('import-gt', vcFile_gt_mda, vcFile_prm); % assume that groundtruth file exists
    irc('validate', vcFile_prm); % assume that groundtruth file exists
catch
    fprintf(2, 'Validation failed\n');
end

% Exit
fprintf('Clustering result wrote to %s\n', firings_out_fname);
exit_deployed_();
end %func


%--------------------------------------------------------------------------
function exit_deployed_()
try
    if isdeployed() || ismcc(), exit(); end
catch
    ;
end
end %func


%--------------------------------------------------------------------------
% 7/21/2018 JJJ: rejecting directories, strictly search for flies
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile, fVerbose)
if nargin<2, fVerbose = 0; end
if isempty(vcFile)
    flag = 0; 
else
    S_dir = dir(vcFile);
    if numel(S_dir) == 1
        flag = ~S_dir.isdir;
    else
        flag = 0;
    end
end
if fVerbose && ~flag
    fprintf(2, 'File does not exist: %s\n', vcFile);
end
end %func

