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

if nargin<2, vcDir_out = ''; end
if nargin<3, vcFile_template = ''; end
if isempty(vcDir_out)
    vcDir_out = strrep(vcDir_in, '/groundtruth/', '/irc/');
end
if ~isdeployed()
    source_path = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(source_path)));
end
if strcmpi(vcDir_in, 'version')
    vcCmd = vcDir_in;
    fprintf('%s\n', irc(vcCmd)); 
    return; 
end

irc('call', 'mkdir', {vcDir_out}); % create temp output directory
    
% inferred from the path
firings_out_fname = fullfile(vcDir_out, 'firings_out.mda');
raw_fname = fullfile(vcDir_in, 'raw.mda');
if ~exist_file_(raw_fname)
    fprintf(2, 'file does not exist: %s\n', raw_fname);
    return;
end
fForceRerun = irc('call', 'read_cfg', {'fForceRerun'});
fClear_mda = irc('call', 'read_cfg', {'fClear_mda'});
if ~exist_file_(firings_out_fname) || fForceRerun
    geom_fname = fullfile(vcDir_in, 'geom.csv');
    prm_fname = fullfile(vcDir_in, 'params.json');
    vcFile_prm = irc('makeprm-mda', raw_fname, geom_fname, prm_fname, vcDir_out, vcFile_template);    
    irc('clear', vcFile_prm);
    irc('run', vcFile_prm);
    irc('export-mda', vcFile_prm, firings_out_fname);
    fprintf('Clustering result wrote to %s\n', firings_out_fname);
    if fClear_mda % clear temp files after exporting mda file
        irc('clear', vcFile_prm); % clear output
    end
end

% create a ground truth
vcFile_gt_mda = fullfile(vcDir_in, 'firings_true.mda');
if exist_file_(vcFile_gt_mda)
    try
        vcFile_prm = fullfile(vcDir_out, 'raw_geom.prm');
        irc('import-gt', vcFile_gt_mda, vcFile_prm);
        irc('validate-mda', vcFile_gt_mda, firings_out_fname, raw_fname); % assume that groundtruth file exists
    catch
        fprintf(2, 'Validation failed\n');
    end
end

% Exit
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

