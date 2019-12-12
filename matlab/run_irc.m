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
version = 2; % default version
t_run_irc = tic;
if ~isdeployed()
    source_path = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(source_path)));
end

if strcmpi(vcDir_in, 'version')
    vcCmd = vcDir_in;
    switch version
        case 1, fprintf('%s\n', irc(vcCmd)); 
        case 2, fprintf('%s\n', irc2(vcCmd)); 
    end
    exit_deployed_();
    return; 
end

vcDir_out = fill_dir_out_(vcDir_in, vcDir_out);
 
% inferred from the path
firings_out_fname = fullfile(vcDir_out, 'firings.mda');
raw_fname = fullfile(vcDir_in, 'raw.mda');
if ~exist_file_(raw_fname)
    fprintf(2, 'file does not exist: %s\n', raw_fname);
    return;
end
fForceRerun = irc('call', 'read_cfg', {'fForceRerun'});
if ~exist_file_(firings_out_fname) || fForceRerun    
    switch version
        case 2
            fParfor = 0;
            irc2(vcDir_in, vcDir_out, vcFile_template, fParfor);
            
        case 1            
            geom_fname = fullfile(vcDir_in, 'geom.csv');  
            if exist_file_(fullfile(vcDir_in, 'argfile.txt'))
                prm_fname = fullfile(vcDir_in, 'argfile.txt');
                version = irc('read-param', prm_fname, 'version');
                if isempty(version), version = 1; end
            else
                prm_fname = fullfile(vcDir_in, 'params.json');
            end
            vcFile_prm = irc('makeprm-mda', raw_fname, geom_fname, prm_fname, vcDir_out, vcFile_template);    
            irc('clear', vcFile_prm);
            irc('run', vcFile_prm);
            irc('export-mda', vcFile_prm, firings_out_fname);
            
            % create a ground truth
            vcFile_gt_mda = fullfile(vcDir_in, 'firings_true.mda');
            vcFile_score = fullfile(vcDir_out, 'raw_geom_score.mat');
            if exist_file_(vcFile_gt_mda)  && ~exist_file_(vcFile_score)
                try
                    irc('validate-mda', vcFile_gt_mda, firings_out_fname, raw_fname); % assume that groundtruth file exists
                catch
                    fprintf(2, 'Validation failed\n');
                end
            end
            fClear_mda = irc('call', 'read_cfg', {'fClear_mda'});
            if fClear_mda % clear temp files after exporting mda file
                irc('clear', fullfile(vcDir_out, 'raw_geom.prm')); % clear output
            end
            fprintf('Clustering result wrote to %s\n', firings_out_fname);
    end %switch    
end
fprintf('#SF-SORTER-RUNTIME#%0.3f#\n', toc(t_run_irc)); % internal time keeping
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


%--------------------------------------------------------------------------
function vcDir_out = fill_dir_out_(vcDir_in, vcDir_out)
if isempty(vcDir_out)
    if ~contains(vcDir_out, 'groundtruth')
        vcDir_out = fullfile(vcDir_in, 'irc2');
    else
        vcDir_out = strrep(vcDir_in, 'groundtruth', 'irc2'); 
    end
end
if ~exist_dir_(vcDir_out), mkdir(vcDir_out); end
end %func


%--------------------------------------------------------------------------
% 8/7/2018 JJJ
function flag = exist_dir_(vcDir)
if isempty(vcDir)
    flag = 0;
else
    S_dir = dir(vcDir);
    if isempty(S_dir)
        flag = 0;
    else
        flag = sum([S_dir.isdir]) > 0;
    end
%     flag = exist(vcDir, 'dir') == 7;
end
end %func