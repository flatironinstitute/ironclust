function p_ironclust(temp_path, raw_fname, geom_fname, prm_fname, vcFile_gt_mda, firings_out_fname, arg_fname)
% cmdstr2 = sprintf("p_ironclust('$(tempdir)','$timeseries$','$geom$','$firings_out$','$(argfile)');");
t1=tic;
if exist(temp_path, 'dir') ~= 7
    mkdir(temp_path);
end
    
version = irc('read-param', arg_fname, 'version');
if isempty(version), version=2; end

disp('===================================================');
fprintf('IronClust Version: %d\n', version);
disp('===================================================');

switch version
    case 1
        vcFile_prm = irc('makeprm-mda', raw_fname, geom_fname, arg_fname, temp_path, prm_fname);
        irc('clear', vcFile_prm); %init 
        irc('run', vcFile_prm);
        irc('export-mda', vcFile_prm, firings_out_fname);
    case 2
        irc2(fileparts(raw_fname), temp_path, arg_fname);
end %switch

% create a ground truth
if exist_file_(vcFile_gt_mda)
    irc('import-gt', vcFile_gt_mda, vcFile_prm); % assume that groundtruth file exists
    irc('validate', vcFile_prm); % assume that groundtruth file exists
end

fprintf('Clustering result wrote to %s\n', firings_out_fname);
fprintf('#SF-SORTER-RUNTIME#%0.3f#\n', toc(t1)); % internal time keeping
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
% 9/26/17 JJJ: Created and tested
function vcFile_new = subs_dir_(vcFile, vcDir_new)
% Substitute dir
[vcDir_new,~,~] = fileparts(vcDir_new);
[~, vcFile, vcExt] = fileparts(vcFile);
vcFile_new = fullfile(vcDir_new, [vcFile, vcExt]);
end % func


%--------------------------------------------------------------------------
% 8/14/18 JJJ: Created and tested
function vcFile_full = subs_file_(vcFile, vcFile_new)
% Substitute dir
[vcDir_new,~,~] = fileparts(vcFile);
[~, vcFile_new1, vcFile_new2] = fileparts(vcFile_new);
vcFile_full = fullfile(vcDir_new, [vcFile_new1, vcFile_new2]);
end % func


%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
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

