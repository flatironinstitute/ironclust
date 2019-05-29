%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function S = meta2struct(vcFile)
% Convert text file to struct
S = struct();
if ~exist_file_(vcFile, 1), return; end

fid = fopen(vcFile, 'r');
mcFileMeta = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', false);
fclose(fid);
csName = mcFileMeta{1};
csValue = mcFileMeta{2};
for i=1:numel(csName)
    vcName1 = csName{i};
    if vcName1(1) == '~', vcName1(1) = []; end
    try         
        eval(sprintf('%s = ''%s'';', vcName1, csValue{i}));
        eval(sprintf('num = str2double(%s);', vcName1));
        if ~isnan(num)
            eval(sprintf('%s = num;', vcName1));
        end
        eval(sprintf('S = setfield(S, ''%s'', %s);', vcName1, vcName1));
    catch
        fprintf('%s = %s error\n', csName{i}, csValue{i});
    end
end
end %func


%--------------------------------------------------------------------------
% 7/21/2018 JJJ: rejecting directories, strictly search for flies
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile, fVerbose)
if nargin<2, fVerbose = 0; end
if isempty(vcFile)
    flag = false; 
elseif iscell(vcFile)
    flag = cellfun(@(x)exist_file_(x, fVerbose), vcFile);
    return;
else
    S_dir = dir(vcFile);
    if numel(S_dir) == 1
        flag = ~S_dir.isdir;
    else
        flag = false;
    end
end
if fVerbose && ~flag
    fprintf(2, 'File does not exist: %s\n', vcFile);
end
end %func


