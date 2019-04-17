%--------------------------------------------------------------------------
function P = file2struct(vcFile_file2struct)
% James Jun 2017 May 23
% Run a text file as .m script and result saved to a struct P
% _prm and _prb can now be called .prm and .prb files
P = []; 
if ~exist_file_(vcFile_file2struct), return; end
    
% load text file. trim and line break. remove comments.  replace 
csLines_file2struct = file2lines_(vcFile_file2struct);
csLines_file2struct = strip_comments_(csLines_file2struct);
if isempty(csLines_file2struct), return; end

P = struct();
try
    eval(cell2mat(csLines_file2struct'));

    S_ws = whos(); 
    csVars = {S_ws.name};
    csVars = setdiff(csVars, {'csLines_file2struct', 'vcFile_file2struct', 'P'});
    for i=1:numel(csVars)
        eval(sprintf('a = %s;', csVars{i}));
        P.(csVars{i}) = a;
    end
catch
    fprintf(2, 'Error in %s:\n\t', vcFile_file2struct);
    fprintf(2, '%s\n', lasterr());
    P=[];
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


%--------------------------------------------------------------------------
% Read a text file and output cell strings separated by new lines
% 7/24/17 JJJ: Code cleanup
function csLines = file2lines_(vcFile_file2struct)
csLines = {};
if ~exist_file_(vcFile_file2struct, 1), return; end

fid = fopen(vcFile_file2struct, 'r');
csLines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

csLines = csLines{1};
end %func


%--------------------------------------------------------------------------
% Strip comments from cell string
% 7/24/17 JJJ: Code cleanup
function csLines = strip_comments_(csLines)
csLines = csLines(cellfun(@(x)~isempty(x), csLines));
csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
csLines = csLines(cellfun(@(x)x(1)~='%', csLines));

% remove comments in the middle
for i=1:numel(csLines)
    vcLine1 = csLines{i};
    iComment = find(vcLine1=='%', 1, 'first');
    if ~isempty(iComment)
        vcLine1 = vcLine1(1:iComment-1);
    end
    vcLine1 = strrep(vcLine1, '...', '');
    if ismember(strsplit(vcLine1), {'for', 'end', 'if'})
        csLines{i} = [strtrim(vcLine1), ', ']; %add blank at the end
    else
        csLines{i} = [strtrim(vcLine1), ' ']; %add blank at the end
    end
end
% csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
csLines = csLines(cellfun(@(x)~isempty(x), csLines));
end %func