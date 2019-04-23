%--------------------------------------------------------------------------
function irc_scoreboard()
irc('addpath');
hMsg = msgbox('Loading files');
[csRecording, csStudy, csStudyset] = load_dataset_();
[csParam] = load_param_();

try, close(hMsg); catch; end

% create GUI
hFig = irc('call', 'create_figure_', {'Fig_scoreboard', [0 0 .5 1], 'IronClust Scoreboard', 0, 1}); 
fh_ax1 = @(x)axes('Parent', hFig, 'OuterPosition', x);

S_fig = set_userdata_(hFig, csRecording);
S_fig = set_userdata_(hFig, csStudy);
S_fig = set_userdata_(hFig, csStudyset);
S_fig = set_userdata_(hFig, csParam);


% create three axes
vhAx = cellfun(@(x)fh_ax1(x), {[0 .1 .15 .3], [.15 .1 .15 .3], [.35 .1 .15 .3], [.5 .1 .15 .3], [.70 .1 .15 .3], [.85 .1 .15 .3]});
[hAx_accu1, hAx_accu2, hAx_prec1, hAx_prec2, hAx_recl1, hAx_recl2] = deal_(vhAx);
xlabel_(vhAx(1:2:end), 'SNR (Vmin/Vrms)');
xlabel_(vhAx(2:2:end), 'Param Rank');
ylabel_(vhAx(1:2), 'Accuracy');
ylabel_(vhAx(3:4), 'Precision');
ylabel_(vhAx(5:6), 'Recall');
arrayfun(@(x)grid(x, 'on'), vhAx);
arrayfun(@(x)ylim(x, [0 1]), vhAx);
S_fig = set_userdata_(hFig, vhAx);
set(hFig, 'UserData', S_fig);

hText1 = uicontrol(hFig, 'Style','text','String','',...
       'Unit', 'Normalized', 'OuterPosition', [.5 .4 .5 .3], ...
       'HorizontalAlignment', 'left', 'BackgroundColor','w', 'Tag', 'text_param');
hText1.String = {'Select a parameter'}; % display meta
S_fig = set_userdata_(hFig, hText1);

hText2  = uicontrol(hFig, 'Style','text','String','',...
       'Unit', 'Normalized', 'OuterPosition', [0 0 1 .1], ...
       'HorizontalAlignment', 'left', 'BackgroundColor','w', 'Tag', 'text_selection');   
hText2.String = {'All selected'}; % display meta
S_fig = set_userdata_(hFig, hText2);

% create table
hTable_studyset = create_table_(hFig, csStudyset, 'Studyset', [0 .8 .5 .2]);
hTable_study = create_table_(hFig, csStudy, 'Study', [0 .6 .5 .2]);
hTable_recording = create_table_(hFig, csRecording, 'Recording', [0 .4 .5 .2]);
hTable_param = create_table_(hFig, csParam, 'Param', [.5 .7 .5 .3]);

% disable when one is selected, click again to enable others
vhTable = [hTable_studyset, hTable_study, hTable_recording, hTable_param];
S_fig = set_userdata_(hFig, vhTable);
arrayfun(@(x)set(x, 'CellSelectionCallback', @cbf_select_table_), vhTable);

% populate table
add_menu_(hFig);
end %func


%--------------------------------------------------------------------------
function run_batch_(h, e)
hFig = h.Parent.Parent;
csInput1 = get_userdata_(hFig, 'csInput1');
csParam1 = get_userdata_(hFig, 'csParam1');
if isempty(csParam1), csParam1 = get_userdata_(hFig, 'csParam'); end
if isempty(csInput1), csInput1 = get_userdata_(hFig, 'csStudyset'); end

sbatch_mda_(csInput1, csParam1);
end %func


%--------------------------------------------------------------------------
function delete_file_(vc)
try
    delete(vc)
catch
end
end %func


%--------------------------------------------------------------------------
% apply paramters to sorting
function sbatch_mda_(csInput, csParam)
% non-blocking call
if ischar(csParam), csParam = {csParam}; end
if ischar(csInput), csInput = {csInput}; end
vcDir_out = locate_('output');
vcDir_in = locate_('input');
fWait = 1;
nParams = numel(csParam);
% run all param in parallel. create a GT if doesn't exist and run
% validation
fh_param = @(x)fullfile(vcDir_out, 'param', x);
for iData = 1:numel(csInput)    
    vcInput1 = csInput{iData};
    fprintf('Processing %s\n', vcInput1); t1=tic;
    csDir_rec1 = locate_('recordings', vcInput1);    
    [csDir_in1, csDir_out1, csFile_param1] = deal({});
    for iRec = 1:numel(csDir_rec1)
        vcDir_rec1 = csDir_rec1{iRec};
        vcRecording1 = strrep(fileparts(vcDir_rec1), vcDir_in, '');
        for iParam = 1:nParams
            vcDir_out1 = fullfile(vcDir_out, csParam{iParam}, vcRecording1);
            [csDir_in1{end+1}, csDir_out1{end+1}, csFile_param1{end+1}] = ...
                deal(vcDir_rec1, vcDir_out1, fh_param(csParam{iParam}));
        end
    end
    irc('sbatch-mda', csDir_in1, csDir_out1, csFile_param1, fWait);
    fprintf('\tFinished %s, took %0.1fs\n', vcInput1, toc(t1));
end

end %func


%--------------------------------------------------------------------------
function vcPath = locate_(vcType, vcKey)
if nargin<2, vcKey = []; end
vcDir_in = irc('call', 'read_cfg', {'path_groundtruth'});
vcDir_out = irc('call', 'read_cfg', {'path_validation'});
vcDir_out = fullfile(vcDir_out, irc('version'));

switch lower(vcType)
    case {'recordings', 'rec'} % return a cell of directories containing raw.mda        
        [~, vcPath] = find_files_(locate_('input', vcKey), '/**/raw.mda');
        
    case 'cache'
        vcPath = fullfile(vcDir_out, 'scoreboard.mat');
    case 'settings'
        vcPath = ircpath_('settings.scoreboard');
    case 'param'
        if isempty(vcKey)
            vcPath = find_files_(vcDir_out, '/param/*_template.prm');
        elseif ischar(vcKey)
            vcPath = fullfile(vcDir_out, 'param', [vcKey, '_template.prm']);
        elseif iscell(vcKey)
            vcPath = cellfun_(@(x)fullfile(vcDir_out, 'param', [x, '_template.prm']), vcKey);
        end
    case {'input', 'in'}
        if ~isempty(vcKey)
            if ischar(vcKey)
                vcPath = fullfile(vcDir_in, vcKey);
            elseif iscell(vcKey)
                vcPath = cellfun_(@(x)fullfile(vcDir_in, x), vcKey);
            end
        else
            vcPath = vcDir_in;
        end
    case {'output', 'out'}
        if ~isempty(vcKey)
            if ischar(vcKey)
                vcPath = fullfile(vcDir_out, vcKey);
            else
                vcPath = cellfun_(@(x)fullfile(vcDir_out, x), vcKey);
            end
        else
            vcPath = vcDir_out;
        end
end %switch
end %func


%--------------------------------------------------------------------------
function varargout = cellfun_(varargin)
if nargout==1
    varargout{1} = cellfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==2
    [varargout{1}, varargout{2}] = cellfun(varargin{:}, 'UniformOutput', 0);    
elseif nargout==3
    [varargout{1}, varargout{2}, varargout{3}] = cellfun(varargin{:}, 'UniformOutput', 0);    
else
    error('cellfun_: nargout exceeds 3');
end   
end %func


%--------------------------------------------------------------------------
function add_menu_(hFig)
drawnow;

set(hFig, 'MenuBar','None');
mh_file = uimenu(hFig,'Label','File'); 
uimenu(mh_file,'Label', 'Edit `settings.scoreboard`', 'Callback', @(h,e)edit_('settings.scoreboard'));
uimenu(mh_file,'Label', 'Edit `default.cfg`', 'Callback', @(h,e)edit_('default.cfg'));

mh_file = uimenu(hFig,'Label','Run'); 
uimenu(mh_file,'Label', 'Generate param files (hashcode_template.prm)', 'Callback', @generate_param_);
uimenu(mh_file,'Label', 'Run batch', 'Callback', @run_batch_);
end %func


%--------------------------------------------------------------------------
% generate *_template.prm files from settings.scoreboard
function generate_param_(h,e)
% load settings.scoreboard
S0 = irc('call', 'file2struct', {locate_('settings')});

% generate struct arrays
csName = fieldnames(S0);
csVal = cellfun_(@(x)S0.(x), csName);
vnParam = zeros(size(csVal));
for iVal = 1:numel(csVal)
    val = csVal{iVal};
    if ~iscell(val), csVal{iVal} = {val}; end
    vnParam(iVal) = numel(csVal{iVal});
end
vcDir_out = locate_('output', 'param');
nParams = prod(vnParam);
vcAns = questdlg(sprintf('Generate %d parameter set?', nParams), 'Confirmation');
if isempty(vcAns), return; end
if ~strcmpi(vcAns, 'Yes'), return; end

hMsg = msgbox('Deleting previous parameters...');
irc('call', 'mkdir', {vcDir_out});
delete_file_(fullfile(vcDir_out, '*_template.prm'));
close_(hMsg);

hMsg = msgbox('Writing files...');
for iParam = 1:nParams
    S1 = new_struct_(csName, csVal, ind2sub_(vnParam, iParam));    
    vcHash1 = DataHash(S1);
    vcFile1 = fullfile(vcDir_out, [vcHash1, '_template.prm']);
    struct2file_(S1, vcFile1);
    fprintf('\tWrote to %s\n', vcFile1);
end
close_(hMsg);

% update parameter list
update_table_('Param', load_param_());
end %func


%--------------------------------------------------------------------------
function S1 = new_struct_(csName, csVal, viVal)
S1 = struct();
for iName = 1:numel(csName)
    [iVal_, name_, val_] = deal(viVal(iName), csName{iName}, csVal{iName});
    try
        S1.(name_) = val_{iVal_};
    catch
        error('new_struct_: invalid index %s(%d)', name_, iVal_);
    end
end
end %func


%--------------------------------------------------------------------------
function vi_sub = ind2sub_(siz, ndx)
assert(prod(siz) >= ndx && ndx>0, 'index out of range');
nout = numel(siz);
siz = double(siz);
lensiz = length(siz);
vi_sub = zeros(numel(siz),1);
k = cumprod(siz);
for i = nout:-1:2,
    vi = rem(ndx-1, k(i-1)) + 1;
    vj = (ndx - vi)/k(i-1) + 1;
    vi_sub(i) = double(vj);
    ndx = vi;
end
vi_sub(1) = ndx;
end %func


%--------------------------------------------------------------------------
function cbf_select_table_(hTbl,e)
hFig = hTbl.Parent;
S_fig = hFig.UserData;
iRow = e.Indices(1,1);
iCol = e.Indices(1,2);
S_tbl = hTbl.UserData;
vhTable = S_fig.vhTable;
csMsg = {};
if iCol == 1
    set(vhTable, 'Enable', 'off');
    set(hTbl, 'Enable', 'on');
    vcType = hTbl.Tag;
    vcKey = hTbl.Data{iRow,1};    
    csMsg{end+1} = sprintf('Selected %s: %s', vcType, vcKey);
elseif iCol > 1 % clear selection
    set(vhTable, 'Enable', 'on');
    csMsg{end+1} = sprintf('All selected');
    vcType = 'All';
    vcKey = [];
%     @TODO: sort by param on all scores
end

% display selected info
S_fig.hText1.String = 'Select a parameter';
csParam1 = get_userdata_(hFig, 'csParam');
csInput1 = get_userdata_(hFig, 'csStudyset');
switch vcType
    case 'Param'
        S_fig.hText1.String = file2cellstr_(locate_('Param', vcKey));
        csParam1 = {vcKey};
    case 'Studyset', csInput1 = {vcKey};
    case 'Study', csInput1 = {vcKey};
    case 'Recording', csInput1 = {vcKey};
end
S_fig.hText2.String = csMsg;
set_userdata_(hFig, csParam1);
set_userdata_(hFig, csInput1);

% @ToDo
%  sort by parameters (# units > .8)
%  update parameter table using the selected values
end %func


%--------------------------------------------------------------------------
% param file: [hashcode].prm. provide full parameter set
% function cbf_select_param_(hTbl,e)
% hFig = hTbl.Parent;
% S_fig = hFig.UserData;
% iRow_ = e.Indices(1,1);
% iCol_ = e.Indices(1,2);
% S_tbl = hTbl.UserData;
% % iRow = get_(S_tbl, 'iRow');
% vhTable = S_fig.vhTable;
% vlEnable = cellfun(@(x)strcmpi(x, 'on'), get(vhTable, 'Enable'));
% assert(all(vlEnable), 'All must be enabled when this is called');
% 
% vcSelection = sprintf('Selected %s: %s', hTbl.Tag, hTbl.Data{iRow_,1});
% S_fig.hText.String = vcSelection;
% 
% %    @ToDo: update table values on the left
% end %func


%--------------------------------------------------------------------------
function hTable = create_table_(hFig, csFiles_raw, vcLabel, pos_table)
% [hTable, hText] = create_table_(hFig, csFiles_raw, vcLabel, pos_table)
% [hTable, hText] = create_table_([], csFiles_raw, vcLabel)

nFiles = numel(csFiles_raw);
csColumnName = irc('call', 'pad_cs', {{vcLabel, 'Frac>.8', 'Accu>8', 'Prec>8', 'Recl>8', '#Units'}, [100 5 5 5 5 5]});
nCols = numel(csColumnName);
vlEditable = false(1, nCols);
cData_rec = cell(nFiles, numel(csColumnName)); 
cData_rec(1:nFiles,1) = csFiles_raw;

% try to find the 
hTable = findobj('Type', 'uitable', 'Tag', vcLabel);
if isempty(hTable)
    switch 1
        case 2
            hTable = uitable(hFig, 'Data', cData_rec, ...
                'ColumnEditable', vlEditable, 'ColumnName', csColumnName, ...
                'Position', pos_table, 'Tag', vcLabel);
        case 1
            hTable = uitable(hFig, 'Data', cData_rec, ...
                'ColumnEditable', vlEditable, 'ColumnName', csColumnName, ...
                'Unit', 'Normalized','OuterPosition', pos_table, 'Tag', vcLabel);
    end
else
    hTable.Data = cData_rec;
end
end %func


%--------------------------------------------------------------------------
function update_table_(vcLabel, csFiles)
create_table_([], csFiles, vcLabel, []);
end %func


%--------------------------------------------------------------------------
function varargout = deal_(vr)
assert(nargout == numel(vr), 'deal_: number should match');
for iArg = 1:nargout
    varargout{iArg} = vr(iArg);
end %for
end %func


%--------------------------------------------------------------------------
function ylabel_(vhAx, vcLabel)
for iAx = 1:numel(vhAx)
    ylabel(vhAx(iAx), vcLabel);
end %for
end %func


%--------------------------------------------------------------------------
function xlabel_(vhAx, vcLabel)
for iAx = 1:numel(vhAx)
    xlabel(vhAx(iAx), vcLabel);
end %for
end %func


%--------------------------------------------------------------------------
function [csRecording, csStudy, csStudyset] = load_dataset_()
vcDir_in = locate_('input');
csFiles_raw = find_files_(vcDir_in, '/**/raw.mda');
% csFiles_raw = flipud(csFiles_raw); % put new files first
if vcDir_in(end) ~= '/', vcDir_in(end+1) = '/'; end
fh_list1 = @(cs)unique(cellfun_(@(x)fileparts(x), cs));
fh_list2 = @(cs)cellfun_(@(x)strrep(x, vcDir_in, ''), cs);
fh_list3 = @(cs)fh_list2(fh_list1(cs));
csRecording = fh_list3(csFiles_raw);
csStudy = fh_list3(csRecording);
csStudyset = fh_list3(csStudy);
end %func


%--------------------------------------------------------------------------
function [csFiles_raw, csDir_raw] = find_files_(vcDir_in, vcFile)
csFiles_raw = arrayfun(@(x)fullfile(x.folder, x.name), ...
    dir(fullfile(vcDir_in, vcFile)), 'UniformOutput', 0);
if nargout>=2
    csDir_raw = cellfun_(@(x)fileparts(x), csFiles_raw);
end
end %func


%--------------------------------------------------------------------------
% load param hash code: vcDir_out\param\paramhash.prm
function csParam = load_param_()
% vcDir_out = locate_('output');
csFiles_param = locate_('Param');
% csFiles_param = arrayfun(@(x)fullfile(x.folder, x.name), ...
%     dir(fullfile(vcDir_out, '/param/*_template.prm')), 'UniformOutput', 0);
[~, csParam, ~] = cellfun_(@fileparts, csFiles_param);
csParam = cellfun_(@(x)strrep(x, '_template', ''), csParam);
end %func


%--------------------------------------------------------------------------
% 10/8/17 JJJ: Created
% 3/20/18 JJJ: captures edit failure (when running "matlab -nodesktop")
function edit_(vcFile)
% vcFile0 = vcFile;
if ~exist_file_(vcFile)
    vcFile = ircpath_(vcFile, 1);
end
if ~exist_file_(vcFile)
    fprintf(2, 'File doesn not exist: %s\n', vcFile);
else
    fprintf('Editing %s\n', vcFile);
    try edit(vcFile); catch, end
end
end %func


%==========================================================================
% call irc.m


%--------------------------------------------------------------------------
function vcFile_full = ircpath_(vcFile, fConditional)
if nargin<2, fConditional=1; end
vcFile_full = irc('call', 'ircpath_', {vcFile, fConditional});
end %func


%--------------------------------------------------------------------------
function vl = exist_file_(vc)
vl = irc('call', 'exist_file_', {vc});
end %func


%--------------------------------------------------------------------------
function vl = file2cellstr_(vc)
vl = irc('call', 'file2cellstr_', {vc});
end %func


%==========================================================================
% from irc.m


%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
end %func


%--------------------------------------------------------------------------
function S_userdata = set_userdata_(varargin)
% Usage
% -----
% S_userdata = set_userdata_(h, val)
% S_userdata = set_userdata_(h, name1, val1, name2, val2, ...)
h = varargin{1};
S_userdata = get(h, 'UserData');
if nargin==2
    [val, vcName] = deal(varargin{2}, inputname(2));
    try
        S_userdata.(vcName) = val;
    catch
        disp(lasterr());
    end
elseif nargin>2
    
    for iArg_in = 2:2:nargin
        [vcName1, val1] = deal(varargin{iArg_in}, varargin{iArg_in+1});
        S_userdata.(vcName1) = val1;
    end %for
else
    return;
end
set(h, 'UserData', S_userdata);
end %func


%--------------------------------------------------------------------------
% 1/31/2019 JJJ: get the field(s) of a struct or index of an array or cell
function varargout = get_(varargin)
% same as struct_get_ function
% retrieve a field. if not exist then return empty
% [val1, val2] = get_(S, field1, field2, ...)
% [val] = get_(cell, index)

if nargin==0, varargout{1} = []; return; end
S = varargin{1};
if isempty(S), varargout{1} = []; return; end

if isstruct(S)
    for i=2:nargin
        vcField = varargin{i};
        try
            varargout{i-1} = S.(vcField);
        catch
            varargout{i-1} = [];
        end
    end
elseif iscell(S)
    try    
        varargout{1} = S{varargin{2:end}};
    catch
        varargout{1} = [];
    end
else
    try    
        varargout{1} = S(varargin{2:end});
    catch
        varargout{1} = [];
    end
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test (.m format)
function struct2file_(S, vcFile)
% Modify the parameter file using the variables in the P struct

csName = fieldnames(S);
% csValue = cellfun(@(vcField)S.(vcField), csName, 'UniformOutput',0);
csLines = cell(size(csName));
for iLine=1:numel(csName)
    vcName_ = csName{iLine}; %find field name with 
    val_ = S.(vcName_);
    if isstruct(val_), continue; end %do not write struct
    csLines{iLine} = sprintf('%s = %s;', vcName_, field2str_(val_));
end
cellstr2file_(vcFile, csLines);
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function cellstr2file_(vcFile, csLines, fVerbose)
% Write a cellstring to a text file
if nargin<3, fVerbose = 0; end
fid = fopen(vcFile, 'w');
for i=1:numel(csLines)
    fprintf(fid, '%s\n', csLines{i});
end
fclose(fid);
if fVerbose
    fprintf('Wrote to %s\n', vcFile);
end
end %func


%---------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function vcStr = field2str_(val, fDoubleQuote)
% convert a value to a strong
if nargin<2, fDoubleQuote = false; end

switch class(val)
    case {'int', 'int16', 'int32', 'uint16', 'uint32'}
        vcFormat = '%d';
    case {'double', 'single'}
        vcFormat = '%g';
        if numel(val)==1
            if mod(val(1),1)==0, vcFormat = '%d'; end
        end
    case 'char'
        if fDoubleQuote
            vcStr = sprintf('"%s"', val);
        else
            vcStr = sprintf('''%s''', val);
        end
        return;
    case 'cell'
        vcStr = '{';
        for i=1:numel(val)
            vcStr = [vcStr, field2str_(val{i})];
            if i<numel(val), vcStr = [vcStr, ', ']; end
        end
        vcStr = [vcStr, '}'];
        return;
    case 'logical'
        vcFormat = '%s';
        if val
            vcStr = '1';
        else
            vcStr = '0';
        end
    otherwise
        vcStr = '';
        fprintf(2, 'field2str_: unsupported format: %s\n', class(val));
        return;
end

if numel(val) == 1
    vcStr = sprintf(vcFormat, val);
else % Handle a matrix or array
    vcStr = '[';
    for iRow=1:size(val,1)
        for iCol=1:size(val,2)
            vcStr = [vcStr, field2str_(val(iRow, iCol))];
            if iCol < size(val,2), vcStr = [vcStr, ', ']; end
        end
        if iRow<size(val,1), vcStr = [vcStr, '; ']; end
    end
    vcStr = [vcStr, ']'];
end
end %func


%--------------------------------------------------------------------------
% 12/11/2018 JJJ: multiple close supported
function close_(varargin)
for iArg = 1:nargin
    vh1 = varargin{iArg};
    for ih1 = 1:numel(vh1)       
        try 
            close(vh1(ih1)); 
        catch
            ;
        end
    end
end
end %func


%--------------------------------------------------------------------------
function val = get_userdata_(h, vcName, fDelete)
% fDelete: delete the field after retrieving
val = [];
if nargin<3, fDelete = 0; end
try
    S_userdata = get(h, 'UserData');
    if ~isfield(S_userdata, vcName), return; end
    val = S_userdata.(vcName);
    if fDelete
        set(h, 'UserData', rmfield(S_userdata, vcName));
    end
catch
    val = [];
end
end %func