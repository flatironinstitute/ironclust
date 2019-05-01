%--------------------------------------------------------------------------
function irc_scoreboard()
irc('addpath');
cd_irc_();

% create figure
hFig = irc('call', 'create_figure_', {'Fig_scoreboard', [0 0 .5 1], 'IronClust Scoreboard', 0, 1}); 
add_menu_(hFig);

% load settings
hMsg = msgbox('Loading files');
vcFile_cache = locate_('cache');
if exist_file_(vcFile_cache)
    disp('Loading from cache');
    S_fig = load(vcFile_cache);
    hFig.UserData = S_fig;
    [csRecording, csStudy, csStudyset, csParam] = ...
        get_userdata_(hFig, 'csRecording', 'csStudy', 'csStudyset', 'csParam');
else
    [csRecording, csStudy, csStudyset] = load_input_();
    [csParam] = load_param_();    
    S_fig = set_userdata_(hFig, csRecording, csStudy, csStudyset, csParam);
end
try, close(hMsg); catch; end

% Add axes
fh_ax1 = @(x)axes('Parent', hFig, 'OuterPosition', x);
%vhAx = cellfun(@(x)fh_ax1(x), {[0 .1 .15 .3], [.15 .1 .15 .3], [.35 .1 .15 .3], [.5 .1 .15 .3], [.70 .1 .15 .3], [.85 .1 .15 .3]});
vhAx = cellfun(@(x)fh_ax1(x), {[0 .1 .33 .3], [0 0 0 0], [.33 .1 .33 .3], [0 0 0 0], [.66 .1 .33 .3], [0 0 0 0]});
[hAx_accu1, hAx_accu2, hAx_prec1, hAx_prec2, hAx_recl1, hAx_recl2] = deal_(vhAx);
xlabel_(vhAx(1:2:end), 'SNR (Vmin/Vrms)');
% xlabel_(vhAx(2:2:end), 'Param Rank');
ylabel_(vhAx(1), 'Accuracy');
ylabel_(vhAx(3), 'Precision');
ylabel_(vhAx(5), 'Recall');
arrayfun(@(x)grid(x, 'on'), vhAx);
arrayfun(@(x)ylim(x, [0 1]), vhAx);
S_fig = set_userdata_(hFig, vhAx);
set(hFig, 'UserData', S_fig);

% add texts
hText1 = uicontrol(hFig, 'Style','text','String','',...
       'Unit', 'Normalized', 'OuterPosition', [.5 .4 .5 .3], ...
       'HorizontalAlignment', 'left', 'BackgroundColor','w', 'Tag', 'text_param');
hText1.String = {'Select a parameter'}; % display meta
hText2  = uicontrol(hFig, 'Style','text','String','',...
       'Unit', 'Normalized', 'OuterPosition', [0 0 1 .1], ...
       'HorizontalAlignment', 'left', 'BackgroundColor','w', 'Tag', 'text_selection');   
hText2.String = {'All selected'}; % display meta
S_fig = set_userdata_(hFig, hText1, hText2);

% create table
hTable_studyset = create_table_(hFig, csStudyset, 'Studyset', [0 .8 .5 .2]);
hTable_study = create_table_(hFig, csStudy, 'Study', [0 .6 .5 .2]);
hTable_recording = create_table_(hFig, csRecording, 'Recording', [0 .4 .5 .2]);
hTable_param = create_table_(hFig, csParam, 'Param', [.5 .7 .5 .3]);

% set table cbf
vhTable = [hTable_studyset, hTable_study, hTable_recording, hTable_param];
S_fig = set_userdata_(hFig, vhTable);
arrayfun(@(x)set(x, 'CellSelectionCallback', @cbf_select_table_), vhTable);
arrayfun(@(x)set(x, 'ButtonDownFcn', @cbf_sort_table_), vhTable);
update_tables_();

[csParam1, csInput1] = deal(csParam, csStudyset);
set_userdata_(hFig, csParam1, csInput1);
end %func


%--------------------------------------------------------------------------
function cd_irc_()
cd(fileparts(mfilename('fullpath')));
end %func


%--------------------------------------------------------------------------
function run_batch_(h, e)
cd_irc_();
hFig = h.Parent.Parent;
csInput1 = get_userdata_(hFig, 'csInput1');
csParam1 = get_userdata_(hFig, 'csParam1');
if isempty(csParam1), csParam1 = get_userdata_(hFig, 'csParam'); end
if isempty(csInput1), csInput1 = get_userdata_(hFig, 'csStudyset'); end

hMsgbox = msgbox('Running batch');
csDir_out = sbatch_mda_(csInput1, csParam1);
assemble_results_(csDir_out)
close_(hMsgbox);
end %func


%--------------------------------------------------------------------------
function cbf_assemble_results_(h,e)
h = msgbox('Assembling results');
[~, csDir_out] = find_files_(locate_('output'), '**/raw_geom_score.mat');
assemble_results_(csDir_out);
close_(h);
end %func


%--------------------------------------------------------------------------
function assemble_results_(csDir_out)
% assemble_results_(csDir_out)
% assemble_results_()
hFig = get_fig_();
vcDir_out = locate_('output');
csParam = get_userdata_(hFig, 'csParam');
csRecording = get_userdata_(hFig, 'csRecording');
[cmSnr, cmPrecision, cmRecall, cmAccuracy] = deal(cell(numel(csRecording), numel(csParam)));
for iDir = 1:numel(csDir_out)
    vcFile_gt1 = fullfile(csDir_out{iDir}, 'raw_geom_score.mat');
    [vcParam1, vcRecording1] = strtok(strrep(csDir_out{iDir}, vcDir_out, ''), '/');
    iParam1 = find(strcmp(vcParam1, csParam));
    iRecording1 = find(strcmp(vcRecording1, csRecording));
    try
        S_score1 = load(vcFile_gt1);
        S_score_clu1 = S_score1.S_score_clu;
        [vrSnr1, vrPrecision1, vrRecall1, vrAccuracy1] = ...
            deal(S_score1.vrSnr_min_gt, 1-S_score_clu1.vrFp, 1-S_score_clu1.vrMiss, ...
                S_score_clu1.vrAccuracy);
        cmSnr{iRecording1,iParam1} = vrSnr1(:);
        cmPrecision{iRecording1,iParam1} = vrPrecision1(:);
        cmRecall{iRecording1,iParam1} = vrRecall1(:);
        cmAccuracy{iRecording1,iParam1} = vrAccuracy1(:);
    catch
        ; %no output generated
    end
end

S_fig = set_userdata_(hFig, cmSnr, cmPrecision, cmRecall, cmAccuracy);
save_cache_(S_fig);
update_tables_();
% close_(hMsgbox);
end %func


%--------------------------------------------------------------------------
function S_fig = save_cache_(S_fig)
if nargin<1
    [hFig, S_fig] = get_fig_();
end
irc('call', 'write_struct', {locate_('cache'), S_fig});
end %func


%--------------------------------------------------------------------------
function [hFig, S_fig] = get_fig_()
persistent hFig_
if ~ishandle(hFig_), hFig_ = []; end
if isempty(hFig_)
    hFig_ = findobj('Type', 'Figure', 'Tag', 'Fig_scoreboard');
end
hFig = hFig_;
S_fig = hFig.UserData;
end %func


%--------------------------------------------------------------------------
function update_tables_(vcType, vcKey)
% arrange the param from the best to worst
if nargin<1, vcType = ''; end
if nargin<2, vcKey = ''; end

[nCols, snr_thresh, accuracy_thresh, iCol_accuracy] = deal(6, 8, .8, 3);

hFig = get_fig_();
vhTable = get_userdata_(hFig, 'vhTable');
[cmSnr, cmPrecision, cmRecall, cmAccuracy, csParam, csParam1] = ...
    get_userdata_(hFig, 'cmSnr', 'cmPrecision', 'cmRecall', 'cmAccuracy', 'csParam', 'csParam1');
[nRecording, nParam] = size(cmSnr);

switch lower(vcType)
    case 'param' % parameter is generated, update table, reset
        csParam = vcKey;
        hTbl_param = vhTable(4);
        hTbl_param.Data = cell(numel(vcKey), nCols);
        hTbl_param.Data(:,1) = csParam;
        [cmSnr, cmPrecision, cmRecall, cmAccuracy] = deal([]);
        set_userdata_(hFig, cmSnr, cmPrecision, cmRecall, cmAccuracy, csParam);
        save_cache_();
        return;
        
    case 'input'
        csParam1 = vcKey;
        if isempty(cmSnr), return; end
        [~, viParam] = ismember(csParam1, csParam);
        fh_ = @(x)cell2mat_(x(:,viParam));
        [mrSnr1, mrPrecision1, mrRecall1, mrAccuracy1] = ...nC_max
            deal(fh_(cmSnr), fh_(cmPrecision), fh_(cmRecall), fh_(cmAccuracy));
        vhAx = get_userdata_(hFig, 'vhAx');
        arrayfun(@(x)hold(x,'off'), vhAx);
        plot(vhAx(1), mrSnr1(:), mrAccuracy1(:), '.'); xylabel_(vhAx(1), 'SNR', 'Accuracy');
        plot(vhAx(3), mrSnr1(:), mrPrecision1(:), '.'); xylabel_(vhAx(3), 'SNR', 'Precision');
        plot(vhAx(5), mrSnr1(:), mrRecall1(:), '.'); xylabel_(vhAx(5), 'SNR', 'Recall');
        arrayfun(@(x)ylim(x,[0,1]), vhAx);
        arrayfun(@(x)grid(x,'on'), vhAx);
end
if isempty(cmSnr), return; end

% param table
[hTbl_studyset, hTbl_study, hTbl_recording, hTbl_param] = deal_(vhTable);
csParam_tbl = hTbl_param.Data(:,1);
% hTbl_param.Data = cell(nParam, nCols);
for iParam_tbl = 1:numel(csParam_tbl)    
    iParam = strcmp(csParam_tbl{iParam_tbl}, csParam);
    fh_ = @(x)cell2mat_(x(:,iParam));
    [vrSnr1, vrPrec1, vrRecall1, vrAccuracy1] = ...
        deal(fh_(cmSnr), fh_(cmPrecision), fh_(cmRecall), fh_(cmAccuracy));   
    
    vl_ = vrSnr1 > snr_thresh;
    fh_ = @(x)nanmean(x(vl_));
    vr_ = [nanmean(vrAccuracy1>accuracy_thresh), fh_(vrAccuracy1), fh_(vrPrec1), fh_(vrRecall1), numel(vl_)];
%     hTbl_param.Data{iParam, 1} = csParam{iParam};
    hTbl_param.Data(iParam_tbl, 2:end) = arrayfun(@num2str, vr_, 'UniformOutput', 0);
end %for
% accuracy rank
% precision rank
% recall rank

% pick the best overall by Accuracy and determine iParam index
if numel(csParam1) == numel(csParam)
    viParam = 1:numel(csParam);
elseif numel(csParam1) >= 1
    [~, viParam] = ismember(csParam1, csParam);
elseif isempty(csParam1)
    vrAccuracy_param = cellfun(@str2num_, hTbl_param.Data(:,iCol_accuracy));
    [~, iParam_] = max(vrAccuracy_param);
    if isnan(iParam_), return; end
    vcParam_best = hTbl_param.Data{iParam_, 1};
    viParam = find(strcmp(vcParam_best, csParam));
end
% clear the recording table (use the best parameter)
% if parameter is selected, use that and pick the best? or average of them?
% if multiple parameters are selected show the averaged results
% if none selected show the best overall. 
% clear the table first
% hTbl_recording.Data(:, 2:end) = {''};
% hTbl_study.Data(:, 2:end) = {''};
% hTbl_studyset.Data(:, 2:end) = {''};

[csRecording, csStudy, csStudyset] = get_userdata_(hFig, 'csRecording', 'csStudy', 'csStudyset');
csRecording_tbl = hTbl_recording.Data(:,1);
% create studyset and study matrix
for iRecording_tbl = 1:nRecording
    try
        iRecording = find(strcmp(csRecording_tbl{iRecording_tbl}, csRecording));
        
        if isempty(iRecording), continue; end
        fh_ = @(x)nanmean(cell2mat_(x(iRecording,viParam)),2);
        [vrSnr1, vrPrec1, vrRecall1, vrAccuracy1] = ...
            deal(fh_(cmSnr), fh_(cmPrecision), fh_(cmRecall), fh_(cmAccuracy));       
        if isempty(vrSnr1), continue; end
        vl_ = vrSnr1 > snr_thresh;
        fh_ = @(x)nanmean(x(vl_));
        vr_ = [nanmean(vrAccuracy1>accuracy_thresh), fh_(vrAccuracy1), fh_(vrPrec1), fh_(vrRecall1), numel(vl_)];
        hTbl_recording.Data(iRecording_tbl, 2:end) = arrayfun(@num2str, vr_, 'UniformOutput', 0);
    catch
        ;
    end
end %for

csStudy_tbl = hTbl_study.Data(:,1);
for iStudy_tbl = 1:numel(csStudy_tbl)
    try
        vcStudy1 = csStudy_tbl{iStudy_tbl};
        viRecording1 = find(contains(csRecording, vcStudy1));
        if isempty(viRecording1), return; end
        fh_ = @(x)nanmean(cell2mat_(x(viRecording1,viParam)),2);
        [vrSnr1, vrPrec1, vrRecall1, vrAccuracy1] = ...
            deal(fh_(cmSnr), fh_(cmPrecision), fh_(cmRecall), fh_(cmAccuracy));       
        if isempty(vrSnr1), continue; end
        vl_ = vrSnr1 > snr_thresh;
        fh_ = @(x)nanmean(x(vl_));
        vr_ = [nanmean(vrAccuracy1>accuracy_thresh), fh_(vrAccuracy1), fh_(vrPrec1), fh_(vrRecall1), numel(vl_)];
        hTbl_study.Data(iStudy_tbl, 2:end) = arrayfun(@num2str, vr_, 'UniformOutput', 0);
    catch
        ;
    end
end %for

csStudyset_tbl = hTbl_studyset.Data(:,1);
for iStudyset_tbl = 1:numel(csStudyset_tbl)
    try
        vcStudyset1 = csStudyset_tbl{iStudyset_tbl};
        viRecording1 = find(contains(csRecording, vcStudyset1));        
        if isempty(viRecording1), return; end
        fh_ = @(x)nanmean(cell2mat_(x(viRecording1,viParam)),2);
        [vrSnr1, vrPrec1, vrRecall1, vrAccuracy1] = ...
            deal(fh_(cmSnr), fh_(cmPrecision), fh_(cmRecall), fh_(cmAccuracy));       
        if isempty(vrSnr1), continue; end
        vl_ = vrSnr1 > snr_thresh;
        fh_ = @(x)nanmean(x(vl_));
        vr_ = [nanmean(vrAccuracy1>accuracy_thresh), fh_(vrAccuracy1), fh_(vrPrec1), fh_(vrRecall1), numel(vl_)];
        hTbl_studyset.Data(iStudyset_tbl, 2:end) = arrayfun(@num2str, vr_, 'UniformOutput', 0);
    catch
        ;
    end
end %for

% accuracy rank
% precision rank
% recall rank

end %func


%--------------------------------------------------------------------------
function vr = cell2mat_(cvr)
% create a matrix that is #vectors x # cells
% remove empty
vi = find(cellfun(@(x)~isempty(x), cvr));
vr = cell2mat(cvr(vi));
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
function csDir_out = sbatch_mda_(csInput, csParam)
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
csDir_out = {};
for iData = 1:numel(csInput)    
    vcInput1 = csInput{iData};
    fprintf('Processing %s\n', vcInput1); t1=tic;
    csDir_rec1 = locate_('recordings', vcInput1);    
    [csDir_in1, csDir_out1, csFile_param1] = deal({});
    for iRec = 1:numel(csDir_rec1)
        vcDir_rec1 = csDir_rec1{iRec};
        vcRecording1 = strrep(vcDir_rec1, vcDir_in, '');
        for iParam = 1:nParams
            vcDir_out1 = fullfile(vcDir_out, csParam{iParam}, vcRecording1);
            [csDir_in1{end+1}, csDir_out1{end+1}, csFile_param1{end+1}] = ...
                deal(vcDir_rec1, vcDir_out1, fh_param(csParam{iParam}));
        end
    end
    irc('sbatch-mda', csDir_in1, csDir_out1, csFile_param1, fWait);
    csDir_out = cat(1, csDir_out, csDir_out1(:));
    fprintf('\tFinished %s, took %0.1fs\n', vcInput1, toc(t1));
end

end %func


%--------------------------------------------------------------------------
function vcPath = locate_(vcType, vcKey)
if nargin<2, vcKey = []; end
persistent vcDir_in vcDir_out
if isempty(vcDir_in) || isempty(vcDir_out)
    vcDir_in = read_cfg_('path_groundtruth');
    vcDir_out = fullfile(read_cfg_('path_validation'), irc('version'));
end

switch lower(vcType)
    case 'disbatch'
        vcPath0 = fullfile(read_cfg_('path_validation'), 'disbatch');
        vcPath = find_files_(vcPath0, '*_status.txt');
        if isempty(vcPath), return; end
        if iscell(vcPath), vcPath = vcPath{1}; end
        
    case {'recordings', 'rec'} % return a cell of directories containing raw.mda        
        [~, vcPath] = find_files_(locate_('input', vcKey), '/**/raw.mda');
        
    case 'cache'
        vcPath = fullfile(vcDir_out, 'scoreboard.mat');
    case 'settings'
        vcPath = ircpath_('settings.scoreboard');
    case 'param'
        if isempty(vcKey)
            vcPath = flipud(find_files_(vcDir_out, '/param/*_template.prm'));
        elseif ischar(vcKey)
            vcPath = fullfile(vcDir_out, 'param', [vcKey, '_template.prm']);
        elseif iscell(vcKey)
            vcPath = cellfun_(@(x)fullfile(vcDir_out, 'param', [x, '_template.prm']), vcKey);
        end
    case {'input', 'in'}
        if ~isempty(vcKey)
            if ischar(vcKey)
                vcPath = fullfile(vcDir_in, vcKey);
                if vcPath(end) ~= filesep(), vcPath(end+1) = filesep(); end
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
                if vcPath(end) ~= filesep(), vcPath(end+1) = filesep(); end
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
uimenu(mh_file,'Label', 'Delete cache', 'Callback', @(h,e)delete_(locate_('cache')));
uimenu(mh_file,'Label', 'Open disbatch status', 'Callback', @(h,e)edit_(locate_('disbatch')));

mh_file = uimenu(hFig,'Label','Run'); 
uimenu(mh_file,'Label', 'Generate param files (hashcode_template.prm)', 'Callback', @generate_param_);
uimenu(mh_file,'Label', 'Run batch', 'Callback', @run_batch_);
uimenu(mh_file, 'Label', 'Assemble batch result', 'Callback', @cbf_assemble_results_);

mh_file = uimenu(hFig,'Label','Help'); 
uimenu(mh_file,'Label', 'Show help', 'Callback', @menu_help_);

% mh_file = uimenu(hFig,'Label','Sort'); 
% uimenu(mh_file,'Label', 'Sort Param table', 'Callback', @cbf_sort_table_);
% uimenu(mh_file,'Label', 'Sort Input tables', 'Callback', @cbf_sort_table_);
end %func


%--------------------------------------------------------------------------
function menu_help_(h,e)
csMsg = {};

csMsg{end+1} = 'To sort the table, left click the cell to sort and right click again';
csMsg{end+1} = 'To apply sorting, select the data on the left column and select `Run` menu -> `Run batch`';

msgbox(csMsg);
end %func


%--------------------------------------------------------------------------
function cbf_sort_table_(hTable,e)
if nargin<2, e=[]; end

if isnumeric(e)
    iCol_sort = e;
else
    Indices = get_userdata_(hTable, 'Indices');
    if isempty(Indices), return; end
    if size(Indices,1)>1, return; end %no multiselection
    iCol_sort = Indices(2);
end
try
    hTable.Data = sortrows(hTable.Data,iCol_sort, 'descend');
catch
    ;
end
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

% update parameter list
update_tables_('Param', load_param_());
close_(hMsg);
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
Indices = e.Indices;
if isempty(Indices), return; end

hMsgbox = msgbox('Updating...');
viRow = e.Indices(:,1);
iCol = e.Indices(1,2);
set_userdata_(hTbl, Indices);
S_tbl = hTbl.UserData;
vhTable = S_fig.vhTable;
csMsg = {};
if iCol == 1
    set(vhTable, 'Enable', 'off');
    set(hTbl, 'Enable', 'on');
    vcType = hTbl.Tag;
    csKey = hTbl.Data(viRow,1);    
    csMsg{end+1} = sprintf('Selected %s:', vcType);
    csMsg = cat(1, csMsg, csKey);
    if numel(csKey)==1
        vcKey = csKey{1};
    else
        vcKey = csKey;
    end
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
        if ischar(vcKey)
            S_fig.hText1.String = file2cellstr_(locate_('Param', vcKey));
        elseif iscell(vcKey)
            cS_ = cellfun_(@file2struct_, locate_('Param', vcKey));
            S_fig.hText1.String = struct2str_(struct_intersect_(cS_{:}));
        end
        csParam1 = csKey;
        update_tables_('Input', csParam1);
    case {'Studyset', 'Study', 'Recording'}, csInput1 = ifeq_(ischar(vcKey), {vcKey}, vcKey);
end
S_fig.hText2.String = csMsg;
set_userdata_(hFig, csParam1, csInput1);
close_(hMsgbox);
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
hTable = create_table_([], csFiles, vcLabel, []);
hFig = hTable.Parent;
switch lower(vcLabel)
    case 'param'
        csParam = csFiles;
        set_userdata_(hFig, csParam);
end
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
function [csRecording, csStudy, csStudyset] = load_input_()
vcDir_in = locate_('input');
csFiles_raw = find_files_(vcDir_in, '/**/raw.mda');
fh_list1 = @(cs)unique(cellfun_(@(x)fileparts(x), cs));
fh_list2 = @(cs)cellfun_(@(x)strrep(x, vcDir_in, ''), cs);
fh_list3 = @(cs)flipud(fh_list2(fh_list1(cs)));
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

function xylabel_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function out1 = title_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = str2num_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = read_cfg_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ircpath_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2cellstr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct2str_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ifeq_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end


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
    for iArg_in = 2:nargin
        [vcName1, val1] = deal(inputname(iArg_in), varargin{iArg_in});
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
function varargout = get_userdata_(varargin)
% varargout = get_userdata_(h, name1, name2, ...)

val = [];
h = varargin{1};
S_userdata = get(h, 'UserData');
for iArg_out = 1:nargout
    vcName_ = varargin{iArg_out + 1};
    if isfield(S_userdata, vcName_)
        varargout{iArg_out} = S_userdata.(vcName_);
    else
        varargout{iArg_out} = [];
    end
end %for
end %func


%--------------------------------------------------------------------------
% 4/23/2019 JJJ: delete either cell of files or multiple arguments
function delete_(varargin)
for iArg = 1:nargin
    csFiles = varargin{iArg};
    if ischar(csFiles), csFiles = {csFiles}; end
    for i=1:numel(csFiles)
        try
            if iscell(csFiles)
                delete(csFiles{i});
            else
                delete(csFiles(i));
            end
        catch
    %         disperr_();
        end
    end
end
end %func


%--------------------------------------------------------------------------
function S_out = struct_intersect_(varargin)

S1 = varargin{1};
for iArg = 2:nargin
    S2 = varargin{iArg};
    csName1 = fieldnames(S1);
    csName2 = fieldnames(S2);
    csName12 = intersect(csName1, csName2);
    S_out = struct();
    if isempty(csName12), return; end

    for iField = 1:numel(csName12)
        vcName_ = csName12{iField};
        val1_ = S1.(vcName_);
        val2_ = S2.(vcName_);
        if ~strcmp(class(val1_), class(val2_))
            fDiff_ = true;
        elseif isempty(val1_) && isempty(val2_)
            fDiff_ = false;
        elseif ischar(val1_)
            fDiff_ = ~strcmpi(val1_, val2_);
        elseif numel(val1_) ~= numel(val2_)
            fDiff_ = true;
        else
            fDiff_ = ~all(val1_ == val2_);
        end
        if ~fDiff_
            S_out.(vcName_) = val1_;
        end
    end %for
    S1 = S_out;
end %arg
end %func