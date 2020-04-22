%--------------------------------------------------------------------------
% 2020/1/9 

function irc2_manual(P)
% global fDebug_ui

S0 = load0_(P.vcFile_prm);

% S_manual
S0.S_manual = load_(strrep(P.vcFile_prm, '.prm', '_manual_irc.mat'));
if ~isempty(S0.S_manual)
    if strcmpi(questdlg_('Load last saved?', 'Confirmation'), 'no')
        S0.S_manual=[];
    end
end
S0.iShank = 1;
if isempty(S0.S_manual), S0.S_manual = create_S_manual_(S0); end

clear mouse_figure;
clear get_fig_cache_ get_tag_ %clear persistent figure handles
close_(get_fig_('FigTrial')); %close previous FigTrial figure
close_(get_fig_('FigTrial_b')); %close previous FigTrial figure
S0 = struct_merge_(S0, ...
    struct('iCluCopy', 1, 'iCluPaste', [], 'hCopy', [], 'hPaste', [], 'nSites', numel(P.viSite2Chan)));
S0 = struct_merge_(S0, figures_manual_(P)); %create figures for manual interface

hMsg = msgbox_('Plotting... (this closes automatically)'); t1=tic;
plot_FigRD_(S0); % ask user before doing so
plot_FigWav_(S0); % hFigWav %do this after for ordering
plot_FigWavCor_(S0);  
S0 = button_CluWav_simulate_(1, [], S0); %select first clu
auto_scale_proj_time_(S0);
keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
set(0, 'UserData', S0);

% logging
% S_log = load_(strrep(P.vcFile_prm, '.prm', '_log.mat'), [], 0);
% if ~isempty(S_log), S0.cS_log = {S_log}; end
% save_log_('start', S0); %crash proof log

% Finish up
close_(hMsg);
fprintf('UI creation took %0.1fs\n', toc(t1));
figure_wait_(0, S0.FigWav);
end %func


%--------------------------------------------------------------------------
function plot_FigWav_(S0)

[P, S_manual] = get_(S0, 'P', 'S_manual');
[hFig, S_fig] = get_fig_cache_('FigWav'); 

% Show number of spikes per clusters
P.LineWidth = 1; %plot a thicker line
P.viSite_clu = S_manual.viSite_clu;
nSites = numel(P.viSite2Chan);
if isempty(S_fig) % initialize
    S_fig.maxAmp = P.maxAmp;
    S_fig.hAx = axes(hFig);
    set(S_fig.hAx, 'Position', [.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual'); 
    grid(S_fig.hAx, 'on');    
    S_fig.vcTitle = '%0.1f uV; [H]elp; (Sft)[Up/Down]; (Sft)[Left/Right]; [M]erge; [S]plit; [D]elete; [A]:Resample; [P]STH; [Z]oom; [Space]:Match';
    xylabel_(S_fig.hAx, 'Cluster #', 'Site #', sprintf(S_fig.vcTitle, S_fig.maxAmp));

    set(hFig, 'KeyPressFcn', @keyPressFcn_FigWav_, 'CloseRequestFcn', @exit_manual_, 'BusyAction', 'cancel');
    axis(S_fig.hAx, [0, S_manual.nClu + 1, 0, nSites + 1]);
    add_menu_(hFig, P);      
%     vrPos_ = get(hFig, 'OuterPosition');
    mouse_figure(hFig, S_fig.hAx, @button_CluWav_);
    S_fig = plot_spkwav_(S_fig, S0); %plot spikes
    S_fig = plot_tnWav_clu_(S_fig, P, S0); %do this after plotSpk_
    S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hSpkAll, S_fig);
%     set(hFig, 'OuterPosition', vrPos_);
else
    S_fig = plot_spkwav_(S_fig, S0); %plot spikes
    try delete(S_fig.vhPlot); catch; end %delete old text
    S_fig = rmfield_(S_fig, 'vhPlot');
    S_fig = plot_tnWav_clu_(S_fig, P, S0); %do this after plotSpk_
    xylabel_(S_fig.hAx, 'Cluster #', 'Site #', sprintf(S_fig.vcTitle, S_fig.maxAmp));
end

% create text
fText = get_set_(S_fig, 'fText', get_set_(P, 'Text', 1));
S_fig = figWav_clu_count_(S_fig, S_manual, fText);
S_fig.csHelp = { ...            
    '[Left-click] Cluter select/unselect (point at blank)', ...
    '[Right-click] Second cluster select (point at blank)', ...
    '[Pan] hold wheel and drag', ...
    '[Zoom] mouse wheel', ...
    '[X + wheel] x-zoom select', ...
    '[Y + wheel] y-zoom select', ...
    '[SPACE] clear zoom', ...
    '[(shift) UP]: increase amplitude scale', ...
    '[(shift) DOWN]: decrease amplitude scale', ...
    '[HOME] Select the first unit', ...
    '[END] Select the last unit', ...
    '------------------', ...
    '[H] Help', ...       
    '[S] Split auto', ...
    '[W] Spike waveforms (toggle)', ...                                  
    '[M] merge cluster', ...
    '[D] delete cluster', ...
    '[A] Resample spikes', ...
    '[Z] zoom selected cluster', ...
    '[R] reset view', ...
    '[1] reset window positions', ...
    '------------------', ...
    '[U] update all', ...  
    '[C] correlation plot', ...            
    '[T] show amp drift vs time', ...            
    '[J] projection view', ...            
    '[V] ISI return map', ...
    '[I] ISI histogram', ...
    '[E] Intensity map', ...
    '[P] PSTH display', ...
    '[O] Overlap average waveforms across sites', ...
    }; 
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function S0 = set0_(varargin)
% Set(0, 'UserData')

S0 = get(0, 'UserData'); 
% set(0, 'UserData', []); %prevent memory copy operation
for i=1:nargin
    try
        S0.(inputname(i)) = varargin{i};
    catch
        disperr_();
    end
end
set(0, 'UserData', S0);
end %func


%--------------------------------------------------------------------------
function S_fig = figWav_clu_count_(S_fig, S_manual, fText)
if nargin==0, [hFig, S_fig] = get_fig_cache_('FigWav'); end
if nargin<3, fText = 1; end

if fText
    csText_clu = arrayfun(@(i)sprintf('%d(%d)', i, S_manual.vnSpk_clu(i)), 1:S_manual.nClu, 'UniformOutput', 0);
else
    csText_clu = arrayfun(@(i)sprintf('%d', i), 1:S_manual.nClu, 'UniformOutput', 0);
end
set(S_fig.hAx, 'Xtick', 1:S_manual.nClu, 'XTickLabel', csText_clu, 'FontSize', 8);
try
    if fText
        xtickangle(S_fig.hAx, -20); 
    else
        xtickangle(S_fig.hAx, 0); 
    end
catch; 
end

S_fig.fText = fText;
if nargout==0
    hFig = get_fig_cache_('FigWav');
    set(hFig, 'UserData', S_fig);
end
end %func


%--------------------------------------------------------------------------
function cvhHide_mouse = mouse_hide_(hFig, hObj_hide, S_fig)
% hide during mouse pan to speed up     
if nargin<3, S_fig = get(hFig, 'UserData'); end
% if nargin<3, S0 = get(0, 'UserData'); end
if nargin == 0 %clear field
%     try S_fig = rmfield(S_fig, 'vhFig_mouse'); catch; end
    try S_fig = rmfield(S_fig, 'cvhHide_mouse'); catch; end
else
    if ~isfield(S_fig, 'vhFig_mouse') && ~isfield(S_fig, 'cvhHide_mouse')
%         S_fig.vhFig_mouse = hFig;
        S_fig.cvhHide_mouse = {hObj_hide};    
    else
%         S_fig.vhFig_mouse(end+1) = hFig;
        S_fig.cvhHide_mouse{end+1} = hObj_hide;
    end
end
cvhHide_mouse = S_fig.cvhHide_mouse;
if nargout==0, set(hFig, 'UserData', S_fig); end
end %func


%--------------------------------------------------------------------------
function S_fig = plot_spkwav_(S_fig, S0)
% fPlot_raw = 0;
if nargin<2, S0 = []; end
if isempty(S0), S0 = get(0, 'UserData'); end
[P, viSite_spk, S_manual] = deal(S0.P, S0.viSite_spk, S0.S_manual);

[cvrX, cvrY, cviSite] = deal(cell(S_manual.nClu, 1));
vnSpk = zeros(S_manual.nClu, 1);
viSites_clu_ = @(x)P.miSites(:, S_manual.viSite_clu(x));
if isfield(S_fig, 'maxAmp')
    maxAmp = S_fig.maxAmp;
else
    maxAmp = P.maxAmp;
end
nSpk_show = get_set_(P, 'nSpk_show', 30);
for iClu = 1:S_manual.nClu        
    try           
        trFet1 = S0.S_manual.ctrFet_sub_clu{iClu};
        viSpk_show = randomSelect_(1:size(trFet1,3), nSpk_show);
        trWav1 = pc2wav_(S0.mrPv_global, trFet1(:,:,viSpk_show));
        viSite_show = viSites_clu_(iClu);
        [cvrY{iClu}, cvrX{iClu}] = tr2plot_(trWav1, iClu, viSite_show, maxAmp, P);
        cviSite{iClu} = viSite_show;
        vnSpk(iClu) = size(trWav1, 3); %subsample 
    catch
        disperr_();
    end
end
S = makeStruct_(cvrY, cviSite, vnSpk);
try
    set(S_fig.hSpkAll, 'XData', cell2mat_(cvrX), 'YData', cell2mat_(cvrY), 'UserData', S);
catch
    S_fig.hSpkAll = plot_(S_fig.hAx, cell2mat_(cvrX), cell2mat_(cvrY), 'Color', [.5 .5 .5], 'LineWidth', .5); %, P.LineStyle); 
    set(S_fig.hSpkAll, 'UserData', S);
end
end %func


%--------------------------------------------------------------------------
function S = rmfield_(S, varargin)
% varargin: list of fields to remove
for i=1:numel(varargin)
    if isfield(S, varargin{i})
        S = rmfield(S, varargin{i});
    end
end
end %func


%--------------------------------------------------------------------------
function auto_split_(fMulti, S0)
% Auto-split feature that calls Hidehiko Inagaki's code
% 20160426
if nargin<1, fMulti = 0; end
if nargin<2, S0 = []; end
if isempty(S0), S0 = get(0, 'UserData'); end
[P, S_manual] = deal(S0.P, S0.S_manual);

if ~isempty(S0.iCluPaste), msgbox_('Select one cluster', 1); return; end
if S_manual.vnSpk_clu(S0.iCluCopy)<3
    msgbox_('At least three spikes required for splitting', 1); return; 
end
    
hMsg = msgbox_('Splitting... (this closes automatically)');
iClu1 = S0.iCluCopy;
iSite1 = S_manual.viSite_clu(iClu1);
if fMulti
    viSites1 = P.miSites(1:P.nSites_fet, iSite1);
else
    viSites1 = iSite1;
end
trSpkWav1 = tnWav2uV_(tnWav_spk_sites_(S_manual.cviSpk_clu{iClu1}, viSites1, S0), P, 0);
[vlSpkIn, mrFet_split, vhAx, hFigTemp] = auto_split_wav_(trSpkWav1, [], 2, viSites1);
[hPoly, hFig_wav] = deal([]);
try 
    drawnow_(); 
    close(hMsg); 
catch
    ;
end
while 1
    vcAns = questdlg_('Split?', 'confirmation', 'Yes', 'No', 'Manual', 'Yes');
    close_(hFig_wav, hPoly); 
    switch lower(vcAns)        
        case 'yes', close_(hFigTemp); break;
        case {'no', ''}, close(hFigTemp); return;            
        case 'manual'
            %vcAns = questdlg_('Select projection', '', 'PC1 vs PC2', 'PC3 vs PC2', 'PC1 vs PC3', 'PC1 vs PC2');
            csAns_opt = {'PC1 vs PC2', 'PC3 vs PC2', 'PC1 vs PC3', 'Waveform'};
            [iAns, fSelected] = listdlg('PromptString', 'Select a projection', ...
                'SelectionMode', 'single', 'ListString', csAns_opt);
            if ~fSelected, close_(hFigTemp); return; end
            vcAns = csAns_opt{iAns};                
            switch vcAns
                case 'PC1 vs PC2', [hAx_, iAx1, iAx2] = deal(vhAx(1), 1, 2);
                case 'PC3 vs PC2', [hAx_, iAx1, iAx2] = deal(vhAx(2), 3, 2);
                case 'PC1 vs PC3', [hAx_, iAx1, iAx2] = deal(vhAx(3), 1, 3);
                case 'Waveform', [vlSpkIn, hFig_wav] = wave_split_manual_(trSpkWav1, viSites1, P);
                otherwise, close_(hFigTemp); return; 
            end              
            if ~strcmpi(vcAns, 'Waveform')                
                axes(hAx_); 
                cla(hAx_);
                [vrX1, vrY1] = deal(mrFet_split(:,iAx1), mrFet_split(:,iAx2));
                plot_(hAx_, vrX1, vrY1, 'k.');                
                hPoly = impoly_();
                if isempty(hPoly), close_(hFigTemp); return; end 
                mrPolyPos = getPosition(hPoly);              
                vlSpkIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
                plot_(hAx_, vrX1(vlSpkIn), vrY1(vlSpkIn), 'b.', vrX1(~vlSpkIn), vrY1(~vlSpkIn), 'r.');
            end            
    end %switch
end
split_clu_(iClu1, vlSpkIn);
end %func


%--------------------------------------------------------------------------
function S0 = keyPressFcn_FigWav_(hObject, event, S0) %amp dist
global fDebug_ui
if isempty(fDebug_ui), fDebug_ui = false; end

if nargin<3, S0 = get(0, 'UserData'); end
[P, S_manual] = get_(S0, 'P', 'S_manual');
P.LineStyle=[];
nSites = numel(P.viSite2Chan);
hFig = hObject;
S_fig = get(hFig, 'UserData');

switch lower(event.Key)
    case {'uparrow', 'downarrow'}
        rescale_FigWav_(event, S0, P);
        clu_info_(S0); %update figpos

    case {'leftarrow', 'rightarrow', 'home', 'end'}
        % switch the current clu
        if strcmpi(event.Key, 'home')
            S0.iCluCopy = 1;
        elseif strcmpi(event.Key, 'end')
            S0.iCluCopy = S_manual.nClu;
        elseif ~key_modifier_(event, 'shift')
            if strcmpi(event.Key, 'leftarrow')
                if S0.iCluCopy == 1, return; end
                S0.iCluCopy = S0.iCluCopy - 1;
            else
                if S0.iCluCopy == S_manual.nClu, return; end
                S0.iCluCopy = S0.iCluCopy + 1;
            end
        else
            if isempty(S0.iCluPaste)
                S0.iCluPaste = S0.iCluCopy;
            end
            if strcmpi(event.Key, 'leftarrow')
                if S0.iCluPaste == 1, return; end
                S0.iCluPaste = S0.iCluPaste - 1;
            else
                if S0.iCluPaste == S_manual.nClu, return; end
                S0.iCluPaste = S0.iCluPaste + 1;
            end
        end
        S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0); %select first clu
        if strcmpi(event.Key, 'home') || strcmpi(event.Key, 'end') %'z' to recenter
            S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); 
        end        
    case 'm', S0 = ui_merge_(S0); % merge clusters                
    case 'space'
        ui_zoom_(S0, hFig);
        % auto-select nearest cluster for black
        mrWavCor = S_manual.mrWavCor;
        mrWavCor(S0.iCluCopy,S0.iCluCopy) = -inf;
        [~,S0.iCluPaste] = max(mrWavCor(:,S0.iCluCopy));
        set(0, 'UserData', S0);
        button_CluWav_simulate_([], S0.iCluPaste);        
    case 's', auto_split_(1, S0);        
    case 'r' %reset view
        hFig_wait = figure_wait_(1);
        axis_([0, S_manual.nClu + 1, 0, numel(P.viSite2Chan) + 1]);
        figure_wait_(0, hFig_wait);        
    case {'d', 'backspace', 'delete'}, S0 = ui_delete_(S0);        
    case 'z', ui_zoom_(S0, hFig);  %zoom
    case 'c', plot_FigCorr_(S0);        
    case 'v', plot_FigIsi_(S0);        
    case 'a', update_spikes_(S0); clu_info_(S0);        
    case 'f', clu_info_(S0);               
    case 'h', msgbox_(S_fig.csHelp, 1);
    case 'w', toggleVisible_(S_fig.hSpkAll); %toggle spike waveforms 
    case 't', plot_FigTime_(S0); % time view        
    case 'j', plot_FigProj_(S0); %projection view        
    case 'n'
        fText = get_set_(S_fig, 'fText', get_set_(P, 'fText', 1));
        figWav_clu_count_(S_fig, S_manual, ~fText);          
    case 'i', plot_FigHist_(S0); %ISI histogram               
    case 'e', plot_FigMap_(S0);        
    case 'u', update_FigCor_(S0);        
    case 'p', plot_psth_trial_(S0, 1);
    case '1', reset_position_();
    otherwise, figure_wait_(0); %stop waiting
end
figure_(hObject); %change the focus back to the current object
end %func


%--------------------------------------------------------------------------
function ui_zoom_(S0, hFig)
iClu = S0.iCluCopy;
S_manual = S0.S_manual;
P = S0.P;
iSiteClu = S_manual.viSite_clu(S0.iCluCopy);
nSites_spk = size(P.miSites,1);
nSites = size(P.miSites,2);
set_axis_(hFig, iClu+[-1,1]*8, iSiteClu+[-1,1]*nSites_spk, [0 S_manual.nClu+1], [0 nSites+1]);
end %func


%--------------------------------------------------------------------------
function set_axis_(hFig, xlim1, ylim1, xlim0, ylim0)
% set the window within the box limit
if nargin <= 3
    % square case
    xlim0 = ylim1;
    ylim1 = xlim1;
    ylim0 = xlim0;
end

hFig_prev = gcf;
figure(hFig); 
dx = diff(xlim1);
dy = diff(ylim1);

if xlim1(1)<xlim0(1), xlim1 = xlim0(1) + [0, dx]; end
if xlim1(2)>xlim0(2), xlim1 = xlim0(2) + [-dx, 0]; end
if ylim1(1)<ylim0(1), ylim1 = ylim0(1) + [0, dy]; end
if ylim1(2)>ylim0(2), ylim1 = ylim0(2) + [-dy, 0]; end

xlim1(1) = max(xlim1(1), xlim0(1));
ylim1(1) = max(ylim1(1), ylim0(1));
xlim1(2) = min(xlim1(2), xlim0(2));
ylim1(2) = min(ylim1(2), ylim0(2));

axis_([xlim1, ylim1]);

figure(hFig_prev);
end %func


%--------------------------------------------------------------------------
function add_menu_(hFig, P)
drawnow_();
posvec = get(hFig, 'OuterPosition');

set(hFig, 'MenuBar','None');
mh_file = uimenu_(hFig,'Label','File'); 
uimenu_(mh_file,'Label', 'Save', 'Callback', @save_manual_);
uimenu_(mh_file,'Label', 'Save figures as .fig', 'Callback', @(h,e)save_figures_('.fig'));
uimenu_(mh_file,'Label', 'Save figures as .png', 'Callback', @(h,e)save_figures_('.png'));

uimenu_(mh_file,'Label', 'Describe', 'Callback', @(h,e)msgbox_(describe_()), 'Separator', 'on');
uimenu_(mh_file,'Label', 'Edit prm file', 'Callback', @edit_prm_);
uimenu_(mh_file,'Label', 'Reload prm file', 'Callback', @reload_prm_);
uimenu_(mh_file,'Label', 'Copy prm file to clipboard', 'Callback', @clipboard_prm_);
uimenu_(mh_file,'Label', 'Open prm folder', 'Callback', @open_prm_folder_);

uimenu_(mh_file,'Label', 'Export units to csv', 'Callback', @export_csv_, 'Separator', 'on');
uimenu_(mh_file,'Label', 'Export unit qualities to csv', 'Callback', @(h,e)export_quality_);
uimenu_(mh_file,'Label', 'Export all mean unit waveforms', 'Callback', @export_tmrWav_clu_);
uimenu_(mh_file,'Label', 'Export selected mean unit waveforms', 'Callback', @(h,e)export_mrWav_clu_);
uimenu_(mh_file,'Label', 'Export all waveforms from the selected unit', 'Callback', @(h,e)export_tnWav_spk_);
uimenu_(mh_file,'Label', 'Export firing rate for all units', 'Callback', @(h,e)export_rate_);
uimenu_(mh_file,'Label', 'Export to NeuroSuite Klusters', 'Callback', @(h,e)export_klusters_, 'Separator', 'on');
uimenu_(mh_file,'Label', 'Export to Phy', 'Callback', @(h,e)export_phy_, 'Separator', 'on');
uimenu_(mh_file,'Label', 'Exit', 'Callback', @exit_manual_, 'Separator', 'on', 'Accelerator', 'Q', 'Separator', 'on');

% shank menu
mh_shank = uimenu_(hFig,'Label','Shank', 'Tag', 'mh_shank'); 
menu_shank_ = @(i)uimenu_(mh_shank, 'Label', sprintf('Shank %d',i), 'Callback', @(h,e)select_shank_(h,e));
viShank_uniq = unique(get_(P, 'viShank_site', []));
if isempty(viShank_uniq), viShank_uniq = 1; end
cell_menu_shank = arrayfun_(@(x)menu_shank_(x), viShank_uniq);
set_userdata_(mh_shank, cell_menu_shank);
set(cell_menu_shank{1}, 'Checked', 'on');

mh_edit = uimenu_(hFig,'Label','Edit'); 
uimenu_(mh_edit,'Label', '[M]erge', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'm'));
uimenu_(mh_edit,'Label', 'Merge auto', 'Callback', @(h,e)merge_auto_());
uimenu_(mh_edit,'Label', '[D]elete', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'd'), 'Separator', 'on');
uimenu_(mh_edit,'Label', 'Delete auto', 'Callback', @(h,e)delete_auto_());
uimenu_(mh_edit,'Label', '[S]plit', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 's'), 'Separator', 'on');
uimenu_(mh_edit,'Label', 'Auto split max-chan', 'Callback', @(h,e)auto_split_(0));
uimenu_(mh_edit,'Label', 'Auto split multi-chan', 'Callback', @(h,e)auto_split_(1));
uimenu_(mh_edit,'Label', 'Annotate', 'Callback', @(h,e)unit_annotate_());

mh_view = uimenu_(hFig,'Label','View'); 
uimenu_(mh_view,'Label', 'Show traces', 'Callback', @(h,e)traces_());
uimenu_(mh_view,'Label', 'View all [R]', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'r'));
uimenu_(mh_view,'Label', '[Z]oom selected', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'z'));
uimenu_(mh_view,'Label', '[W]aveform (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'w'));
uimenu_(mh_view,'Label', '[N]umbers (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'n'));
% uimenu_(mh_view,'Label', 'Show raw waveform', 'Callback', @(h,e)raw_waveform_(h), ...
%     'Checked', ifeq_(get_(P, 'fWav_raw_show'), 'on', 'off'));
%uimenu_(mh_view,'Label', 'Threshold by sites', 'Callback', @(h,e)keyPressFcn_thresh_(hFig, 'n'));
% uimenu_(mh_view,'Label', '.prm file', 'Callback', @edit_prm_);
uimenu_(mh_view,'Label', 'Show averaged waveforms on all channels','Callback', @(h,e)ui_show_all_chan_(1,h));
uimenu_(mh_view,'Label', 'Show global drift','Callback', @(h,e)plot_drift_());
uimenu_(mh_view,'Label', 'Show drift view','Callback', @(h,e)ui_show_drift_view_(1,h));
uimenu_(mh_view,'Label', 'Reset window positions[1]', 'Callback', @reset_position_);

mh_proj = uimenu_(hFig,'Label','Projection'); 
uimenu_(mh_proj, 'Label', 'vpp', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.vcFet_show, {'vpp', 'vmin'}));
uimenu_(mh_proj, 'Label', 'pca', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.vcFet_show, {'pca'}));
uimenu_(mh_proj, 'Label', 'ppca', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.vcFet_show, {'ppca', 'private pca'}));
% uimenu_(mh_proj, 'Label', 'cov', 'Callback', @(h,e)proj_view_(h), ...
%     'Checked', if_on_off_(P.vcFet_show, {'cov', 'spacetime'}));

mh_trials = uimenu_(hFig,'Label','Trials', 'Tag', 'mh_trials');
set_userdata_(mh_trials, P);
update_menu_trials_(mh_trials);

mh_info = uimenu_(hFig,'Label','','Tag', 'mh_info'); 
uimenu_(mh_info, 'Label', 'Annotate unit', 'Callback', @unit_annotate_);
uimenu_(mh_info, 'Label', 'single', 'Callback', @(h,e)unit_annotate_(h,e,'single'));
uimenu_(mh_info, 'Label', 'multi', 'Callback', @(h,e)unit_annotate_(h,e,'multi'));
uimenu_(mh_info, 'Label', 'noise', 'Callback', @(h,e)unit_annotate_(h,e,'noise'));
uimenu_(mh_info, 'Label', 'clear annotation', 'Callback', @(h,e)unit_annotate_(h,e,''));
uimenu_(mh_info, 'Label', 'equal to', 'Callback', @(h,e)unit_annotate_(h,e,'=%d'));

% mh_history = uimenu_(hFig, 'Label', 'History', 'Tag', 'mh_history'); 

mh_help = uimenu_(hFig,'Label','Help'); 
uimenu_(mh_help, 'Label', '[H]elp', 'Callback', @help_FigWav_);
uimenu_(mh_help, 'Label', 'Wiki on GitHub', 'Callback', @(h,e)wiki_());
uimenu_(mh_help, 'Label', 'About', 'Callback', @(h,e)msgbox_(about_()));
uimenu_(mh_help, 'Label', 'Post an issue on GitHub', 'Callback', @(h,e)issue_('search'));
uimenu_(mh_help, 'Label', 'Search issues on GitHub', 'Callback', @(h,e)issue_('post'));

drawnow_();
set(hFig, 'OuterPosition', posvec);
end %func


%--------------------------------------------------------------------------
function exit_manual_(src, event)
try    
    if ~ishandle(src), return; end
    if ~isvalid(src), return; end
    S0 = get(0, 'UserData'); 
    fExit = save_manual_();
    if ~fExit, return; end 
    % These figures get closed when the main window gets closed.    
    csFig_close = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', ...
        'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist', 'FigClust', ...
        'FigAux', 'FigDrift', 'FigDriftAll', 'FigImport_trial'};
    if ~isfield(S0, 'csFig')
        S0.csFig = csFig_close;
    else
        S0.csFig = union(S0.csFig, csFig_close);
    end
    delete_multi_(get_fig_all_(S0.csFig), src);
    close_(get_fig_('FigTrial'), get_fig_('FigTrial_b'), get_fig_('FigAux'));    
catch
    disperr_();
    close(src);
end
set(0, 'UserData', []); % clear previous
end %func


%--------------------------------------------------------------------------
function plot_FigRD_(S0)
% P = funcDefStr_(P, ...
%     'delta1_cut', .5, 'rho_cut', -2.5, 'fDetrend', 0, 'fAskUser', 1, ...
%     'min_count', 50, 'y_max', 2, 'rhoNoise', 0, 'minGamma', -1, 'fNormDelta', 0, 'fExclMaxRho', 0, 'fLabelClu', 1);
% P.y_max = 1;
% P.fAskUser = 1;

[hFig, S_fig] = get_fig_cache_('FigRD');
hAx = get_(S_fig, 'hAx'); 
if isempty(hAx)
    hAx = axes(hFig);
    S_fig.hAx = hAx;
else    
    cla(hAx);     
end
hold(hAx, 'on'); 

[rho, delta, icl] = get_(S0.S_clu, 'rho', 'delta', 'icl');
delta_cut = get_set_(S0.P, 'delta_cut', 1);
plot(hAx, rho, delta, 'k.', rho(icl), delta(icl), 'r.');
axis(hAx, [0 .1 0 10]);
plot(hAx, get(hAx,'XLim'), delta_cut*[1,1], 'r-');
xlabel(hAx, 'rho'); ylabel(hAx, 'delta');
grid(hAx, 'on'); 
set_fig_(hFig, S_fig);
end %func


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Type Figure added
function hFig = set_fig_(vcTag, S_fig)
% return figure handle based on the tag
% hFig = set_fig_(vcTag, S_fig)
% hFig = set_fig_(hFig, S_fig)
hFig = [];
try
    if ischar(vcTag)
        hFig = findobj('Tag', vcTag, 'Type', 'Figure');
    else
        hFig = vcTag;
    end
    set(hFig, 'UserData', S_fig); %figure property
catch
     disperr_();
end
end %end


%--------------------------------------------------------------------------
function S_manual = create_S_manual_(S0)
% S_auto is copied to S_manual and will change. S_auto not referred anymore

% subsample fet
[ctrFet_sub_clu, cviSpk_sub_clu] = clu_fet_(S0);
[trWav_clu, trWav_full_clu] = clu_mean_wav_(S0, ctrFet_sub_clu);
mrWavCor = clu_wavcor_(S0, trWav_full_clu);
csNote_clu = cell(S0.S_auto.nClu, 1);
S_manual = S0.S_auto;
S_manual = struct_add_(S_manual, ctrFet_sub_clu, cviSpk_sub_clu, ...
    trWav_clu, trWav_full_clu, mrWavCor, csNote_clu);
end %func


%--------------------------------------------------------------------------
function S = struct_add_(varargin)
% S = struct_add_(S, var1, var2, ...)
% output
% S.var1=var1; S.var2=var2; ...

S = varargin{1};
for i=2:numel(varargin)
    try
        S.(inputname(i)) = varargin{i};
    catch
        disperr_();
    end
end
end %func


%--------------------------------------------------------------------------
function [ctrFet_sub_clu, cviSpk_sub_clu] = clu_fet_(S0, max_sample)
if nargin<2, max_sample = []; end
if isempty(max_sample), max_sample = 2^12; end

fUseSecondSite = 0;
t_fun = tic;
[S_auto, P] = get_(S0, 'S_auto', 'P');
nSites = numel(P.viSite2Chan);
[viClu_spk, viSite_clu, nClu] = get_(S_auto, 'viClu', 'viSite_clu', 'nClu');
[viSite_spk, viSite2_spk] = get_(S0, 'viSite_spk', 'viSite2_spk');
% nSites = numel(P.viSite2Chan);
% for each cluster create cluster average waveform and subsample features
if isempty(viSite2_spk) || ~fUseSecondSite
    cviSite_spk_fet = {viSite_spk};
else
    cviSite_spk_fet = {viSite_spk, viSite2_spk};
end
[ctrFet_sub_clu, cviSpk_sub_clu] = deal(cell(1, nClu));

fprintf('clu_fet_: ');
for iFet = 1:numel(cviSite_spk_fet)
    for iSite = 1:nSites  
        [trPc_spk1, viSpk1] = load_fet_site_(S0, iFet, iSite);
        viClu_site1 = find(viSite_clu == iSite);
        if isempty(viClu_site1), continue; end % no cluster found
        viClu1 = viClu_spk(viSpk1);              
        cviiSpk1_clu1 = arrayfun_(@(x)subsample_vr_(find(viClu1(:)==x), max_sample), viClu_site1);
        cviSpk1_clu1 = cellfun_(@(x)viSpk1(x), cviiSpk1_clu1);
        ctrFet_clu1 = cellfun_(@(x)trPc_spk1(:,:,x), cviiSpk1_clu1);
        ctrFet_sub_clu(viClu_site1) = join_cell_(ctrFet_sub_clu(viClu_site1), ctrFet_clu1, 3);
        cviSpk_sub_clu(viClu_site1) = join_cell_(cviSpk_sub_clu(viClu_site1), cviSpk1_clu1, 1);
        fprintf('.');
    end
end
fprintf(' took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
% join cells element-wise
function cellC = join_cell_(cellA, cellB, dimm)
assert(numel(cellA)==numel(cellB), 'join_cell_: number of elements must match');
cellC = cellfun_(@(x,y)cat(dimm, x, y), cellA, cellB);
end %func


%--------------------------------------------------------------------------
function S = load_(vcFile_mat)
S = [];
if exist_file_(vcFile_mat)
    try
        S = load(vcFile_mat);
    catch
    end
end
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
function S_fig = figures_manual_(P)
% 'iFig', [], 'name', '', 'pos', [], 'fToolbar', 0, 'fMenubar', 0);
create_figure_('FigMap', [0 .25 .15 .75], ['Probe map; ', P.vcFile_prm], 1, 0);  
create_figure_('FigInfo', [0 0 .15 .25], ['Unit info; ', P.vcFile_prm], 1, 0);

create_figure_('FigWav', [.15 .25 .35 .75],['Averaged waveform: ', P.vcFile_prm], 0, 1);
create_figure_('FigTime', [.15 0 .7 .25], ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);

create_figure_('FigProj', [.5 .25 .35 .5], ['Feature projection: ', P.vcFile_prm]);
create_figure_('FigWavCor', [.5 .75 .35 .25], ['Waveform correlation (click): ', P.vcFile_prm]);

create_figure_('FigHist', [.85 .75 .15 .25], ['ISI Histogram: ', P.vcFile_prm]);
create_figure_('FigIsi', [.85 .5 .15 .25], ['Return map: ', P.vcFile_prm]);
create_figure_('FigCorr', [.85 .25 .15 .25], ['Time correlation: ', P.vcFile_prm]);
create_figure_('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.vcFile_prm]);

% drawnow_();
csFig = {'FigMap', 'FigInfo', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
cvrFigPos0 = cellfun(@(vc)get(get_fig_(vc), 'OuterPosition'), csFig, 'UniformOutput', 0);
S_fig = cell2struct(cvrFigPos0, csFig, 2);
end %func


%--------------------------------------------------------------------------
function varargout = cellfun_(varargin)
if nargout == 0
    cellfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==1
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
function varargout = arrayfun_(varargin)
if nargout == 0
    arrayfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==1
    varargout{1} = arrayfun(varargin{:}, 'UniformOutput', 0);
elseif nargout==2
    [varargout{1}, varargout{2}] = arrayfun(varargin{:}, 'UniformOutput', 0);    
elseif nargout==3
    [varargout{1}, varargout{2}, varargout{3}] = arrayfun(varargin{:}, 'UniformOutput', 0);    
else
    error('arrayfun_: nargout exceeds 3');
end   
end %func


%--------------------------------------------------------------------------
function drawnow_()
try
    drawnow limitrate;
catch
    drawnow();
end
end %func


%--------------------------------------------------------------------------
function auto_scale_proj_time_(S0, fPlot)
% auto-scale and refgresh


if nargin<2, fPlot = 0; end

autoscale_pct = get_set_(S0.P, 'autoscale_pct', 99.5);
[hFig_proj, S_fig_proj] = get_fig_cache_('FigProj');
[mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, S_fig_proj.viSites_show);
if isempty(mrMin2) || isempty(mrMax2)
    cmrAmp = {mrMin1, mrMax1};
else
    cmrAmp = {mrMin1, mrMax1, mrMin2, mrMax2};
end
S_fig_proj.maxAmp = max(cellfun(@(x)quantile(x(:), autoscale_pct/100), cmrAmp));
set(hFig_proj, 'UserData', S_fig_proj);


% Update time
[hFig_time, S_fig_time] = get_fig_cache_('FigTime');
iSite = S0.S_manual.viSite_clu(S0.iCluCopy);
[vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, S0.iCluCopy, S0); % plot iCluCopy
if isempty(S0.iCluPaste)
    cvrFet = {vrFet1};
else
    [vrFet2, vrTime2, vcYlabel, viSpk2] = getFet_site_(iSite, S0.iCluPaste, S0); % plot iCluCopy
    cvrFet = {vrFet1, vrFet2};
end
% S_fig_time.maxAmp = quantile(vrFet, autoscale_pct/100);
S_fig_time.maxAmp = max(cellfun(@(x)quantile(x(:), autoscale_pct/100), cvrFet));
set(hFig_time, 'UserData', S_fig_time);

% plot
if fPlot
    keyPressFcn_cell_(get_fig_cache_('FigWav'), {'j', 't'}, S0); 
else
    rescale_FigProj_(S_fig_proj.maxAmp, hFig_proj, S_fig_proj, S0);    
    rescale_FigTime_(S_fig_time.maxAmp, S0);
end
end %func


%--------------------------------------------------------------------------
function [hFig, S_fig] = plot_FigWavCor_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P;
nClu = S0.S_manual.nClu;
[hFig, S_fig] = get_fig_cache_('FigWavCor'); 
mrWavCor = S0.S_manual.mrWavCor;

figure_wait_(1, hFig);
% Plot
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    set(S_fig.hAx, 'Position', [.1 .1 .8 .8], 'XLimMode', 'manual', 'YLimMode', 'manual', 'Layer', 'top');
    set(S_fig.hAx, {'XTick', 'YTick'}, {1:nClu, 1:nClu});
    axis_(S_fig.hAx, [0 nClu 0 nClu]+.5);
    axis(S_fig.hAx, 'xy');
    grid(S_fig.hAx, 'on');
    xlabel(S_fig.hAx, 'Clu#'); 
    ylabel(S_fig.hAx, 'Clu#');
    S_fig.hImWavCor = imagesc(mrWavCor, P.corrLim);  %clears title and current figure
    S_fig.hCursorV = line([1 1], [.5 nClu+.5], 'Color', [0 0 0], 'LineWidth', 1.5); 
    S_fig.hCursorH = line([.5 nClu+.5], [1 1], 'Color', [1 0 0], 'LineWidth', 1.5);             
    colorbar(S_fig.hAx);
    S_fig.vcTitle = '[S]plit; [M]erge; [D]elete';
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigWavCor_);
    mouse_figure(hFig, S_fig.hAx, @button_FigWavCor_);
    S_fig.hDiag = plotDiag_([0, nClu, .5], 'Color', [0 0 0], 'LineWidth', 1.5);
else
    set(S_fig.hImWavCor, 'CData', mrWavCor);
    set(S_fig.hCursorV, 'xdata', [1 1], 'ydata', [.5 nClu+.5]);
    set(S_fig.hCursorH, 'xdata', .5+[0 nClu], 'ydata', [1 1]);
end
% output
set(hFig, 'UserData', S_fig);
figure_wait_(0, hFig);
end %func


%--------------------------------------------------------------------------
function S0 = button_CluWav_simulate_(iCluCopy, iCluPaste, S0)
if nargin<3,  S0 = get(0, 'UserData'); end
if nargin<2, iCluPaste = []; end
if iCluCopy == iCluPaste, iCluPaste = []; end
hFig_wait = figure_wait_(1);

S0 = update_cursor_(S0, iCluCopy, 0);
S0 = update_cursor_(S0, iCluPaste, 1);
S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'j','t','c','i','v','e','f'}, S0); %'z' to recenter

auto_scale_proj_time_(S0);
figure_wait_(0, hFig_wait);

plot_raster_(S0); %psth
ui_show_elective_();
if nargout==0, set(0, 'UserData', S0); end
end


%--------------------------------------------------------------------------
function button_CluWav_(xyPos, vcButton)
if strcmpi(vcButton, 'normal')
    event.Button = 1;
elseif strcmpi(vcButton, 'alt')
    event.Button = 3;
else
    return;
end
xPos = round(xyPos(1));
S0 = get(0, 'UserData');
switch(event.Button)
    case 1 %left click. copy clu and delete existing one
        S0 = update_cursor_(S0, xPos, 0);
    case 2 %middle, ignore
        return; 
    case 3 %right click. paste clu
        S0 = update_cursor_(S0, xPos, 1);
end
hFig_wait = figure_wait_(1);
S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'j','t','c','i','v','e','f'}, S0); %'z'
auto_scale_proj_time_(S0);
set(0, 'UserData', S0);
plot_raster_(S0);
ui_show_elective_();
figure_wait_(0, hFig_wait);
end %func


%--------------------------------------------------------------------------
function [trWav_clu, trWav_full_clu] = clu_mean_wav_(S0, ctrFet_clu)

trWav_clu = cellfun_(@(x)pc2wav_(S0.mrPv_global, mean(x,3)), ctrFet_clu);
trWav_clu = cat(3, trWav_clu{:});  

if nargout>=2
    [S_auto, P] = get_(S0, 'S_auto', 'P');
    dimm = size(trWav_clu); dimm(2) = numel(P.viSite2Chan);
    trWav_full_clu = zeros(dimm, 'like', trWav_clu);
    nClu = dimm(3);
    sites_clu_ = @(x)P.miSites(:, S_auto.viSite_clu(x));
    for iClu = 1:nClu
        trWav_full_clu(:,sites_clu_(iClu),iClu) = trWav_clu(:,:,iClu);
    end
end
end %func


%--------------------------------------------------------------------------
function mrWavCor = clu_wavcor_(S0, trWav_full_clu)

[S_auto, P] = get_(S0, 'S_auto', 'P');
mrWavCor = zeros(S_auto.nClu, 'single');
nShift = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
for iClu=1:S_auto.nClu
    mrWavCor(:,iClu) = wavcor_(trWav_full_clu(:,:,iClu), trWav_full_clu, nShift);
end
mrWavCor = max(mrWavCor, mrWavCor');
end %func


%--------------------------------------------------------------------------
function S_fig = plot_tnWav_clu_(S_fig, P, S0)
% Substituting plot_spk_

if ~isfield(P, 'LineWidth'), P.LineWidth=1; end
trWav_clu = S0.S_manual.trWav_clu;
[nSamples, nSites, nClu] = size(trWav_clu);
nChans_show = size(P.miSites, 1);
sites_clu_ = @(x)P.miSites(:, S0.S_manual.viSite_clu(x));

% determine x
x_offset = P.spkLim(2) / (diff(P.spkLim)+1); %same for raw and filt
vrX = (1:nSamples*nClu)/nSamples + x_offset;
vrX(1:nSamples:end) = nan;
vrX(nSamples:nSamples:end) = nan;
trWav_clu = trWav_clu / S_fig.maxAmp;

% nChans_show = size(P.miSites,1);
mrX = repmat(vrX(:), [1, nChans_show]);
mrX = reshape(mrX, [nSamples, nClu, nChans_show]);
mrX = reshape(permute(mrX, [1 3 2]), [nSamples*nChans_show, nClu]);

mrY = zeros(nSamples * nChans_show, nClu, 'single');
for iClu=1:nClu
    viSites1 = sites_clu_(iClu);
    mrY1 = trWav_clu(:,:,iClu);
    mrY1 = bsxfun(@plus, mrY1, single(viSites1'));
    mrY(:,iClu) = mrY1(:);
end
    
% if isempty(P.LineStyle)
if isfield(S_fig, 'vhPlot')
    plot_update_(S_fig.vhPlot, mrX, mrY); 
else
    S_fig.vhPlot = plot_group_(S_fig.hAx, mrX, mrY, 'LineWidth', P.LineWidth); 
end
end %func


%--------------------------------------------------------------------------
% 11/19/2018 JJJ: version-neutral uimenu
function h = uimenu_(varargin)
cell_args = varargin;
hFig = varargin{1};
position0 = [];
try
    if isa(hFig, 'matlab.ui.Figure')
        position0 = hFig.OuterPosition;
    end
catch
    ;
end
% csFields = varargin(2:2:end);
if version_matlab_() >= 2017.5 % new version calls Label as Text
    cell_args = cellstr_subs_(cell_args, 'Label', 'Text', 1);
else
    cell_args = cellstr_subs_(cell_args, 'Text', 'Label', 1);
end
h = uimenu(cell_args{:});
if ~isempty(position0)
    drawnow;
    hFig.OuterPosition = position0;
end
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
function vc = if_on_off_(vc, cs)
if ischar(cs), cs = {cs}; end
vc = ifeq_(ismember(vc, cs), 'on', 'off');
end %func


%---------------------------------------------------------------------------
function out = ifeq_(if_, true_, false_)
if (if_)
    out = true_;
else
    out = false_;
end
end %func


%--------------------------------------------------------------------------
% 12/21/17 JJJ: Get the tag by name which is cached (like hash table)
function hObj = get_tag_(vcTag, vcType)
% clear before starting manual
% Return from persistent cache
% Create a new figure if Tag doesn't exist

persistent S_tag_cache_
if nargin<2, vcType = []; end
if isempty(S_tag_cache_)
    S_tag_cache_ = struct(); 
else
    if isfield(S_tag_cache_, vcTag)
        hObj = S_tag_cache_.(vcTag);
        if isvalid_(hObj), return; end
    end
end
hObj = findobj('Tag', vcTag, 'Type', vcType);
S_tag_cache_.(vcTag) = hObj;
end %func


%--------------------------------------------------------------------------
function cell_out = call_irc2_(dbstack1, cell_input, nargout)
vcFunc = dbstack1(1).name;
try
    switch nargout
        case 0, cell_out{1} = []; irc2('call', vcFunc, cell_input);
        case 1, cell_out{1} = irc2('call', vcFunc, cell_input);
        case 2, [cell_out{1}, cell_out{2}] = irc2('call', vcFunc, cell_input);
        case 3, [cell_out{1}, cell_out{2}, cell_out{3}] = irc2('call', vcFunc, cell_input);
        case 4, [cell_out{1}, cell_out{2}, cell_out{3}, cell_out{4}] = irc2('call', vcFunc, cell_input);
        otherwise, error('call_irc2_: undefined func: %s', vcFunc);
    end
catch ME
    fprintf(2, 'call_irc2_: %s\n', ME.message);
    rethrow(ME);
end
end %func


%--------------------------------------------------------------------------
function cell_out = call_irc_(dbstack1, cell_input, nargout)
vcFunc = dbstack1(1).name;
try
    switch nargout
        case 0, cell_out{1} = []; irc('call', vcFunc, cell_input);
        case 1, cell_out{1} = irc('call', vcFunc, cell_input);
        case 2, [cell_out{1}, cell_out{2}] = irc('call', vcFunc, cell_input);
        case 3, [cell_out{1}, cell_out{2}, cell_out{3}] = irc('call', vcFunc, cell_input);
        case 4, [cell_out{1}, cell_out{2}, cell_out{3}, cell_out{4}] = irc('call', vcFunc, cell_input);
        otherwise, error('call_irc_: undefined func: %s', vcFunc);
    end
catch ME
    fprintf(2, 'call_irc_: %s\n', ME.message);
    rethrow(ME);
end
end %func


%--------------------------------------------------------------------------
function viTime1 = randomSelect_(viTime1, nShow)
if isempty(viTime1), return; end
if numel(viTime1) > nShow
    viTime1 = viTime1(randperm(numel(viTime1), nShow));
end
end %func


%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
end %func


%--------------------------------------------------------------------------
function vi = subsample_vr_(vi, nMax)
if numel(vi)>nMax
    nSkip = floor(numel(vi)/nMax);
    if nSkip>1, vi = vi(1:nSkip:end); end
    if numel(vi)>nMax
        try
            nRemove = numel(vi) - nMax;
            viRemove = round(linspace(1, numel(vi), nRemove));
            viRemove = min(max(viRemove, 1), numel(vi));
            vi(viRemove) = [];
        catch
            vi = vi(1:nMax);
        end
    end
end
end %func


%--------------------------------------------------------------------------
function vhFig = figure_wait_(fWait, vhFig)
% set all figures pointers to watch
if nargin<2, vhFig = gcf; end
fWait_prev = strcmpi(get_(vhFig(1), 'Pointer'), 'watch');
if fWait && ~fWait_prev      
    set_(vhFig, 'Pointer', 'watch'); 
    drawnow_();
else
    set_(vhFig, 'Pointer', 'arrow');
end
end %func


%--------------------------------------------------------------------------
function vc = set_(vc, varargin)
% Set handle to certain values
% set_(S, name1, val1, name2, val2)

if isempty(vc), return; end
if isstruct(vc)
    for i=1:2:numel(varargin)        
        vc.(varargin{i}) = varargin{i+1};
    end
    return;
end
if iscell(vc)
    for i=1:numel(vc)
        try
            set(vc{i}, varargin{:});
        catch
        end
    end
elseif numel(vc)>1
    for i=1:numel(vc)
        try
            set(vc(i), varargin{:});
        catch
        end
    end
else
    try
        set(vc, varargin{:});
    catch
    end 
end
end %func


%--------------------------------------------------------------------------
function axis_(arg1, arg2)
if nargin==1
    [hAx_, lim_] = deal(gca, arg1);
else
    [hAx_, lim_] = deal(arg1, arg2);
end
if ischar(lim_)
    axis(hAx_, lim_);
    return;
end
if any(isnan(lim_)), return; end
try
    axis(hAx_, [sort(lim_(1:2)), sort(lim_(3:4))]);
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
% 8/9/17 JJJ: Generalized to any figure objects
function S0 = keyPressFcn_cell_(hObject, csKey, S0)
% Simulate key press function

if nargin<3, S0 = get(0, 'UserData'); end
% figure_wait_(1); 
event1.Key = '';
if ischar(csKey), csKey = {csKey}; end
nKeys = numel(csKey);
keyPressFcn_ = get(hObject, 'KeyPressFcn');
for i=1:nKeys
    event1.Key = csKey{i};
    S0 = keyPressFcn_(hObject, event1, S0);
end
if nargout==0, set(0, 'UserData', S0); end
end %func


%--------------------------------------------------------------------------
function S0 = update_cursor_(S0, iClu, fPaste)
if isempty(iClu), return; end
if isempty(S0), S0 = get(0, 'UserData'); end
nClu = S0.S_manual.nClu;
% [hFig, S_fig] = get_fig_cache_('FigWav');

if ~isfield(S0, 'hCopy'), S0.hCopy = []; end
if ~isfield(S0, 'hPaste'), S0.hPaste = []; end

if ~fPaste
    iCluCopy = iClu;    
    if iCluCopy <1 || iCluCopy > nClu, return; end
    update_plot_(S0.hPaste, nan, nan); %hide paste
    S0.iCluPaste = []; 
    [S0.iCluCopy, S0.hCopy] = plot_tmrWav_clu_(S0, iCluCopy, S0.hCopy, [0 0 0]);
else
    iCluPaste = iClu;    
    if iCluPaste < 1 || iCluPaste > nClu || S0.iCluCopy == iCluPaste, return; end
    [S0.iCluPaste, S0.hPaste] = plot_tmrWav_clu_(S0, iCluPaste, S0.hPaste, [1 0 0]);
end
% set(hFig, 'UserData', S_fig);
cursor_FigWavCor_(S0);
if nargout==0, set(0, 'UserData', S0); end
end %func


%--------------------------------------------------------------------------
function cursor_FigWavCor_(S0)
if nargin==0, S0 = get(0, 'UserData'); end
[mrWavCor, nClu] = get_(S0.S_manual ,'mrWavCor', 'nClu');

[hFig, S_fig] = get_fig_cache_('FigWavCor');
if isempty(S_fig)
    [hFig, S_fig] = plot_FigWavCor_(S0);
end
iClu1 = S0.iCluCopy;
if isempty(S0.iCluPaste)
    iClu2 = S0.iCluCopy;
else
    iClu2 = S0.iCluPaste;
end

cor12 = mrWavCor(iClu1, iClu2);
set(S_fig.hCursorV, 'XData', iClu1*[1,1], 'YData', [.5, nClu+.5]);
title_(S_fig.hAx, sprintf('Clu%d vs. Clu%d: %0.3f; %s', iClu1, iClu2, cor12, S_fig.vcTitle));    
if iClu1==iClu2, color_H = [0 0 0]; else color_H = [1 0 0]; end
set(S_fig.hCursorH, 'YData', iClu2*[1,1], 'XData', [.5, nClu+.5], 'Color', color_H);
xlim_(S_fig.hAx, trim_lim_(iClu1 + [-6,6], [.5, nClu+.5]));
ylim_(S_fig.hAx, trim_lim_(iClu2 + [-6,6], [.5, nClu+.5]));
end %func


%--------------------------------------------------------------------------
function update_plot_(hPlot, vrX, vrY, S_plot)
% update the plot with new x and y 

if nargin<4, S_plot = []; end
if isempty(hPlot), return; end
% selective plot to speed up plotting speed
if isempty(vrY) || isempty(vrX)
    hide_plot_(hPlot);
%     set(hPlot, 'XData', nan, 'YData', nan);  %visible off
    return;
end

% only update if both x and y are changed
vrX1 = get(hPlot, 'XData');
vrY1 = get(hPlot, 'YData');
fUpdate = 1;
if (numel(vrX1) == numel(vrX)) && (numel(vrY1) == numel(vrY))
    if (std(vrX1(:) - vrX(:)) == 0) && (std(vrY1(:) - vrY(:)) == 0)
        fUpdate = 0; 
    end
end
if fUpdate, set(hPlot, 'xdata', vrX, 'ydata', vrY); end
if ~isempty(S_plot), set(hPlot, 'UserData', S_plot); end
end %func


%--------------------------------------------------------------------------
function [iClu, hPlot] = plot_tmrWav_clu_(S0, iClu, hPlot, vrColor)
[S_manual, P] = get_(S0, 'S_manual', 'P');
[hFig, S_fig] = get_fig_cache_('FigWav');
if ~isvalid_(hPlot)
    hPlot = plot_(nan, nan, 'Color', vrColor, 'LineWidth', 2, 'Parent', S_fig.hAx);
end
mrWav_clu1 = S_manual.trWav_full_clu(:,:, iClu);
multiplot(hPlot, S_fig.maxAmp, wav_clu_x_(iClu, P), mrWav_clu1);
uistack_(hPlot, 'top');
end %func


%--------------------------------------------------------------------------
function vrX = wav_clu_x_(iClu, P)
% determine x range of a cluster
if P.fWav_raw_show    
    spkLim = P.spkLim_raw;
    dimm_raw = get0_('dimm_raw');
    if dimm_raw(1) ~= diff(spkLim)+1, spkLim = P.spkLim * 2; end %old format
else
    spkLim = P.spkLim;
end
nSamples = diff(spkLim) + 1;
x_offset = spkLim(2) / nSamples + iClu - 1;

vrX = (1:nSamples) / nSamples + x_offset;
vrX([1,end]) = nan;
vrX = single(vrX(:));
end %func


%--------------------------------------------------------------------------
function uistack_(h, vc)
try
    uistack(h, vc);
catch
end
end %func


%--------------------------------------------------------------------------
function flag = isvalid_(h)
if isempty(h), flag = 0; return; end
try
    flag = isvalid(h);
catch
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
function hTitle = title_(hAx, vc)
% title_(vc)
% title_(hAx, vc)

if nargin==1, vc=hAx; hAx=[]; end
% Set figure title

if isempty(hAx), hAx = gca; end
hTitle = get_(hAx, 'Title');
if isempty(hTitle)
    hTitle = title(hAx, vc, 'Interpreter', 'none', 'FontWeight', 'normal');
else
    set_(hTitle, 'String', vc, 'Interpreter', 'none', 'FontWeight', 'normal');
end
end %func


%--------------------------------------------------------------------------
function ylim_(arg1, arg2)
% ylim function
% ylim_(lim_)
% ylim_(hAx, lim_)
if nargin==1
    [hAx_, lim_] = deal(gca, arg1);
else
    [hAx_, lim_] = deal(arg1, arg2);
end
if any(isnan(lim_)), return; end
try
    ylim(hAx_, sort(lim_));
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
function xlim_(arg1, arg2)
% ylim function
% ylim_(lim_)
% ylim_(hAx, lim_)
if nargin==1
    [hAx_, lim_] = deal(gca, arg1);
else
    [hAx_, lim_] = deal(arg1, arg2);
end
if any(isnan(lim_)), return; end
try
    xlim(hAx_, sort(lim_));
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
% find intersection of two limit ranges
function xlim1 = trim_lim_(xlim1, xlim0)
dx = diff(xlim1);

if xlim1(1)<xlim0(1), xlim1 = xlim0(1) + [0, dx]; end
if xlim1(2)>xlim0(2), xlim1 = xlim0(2) + [-dx, 0]; end
xlim1(1) = max(xlim1(1), xlim0(1));
xlim1(2) = min(xlim1(2), xlim0(2));
end %func


%--------------------------------------------------------------------------
function hide_plot_(vhPlot)
for i=1:numel(vhPlot), set(vhPlot(i), 'XData', nan, 'YData', nan); end
end %func


%--------------------------------------------------------------------------
function plot_FigMap_(S0)
if nargin<1, S0 = get(0, 'UserData'); end 
P = S0.P; 
trWav_full_clu = S0.S_manual.trWav_full_clu;
[hFig, S_fig] = get_fig_cache_('FigMap');

mrWav1 = trWav_full_clu(:,:,S0.iCluCopy);
vrVpp = squeeze_(max(mrWav1) - min(mrWav1));
mrVpp = repmat(vrVpp(:)', [4, 1]);
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    [S_fig.mrPatchX, S_fig.mrPatchY] = probe_map_(P);
    S_fig.hPatch = patch(S_fig.mrPatchX, S_fig.mrPatchY, mrVpp, ...
        'EdgeColor', 'k', 'FaceColor', 'flat');
    S_fig.alim = [min(S_fig.mrPatchX(:)), max(S_fig.mrPatchX(:)), min(S_fig.mrPatchY(:)), max(S_fig.mrPatchY(:))];
    S_fig.cell_alim = get_lim_shank_(P);
    colormap jet;
    mouse_figure(hFig);   
    nSites = size(P.mrSiteXY,1);    
    csText = arrayfun(@(i)sprintf('%d', i), 1:nSites, 'UniformOutput', 0);
    S_fig.hText = text(P.mrSiteXY(:,1), P.mrSiteXY(:,2), csText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');    
    xlabel('X Position (\mum)');
    ylabel('Y Position (\mum)');
else    
    set(S_fig.hPatch, 'CData', mrVpp);    
end
try
    iShank1 = P.viShank_site(S_manual.viSite_clu(S0.iCluCopy));
    axis_(S_fig.hAx, S_fig.cell_alim{iShank1});
    vcTitle = sprintf('Max: %0.1f uVpp (Shank %d)', max(vrVpp), iShank1);
catch
    axis_(S_fig.hAx, S_fig.alim);
    vcTitle = sprintf('Max: %0.1f uVpp', max(vrVpp));
end
title_(S_fig.hAx, vcTitle);
caxis(S_fig.hAx, [0, max(vrVpp)]);

set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
% Remove leading singular dimension
% 12/15/17 JJJ: squeeze out specific dimension
% 7/26/17 JJJ: code cleanup and testing
function val = squeeze_(val, idimm)
% val = squeeze_(val) : when squeezeing matrix, transpose if leading dimm is 1
% val = squeeze_(val, idimm): permute specified dimension out
size_ = size(val);
if nargin>=2
    dimm_ = [setdiff(1:ndims(val), idimm), idimm];
    val = permute(val, dimm_);
elseif numel(size_)==2 && size_(1) == 1
    val = val';
else
    val = squeeze(val);
end
end


%--------------------------------------------------------------------------
function cell_alim = get_lim_shank_(P)
vrSiteHW = get_set_(P, 'vrSiteHW', [12,12]);
nSites = size(P.mrSiteXY,1);
viShank_site = get_set_(P, 'viShank_site', ones(nSites,1));
if isempty(viShank_site), viShank_site = ones(nSites,1); end

viShank_unique = 1:max(viShank_site);
cell_alim = cell(size(viShank_unique));
[dx, dy] = deal(abs(vrSiteHW(2)), abs(vrSiteHW(1)));
for iShank = 1:numel(viShank_unique)
    viSite1 = find(viShank_site == iShank);
    vrX1 = P.mrSiteXY(viSite1,1);
    vrY1 = P.mrSiteXY(viSite1,2);    
    xlim1 = [min(vrX1), max(vrX1)] + [-1,2] * dx;
    ylim1 = [min(vrY1), max(vrY1)] + [-1,2] * dy;
    cell_alim{iShank} = [xlim1, ylim1];
end
end %func


%--------------------------------------------------------------------------
function [mrPatchX, mrPatchY] = probe_map_(P)
vrSiteHW = get_set_(P, 'vrSiteHW', [12, 12]);
vrX = [0 0 1 1] * vrSiteHW(2); 
vrY = [0 1 1 0] * vrSiteHW(1);
mrPatchX = bsxfun(@plus, P.mrSiteXY(:,1)', vrX(:));
mrPatchY = bsxfun(@plus, P.mrSiteXY(:,2)', vrY(:));
end %func


%--------------------------------------------------------------------------
function plot_FigHist_(S0)

if nargin<1, S0 = get(0, 'UserData'); end
[hFig, S_fig] = get_fig_cache_('FigHist');

nBins_hist = 50; % @TODO: put this in param file
P = S0.P;
vrX = logspace(0, 4, nBins_hist);
vrY1 = isi_hist_(S0, S0.iCluCopy, vrX); 
vcTitle = sprintf('Cluster %d', S0.iCluCopy);

% draw
if isempty(S_fig) %first time the iCluPaste is always empty
    S_fig.hAx = axes_new_(hFig);
    S_fig.hPlot1 = stairs(S_fig.hAx, nan, nan, 'k'); 
    S_fig.hPlot2 = stairs(S_fig.hAx, nan, nan, 'r');     
    xlim_(S_fig.hAx, [1 10000]); %in msec
    grid(S_fig.hAx, 'on');
    xlabel(S_fig.hAx, 'ISI (ms)');
    ylabel(S_fig.hAx, 'Prob. Density');
    set(S_fig.hAx, 'XScale', 'log');
end
update_plot_(S_fig.hPlot1, vrX, vrY1);
if ~isempty(S0.iCluPaste)
    vrY2 = isi_hist_(S0, S0.iCluPaste, vrX);
    vcTitle = sprintf('Cluster %d (black) vs %d (red)', S0.iCluCopy, S0.iCluPaste);
    update_plot_(S_fig.hPlot2, vrX, vrY2);
else
    update_plot_(S_fig.hPlot2, nan, nan);
end
title_(S_fig.hAx, vcTitle);

set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function vnHist = isi_hist_(S0, iClu1, vrX)

vrTime1 = double(clu_time_(S0, iClu1)) / S0.P.sRateHz;
vnHist = hist(diff(vrTime1)*1000, vrX);
vnHist(end)=0;
vnHist = vnHist ./ sum(vnHist);
end


%--------------------------------------------------------------------------
function plot_FigIsi_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P; 
[hFig, S_fig] = get_fig_cache_('FigIsi');

[vrX1, vrY1] = get_returnMap_(S0, S0.iCluCopy);
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    S_fig.hPlot1 = plot_(S_fig.hAx, nan, nan, 'ko');
    S_fig.hPlot2 = plot_(S_fig.hAx, nan, nan, 'ro');
    set(S_fig.hAx, 'XScale','log', 'YScale','log');   
    xlabel('ISI_{k} (ms)'); ylabel('ISI_{k+1} (ms)');
    axis_(S_fig.hAx, [1 10000 1 10000]);
    grid(S_fig.hAx, 'on');
    % show refractory line
    line(get(S_fig.hAx,'XLim'), P.spkRefrac_ms*[1 1], 'Color', [1 0 0]);
    line(P.spkRefrac_ms*[1 1], get(S_fig.hAx,'YLim'), 'Color', [1 0 0]);
end  
update_plot_(S_fig.hPlot1, vrX1, vrY1);
if ~isempty(S0.iCluPaste)    
    [vrX2, vrY2] = get_returnMap_(S0, S0.iCluPaste);
    update_plot_(S_fig.hPlot2, vrX2, vrY2);
else
    update_plot_(S_fig.hPlot2, nan, nan);
end

set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function [vrX, vrY] = get_returnMap_(S0, iClu)
P = S0.P;
vrTime1 = double(clu_time_(S0, iClu)) / P.sRateHz;
vrIsi1 = diff(vrTime1 * 1000); % in msec
vrX = vrIsi1(1:end-1);
vrY = vrIsi1(2:end);
viShow = randperm(numel(vrX), min(P.nShow, numel(vrX)));
vrX = vrX(viShow);
vrY = vrY(viShow);
end


%--------------------------------------------------------------------------
function [viTime_clu1, viSpk_clu1] = clu_time_(S0, iClu1)
% returns time in sec
[S_manual, viTime_spk] = get_(S0, 'S_manual', 'viTime_spk');
viSpk_clu1 = S_manual.cviSpk_clu{iClu1};
viTime_clu1 = viTime_spk(S_manual.cviSpk_clu{iClu1});
end %func


%--------------------------------------------------------------------------
function plot_FigProj_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P;
viSite_clu = S0.S_manual.viSite_clu;
[hFig, S_fig] = get_fig_cache_('FigProj');

iClu1 = S0.iCluCopy;
iClu2 = S0.iCluPaste;
update_plot2_proj_(); %erase prev objects

%---------------
% Compute
iSite1 = viSite_clu(iClu1);
% miSites = P.miSites;
if ~isfield(P, 'viSites_show')
    P.viSites_show = sort(P.miSites(1:P.nSites_fet, iSite1), 'ascend');
end
viSites_show = P.viSites_show;
nSites = numel(P.viSites_show);
cell_plot = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};
switch lower(P.vcFet_show)
    case {'vpp', 'vmin', 'vmax'}
        vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        vcYLabel = 'Site # (%0.0f \\muV_{min})';
    otherwise
        vcXLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', P.vcFet_show, P.vcFet_show, P.vcFet_show);
        vcYLabel = sprintf('Site # (%%0.0f %s)', P.vcFet_show);    
end
vcTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

%----------------
% display
if isempty(S_fig)
    S_fig.maxAmp = P.maxAmp;    
    S_fig.hAx = axes_new_(hFig);
    set(S_fig.hAx, 'Position', [.1 .1 .85 .85], 'XLimMode', 'manual', 'YLimMode', 'manual');
    S_fig.hPlot0 = line(nan, nan, 'Color', P.mrColor_proj(1,:), 'Parent', S_fig.hAx);
    S_fig.hPlot1 = line(nan, nan, 'Color', P.mrColor_proj(2,:), 'Parent', S_fig.hAx); %place holder
    S_fig.hPlot2 = line(nan, nan, 'Color', P.mrColor_proj(3,:), 'Parent', S_fig.hAx); %place holder
    set([S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2], cell_plot{:}); %common style
    S_fig.viSites_show = []; %so that it can update
    S_fig.vcFet_show = 'vpp';
    % plot boundary
    plotTable_([0, nSites], '-', 'Color', [.5 .5 .5]); %plot in one scoop
    plotDiag_([0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5); %plot in one scoop
    mouse_figure(hFig);
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigProj_);
    S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlot0, S_fig);
    set_fig_(hFig, S_fig);
end
% get features for x0,y0,S_plot0 in one go
%[mrMin, mrMax, vi0, vi1, vi2] = fet2proj_(S0, P.viSites_show);
[mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, P.viSites_show);
% S_fig.maxAmp %debug
if ~isfield(S_fig, 'viSites_show'), S_fig.viSites_show = []; end
if ~equal_vr_(S_fig.viSites_show, P.viSites_show) || ...
    ~equal_vr_(S_fig.vcFet_show, P.viSites_show)
    plot_proj_(S_fig.hPlot0, mrMin0, mrMax0, P, S_fig.maxAmp);
end

plot_proj_(S_fig.hPlot1, mrMin1, mrMax1, P, S_fig.maxAmp);
if ~isempty(iClu2)
    plot_proj_(S_fig.hPlot2, mrMin2, mrMax2, P, S_fig.maxAmp);
    vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', iClu1, iClu2, vcTitle);
else
    update_plot_(S_fig.hPlot2, nan, nan);
    vcTitle = sprintf('Clu%d (black); %s', iClu1, vcTitle);
end

% Annotate axes
axis_(S_fig.hAx, [0 nSites 0 nSites]);
set(S_fig.hAx,'XTick',.5:1:nSites,'YTick',.5:1:nSites, 'XTickLabel', P.viSites_show, 'YTickLabel', P.viSites_show, 'Box', 'off');
xlabel(S_fig.hAx, sprintf(vcXLabel, S_fig.maxAmp));   
ylabel(S_fig.hAx, sprintf(vcYLabel, S_fig.maxAmp));  
title_(S_fig.hAx, vcTitle);
vcFet_show = P.vcFet_show;
S_fig = struct_merge_(S_fig, ...
    makeStruct_(vcTitle, iClu1, iClu2, viSites_show, vcXLabel, vcYLabel, vcFet_show));
S_fig.csHelp = { ...
    '[D]raw polygon', ...
    '[S]plit cluster', ...
    '(shift)+Up/Down: change scale', ...
    '[R]eset scale', ...
    'Zoom: mouse wheel', ...
    'Drag while pressing wheel: pan'};
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function keyPressFcn_FigProj_(hFig, event)
S0 = get(0, 'UserData');
[P, S_clu] = deal(S0.P, S0.S_clu);
[hFig, S_fig] = get_fig_cache_('FigProj');
S_plot1 = get(S_fig.hPlot1, 'UserData');
viSites_show = S_plot1.viSites_show;
hFig_wait = figure_wait_(1);
switch lower(event.Key)
    case {'uparrow', 'downarrow'}
        rescale_FigProj_(event, hFig, S_fig, S0);

    case {'leftarrow', 'rightarrow'} % change channels
        fPlot = 0;
        if strcmpi(event.Key, 'leftarrow')
            if min(S_fig.viSites_show)>1
                S_fig.viSites_show=S_fig.viSites_show-1; 
                fPlot = 1;
            end
        else
            if max(S_fig.viSites_show) < max(P.viSite2Chan)
                S_fig.viSites_show=S_fig.viSites_show+1;                 
                fPlot = 1;
            end
        end
        if fPlot
            set(hFig, 'UserData', S_fig);
            S0.P.viSites_show = S_fig.viSites_show;
            plot_FigProj_(S0);
        end
        
    case 'r' %reset view
        axis_([0 numel(viSites_show) 0 numel(viSites_show)]);

    case 's' %split
        figure_wait_(0, hFig_wait);
        if ~isempty(S0.iCluPaste)
            msgbox_('Select one cluster to split'); return;
        end
        S_plot1 = select_polygon_(S_fig.hPlot1); 
        if ~isempty(S_plot1)
            [fSplit, vlIn] = plot_split_(S_plot1);
            if fSplit
                S_clu = split_clu_(S0.iCluCopy, vlIn);
            else
                update_plot2_proj_();
            end
        end
        
    case 'm'
        ui_merge_(S0);
        
    case 'f'
        disp('keyPressFcn_FigProj_: ''f'': not implemented yet');
        
    case 'b' %background spikes
        toggleVisible_(S_fig.hPlot0);

    case 'h' %help
        msgbox_(S_fig.csHelp, 1);
end %switch
figure_wait_(0, hFig_wait);
end %func


%--------------------------------------------------------------------------
function flag = equal_vr_(vr1, vr2)
if all(size(vr1) == size(vr2))
    ml = vr1 == vr2;
    flag = all(ml(:));
else
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
function update_plot2_proj_(vrX, vrY)
if nargin==0, vrX=nan; vrY=nan; end
[hFig, S_fig] = get_fig_cache_('FigProj');
% erase polygon
if nargin==0
    try
        update_plot_(S_fig.hPlot2, vrX, vrY);
        delete(findobj(get(S_fig.hAx, 'Child'), 'Type', 'hggroup'));
    catch
        ;
    end
end
end


%--------------------------------------------------------------------------
function plot_proj_(hPlot, mrMin, mrMax, P, maxAmp)
if nargin<5
    [hFig, S_fig] = get_fig_cache_('FigProj');
    maxAmp = S_fig.maxAmp;
end
[vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, maxAmp, P.maxSite_show, P);

% make struct
maxPair = P.maxSite_show;
viSites_show = P.viSites_show;
S_plot = makeStruct_(mrMax, mrMin, viSites_show, viPlot, tr_dim, maxPair, maxAmp);

update_plot_(hPlot, vrX, vrY, S_plot);
end %func


%--------------------------------------------------------------------------
function [vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, maxAmp, maxPair, P)
if nargin<4, maxPair = []; end
if nargin<5, P = get0_('P'); end
% switch lower(P.vcFet_show)
%     case {'vpp', 'vmin', 'vmax'}
%         mrMax = linmap_(mrMax', [0, maxAmp/2], [0,1], 1);
%         mrMin = linmap_(mrMin', [0, maxAmp], [0,1], 1);
%     otherwise
mrMax = linmap_(mrMax', [0, 1] * maxAmp, [0,1], 1);
mrMin = linmap_(mrMin', [0, 1] * maxAmp, [0,1], 1);            
% end
[nEvt, nChans] = size(mrMin);
if isempty(maxPair), maxPair = nChans; end
[trX, trY] = deal(nan([nEvt, nChans, nChans], 'single'));
for chY = 1:nChans
    vrY1 = mrMin(:,chY);
    vlY1 = vrY1>0 & vrY1<1;
    for chX = 1:nChans
        if abs(chX-chY) > maxPair, continue; end
        if chY > chX
            vrX1 = mrMin(:,chX);
        else
            vrX1 = mrMax(:,chX);
        end
        viPlot1 = find(vrX1>0 & vrX1<1 & vlY1);
        trX(viPlot1,chY,chX) = vrX1(viPlot1) + chX - 1;
        trY(viPlot1,chY,chX) = vrY1(viPlot1) + chY - 1;
    end
end
% plot projection
viPlot = find(~isnan(trX) & ~isnan(trY));
vrX = trX(viPlot);  vrX=vrX(:);
vrY = trY(viPlot);  vrY=vrY(:);
tr_dim = size(trX);
end %func


%--------------------------------------------------------------------------
function vr = linmap_(vr, lim1, lim2, fSat)
if nargin< 4
    fSat = 0;
end
if numel(lim1) == 1, lim1 = [-abs(lim1), abs(lim1)]; end
if numel(lim2) == 1, lim2 = [-abs(lim2), abs(lim2)]; end

if fSat
    vr(vr>lim1(2)) = lim1(2);
    vr(vr<lim1(1)) = lim1(1);
end
if lim1(1)==lim1(2)
    vr = vr / lim1(1);
else
    vr = interp1(lim1, lim2, vr, 'linear', 'extrap');
end
end %func


%--------------------------------------------------------------------------
function [trFet1, viSpk1] = get_fet_clu_(S_manual, iClu1, viSites2)
% viSites: optional site # to project on
if nargin<3, viSites2=[]; end

if isempty(iClu1)
    viClu1 = find(ismember(S_manual.viSite_clu, viSites2));
    [ctrFet1, cviSpk1] = arrayfun_(@(x)get_fet_clu_(S_manual, x, viSites2), viClu1);
    trFet1 = cat_(3, ctrFet1);
    viSpk1 = cat_([], cviSpk1);
    return;
end
iSite1 = S_manual.viSite_clu(iClu1);
viSite1 = S_manual.P.miSites(:,iSite1);
viSpk1 = S_manual.cviSpk_sub_clu{iClu1};
trFet1 = S_manual.ctrFet_sub_clu{iClu1};
if ~isempty(viSites2)
    dimm2 = size(trFet1); dimm2(2) = numel(viSites2);
    trFet2 = zeros(dimm2, 'like', trFet1);
    [~, vii1, vii2] = intersect(viSite1, viSites2);   
    if isempty(vii1), return; end
    for i1 = 1:numel(vii1)
        trFet2(:,vii2(i1),:) = trFet1(:,vii1(i1),:);
    end
    trFet1 = trFet2;
end
end %func


%--------------------------------------------------------------------------
function vr1 = cat_(dimm, cell1)
if isempty(dimm)
    try
        vr1 = cat(1, cell1{:});
    catch
        vr1 = cat(2, cell1{:});
    end
else
    vr1 = cat(dimm, cell1{:});
end
end %func


%--------------------------------------------------------------------------
function [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, viSites0)
% show spikes excluding the clusters excluding clu1 and 2

[S_manual, iClu1, iClu2] = get_(S0, 'S_manual', 'iCluCopy', 'iCluPaste');
trFet_clu0 = get_fet_clu_(S_manual, [], viSites0);
trFet_clu1 = get_fet_clu_(S_manual, iClu1, viSites0);
if ~isempty(iClu2)
    trFet_clu2 = get_fet_clu_(S_manual, iClu2, viSites0);
else
    trFet_clu2 = [];
end
fet2pc_ = @(x,i)permute(x(i,:,:),[2,3,1]);

% put these sites in order
[mrMin0, mrMax0] = deal(fet2pc_(trFet_clu0,1), fet2pc_(trFet_clu0,2)); 
[mrMin1, mrMax1] = deal(fet2pc_(trFet_clu1,1), fet2pc_(trFet_clu1,2)); 
if ~isempty(iClu2)  
    [mrMin2, mrMax2] = deal(fet2pc_(trFet_clu2,1), fet2pc_(trFet_clu2,2)); 
else
    [mrMin2, mrMax2] = deal([]);
end
        
% todo
% switch lower(P.vcFet_show)
%     case 'vpp'
%     case 'vmin'
%         [mrMin0, mrMax0] = getFet_spk_(viSpk00, viSites0, S0); %getall spikes whose center lies in certain range
%         [mrMin1, mrMax1] = getFet_spk_(viSpk01, viSites0, S0); %getall spikes whose center lies in certain range
%         if ~isempty(iClu2)  
%             [mrMin2, mrMax2] = getFet_spk_(viSpk02, viSites0, S0);
%         end            
% end %switch
% [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = ...
%     multifun_(@(x)abs(x), mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2);
end %func


%--------------------------------------------------------------------------
function hPlot = plotTable_(lim, varargin)

vrX = floor((lim(1)*2:lim(2)*2+1)/2);
vrY = repmat([lim(1), lim(2), lim(2), lim(1)], [1, ceil(numel(vrX)/4)]);
vrY = vrY(1:numel(vrX));
hPlot = plot_([vrX(1:end-1), fliplr(vrY)], [vrY(1:end-1), fliplr(vrX)], varargin{:});
end %func


%--------------------------------------------------------------------------
function hPlot = plotDiag_(lim, varargin)
[vrX, vrY] = plotDiag__(lim);
% vrY = floor((lim(1)*2:lim(2)*2+1)/2);
% vrX = [vrY(2:end), lim(end)];
% hPlot = plot_([vrX(1:end-1), fliplr(vrY)], [vrY(1:end-1), fliplr(vrX)], varargin{:});
hPlot = plot_(vrX, vrY, varargin{:});
end %func


%--------------------------------------------------------------------------
function [vrX, vrY] = plotDiag__(lim)
% lim: [start, end] or [start, end, offset]

vrY0 = floor((lim(1)*2:lim(2)*2+1)/2);
% vrY0 = lim(1):lim(end);
vrX0 = [vrY0(2:end), lim(2)];
vrX = [vrX0(1:end-1), fliplr(vrY0)];
vrY = [vrY0(1:end-1), fliplr(vrX0)];
if numel(lim)>=3
    vrX = vrX + lim(3);
    vrY = vrY + lim(3);
end
end %func


%--------------------------------------------------------------------------
function [mrMin, mrMax] = getFet_spk_(viSpk1, viSites1, S0)
% get feature for the spikes of interest

if nargin<3, S0 = get(0, 'UserData'); end
P = S0.P;

switch lower(P.vcFet_show)
    case {'vmin', 'vpp'}   
        tnWav_spk1 = tnWav2uV_(tnWav_spk_sites_(viSpk1, viSites1, S0), P);
        [mrMin, mrMax] = multifun_(@(x)abs(permute(x,[2,3,1])), min(tnWav_spk1), max(tnWav_spk1));
    case {'cov', 'spacetime'}
        [mrMin, mrMax] = calc_cov_spk_(viSpk1, viSites1);
    case 'pca'
        [mrMin, mrMax] = pca_pc_spk_(viSpk1, viSites1); %getall spikes whose center lies in certain range
    otherwise
        error('not implemented yet');
end
end %func


%--------------------------------------------------------------------------
function varargout = multifun_(hFun, varargin)
% apply same function to the input, unary function only

if nargout ~= numel(varargin), error('n arg mismatch'); end
for i=1:nargout
    try
        varargout{i} = hFun(varargin{i});
    catch
        varargout{i} = varargin{i};
    end
end
end %func


%--------------------------------------------------------------------------
% 171201 JJJ: Unique sites handling for diagonal plotting
function tnWav_spk1 = tnWav_spk_sites_(viSpk1, viSites1, S0, fWav_raw_show)
% reorder tnWav1 to viSites1
% P = get0_('P');
% if nargin<3, fWav_raw_show = P.fWav_raw_show; end
if nargin<3, S0 = []; end
if isempty(S0), S0 = get(0, 'UserData'); end
if nargin<4, fWav_raw_show = get_set_(S0.P, 'fWav_raw_show', 0); end

% unique exception handling %171201 JJJ
[viSites1_uniq, ~, viiSites1_uniq] = unique(viSites1);
if numel(viSites1_uniq) ~= numel(viSites1)
    tnWav_spk11 = tnWav_spk_sites_(viSpk1, viSites1_uniq, S0, fWav_raw_show);
    tnWav_spk1 = tnWav_spk11(:,viiSites1_uniq,:);
    return;
end

[viSite_spk, P] = deal(S0.viSite_spk, S0.P);
try
    tnWav = get_spkwav_(P, fWav_raw_show);
    tnWav1 = tnWav(:,:,viSpk1);
    nT_spk = size(tnWav, 1);
catch % decompress from pc
    tnWav1 = pc2wav_(S0.mrPv_global, S0.trPc_spk(:,:,viSpk1));
    nT_spk = size(S0.mrPv_global,1);
end
nSpk1 = numel(viSpk1);
viSites_spk1 = viSite_spk(viSpk1);
viSites_spk_unique = unique(viSites_spk1);
tnWav_spk1 = zeros([nT_spk, numel(viSites1), nSpk1], 'like', tnWav1);
for iSite1 = 1:numel(viSites_spk_unique) %only care about the first site
    iSite11 = viSites_spk_unique(iSite1); %center sites group
    viSpk11 = find(viSites_spk1 == iSite11); %dangerous error
    viSites11 = P.miSites(:, iSite11);        
    [vlA11, viiB11] = ismember(viSites11, viSites1);
    tnWav_spk1(:,viiB11(vlA11),viSpk11) = tnWav1(:,vlA11,viSpk11);
end    
end %func


%--------------------------------------------------------------------------
function clu_info_(S0)
% This also plots cluster position
if nargin<1, S0 = get(0, 'UserData'); end
csNote_clu = S0.S_manual.csNote_clu;
mh_info = get_tag_('mh_info', 'uimenu');
S_clu1 = get_cluInfo_(S0, S0.iCluCopy);
if ~isempty(S0.iCluPaste)
    S_clu2 = get_cluInfo_(S0, S0.iCluPaste);
    vcLabel = sprintf('Unit %d "%s" vs. Unit %d "%s"', ...
        S0.iCluCopy, csNote_clu{S0.iCluCopy}, ...
        S0.iCluPaste, csNote_clu{S0.iCluPaste});
    set(mh_info, 'Label', vcLabel);
else
    S_clu2 = [];
    vcLabel = sprintf('Unit %d "%s"', S0.iCluCopy, csNote_clu{S0.iCluCopy});
    set(mh_info, 'Label', vcLabel);
end
plot_FigPos_(S0, S_clu1, S_clu2);
end %func


%--------------------------------------------------------------------------
function xyPos = clu_pos_(S0, iClu1)
% compute cluster centroid based on the subsampled spikes

trFet1 = get_fet_clu_(S0.S_manual, iClu1);
iSite1 = S0.S_manual.viSite_clu(iClu1);
nSpk1 = size(trFet1,3);
mrPos1 = calc_pos_spk_(trFet1, repmat(iSite1, [nSpk1,1]), S0.P);
xyPos = mean(mrPos1,1);
end %func


%--------------------------------------------------------------------------
function S_cluInfo = get_cluInfo_(S0, iClu)

% determine cluster position
if isempty(iClu), S_cluInfo=[]; return; end
[P, S_manual] = get_(S0, 'P', 'S_manual');

iSite1 = S_manual.viSite_clu(iClu);
viSite = P.miSites(:, iSite1);

xyPos = clu_pos_(S0, iClu);
vcPos = sprintf('Unit %d (x,y):(%0.1f, %0.1f)[pix]', iClu, xyPos/P.um_per_pix);

trFet1 = S0.S_manual.ctrFet_sub_clu{iClu};
nSpk_show = get_set_(P, 'nSpk_show', 30);
viSpk_show = randomSelect_(1:size(trFet1,3), nSpk_show);
trWav = pc2wav_(S0.mrPv_global, trFet1(:,:,viSpk_show));
mrWav_clu = mean(trWav, 3);

% if P.fWav_raw_show
%     trWav = fft_lowpass_(trWav, get_set_(P, 'fc_spkwav_show', []), P.sRateHz);
% end
S_cluInfo = makeStruct_(xyPos, iClu, mrWav_clu, viSite, vcPos, trWav, P);
try
    S_cluInfo.l_ratio = S_manual.vrLRatio_clu(iClu);
    S_cluInfo.isi_ratio = S_manual.vrIsiRatio_clu(iClu);
    S_cluInfo.iso_dist = S_manual.vrIsoDist_clu(iClu);
    S_cluInfo.snr = S_manual.vrSnr_clu(iClu);
    S_cluInfo.uVmin = S_manual.vrVmin_uv_clu(iClu);
    S_cluInfo.uVpp = S_manual.vrVpp_uv_clu(iClu);
catch    
end
end %func


%--------------------------------------------------------------------------
function plot_FigPos_(S0, S_clu1, S_clu2)
[hFig, S_fig] = get_fig_cache_('FigPos');
S_manual = S0.S_manual;
if nargin<2, S_clu2 = []; end

% plot waveform in space
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
else
    cla(S_fig.hAx); hold(S_fig.hAx, 'on');
end
fPlot_spk = isempty(S_clu2);
plot_unit_(S_clu1, S_fig.hAx, [0 0 0], fPlot_spk);
vrPosXY1 = S_clu1.xyPos;
nSpk1 = S_manual.vnSpk_clu(S_clu1.iClu);
try
    if isempty(S_clu2)        
        vcTitle = sprintf('Unit %d: #spikes:%d; x:%0.1fum; y:%0.1fum', S_clu1.iClu, nSpk1, vrPosXY1);
        try
            vcTitle = sprintf('%s\nSNR:%0.1f; %0.1fuVmin; %0.1fuVpp\nIsoD:%0.1f; ISIr:%0.2f; Lrat:%0.2f', ...
                vcTitle, S_clu1.snr, S_clu1.uVmin, S_clu1.uVpp, S_clu1.iso_dist, S_clu1.isi_ratio, S_clu1.l_ratio);
        catch
        end
    else
        nSpk2 = S_manual.vnSpk_clu(S_clu2.iClu);
        vrPosXY2 = S_clu2.xyPos;
        plot_unit_(S_clu2, S_fig.hAx, [1 0 0], fPlot_spk);
        vcTitle = sprintf('Units %d/%d (black/red); (%d/%d) spikes\nSNR=%0.1f/%0.1f; (X,Y)=(%0.1f/%0.1f, %0.1f/%0.1f)um', ...
            S_clu1.iClu, S_clu2.iClu, nSpk1, nSpk2, S_clu1.snr, S_clu2.snr, ...
            vrPosXY1(1), vrPosXY2(1), vrPosXY1(2), vrPosXY2(2));
    end
    title_(S_fig.hAx, vcTitle);
catch
    ;
end
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function [hSpk, hSpkAll] = plot_unit_(S_clu1, hAx, vcColor0, fPlot_spk)
if isempty(S_clu1), return; end
if nargin<2, hAx = axes_new_('FigWav'); end
if nargin<3, vcColor0 = [0 0 0]; end
if nargin<4, fPlot_spk = true; end

P = S_clu1.P;
[~, S_figWav] = get_fig_cache_('FigWav');
maxAmp = S_figWav.maxAmp;
% plot individual unit
nSamples = size(S_clu1.mrWav_clu,1);
vrX = (1:nSamples)'/nSamples;
vrX([1,end])=nan; % line break

if ~fPlot_spk %~isequal(vcColor0, [0 0 0])
    trWav1 = zeros(1,1,0);
else
    trWav1 = S_clu1.trWav;
end

for iWav = size(trWav1,3):-1:0
    if iWav==0
        mrY1 = S_clu1.mrWav_clu / maxAmp;
        lineWidth=1.5;
        vcColor = vcColor0;
    else
        mrY1 = trWav1(:,:,iWav) / maxAmp;
        lineWidth=.5;
        vcColor = .5*[1,1,1];
    end
    vrX1_site = P.mrSiteXY(S_clu1.viSite, 1) / P.um_per_pix;
    vrY1_site = P.mrSiteXY(S_clu1.viSite, 2) / P.um_per_pix;
    mrY1 = bsxfun(@plus, mrY1, vrY1_site');
    mrX1 = bsxfun(@plus, repmat(vrX, [1, size(mrY1, 2)]), vrX1_site');
    line(mrX1(:), mrY1(:), 'Color', vcColor, 'Parent', hAx, 'LineWidth', lineWidth);
end

xlabel(hAx, 'X pos [pix]');
ylabel(hAx, 'Z pos [pix]');
grid(hAx, 'on');
xlim_(hAx, [min(mrX1(:)), max(mrX1(:))]);
ylim_(hAx, [floor(min(mrY1(:))-1), ceil(max(mrY1(:))+1)]);
end %func


%--------------------------------------------------------------------------
function rescale_FigWav_(event, S0, P)
% set(0, 'UserData', S0);

[S_fig, maxAmp_prev, hFigWav] = set_fig_maxAmp_('FigWav', event);                
set_fig_(hFigWav, plot_tnWav_clu_(S_fig, P, S0));
multiplot(S0.hCopy, S_fig.maxAmp);
if ~isempty(S0.iCluPaste)
    multiplot(S0.hPaste, S_fig.maxAmp);
end
rescale_spikes_(S_fig.hSpkAll, maxAmp_prev, P);
title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp)); %update scale
end


%--------------------------------------------------------------------------
function plot_FigTime_(S0)
% plot FigTime window. Uses subsampled data

S_manual = S0.S_manual; 
P = S0.P; 
[hFig, S_fig] = get_fig_cache_('FigTime');

%----------------
% collect info
iSite = S_manual.viSite_clu(S0.iCluCopy);
[vrFet0, vrTime0] = getFet_site_(iSite, [], S0);    % plot background    
[vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, S0.iCluCopy, S0); % plot iCluCopy

vcTitle = '[H]elp; (Sft)[Left/Right]:Sites/Features; (Sft)[Up/Down]:Scale; [B]ackground; [S]plit; [R]eset view; [P]roject; [M]erge; (sft)[Z] pos; [E]xport selected; [C]hannel PCA';
if ~isempty(S0.iCluPaste)
    [vrFet2, vrTime2] = getFet_site_(iSite, S0.iCluPaste, S0);
    vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', S0.iCluCopy, S0.iCluPaste, vcTitle);
else
    vrFet2 = [];
    vrTime2 = [];
    vcTitle = sprintf('Clu%d (black); %s', S0.iCluCopy, vcTitle);
end
time_lim = double([0, abs(S0.viTime_spk(end))] / P.sRateHz);

%------------
% draw
if isempty(S_fig)
    S_fig.maxAmp = P.maxAmp;
    S_fig.hAx = axes_new_(hFig);
    set(S_fig.hAx, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');
    grid(S_fig.hAx, 'on');
    
    % first time
    S_fig.hPlot0 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(1,:), 'MarkerSize', 3, 'LineStyle', 'none');
    S_fig.hPlot1 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(2,:), 'MarkerSize', 5, 'LineStyle', 'none');
    S_fig.hPlot2 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(3,:), 'MarkerSize', 5, 'LineStyle', 'none');   %place holder  
    xlabel(S_fig.hAx, 'Time (s)');         
    
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigTime_);
    S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlot0, S_fig);
    if ~isempty(P.time_tick_show) %tick mark
        set(S_fig.hAx, 'XTick', time_lim(1):P.time_tick_show:time_lim(end));
    end
end
vpp_lim = [0, abs(S_fig.maxAmp)];
% iFet = S_fig.iFet;
% iFet = 1;
if ~isfield(S_fig, 'iSite'), S_fig.iSite = []; end
update_plot_(S_fig.hPlot0, vrTime0, vrFet0);
update_plot_(S_fig.hPlot1, vrTime1, vrFet1);
update_plot_(S_fig.hPlot2, vrTime2, vrFet2);
mouse_figure(hFig, S_fig.hAx); % allow zoom using wheel
% button click function to select individual spikes, all spikes plotted

if isfield(S_fig, 'vhAx_track')
    toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 0);
end

if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
toggleVisible_(S_fig.hPlot0, S_fig.fPlot0);

axis_(S_fig.hAx, [time_lim, vpp_lim]);
title_(S_fig.hAx, vcTitle);    
ylabel(S_fig.hAx, vcYlabel);

S_fig = struct_merge_(S_fig, makeStruct_(iSite, time_lim, P, vpp_lim, viSpk1));
S_fig.csHelp = {...
    'Up/Down: change channel', ...
    'Left/Right: Change sites', ...
    'Shift + Left/Right: Show different features', ...
    'r: reset scale', ...
    'a: auto-scale', ...
    'c: show pca across sites', ...
    'e: export cluster info', ...
    'f: export cluster feature', ...
    'Zoom: mouse wheel', ...
    'H-Zoom: press x and wheel. space to reset', ...
    'V-Zoom: press y and wheel. space to reset', ...
    'Drag while pressing wheel: pan'};
        
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function [vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, iClu, S0)
% just specify iSite to obtain background info
% 2016 07 07 JJJ
% return feature correspojnding to a site and cluster
% requiring subsampled info: cvrVpp_site and cmrFet_site. store in S0

if nargin < 2, iClu = []; end
P = S0.P;


if ~isfield(P, 'vcFet_show'), P.vcFet_show = 'vpp'; end
if isempty(iClu)
    viSite = S0.P.miSites(1:P.nSites_fet,iSite);
    [vrFet1, viSpk1] = get_fet_clu_(S0.S_manual, iClu, viSite);
    vrFet1 = abs(squeeze(vrFet1(1,1,:)));
else
    [vrFet1, viSpk1] = get_fet_clu_(S0.S_manual, iClu, iSite);
    vrFet1 = abs(squeeze(vrFet1(1,1,:)));
end
vrTime1 = double(S0.viTime_spk(viSpk1)) / P.sRateHz;
vcYlabel = sprintf('Site %d (PC)', iSite);
% % label
% switch lower(P.vcFet_show)
%     case {'vpp', 'vmin'} %voltage feature
%         vcYlabel = sprintf('Site %d (\\mu%s)', iSite, P.vcFet_show);
%     otherwise %other feature options
%         vcYlabel = sprintf('Site %d (%s)', iSite, P.vcFet_show);
% end

end %func 


%--------------------------------------------------------------------------
% 122917 JJJ: modified
function S_fig = rescale_FigProj_(event, hFig, S_fig, S0)
% S_fig = rescale_FigProj_(event, hFig, S_fig, S0)
% S_fig = rescale_FigProj_(maxAmp)

if nargin<2, hFig = []; end
if nargin<3, S_fig = []; end
if nargin<4, S0 = []; end
if isempty(hFig) || isempty(S_fig), [hFig, S_fig] = get_fig_cache_('FigProj'); end
if isempty(S0), S0 = get(0, 'UserData'); end
P = S0.P;

if isnumeric(event)
    S_fig.maxAmp = event;
else
    S_fig.maxAmp = change_amp_(event, S_fig.maxAmp);     
end
vhPlot = [S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2];
if isempty(S0.iCluPaste), vhPlot(end) = []; end
rescaleProj_(vhPlot, S_fig.maxAmp, S0.P);
switch lower(P.vcFet_show)
    case {'vpp', 'vmin', 'vmax'}
        S_fig.vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        S_fig.vcYLabel = 'Site # (%0.0f \\muV_{min})';
    otherwise
        S_fig.vcXLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', P.vcFet_show, P.vcFet_show, P.vcFet_show);
        S_fig.vcYLabel = sprintf('Site # (%%0.0f %s)', P.vcFet_show);    
end
xlabel(S_fig.hAx, sprintf(S_fig.vcXLabel, S_fig.maxAmp));   
ylabel(S_fig.hAx, sprintf(S_fig.vcYLabel, S_fig.maxAmp));  
if nargout==0
    set(hFig, 'UserData', S_fig);
end
end


%--------------------------------------------------------------------------
function rescale_FigTime_(event, S0)
% rescale_FigTime_(event, S0, P)
% rescale_FigTime_(maxAmp, S0, P)
if nargin<2, S0 = []; end

if isempty(S0), S0 = get0_(); end

[S_fig, maxAmp_prev] = set_fig_maxAmp_('FigTime', event);
ylim_(S_fig.hAx, [0, 1] * S_fig.maxAmp);
% imrect_set_(S_fig.hRect, [], [0, S_fig.maxAmp]);
iSite = S0.S_manual.viSite_clu(S0.iCluCopy);
end %func


%--------------------------------------------------------------------------
function [S_fig, maxAmp_prev, hFig] = set_fig_maxAmp_(vcFig, event)
[hFig, S_fig] = get_fig_cache_(vcFig);
if isempty(S_fig)
    P = get0_('P');
    S_fig.maxAmp = P.maxAmp;
end
maxAmp_prev = S_fig.maxAmp;
if isnumeric(event)
    S_fig.maxAmp = event;
else
    S_fig.maxAmp = change_amp_(event, maxAmp_prev);
end
set(hFig, 'UserData', S_fig);
end


%--------------------------------------------------------------------------
function rescaleProj_(vhPlot1, maxAmp, P)
if nargin<3, P = get0_('P'); end
for iPlot=1:numel(vhPlot1)
    hPlot1 = vhPlot1(iPlot);
    S_plot1 = get(hPlot1, 'UserData');
    update_plot2_proj_();
    S_plot1 = struct_delete_(S_plot1, 'hPoly'); %, 'hPlot_split'
    [vrX, vrY, viPlot, tr_dim] = amp2proj_(S_plot1.mrMin, S_plot1.mrMax, maxAmp, P.maxSite_show, P);
    S_plot1 = struct_add_(S_plot1, viPlot, vrX, vrY, maxAmp);
    set(hPlot1, 'XData', vrX, 'YData', vrY, 'UserData', S_plot1);
end
end %func


%--------------------------------------------------------------------------
function S = struct_delete_(S, varargin)
% delete and set to empty

for i=1:numel(varargin)
    try 
        delete(S.(varargin{i}));
        S.(varargin{i}) = [];
    catch
        ;
    end
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function varargout = get0_(varargin)
% returns get(0, 'UserData') to the workspace
% [S0, P] = get0_();
% [S0, P, S_clu] = get0_();
% [var1, var2] = get0_('var1', 'var2'); %sets [] if doesn't exist
S0 = get(0, 'UserData'); 
if ~isfield(S0, 'S_clu'), S0.S_clu = []; end
if nargin==0
    varargout{1} = S0; 
    if nargout==0, assignWorkspace_(S0); return; end
    if nargout>=1, varargout{1} = S0; end
    if nargout>=2, varargout{2} = S0.P; end
    if nargout>=3, varargout{3} = S0.S_clu; end
    return;
end
for i=1:nargin
    try                
        eval(sprintf('%s = S0.%s;', varargin{i}, varargin{i}));
        varargout{i} = S0.(varargin{i});
    catch
        varargout{i} = [];
    end
end
end %func


%--------------------------------------------------------------------------
function [maxAmp, mrAmp_prev] = change_amp_(event, maxAmp, varargin)
% varargin: plot object to rescale
% Change amplitude scaling 
% change_amp_(event, maxAmp, varargin)
% change_amp_(event) % directly set
% if nargin<3, hPlot=[]; end
factor = sqrt(2);
if key_modifier_(event, 'shift'), factor = factor ^ 4; end
if isempty(maxAmp)
    vhLine = varargin(:);
    vhLine = vhLine(~isempty_(vhLine));
    maxAmp = get_userdata_(vhLine{1}, 'scale');
end
mrAmp_prev = maxAmp;
if strcmpi(event.Key, 'uparrow')
    maxAmp = maxAmp / factor;
elseif strcmpi(event.Key, 'downarrow')
    maxAmp = maxAmp * factor;
end
for iPlot = 1:numel(varargin)
    try
        multiplot(varargin{iPlot}, maxAmp);
    catch            
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


%--------------------------------------------------------------------------
function plot_FigCorr_(S0)
% hFigCorr plot
jitter_ms = .5; % bin size for correlation plot
nLags_ms = 25; %show 25 msec

if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P; 
S_manual = S0.S_manual;
P.jitter_ms = jitter_ms;
P.nLags_ms = nLags_ms;

[hFig, S_fig] = get_fig_cache_('FigCorr'); 
iClu1 = get_set_(S0, 'iCluCopy', 1);
iClu2 = get_(S0, 'iCluPaste');
if isempty(iClu2), iClu2 = iClu1; end

jitter = round(P.sRateHz / 1000 * P.jitter_ms); %0.5 ms
nLags = round(P.nLags_ms / P.jitter_ms);

vi1 = int32(double(clu_time_(S0, iClu1)) /jitter);

if iClu1~=iClu2
    vi1 = [vi1, vi1-1, vi1+1]; %allow missing one
end
vi2 = int32(double(clu_time_(S0, iClu2)) /jitter);
viLag = -nLags:nLags;
vnCnt = zeros(size(viLag));
for iLag=1:numel(viLag)
    if iClu1 == iClu2 && viLag(iLag)==0, continue; end
    vnCnt(iLag) = numel(intersect(vi1, vi2+viLag(iLag)));
end
vrTime_lag = viLag * P.jitter_ms;

%--------------
% draw
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    S_fig.hBar = bar_(vrTime_lag, vnCnt, 1);     
    xlabel('Time (ms)'); 
    ylabel('Counts');
    grid on; 
    set(S_fig.hAx, 'YScale', 'log');
else
    set(S_fig.hBar, 'XData', vrTime_lag, 'YData', vnCnt);
end
if iClu1==iClu2
    title_(S_fig.hAx, sprintf('Clu%d', iClu1));
else
    title_(S_fig.hAx, sprintf('Clu%d vs Clu%d', iClu1, iClu2));
end
xlim_(S_fig.hAx, [-nLags, nLags] * P.jitter_ms);
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function hPlot = bar_(varargin)
try hPlot = bar(varargin{:}); catch, hPlot = []; end
end %func


%--------------------------------------------------------------------------
function hPlot = stem_(varargin)
try hPlot = stem(varargin{:}); catch, hPlot = []; end
end %func


%--------------------------------------------------------------------------
function S_txt = loadjson_(vcArg_txt)
S_txt=[]; 
if ~exist_file_(vcArg_txt), return; end
addpath_('jsonlab-1.5/');
S_txt = loadjson(vcArg_txt);
end %func


%--------------------------------------------------------------------------
function vl = isempty_(cvr)
if iscell(cvr)
    vl = cellfun(@isempty, cvr);
else
    vl = isempty(cvr);
end
end %func


%--------------------------------------------------------------------------
function keyPressFcn_FigTime_(hObject, event, S0)

if nargin<3, S0 = get(0, 'UserData'); end
[P, S_manual, hFig] = deal(S0.P, S0.S_manual, hObject);
S_fig = get(hFig, 'UserData');

nSites = numel(P.viSite2Chan);
switch lower(event.Key)
    case {'leftarrow', 'rightarrow'}
        if ~isVisible_(S_fig.hAx)
            msgbox_('Channel switching is disabled in the position view'); return; 
        end
        factor = key_modifier_(event, 'shift')*3 + 1;
        if strcmpi(event.Key, 'rightarrow')
            S_fig.iSite = min(S_fig.iSite + factor, nSites);
        else
            S_fig.iSite = max(S_fig.iSite - factor, 1);
        end
        set(hFig, 'UserData', S_fig);        
        update_FigTime_(S0);                

    case {'uparrow', 'downarrow'} %change ampl
        if ~isVisible_(S_fig.hAx)
            msgbox_('Zoom is disabled in the position view'); return; 
        end
        rescale_FigTime_(event, S0, P);
        
    case 'r' %reset view
        if ~isVisible_(S_fig.hAx), return; end
        axis_(S_fig.hAx, [S_fig.time_lim, S_fig.vpp_lim]);
        imrect_set_(S_fig.hRect, S_fig.time_lim, S_fig.vpp_lim);

    case 'm' %merge
        ui_merge_(S0);
        
    case 'h' %help
        msgbox_(S_fig.csHelp, 1);
        
    case 'b' %background spike toggle
        if isVisible_(S_fig.hAx)  
            S_fig.fPlot0 = toggleVisible_(S_fig.hPlot0);
        else
            S_fig.fPlot0 = toggleVisible_(S_fig.hPlot0_track);
        end
        set(hFig, 'UserData', S_fig);
        
    case 't'
        plot_FigTime_(S0);
        
    case 's' %split. draw a polygon
        if ~isempty(S0.iCluPaste)
            msgbox_('Select one cluster'); return;
        end
        try
            hPoly = impoly_();
            if isempty(hPoly); return; end
            mrPolyPos = getPosition(hPoly);
            vrX1 = double(get(S_fig.hPlot1, 'XData'));
            vrY1 = double(get(S_fig.hPlot1, 'YData'));
            vlIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
            hSplit = line(vrX1(vlIn), vrY1(vlIn), 'Color', [1 0 0], 'Marker', '.', 'LineStyle', 'none');
            if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes'), 'yes')
                split_clu_(S0.iCluCopy, vlIn);
            end
            delete_multi_(hPoly, hSplit);
        catch
            disp(lasterror());
        end
        
    case 'p' %update projection view
        vrPos = getPosition(S_fig.hRect);
        tlim_proj = [vrPos(1), sum(vrPos([1,3]))];
        P.tlim_proj = tlim_proj;
        plot_FigProj_(S0);        
end %switch
end %func


%--------------------------------------------------------------------------
function update_FigTime_(S0)
% display features in a new site

[hFig, S_fig] = get_fig_cache_('FigTime');
if ~isVisible_(S_fig.hAx), return; end
[vrFet0, vrTime0, vcYlabel] = getFet_site_(S_fig.iSite, [], S0);
if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
toggleVisible_(S_fig.hPlot0, S_fig.fPlot0);
update_plot_(S_fig.hPlot0, vrTime0, vrFet0);
set(S_fig.hPlot1, 'YData', getFet_site_(S_fig.iSite, S0.iCluCopy, S0));
if ~isempty(S0.iCluPaste)
    set(S_fig.hPlot2, 'YData', getFet_site_(S_fig.iSite, S0.iCluPaste, S0));
else
    hide_plot_(S_fig.hPlot2);
end
ylim_(S_fig.hAx, [0, 1] * S_fig.maxAmp);
% imrect_set_(S_fig.hRect, [], [0, 1] * S_fig.maxAmp);    
grid(S_fig.hAx, 'on');
ylabel(S_fig.hAx, vcYlabel);
end %func


%--------------------------------------------------------------------------
function flag = isVisible_(hObj)
flag = strcmpi(get(hObj, 'Visible'), 'on');
end %func


%--------------------------------------------------------------------------
function update_spikes_(S0)

hMsg = msgbox_open_('Updating spikes');
fig_prev = gcf;
hFig_wait = figure_wait_(1);
[~, S_fig] = get_fig_cache_('FigWav');
plot_spkwav_(S_fig, S0);
close_(hMsg);
figure_wait_(0, hFig_wait);
figure(fig_prev);
end %func


%--------------------------------------------------------------------------
function fExit = save_manual_(varargin)
% just save S_manual

S0 = get0_();
P = S0.P;
vcFile_manual = subsFileExt_(P.vcFile_prm, '_manual_irc.mat');
fExit = 1;
switch lower(questdlg_(['Save to ', vcFile_manual, ' ?'], 'Confirmation', 'Yes'))
    case 'yes'
        hMsg = msgbox_('Saving... (this closes automatically)');
        struct_save_(S0.S_manual, vcFile_manual); % 1 will skip figure saving
        fExit = 1;
        close_(hMsg);
    case 'no'
        fExit = 1;
    case 'cancel' 
        fExit = 0;
        return;
end
end %func;


%--------------------------------------------------------------------------
function button_FigWavCor_(xyPos, vcButton)
S0 = get(0, 'UserData');
xyPos = round(xyPos);
switch lower(vcButton)
    case 'normal' %left click        
        S0.iCluCopy = xyPos(1);
        if diff(xyPos) == 0
            S0.iCluPaste = [];
        else
            S0.iCluPaste = xyPos(2);
        end
        S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
        S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom        
end %switch
end %func


%--------------------------------------------------------------------------
function keyPressFcn_FigWavCor_(hObject, event)
S0 = get(0, 'UserData');
switch lower(event.Key)
    case 'm' %merge
        ui_merge_(S0);
    case 's' %split
        auto_split_(1); %multi
    case {'d', 'backspace', 'delete'} %delete
        ui_delete_(S0);        
end %switch
end %func


%--------------------------------------------------------------------------
function help_FigWav_(hObject, event)
[~, S_fig] = get_fig_cache_('FigWav');
msgbox_(S_fig.csHelp, 1);
end %func


%--------------------------------------------------------------------------
function select_shank_(hObject, event)
% select active shank
iShank = str2num(strrep(hObject.Text, 'Shank ', ''));
disp(hObject.Text);
cell_menu_shank = get_userdata_(hObject.Parent, 'cell_menu_shank');
for iShank_ = 1:numel(cell_menu_shank)
    set(cell_menu_shank{iShank_}, 'Checked', ifeq_(iShank_==iShank, 'on', 'off'));
end
set0_(iShank);
end %func


%--------------------------------------------------------------------------
% irc2.m
function varargout = load0_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = is_sorted_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_fet_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_set_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = xylabel_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = cell2mat_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = version_matlab_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = load_fet_site_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = validate_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = pc2wav_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = wavcor_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = describe_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = calc_pos_spk_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_save_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = subsFileExt_(varargin), cell_out = call_irc2_(dbstack(), varargin, nargout); varargout = cell_out; end

%--------------------------------------------------------------------------
% irc.m

function varargout = msgbox_open_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_raster_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_split_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = select_polygon_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = calc_cov_spk_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = pca_pc_spk_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_spkwav_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = rescale_spikes_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = imrect_set_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = assignWorkspace_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = addpath_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = ui_show_elective_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end

function varargout = export_csv_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = export_klusters_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = export_mrWav_clu_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = export_phy_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = export_quality_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = export_rate_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = export_tmrWav_clu_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = export_tnWav_spk_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = help_FigWav_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = issue_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = merge_auto_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = open_prm_folder_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_drift_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = proj_view_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = reload_prm_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = save_figures_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = traces_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = ui_show_all_chan_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = ui_show_drift_view_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = unit_annotate_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = wiki_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = delete_multi_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_fig_all_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end

% function varargout = clu_info_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = figure_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = key_modifier_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = plot_FigCorr_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_psth_trial_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = reset_position_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = toggleVisible_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = ui_delete_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = ui_merge_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = update_FigCor_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = update_spikes_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = about_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = clipboard_prm_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = delete_auto_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = edit_prm_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = keyPressFcn_FigWavCor_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = axes_new_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_update_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_group_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = impoly_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = split_clu_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = tnWav2uV_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = wave_split_manual_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = auto_split_wav_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = tr2plot_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = disperr_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = update_menu_trials_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = create_figure_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = close_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = save_log_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = text3_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = text_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = cellstr_subs_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_probe_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = S_clu_position_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = S_clu_quality_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = exist_file_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = questdlg_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_fig_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = struct_merge_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = delete_clu_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = msgbox_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = correlogram_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
% function varargout = keyPressFcn_cell_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = S_clu_prune_icl_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_fig_cache_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end