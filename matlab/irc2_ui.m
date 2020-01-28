%--------------------------------------------------------------------------
% 2020/1/9 

function irc2_ui(P)
% global fDebug_ui

if ~is_sorted_(P)
    fprintf(2, 'File must to be sorted first (run "irc spikesort %s")\n', P.vcFile_prm); 
    return; 
end
S0 = load0_(P.vcFile_prm);
S0.trPc_spk = load_fet_(S0, P, 1); % 
S0.S_ui = load_(strrep(P.vcFile_prm, '.prm', '_ui_irc.mat'));
if ~isempty(S0.S_ui)
    % ask if you want to restart
    if strcmpi(questdlg_('Load last saved?', 'Confirmation'), 'no')
        S0.S_ui=[];
    end
end
if isempty(S0.S_ui)
    S0.S_ui = create_S_ui_(S0); 
end

clear mouse_figure;
clear get_fig_cache_ get_tag_ %clear persistent figure handles
close_(get_fig_('FigTrial')); %close previous FigTrial figure
close_(get_fig_('FigTrial_b')); %close previous FigTrial figure
S0.S_ui = struct_merge_(S0.S_ui, ...
    struct('iCluCopy', 1, 'iCluPaste', [], 'hCopy', [], 'hPaste', [], 'nSites', numel(P.viSite2Chan)));
hMsg = msgbox_('Plotting... (this closes automatically)'); t1=tic;

S0.S_ui = struct_merge_(S0.S_ui, figures_manual_(P)); %create figures for manual interface
plot_FigRD_(S0); % ask user before doing so
plot_FigWav_(S0); % hFigWav %do this after for ordering
plot_FigWavCor_(S0);  % need `S_clu.mrWavCor`
button_CluWav_simulate_(1, [], S0); %select first clu
auto_scale_proj_time_(S0);
S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
% S_log = load_(strrep(P.vcFile_prm, '.prm', '_log.mat'), [], 0);
% if ~isempty(S_log), S0.cS_log = {S_log}; end
% save_log_('start', S0); %crash proof log

% Finish up
close_(hMsg);
fprintf('UI creation took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function plot_FigWav_(S0)

[P, S_auto] = get_(S0, 'P', 'S_auto');
[hFig, S_fig] = get_fig_cache_('FigWav'); 

% Show number of spikes per clusters
P.LineWidth = 1; %plot a thicker line
P.viSite_clu = S_auto.viSite_clu;
nSites = numel(P.viSite2Chan);
if isempty(S_fig) % initialize
    S_fig.maxAmp = P.maxAmp;
    S_fig.hAx = axes(hFig);
    set(S_fig.hAx, 'Position', [.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual'); 
    grid(S_fig.hAx, 'on');    
    S_fig.vcTitle = '%0.1f uV; [H]elp; (Sft)[Up/Down]; (Sft)[Left/Right]; [M]erge; [S]plit; [D]elete; [A]:Resample; [P]STH; [Z]oom; [Space]:Match';
    xylabel_(S_fig.hAx, 'Cluster #', 'Site #', sprintf(S_fig.vcTitle, S_fig.maxAmp));

    set(hFig, 'KeyPressFcn', @keyPressFcn_FigWav_, 'CloseRequestFcn', @exit_manual_, 'BusyAction', 'cancel');
    axis(S_fig.hAx, [0, S_auto.nClu + 1, 0, nSites + 1]);
    add_menu_(hFig, P);      
%     vrPos_ = get(hFig, 'OuterPosition');
    mouse_figure(hFig, S_fig.hAx, @button_CluWav_);
    S_fig = plot_spkwav_(S_fig, S0); %plot spikes
    S_fig = plot_tnWav_clu_(S_fig, P, S0); %do this after plotSpk_
    S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hSpkAll, S_fig);
%     set(hFig, 'OuterPosition', vrPos_);
else
%     mh_info = [];
    S_fig = plot_spkwav_(S_fig, S0); %plot spikes
    try delete(S_fig.vhPlot); catch; end %delete old text
    S_fig = rmfield_(S_fig, 'vhPlot');
    S_fig = plot_tnWav_clu_(S_fig, P, S0); %do this after plotSpk_
%     xylabel_(S_fig.hAx, 'Cluster #', 'Site #');
end

% create text
% S0 = set0_(mh_info);
fText = get_set_(S_fig, 'fText', get_set_(P, 'Text', 1));
S_fig = figWav_clu_count_(S_fig, S_auto, fText);
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
function S_fig = figWav_clu_count_(S_fig, S_clu, fText)
if nargin==0, [hFig, S_fig] = get_fig_cache_('FigWav'); end
if nargin<3, fText = 1; end

if fText
    csText_clu = arrayfun(@(i)sprintf('%d(%d)', i, S_clu.vnSpk_clu(i)), 1:S_clu.nClu, 'UniformOutput', 0);
else
    csText_clu = arrayfun(@(i)sprintf('%d', i), 1:S_clu.nClu, 'UniformOutput', 0);
end
set(S_fig.hAx, 'Xtick', 1:S_clu.nClu, 'XTickLabel', csText_clu, 'FontSize', 8);
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
[P, viSite_spk, S_auto] = deal(S0.P, S0.viSite_spk, S0.S_auto);

[cvrX, cvrY, cviSite] = deal(cell(S_auto.nClu, 1));
vnSpk = zeros(S_auto.nClu, 1);
miSites_clu = P.miSites(:, S_auto.viSite_clu);
if isfield(S_fig, 'maxAmp')
    maxAmp = S_fig.maxAmp;
else
    maxAmp = P.maxAmp;
end

for iClu = 1:S_auto.nClu        
    try   
        viSpk_show = randomSelect_(S_clu_viSpk_(S_auto, iClu, viSite_spk), P.nSpk_show);
        trWav1 = pc2wav_(S0.mrPv_global, S0.trPc_spk(:,:,viSpk_show));
        viSite_show = miSites_clu(:, iClu);
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
[P, S_clu] = deal(S0.P, S0.S_clu);

if ~isempty(S0.iCluPaste), msgbox_('Select one cluster', 1); return; end
if S_clu.vnSpk_clu(S0.iCluCopy)<3
    msgbox_('At least three spikes required for splitting', 1); return; 
end
    
hMsg = msgbox_('Splitting... (this closes automatically)');
iClu1 = S0.iCluCopy;
iSite1 = S_clu.viSite_clu(iClu1);
if fMulti
    viSites1 = P.miSites(1:P.nSites_fet, iSite1);
else
    viSites1 = iSite1;
end
% mrSpkWav1 = tnWav2uV_(tnWav_sites_(tnWav_spk, S_clu.cviSpk_clu{iClu1}, viSites1));
trSpkWav1 = tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu1}, viSites1, S0), P, 0);
% mrSpkWav1 = tnWav2uV_(tnWav_spk_sites_(find(S_clu.viClu==iClu1), viSites1, S0), P);
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
[P, S_auto] = get_(S0, 'P', 'S_clu');
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
            S0.iCluCopy = S_auto.nClu;
        elseif ~key_modifier_(event, 'shift');
            if strcmpi(event.Key, 'leftarrow')
                if S0.iCluCopy == 1, return; end
                S0.iCluCopy = S0.iCluCopy - 1;
            else
                if S0.iCluCopy == S_auto.nClu, return; end
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
                if S0.iCluPaste == S_auto.nClu, return; end
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
        mrWavCor = S_auto.mrWavCor;
        mrWavCor(S0.iCluCopy,S0.iCluCopy) = -inf;
        [~,S0.iCluPaste] = max(mrWavCor(:,S0.iCluCopy));
        set(0, 'UserData', S0);
        button_CluWav_simulate_([], S0.iCluPaste);        
    case 's', auto_split_(1, S0);        
    case 'r' %reset view
        hFig_wait = figure_wait_(1);
        axis_([0, S0.S_clu.nClu + 1, 0, numel(P.viSite2Chan) + 1]);
        figure_wait_(0, hFig_wait);        
    case {'d', 'backspace', 'delete'}, S0 = ui_delete_(S0);        
    case 'z' %zoom
        ui_zoom_(S0, hFig);
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
        figWav_clu_count_(S_fig, S_auto, ~fText);          
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

mh_history = uimenu_(hFig, 'Label', 'History', 'Tag', 'mh_history'); 

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
    P = S0.P;
    fExit = save_manual_(P);
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
function S_ui = create_S_ui_(S0)
% P.fParfor = 0; % disable parfor for GUI
% fDebug_ui = false;

% keep trPc_spk loaded for manual

% compute cluster average waveform and waveform similarity matrix
% similar auto-merging
fprintf(2, 'create_S_ui_: todo\n');
S_ui = [];
% if isempty(get_(S0.S_clu, 'trWav_spk_clu'))
%     S0.S_clu = calc_clu_wav_(S0, P);    
%     S0.S_clu = S_clu_quality_(S0.S_clu, P);
% end
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
function [S_clu, viClu_delete] = calc_clu_wav_(S0, P)

S_clu = S0.S_clu;
S_clu.nClu = max(S_clu.viClu);
cviSpk_clu = arrayfun_(@(x)find(S_clu.viClu==x), (1:S_clu.nClu)');
viSite_clu = cellfun(@(x)mode(S0.viSite_spk(x)), cviSpk_clu);
nT = size(S0.mrPv_global,1);
if isempty(get_(S0, 'trPc_spk'))
    trPc_spk = load_fet_(S0, P, 1);
else
    trPc_spk = S0.trPc_spk;
end
nSites_spk = size(trPc_spk,2);
nSites = numel(P.viSite2Chan);
trWav_clu = zeros(nT, nSites_spk, S_clu.nClu, 'single');
tmrWav_clu = zeros(nT, nSites, S_clu.nClu, 'single'); %global waveform
for iClu = 1:S_clu.nClu
    viSpk1 = cviSpk_clu{iClu};
    if isempty(viSpk1), continue; end
    viSite1 = S0.viSite_spk(viSpk1);
    iSite1 = viSite_clu(iClu);
    trPc1 = trPc_spk(:,:,viSpk1);
    mrPc1 = mean(trPc1(:,:,viSite1==iSite1),3);
    mrWav1 = S0.mrPv_global * mrPc1;
    trWav_clu(:,:,iClu) = mrWav1;
    viSite_clu1 = P.miSites(:,iSite1);
    tmrWav_clu(:,viSite_clu1,iClu) = mrWav1;
end
[trWav_raw_clu, tmrWav_raw_clu] = deal(trWav_clu, tmrWav_clu); % for now raw is not saved
[trWav_spk_clu, tmrWav_spk_clu] = deal(trWav_clu, tmrWav_clu); % for now raw is not saved
S_clu = struct_add_(S_clu, trWav_clu, trWav_raw_clu, tmrWav_clu, viSite_clu, tmrWav_raw_clu, trWav_spk_clu, tmrWav_spk_clu);

% compute SNR
mrWav1_clu = squeeze_(trWav_clu(:,1,:),2); 
vrVmax_clu = max(mrWav1_clu)';
vrVmin_clu = min(mrWav1_clu)';
vrVpp_clu = vrVmax_clu - vrVmin_clu;
vrRms_site = S0.vrThresh_site(:) / S0.P.qqFactor;
S_clu.vrSnr_clu = abs(vrVmin_clu(:)) ./ vrRms_site(viSite_clu);
S_clu.vrSnr2_clu = abs(vrVpp_clu(:)) ./ vrRms_site(viSite_clu);

viClu_delete = [];

% update similarity
S0.S_clu = S_clu;
mrWavCor = wave_similarity_clu_(S0, P);
S_clu.mrWavCor = mrWavCor + mrWavCor'; % make it symmetric

S_clu = S_clu_position_(S_clu);
S_clu = S_clu_refresh_(S_clu, 1, S0.viSite_spk);
if ~isfield(S_clu, 'csNote_clu'), S_clu.csNote_clu = cell(S_clu.nClu, 1); end
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
create_figure_('FigPos', [0 0 .15 .5], ['Unit position; ', P.vcFile_prm], 1, 0);
create_figure_('FigMap', [0 .5 .15 .5], ['Probe map; ', P.vcFile_prm], 1, 0);  

create_figure_('FigWav', [.15 .25 .35 .75],['Averaged waveform: ', P.vcFile_prm], 0, 1);
create_figure_('FigTime', [.15 0 .7 .25], ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);

create_figure_('FigProj', [.5 .25 .35 .5], ['Feature projection: ', P.vcFile_prm]);
create_figure_('FigWavCor', [.5 .75 .35 .25], ['Waveform correlation (click): ', P.vcFile_prm]);

create_figure_('FigHist', [.85 .75 .15 .25], ['ISI Histogram: ', P.vcFile_prm]);
create_figure_('FigIsi', [.85 .5 .15 .25], ['Return map: ', P.vcFile_prm]);
create_figure_('FigCorr', [.85 .25 .15 .25], ['Time correlation: ', P.vcFile_prm]);
create_figure_('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.vcFile_prm]);

% drawnow_();
csFig = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
cvrFigPos0 = cellfun(@(vc)get(get_fig_(vc), 'OuterPosition'), csFig, 'UniformOutput', 0);
S_fig = cell2struct(cvrFigPos0, csFig, 2);
end %func





%--------------------------------------------------------------------------
% keeps trPc in the main memory, not sent out to workers
function [trWav_clu, mrDist_clu] = mean_clu_(S0)
% use 1000 waveforms to compute average, store 1000 subsamples

P=S0.P;
S_clu = S0.S_clu;
t_fun = tic;
% fprintf('\tAutomated merging based on waveform similarity...\n'); t_template=tic;
[viClu, vrRho] = struct_get_(S_clu, 'viClu', 'rho');
[ccviSpk_site_load, ccviSpk_site2_load, type_fet, dimm_fet, mrPv, vrThresh_site] = ...
    get_(S0, 'ccviSpk_site_load', 'ccviSpk_site2_load', 'type_fet', ...
        'dimm_fet', 'mrPv_global', 'vrThresh_site');
nShift_max = ceil(diff(P.spkLim) * P.frac_shift_merge / 2);
viShift = -nShift_max:nShift_max;
[knn, nSpk_min, vcFile_prm] = get_(P, 'knn', 'knn', 'vcFile_prm');
nClu = S_clu.nClu;

% try loading the entire miKnn to RAM
nSites = size(P.miSites,2);
[viLim_drift, mlDrift] = get_(S_clu.S_drift, 'viLim_drift', 'mlDrift');
% nDrift = size(mlDrift,1);
S_auto = makeStruct_(nClu, nSites, knn, vcFile_prm, nSpk_min, ...
    vrRho, viClu, viLim_drift, ccviSpk_site_load, ccviSpk_site2_load, ...
    type_fet, dimm_fet, mrPv, vrThresh_site, viShift, mlDrift, maxWavCor);

% compute pairwise distance in parallel by sites
[cmlDist_site, cvlExist_site] = deal(cell(nSites, 1));
fParfor = get_set_(P, 'fParfor', 1) && nSites > 1;
if fParfor %&& ~isLargeRecording_(P)
    try
        parfor iSite = 1:nSites
            try
                [cmlDist_site{iSite}, cvlExist_site{iSite}] = wave_drift_site_(iSite, S_auto);
            catch
            end
        end
    catch
    end
end

% merge cluster pairwise distance
mlDist_clu = false(nClu);
mlDist_clu(sub2ind([nClu,nClu], 1:nClu, 1:nClu)) = true; % self join
vlExist_clu = false(1, nClu);
for iSite = 1:nSites
    if isempty(cmlDist_site{iSite})
        [cmlDist_site{iSite}, cvlExist_site{iSite}] = wave_drift_site_(iSite, S_auto);
    end
    mlDist_clu = mlDist_clu | cmlDist_site{iSite};
    vlExist_clu = vlExist_clu | cvlExist_site{iSite};
    [cmlDist_site{iSite}, cvlExist_site{iSite}] = deal([]); % clear memory
end
viClu_remove = find(~vlExist_clu);
fprintf('\n\twave_similarity_: took %0.1fs\n', toc(t_fun));
end %func


%--------------------------------------------------------------------------
function [ctrPc_drift, cviClu_drift, mrDist_clu] = wave_drift_site_(iSite1, S_auto)

NUM_KNN = 10;
fUseSecondSite = 1;

import_struct_(S_auto);
nDrift = size(mlDrift, 1);
% fprintf('\twave_similarity_site_pre_: Site%d... ', iSite1); t_fun=tic;
miKnn1 = load_miKnn_site_(S_auto, iSite1);
miKnn1 = miKnn1(1:min(NUM_KNN, size(miKnn1,1)), :);
% [viLim_drift, nDrift, viClu, nClu, nSpk_min, vrThresh_site, mrPv, mlDrift] = ...
%     get_(S_auto, 'viLim_drift', 'nDrift', 'viClu', 'nClu', 'nSpk_min', 'vrThresh_site', 'mrPv', 'mlDrift');
thresh1 = vrThresh_site(iSite1);
[trPc1, viSpk1] = load_fet_site_(S_auto, 1, iSite1);
% cc1_drift_clu = cell(nDrift, nClu);
cvii1_drift = vi2cell_(discretize(viSpk1, viLim_drift), nDrift);
[vrRho1, viClu1] = deal(S_auto.vrRho(viSpk1), viClu(viSpk1));
vrRho_1 = copy_mask_(S_auto.vrRho, viSpk1);
[~, miiKnn1] = ismember(miKnn1, viSpk1);

if fUseSecondSite
    [trPc2, viSpk2] = load_fet_site_(S_auto, 2, iSite1);
else
    trPc2 = [];
end
fSecondSite = ~isempty(trPc2);
if fSecondSite   
    vrRho_2 = copy_mask_(S_auto.vrRho, viSpk2);
    [~, miiKnn2] = ismember(miKnn1, viSpk2);
end
% cc2_drift_clu = cell(nDrift, nClu);
cviClu_drift = cell(nDrift,1);
ctrPc_drift = cell(nDrift,1);

for iDrift = 1:nDrift    
    vii1 = cvii1_drift{iDrift};
    if isempty(vii1), continue; end
    [vrRho11, viClu11, miKnn11, miiKnn11] = ...
        deal(vrRho1(vii1), viClu1(vii1), miKnn1(:,vii1), miiKnn1(:,vii1));
    [cviiSpk_clu_, ~, viClu_uniq] = vi2cell_(viClu11, nClu);
    if fSecondSite, miiKnn21 = miiKnn2(:,vii1); end
    [viClu_drift1, cmrPc_drift1] = deal([], {}); 
    for iClu = viClu_uniq
        vii_ = cviiSpk_clu_{iClu};
        [miKnn11_, miiKnn11_] = deal(miKnn11(:,vii_), miiKnn11(:,vii_));
        vrRho11_T = vrRho11(vii_)';
        vii1_ = miiKnn11_(vrRho_1(miKnn11_) >= vrRho11_T);  
        mrPc1 = mean_conditional_(trPc1, vii1_, nSpk_min, mrPv, thresh1);
        if ~isempty(mrPc1)
            viClu_drift1(end+1) = iClu;
            cmrPc_drift1{end+1} = mrPc1;
        end
        if fSecondSite
            miiKnn21_ = miiKnn21(:,vii_);
            vii2_ = miiKnn21_(vrRho_2(miKnn11_) >= vrRho11_T);
            mrPc2 = mean_conditional_(trPc2, vii2_, nSpk_min, mrPv, thresh1);
            if ~isempty(mrPc2)
                viClu_drift1(end+1) = iClu;
                cmrPc_drift1{end+1} = mrPc2;
            end
        end         
    end
    cviClu_drift{iDrift} = viClu_drift1;
    ctrPc_drift{iDrift} = cat(3, cmrPc_drift1{:});
end
[trPc1, trPc2, miKnn1, miiKnn1, miiKnn2] = deal([]); % clear memory
vlExist_clu = false(1, nClu);
vlExist_clu([cviClu_drift{:}]) = true;
mrDist_clu = zeros(nClu, 'single');

norm_mr_ = @(mr)mr ./ sqrt(sum(mr.^2,1)); 
tr2mr_pv_norm_ = @(tr,mr)norm_mr_(reshape(mr*reshape(tr,size(tr,1),[]),[],size(tr,3))); 

% distance calculation
for iDrift = 1:nDrift
    viDrift1 = find(mlDrift(:,iDrift));
    viClu1 = cviClu_drift{iDrift};
    if isempty(viClu1), continue; end
    trPc_clu1 = cat(3, ctrPc_drift{iDrift});
    if isempty(trPc_clu1), continue; end
    viClu2 = [cviClu_drift{viDrift1}];
    trPc_clu2  = cat(3, ctrPc_drift{viDrift1});
    if isempty(trPc_clu2), continue; end
    mrWav_clu2 = tr2mr_pv_norm_(trPc_clu2, mrPv);
    for iiClu1 = 1:numel(viClu1)
        iClu1 = viClu1(iiClu1);
        mrWav11 = pc2wav_shift_(trPc_clu1(:,:,iiClu1), mrPv, viShift);
        for iiClu2 = 1:numel(viClu2)
            iClu2 = viClu2(iiClu2);         
            if iClu2 > iClu1 % symmetric
                dist12 = max(mrWav_clu2(:,iiClu2)' * mrWav11);
                mrDist_clu(iClu2, iClu1) = max(mrDist_clu(iClu2, iClu1), dist12);
            end
        end
    end
end

fprintf('.');
end %func


%--------------------------------------------------------------------------
function mrWav = pc2wav_shift_(trPc, mrPv, viShift)

if isempty(trPc), mrWav=[]; return; end
[nPc, nSites, nSpk, nT] = deal(size(trPc,1), size(trPc,2), size(trPc,3), size(mrPv,1));
trWav = reshape(mrPv*reshape(trPc, nPc,[]), [nT, nSites, nSpk]);
mrWav = reshape(shift_trWav_(trWav, viShift), nT*nSites, []);
mrWav = mrWav ./ sqrt(sum(mrWav.^2,1)); 
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
if nargin<1, S0 = get(0, 'UserData'); end
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
iSite = S0.S_clu.viSite_clu(S0.iCluCopy);
% [vrFet0, vrTime0] = getFet_site_(iSite, [], S0);    % plot background    
[vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, S0.iCluCopy, S0); % plot iCluCopy
if isempty(S0.iCluPaste)
%     vrFet = [vrFet0(:); vrFet1(:)];
    cvrFet = {vrFet1};
else
    [vrFet2, vrTime2, vcYlabel, viSpk2] = getFet_site_(iSite, S0.iCluPaste, S0); % plot iCluCopy
%     vrFet = [vrFet0(:); vrFet1(:); vrFet2(:)];
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
    rescale_FigTime_(S_fig_time.maxAmp, S0, S0.P);
end
end %func


%--------------------------------------------------------------------------
function [hFig, S_fig] = plot_FigWavCor_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
S_clu = S0.S_clu; P = S0.P;
[hFig, S_fig] = get_fig_cache_('FigWavCor'); 

figure_wait_(1, hFig);
nClu = S_clu.nClu;
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
    S_fig.hImWavCor = imagesc(S_clu.mrWavCor, P.corrLim);  %clears title and current figure
    S_fig.hCursorV = line([1 1], [.5 nClu+.5], 'Color', [0 0 0], 'LineWidth', 1.5); 
    S_fig.hCursorH = line([.5 nClu+.5], [1 1], 'Color', [1 0 0], 'LineWidth', 1.5);             
    colorbar(S_fig.hAx);
    S_fig.vcTitle = '[S]plit; [M]erge; [D]elete';
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigWavCor_);
    mouse_figure(hFig, S_fig.hAx, @button_FigWavCor_);
    S_fig.hDiag = plotDiag_([0, nClu, .5], 'Color', [0 0 0], 'LineWidth', 1.5);
else
    set(S_fig.hImWavCor, 'CData', S_clu.mrWavCor);
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
set(0, 'UserData', S0);

auto_scale_proj_time_(S0);
figure_wait_(0, hFig_wait);

plot_raster_(S0); %psth
ui_show_elective_();
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
function trWav_clu = S0_trWav_clu_(S0)
nSamples_max = 2^10;

cviSite_spk_clu = cellfun_(@(x)S0.viSite_spk(x), S0.S_auto.cviSpk_clu);
cviSite_clu = arrayfun_(@(x)x,S0.S_auto.viSite_clu);
cviSpk_clu = cellfun_(@(x,y,z)x(y==z), S0.S_auto.cviSpk_clu, cviSite_spk_clu, cviSite_clu);

trWav_clu = cellfun_(@(x)mean(...
    pc2wav_(S0.mrPv_global, S0.trPc_spk(:,:,subsample_vr_(x,nSamples_max))),3), ...
        cviSpk_clu);
trWav_clu = cat(3, trWav_clu{:});
end %func


%--------------------------------------------------------------------------
function S_fig = plot_tnWav_clu_(S_fig, P, S0)
% Substituting plot_spk_
% S0 = get(0, 'UserData'); 
S_auto = S0.S_auto;
if ~isfield(P, 'LineWidth'), P.LineWidth=1; end
% trWav_clu = ifeq_(P.fWav_raw_show, S_auto.trWav_raw_clu, S_auto.trWav_spk_clu);
trWav_clu = S0_trWav_clu_(S0);
[nSamples, nSites, nClu] = size(trWav_clu);
nChans_show = size(P.miSites, 1);
miSites_clu = P.miSites(:, S_auto.viSite_clu);
% nSites = numel(P.viSite2Chan);

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
    viSites1 = miSites_clu(:,iClu);
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
function [viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S_clu, iClu1, viSite_spk)
% get a subset of cluster that is centered
% return only centered spikes
% if nargin<2, S0 = get(0, 'UserData'); end
% S_clu = S0.S_clu;
if nargin<3, viSite_spk = get0_('viSite_spk'); end
iSite_clu1 = S_clu.viSite_clu(iClu1);
viSpk_clu1 = S_clu.cviSpk_clu{iClu1};
viSite_clu1 = viSite_spk(viSpk_clu1);
viiSpk_clu1 = find(viSite_clu1 == iSite_clu1);
viSpk_clu1 = viSpk_clu1(viiSpk_clu1);
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
    rethrow ME;
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
    rethrow ME;
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


%--------------------------------------------------------------------------
% irc.m
function varargout = plot_update_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = plot_group_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = impoly_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = split_clu_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = tnWav2uV_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = tnWav_spk_sites_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
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
function varargout = keyPressFcn_cell_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = S_clu_prune_icl_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end
function varargout = get_fig_cache_(varargin), cell_out = call_irc_(dbstack(), varargin, nargout); varargout = cell_out; end