function convert_mda_ui(vcMode, vcDir_in, vcDir_out)

switch vcMode
    case 'english', convert_mda_english_(vcDir_in, vcDir_out);     
    case 'crcns', convert_mda_crcns_(vcDir_in, vcDir_out);            
    case 'buzsaki', convert_mda_buzsaki_intra_(vcDir_in, vcDir_out); % no GUI
    case 'juxta', convert_mda_juxta_(vcDir_in, vcDir_out);  % no GUI
    case 'summarize_study', summarize_study_(vcDir_in);
    otherwise, error('invalid mode');
end %switch
end %func


%--------------------------------------------------------------------------
% 1/11/2019 parse crcns buzsaki lab dataset
% export to mda
% @TODO: merge *.### multi session recordings from 
% @TODO: select larger shank
function convert_mda_english_(vcDir_in, vcDir_out)
if isempty(vcDir_in)
    S_cfg = file2struct_('dan_english.cfg');
    vcDir_in = S_cfg.vcDir_in;
end
if isempty(vcDir_out)
    S_cfg = file2struct_('dan_english.cfg');
    vcDir_out = S_cfg.vcDir_out;
end
dir_ext_ = @(ext)arrayfun_(@(x)fullfile(x.folder, x.name), dir(fullfile(vcDir_in, ext)));
csFiles_xml = dir_ext_(fullfile('*','*.xml'));
vcFile_mat = fullfile(vcDir_out, 'settings_english.mat');
table_data = load_table_paired_(vcFile_mat, csFiles_xml);

S_fig = struct('vcFile_mat', vcFile_mat, 'vcDir_out', vcDir_out);
hFig = create_gui_crcns_(csFiles_xml, table_data, S_fig);   
end %func


%--------------------------------------------------------------------------
% 3/29/2019 JJJ: convert to english format
function S_mda = export_spikeforest_english_(vcFile_rhd1, vcDir_out1, mrLim_incl1, viSite_ext1, vcFile_prb)
% S_mda = export_spikeforest_english_(vcFile_rhd1, vcDir_out1, mrLim_incl1, viChan_ext1)
% S_mda = export_spikeforest_english_(vcFile_rhd1, vcDir_out1, mrLim_incl1, vcFile_prb)

if nargin<2, vcDir_out1 = []; end
if nargin<3, mrLim_incl1 = []; end
if nargin<4, viSite_ext1 = []; end
if nargin<5, vcFile_prb = []; end

S_cfg = file2struct_(ircpath_('dan_english.cfg'));
if isempty(vcFile_prb)    
    vcFile_prb = S_cfg.vcFile_probe;
end
S_prb = load_prb_(vcFile_prb);
mrSiteXY = S_prb.mrSiteXY;

if matchFileExt_(vcFile_rhd1, '.rhd')
    [vcFile_dat1, vcFile_meta, S_meta1] = rhd2bin_(vcFile_rhd1);
    [nChans1, sRateHz1, vcDataType, uV_per_bit] = get_(S_meta1, 'nChans', 'sRateHz', 'vcDataType', 'uV_per_bit');
elseif matchFileExt_(vcFile_rhd1, '.xml')
    vcFile_dat1 = subsFileExt_(vcFile_rhd1, '.dat');
    S_xml1 = xml2struct(vcFile_rhd1);
    S_ = S_xml1.parameters.acquisitionSystem;
    sRateHz1 = str2num_(S_.samplingRate.Text);
    nChans1 = str2num_(S_.nChannels.Text);
    nBits1 = str2num_(S_.nBits.Text);
    voltageRange1 = str2num_(S_.voltageRange.Text);
    amplification1 = str2num_(S_.amplification.Text);
    uV_per_bit = (voltageRange1 / amplification1 * 1e6) / 2 ^ nBits1;
    vcDataType = sprintf('int%d', nBits1);
    cviChan_group = parse_channelGroups_(S_xml1.parameters.anatomicalDescription.channelGroups);
end

mnWav_ext = load_bin_(vcFile_dat1, vcDataType, nChans1);
if ~isempty(mrLim_incl1)
    mnWav_ext = mnWav_ext(:, lim2range_(size(mnWav_ext,2), mrLim_incl1));
end

iChan_intra = get_set_(S_cfg, 'iChan_intra', nChans1);
vnWav_int = mnWav_ext(iChan_intra,:)';
mnWav_ext = mnWav_ext(S_prb.viSite2Chan,:);
if ~isempty(viSite_ext1), mnWav_ext = mnWav_ext(viSite_ext1,:); end

% load stim
viChan_stim = get_set_(S_cfg, 'viChan_stim', 1);
vcLabel_stim = get_set_(S_cfg, 'vcLabel_stim', 'current injection');
try
    nChans_aux = get_set_(S_cfg, 'nChans_aux', 8);
    vcFile_stim1 = subs_file_(vcFile_dat1, 'analogin.dat');
    if exist_file_(vcFile_stim1)
        mnWav_aux = load_bin_(vcFile_stim1, vcDataType, nChans_aux);
        vnWav_stim = sum(mnWav_aux(viChan_stim,:), 1)';
    else
        vnWav_stim = [];
    end
catch
    fprintf(2, 'export_spikeforest_english_: loading `analogin.dat` error: \n\t%s\n', lasterr());
    vnWav_stim= [];
end

S_json = struct('spike_sign', -1, 'samplerate', sRateHz1, 'scale_factor', uV_per_bit);
t_dur1 = filesize_(vcFile_dat1) / nChans1 / sRateHz1 / bytesPerSample_(vcDataType); 
if isa(vnWav_int, 'int16')
    rSaturation = mean(vnWav_int == int16(2^15-1) | vnWav_int == int16(-2^15));
elseif isa(vnWav_int, 'uint16')
    rSaturation = mean(vnWav_int == uint16(2^16-1));
else
    rSaturation = 0;
end

% build message
csMsg = {
    sprintf('channel count: %d', nChans1);
    sprintf('site count: %d', size(mnWav_ext,1));
	sprintf('recording duration: %0.1f sec', t_dur1);
	sprintf('ADC sampling rate: %0.1f', sRateHz1);
	sprintf('scale factor: %0.4f uV/bit', uV_per_bit);
	sprintf('probe file: %s', vcFile_prb);
	sprintf('ADC saturation: %0.2f %%', rSaturation * 100)
    };
csMsg = {cellstr2vc_(csMsg)};

spkLim_ms_intra = get_set_(S_cfg, 'spkLim_ms_intra', [-1,1]);
spkLim1 = round(sRateHz1*spkLim_ms_intra / 1000);
fJuxta = get_set_(S_cfg, 'fJuxta', 1);
[viSpk_gt, S_intra] = detect_spikes_intra_(vnWav_int, vnWav_stim, spkLim1, fJuxta);            

% output
S_mda = makeStruct_(mnWav_ext, vnWav_int, vnWav_stim, mrSiteXY, S_json, viSpk_gt, csMsg, S_intra, vcLabel_stim);
if ~isempty(vcDir_out1) % write to file
    mkdir_(vcDir_out1);
    export_spikeforest_intra_(vcDir_out1, S_intra, S_json);    
    export_spikeforest_(vcDir_out1, mnWav_ext, mrSiteXY, S_json, viSpk_gt);     
end
end %func


%--------------------------------------------------------------------------
function S_mda = export_spikeforest_crcns_(vcFile_dat1, vcDir_out1, mrLim_incl1, viChan_ext1, vcFile_prb)
% S_mda = export_spikeforest_crcns_(vcFile_dat1, vcDir_out1, mrLim_incl1, viChan_ext1)

if nargin<2, vcDir_out1 = []; end
if nargin<3, mrLim_incl1 = []; end
if nargin<4, viChan_ext1 = []; end
if nargin<5, vcFile_prb = []; end

if matchFileExt_(vcFile_dat1, {'.rhd', '.xml'})
    % dan english format
    S_mda = export_spikeforest_english_(vcFile_dat1, vcDir_out1, mrLim_incl1, viChan_ext1, vcFile_prb);
    return;    
end

vcDataType = 'int16';
[site_dist, shank_dist] = deal(25, 250);
% params.json
vcFile_xml = subsFileExt_(vcFile_dat1, '.xml');
S_xml1 = load_xml_neuroscope_(vcFile_xml);
[sRateHz1, nChans1, scale_factor1, viSite2chan1, viSite2chan2, viSite2chan3, nBits1] = ...
    struct_get_(S_xml1, 'sRateHz', 'nChans', 'uV_per_bit', ...
        'viSite2chan1', 'viSite2chan2', 'viSite2chan3', 'nBits');
S_json = struct('spike_sign', -1, 'samplerate', sRateHz1, 'scale_factor', scale_factor1);

% filter out legacy format
if nBits1 ~= 16, return; end
if isempty(viSite2chan2) % find intracellular channel
    viSite2chan2 = viSite2chan1(end-2:end);
    viSite2chan1 = viSite2chan1(1:end-3);
end     

% raw.mda
% [mnWav1, vnWav_int1, viSpk_gt1] = deal([]);
% if ~exist_file_(fullfile(vcDir_out1, 'raw.mda')) || ~exist_file_(fullfile(vcDir_out1, 'firings_true.mda'))
mnWav1 = reshape_(load_bin_(vcFile_dat1, vcDataType), nChans1);
if ~isempty(mrLim_incl1)
    mnWav1 = mnWav1(:, lim2range_(size(mnWav1,2), mrLim_incl1));
end
if ~isempty(viSite2chan3) 
    viSite2chan_int1 = viSite2chan3;
else
    viSite2chan_int1 = viSite2chan2;
end
if numel(viSite2chan_int1) > 1
    vnWav_stim = mnWav1(viSite2chan_int1(end-1),:)'; % load intracelluar recordings
else
    vnWav_stim = [];
end
vnWav_int = mnWav1(viSite2chan_int1(end),:)'; % load intracelluar recordings

% firings_true.mda
if isempty(mnWav1), return; end

spkLim1 = round(sRateHz1*[-1,1] / 1000);
[viSpk_gt, S_intra] = detect_spikes_intra_(vnWav_int, vnWav_stim, spkLim1);            

% shank selection
if isempty(viSite2chan3)
    mnWav_ext = mnWav1(viSite2chan1,:); % use three column map and find first and last
else
    mnWav_ext = mnWav1([viSite2chan1, viSite2chan2],:);
end
if ~isempty(viChan_ext1), mnWav_ext = mnWav_ext(viChan_ext1,:); end

% geom.csv
nSites1 = size(mnWav_ext,1);
switch nSites1
    case 32
        S_prb1 = load_prb_('poly3nn.prb');
        mrSiteXY = S_prb1.mrSiteXY;   
    otherwise
        mrSiteXY = [zeros(1, nSites1); (0:nSites1-1)]' * site_dist;
        if ~isempty(viSite2chan3) % two shanks
            viSite_shank2 = (1:numel(viSite2chan2)) + numel(viSite2chan1);
            mrSiteXY(viSite_shank2,1) = shank_dist;  % shank separation distance
        end
end %switch          

% set the 
t_dur1 = filesize_(vcFile_dat1) / nChans1 / sRateHz1 / bytesPerSample_(vcDataType); 
csMsg = {
    sprintf('channel_count: %d\nt_dur1: %0.1fs\nsRateHz: %0.1f\nscale_factor: %0.4f', ...
        nChans1, t_dur1, sRateHz1, scale_factor1),
    sprintf('channel group 1:\n    %s', sprintf('%d,', viSite2chan1)),
    sprintf('channel group 2:\n    %s', sprintf('%d,', viSite2chan2)),        
	sprintf('channel group 3:\n    %s', sprintf('%d,', viSite2chan3)),
    };
% output
vcLabel_stim = 'Current stim';
S_mda = makeStruct_(mnWav_ext, vnWav_int, vnWav_stim, mrSiteXY, S_json, viSpk_gt, csMsg, S_intra, vcLabel_stim);
if ~isempty(vcDir_out1) % write to file
    mkdir_(vcDir_out1);
    export_spikeforest_intra_(vcDir_out1, vnWav_int, S_intra);
    export_spikeforest_(vcDir_out1, mnWav_ext, mrSiteXY, S_json, viSpk_gt);     
end
end %func


%--------------------------------------------------------------------------
% convert crcns files to mda, only include selected range
% pay attention to geometry files
function cbf_menu_crcns_(h, e)
% read table and convert file
hFig = h.Parent.Parent;
S_fig = hFig.UserData;
[hTbl, iFile, vcDir_out] = struct_get_(S_fig, 'hTbl', 'iFile', 'vcDir_out');
csFiles_dat = hTbl.Data(:,1);
cLim_incl = hTbl.Data(:,3);
cChan_ext = hTbl.Data(:,4);
try
    vcFile_prb = hTbl.Data{iFile,5};
catch
    vcFile_prb = '';
end

vcFile_dat1 = csFiles_dat{iFile};
[vcDir_, vcDir12] = fileparts(vcFile_dat1);
[~, vcDir11] = fileparts(vcDir_);
vcDir_out1 = fullfile(vcDir_out, sprintf('%s_%s', vcDir11, vcDir12));
S_cfg = file2struct_('dan_english.cfg');
switch lower(h.Label)
    case 'edit settings'
        edit_('dan_english.cfg');
    case 'load'
        
    case {'convert and next', 'convert', 'convert and summarize recording'}
        figure_wait_(1, hFig);
        try
            mrLim_incl1 = str2num_(cLim_incl{iFile});
            viChan_ext1 = str2num_(cChan_ext{iFile});            
            export_spikeforest_crcns_(vcFile_dat1, vcDir_out1, mrLim_incl1, viChan_ext1, vcFile_prb);                               
        catch
            errordlg(lasterr());
            return;
        end
        figure_wait_(0, hFig);
        if strcmpi(h.Label, 'convert') || strcmpi(h.Label, 'convert and summarize recording')
            msgbox_(sprintf('Exported to %s', vcFile_dat1), 1);
            summarize_recording_(S_cfg.vcDir_out, vcFile_dat1);
            return; 
        else
            msgbox_(sprintf('Exported to %s\n Advancing to next', vcFile_dat1), 1);
        end
        
    case 'skip and next'
        set_table_(hTbl, iFile, 2, false);
        cbf_table_save_(hTbl);
        
    case 'open output'
        if exist_dir_(vcDir_out1)
            winopen_(vcDir_out1);
        else
            errordlg({'Folder does not exist', vcDir_out1});
        end
        return;
        
    case 'summarize study'
        vcFile_summary = summarize_study_(S_cfg.vcDir_out);
        winopen_(vcFile_summary);
        return;
        
    case 'summarize recording'
        summarize_recording_(S_cfg.vcDir_out, vcFile_dat1);
        return;        
        
    otherwise, return;
end %switch

if iFile < numel(csFiles_dat)
    cbf_table_crcns_(hTbl, iFile+1);
else
    msgbox_('all done', 1);
end
end %func


%--------------------------------------------------------------------------
% 1/11/2019 parse crcns buzsaki lab dataset
% export to mda
% @TODO: merge *.### multi session recordings from 
% @TODO: select larger shank
function convert_mda_crcns_(vcDir_in, vcDir_out)
if isempty(vcDir_in) && isempty(vcDir_out)
    vcDir_in = 'K:\BuzsakiLab\groundtruth_tetrode\';
    vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\crcns';    
    [csFiles_xml, csFiles_dat] = dir_set_(vcDir_in, '.xml', '.dat');
elseif isempty(vcDir_out)
    csFiles_dat = {vcDir_in};
    csFiles_xml = {strrep(vcDir_in, '.dat', '.xml')};
    vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\crcns';    
end
% fPlot = 1; 
min_spikes = 30;

% load timestamps and parse out non-noise clusters
% csCell_merge = {'D14521', 'D14921', 'D15711', 'D15712', 'D16311', 'D16613', ...
%     'D17012', 'D17013', 'D17014', 'D17111', 'D17211', 'D17212', 'D18011', ...
%     'D18021', 'D18711', 'D18712', 'D18811', 'D18911'};
% S_json = struct('spike_sign', -1); % base setting
% vcDataType = 'int16';
vcFile_mat = fullfile(vcDir_out, 'settings_crcns.mat');
table_data = load_table_paired_(vcFile_mat, csFiles_dat, min_spikes);
S_fig = struct('vcFile_mat', vcFile_mat, 'vcDir_out', vcDir_out);
hFig = create_gui_crcns_(csFiles_dat, table_data, S_fig);   
end %func


%--------------------------------------------------------------------------
function [vnSpikes_crcns, cS_crcns] = count_spikes_crcns_(csFiles_dat)
vnSpikes_crnc = zeros(size(csFiles_dat));
cS_crcns = cell(size(csFiles_dat));
for iFile = 1:numel(csFiles_dat)
    try
       S_mda1 = export_spikeforest_crcns_(csFiles_dat{iFile});
       [mrSiteXY, S_json, viSpk_gt, csMsg] = ...
           struct_get_(S_mda1, 'mrSiteXY', 'S_json', 'viSpk_gt', 'csMsg');
       vnSpikes_crcns(iFile) = numel(viSpk_gt);
       cS_crcns{iFile} = makeStruct_(mrSiteXY, S_json, viSpk_gt, csMsg);
    catch
        ;
    end
end %for
end %func


%--------------------------------------------------------------------------
% append new files to the table
function [table_data, vnSpikes_crcns] = load_table_paired_(vcFile_mat, csFiles_xml, min_spikes)
% # Usages
% [table_data] = load_table_paired_(vcFile_mat, csFiles_rhd)
%    for dan english dataset
% [table_data, vnSpikes_crcns] = load_table_paired_(vcFile_mat, csFiles_dat, min_spikes)
%    for crcns dataset
% see also create_gui_crcns_

if nargin<3, min_spikes = []; end
nCols = 5;

% initialize table_data
table_data = cell(numel(csFiles_xml), nCols);
table_data(:,1) = csFiles_xml(:);
table_data(:,2) = {true};
table_data(:,3) = {''};
table_data(:,4) = {''};
table_data(:,5) = {''};
get_end_ = @(x,i)x{end-i};
dir2id_ = @(x)cellfun_(@(x1)get_end_(strsplit(x1,'/'),1), x);
if exist_file_(vcFile_mat)
    S_mat = load(vcFile_mat);
    table_data_old = struct_get_(S_mat, 'table_data');
    [~,viA] = ismember(dir2id_(table_data_old(:,1)), dir2id_(csFiles_xml));
    table_data(viA,1:size(table_data_old,2)) = table_data_old;    
%     vnSpikes_crcns = get_(S_mat, 'vnSpikes_crcns');
    if isempty(min_spikes), return; end
else
    vnSpikes_crcns = [];
end

if ~isempty(min_spikes)
    if isempty(vnSpikes_crcns)
        [vnSpikes_crcns, ~] = count_spikes_crcns_(csFiles_xml);
    end
    table_data = table_data(vnSpikes_crcns >= min_spikes,:);
end

if ~exist_file_(vcFile_mat)
    [vcDir_out,~,~] = fileparts(vcFile_mat);
    mkdir_(vcDir_out); % create a directory if doesn't exist
    if isempty(vnSpikes_crcns)
        save(vcFile_mat, 'table_data');
    else
        save(vcFile_mat, 'table_data', 'vnSpikes_crcns');
    end
end
end % func


%--------------------------------------------------------------------------
% no GUI in this
function convert_mda_buzsaki_intra_(vcDir_in, vcDir_out)
if isempty(vcDir_in)
    vcDir_in = 'K:\globus\Buzsaki\juxtaExtra';
end
if isempty(vcDir_out)
    vcDir_out = 'K:\spikeforest\groundtruth\paired_recordings\neuronexus32c';
end
fMatchProbe = 0; 

[csFiles_JC, csFiles_clu, csFiles_res, csFiles_EC, csFiles_EC_XML] = ...
    dir_set_(vcDir_in, '_JC.dat', '_JC.clu.1', '_JC.res.1', '_EC.dat', '_EC.xml');

% load timestamps and parse out non-noise clusters
S_prb_poly2 = load_prb_('poly2nn.prb');
S_prb_poly3 = load_prb_('poly3nn.prb');
[mrSiteXY_poly2, mrSiteXY_poly3] = deal(S_prb_poly2.mrSiteXY, S_prb_poly3.mrSiteXY);
if fMatchProbe
    [miSites_poly2, mrDist_site_poly2] = mrSiteXY_to_miSites_(mrSiteXY_poly2);
    [miSites_poly3, mrDist_site_poly3] = mrSiteXY_to_miSites_(mrSiteXY_poly3);
    P_bandpass = struct('freqLim', [500 6000], 'freqLim_width', [100 1000], 'fGpu', 0);
    miRank_poly2 = rankorder_mr_(1./mrDist_site_poly2);
    miRank_poly3 = rankorder_mr_(1./mrDist_site_poly3);
end
S_json = struct('spike_sign', -1);

for iFile = 1:numel(csFiles_JC)
    try
        % firings_true.mda
        viTime1 = load(csFiles_res{iFile});
        viClu1 = load(csFiles_clu{iFile});
        viClu1 = viClu1(2:end); % fist value is the total number of clusteres
        viSpk_gt1 = viTime1(viClu1~=0);        
        
        % params.json
        S_xml1 = xml2struct(csFiles_EC_XML{iFile}); % depends on an external file
        sRateHz1 = str2num_(S_xml1.parameters.acquisitionSystem.samplingRate.Text);
        nChans1 = str2num_(S_xml1.parameters.acquisitionSystem.nChannels.Text);
        nBits1 = str2num_(S_xml1.parameters.acquisitionSystem.nBits.Text);
        voltageRange1 = str2num_(S_xml1.parameters.acquisitionSystem.voltageRange.Text);
        amplification1 = str2num_(S_xml1.parameters.acquisitionSystem.amplification.Text);
        scale_factor1 = (voltageRange1 / amplification1 * 1e6) / 2 ^ nBits1;
        cS_channels1 = S_xml1.parameters.anatomicalDescription.channelGroups.group;
        if numel(cS_channels1) > 1, cS_channels1 = cS_channels1{1}; end
        viSite2chan1 = fliplr(cellfun(@(x)str2num_(x.Text), ... % top to bottom order
            cS_channels1.channel)) + 1;
        S_json1 = struct_set_(S_json, 'samplerate', sRateHz1, 'scale_factor', scale_factor1);
        
        % filter out lower channel counts
        assert(nBits1==16, 'convert_mda_buzsaki_: nBits==16');
        if nChans1 ~= 33, continue; end                
        
        % output directory
        vcFile_EC1 = csFiles_EC{iFile}; % read binary file and export as mda. needs nChans        
        [vcDir_, vcDir12] = fileparts(fileparts(vcFile_EC1));
        [~, vcDir11] = fileparts(vcDir_);
        vcDir_out1 = fullfile(vcDir_out, sprintf('%s_%s', vcDir11, vcDir12));
          
        fprintf('\n%d/%d: %s, #ch=%d, sRate=%0.1f, scale_factor1=%0.4f\n\t%s\n', ...
            iFile, numel(csFiles_JC), vcFile_EC1, nChans1, sRateHz1, scale_factor1, ...
                sprintf('%d,', viSite2chan1));        
            
        % raw.mda
        if ~exist_file_(fullfile(vcDir_out1, 'raw.mda'))
            mnWav1 = reshape_(load_bin_(vcFile_EC1, 'int16'), nChans1);
            mnWav1 = mnWav1(viSite2chan1,:); % use three column map and find first and last
        else
            mnWav1 = [];
        end
        export_spikeforest_(vcDir_out1, mnWav1, mrSiteXY_poly3, S_json1, viSpk_gt1);        
    catch
        disp(lasterr());
    end
    
    % plot pairwise signal correlation and compare with two column vs.
    if fMatchProbe
        nSamples11 = min(size(mnWav1,2), round(sRateHz1 * 10)); % look at the first 10s
        mrWav11 = single(mnWav1(:,1:nSamples11)');
        mrWav11 = ms_bandpass_filter_(mrWav11, setfield(P_bandpass, 'sRateHz', sRateHz1));
        mrCorr_site1 = corr_(mrWav11);
        miSort_site1 = rankorder_mr_(mrCorr_site1);
    end    
end %for

end %func


%--------------------------------------------------------------------------
function convert_mda_juxta_(vcFile_dat1, vcRange1)
vcDir_in = fileparts(vcFile_dat1);
vcDir_out1 = fullfile(vcDir_in, 'mda');
nlim1 = str2num_(vcRange1);
vcDataType = 'int16';

S_xml1 = load_xml_neuroscope_(subsFileExt_(vcFile_dat1, '.xml'));
[sRateHz1, nChans1, scale_factor1, viSite2chan1, viSite2chan2, viSite2chan3, nBits1] = ...
    struct_get_(S_xml1, 'sRateHz', 'nChans', 'uV_per_bit', ...
        'viSite2chan1', 'viSite2chan2', 'viSite2chan3', 'nBits');
S_json1 = struct('spike_sign', -1, 'samplerate', sRateHz1, 'scale_factor', scale_factor1);

% filter out legacy format
assert(nBits1 == 16, 'convert_mda_juxta_: nBit == 16');
assert(~isempty(viSite2chan2), 'convert_mda_juxta_: no intracellular spikes detected');

% display info
t_dur1 = filesize_(vcFile_dat1) / nChans1 / sRateHz1 / bytesPerSample_(vcDataType); 
fprintf('\n%s: #ch=%d, t_dur1=%0.1fs, sRateHz=%0.1f, scale_factor1=%0.4f\n', ...
    vcFile_dat1, nChans1, t_dur1, sRateHz1, scale_factor1);
fprintf('\tchannel group 1: %s\n', sprintf('%d,', viSite2chan1));
fprintf('\tchannel group 2: %s\n', sprintf('%d,', viSite2chan2));        
fprintf('\tchannel group 3: %s\n', sprintf('%d,', viSite2chan3));        

% raw.mda
mnWav1 = reshape_(load_bin_(vcFile_dat1, vcDataType), nChans1);
if ~isempty(viSite2chan3) 
    viSite2chan_int1 = viSite2chan3;
else
    viSite2chan_int1 = viSite2chan2;
end
if isempty(nlim1)
    mnWav_ext1 = mnWav1(viSite2chan1,:);
else
    mnWav_ext1 = mnWav1(viSite2chan1,nlim1(1):nlim1(2));
end

% firings_true.mda
if numel(viSite2chan_int1) > 1
    vnWav_stim1 = mnWav1(viSite2chan_int1(end-1),:)'; % load intracelluar recordings
else
    vnWav_stim1 = [];
end
if isempty(nlim1)
    vnWav_int1 = mnWav1(viSite2chan_int1(end),:)'; % load intracelluar recordings
else
    vnWav_int1 = mnWav1(viSite2chan_int1(end),nlim1(1):nlim1(2))'; % load intracelluar recordings
end
spkLim1 = round(sRateHz1*[-1,1] / 1000);
fJuxta = 1;
[viSpk_gt1, S_intra] = detect_spikes_intra_(vnWav_int1, vnWav_stim1, spkLim1, fJuxta);    

% geom.csv
nSites1 = numel(viSite2chan1);
switch nSites1
    case 32
        S_prb1 = load_prb_('poly3nn.prb');
        mrSiteXY1 = S_prb1.mrSiteXY;   
    otherwise
        mrSiteXY1 = [zeros(1, nSites1); (0:nSites1-1)]' * 25;
end %switch

% output 
export_spikeforest_(vcDir_out1, mnWav_ext1, mrSiteXY1, S_json1, viSpk_gt1);        
end %func


%--------------------------------------------------------------------------
function hFig = create_gui_crcns_(csFiles_dat, table_data, S_fig)
if nargin<3, S_fig = []; end

hFig = create_figure_('Fig_crcns', [0 0 .5 1], 'Paired recording selector', 0, 1); 

% create table
nFiles = numel(csFiles_dat);
csColumnName = pad_cs_({'File', 'Convert', 'Include time', 'Include Chan (Vext)', 'Probe file'}, [120, 0, 100, 20, 20]);
vlEditable = [false, true, false, true, true];
cData = cell(nFiles, numel(csColumnName)); 
nCol = min(size(cData,2), size(table_data,2));
nRow = min(size(cData,1), size(table_data,1));
cData(1:nRow,1:nCol) = table_data(1:nRow,1:nCol);

hTbl = uitable(hFig, 'Data', cData, ...
    'ColumnEditable', vlEditable, 'ColumnName', csColumnName, ...
    'Unit', 'Normalized','OuterPosition',[0 .8 .8 .2]);
hText  = uicontrol(hFig, 'Style','text','String','',...
           'Unit', 'Normalized', 'OuterPosition', [.8 .8 .2 .2], ...
           'HorizontalAlignment', 'left', 'BackgroundColor','w');
hText.String = {'Click on a file name to plot'}; % display meta

% create axes
hAx1 = axes(hFig, 'OuterPosition', [0 .7 1 .1]); xylabel_(hAx1, 't','V_int','V_int'); % all subsampled
hAx2 = axes(hFig, 'OuterPosition', [0 .6 1 .1]); xylabel_(hAx2, 't','dV_int/dt','dV_int/dt');    
hAx3 = axes(hFig, 'OuterPosition', [0 .5 1 .1]); xylabel_(hAx3, 't','I_int', 'I_int');
hAx4 = axes(hFig, 'OuterPosition', [0 0 1 .5]); xylabel_(hAx4, 't','V_ext','V_ext (Change scale using UP/DOWN arrows)');
S_fig = struct_append_(S_fig, makeStruct_(hTbl, hAx1, hAx2, hAx3, hAx4, hText));
hFig.UserData = S_fig;

% create menu
set(hFig, 'MenuBar', 'none'); 
mh_convert = uimenu_(hFig,'Label','Action'); 
uimenu_(mh_convert,'Label', 'convert and next', 'Callback', @cbf_menu_crcns_);
uimenu_(mh_convert,'Label', 'convert', 'Callback', @cbf_menu_crcns_);
uimenu_(mh_convert,'Label', 'skip and next', 'Callback', @cbf_menu_crcns_);
uimenu_(mh_convert,'Label', 'open output', 'Callback', @cbf_menu_crcns_);
uimenu_(mh_convert,'Label', 'summarize study', 'Callback', @cbf_menu_crcns_);
uimenu_(mh_convert,'Label', 'summarize recording', 'Callback', @cbf_menu_crcns_);

% add call back functions 
hTbl.CellSelectionCallback = @cbf_table_crcns_;
hTbl.CellEditCallback = @cbf_table_save_;
hTbl.BusyAction  = 'cancel';
cbf_table_crcns_(hTbl, 1);
end %func


%--------------------------------------------------------------------------
% callback function for selecting a table on crcns figure
function cbf_table_crcns_(hTbl, arg2)
% cbf_table_crcns_(hTbl, event)
% cbf_table_crcns_(hTbl, iFile)

nPlot_max = 1000000;
hFig = hTbl.Parent;

if nargin<2, arg2 = struct_get_(hFig.UserData, 'iFile'); end
try
    [iFile, iCol] = deal(arg2.Indices(1), arg2.Indices(2));
catch
    [iFile, iCol] = deal(arg2, 1);
end
if iCol~=1, return; end
try
    vcFile_dat1 = hTbl.Data{iFile,1};
catch
    return;
end

figure_wait_(1, hFig);
S_fig1 = hFig.UserData;
[hAx1, hAx2, hAx3, hAx4, hText] = ...
    struct_get_(S_fig1, 'hAx1', 'hAx2', 'hAx3', 'hAx4', 'hText');
cellfun(@(x)cla(x), {hAx1, hAx2, hAx3, hAx4});
cellfun(@(x)hold(x, 'on'), {hAx1, hAx2, hAx3, hAx4});
try
    vcFile_prb = hTbl.Data{iFile,5};
catch
    vcFile_prb = []; 
end
if isempty(vcFile_prb)
    S_cfg = file2struct_(ircpath_('dan_english.cfg'));
    vcFile_prb = S_cfg.vcFile_probe;
    hTbl.Data{iFile,5} = vcFile_prb;
end
S_mda = export_spikeforest_crcns_(vcFile_dat1, [], [], [], vcFile_prb);
% set_userdata_(hFig, S_mda);

[mnWav_ext, vnWav_int, vnWav_stim, mrSiteXY, S_json, viSpk_gt, csMsg, S_intra] = ...
    struct_get_(S_mda, 'mnWav_ext', 'vnWav_int', 'vnWav_stim', 'mrSiteXY', 'S_json', ...
        'viSpk_gt', 'csMsg', 'S_intra');
nTime = numel(vnWav_int);
if isempty(vnWav_stim), vnWav_stim = zeros(size(vnWav_int)); end
[vrFilt_intra, thresh_intra] = struct_get_(S_intra, 'vrWav_filt', 'thresh');
S_fig1.hText.String = csMsg;
nSkip_plot = min(4, max(1, floor(numel(vnWav_int)/nPlot_max)));
switch 2
    case 2
        nSkip_plot = 1;
        [vrX_plot, vrY1_plot, vrY2_plot, vrY3_plot] = deal((1:numel(vnWav_int))', vnWav_int, vrFilt_intra, vnWav_stim);        
        plot__ = @fastplot;
    case 1
        [vrX_plot, vrY1_plot, vrY2_plot, vrY3_plot] = plot_subsample_(nSkip_plot, viSpk_gt, vnWav_int, vrFilt_intra, vnWav_stim);
        plot__ = @plot;
end
xlim0 = [0, vrX_plot(end)];    

sRateHz = S_json.samplerate;
fh_title = @(x,y)title(hAx1, sprintf('%d spikes, %0.3f Hz', numel(x), numel(x)/lim2duration_(y)*sRateHz));
fh_axis = @(h,y)set(h, 'XLim', xlim0, 'YLim', [min(y), max(y)+eps()]);

vhPlot1 = plot__(hAx1, vrX_plot, vrY1_plot, 'k-', viSpk_gt, vnWav_int(viSpk_gt), 'k.');
xylabel_(hAx1, '', 'V_int');
fh_axis(hAx1, vrY1_plot);

vhPlot2 = plot__(hAx2, vrX_plot, vrY2_plot, 'k-', viSpk_gt, vrFilt_intra(viSpk_gt), 'k.');
plot(hAx2, [1, numel(vrFilt_intra)], -double(thresh_intra) * [1,1], 'r-');
vcLabel_filt = get_set_(S_intra, 'vcLabel', 'd/dt(V_int)');
xylabel_(hAx2, '', vcLabel_filt);
fh_axis(hAx2, vrY2_plot);

vhPlot3 = plot__(hAx3, vrX_plot, vrY3_plot, 'k-', viSpk_gt, vnWav_stim(viSpk_gt), 'k.'); 
xylabel_(hAx3, 'Time (adc)', S_mda.vcLabel_stim);
fh_axis(hAx3, vrY3_plot);

mrLim_incl = str2num_(hTbl.Data{iFile,3});
fSave = 0;
if isempty(mrLim_incl)
    mrLim_incl = [1, nTime];
    set_table_(hTbl, iFile, 3, lim2str_(mrLim_incl));
    viSpk_plot = 1:numel(viSpk_gt);
    fSave = 1;
else
    viSpk_plot = vr_select_lim_(viSpk_gt, mrLim_incl);
end
fh_title(viSpk_gt(viSpk_plot), mrLim_incl); % update title

viChan_ext = str2num_(hTbl.Data{iFile,4});
if isempty(viChan_ext)
    set_table_(hTbl, iFile, 4, num2str_(1:size(mnWav_ext,1)));
    fSave = 1;
end

fMeanSubt_ext = size(mnWav_ext,1) >= 16;
mnWav_ext = filt_car_gpu_(mnWav_ext', fMeanSubt_ext, 1);

viX_ext = 1:nSkip_plot:size(mnWav_ext,1);
if nSkip_plot>1
    scale_ext = double(max(max(mnWav_ext) - min(mnWav_ext))) / 2;
    hPlot = multiplot(hAx4, scale_ext, viX_ext, mnWav_ext(viX_ext,:));
else
    hPlot = plot__(hAx4, viX_ext, mnWav_ext(viX_ext,:));
end
set(hAx4, 'YLim', [0, size(mnWav_ext,2)+1], 'YTick', 1:size(mnWav_ext,2), 'XLim', xlim0);
hFig.KeyPressFcn = @(h,e)multifun_cbf_(h,e,hPlot);

% add context menu
c_ = uicontextmenu(hFig);
cPlot_spk = plot2plot_([vhPlot1(2), vhPlot2(2), vhPlot3(2)], 'r.', viSpk_plot);
c_.UserData = makeStruct_(cPlot_spk, fh_title);
uimenu_(c_,'Label','include select','Callback',@cbf_axes_box_);
uimenu_(c_,'Label','exclude select','Callback',@cbf_axes_box_);
uimenu_(c_,'Label','include all','Callback',@cbf_axes_box_);
uimenu_(c_,'Label','exclude all','Callback',@cbf_axes_box_);
% uimenu_(c_,'Label','exclude spikes','Callback',@cbf_axes_box_,'Separator', 'on');
uimenu_(c_,'Label','zoom in x','Callback',@reset_view_,'Separator', 'on');
uimenu_(c_,'Label','set x lim','Callback',@reset_view_,'Separator', 'on');
uimenu_(c_,'Label','zoom in y','Callback',@reset_view_);
uimenu_(c_,'Label','zoom in xy','Callback',@reset_view_);
uimenu_(c_,'Label','center x','Callback',@reset_view_,'Separator', 'on');
uimenu_(c_,'Label','zoom out x','Callback',@reset_view_);
uimenu_(c_,'Label','reset x','Callback',@reset_view_,'Separator', 'on');
uimenu_(c_,'Label','reset y','Callback',@reset_view_);
uimenu_(c_,'Label','reset xy','Callback',@reset_view_);
uimenu_(c_,'Label','redraw','Callback',@(h,e)cbf_table_crcns_(hTbl), 'Separator', 'on');
arrayfun(@(x)set(x, 'UIContextMenu', c_, 'BusyAction', 'cancel'), [hAx1, hAx2, hAx3, hAx4]);

% link plots
cellfun(@(x)grid(x,'on'), {hAx1, hAx2, hAx3, hAx4});
cellfun(@(x)set_userdata_(x, 'axis0', [double(x.XLim), double(x.YLim)]), {hAx1, hAx2, hAx3, hAx4}, 'UniformOutput', 0);
set_userdata_(hFig, 'iFile', iFile, 'nTime', nTime, 'sRateHz', sRateHz);
set(hFig, 'Name', sprintf('(#%d) crcns data selector: %s', iFile, vcFile_dat1));
if fSave, cbf_table_save_(hTbl); end
linkaxes([hAx1, hAx2, hAx3, hAx4], 'x');
drawnow_();
figure_wait_(0, hFig);
end %func


%--------------------------------------------------------------------------
function cbf_table_save_(hTbl, callbackdata)
% update selected range on the table
hFig1 = hTbl.Parent;
table_data = hTbl.Data;
save(get_userdata_(hFig1, 'vcFile_mat'), 'table_data', '-append');
end %func


%--------------------------------------------------------------------------
% draw a box and restrict the plot selection
function cbf_axes_box_(source, callbackdata)
[cPlot_spk, fh_title] = struct_get_(source.Parent.UserData, 'cPlot_spk', 'fh_title');

switch source.Label
    case 'include select', fInclude = 1; xlim=[];
    case 'exclude select', fInclude = 0; xlim=[];
    case 'include all', fInclude = 1; xlim=[1, inf];
    case 'exclude all', fInclude = 0; xlim=[1, inf];
    otherwise, return;
end %switch

% msgbox('Draw a rectangle and double click', 'modal');
if isempty(xlim)
    [xlim, ~] = uirect_();
    if isempty(xlim), return; end
end

% update selected range on the table
hFig1 = source.Parent.Parent;
S_fig1 = hFig1.UserData;
[iFile, hTbl, hAx1, nTime] = struct_get_(S_fig1, 'iFile', 'hTbl', 'hAx1', 'nTime');
mrLim_incl = str2num_(hTbl.Data{iFile, 3});
mrLim_incl = update_lim_(nTime, mrLim_incl, xlim, fInclude);
set_table_(hTbl, iFile, 3, lim2str_(mrLim_incl));

% update display
figure_wait_(1, hFig1);
viPlot = [];
for iPlot = 1:numel(cPlot_spk)
    hPlot_ = cPlot_spk{iPlot};
    cXY_ = hPlot_.UserData;
    [vx_, vy_] = deal(cXY_{1}, cXY_{2});
    if isempty(viPlot)
        viPlot = vr_select_lim_(vx_, mrLim_incl);
        fh_title(vx_(viPlot), mrLim_incl);
    end
    set(hPlot_, 'XData', vx_(viPlot), 'YData', vy_(viPlot));
end %for

% save table data
cbf_table_save_(hTbl);
drawnow_(); 
figure_wait_(0, hFig1);
end % func


%--------------------------------------------------------------------------
% 1/18/2019 detect intracellular spikes
function [viSpk1, S_intra] = detect_spikes_intra_(vnWav_int1, vnWav_stim1, spkLim, fJuxta)
if nargin<3, spkLim = []; end
if nargin<4, fJuxta = 0; end
if isempty(spkLim), spkLim = [-30 30]; end
[min_stim_step, thresh_stim_mad] = deal(10, 10);

vrWav_int1 = single(vnWav_int1);
% mad_int1 = mad_nonzero_(vrWav_int1);

if fJuxta
    qqFactor = 8;
    fft_thresh = 8;
    vrWav_int1 = fft_clean_(vrWav_int1, struct('fft_thresh', fft_thresh, 'fGpu', 1));   
    vrWav_filt = -ms_bandpass_filter_(vrWav_int1, ...
        struct('freqLim', [300 6000], 'sRateHz', 30000, 'freqLim_width', [100 1000]));
    vcLabel = '-bandpass(V_int)';
else
    qqFactor = 12;
    vrWav_filt = -ndiff_(vrWav_int1, 4);
    vcLabel = 'd/dt(V_int)';
end
nRefrac = diff(spkLim);
[viSpk1, vrSpk1, thresh] = spikeDetectSingle_fast_(vrWav_filt, struct('qqFactor', qqFactor));
if 0
    vrAmp_int_spk1 = abs(vrWav_int1(viSpk1)) / mad_nonzero_(vrWav_int1,0);
%     vrAmp_filt_spk = abs(vrWav_filt(viSpk1)) / mad_nonzero_(vrWav_filt,0);
    viiKeep_spk1 = find(vrAmp_int_spk1 >= qqFactor / .6745);
    [viSpk1, vrSpk1] = deal(viSpk1(viiKeep_spk1), vrSpk1(viiKeep_spk1));
%     figure; plot(vrAmp_int_spk, vrAmp_filt_spk, '.');
end

[viSpk1, vrSpk1] = spike_refrac_(viSpk1, vrSpk1, [], nRefrac); %same site spikes

% remove step artifact
vlKeep_spk = true(size(viSpk1));
if ~isempty(vnWav_stim1) && ~isempty(viSpk1)
    max_step_stim1 = max(abs(diff(vnWav_stim1)));
    if max_step_stim1 > min_stim_step % if significant current stimulation is given remove artifact
        vrPre_spk1 = median(vr2mr3_(vnWav_stim1,viSpk1, [spkLim(1), 0]));
        vrPost_spk1 = median(vr2mr3_(vnWav_stim1,viSpk1, [0, spkLim(2)]));
        mad_stim = single(mad_nonzero_(vnWav_stim1));
        vrStim_mad_spk1 = abs(single(vrPre_spk1) - single(vrPost_spk1)) / mad_stim;
        vlKeep_spk(vrStim_mad_spk1 >= thresh_stim_mad) = 0;
    end
end

% remove negative intracelluar spike artifact
if ~fJuxta && ~isempty(viSpk1)
    viSpk_pre1 = max(viSpk1 - spkLim(2), 1);
    vr_dAmp = vnWav_int1(viSpk1) - vnWav_int1(viSpk_pre1); % should be positive
    vlKeep_spk(vr_dAmp<0) = 0;
    viSpk1 = viSpk1(vlKeep_spk);
end

% get amplitudes
vrAmp_filt_spk = vrWav_filt(viSpk1);
vrAmp_int_spk = vrWav_int1(viSpk1);
S_intra = makeStruct_(vrWav_filt, thresh, vcLabel, vrAmp_filt_spk, vrAmp_int_spk);
end %func


%--------------------------------------------------------------------------
function [vcFile_true, vcFile_json] = export_spikeforest_intra_(vcDir_out1, S_intra, S_json)
% create raw_true.mda file 
mkdir_(vcDir_out1);
vcFile_true = fullfile(vcDir_out1, 'raw_true.mda');
vcFile_json = fullfile(vcDir_out1, 'params_true.json');

% write meta
S_true_json = struct('samplerate', S_json.samplerate, 'scale_factor', 1, 'filter', S_intra.vcLabel, 'unit', 'uV');
struct2json_(S_true_json, vcFile_json);

vrWav_intra = single(S_intra.vrWav_filt(:)') * single(-1 * S_json.scale_factor);
% write firings_true.mda
writemda_(vrWav_intra, vcFile_true);
end %func


%--------------------------------------------------------------------------
% 1/9/2019 JJJ: export to spikeforest format
% save raw.mda, geom.csv, params.json, firings_true.mda
function export_spikeforest_(vcDir_out1, mnWav_T1, mrSiteXY1, S_json1, mrGt1)
% export_spikeforest_paired_(vcDir_out1, mnWav1, mrSiteXY1, S_json1, mrGt1)
% S_json: {'samplerate', 'spike_sign', 'scale_factor'}
% mrGt1: 3 x nSpk double matrix. elec#, time index, cluster number

% make directory
mkdir_(vcDir_out1);

% write raw.mda
writemda_(mnWav_T1, fullfile(vcDir_out1, 'raw.mda'));

% write geom.csv
csvwrite(fullfile(vcDir_out1, 'geom.csv'), mrSiteXY1);

% Write to params.json
struct2json_(S_json1, fullfile(vcDir_out1, 'params.json'));

% write firings_true.mda
if ~isempty(mrGt1)
    if size(mrGt1,1) ~= 3
        viTime1 = mrGt1;
        mrGt1 = ones(3, numel(viTime1), 'double');
        mrGt1(2,:) = viTime1;
    end
    writemda_(mrGt1, fullfile(vcDir_out1, 'firings_true.mda'));
end %if
end %func


%--------------------------------------------------------------------------
% assume one GT per recording
function summarize_recording_(vcDir_study, vcFile_recording)
% K:\spikeforest\irc_out\paired_recordings\neuronexus32c
% todo: support firings_true.mda files 
fUseCache = 0;

[~, vcRec] = fileparts(fileparts(vcFile_recording));
[csDir_full, csDir_rec] = subdir_(vcDir_study);
iRec = find(cellfun(@(x)contains(x, vcRec), csDir_rec));
if isempty(iRec), errordlg('mda file is not found. Convert first'); return; end
iRec = iRec(1);

hFig_wait = figure_wait_(1);

% load S_gt1
vcDir1 = csDir_full{iRec};
vcFile_gt1 = fullfile(vcDir1, 'raw_geom_gt1.mat');
vcFile_gt_mda1 = fullfile(vcDir1, 'firings_true.mda');
if exist_file_(vcFile_gt1) && fUseCache
    S_gt1 = load(vcFile_gt1);
elseif exist_file_(vcFile_gt_mda1)        
    S_gt1 = mda2gt1_(vcFile_gt_mda1);
    struct_save_(S_gt1, vcFile_gt1);
else
    disperr_('groundtruth file not found');
end

% display info
[vpp1, vmin1, snr_rms1, snr_min1, nSpk1, iSite1] = ...
    struct_get_(S_gt1, 'vrVpp_clu', 'vrVmin_clu', 'vrSnr_sd_clu', 'vrSnr_min_clu', 'vnSpk_clu', 'viSite_clu');
snr_pp1 = snr_min1 * vpp1/vmin1;    
vnoise1 = mean(S_gt1.vrVrms_site); 
dimm1 = size(S_gt1.trWav_clu);
[nSamples1, nSites1] = deal(dimm1(1), dimm1(2));

% load intracellular traces
vrWav_int = readmda_(fullfile(vcDir1, 'raw_true.mda'));
S_true = loadjson_(fullfile(vcDir1, 'params_true.json'));
[sRateHz, uV_per_bit, vcFilter_intra] = struct_get_(S_true, 'samplerate', 'scale_factor', 'filter');
spkLim = round(read_cfg_('spkLim_ms_gt')/1000 * sRateHz);
mrWav_int_spk = vr2mr3_(vrWav_int, S_gt1.viTime, spkLim);
t_dur1 = numel(vrWav_int) / sRateHz;
vrWav_int_spk = mean(mrWav_int_spk,2);
noise_int = mad_nonzero_(vrWav_int, 0) / .6745;
peak_int = max(abs(vrWav_int_spk));
snr_int = peak_int / noise_int;

% build message
csMsg = {
    sprintf('# Recording info');
    sprintf('  File name: %s', vcRec);
    sprintf('  site count: %d', nSites1);
    sprintf('  ADC sampling rate: %0.1f Hz', sRateHz);
    sprintf('  scale factor: %0.4f uV/bit', uV_per_bit);
	sprintf('  exported duration: %0.1f sec', t_dur1);
	sprintf('  Number of spikes: %d', nSpk1);
	sprintf('# Extracellular info');
	sprintf('  SNR Vpp: %0.4f', snr_pp1);
	sprintf('  SNR Vp: %0.4f', snr_min1);
	sprintf('  noise: %0.1f uV', vnoise1);
	sprintf('  peak: %0.1f uV', vpp1);
	sprintf('  peak-to-peak: %0.1f uV', vmin1);
    sprintf('# Intracellular info');
    sprintf('  SNR Vp: %0.4f', snr_int);
	sprintf('  noise (site average): %0.1f uV', noise_int);
    sprintf('  peak (spike average): %0.1f uV', peak_int);
	sprintf('  Filter: %s', vcFilter_intra);
    };
csMsg = {cellstr2vc_(csMsg)};
disp_cs_(csMsg);

% create a window and display info
vcTitle_fig = ['Paired recording viewer: ', vcRec];
hFig = create_figure_('Fig_paired_recording', [.5 0 .5 1], vcTitle_fig, 0, 1); 
hText  = uicontrol(hFig, 'Style','text','String','',...
           'Unit', 'Normalized', 'OuterPosition', [.8 .5 .2 .5], ...
           'HorizontalAlignment', 'left', 'BackgroundColor','w');
hText.String = csMsg;
hAx1 = axes(hFig, 'OuterPosition', [0 .8 .8 .2]); 
hAx2 = axes(hFig, 'OuterPosition', [0 0 .8 .8]); 
vrT_plot = (spkLim(1):spkLim(2)) / sRateHz * 1000;
plot(hAx1, vrT_plot, vrWav_int_spk); 
plot(hAx2, vrT_plot, S_gt1.trWav_clu);
legend(hAx2, arrayfun(@(x)sprintf('Chan %d', x), 1:nSites1, 'UniformOutput', 0), 'Location', 'southeast');
grid_([hAx1, hAx2], 'on');
linkaxes([hAx1, hAx2], 'x');
xylabel_(hAx1, 'time (ms)','V_int (uV)','Intracellular average spike waveform'); % all subsampled
xylabel_(hAx2, 'time (ms)','V_ext (uV)','Extracellular average spike waveforms');    

% save recording summary
drawnow_();
vcFile_png = fullfile(vcDir1, sprintf('%s_summary.png', vcRec));
saveas(hFig, vcFile_png);
msgbox_(['Saved to ', vcFile_png]);
disp(['Saved to ', vcFile_png]);
figure_wait_(0, hFig_wait);
end %func


%--------------------------------------------------------------------------
% assume one GT per recording
function vcFile_tbl = summarize_study_(vcDir_study)
% K:\spikeforest\irc_out\paired_recordings\neuronexus32c
% todo: support firings_true.mda files 

[csDir_full, csDir_rec] = subdir_(vcDir_study);

nRec = numel(csDir_full);
[vnSpikes, vr_snr_pp, vr_snr_min, vr_snr_rms, vrVpp, vrVmin, vrNoise, viSite] = deal(zeros(nRec,1));
% generate a table output
for iRec = 1:nRec
    vcDir1 = csDir_full{iRec};
    vcFile_gt1 = fullfile(vcDir1, 'raw_geom_gt1.mat');
    vcFile_gt_mda1 = fullfile(vcDir1, 'firings_true.mda');
    if exist_file_(vcFile_gt1)
        S_gt1 = load(vcFile_gt1);
    elseif exist_file_(vcFile_gt_mda1)        
        S_gt1 = mda2gt1_(vcFile_gt_mda1);
        struct_save_(S_gt1, vcFile_gt1);
    else
        disperr_('groundtruth file not found');
    end
    [vpp1, vmin1, snr_rms1, snr_min1, nSpk1, iSite1] = ...
        struct_get_(S_gt1, 'vrVpp_clu', 'vrVmin_clu', 'vrSnr_sd_clu', 'vrSnr_min_clu', 'vnSpk_clu', 'viSite_clu');
    snr_pp1 = snr_min1 * vpp1/vmin1;    
    vnoise1 = mean(S_gt1.vrVrms_site); 
    [vnSpikes(iRec), vr_snr_pp(iRec), vr_snr_min(iRec), vr_snr_rms(iRec), vrVpp(iRec), vrVmin(iRec), vrNoise(iRec), viSite(iRec)] = ...
        deal(nSpk1, snr_pp1, snr_min1, snr_rms1, vpp1, vmin1, vnoise1, iSite1);
end %for
tbl_study = table(csDir_rec(:), vnSpikes, vr_snr_pp, vr_snr_min, vr_snr_rms, vrVpp, vrVmin, vrNoise, viSite, ...
    'VariableNames', {'Recording', 'nSpikes', 'snr_pp', 'snr_min', 'snr_rms', 'Vpp_uV', 'Vmin_uV', 'noise_uV', 'site_peak'});
tbl_study = sortrows(tbl_study, 'snr_pp', 'ascend');
disp(tbl_study);
vcFile_tbl = fullfile(vcDir_study, 'summary_study.xlsx');
writetable(tbl_study, vcFile_tbl);
fprintf('Table saved as %s\n', vcFile_tbl);
end %func


%--------------------------------------------------------------------------
function [miSites, mrDist_site] = mrSiteXY_to_miSites_(mrSiteXY, radius_um)
if nargin<2, radius_um = []; end
mrDist_site = pdist2(mrSiteXY, mrSiteXY);
miSites = findNearSites_(mrSiteXY);
if isempty(radius_um), return; end

nSites_spk = max(sum(mrDist_site <= radius_um));
miSites = miSites(1:nSites_spk,:);
mrDist_site = mrDist_site(1:nSites_spk,:);
end %func


%--------------------------------------------------------------------------
function vc = lim2str_(mrLim)
% mrLIm: nx2
cs = cell(1, size(mrLim,1));
for iRow = 1:size(mrLim,1)
    cs{iRow} = sprintf('[%g,%g]; ', mrLim(iRow,1), mrLim(iRow,2));
end %for
vc = cell2mat(cs);
end %func


%--------------------------------------------------------------------------
% convert limit matrix ([a1,b1;a2,b2;a3,b3]) to sum of durations
% include both edges
function duration = lim2duration_(mrLim_incl)
duration = sum(diff(mrLim_incl')+1);
end %func


%--------------------------------------------------------------------------
function multifun_cbf_(hFig, event, hPlot)
% varargin: plot object to rescale
% Change amplitude scaling 
% change_amp_(event, maxAmp, varargin)
% change_amp_(event) % directly set
% if nargin<3, hPlot=[]; end
factor = sqrt(2);
if key_modifier_(event, 'shift'), factor = factor ^ 4; end
mrAmp_prev = get_userdata_(hPlot, 'scale');
if strcmpi(event.Key, 'uparrow')
    maxAmp = mrAmp_prev / factor;
elseif strcmpi(event.Key, 'downarrow')
    maxAmp = mrAmp_prev * factor;
else
    return;
end
multiplot(hPlot, maxAmp);
end %func


%--------------------------------------------------------------------------
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
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
end %func


%--------------------------------------------------------------------------
% 8/22/18 JJJ: changed from the cell output to varargout
% 9/26/17 JJJ: Created and tested
function varargout = struct_get_(varargin)
% Obtain a member of struct
% cvr = cell(size(varargin));
% if varargin is given as cell output is also cell
S = varargin{1};
for iArg=1:nargout
    vcName = varargin{iArg+1};
    if iscell(vcName)
        csName_ = vcName;
        cell_ = cell(size(csName_));
        for iCell = 1:numel(csName_)
            vcName_ = csName_{iCell};
            if isfield(S, vcName_)
                cell_{iCell} = S.(vcName_);
            end
        end %for
        varargout{iArg} = cell_;
    elseif ischar(vcName)
        if isfield(S, vcName)
            varargout{iArg} = S.(vcName);
        else
            varargout{iArg} = [];
        end
    else
        varargout{iArg} = [];
    end
end %for
end %func


%--------------------------------------------------------------------------
% 1/8/2019 JJJ: returns a set of path that contains all the sets
function varargout = dir_set_(varargin)
% [] = dir_set_(vcDir_in, pattern, replace1, replace2, replace3, ...)
[vcDir_in, vcPostfix1] = deal(varargin{1}, varargin{2});
csFiles1 = dir_(fullfile(vcDir_in, ['**', filesep(), '*', vcPostfix1]));
csPostfix = varargin(2:end);
% cm_files = zeros(numel(csFiles1), numel(csPostfix));
[xx,yy] = meshgrid(csFiles1, csPostfix');
cc = cellfun(@(xx1,yy1)strrep(xx1, vcPostfix1, yy1), xx, yy, 'UniformOutput', 0);
cc = cc(:, all(cellfun(@exist_file_, cc)));
for iArg_out = 1:nargout
    varargout{iArg_out} = cc(iArg_out,:);
end
end %func


%--------------------------------------------------------------------------
function frewind_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function disperr_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function struct_save_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function edit_prm_file_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function mkdir_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function edit_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function figure_wait_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function msgbox_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function set_table_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function winopen_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function xylabel_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end
function drawnow_(varargin), fn=dbstack(); irc('call', fn(1).name, varargin); end


function out1 = meta2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_merge_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ircpath_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = read_cfg_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = file2struct_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = loadjson_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = filesize_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = bytesPerSample_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = fread_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gather_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2ref_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = car_reject_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_copy_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cast_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mr2tr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsample_vr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cell2mat_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_default_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = get_filter_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = set0_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = postCluster_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = gt2mda_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = exist_dir_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = meanSubt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
% function out1 = zscore_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = templateMatch_post_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = memory_matlab_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = calc_drift_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = recording_duration_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = squeeze_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_set_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = S_clu_refresh_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = find_site_spk23_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = mn2tn_wav_spk2_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = matchFileExt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = arrayfun_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cellstr2vc_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = lim2range_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_bin_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_prb_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = load_xml_neuroscope_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = reshape_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subsFileExt_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = str2num_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = corr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = subs_file_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = cellfun_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = ms_bandpass_filter_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = rankorder_mr_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = create_figure_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = struct_append_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = uimenu_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = filt_car_gpu_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end
function out1 = pad_cs_(varargin), fn=dbstack(); out1 = irc('call', fn(1).name, varargin); end


function [out1, out2] = readmda_header_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = mr2thresh_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = gpuArray_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = filt_car_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = findNearSites_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = spatialMask_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end
function [out1, out2] = shift_range_(varargin), fn=dbstack(); [out1, out2] = irc('call', fn(1).name, varargin); end


function [out1, out2, out3] = plan_load_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = detect_spikes_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = cuda_delta_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = cuda_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
% function [out1, out2, out3] = unique_count_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = rho_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = delta_drift_knn_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end
function [out1, out2, out3] = rhd2bin_(varargin), fn=dbstack(); [out1, out2, out3] = irc('call', fn(1).name, varargin); end

