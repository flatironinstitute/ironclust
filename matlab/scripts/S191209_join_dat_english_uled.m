% goal: join .dat files and export as .mda files and geom.csv
% define constants

csDir = {'m90_190915_125303', 'm90_190916_105341', 'm90_190917_112752', 'm90_190918_105939'};
%vcFile_xml0 = '/mnt/ceph/users/jjun/DanEnglish/PV_ChR2_Chronic/%s/amplifier.xml';
vcFile_xml0 = 'D:\\globus\\DanEnglish\\PV_ChR2_Chronic\\%s\\amplifier.xml';
vcDir_out = 'D:\globus\DanEnglish\PV_ChR2_Chronic';
csFile_mda_shank = arrayfun(@(x)fullfile(vcDir_out, sprintf('shank%d', x), 'raw.mda'), ...
        1:4, 'UniformOutput', 0);
% 8.9, 6.6, 5.7, 3.4 hours

P = struct('nChans_ain', 8, 'viChans_ain', [6,7], 'uV_per_bit', .195, ...
    'MAX_BYTES_LOAD', .5e9, 'vcDataType', 'int16');

% Total duration
% m90_190915_125303
% # stim
%            3        3331        3320
% m90_190917_112752
% # stim
%            0        1062        1060
% m90_190916_105341
% # stim
%                     5266        5240
% m90_190918_105939
% # stim
%            7        3550        3533

% last 1 hour (1e8 samples) duration
% m90_190918_105939
% # stim
%         1495        1496
% m90_190916_105341
% # stim
%         1675        1675
% m90_190917_112752
% # stim
%    999   999
% m90_190918_105939
% # stim
%         1495        1496


%% load and save analysis


% iFile = 3;
v_fid_shank = cellfun(@(x)fopen(x,'w+'), csFile_mda_shank);
for iFile = 1:numel(csDir)
    vcDataID1 = csDir{iFile};
    disp(vcDataID1); % show dataID
    vcFile_xml1 = sprintf(vcFile_xml0, csDir{iFile});
    S_xml1 = irc('call', 'load_xml_neuroscope_', {vcFile_xml1});
    cviSite2chan = S_xml1.cviSite2chan; %{S_xml1.viSite2chan1, S_xml1.viSite2chan2, S_xml1.viSite2chan3, S_xml1.viSite2chan4};
    vcFile_dat1 = S_xml1.vcFile_dat;
    fid_r1 = fopen(vcFile_dat1, 'r');
    S_dir1 = dir(vcFile_dat1);
    nBytes1 = S_dir1.bytes;
    bytes_sample1 = 2 * S_xml1.nChans;
    nSamples1 = floor(nBytes1 / bytes_sample1);
    nSamples_load = floor(min(P.MAX_BYTES_LOAD, nBytes1) / bytes_sample1);
    nSamples_loaded = 0;
    while nSamples_loaded < nSamples1
        if nSamples_loaded + nSamples_load > nSamples1
            nSamples_load = nSamples1 - nSamples_loaded;
        end
        mnWav_T1 = fread(fid_r1, [S_xml1.nChans, nSamples_load], ['*', P.vcDataType]);
        for iShank = 1:numel(cviSite2chan)
            mnWav_T11 = mnWav_T1(cviSite2chan{iShank},:); % temporal downlsample
            writemda_fid(v_fid_shank(iShank), mnWav_T11);
            fprintf('.');
        end
        fprintf(' ');
        nSamples_loaded = nSamples_loaded + nSamples_load;
    end
    fclose(fid_r1);
    fprintf('\n');
    mnWav_T1 = [];    
end
% arrayfun(@(x)fclose(x), v_fid_shank);
arrayfun(@(x)writemda_fid(x, 'close'), v_fid_shank);

%%
for i=1:4 
    irc2('auto', sprintf('D:\\Globus\\DanEnglish\\PV_ChR2_Chronic\\shank%d\\irc2\\raw_geom.prm',i));
end

%% run ironclust
for i=1:4 
    irc2(sprintf('D:\\Globus\\DanEnglish\\PV_ChR2_Chronic\\shank%d',i));
end
%
irc2 D:\Globus\recordings\kf19_nt27\irc2\raw_geom.prm

%% dan english dataset
irc2 