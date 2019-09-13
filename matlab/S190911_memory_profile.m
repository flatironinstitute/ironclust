% memory profile loop
% paramteric test
%   fSave_spkwav {0, 1}
%   nChans, timeDuration
%
%   ToDo: extend time and channels, add noise
%       plan: use compiled matlab (./run_irc myfileloc)

%% 1. generate a new recording by tiling
vcDir0 = '~/raid/groundtruth/hybrid_synth/static_siprobe_bench';
vcFile0 = 'rec_64c_1200s';
vcDir_in = fullfile(vcDir0, vcFile0);
nChans0 = 64;
duration0 = 1200; % 20 min
shank_spacing = 250; % shank spacing in um

% 1. load the original (vcDir0) to memory (~4.6GB)
readmda_header_ = @(x)irc('call','readmda_header_',{x});
readmda_ = @(x)irc('call','readmda',{x});
writemda_ = @(x,vcFile)irc('call','writemda',{x,vcFile});
% [S_mda, fid_r] = readmda_header_(fullfile(vcDir0, 'raw.mda'));
fwrite_mda_ = @(x,y,z)irc('call', 'fwrite_mda', {x,y,z});
mkdir_ = @(x)irc('call', 'mkdir', {x});
exist_dir_ = @(x)irc('call', 'exist_dir', {x});


%% 2. expand or shrink channels [8, 16, 32, *64*, 128, 256, 512] and times [300, 600, *1200*, 2400, 4800, 9600]

viChan_pow = -3:3;
viTime_pow = -2:3;
mnWav0 = readmda_(fullfile(vcDir_in, 'raw.mda'));
mrSiteXY0 = csvread(fullfile(vcDir_in, 'geom.csv'));
mrFirings0 = readmda_(fullfile(vcDir_in, 'firings_true.mda'));
nSamples0 = size(mnWav0,2);
nSpk_gt0 = size(mrFirings0, 2);

for iChan_pow = viChan_pow
    nChans1 = nChans0 * (2^iChan_pow);
    % create a channel map
    if iChan_pow==0
        mnWav1 = mnWav0;
        mrSiteXY1 = mrSiteXY0;
    elseif iChan_pow>0
        mnWav1 = [mnWav1; mnWav1];
        mrSiteXY2 = [mrSiteXY1(:,1) + shank_spacing * iChan_pow, mrSiteXY1(:,2)];
        mrSiteXY1 = [mrSiteXY1; mrSiteXY2];
    elseif iChan_pow<0
        mnWav1 = mnWav0(1:nChans1,:);
        mrSiteXY1 = mrSiteXY0(1:nChans1,:);
    end
    assert(size(mnWav1,1) == nChans1, 'nChans must match');
    
    % 3. repeat time
    for iTime_pow = viTime_pow
        t1=tic;
        nRepeat = (2^iTime_pow);
        duration1 = duration0 * nRepeat;
        nSamples1 = int32(size(mnWav1,2)) * nRepeat;        
        vcDir_out1 = fullfile(vcDir0, sprintf('rec_%dc_%ds', nChans1, duration1));
        if exist_dir_(vcDir_out1), continue; end
        mkdir_(vcDir_out1);
        
        % write raw.mda
        fid1 = fwrite_mda_(fullfile(vcDir_out1, 'raw.mda'), [nChans1, nSamples1], class(mnWav0));
        if nRepeat >= 1
            arrayfun(@(x)fwrite(fid1, mnWav1, class(mnWav1)), 1:nRepeat);
        else
            fwrite(fid1, mnWav1(:,1:nSamples1), class(mnWav1));
        end
        fclose(fid1);
        
        % write geom.csv
        csvwrite(fullfile(vcDir_out1, 'geom.csv'), mrSiteXY1);
        
        % copy params.json
        copyfile(fullfile(vcDir_in, 'params.json'), fullfile(vcDir_out1, 'params.json'));        
        
        % write groundtruth        
        if nRepeat >= 1
            fid2 = fwrite_mda_(fullfile(vcDir_out1, 'firings_true.mda'), [3, nSpk_gt0*nRepeat], class(mrFirings0));            
            arrayfun(@(x)fwrite(fid2, mrFirings0+[0;nSamples0*(x-1);0], class(mrFirings0)), 1:nRepeat);
            fclose(fid2);
        else
            iLast1 = find(mrFirings0(2,:) <= nSamples1, 1, 'last');
            writemda_(mrFirings0(:, 1:iLast1), fullfile(vcDir_out1, 'firings_true.mda'));
        end              
        fprintf('Wrote to %s (took %0.1fs)\n', vcDir_out1, toc(t1));
    end
end



%% 3. loop over the files, exract values

% change paramter

% loop over the files
vS_dir1 = dir(fullfile(vcDir0, 'rec_*')); 
csFiles_batch = {vS_dir1.name};
% csFiles_batch = {'rec_16c_1200s_11', 'rec_16c_600s_11', 'rec_32c_1200s_11', ...
%     'rec_32c_600s_11', 'rec_64c_600s_11', ...
%     'rec_64c_1200s_11', 'rec_64c_2400s_11', 'rec_64c_4800s_11', ...
%     'rec_128c_1200s_11', 'rec_128c_2400s_11', 'rec_128c_4800s_11', ...
%     'rec_256c_1200s_11', 'rec_256c_2400s_11', 'rec_256c_4800s_11'};

csParams = {'param_set1.prm', 'param_set2.prm'};

[xx,yy] = meshgrid(1:numel(csFiles_batch), 1:numel(csParams));

fh_bench = @(x,y)irc('benchmark', fullfile(vcDir0, csFiles_batch{x}), csParams{y});
cS_bench = arrayfun(@(x,y)fh_bench(x,y), xx, yy, 'UniformOutput', 0);
vS_bench = cell2mat(cS_bench);
mrPeakMem_batch = [vS_bench.memory_gb];
mrRuntime_batch = [vS_bench.runtime_sec];

cell_index_ = @(x,y)x{y};
str2num_strip_ = @(x)str2double(x(x>='0' & x<='9'));
strsplit2num_ = @(x,y,z)str2num_strip_(cell_index_(strsplit(x, y), z));
vnChans_batch = cellfun(@(x)strsplit2num_(x, '_', 2), csFiles_batch);
vrMinutes_batch = cellfun(@(x)strsplit2num_(x, '_', 3), csFiles_batch) / 60;

%% unshuffle the data being processed


%% 4. plot result (create two tables), also consider creating a bar plot
for iParam = 1:numel(csParams)
    vcParam1 = csParams{iParam};
    
    fprintf('Paramer set: %s\n', vcParam1);
    peakMemory_GiB = mrPeakMem_batch(iParam,:)'; 
    nChans = vnChans_batch'; 
    duration_min = vrMinutes_batch';
    MB_per_chan_min = peakMemory_GiB ./ nChans ./ duration_min * 1024;
    MB_per_chan = peakMemory_GiB ./ nChans  * 1024;
    MB_per_min = peakMemory_GiB ./ duration_min  * 1024;
    table(peakMemory_GiB, nChans, duration_min, MB_per_chan_min, MB_per_chan, MB_per_min, 'rownames', csFiles_batch) %{'PeakMem_GiB', 'nChans', 'duration_sec'})

    img_table1 = peakMemory_GiB;
    img_table1(sub2ind(size(img_table1), vnChans_batch, vrMinutes_batch) = img_table1;

    figure; 
    imagesc(reshape(MB_per_chan_min, 2,3), 'xdata', unique(vnChans_batch), 'ydata', unique(vrMinutes_batch)); % may need to unravel
    xlabel('nChans'); ylabel('min');
    set(gca,'XTick', [16 32 64], 'YTick', [10 20]);
    title(sprintf('Normalized peak memory (MB/chan/min) for %s', vcParam1));
end

% [param_set1.prm: fSave_spkwav=1]
%                         peakMemory_GiB    nChans    duration_min    MB_per_chan_min    MB_per_chan    MB_per_min
%                         ______________    ______    ____________    _______________    ___________    __________
% 
%     rec_16c_1200s_11        0.777           16           20             2.4864           49.728         39.782  
%     rec_16c_600s_11          0.64           16           10              4.096            40.96         65.536  
%     rec_32c_1200s_11        1.003           32           20             1.6048           32.096         51.354  
%     rec_32c_600s_11         0.761           32           10             2.4352           24.352         77.926  
%     rec_64c_1200s_11        1.483           64           20             1.1864           23.728          75.93  
%     rec_64c_600s_11         1.027           64           10             1.6432           16.432         105.16  
%
%
% [param_set1.prm: fSave_spkwav=0]
%                         peakMemory_GiB    nChans    duration_min    MB_per_chan_min    MB_per_chan    MB_per_min
%                         ______________    ______    ____________    _______________    ___________    __________
% 
%     rec_16c_1200s_11        0.541           16           20             1.7312           34.624         27.699  
%     rec_16c_600s_11         0.478           16           10             3.0592           30.592         48.947  
%     rec_32c_1200s_11        0.615           32           20              0.984            19.68         31.488  
%     rec_32c_600s_11         0.559           32           10             1.7888           17.888         57.242  
%     rec_64c_1200s_11        0.693           64           20             0.5544           11.088         35.482  
%     rec_64c_600s_11         0.627           64           10             1.0032           10.032         64.205  


%% parse the output and plot 