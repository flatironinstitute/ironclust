% memory profile loop
% paramteric test
%   fSave_spkwav {0, 1}
%   nChans, timeDuration
%
%   ToDo: extend time and channels, add noise
%       plan: use compiled matlab (./run_irc myfileloc)

%% 1. Define functions and data

vcDir_from_original = '~/ceph/groundtruth/hybrid_synth/static_siprobe/rec_64c_1200s_11';
vcDir_to_original = '~/raid/groundtruth/hybrid_synth/static_siprobe_bench/rec_64c_1200s';
if ~exist(vcDir_to_original, 'dir')
    system(sprintf('mkdir -p %s', vcDir_to_original));
    vcCmd = sprintf('cp %s %s', fullfile(vcDir_from_original, '*.*'), vcDir_to_original);
    system(vcCmd);
end

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


%% 2. Generate data. Expand or shrink channels 
% Channels: [8, 16, 32, *64*, 128, 256, 512] and times: [300, 600, *1200*, 2400, 4800, 9600]

viChan_pow = -3:2;
viTime_pow = -2:2;
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
        if exist_dir_(vcDir_out1)
            fprintf('%s already exists\n', vcDir_out1); 
            continue; 
        end
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
% setting: Oct 3 2019
% CentOS Workstation, 12 woorkers + GPU (Quadro P4000, 8GB)
% pause(3600*4); % start delay of 4 hours

vnChans_uniq = 64 * 2.^[-3:4];
vrDuration_uniq = 1200 * 2.^[-2:4];

csParam = {};
csParam{1} = 'param1.prm'; % fGpu=0, fParfor=0
csParam{2} = 'param2.prm'; % fGpu=0, fParfor=1
csParam{3} = 'param3.prm'; % fGpu=1, fParfor=0
csParam{4} = 'param4.prm'; % fGpu=1, fParfor=1
switch 1
    case 1, vcVersion = irc2('version');
    case 2, vcVersion = 'v5.0.3';
end

[xx1,yy1] = meshgrid(1:numel(vnChans_uniq), 1:numel(vrDuration_uniq));
vnChans_batch = vnChans_uniq(xx1(:));
vrDuration_batch = vrDuration_uniq(yy1(:));
cS_bench_param = cell(size(csParam));
%for iParam = 1:numel(csParam)
for iParam = 3:3
    try
        vcParam1 = csParam{iParam};
        csFiles_batch_in = arrayfun(@(x,y)...
            fullfile(vcDir0, sprintf('rec_%dc_%ds', vnChans_uniq(x), vrDuration_uniq(y))), ...
                xx1(:), yy1(:), 'UniformOutput', 0);
        [~, vcParam_name1] = fileparts(vcParam1);
        vcPath_out1 = fullfile(sprintf('irc2_%s', vcVersion), vcParam_name1);
        csFiles_batch_out = cellfun(@(x)strrep(x, 'groundtruth', vcPath_out1), csFiles_batch_in, 'UniformOutput', 0);
        cS_bench_param{iParam} = cellfun(@(x,y)irc2('benchmark',x,y,vcParam1), csFiles_batch_in, csFiles_batch_out, 'UniformOutput', 0); % must be transposed for accurate cache result
    catch
        disp(lasterr());
    end
end



% 4. plot result (create two tables), also consider creating a bar plot
% for iParam = 1:numel(csParam)
for iParam = 3:3
    cS_bench = cS_bench_param{iParam};
    % parameter select and plot
    vS_bench = cell2mat(cS_bench);
    vrPeakMem_bench = [vS_bench.memory_gb]; vrPeakMem_bench = vrPeakMem_bench(:);
    vrRuntime_bench = [vS_bench.runtime_sec]; vrRuntime_bench = vrRuntime_bench(:);

    nunique_ = @(x)numel(unique(x));
    title_ = @(x)irc('call','title',{x},1);
    lg = @(x)log(x)/log(2);
    for iMode = 1:2
        switch iMode
            case 1, vrPlot = vrPeakMem_bench(:); vcMode = 'Peak memory (GB)'; vrPlot = log(vrPlot)/log(2);
            case 2, vrPlot = vrRuntime_bench(:); vcMode = 'Runtime (s)'; vrPlot = log(vrPlot)/log(2);
        end
        nChans = vnChans_batch(:); 
        duration_sec = vrDuration_batch(:);
        MB_per_chan_min = vrPlot ./ nChans ./ duration_sec * 1024;
        MB_per_chan = vrPlot ./ nChans  * 1024;
        MB_per_min = vrPlot ./ duration_sec  * 1024;
        table(nChans, duration_sec, vrPlot, MB_per_chan_min, MB_per_chan, MB_per_min, 'rownames', csFiles_batch_in) %{'PeakMem_GiB', 'nChans', 'duration_sec'})

        figure; 
        subplot(2,2,1:2);
        img = reshape(vrPlot, nunique_(duration_sec), nunique_(nChans));
        imagesc(img);
        set(gca, 'XTickLabel', unique(nChans), 'YTickLabel', unique(duration_sec)); % may need to unravel
        xlabel('nChans'); ylabel('min');
    %     set(gca,'XTick', [16 32 64], 'YTick', [10 20]);
        title_(sprintf('%s: %s (irc2 %s)', vcMode, csParam{iParam}, vcVersion));

        subplot 223; plot(img,'b'); 
        xlabel('Duration (s)'); 
        set(gca,'XTickLabel', unique(duration_sec), 'XTick', 1:nunique_(duration_sec));         
        hold on; plot([1, size(img,1)], [1, size(img,1)], 'r');
        set(gca,'YTickLabel', 2.^get(gca,'YTick'), 'YTick', get(gca,'YTick'));
        ylabel(vcMode); grid on; axis tight;

        subplot 224; plot(img','k'); 
        xlabel('#Chans'); 
        set(gca,'XTickLabel', unique(nChans), 'XTick', 1:nunique_(nChans));    
        hold on; plot([1, size(img,2)], [1, size(img,2)], 'r'); 
        set(gca,'YTickLabel', 2.^get(gca,'YTick'), 'YTick', get(gca,'YTick'));
        ylabel(vcMode); grid on; axis tight;       
    end
end


%%
figure('Color','w'); 
vrPlot1 = vrPeakMem_bench(:); vcMode1 = 'Peak memory (GB)'; vrPlot1 = log(vrPlot1)/log(2);
img1 = reshape(vrPlot1, nunique_(duration_sec), nunique_(nChans));

img11 = img1(:,4:end);
ax1=subplot(221); plot(img11,'b'); axis([1,6,-2,8]); grid on; hold on;  
set(gca,'XTickLabel', {}, 'XTick', 1:nunique_(duration_sec));         
plot([1, size(img11,1)], [1, size(img11,1)]+2, 'r');
set(gca,'YTickLabel', 2.^get(gca,'YTick'), 'YTick', get(gca,'YTick'));
ylabel(vcMode1);
ax1.Position = [0.1300    0.5838    0.3262    0.3412] + [0 0 .1 .1];

img12 = img1(3:end,:);
ax2=subplot(222); plot(img12','b'); axis([2,7,-2,8]); grid on; hold on; 
set(gca,'XTickLabel', {}, 'XTick', 1:nunique_(nChans));  
plot([1, size(img12,2)], [1, size(img12,2)]+1, 'r'); 
set(gca,'YTickLabel', {}, 'YTick', get(gca,'YTick'));
ax2.Position = [0.5703    0.5838    0.3262    0.3412] + [0 0 .05 .05];


vrPlot2 = vrRuntime_bench(:); vcMode2 = 'Runtime (s)'; vrPlot2 = log(vrPlot2)/log(2);
img2 = reshape(vrPlot2, nunique_(duration_sec), nunique_(nChans));

img21 = img2(:,4:end);
ax3=subplot(223); plot(img21,'k'); axis([1,6,5,14]); grid on; hold on;  
xlabel('Duration (s)'); 
set(gca,'XTickLabel', unique(duration_sec), 'XTick', 1:nunique_(duration_sec));         
plot([1, size(img21,1)], [1, size(img21,1)]+8, 'r');
set(gca,'YTickLabel', 2.^get(gca,'YTick'), 'YTick', get(gca,'YTick'));
ylabel(vcMode2); 
ax3.Position = [0.1300    0.1100    0.3262    0.3412] + [0 0 .05 .05];

img22 = img2(3:end,:);
ax4=subplot(224); plot(img22','k'); axis([2,7,5,14]); grid on; hold on; 
xlabel('#Chans'); 
set(gca,'XTickLabel', unique(nChans), 'XTick', 1:nunique_(nChans));    
plot([1, size(img22,2)], [1, size(img22,2)]+7, 'r'); 
set(gca,'YTickLabel', {}, 'YTick', get(gca,'YTick'));
ax4.Position = [0.5703    0.1100    0.3262    0.3412] + [0 0 .05 .05];
