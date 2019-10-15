% data loader
switch 5
    case 1, vcFile = '/mnt/ceph/users/jjun/groundtruth/paired_recordings/english/m14_190326_155432/raw_geom_gt1.mat';        
    case 2, vcFile = '/mnt/ceph/users/jjun/groundtruth/paired_recordings/english/m14_190326_160710_cell1/raw_geom_gt1.mat';
    case 3, vcFile = '/mnt/ceph/users/jjun/groundtruth/paired_recordings/english/m15_190315_142052_cell1/raw_geom_gt1.mat';
    case 4, vcFile = '/mnt/ceph/users/jjun/groundtruth/paired_recordings/english/m15_190315_145422/raw_geom_gt1.mat';
    case 5, vcFile = '/mnt/ceph/users/jjun/groundtruth/paired_recordings/english/m15_190315_150831_cell1/raw_geom_gt1.mat';
    case 6, vcFile = '/mnt/ceph/users/jjun/groundtruth/paired_recordings/english/m15_190315_152315_cell1/raw_geom_gt1.mat';
    case 7, vcFile = '/mnt/ceph/users/jjun/groundtruth/paired_recordings/english/m52_190731_145204_cell3/raw_geom_gt1.mat';
end
S_gt = load(vcFile);    
fLoad_spk = true;
if isfield(S_gt, 'tnWav_spk')
    if ~isempty(S_gt.tnWav_spk)
        fLoad_spk = false;
    end
end
if fLoad_spk
    if ~exist(S_gt.vcFile_spk)
        [~,vcFile1,vcExt1] = fileparts(S_gt.vcFile_spk);
        S_gt.vcFile_spk = fullfile(fileparts(vcFile), [vcFile1, vcExt1]);
    end
    S_gt.tnWav_spk = irc('call', 'load_bin', {S_gt.vcFile_spk, S_gt.type_spk, S_gt.dimm_spk});
else
    vcFile = S_gt.P.vcFile;
end
[~,vcDataID] = fileparts(fileparts(vcFile));

%

sRateHz = 30000;
mrSiteXY = csvread(fullfile(fileparts(vcFile), 'geom.csv'));
vrTime_spk = double(S_gt.viTime) / sRateHz;
nSpk = numel(vrTime_spk);
% trWav_raw = single(S_gt.tnWav_raw);
trWav = single(S_gt.tnWav_spk);
xpos_lim = [min(mrSiteXY(:,1)), max(mrSiteXY(:,1))];
ypos_lim = [min(mrSiteXY(:,2)), max(mrSiteXY(:,2))];
time_lim = [min(vrTime_spk), max(vrTime_spk)]; 

diff_1d_ = @(x)(x([3,3:end,end]) - x([1,1:end-2,end-2]))/2;
diff_2d_ = @(x)(x([3,3:end,end],:) - x([1,1:end-2,end-2],:))/2;
xypos_ = @(x)([sum(x.^2 .* mrSiteXY(:,1),1); sum(x.^2 .* mrSiteXY(:,2),1)] ./ sum(x.^2,1))';
title_ = @(x)title(x,'Interpreter', 'None');

% %% extract pc1 and plot position
% trWav_d = (trWav([3,3:end,end],:,:) - trWav([1,1:end-2,end-2],:,:))/2;
% 
% mrWav_mean = squeeze(mean(trWav,3));
% mrWav_d_mean = squeeze(mean(trWav_d,3));
% figure; plot(mrWav_mean);
% figure; plot(mrWav_d_mean);
% figure; plot(mrWav_mean, mrWav_d_mean);

% track the range
% trWav1 = trWav - mean(trWav,2);
switch 2
    case 1, trWav1 = medfilt1(trWav, 30, [], 3,'truncate');
    case 2, trWav1 = medfilt1(trWav-median(trWav,2), 30, [], 3,'truncate'); % subtract median across sites
end

figure('Color','w','Name', vcFile); 

switch 1
    case 1, mrVpp_spk = squeeze(range(trWav1));
    case 2, trWav1(trWav1>0)=0; mrVpp_spk = squeeze(-min(trWav1));
end
nPix = min(nSpk, 1024);
mrVpp_time = interp1(vrTime_spk, mrVpp_spk', linspace(min(vrTime_spk),max(vrTime_spk),nPix))';
mrDist_time = squareform(pdist(mrVpp_time'));
subplot(3,2,[2,4]); cla; imagesc(log(mrDist_time), [0,8], 'XData', vrTime_spk, 'YData', vrTime_spk); 
colormap jet; title('Extracellular spike similarity'); ylabel('Time (s)');


% [vrVpp_spk, viChMax_spk] = max(mrVpp_spk);
% figure; plot(viChMax_spk, '.');

mrXY_spk = xypos_(mrVpp_spk);
vrA_spk = log(sum(mrVpp_spk.^2))';
ax1=subplot(3,2,1); plot(vrTime_spk, mrXY_spk(:,1), 'k-');  grid on;  ylabel('x pos (um)'); title_(vcDataID);
ax2=subplot(3,2,3); plot(vrTime_spk, mrXY_spk(:,2), 'k-'); grid on; ylabel('y pos (um)');
ax3=subplot(3,2,5); plot(vrTime_spk, vrA_spk, 'k-'); grid on; ylabel('log power'); xlabel('Time (s)');
linkaxes([ax1,ax2,ax3],'x');

vrDist_spk = sqrt(sum((mrXY_spk-median(mrXY_spk)).^2,2));
ax4=subplot(3,2,6); 
switch 2
    case 1, plot(vrTime_spk, cumsum(sqrt(sum(diff_2d_(mrXY_spk).^2,2)))); grid on; ylabel('Accumulated movement (um)');
    case 2, plot(vrTime_spk, vrDist_spk, 'k-'); grid on; ylabel('Displacement (um)');
end
xlabel('Time (s)'); 
% figure; scatter(mrXY_spk(:,1), mrXY_spk(:,2), 1, S_gt.viTime_spk); 
% grid on; axis equal; xlim(xpos_lim); ylim(ypos_lim);
linkaxes([ax1,ax2,ax3,ax4],'x'); xlim(time_lim);


return;

% 
% 
% %% raw peak to peak to plot spatial variable
% nPc_spatial = 3;
% figure; plot(mean(trWav_raw,3),'k'); % spatial averaging reference
% 
% % denoised by subtracting site average and adding back ensemble average
% trWav_raw1 = trWav_raw - mean(trWav_raw,2) + mean(mean(trWav_raw,2),3);
% hold on; plot(mean(trWav_raw1,3),'b'); % spatial averaging reference
% 
% mrVpp_raw1 = squeeze(range(trWav_raw1));
% [mrPc1, mrPv1, c] = pca(mrVpp_raw1, 'NumComponents', nPc_spatial, 'Centered', false);    
% figure; imagesc(mrPv1);
% figure; plot(1-cumsum(c)/sum(c));
% 
% % plot image intensity vs time, keep three pc
% mrVpp_raw2 = mrPv1 * mrPc1'; % denoised vpp
% Y = tsne(mrVpp_raw2');
% figure; scatter(Y(:,1),Y(:,2), 3, 1:size(Y,1));
% 
% %% track centroid movement vs time
% N_MED_FILT = 16;
% 
% [a,b,c] = pca(reshape(trWav, size(trWav,1),[])', 'NumComponents', 1, 'Centered', true);
% figure; plot(1-cumsum(c)/sum(c));
% mrPc = reshape(b, [size(trWav,3), size(trWav,2)]);
% mrPc = medfilt1(mrPc,N_MED_FILT,[],1);
% figure; imagesc(abs(mrPc));
% vrX_spk = sum(mrPc.^2 .* mrSiteXY(:,1)',2) ./ sum(mrPc.^2,2);
% vrY_spk = sum(mrPc.^2 .* mrSiteXY(:,2)',2) ./ sum(mrPc.^2,2);
% %figure; scatter(vrX_spk, vrY_spk, 3, S_gt.viTime_spk);
% vrA_spk = log(sum(mrPc.^2,2));
% figure; 
% ax1=subplot(311); plot(vrTime_spk, vrX_spk, 'b.-');  grid on; 
% ax2=subplot(312); plot(vrTime_spk, vrY_spk, 'r.-'); grid on;
% ax3=subplot(313); plot(vrTime_spk, vrA_spk, 'k.-'); grid on;
% 
% linkaxes([ax1,ax2,ax3],'x');
% xlabel('Time (s)');
% 
% figure; plot3(vrX_spk, vrY_spk, vrA_spk, '.');
% 
% figure; plot(vrTime_spk, cumsum(sqrt(diff_2d_(vrX_spk).^2 + diff_2d_(vrY_spk).^2)), 'k.-');
% xlabel('Time (s)'); ylabel('Cumulative electrode movement (um)');
% 
% %% denoise waveform up to 6 components
% nPc = 6;
% viT = 41:80; %2ms/3 range
% trWav1 = trWav(viT,:,:); 
% if 1
%     trWav1 = trWav1 - mean(trWav1,2); % spatial subtract
% end
% % [~,iSite_max] = max(range(mrWav_mean)) % ch10 max
% [nSamples, nSites, nSpk] = size(trWav);
% trWav_denoise1 = zeros(size(trWav1), 'single');
% for iSite = 1:nSites
%     [mrPc1, mrPv1, c] = pca(squeeze(trWav1(:,iSite,:)), 'NumComponents', nPc, 'Centered', false);    
%     trWav_denoise1(:,iSite,:) = mrPv1 * mrPc1';
% end
% % figure; plot(1-cumsum(c)/sum(c));
% figure; hold on;
% iSpk = 5;
% plot(trWav1(:,:,iSpk),'k');
% plot(trWav_denoise1(:,:,iSpk),'b');
% 
% %% time range
% [~,iSite_max] = max(range(mrWav_d_mean)) % ch10 max
% [mrPc,mrPv_d,c] = pca(squeeze(trWav_d(viT,iSite_max,:)), 'NumComponents', nPc, 'Centered', false);
% figure; plot(1-cumsum(c)/sum(c));
% figure; plot(S_gt.viTime_spk, mrPc(:,1), '.');
% figure; plot(mrPv_d);
% 
% %%
% mrA_spk = squeeze(range(tnWav_spk));
% [vrA_spk, viMaxCh_spk] = max(mrA_spk);
% 
% mrWav_spk = reshape(tnWav_spk, size(tnWav_spk,1), []);
% [mrPc,b,c]=pca(single(mrWav_spk), 'NumComponents', 6);