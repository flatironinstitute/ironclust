% test pca compression after realignment

load tnWav_spk_

%%

trWav_spk = tnWav_spk_; % nSamples x nSites_spk x nSpk
dimm_spk = size(trWav_spk);

%% determine peak
nSites_peak = ceil(size(trWav_spk,2)/2);
[~, viMax_spk] = max(sum(trWav_spk(:,1:nSites_peak,:).^2,2));
viMax_spk = viMax_spk(:);

figure; hist(viMax_spk, min(viMax_spk):max(viMax_spk));

%% compute pca VAF: 0.1036 remaining
nPc = 5;
% [a,b,c]=pca(reshape(trWav_spk, dimm_spk(1), [])', 'Centered', false, 'NumComponents', nPc); %mr2=a*b';
% vaf_error = 1-sum(c(1:nPc))/sum(c)
% figure; plot(1-cumsum(c)/sum(c), 'k.-')

% compute pca VAF of central peaks: .088 remaining
iPeak = mode(viMax_spk);
trWav_spk2 = trWav_spk(2:end-1,:, viMax_spk==iPeak);
trWav_spk3 = trWav_spk(2:end-1,:, viMax_spk==iPeak-1);
trWav_spk4 = trWav_spk(2:end-1,:, viMax_spk==iPeak+1);
trWav_spk234 = cat(3, trWav_spk2, trWav_spk3, trWav_spk4);
[a,b,c]=pca(reshape(trWav_spk234, size(trWav_spk234,1), [])', 'Centered', false, 'NumComponents', nPc); %mr2=a*b';
vaf_error = 1-sum(c(1:nPc))/sum(c)
figure; plot(1-cumsum(c)/sum(c), 'k.-')

% compute pca VAF of central peaks: .088 remaining
trWav_spk2 = trWav_spk(2:end-1,:, viMax_spk==iPeak);
trWav_spk3 = trWav_spk(1:end-2,:, viMax_spk==iPeak-1);
trWav_spk4 = trWav_spk(3:end,:, viMax_spk==iPeak+1);
trWav_spk234 = cat(3, trWav_spk2, trWav_spk3, trWav_spk4);
[a,b,c]=pca(reshape(trWav_spk234, size(trWav_spk234,1), [])', 'Centered', false, 'NumComponents', nPc); %mr2=a*b';
vaf_error = 1-sum(c(1:nPc))/sum(c)
hold on; plot(1-cumsum(c)/sum(c), 'b.-')

%% fit model by channel
nPc = 6;
viiSpk2 = find(viMax_spk==mode(viMax_spk));
[a,b,c]=pca(reshape(trWav_spk2, [], size(trWav_spk2,3))', 'Centered', false, 'NumComponents', nPc); %mr2=a*b';
vaf_error = 1-sum(c(1:nPc))/sum(c)
figure; plot(1-cumsum(c)/sum(c), '.-')

%% power pca only
trFft_spk = fft(trWav_spk);
trFft_spk = abs(trFft_spk(2:end/2,:,:));
[a,b,c]=pca(reshape(trFft_spk, size(trFft_spk,1), [])', 'Centered', false, 'NumComponents', nPc); %mr2=a*b';
vaf_error = 1-sum(c(1:nPc))/sum(c)



