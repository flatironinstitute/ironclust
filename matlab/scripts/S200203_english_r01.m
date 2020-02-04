if ispc()
    vcDir = 'D:\Globus\DanEnglish\juxta_cell_curated';
else
    vcDir = '/mnt/ceph/users/jjun/DanEnglish/juxta_cell_curated';
end
S_irc=load(fullfile(vcDir, 'scores_prmset_ironclust.mat'));
S_ms4=load(fullfile(vcDir, 'scores_prmset_mountainsort4.mat'));
S_ks2=load(fullfile(vcDir, 'scores_prmset_kilosort2.mat'));

%% plot and compares
cS_sorter = {S_irc, S_ks2, S_ms4};
ccVar = {'vrF1_gt', 'vrPrecision_gt', 'vrRecall_gt'};
[cvrF1_sorter, cvrPrecision_sorter, cvrRecall_sorter] = deal(cell(size(cS_sorter)));
for iSorter = 1:numel(cS_sorter)
    S_sorter = cS_sorter{iSorter};
    cmr_var = cellfun(@(y)cellfun(@(x)get_nan(x,y), S_sorter.ccScore_prmset_rec), ccVar, 'UniformOutput',0);
    cvrMean_var = cellfun(@(x)nanmean(x,2), cmr_var, 'UniformOutput', 0);
    [~, iPrmset_best] = max(cvrMean_var{1});
    [~, iPrmset_worst] = min(cvrMean_var{1});
    cvrF1_sorter{iSorter} = cmr_var{1}(iPrmset_best,:);
    cvrPrecision_sorter{iSorter} = cmr_var{2}(iPrmset_best,:);
    cvrRecall_sorter{iSorter} = cmr_var{3}(iPrmset_best,:);
end

%%
SNR_THRESH = 6;
mrF1_sorter = cell2mat(cvrF1_sorter')';
mrPrecision_sorter = cell2mat(cvrPrecision_sorter')';
mrRecall_sorter = cell2mat(cvrRecall_sorter')';
vrSnr_rec = [3.5047    2.9944    4.4937    9.8666    4.7681    3.1651    3.6722    4.7454    5.2046    6.8384   19.8035   44.3852    4.6029    7.7731    4.4694 5.8545    2.0890    7.1567    4.3370    1.5844    1.8871    3.8094    2.7723    3.3634    4.4422    0.6719    3.6551    6.9088    4.1878    3.2001    2.8092    3.8817    5.4892];
ml = repmat(vrSnr_rec(:)>=SNR_THRESH, [1,3]);

%
xlabel_ = @()set(gca,'XTick',1:3,'XTickLabel',{'IronClust','KiloSort2','MountainSort4'});
angle_ = @()xtickangle(45);
ylim_ = @()ylim([0 110]);
vcSnr = sprintf(' | SNR>=%0.1f', SNR_THRESH);

figure('Color','w'); 
subplot(1,3,1);
yyaxis left; bar_error(sum(~isnan(mrF1_sorter)),-1/4); ylabel('# Units sorted'); ylim([0 33]); grid on; ylim([0 35]);
yyaxis right; bar_error(sum(mrF1_sorter>=80),1/4); ylabel('# Units | F1>=80'); xlabel_(); angle_(); ylim([0 35]);

mrF1_sorter(~ml)=nan;
mrPrecision_sorter(~ml)=nan;
mrRecall_sorter(~ml)=nan;

subplot(1,3,2); bar_error(mrF1_sorter); ylabel(['F1-score', vcSnr]); xlabel_(); angle_();  ylim_();
subplot(1,3,3); 
yyaxis left; bar_error(mrPrecision_sorter,-1/4); ylabel(['Precision', vcSnr]); ylim_();
yyaxis right; bar_error(mrRecall_sorter,1/4); ylabel(['Recall', vcSnr]); xlabel_(); angle_();  ylim_();


%% do stats
viCompare = find(all(~isnan(mrF1_sorter),2));
[~,vpF1(1)] = ttest2(mrF1_sorter(viCompare,1), mrF1_sorter(viCompare,2));
[~,vpF1(2)] = ttest2(mrF1_sorter(viCompare,1), mrF1_sorter(viCompare,3));
[~,vpF1(3)] = ttest2(mrF1_sorter(viCompare,2), mrF1_sorter(viCompare,3));
disp(vpF1)

[~,vpPrec(1)] = ttest2(mrPrecision_sorter(viCompare,1), mrPrecision_sorter(viCompare,2));
[~,vpPrec(2)] = ttest2(mrPrecision_sorter(viCompare,1), mrPrecision_sorter(viCompare,3));
[~,vpPrec(3)] = ttest2(mrPrecision_sorter(viCompare,2), mrPrecision_sorter(viCompare,3));
disp(vpPrec)

[~,vpRecl(1)] = ttest2(mrRecall_sorter(viCompare,1), mrRecall_sorter(viCompare,2));
[~,vpRecl(2)] = ttest2(mrRecall_sorter(viCompare,1), mrRecall_sorter(viCompare,3));
[~,vpRecl(3)] = ttest2(mrRecall_sorter(viCompare,2), mrRecall_sorter(viCompare,3));
disp(vpRecl)

mean_ = @(mr,x)mean(mr(viCompare,x));
mean_std_ = @(mr,x)[mean(mr(viCompare,x)), std(mr(viCompare,x))];
samplesize_ = @(x,y)sampsizepwr('t2', mean_std_(mrF1_sorter,x), mean_(mrF1_sorter,y));
n12 = min(samplesize_(1,2), samplesize_(2,1))
n13 = min(samplesize_(1,3), samplesize_(3,1))
n23 = min(samplesize_(2,3), samplesize_(3,2))

% vpF1 =
% 
%     0.5659    0.8534    0.4492
% 
% 
% vpPrec =
% 
%     0.9499    0.9648    0.9869
% 
% 
% vpRecl =
% 
%     0.3158    0.9092    0.2253