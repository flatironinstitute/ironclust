% fGpu=0
% Runtime:
%     Detect + feature        40.4s
%     Cluster + auto-merge    49.7s
%         Cluster             26.8s
%         auto-merge          22.9s
%     Total                   90.0s
%     Runtime speed           x3.2 realtime
% Execution:
%     fGpu (GPU use):         0
%     fParfor (parfor use):   0

% fGpu=1
% Runtime:
%     Detect + feature        17.6s
%     Cluster + auto-merge    25.6s
%         Cluster             2.4s
%         auto-merge          23.2s
%     Total                   43.3s
%     Runtime speed           x6.7 realtime
% Execution:
%     fGpu (GPU use):         1
%     fParfor (parfor use):   0

vr_ratio = [40.4/17.6, 26.8/2.4, 90/43.3];
hFig=figure('Color','w','Unit','Normalized','OuterPosition',[.75 0 .25 .5]); 
h=bar(vr_ratio, .5); 
text_ = @(x,y)text(x,y,sprintf('%0.1f',y),...
    'HorizontalAlignment','center','VerticalAlignment','top','Color','w','FontSize',16);
arrayfun(@(x,y)text_(x,y), 1:3, vr_ratio);
h.FaceColor=[0,1,0];
h.LineStyle='none';
set(gca,'XTick',1:3,'XTickLabel', {'Filter', 'Cluster', 'Overall'}, 'XLim', [.5 3.5],'FontSize',16);
grid on; 
ylabel('GPU/CPU speed ratio','FontSize',16);