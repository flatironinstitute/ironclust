function [hbar, herr] = bar_error(mr, width)
if nargin<2, width=0; end

if size(mr,1)>1
    vrY = nanmean(mr);
    vrE = nanstd(mr);
    vrX = 1:size(mr,2);    
else
    vrY = mr;
    vrE = [];
    vrX = 1:numel(vrY);
end
if width==0
    hbar = bar(vrX, vrY, .5);
else
    hbar = bar(vrX+width/2, vrY, width);
end
if ~isempty(vrE)
    hold on;
    herr = errorbar(vrX+width/2, vrY, zeros(size(vrY)), vrE);
    herr.LineStyle = 'none';  
    herr.Color = [0 0 0];                            
    grid on;
end
axis([min(vrX)-.5, max(vrX)+.5, 0, max(vrY)]);
end %func