function mrC2 = ranksort(mrC, fKillDiag)
if nargin<2, fKillDiag=0; end
[~,mrC1]=sort(mrC); [~,mrC2]=sort(mrC1); 
if fKillDiag
    dimm = size(mrC);
    mrC2(sub2ind(dimm, 1:dimm(1), 1:dimm(2)))=nan; 
end
end 