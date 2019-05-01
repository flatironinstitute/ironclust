function [S_d1, S_d2] = struct_compare(S1,S2)

csName1 = fieldnames(S1);
csName2 = fieldnames(S2);
csName12 = intersect(csName1, csName2);
[S_d1,S_d2] = deal(struct());
if isempty(csName12), return; end

for iField = 1:numel(csName12)
    vcName_ = csName12{iField};
    val1_ = S1.(vcName_);
    val2_ = S2.(vcName_);
    if ~strcmp(class(val1_), class(val2_))
        fDiff_ = true;
    elseif isempty(val1_) && isempty(val2_)
        fDiff_ = false;
    elseif ischar(val1_)
        fDiff_ = ~strcmpi(val1_, val2_);
    elseif numel(val1_) ~= numel(val2_)
        fDiff_ = true;
    elseif iscell(val1_)
        if numel(val1_) ~= numel(val2_)
            fDiff_ = true;
        else
            for iCell = 1:numel(val1_)
                if numel(val1_{iCell}) ~= numel(val2_{iCell})
                    fDiff_ = true;
                elseif ~all(val1_{iCell} == val2_{iCell})
                    fDiff_ = true;
                end
            end
        end
    else
        fDiff_ = ~all(val1_ == val2_);
    end
    if fDiff_
        S_d1.(vcName_) = val1_;
        S_d2.(vcName_) = val2_;
    end
end %for
end %func