function mat_ = neCell2mat(cell_)
    %NECELL2MAT Like cell2mat, but keeps only nonempty cells
    nonempty = cellfun(@(x) ~isempty(x), cell_);
    cell_ = cell_(nonempty);
    try
%         mat_ = cell2mat(cell_(nonempty));
        mat_ = cat(1, cell_{:});
    catch
        mat_ = cat(2, cell_{:});
    end
end
