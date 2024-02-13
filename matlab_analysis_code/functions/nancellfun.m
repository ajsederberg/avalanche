function y = nancellfun(cell_arr, field_name)
% extracts the specified field (must be a scalar) from the cell array
% cell_arr
    y = nan(size(cell_arr));
    isnotempty = cellfun(@(x) ~isempty(x), cell_arr);

    y(isnotempty) = cellfun(@(x) x.(field_name), cell_arr(isnotempty));
end