function A_out = applyCellfunskipEmpty(x_struct, var_field_name, sub_field_name, fn_x, ii_inds, jj_inds)
% Applies the function fn_x to contents of x_struct.(var_field_name).(sub_field_name)
% when x_struct.(var_field_name) is non-empty. fn_x must return a scalar
% value. 
% ii_inds and jj_inds allows selecting certain row (ii_inds) and column
% (jj_inds) from the x_struct.(var_field_name) array. 

    if nargin < 6
        ii_inds = 1:size(x_struct.(var_field_name));
        jj_inds = 1:size(x_struct.(var_field_name));
    end

    use_cell_arr = x_struct.(var_field_name)(ii_inds, jj_inds);
    has_entry = cellfun(@(x) ~isempty(x), use_cell_arr);
    
    A_out = nan(size(has_entry));
    A_out(has_entry) = cellfun(@(x) fn_x(x.(sub_field_name)), use_cell_arr(has_entry));
    
end

