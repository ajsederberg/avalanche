function A_out = pullFieldskipEmpty(x_struct, var_field_name, sub_field_name, ii_inds, jj_inds)
% Pulls the scalar value out of x_struct.(var_field_name).(sub_field_name)
% when x_struct.(var_field_name) is non-empty.

    use_cell_arr = x_struct.(var_field_name)(ii_inds, jj_inds);
    has_entry = cellfun(@(x) ~isempty(x), use_cell_arr);
    
    A_out = nan(size(has_entry));
    A_out(has_entry) = cellfun(@(x) x.(sub_field_name), use_cell_arr(has_entry));
    
end