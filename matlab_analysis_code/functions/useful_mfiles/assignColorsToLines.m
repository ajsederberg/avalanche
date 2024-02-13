function plot_handles = assignColorsToLines(plot_handles, color_map_matrix, color_field_name)
% this will work with bar handles if the color field name is provided
% ('FaceColor')
if size(color_map_matrix, 1) < length(plot_handles)
    error('color map matrix should be a Nx3 matrix of rgb colors, with N >= length(plot_handles)')
end

if nargin < 3 
    color_field_name = 'Color';
end
for ii = 1:length(plot_handles)
    try
        plot_handles(ii).(color_field_name) = color_map_matrix(ii, :);
    catch
        set(plot_handles(ii), color_field_name, color_map_matrix(ii, :));
    end
end