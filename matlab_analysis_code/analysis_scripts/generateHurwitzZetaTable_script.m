% Script to generate and save look-up table for the Hurwitz zeta functions

%%%%%% 2022-01-18 versions
% a_min = 1.6;
% a_max = 2.5;
% a_res = 0.005;
% grid_a_vals = a_min:a_res:a_max;
% 
% %     cdf_Paxm = hurwitzZeta(a_hat, x_vals)/hurwitzZeta(a_hat, x_min);
% grid_x_vals = [1:99 unique(round(logspace(2, 4, 101)))];

%%%%%% 2022-10-29 versions
a_min = 1.3;
a_max = 2.9;
a_res = 0.005;
grid_a_vals = a_min:a_res:a_max;

%     cdf_Paxm = hurwitzZeta(a_hat, x_vals)/hurwitzZeta(a_hat, x_min);
grid_x_vals = [1:99 unique(round(logspace(2, 5, 151)))];

grid_HZ = zeros(length(grid_a_vals), length(grid_x_vals));
tic
for ii = 1:length(grid_a_vals)
    this_a = grid_a_vals(ii);
    parfor jj = 1:length(grid_x_vals)
        grid_HZ(ii, jj) = hurwitzZeta(this_a, grid_x_vals(jj));
    end
    toc
end

%% generate for the derivative as well
figure();
colormap jet
grid_dHZ = zeros(length(grid_a_vals), length(grid_x_vals));

for ii = 1:length(grid_a_vals)
    for jj = 1:length(grid_x_vals)
        grid_dHZ(ii, jj) = hurwitzZeta(1, grid_a_vals(ii), grid_x_vals(jj));
    end
    imagesc(grid_dHZ)
    pause(0.01)
end
%%
spl_dhurwitzZeta = @(a, x) interp2(grid_x_vals, grid_a_vals, grid_dHZ, x, a, 'spline');
lin_dhurwitzZeta = @(a, x) interp2(grid_x_vals, grid_a_vals, grid_dHZ, x, a, 'linear');
mak_dhurwitzZeta = @(a, x) interp2(grid_x_vals, grid_a_vals, grid_dHZ, x, a, 'makima');

%%
% 
spl_hurwitzZeta = @(a, x) interp2(grid_x_vals, grid_a_vals, grid_HZ, x, a, 'spline');
lin_hurwitzZeta = @(a, x) interp2(grid_x_vals, grid_a_vals, grid_HZ, x, a, 'linear');
mak_hurwitzZeta = @(a, x) interp2(grid_x_vals, grid_a_vals, grid_HZ, x, a, 'makima');

test_a_vals = [1.944 1.8717 2.319];
test_x_vals = [1 2 3 10 30 105 148 192 401 2e4];

tic
test_spl = spl_hurwitzZeta(test_a_vals', test_x_vals);
t_spl = toc;
tic
test_lin = lin_hurwitzZeta(test_a_vals', test_x_vals);
t_lin = toc;
tic
test_mak = mak_hurwitzZeta(test_a_vals', test_x_vals);
t_mak = toc;
tic
test_exact = bsxfun(@(ax,xy) hurwitzZeta(ax, xy), test_a_vals', test_x_vals);
t_exact = toc;
err_spl = test_spl - test_exact;
err_lin = test_lin - test_exact;
err_mak = test_mak - test_exact;
figure();
hold on
histogram(log10(abs(err_spl)))
histogram(log10(abs(err_lin)))
histogram(log10(abs(err_mak)))
legend({['spline (' num2str(t_exact/t_spl, '%1.0f') 'x)'], ...
    ['linear (' num2str(t_exact/t_lin, '%1.0f') 'x)'], ...
    ['mod-akima (' num2str(t_exact/t_mak, '%1.0f') 'x)']})
xlabel('log_{10} error')
%%

save(['avalanches_matlab_code/analysis_results/hurwitzLookup_' ...
    datestr(now, 'yyyymmdd')], 'grid_a_vals', 'grid_x_vals', 'grid_HZ', ... 
    'grid_dHZ', 'spl_dhurwitzZeta', 'lin_dhurwitzZeta', 'mak_dhurwitzZeta', ...
    'spl_hurwitzZeta', 'lin_hurwitzZeta', 'mak_hurwitzZeta');