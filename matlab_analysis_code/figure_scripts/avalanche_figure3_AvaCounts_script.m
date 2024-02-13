% Avalanche count plots


rep_list = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};
eta_list_sf1 = 1:10;
eps_list_sf1 = -2:-2:-14; % [ -6.0, -8.0, -10.0, -12.0, -14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_sf1, num_ava_sf1] = loadSimulationRunResults('run_f1smallfinesweep', rep_list);
N_samples_sf1 = double(cellfun(@(x) x.num_field_samples, x_info_sf1));
%% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/fields1smallfinesweep/';
x_reps_sf1 = cell(length(rep_list), 1);
for ii = 1:length(x_reps_sf1)
    rep_dir = [results_dir 'rep' rep_list{ii} '/' ];
    rep_fn = dir([rep_dir 'ava_decade_analysis*.mat']);
    [~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');
    x_reps_sf1{ii} = load([rep_dir rep_fn(date_ord(1)).name]);
end
sum_stats_sf1 = pullSummaryFitInfo(x_reps_sf1, eta_list_sf1, eps_list_sf1);

%% Make plot for this simulation, pull in some examples
% plotSizeDurPDFScaling(x_reps_T5, eta_T5_infTC, eps_T5_infTC, dyn_tau, 'B', h_Fig, h_sp5, f1f5_colors(2, :));
i_eps_ex2 = 3;
j_eta_ex2 = 7;
i_eps_ex1 = 2;
j_eta_ex1 = 4;
i_rep = 4;
ex_color2 = [0 0.8 0.2];
ex_color1 = [1 0.5 0];


%% Compute p_silence: N = 128

eta_list_ultra = linspace(0, 10, 101);
eps_list_ultra = linspace(0, -14, 141);
[p_sil_vs_ee, p_ava_vs_ee, ave_p_silence, ave_p_ava] = computeProbSilenceAvalanche(eta_list_ultra, eps_list_ultra, 128);
%% Compute p_silence and P_ava for 8-cell groups
num_draws = 2000;
[p_sil_vs_ee8, p_ava_vs_ee8, ave_p_silence8, ave_p_ava8] = computeProbSilenceAvalanche(eta_list_ultra, eps_list_ultra, 8, num_draws);
%% Compute p_silence and P_ava for 1024-cell groups
num_draws1k = 8;
[p_sil_vs_ee1024, p_ava_vs_ee1024, ave_p_silence1024, ave_p_ava1024] = computeProbSilenceAvalanche(eta_list_ultra, eps_list_ultra, 1024, num_draws1k);

%% interpolate
p_silence_ee = interp2(eta_list_ultra', eps_list_ultra, ave_p_silence, eta_list_sf1', eps_list_sf1);
p_ava_ee = interp2(eta_list_ultra', eps_list_ultra, ave_p_ava, eta_list_sf1', eps_list_sf1);

%%
for jj = 1:length(eta_list_sf1)
        nexttile
    imagesc(h_vals, eps_list_ultra, squeeze(p_sil_vs_ee(:, eta_list_ultra == eta_list_sf1(jj), :)), [0 1])
    title(['\eta = ' num2str(eta_list_sf1(jj))])
    yyaxis right
    plot(h_vals, normpdf(h_vals), 'linewidth', 1.5)
end

%%
[eta_list_overlap, eta_inds] = intersect(eta_list_ultra, eta_list_sf1);
figure()
nexttile
hold on
ph = plot(h_vals, squeeze(p_sil_vs_ee(eps_list_ultra == -6, eta_inds, :)));
assignColorsToLines(ph, parula(length(ph)))
legend(ph, num2str(eta_list_overlap(:), 'eta = %1.0f'))
plot(h_vals, normpdf(h_vals), 'k', 'LineWidth', 1.5)
title('\epsilon = -6')
nexttile
hold on
ph = plot(h_vals, squeeze(p_sil_vs_ee(1:10:end, eta_list_sf1 == 2, :)));
assignColorsToLines(ph, cool(length(ph)))
legend(ph, num2str(eps_list_ultra(1:10:end)', 'eps = %1.0f'))

ylabel('probability of silence at any time point')
yyaxis right
plot(h_vals, normpdf(h_vals), 'k', 'LineWidth', 1.5)
ylabel('p(h)')
xlabel('h')
%%
figure()
for i_h = 1:25
nexttile
imagesc(eta_list_sf1, eps_list_ultra,  squeeze(p_sil_vs_ee(:, :, 2*i_h)))
title(['h = ' num2str(h_vals(2*i_h))])
xlabel('\eta')
ylabel('\epsilon')
end
%% figure();
plot(h_vals, p_silence_ex1, h_vals, p_silence_ex2)

[scaled_Nava_SF1, scaled_rate_SF1, ave_rate_SF1, std_SF1, scaled_std_SF1] = scaleFRAvalancheCounts(x_info_sf1, num_ava_sf1, eta_list_sf1);
% [scaled_Nava_UF1, scaled_rate_UF1, ave_rate_UF1, std_UF1, scaled_std_UF1] = scaleFRAvalancheCounts(x_info_uf1, num_ava_uf1, eta_list_uf1);
%% MAKE THE FIGURE. 
% h_Fig = makeMyFigure(30, 15);
% tiledlayout(2, 3, 'TileSpacing', 'Compact', 'TileIndexing','columnmajor')
h_Fig = makeMyFigure(28, 14);
tiledlayout(2, 3, 'TileSpacing', 'Compact', 'TileIndexing','rowmajor')

%%%%%%%% plot of avalanche count again p_ava
nexttile
x = 1e7*p_ava_ee/1e6;
y = squeeze(mean(num_ava_sf1, 1))/1e6;
dy = squeeze(std(num_ava_sf1,[], 1))/sqrt(length(rep_list)-1)/1e6;

hold on
eh = errorbar(x', y', dy', '.', 'markersize', 8, 'Color', 0.5*[1 1 1], 'capsize', 0);
% assignColorsToLines(eh, parula(length(eh)));
% plot(ave_rate_SF1(:), num_ava_sf1(:), '.', 'markersize', 8, 'color', 0.5*[1 1 1])
% plot(ave_rate_UF1(:), num_ava_uf1(:), '.', 'color', [.5 0.5 1])
xlabel('P_{silence}(1-P_{silence}) N_T (\times 10^6)')
ylabel('Avalanche count (N_{ava})')
eqline
x_lims = xlim;
xlim([0 x_lims(2)])
y_lims = ylim;
ylim([0 y_lims(2)])
% axis square
set(gca, 'color', 'none', 'fontsize', 12)
plot(x(i_eps_ex1, j_eta_ex1), ...
    num_ava_sf1(i_rep, i_eps_ex1, j_eta_ex1), '*' , 'color', ex_color1,'linewidth', 1.5)
plot(x(i_eps_ex2, j_eta_ex2), ...
    num_ava_sf1(i_rep, i_eps_ex2, j_eta_ex2), '*', 'color', ex_color2, 'linewidth', 1.5)


text(-0.4, 1.15, 'A', 'Units', 'normalized', 'FontSize',16)



sh = nexttile;
imagesc(h_vals, -eps_list_ultra, squeeze(p_ava_vs_ee(:, eta_list_ultra == 2, :)), [0 0.25])
colormap(sh, flipud(gray));
xlabel('h')
ylabel('\epsilon (bias)')
set(gca, 'ytick', 2:6:14, 'ydir', 'reverse', 'fontsize', 12)
ch = colorbar;
ch.Label.String = 'P_{silence}(1 - P_{silence})';
ch.Ticks = [0 0.25];
% ch.Location = 'northoutside';
% ch.Position(3) = ch.Position(3)+.01;
% % yyaxis right
% % % nexttile
% % bar(h_vals, normpdf(h_vals), 'BarWidth',1, 'EdgeColor','none', 'FaceAlpha', 0.5)
% % xlabel('h')
% % ylabel('p(h)')
title(['\eta = ' num2str(eta_list_sf1(2))])
% set(gca, 'box', 'off', 'color', 'none', 'fontsize', 12)
text(-0.4, 1.15, 'B', 'Units', 'normalized', 'FontSize',16)

% % ava_prob_dot_h_prob = zeros([size(squeeze(p_ava_vs_ee(:, eta_list_ultra == 2, :))) 3]);
% % h_prob = sqrt(2*pi)*repmat(normpdf(h_vals), [size(ava_prob_dot_h_prob, 1) 1 3]);
% % all_shades = lines(3);
% % shade_p_h = 1 - all_shades(1, :);
% % h_prob_cmap(:, :, 1) = shade_p_h(1)*h_prob(:, :, 1);
% % h_prob_cmap(:, :, 2) = shade_p_h(2)*h_prob(:, :, 1);
% % h_prob_cmap(:, :, 3) = shade_p_h(3)*h_prob(:, :, 1);
% % 
% % 
% % p_ava_e2 = squeeze(p_ava_vs_ee(:, eta_list_ultra == 2, :))/max(p_ava_e2(:));
% % shade_p_ava = 1 - 0*all_shades(2, :);
% % p_ava_cmap(:, :, 1) = shade_p_ava(1)*p_ava_e2;
% % p_ava_cmap(:, :, 2) = shade_p_ava(2)*p_ava_e2;
% % p_ava_cmap(:, :, 3) = shade_p_ava(3)*p_ava_e2;
% % 
% % % imagesc(h_vals, -eps_list_ultra, squeeze(p_ava_vs_ee(:, 2, :)), [0 0.25])
% % % imagesc(h_vals, -eps_list_ultra, (1-p_ava_cmap).*(h_prob_cmap + 1)/2)
% % imagesc(h_vals, -eps_list_ultra, (2- h_prob_cmap - p_ava_cmap) /2)
% % 

%%%%%%% same as the first, but showing eta-dependence
sh = nexttile;
imagesc(h_vals, eta_list_sf1, squeeze(p_ava_vs_ee(eps_list_ultra == -8, :, :)), [0 0.25])
colormap(sh, flipud(gray));
xlabel('h')
ylabel('\eta (gain)')
set(gca, 'ytick', 2:2:10, 'ydir', 'normal', 'fontsize', 12)
title('\epsilon = 8')

ch = colorbar;
ch.Label.String = 'P_{silence}(1 - P_{silence})';
ch.Ticks = [0 0.25];
text(-0.4, 1.15, 'C', 'Units', 'normalized', 'FontSize',16)

% ch = colorbar;
% ch.Label.String = ['P_{silence}(1 - P_{silence}), \epsilon = ' num2str(-8)];
% ch.Ticks = [0 0.25];
% ch.Location = 'northoutside';
% ch.Label.Rotation = -90;
% % yyaxis right
% % % plot(h_vals, normpdf(h_vals), '-.', 'linewidth', 1.5)
% % bar(h_vals, normpdf(h_vals), 'BarWidth',1, 'EdgeColor','none', 'FaceAlpha', 0.5)

% % ylabel('p(h)')




% ch.Label.Rotation = -90;
% yyaxis right
% plot(h_vals, normpdf(h_vals), '-.', 'linewidth', 1.5)





% nexttile
% size_min_all = sum_stats_sf1.tau_x_min;
% ave_min_size = squeeze(exp(mean(log(size_min_all), 3, 'omitnan')));
% plot(eps_list_sf1, ave_min_size, '.', 'markersize', 10, 'color', 0.5*[1 1 1])

% ylabel('Power law cutoff S_{min}')
% xlabel('\epsilon')
% plot(h_vals, squeeze(p_ava_vs_ee(eps_list_ultra == -8, 2:4:10, :)))
% imagesc(h_vals, eps_list_ultra, squeeze(p_sil_vs_ee(:, 2, :).*(1 - p_sil_vs_ee(:, 2, :))), [0 0.25])
% colormap(sh, flipud(gray));
% xlabel('h')
% ylabel('\epsilon')
% set(gca, 'ytick', -14:6:-2, 'ydir', 'normal', 'fontsize', 12)
% ch = colorbar;
% ch.Label.String = ['P_{silence}(1 - P_{silence}), \eta = ' num2str(eta_list_sf1(2))];
% ch.Ticks = [0 0.25];
% ch.Location = 'northoutside';
% % ch.Label.Rotation = -90;
% yyaxis right
% % plot(h_vals, normpdf(h_vals), '-.', 'linewidth', 1.5)
% bar(h_vals, normpdf(h_vals), 'BarWidth',1, 'EdgeColor','none', 'FaceAlpha', 0.5)

%figure(), ph = plot(p_silence_ee, squeeze(mean(scaled_Nava_SF1, 1)), '.', 'markersize', 8); assignColorsToLines(ph, parula(length(ph)))


% xlabel('Psilence')
% ylabel('\eta Nava')
% % xlim([0 0.25])
% set(gca, 'color', 'none', 'fontsize', 12)
% 
% sh = nexttile;
% size_min_all = sum_stats_sf1.tau_x_min;
% ave_min_size = squeeze(exp(mean(log(size_min_all), 3, 'omitnan')));
% imagesc(eta_list_sf1, eps_list_sf1, log10(ave_min_size), log10([1 100]))
% hold on
% rms_rate = squeeze(sqrt(mean(scaled_rate_SF1, 1)));
% contour(eta_list_sf1, eps_list_ultra, ave_p_silence, 0.05:0.1:1, ...
%     'LineColor', [0 0 0], 'ShowText','on')
% nan_parula = [0.8 0.8 0.8; parula];
% colormap(sh, nan_parula);
% ch = colorbar;
% ch.Ticks = log10([1 10 100]);
% ch.TickLabels = {'1', '10', '100'};
% xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
% set(gca, 'YDir', 'normal', 'FontSize', 12)
% axis square





% nexttile

% hold on
% plot(sqrt(scaled_std_SF1(:)), scaled_Nava_SF1(:), '.', 'color', 0.5*[1 1 1])
% plot(sqrt(scaled_std_UF1(:)), scaled_Nava_UF1(:), '.', 'color', [.5 .5 1])




% 
% yyaxis right
nexttile
plot(-eps_list_ultra, p_ava_vs_ee(:, 1, h_vals == 0), 'k', 'linewidth', 1.5)
x_int1 = -eps_list_sf1([i_eps_ex1 i_eps_ex2]);
y_int1 = interp1(-eps_list_ultra', p_ava_vs_ee(:, 1, h_vals == 0), x_int1);
hold on
plot(x_int1(1), y_int1(1), 'x', 'color', ex_color1,'linewidth', 3, 'markersize', 10);
plot(x_int1(2), y_int1(2), 'x', 'color', ex_color2, 'linewidth', 3, 'markersize', 10);
xlabel('\epsilon (bias)')
ylabel('P_S(1-P_S) at h = 0', 'Color', [0 0 0])
set(gca, 'YColor', [0 0 0], 'color', 'none')
text(-0.4, 1.2, 'D', 'Units', 'normalized', 'FontSize',16)

h_ex1 = nexttile;
sp_list = [h_ex1 h_ex1 h_ex1 h_ex1];
make_subplots = [true false false false];
plotSizeDurPDFScaling(x_reps_sf1{i_rep}, eta_list_sf1(j_eta_ex1), eps_list_sf1(i_eps_ex1), [], [], h_Fig, sp_list, ex_color1, make_subplots);
% axis normal
text(0.7, 0.8, '\epsilon < \epsilon_0', 'Units','normalized', 'FontSize',12)
text(-0.4, 1.2, 'E', 'Units', 'normalized', 'FontSize',16)


h_ex2 = nexttile;
sp_list(1) = h_ex2;
plotSizeDurPDFScaling(x_reps_sf1{i_rep}, eta_list_sf1(j_eta_ex2), eps_list_sf1(i_eps_ex2), [], [], h_Fig, sp_list, ex_color2, make_subplots);
% axis normal
text(0.7, 0.8, '\epsilon > \epsilon_0', 'Units','normalized', 'FontSize',12)
text(-0.2, 1.2, 'F', 'Units', 'normalized', 'FontSize',16)


set(h_ex1, 'FontSize', 13);
set(h_ex2, 'FontSize', 13);
%%
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/fig3_ava_counts_' datestr(now, 'yyyymmdd')])
print(gcf, '-dpng', 'avalanches_matlab_code/plots/paper_figures/fig3_ava_counts.png')
%% Supporting: show size cutoff figures

%%%%%%%%%%%% size cutoffs
% sh= nexttile;
% imagesc(-eps_list_sf1, eta_list_sf1, log10(ave_min_size)', [0 2]);
% ylabel('\eta')
% 
% ch = colorbar;
% ch.Ticks = log10([1 10 100]);
% ch.TickLabels = {'1', '10', '100'};
% ch.Label.String = 'Power Law Size Cutoff';
% ch.FontSize = 11;
% % ch.Location = 'southoutside';
% % ch.Position(3) = ch.Position(3) - .1 ;
% cmap_size = flipud(pink);
% cmap_size = cmap_size(1:end - 20, :);
% % cmap_size(:, 1) = 1;
% cmap_size = [0.5*[1 1 1]; cmap_size; 0.5*[1 1 1]];
% colormap(sh, cmap_size);


%% Make the Information Figure: load ave_MLE_Info

x_info = load('average_info_calculations_20221214.mat');
eta_vals = x_info.eta_vals;
eps_vals = x_info.eps_vals;
ave_MLE_Info128 = x_info.ave_MLE_Info128;
ave_MLE_Info1024 = x_info.ave_MLE_Info1024;



%%
fracNaN_MLE_Info128 = mean((x_info.mle_I_vals_128 < 0), 3);

exact_Info = x_info.all_I_vals;
%% 
ave_Info8 = squeeze(mean(exact_Info(:, :, end, :), 4));
% 
%% Information plot with P_ava contours, for N = 1024 and N = 8 (suppl)
% note : must run ave_p_ava calculations for N = 1024 still , currently
% showing N = 128 calculations (12/23/22)
makeMyFigure(28, 24);
tiledlayout(4, 4, 'Padding','tight')
x_Loc_L = -0.6;
y_Loc_L = 1.3;
x_Loc_L0 = -1;
y_Loc_L0 = 0;

nexttile([2 2])
imagesc(eta_vals, eps_vals, ave_MLE_Info1024, [0 3.5]);
ch = colorbar;
ch.Label.String = 'I(h^*, h) for N = 1024 cells (MLE)';
hold on
% contour(eta_vals, eps_vals, ave_MLE_Info128, ...
%     0.5:0.5:max(ave_MLE_Info128(:)), 'color', [1 0.5 0.5], 'linewidth', 2)
ylim([0 14])
set(gca, 'fontsize', 13)
% nexttile
% imagesc(eta_list_ultra, -eps_list_ultra, ave_p_ava)

% smooth p_ava for contours
smooth_p_ava1024 = imfilter(ave_p_ava1024, fspecial('gaussian', 10, 2), 'replicate');
% compute max along epsilon at each value of eta
[~, inds_max1024] = max(smooth_p_ava1024, [], 1);
eps_vals_max = eps_list_ultra(inds_max1024);

[~, inds_maxI] = max(ave_MLE_Info1024, [], 1);
eps_vals_maxI = eps_vals(inds_maxI);


c_lines = [1e-3 3e-2 1e-1 2e-1];
% hold on
contour(eta_list_ultra, -eps_list_ultra, smooth_p_ava1024, c_lines, 'color', [0 0 0],...
    'linewidth', 2, 'ShowText', 'on')
plot(eta_list_ultra, -eps_vals_max, 'm', 'linewidth', 1.5)
plot(eta_vals, eps_vals_maxI, 'w', 'linewidth', 1.5)

% just plot where the local maxes are by d^2/deps^2 
legend({'P_{ava} contours', 'P_{ava} ridge', 'I_{max}(\eta)'},...
    'Color', 0.8*[ 1 1 1], 'Location','Southwest')
xlabel('\eta (input scaling)')
ylabel('\epsilon (bias)')
title('Ensembles of N = 1024 cells')

text(x_Loc_L0, y_Loc_L0, 'A', 'Fontsize', 16)
% axis square
% 
% hold on
% eps0 = log(2^(1/128)-1);
% plot(eta_vals, 0*eta_vals - eps0, 'r')

% plot exact for 8-cells
nexttile([2 2])
imagesc(eta_vals, eps_vals, ave_Info8, [0 1.2]);
% imagesc(eta_vals, eps_vals, ave_Info8./(1e-3+r_i));
hold on
% contour(eta_vals, eps_vals, ave_Info8, ...
%     0.6:0.2:max(ave_Info8(:)), 'color', [1 0.5 0.5], 'linewidth', 2)
ylim([0 14])
ch = colorbar;
ch.Label.String = 'I({s_i}, h) for ensembles of N = 8 cells (exact)';

% interpolate value of MI
interp_MLE_Info8 = interp2(eta_vals, eps_vals', ave_Info8, eta_list_ultra, -eps_list_ultra', 'spline');
[~, inds_maxI8] = max(interp_MLE_Info8, [], 1);
eps_vals_maxI8 = -eps_list_ultra(inds_maxI8);


% [~, inds_maxI8] = max(ave_Info8, [], 1);
% eps_vals_maxI8 = eps_vals(inds_maxI8);
set(gca, 'fontsize', 13)

% nexttile
% imagesc(eta_list_ultra, -eps_list_ultra, ave_p_ava)
% smooth p_ava for contours
smooth_p_ava8 = imfilter(ave_p_ava8, fspecial('gaussian', 10, 2), 'replicate');
% compute max along epsilon at each value of eta
[~, inds_max8] = max(smooth_p_ava8, [], 1);
eps_vals8_max = eps_list_ultra(inds_max8);


c_lines = [1e-3 3e-2 1e-1 2e-1];
% hold on
contour(eta_list_ultra, -eps_list_ultra, smooth_p_ava8, c_lines, 'color', [0 0 0],...
    'linewidth', 2, 'ShowText', 'on')
plot(eta_list_ultra, -eps_vals8_max, 'm', 'linewidth', 1.5)
plot(eta_list_ultra, eps_vals_maxI8, 'w', 'linewidth', 1.5)

% just plot where the local maxes are by d^2/deps^2 
legend({'P_{ava} contours', 'P_{ava} ridge', 'I_{max}(\eta)'}, ...
    'Color', 0.8*[ 1 1 1], 'Location','Southwest')
xlabel('\eta (input scaling)')
ylabel('\epsilon (bias)')
title('Ensembles of N = 8 cells')
text(x_Loc_L0, y_Loc_L0, 'B', 'Fontsize', 16)




eta_list = [2 5 7 9];
label_list = {'C', 'D', 'E', 'F'};
for ii = 1:length(eta_list)
nexttile
select_eta = eta_list(ii);
plot(eps_vals, ave_MLE_Info1024(:, eta_vals == select_eta), 'linewidth', 2)
xlim([0 14])
ylabel('I(h, h*), N = 128')
y_lims = ylim;
ylim([0 y_lims(2)])
yyaxis right
plot(-eps_list_ultra, ave_p_ava1024(:, eta_list_ultra == select_eta), 'linewidth', 2)
title(['\eta = ' num2str(select_eta)])
xlabel('\epsilon (bias)')
ylabel('P_{ava}')
axis square
set(gca, 'color', 'none', 'fontsize', 12)
text(x_Loc_L, y_Loc_L, label_list{ii}, 'units', 'normalized', 'fontsize', 16)
end

label_list = {'G', 'H', 'I', 'J'};

eta_list = [2 5 7 9];
for ii = 1:length(eta_list)
nexttile
select_eta = eta_list(ii);
plot(eps_vals, ave_Info8(:, eta_vals == select_eta), 'linewidth', 2)
xlim([0 14])
ylabel('I(h, {s_i}), N = 8')
y_lims = ylim;
ylim([0 y_lims(2)])

yyaxis right
plot(-eps_list_ultra, smooth_p_ava8(:, eta_list_ultra == select_eta), 'linewidth', 2)
title(['\eta = ' num2str(select_eta)])
xlabel('\epsilon (bias)')
ylabel('P_{ava}')
axis square
set(gca, 'color', 'none', 'fontsize', 12)
text(x_Loc_L, y_Loc_L, label_list{ii}, 'Units','normalized', 'fontsize', 16)
end

print(gcf, '-dpng', ['avalanches_matlab_code/plots/paper_figures/suppl_info_avalanches_' datestr(now, 'yyyymmdd') '.png'])
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/suppl_info_avalanches_' datestr(now, 'yyyymmdd') '.pdf'])

%% Simplified Info contours - N = 128 only

makeMyFigure(28, 18);
tiledlayout(5, 6)
x_Loc_L = -0.4;
y_Loc_L = 1.1;
x_Loc_L0 = -1.5;
y_Loc_L0 = 0;

nexttile([3 3])
imagesc(eta_list_ultra, -eps_list_ultra, ave_p_ava, [0 0.25]);
ch = colorbar;
ch.Label.String = 'Probability of Avalanche';
hold on
% contour(eta_vals, eps_vals, ave_MLE_Info128, ...
%     0.5:0.5:max(ave_MLE_Info128(:)), 'color', [1 0.5 0.5], 'linewidth', 2)
ylim([0 14])
set(gca, 'fontsize', 13)
% nexttile
% imagesc(eta_list_ultra, -eps_list_ultra, ave_p_ava)

% smooth p_ava for contours
smooth_p_ava = imfilter(ave_p_ava, fspecial('gaussian', 10, 2), 'replicate');
% compute max along epsilon at each value of eta
[~, inds_max] = max(smooth_p_ava, [], 1);
eps_vals_max = eps_list_ultra(inds_max);


% interpolate value of MI
interp_MLE_Info128 = interp2(xSCI.eta_vals, xSCI.eps_vals', info2_vals, eta_list_ultra, -eps_list_ultra', 'spline');
% interp_MLE_Info128 = interp2(eta_vals, eps_vals', ave_MLE_Info128, eta_list_ultra, -eps_list_ultra', 'spline');
[~, inds_maxI] = max(interp_MLE_Info128, [], 1);
eps_vals_maxI = -eps_list_ultra(inds_maxI);

% interpolate validity region (fraction of info cacluations at eache ps,
% eta value that generated negative values -> means calculation failed)
interp_Valid_Info128 = interp2(eta_vals, eps_vals', fracNaN_MLE_Info128, eta_list_ultra, -eps_list_ultra', 'linear');


c_lines = [1e-3 2e-2 4e-2 6e-2 8e-2 1e-1 2e-1];
% hold on
contour(eta_list_ultra, -eps_list_ultra, smooth_p_ava, c_lines, 'color', [0 0 0],...
    'linewidth', 2, 'ShowText', 'on');
plot(eta_list_ultra, -log(2^(1/128)-1) + 0*eta_list_ultra, '-.',...
    'color', [1 0.5 1], 'linewidth', 2)

plot(eta_list_ultra, -eps_vals_max, 'm', 'linewidth', 1.5)
% plot(eta_vals, eps_vals_maxI, 'w', 'linewidth', 1.5)
xlabel('\eta (input scaling)')
ylabel('\epsilon (bias)')
text(x_Loc_L0, y_Loc_L0, 'A', 'Fontsize', 16)


nexttile([3 3])
imagesc(eta_list_ultra, -eps_list_ultra, interp_MLE_Info128, [0 2.5])
% imagesc(eta_vals, eps_vals, ave_MLE_Info128, [0 2.5]);
ch = colorbar;
ch.Label.String = 'I(h^*, h) for N = 128 cells (MLE)';
hold on
% contour(eta_vals, eps_vals, ave_MLE_Info128, ...
%     0.5:0.5:max(ave_MLE_Info128(:)), 'color', [0 0 0], 'linewidth', 2)
[cntr_mat, ph4] = contour(eta_list_ultra, -eps_list_ultra, smooth_p_ava, [0.001 0.001], 'color', [0 0 0], ...
    'linestyle', '-.', 'linewidth', 2);

[~, ph3] = contour(eta_list_ultra, -eps_list_ultra, interp_Valid_Info128, [0.1 0.1], 'color', [0.5 0.5 0.5], ...
    'linewidth', 2);
ph5 = plot(eta_list_ultra, -log(2^(1/128)-1) + 0*eta_list_ultra, '-.',...
    'color', [1 0.5 1], 'linewidth', 2);
set(gca, 'fontsize', 13)
ph1 = plot(eta_list_ultra, -eps_vals_max, 'm', 'linewidth', 1.5);
% ph2 = plot(eta_list_ultra, eps_vals_maxI, 'k', 'linewidth', 1.5);

% just plot where the local maxes are by d^2/deps^2 
lh = legend([ph3, ph4, ph5, ph1], {'I_{MLE} Valid','P_{ava} = 0.001', ... 
    '\epsilon_0', '\epsilon^*(\eta)'},...
    'Color', 0.8*[ 1 1 1], 'Location','Southeast');

ylim([0 14])
xlabel('\eta (input scaling)')
ylabel('\epsilon (bias)')
% title('Ensembles of N = 128 cells')
text(x_Loc_L0, y_Loc_L0, 'B', 'Fontsize', 16)

% axis square
% 
% hold on
% eps0 = log(2^(1/128)-1);
% plot(eta_vals, 0*eta_vals - eps0, 'r')
eta_list = [2 5 9];
label_list = {'C', 'D', 'E'};
for ii = 1:length(eta_list)
nexttile([2 2])
select_eta = eta_list(ii);
plot(eps_vals, ave_MLE_Info128(:, eta_vals == select_eta), 'k', 'linewidth', 2)
xlim([0 14])
ylabel('I(h, h*), N = 128')
y_lims = ylim;
y_lims(1) = 0;
ylim(y_lims)
hold on
plot(-eps_vals_max(eta_list_ultra == select_eta)*[1 1], y_lims, 'm', 'linewidth', 2)

p_ava_at_eta = ave_p_ava(:, eta_list_ultra == select_eta);
p0001 = find(p_ava_at_eta > 0.001, 1);
plot(-eps_list_ultra(p0001)*[1 1], y_lims, 'k-.', 'linewidth', 2);
% yyaxis right
% plot(-eps_list_ultra, ave_p_ava(:, eta_list_ultra == select_eta), 'linewidth', 2)
title(['\eta = ' num2str(select_eta)])
xlabel('\epsilon (bias)')
% ylabel('P_{ava}')
axis square
set(gca, 'color', 'none', 'fontsize', 12)
text(x_Loc_L, y_Loc_L, label_list{ii}, 'units', 'normalized', 'fontsize', 16)
end


print(gcf, '-dpng', ['avalanches_matlab_code/plots/paper_figures/fig5_info_avalanches_' datestr(now, 'yyyymmdd') '.png'])
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/fig5_info_avalanches_' datestr(now, 'yyyymmdd') '.pdf'])


%% Below this point: aggregating other simulation runs
rep_list = {'A', 'B', 'C', 'D', 'E'};
eta_list_uf1 = [4.0, 6.0, 8.0];
eps_list_uf1 = [ -6.0, -8.0, -10.0, -12.0, -14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_uf1, num_ava_uf1] = loadSimulationRunResults('run_f1ultrafinesweep', rep_list);

%%

[x_info_f1, num_ava_f1] = loadSimulationRunResults('run_f1finesweep');

%% time sweeps
[x_info_T1, num_ava_T1] = loadSimulationRunResults('run_f1_e1204_timesweep');
[x_info_T5, num_ava_T5] = loadSimulationRunResults('run_f5_e1204_timesweep');
%% goal - plot eta*Nava vs. 1/eta <fr> 

%% do the SF results still match? 

ave_rate_SF = cellfun(@(x) mean(x.cell_FR), x_info_sf1);
eta_arr_SF = repmat(eta_list_sf1(:), [1 size(ave_rate_SF, [2 1])]);
eta_arr_SF = permute(eta_arr_SF, [3 2 1]);
scaled_Nava = num_ava_sf1.*eta_arr_SF./N_samples_sf1;
scaled_rate = ave_rate_SF./eta_arr_SF;

rate_bins = linspace(0, 5, 51);
count_bins = linspace(0, 1e7, 51);
%% SF 
figure()
nexttile
plot(ave_rate_SF(:),num_ava_sf1(:), '.' )
nexttile
plot(scaled_rate(:), scaled_Nava(:), '.')
nexttile
h_cts = histcounts2(scaled_Nava(:), scaled_rate(:), count_bins, rate_bins);
imagesc(rate_bins, count_bins, h_cts)
% plot(scaled_rate(:), (scaled_Nava(:)), '.')
set(gca, 'ydir', 'normal')
nexttile
x = reshape(scaled_rate, [numel(scaled_rate(:, :, 1)) size(scaled_rate, 3)]);
y = reshape(scaled_Nava, [numel(scaled_Nava(:, :, 1)) size(scaled_Nava, 3)]);
ph = plot(x, y, '.');
assignColorsToLines(ph, parula(length(ph)))

%% now for the ... UF1 results
[scaled_Nava_UF1, scaled_rate_UF1, ave_rate_UF1] = scaleFRAvalancheCounts(x_info_uf1, num_ava_uf1, eta_list_uf1);
N_samples_uf1 = nan(size(scaled_Nava_UF1));
N_samples_uf1(scaled_Nava_UF1 > 0) = cellfun(@(x) x.num_field_samples, x_info_uf1((scaled_Nava_UF1 > 0)));
scaled_Nava_UF1 = scaled_Nava_UF1./N_samples_uf1;
%% t1 and t5
eta_T1 = x_info_T1{1}.eta;
[scaled_Nava_T1, scaled_rate_T1, ave_rate_T1] = scaleFRAvalancheCounts(x_info_T1, num_ava_T1, eta_T1);
N_samples_T1 = cellfun(@(x) double(x.num_field_samples), x_info_T1);

scaled_Nava_T1 = scaled_Nava_T1./N_samples_T1;
%%
figure()
nexttile
plot(ave_rate_UF1(:),num_ava_uf1(:), '.' )
nexttile
plot(scaled_rate_UF1(:), scaled_Nava_UF1(:), '.')
nexttile
h_cts = histcounts2(scaled_Nava_UF1(:), scaled_rate_UF1(:), count_bins, rate_bins);
imagesc(rate_bins, count_bins, h_cts)
% plot(scaled_rate(:), (scaled_Nava(:)), '.')
set(gca, 'ydir', 'normal')

nexttile
x = reshape(scaled_rate_UF1, [numel(scaled_rate_UF1(:, :, 1)) size(scaled_rate_UF1, 3)]);
y = reshape(scaled_Nava_UF1, [numel(scaled_Nava_UF1(:, :, 1)) size(scaled_Nava_UF1, 3)]);
ph = plot(x, log10(y), '.');
assignColorsToLines(ph, parula(length(ph)))

%%
figure()
nexttile
hold on
plot(ave_rate_SF(:),num_ava_sf1(:)./N_samples_sf1(:), '.', 'markersize', 12 )

plot(ave_rate_UF1(:),num_ava_uf1(:)./N_samples_uf1(:), '.' , 'markersize', 12 )
plot(ave_rate_T1(:), num_ava_T1(:)./N_samples_T1(:), '.', 'markersize', 12)
xlabel("<s_i>")
ylabel('N_{ava}')
set(gca, 'fontsize', 12, 'color', 'none')
legend({'N = 128', 'N = 1024', 'Tsweep, N = 1024'})

xlim
nexttile
hold on


plot(scaled_rate(:), scaled_Nava(:), '.', 'markersize', 12 )
plot(scaled_rate_UF1(:), scaled_Nava_UF1(:), '.', 'markersize', 12 )
plot(scaled_rate_T1(:), scaled_Nava_T1(:), '.', 'markersize', 12)

ylabel('eta * N_{ava}')
xlabel('<s_i> / eta')
set(gca, 'fontsize', 12, 'color', 'none')
legend({'N = 128', 'N = 1024'})

%% at say value of eta, eps, how do firing rate distributions differ? 
sf_eta = 4;
sf_eps = 3; 
uf_eta = 1;
uf_eps = 1;

fr_sf1_ij = cell2mat(cellfun(@(x) sort(x.cell_FR), x_info_sf1(:, sf_eps, sf_eta), 'UniformOutput',false)');
fr_uf1_ij = cell2mat(cellfun(@(x) sort(x.cell_FR), x_info_uf1(:, uf_eps, uf_eta), 'UniformOutput',false)');
figure()
hold on
plot((1:128)/128, log(fr_sf1_ij), 'k')
plot((1:1024)/1024, log(fr_uf1_ij), 'b')

%% is teh percentile at which <FR> = 0.1 (or whatever arb. threshold) determinative
% of avalanche count? 
thr = 0.1;

fr_pctl_thr_sf1 = cellfun(@(x) mean(x.cell_FR > thr), x_info_sf1);
uf_has_entry = cellfun(@(x) ~isempty(x), x_info_uf1);
fr_pctl_thr_uf1 = nan(size(fr_pctl_thr_uf1));
fr_pctl_thr_uf1(uf_has_entry) = cellfun(@(x) mean(x.cell_FR > thr), x_info_uf1(uf_has_entry));
figure()
hold on
plot(fr_pctl_thr_uf1(:), num_ava_uf1(:)./N_samples_uf1(:), 'o')
plot(fr_pctl_thr_sf1(:), num_ava_sf1(:)./N_samples_sf1(:), 'o')
%%
figure()
for ii = 1:5
    nexttile
    imagesc(squeeze(num_ava_sf1(ii, :, :)))
    
    nexttile
    imagesc(squeeze(fr_pctl_thr_sf1(ii, :, :)), [0 1])

    nexttile
    imagesc(squeeze(num_ava_uf1(ii, :, :)))
    nexttile
    imagesc(squeeze(fr_pctl_thr_uf1(ii, :, :)), [0 1])
end

%% variance of firing rate == eta?
[scaled_Nava_SF1, scaled_rate_SF1, ave_rate_SF1, std_SF1, scaled_std_SF1] = scaleFRAvalancheCounts(x_info_sf1, num_ava_sf1, eta_list_sf1);

% std_fr_sf1 = cellfun(@(x) std(x.cell_FR), x_info_sf1);
% ave_fr_sf1 = cellfun(@(x) mean(x.cell_FR), x_info_sf1);

figure(),
nexttile
plot(scaled_rate(:), (scaled_Nava(:)), '.')

nexttile
plot(scaled_std_SF1(:), scaled_Nava(:), '.')

ave_std_SF1 = squeeze(sqrt(mean(scaled_std_SF1.^2, 1)));
mult_ave_Nava = squeeze(exp(mean(log(num_ava_sf1), 1)));

mult_ave_scaledNava = squeeze(exp(mean(log(scaled_Nava_SF1), 1)));
nexttile
ph = plot(ave_std_SF1, mult_ave_Nava, '.');
assignColorsToLines(ph, parula(length(ph)))

nexttile
ph = plot(ave_std_SF1, mult_ave_scaledNava, '.');
assignColorsToLines(ph, parula(length(ph)))
%% 

%%

function [scaled_Nava_UF1, scaled_rate_UF1, ave_rate_UF1, std_rate_UF1, scaled_rate_std] = scaleFRAvalancheCounts(x_info_uf1, num_ava_uf1, eta_list_uf1)

has_entry_UF1 = cellfun(@(x) ~isempty(x), x_info_uf1);
ave_rate_UF1 = nan(size(has_entry_UF1));
ave_rate_UF1(has_entry_UF1) = cellfun(@(x) mean(x.cell_FR), x_info_uf1(has_entry_UF1));

% same, compute std(rate)
std_rate_UF1 = nan(size(has_entry_UF1));
std_rate_UF1(has_entry_UF1) = cellfun(@(x) std(x.cell_FR), x_info_uf1(has_entry_UF1));

% put the eta values in the correct array size
eta_arr_UF1 = repmat(eta_list_uf1(:), [1 size(ave_rate_UF1, [2 1])]);
eta_arr_UF1 = permute(eta_arr_UF1, [3 2 1]);
scaled_Nava_UF1 = num_ava_uf1.*eta_arr_UF1;
scaled_rate_UF1 = ave_rate_UF1./eta_arr_UF1;
scaled_rate_std = std_rate_UF1./eta_arr_UF1;

end

