% Script to plot basic parameters from each run (like the coupling
% matrices)

rep_list = {'A', 'B', 'C', 'D', 'E'};

x_info = cell(length(rep_list),1 );
% file_name_rep =  @(str) ['avalanches_code/data/fields5ultrafinesweep/sweep_ultrafine_rep' str ...
%     '_av_stim5e-10.0et6.0ph1.0p1.0tstatfr_stats.mat'];

file_name_rep =  @(str) ['avalanches_code/data/fields5timesweep/time' str ...
    '_stim5e-12.0et4.0ph1.0p1.0t0.1fr_stats.mat'];


for ii = 1:length(rep_list)
    x_info{ii} = load(file_name_rep(rep_list{ii}));
end

%%
J_mats = cellfun(@(x) x.J, x_info, 'UniformOutput', false);
[~, ord1] = sort(J_mats{1}(1, :));
figure()
for ii = 1:length(J_mats)
    subplot(length(J_mats), 2, 2*ii - 1)
    imagesc(J_mats{ii}(:, ord1))
    
    subplot(length(J_mats), 2, 2*ii )
    histogram(J_mats{ii})
end

%% Firing rate distribution
eff_bin_size = 0.003;    % bin size in seconds
fr_mats = cellfun(@(x) x.cell_FR/eff_bin_size, x_info, 'UniformOutput',false);
figure()
tiledlayout(3, 2)
for ii = 1:length(fr_mats)
    nexttile
    histogram(fr_mats{ii})
end

%% LOAD: 5-field Ultrafine sweep with replicates (5USR): plotting firing rate, avalanche counts, and field values

rep_list = {'A', 'B', 'C', 'D', 'E'};
eta_list_uf5 = [4.0, 6.0, 8.0];
eps_list_uf5 = [ -6.0, -7.0, -8.0, -9.0, -10.0, -12.0] ;%-14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_uf5, num_ava_uf5] = loadSimulationRunResults('run_f5ultrafinesweep', rep_list);

% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/fields5ultrafinesweep/';
x_reps_uf5 = cell(length(rep_list), 1);
for ii = 1:length(x_reps_uf5)
    rep_dir = [results_dir 'rep' rep_list{ii} '/' ];
    rep_fn = dir([rep_dir '*.mat']);
    x_reps_uf5{ii} = load([rep_dir rep_fn(1).name]);
end
sum_stats_uf5 = pullSummaryFitInfo(x_reps_uf5, eta_list_uf5, eps_list_uf5);

%% LOAD: 1-field Ultrafine sweep with replicates (5USR): plotting firing rate, avalanche counts, and field values

rep_list = {'A', 'B', 'C', 'D', 'E'};
eta_list_uf1 = [4.0, 6.0, 8.0];
eps_list_uf1 = [ -6.0, -8.0, -10.0, -12.0, -14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_uf1, num_ava_uf1] = loadSimulationRunResults('run_f1ultrafinesweep', rep_list);

% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/fields1ultrafinesweep/';
x_reps_uf1 = cell(length(rep_list), 1);
for ii = 1:length(x_reps_uf1)
    rep_dir = [results_dir 'rep' rep_list{ii} '/' ];
    rep_fn = dir([rep_dir 'ava_decade_analysis*.mat']);
    [~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');
    x_reps_uf1{ii} = load([rep_dir rep_fn(date_ord(1)).name]);
end
sum_stats_uf1 = pullSummaryFitInfo(x_reps_uf1, eta_list_uf1, eps_list_uf1);


%% LOAD: 1-field time sweep with replicates (T1): plotting firing rate, avalanche counts, and field values

% rep_list = {'A', 'B', 'C', 'D', 'E'};
% eta_list_uf5 = [4.0, 6.0, 8.0];
% eps_list_uf5 = [ -6.0, -7.0, -8.0, -9.0, -10.0, -12.0] ;%-14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_T1, num_ava_T1] = loadSimulationRunResults('run_f1_e1204_timesweep');
% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/timesweep_stim1_e12et04/';
rep_fn = dir([results_dir 'ava_decade_analysis*.mat']);
[~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');

x_reps_T1 = load([results_dir rep_fn(date_ord(1)).name]);
sum_stats_T1 = pullSummaryFitInfo({x_reps_T1}, x_reps_T1.x_var_list, x_reps_T1.y_var_list);
%% LOAD: 5-field time sweep with replicates (T5): plotting firing rate, avalanche counts, and field values

% rep_list = {'A', 'B', 'C', 'D', 'E'};
% eta_list_uf5 = [4.0, 6.0, 8.0];
% eps_list_uf5 = [ -6.0, -7.0, -8.0, -9.0, -10.0, -12.0] ;%-14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_T5, num_ava_T5] = loadSimulationRunResults('run_f5_e1204_timesweep');
% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/timesweep_stim5_e12et04/';
rep_fn = dir([results_dir 'ava_decade_analysis*.mat']);
[~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');

x_reps_T5 = load([results_dir rep_fn(date_ord(1)).name]);
sum_stats_T5 = pullSummaryFitInfo({x_reps_T5}, x_reps_T5.x_var_list, x_reps_T5.y_var_list);

%% Also want the infinite-time constant simulations, from the UF run. 

% for T1 runs, x_info_T1 entries indicate epsilon = -8 and eta = 4. (These
% were run labeling parameters with epsilon = epsilon'*eta, so the param
% struct says epsilon = -2. 
eta_T1_infTC = 4;
eps_T1_infTC = -12;
sum_stats_T1_infTC = pullSummaryFitInfo(x_reps_uf1, eta_T1_infTC, eps_T1_infTC);

% for T5 runs, epsilon = -12 and eta = 4
eta_T5_infTC = 4;
eps_T5_infTC = -12;
sum_stats_T5_infTC = pullSummaryFitInfo(x_reps_uf5, eta_T5_infTC, eps_T5_infTC);

%% Figure 2: want a panel showing scaling, then a summary of goodness of fit
h_Fig = makeMyFigure(2.54*6.5*1.6, 2.54*6.5*1.6);
tiledlayout(h_Fig, 4, 4, 'TileSpacing', 'tight', 'Padding', 'tight');
h_sp1(1) = nexttile;%([2 1]);
h_sp1(2) = nexttile;%([2 1]);
h_sp1(3) = nexttile;%([2 1]);
h_sp1(4) = nexttile;%([2 1]);
% h_sp1(5) = nexttile();

h_sp5(1) = nexttile;%([2 1]);
h_sp5(2) = nexttile;%([2 1]);
h_sp5(3) = nexttile;%([2 1]);
h_sp5(4) = nexttile;%([2 1]);
% h_sp5(5) = nexttile;


%     h_sp1(1) = subplot(4, 7, [1 2]);
%     h_sp1(2) = subplot(4, 7, [3 4]);
%     h_sp1(3) = subplot(4, 7, [5 6]);
%     h_sp1(4) = subplot(4, 7, 7);
%    
%     
%     h_sp5(1) = subplot(4, 7, 7+[1 2]);
%     h_sp5(2) = subplot(4, 7, 7+[3 4]);
%     h_sp5(3) = subplot(4, 7, 7+[5 6]);
%     h_sp5(4) = subplot(4, 7, 14);
% set time scale for the plot
dyn_tau = 1;
% set colors: 
f1f5_colors = lines(2);


% the string 'A', 'B', etc. selects one of five repeated simulations  
plotSizeDurPDFScaling(x_reps_T1, eta_T1_infTC, eps_T1_infTC, dyn_tau, 'B', h_Fig, h_sp1, f1f5_colors(1, :));

plotSizeDurPDFScaling(x_reps_T5, eta_T5_infTC, eps_T5_infTC, dyn_tau, 'B', h_Fig, h_sp5, f1f5_colors(2, :));

% set layout
% these plots show tau, alpha, and gamma across multiple dynamical
% timescale values
h_spG1(1) = nexttile;
h_spG1(2) = nexttile;   % plot for alpha values
h_spG1(3) = nexttile;   % plot for gamma values
h_spG1(5) = nexttile;   % plot for gamma fit quality

% These plots show the alpha and tau fit quality statistics
h_spQ(1) = nexttile;%([2 1]);
h_spQ(2) = nexttile;
% create handle to the next tile for plot for gamma fit minus gamma prediction  
h_spG1(4) = nexttile;

% now make plots
set(gcf, 'CurrentAxes', h_spQ(1));
hold on
plot(log10(sum_stats_T1.rep_eta_vals), sum_stats_T1.tau_KS_min, 'o', 'color', f1f5_colors(1, :))
plot(log10(sum_stats_T5.rep_eta_vals), sum_stats_T5.tau_KS_min, 'o', 'color', f1f5_colors(2, :))

plot(1.4, squeeze(sum_stats_T1_infTC.tau_KS_min), 'o', 'color', f1f5_colors(1, :))
plot(1.4, squeeze(sum_stats_T5_infTC.tau_KS_min), 'o', 'color', f1f5_colors(2, :))


ylabel('KS statistic: \tau')
xlabel('\tau_F')
axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')


set(gcf, 'currentaxes', h_spQ(2));
hold on
plot(log10(sum_stats_T1.rep_eta_vals), sum_stats_T1.alpha_KS_min, 'o', 'color', f1f5_colors(1, :))
plot(log10(sum_stats_T5.rep_eta_vals), sum_stats_T5.alpha_KS_min, 'o', 'color', f1f5_colors(2, :))

plot(1.4, squeeze(sum_stats_T1_infTC.tau_KS_min), 'o', 'color', f1f5_colors(1, :))
plot(1.4, squeeze(sum_stats_T5_infTC.tau_KS_min), 'o', 'color', f1f5_colors(2, :))


ylabel('KS statistic: \alpha')
xlabel('\tau_F')
axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')


% put infinite-time constant point on the gamma-gamma_pred plot
set(gcf, 'CurrentAxes', h_spG1(5))
hold on
plot(1.4, squeeze(sum_stats_T1_infTC.gamma_fit_values-sum_stats_T1_infTC.gamma_pred_values), 'o', 'color', f1f5_colors(1, :))
plot(1.4, squeeze(sum_stats_T5_infTC.gamma_fit_values -sum_stats_T5_infTC.gamma_pred_values), 'o', 'color', f1f5_colors(2, :))

% this plots in the first 5 subplot listed in h_spG1
plotTauAlphaGammaSummary({x_reps_T1, x_reps_T5}, {'1 field', '5 fields'}, ...
    h_Fig, h_spG1);
%% labels

text(h_sp1(1), -.4, 1, 'A', 'FontSize',16, 'Units','normalized')
text(h_sp1(2), -.4, 1, 'B', 'FontSize',16, 'Units','normalized')
text(h_sp1(3), -.4, 1, 'C', 'FontSize',16, 'Units','normalized')
text(h_sp1(4), -.4, 1, 'D', 'FontSize',16, 'Units','normalized')


text(h_sp5(1), -.4, 1, 'E', 'FontSize',16, 'Units','normalized')
text(h_sp5(2), -.4, 1, 'F', 'FontSize',16, 'Units','normalized')
text(h_sp5(3), -.4, 1, 'G', 'FontSize',16, 'Units','normalized')
text(h_sp5(4), -.4, 1, 'H', 'FontSize',16, 'Units','normalized')


text(h_spG1(1), -.4, 1, 'I', 'FontSize',16, 'Units','normalized')
text(h_spG1(2), -.4, 1, 'J', 'FontSize',16, 'Units','normalized')
text(h_spG1(3), -.4, 1, 'K', 'FontSize',16, 'Units','normalized')
text(h_spG1(4), -.4, 1, 'O', 'FontSize',16, 'Units','normalized')
text(h_spQ(1), -.4, 1, 'M', 'FontSize',16, 'Units','normalized')
text(h_spQ(2), -.4, 1, 'N', 'FontSize',16, 'Units','normalized')


text(h_spG1(5), -.4, 1, 'L', 'FontSize',16, 'Units','normalized')


% set unused plots to be invisible
% set(h_sp1(5), 'Visible', 'off');
% set(h_sp5(5), 'Visible', 'off');
%%
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/figure2_timescales_' datestr(now, 'yyyymmdd')])
print(gcf, '-dsvg', ['avalanches_matlab_code/plots/paper_figures/figure2_timescales_' datestr(now, 'yyyymmdd')])
% %% Plots of KS stats and distance from gamma_pred = gamma_true
% 
% makeMyFigure(20, 14);
% 
% subplot(2, 3, 1)
% 
% subplot(2, 3, 2)
% 
% % subplot(2, 3, 3)
% % hold on
% % h_eb = plotXYErrorbars(sum_stats_T1.gamma_fit_values, sum_stats_T1.gamma_pred_values, sum_stats_T1.se_gamma_fit_values, sum_stats_T1.se_gamma_pred_values);
% % assignColorsToLines(h_eb, 0.5 - 0.5*gray(length(h_eb)))
% % eqline
% % xlabel('\gamma fit')
% % ylabel('\gamma prediction')
% % xlim([1 1.8])
% % ylim([1 1.8])
% % axis square
% % set(gca, 'color', 'none')
% 
% 
% 
% 
% % subplot(2, 3, 5)
% % hold on
% % plot(log10(sum_stats_T1.rep_eta_vals), log10(sum_stats_T1.alpha_x_min), 'o', 'color', [0 0 0])
% % plot(log10(sum_stats_T5.rep_eta_vals), log10(sum_stats_T5.alpha_x_min), 's', 'color', [0 0 1])
% % ylabel('lower cutoff: \alpha')
% % xlabel('\tau_F')
% % axis square
% % set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')
% 
% 
% subplot(2, 3, 3)
% hold on
% plot(log10(sum_stats_T1.rep_eta_vals), sum_stats_T1.gamma_fit_values - sum_stats_T1.gamma_pred_values, 'o', 'color', [0 0 0])
% plot(log10(sum_stats_T5.rep_eta_vals), sum_stats_T5.gamma_fit_values - sum_stats_T5.gamma_pred_values, 'o', 'color', [0 0 1])
% % h_eb = plotXYErrorbars(sum_stats_T5.gamma_fit_values, sum_stats_T5.gamma_pred_values, sum_stats_T5.se_gamma_fit_values, sum_stats_T5.se_gamma_pred_values);
% % assignColorsToLines(h_eb, 0.5 - 0.5*gray(length(h_eb)))
% 
% plot(1.4, squeeze(sum_stats_T1_infTC.gamma_fit_values-sum_stats_T1_infTC.gamma_pred_values), 'o-', 'color', [0 0 0])
% plot(1.4, squeeze(sum_stats_T5_infTC.gamma_fit_values -sum_stats_T5_infTC.gamma_pred_values), 'o-', 'color', [0 0 1])
% 
% 
% 
% 
% xlabel('\tau_F')
% ylabel('\gamma (fit - prediction)')
% % xlim([1 1.8])
% set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')
% 
% ylim([-0.25 0.25])
% 
% axis square
% set(gca, 'color', 'none')
% 
% print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/timesweep_1v5_' datestr(now, 'yyyymmdd')])
%% Plot lower cut offs 
makeMyFigure(20, 10);
x_min_bins = [3 5 9 16 28 50];
x_min_logbins = log10(x_min_bins); %[0.45 0.7 0.95 1.2 1.45 1.7];
n_xmins = length(x_min_logbins) - 1;
n_tauF = length(sum_stats_T1.rep_eta_vals);
xmin_cts_T1_tau = zeros(n_xmins, n_tauF);
xmin_cts_T5_tau = zeros(n_xmins, n_tauF);
xmin_cts_T1_alpha = zeros(n_xmins, n_tauF);
xmin_cts_T5_alpha = zeros(n_xmins, n_tauF);
for i_et = 1:size(sum_stats_T1.tau_x_min, 2)
    xmin_cts_T1_tau(:, i_et) = histcounts(log10(sum_stats_T1.tau_x_min(:, i_et)), x_min_logbins);
    xmin_cts_T5_tau(:, i_et) = histcounts(log10(sum_stats_T5.tau_x_min(:, i_et)), x_min_logbins);
    xmin_cts_T1_alpha(:, i_et) = histcounts(log10(sum_stats_T1.alpha_x_min(:, i_et)), x_min_logbins);
    xmin_cts_T5_alpha(:, i_et) = histcounts(log10(sum_stats_T5.alpha_x_min(:, i_et)), x_min_logbins);
    
end
sh = subplot(2, 3, 1);
imagesc(log10(sum_stats_T1.rep_eta_vals), x_min_logbins, xmin_cts_T1_tau)
set(gca, 'ytick', x_min_logbins, 'yticklabel', x_min_bins)
axis image
ylabel('S_{min}, N_F = 1')
xlabel('\tau_F')
% axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')
colormap(sh, 1-gray);


sh = subplot(2, 3, 2);
imagesc(log10(sum_stats_T1.rep_eta_vals), x_min_logbins, xmin_cts_T1_alpha)
set(gca, 'ytick', x_min_logbins, 'yticklabel', x_min_bins)
axis image
ylabel('D_{min}, N_F = 1')
xlabel('\tau_F')
% axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')
colormap(sh, 1-gray);

sh = subplot(2, 3, 4);
imagesc(log10(sum_stats_T5.rep_eta_vals), x_min_logbins, xmin_cts_T5_tau)
set(gca, 'ytick', x_min_logbins, 'yticklabel', x_min_bins)
axis image
ylabel('S_{min}, N_F = 5')
xlabel('\tau_F')
% axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')
colormap(sh, 1-gray);



sh = subplot(2, 3, 5);
imagesc(log10(sum_stats_T5.rep_eta_vals), x_min_logbins, xmin_cts_T5_alpha)
set(gca, 'ytick', x_min_logbins, 'yticklabel', x_min_bins)
axis image
ylabel('D_{min}, N_F = 5')
xlabel('\tau_F')
% axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(10.^(-1:1:1)', '%1.1f'), 'color', 'none')
colormap(sh, 1-gray);




%%
ava_mask_uf1 = double(squeeze(nanmean(num_ava_uf1, 1)) > 1e3);
hFig = plotSumStats(sum_stats_uf1, ava_mask_uf1);

% %%
% % Plot cdfs of firing rates
% hFR = plotFRSummaries('run_f5ultrafinesweep','A', x_info, num_ava_uf5);
%     
% hFig = plotSumStats(sum_stats_uf5, num_ava_uf5);

% %% number of avalanches
% 
% figure()
% eta_color = lines(length(eta_list));
% lh = 0*eta_list;
% hold on
% for ii = 1:length(eta_list)
% %     subplot(1, 3, ii)
%     ph = plot(eps_list, log10(squeeze(num_ava(:, :, ii))), ...
%         'color', eta_color(ii, :));
%     lh(ii) = ph(1);
% end
% legend(lh, num2str(eta_list', 'eta = %1.0f'))
% ylabel('log_{10} N_{ava}')
% 
% xlabel('\epsilon')
% y_lims = ylim;
% 
% 
% set(gca, 'ytick', ceil(2*y_lims(1))/2:0.5:floor(2*y_lims(2))/2, 'fontsize', 15)
% 
% % end
%%


%% pull useful data: ULTRAFINE RUN
num_eta = length(eta_list);
num_eps = length(eps_list);
num_reps = length(x_reps);
tau_values = nan(num_eps, num_eta, num_reps);
tau_KS_min = nan(num_eps, num_eta, num_reps);
tau_maxSurrKS = nan(num_eps, num_eta, num_reps);
se_tau_values = nan(num_eps, num_eta, num_reps);

alpha_KS_min = nan(num_eps, num_eta, num_reps);
alpha_maxSurrKS = nan(num_eps, num_eta, num_reps);
alpha_values = nan(num_eps, num_eta, num_reps);
se_alpha_values = nan(num_eps, num_eta, num_reps);

gamma_values = nan(num_eps, num_eta, num_reps);
se_gamma_values = nan(num_eps, num_eta, num_reps);

get_g_val = @(g_arr) cellfun(@(x) x.mhat, g_arr);



for ii = 1:num_reps
    
    alpha_values(:, :, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat');
    se_alpha_values(:, :, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'se_a_hat');
    alpha_KS_min(:, :, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'min_KS_stat');
    alpha_maxSurrKS(:, :, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'ks_surrogate', @(x) max(x));
   
    tau_values(:, :, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_tau_pfit', 'a_hat');
    se_tau_values(:, :, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_tau_pfit', 'a_hat');
    tau_KS_min(:, :, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_tau_pfit', 'min_KS_stat');
    tau_maxSurrKS(:, :, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_tau_pfit', 'ks_surrogate', @(x) max(x));

%     alpha_vals = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat');
%     alpha_vals = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat');
%     alpha_vals = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat');

    
    
    has_entry = cellfun(@(x) ~isempty(x), x_reps{ii}.all_gamma_pfit);

    gamma_arr = cell(num_eps, num_eta);
    gamma_arr(has_entry) = cellfun(@(x) get_g_val(x), x_reps{ii}.all_gamma_pfit(has_entry), ... 
        'UniformOutput', false);
    gamma_vals = nan(num_eps, num_eta);
    se_gamma_vals = nan(num_eps, num_eta);
    gamma_vals(has_entry) = cellfun(@(x) mean(x(2:5)), gamma_arr(has_entry) );
    se_gamma_vals(has_entry) = cellfun(@(x) std(x(2:5)), gamma_arr(has_entry) );


    gamma_values(:, :, ii) = gamma_vals;
    se_gamma_values(:, :, ii) = se_gamma_vals;
end

%% Plot summary data


%% 
figure(), 
for i_eta = 1:num_eta
    subplot(1, num_eta, i_eta)
    
    errorbar(eps_list, ave_gammas(:, i_eta), se_gammas(:, i_eta)) 
    hold on, 
    errorbar(eps_list, gamma_predicted(:, i_eta), gamma_predicted_se(:, i_eta), 'o-')
    xlabel('\epsilon')
    ylabel('\gamma')
    title(['\eta = ' num2str(eta_list(i_eta))])
    
end

legend({'fit', 'predicted'})
%% 5-field sweep (5FS): plotting firing rate, avalanche counts, and field values
[x5_info, x5_num_ava] = loadSimulationRunResults('run_f5finesweep');


%% Summaries from avalanche analysis
results_dir5 = 'avalanches_matlab_code/analysis_results/fields5finesweep/';
x5_reps = cell(1, 1);
for ii = 1:1
    rep_fn = dir([results_dir5 '*.mat']);
    [~, ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');
    
    x5_reps{ii} = load([results_dir5 rep_fn(ord(1)).name]);
end
%% Pull useful info, 5-field sweep (not ultra fine)
sum_stats5 = pullSummaryFitInfo(x5_reps, x5_reps{1}.x_var_list, x5_reps{1}.y_var_list);

%% plot 
hFig5 = plotSumStats(sum_stats5, x5_num_ava);

%%  Plot sum_stats and sum_stats5 on same axes for certain params
% positive: these are consistent, even though the replicated one has 5x the
% oberservations 

plot_fields = {'alpha_values', 'tau_values', 'gamma_values', ...
    'alpha_KS_min', 'alpha_medSurrKS', 'tau_KS_min', 'tau_medSurrKS'};

eta_val = 6;
nC = ceil(length(plot_fields)/2);

makeMyFigure(12, 30);
for ii = 1:length(plot_fields)
    subplot(nC, 2, ii)
    hold on
    plot(sum_stats.eps_list, ...
        squeeze(sum_stats.(plot_fields{ii})(:, sum_stats.eta_list == eta_val, :)),...
        'o-', 'linewidth', 1)
    
    plot(sum_stats5.eps_list, ...
        sum_stats5.(plot_fields{ii})(:, sum_stats5.eta_list == eta_val),...
        'ko-', 'linewidth', 1)
   ylabel(plot_fields{ii}, 'Interpreter', 'none')
   xlabel('\epsilon')
   
end

subplot(nC, 2, 2*nC)
deltag_pred = (sum_stats.alpha_values - 1)./(sum_stats.tau_values - 1) - sum_stats.gamma_values;
delta_g_pred5 = (sum_stats5.alpha_values - 1)./(sum_stats5.tau_values - 1) - sum_stats5.gamma_values;
hold on
plot(sum_stats.eps_list, squeeze(deltag_pred(:, sum_stats.eta_list == eta_val, :)),...
        'o-', 'linewidth', 1)
plot(sum_stats5.eps_list, delta_g_pred5(:, sum_stats.eta_list == eta_val),...
        'ko-', 'linewidth', 1)
ylabel('\gamma_{pred} - \gamma_{fit}')    

suptitle(['\eta = ' num2str(eta_val)])

%% 1-field Ultrafine sweep with replicates (1USR): plotting firing rate, avalanche counts, and field values

rep_list = {'A', 'B', 'C', 'D'};
eta_list = [4.0, 6.0, 8.0];
eps_list = [ -6.0,  -8.0, -10.0, -12.0];
% eps_list = [ -10.0, -12.0];
[xu1_info, num_ava_u1] = loadSimulationRunResults('run_f1ultrafinesweep', rep_list);



%% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/fields1ultrafinesweep/';
x_reps = cell(length(rep_list), 1);
for ii = 1:length(x_reps)
    rep_dir = [results_dir 'rep' rep_list{ii} '/' ];
    rep_fn = dir([rep_dir '*.mat']);
    x_reps{ii} = load([rep_dir rep_fn(1).name]);
end

sum_stats_uf1 = pullSummaryFitInfo(x_reps, eta_list, eps_list);
hFig = plotSumStats(sum_stats_uf1, num_ava);

%% Plot cdfs of firing rates
hFR = plotFRSummaries('run_f1ultrafinesweep','A', xu1_info, num_ava_u1);
    

%% number of avalanches

figure()
eta_color = lines(length(eta_list));
lh = 0*eta_list;
hold on
for ii = 1:length(eta_list)
%     subplot(1, 3, ii)
    ph = plot(eps_list, log10(squeeze(num_ava(:, :, ii))), ...
        'color', eta_color(ii, :));
    lh(ii) = ph(1);
end
legend(lh, num2str(eta_list', 'eta = %1.0f'))
ylabel('log_{10} N_{ava}')

xlabel('\epsilon')
y_lims = ylim;


set(gca, 'ytick', ceil(2*y_lims(1))/2:0.5:floor(2*y_lims(2))/2, 'fontsize', 15)

% end


%% Summaries from avalanche analysis
results_dir1 = 'avalanches_matlab_code/analysis_results/fields1finesweep/';
x1_reps = cell(1, 1);
for ii = 1:1
    rep_fn = dir([results_dir1 '*.mat']);
    [~, ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');
    
    x1_reps{ii} = load([results_dir1 rep_fn(ord(1)).name]);
end

%% 1-field sweep (1FS): plotting firing rate, avalanche counts, and field values
[x1_info, x1_num_ava] = loadSimulationRunResults('run_f1finesweep');


%% Pull useful info, 1-field sweep (not ultra fine)
sum_stats1 = pullSummaryFitInfo(x1_reps, x1_reps{1}.x_var_list, x1_reps{1}.y_var_list);

%% plot 
% mask by minimum FR or # avalanches
x1_mask = double(squeeze(mean(x1_num_ava, 1)) > 1e3);
hFig1 = plotSumStats(sum_stats1, x1_mask);
print(hFig1, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/f1finesweep_summary_' datestr(now, 'yyyymmdd')])

x5_mask = double(squeeze(mean(x5_num_ava, 1)) > 1e3);
hFig5 = plotSumStats(sum_stats5, x5_mask);
print(hFig5, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/f5finesweep_summary_' datestr(now, 'yyyymmdd')])

%% Plot gamma - gamma_pred for fixed eta again epsilon, for 1-field and 5-field
for eta_val = 4:2:8
    hFig = plotSumStat_deltaG({sum_stats1, sum_stats_uf1}, eta_val);
    
    set(gca, 'fontsize', 13, 'color','none')
    xlabel('\epsilon')
    legend({'1 field', 'long sample 1 fields'})
    x_lims = xlim;
%     xlim([-12 x_lims(2)])
    print(hFig, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/f1finesweep_deltaG_eta' num2str(eta_val, '%1.0f') '_' datestr(now, 'yyyymmdd')])

%     hFig = plotSumStat_deltaG({sum_stats5, sum_stats_uf5}, eta_val);
%     
%     set(gca, 'fontsize', 13, 'color','none')
%     xlabel('\epsilon')
%     legend({'5 field', 'long sample 5 fields'})
%     x_lims = xlim;
%     xlim([-12 x_lims(2)])
%     print(hFig, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/f5finesweep_deltaG_eta' num2str(eta_val, '%1.0f') '_' datestr(now, 'yyyymmdd')])

%     hFig = plotSumStat_deltaG({sum_stats_uf1}, eta_val_uf1);
%     
%     set(gca, 'fontsize', 13, 'color','none')
%     xlabel('\epsilon')
%     legend({'1 field', '5 fields'})
%     x_lims = xlim;
%     xlim([-12 x_lims(2)])
%     print(hFig, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/f1f5finesweep_deltaG_eta' num2str(eta_val, '%1.0f') '_' datestr(now, 'yyyymmdd')])

end

%%
eta_val = 8;
hFig = plotSumStat_deltaG({sum_stats_uf1, sum_stats_uf5}, eta_val);

set(gca, 'fontsize', 13, 'color','none')
xlabel('\epsilon')
legend({'1 field', '5 fields'})

%%
eta_val = 10;
hFig = plotSumStat_deltaG({sum_stats1, sum_stats5}, eta_val);

set(gca, 'fontsize', 13, 'color','none')
xlabel('\epsilon')
legend({'1 field', '5 fields'})

%% plot delta-gamma
hFig = plotSumStat_deltaG({sum_stats_uf1}, eta_list_uf1);
print(hFig, '-dpdf', ...
            ['avalanches_matlab_code/plots/paper_figures/ultrafine1_deltaG_alletas_' datestr(now, 'yyyymmdd')])

%% plot other coefficients
sel_fields = {'tau_values', 'alpha_values', 'gamma_values'};

% for eta_val = 4:2:8
    for i_fn = 1:length(sel_fields)

        sel_field_name = sel_fields{i_fn};
        hFig = plotSumStat_field({sum_stats_uf1}, eta_list_uf1, sel_field_name);
        print(hFig, '-dpdf', ...
            ['avalanches_matlab_code/plots/paper_figures/ultrafine1_' ...
            sel_field_name '_alletas_' datestr(now, 'yyyymmdd')])

    end
% end
%%

% for eta_val = 4:2:8
    for i_fn = 1:length(sel_fields)

        sel_field_name = sel_fields{i_fn};
        hFig = plotSumStat_field({sum_stats_uf1, sum_stats_uf5}, eta_val, sel_field_name);
        print(hFig, '-dpdf', ...
            ['avalanches_matlab_code/plots/paper_figures/uf15_ultrafinesweep_' ...
            sel_field_name '_eta' num2str(eta_val, '%1.0f') '_' datestr(now, 'yyyymmdd')])

    end
% end

%%
% NEXT : find folder for figure supblot, print out the various figures made
% in this script, add to the sketch file. 
% print(gcf, '-dpdf', 

%% helper functions


function hFR = plotFRSummaries(sim_string, rep_string, x_info, num_ava)
    
    ss_str = selectSimulation(sim_string, rep_string);
    
    setLoadingFunctionNames
    plot_field = 'cell_FR'; %'field_samples';%
    num_rows = length(y_var_list);
    num_cols = length(x_var_list);
    fr_bins = linspace(0, 1, 501); %linspace(-3, 3, 101); %
    n_samp = zeros(size(x_info));
    hFR = figure();
    for i_et = 1:length(x_var_list)
        for i_ep = 1:length(y_var_list)
            try

                fr_vals = cellfun(@(x) x.(plot_field), x_info(:, i_ep, i_et), 'UniformOutput', false);
                subplot(num_rows, num_cols, i_et + (i_ep -1)*num_cols)
                hold on
                for ii = 1:length(fr_vals)
                    n_samp(ii, i_ep, i_et) = size(x_info{ii, i_ep, i_et}.field_samples, 2);
                    histogram(fr_vals{ii}, fr_bins, 'normalization', 'cdf', 'DisplayStyle', 'stairs')
                end
                legend(num2str(num_ava(:, i_ep, i_et), 'n_{ava} = %1.0f'))
                xlabel(plot_field, 'Interpreter', 'none')
                ylabel('cdf over cells')
                title([y_lab_sum ' = ' num2str(y_var_list(i_ep)) ', ' x_lab_sum ' = ' num2str(x_var_list(i_et))])
            end
        end

    end

    % n_samp = unique(cellfun(@(x) size(x.field_samples, 2), x_info));
    suph = suptitle({[plot_field ', number of field samples is ' num2str(unique(n_samp)')]; ... 
        data_dir});
    suph.Interpreter = 'none';

end



function A_out = applyCellfunskipEmpty(x_struct, var_field_name, sub_field_name, fn_x)
% Applies the function fn_x to contents of x_struct.(var_field_name).(sub_field_name)
% when x_struct.(var_field_name) is non-empty. fn_x must return a scalar
% value. 
    has_entry = cellfun(@(x) ~isempty(x), x_struct.(var_field_name));
    
    A_out = nan(size(has_entry));
    A_out(has_entry) = cellfun(@(x) fn_x(x.(sub_field_name)), x_struct.(var_field_name)(has_entry));
    
end


%% plot summary data
function hFig = plotSumStats(sum_stats, ava_mask)

% Set zero-values of ava_mask to nan
ava_mask(ava_mask == 0) = nan;
% load needed variables
%%
tau_values = sum_stats.tau_values;
alpha_values = sum_stats.alpha_values;
gamma_values = sum_stats.gamma_fit_values;
alpha_KS_min = sum_stats.alpha_KS_min;
tau_KS_min = sum_stats.tau_KS_min;
alpha_medSurrKS = sum_stats.alpha_medSurrKS;
tau_medSurrKS = sum_stats.tau_medSurrKS;

eta_list = sum_stats.eta_list;
eps_list = sum_stats.eps_list;

gamma_pred_values = (alpha_values - 1)./(tau_values - 1);

if ndims(tau_values) == 3
    ave_taus = nanmean(tau_values, 3);
    ave_alphas = nanmean(alpha_values, 3);
    
    gamma_predicted = nanmean(gamma_pred_values, 3);
    ave_gammas = nanmean(gamma_values, 3);
    ave_alpha_KS = nanmedian(alpha_KS_min, 3);
    ave_tau_KS = nanmedian(tau_KS_min, 3);
    surr_alpha_KS = nanmedian(alpha_medSurrKS, 3);
    surr_tau_KS = nanmedian(tau_medSurrKS, 3);

else
    ave_taus = tau_values;
    ave_alphas = alpha_values;
    
    gamma_predicted = gamma_pred_values;
    ave_gammas = gamma_values;
    ave_alpha_KS = alpha_KS_min;
    ave_tau_KS = tau_KS_min;
    surr_alpha_KS = alpha_medSurrKS;
    surr_tau_KS = tau_medSurrKS;
end

%% now for the plots

hFig = makeMyFigure(30, 10);
colormap([ 1 1 1; parula])

subplot(2, 4, 1)

imagesc(eta_list, 1:length(eps_list), ave_taus.*ava_mask, [1.5 2.2])
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)
title('\tau (size scaling)')

subplot(2, 4, 2)
imagesc(eta_list, 1:length(eps_list), ave_alphas.*ava_mask, [1.5 2.2])
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)
title('\alpha (duration scaling)')


subplot(2, 4, 3)
imagesc(eta_list, 1:length(eps_list), gamma_predicted.*ava_mask, [1.1 1.3])
title('\gamma predicted')
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)


subplot(2, 4, 4)
imagesc(eta_list, 1:length(eps_list), ave_gammas.*ava_mask, [1.1 1.3])
title('\gamma from fit')
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)

% plot the KS values
ks_lim = [0 0.03];
subplot(2, 4, 5)
imagesc(eta_list, 1:length(eps_list), ave_alpha_KS.*ava_mask, ks_lim)
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)
title('\alpha (KS value)')


subplot(2, 4, 6)
imagesc(eta_list, 1:length(eps_list), ave_tau_KS.*ava_mask, ks_lim)
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)
title('\tau (KS value)')


subplot(2, 4, 7)
imagesc(eta_list, 1:length(eps_list), surr_alpha_KS.*ava_mask, ks_lim)
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)
title('\alpha (median surrogate KS)')


subplot(2, 4, 8)
imagesc(eta_list, 1:length(eps_list), surr_tau_KS.*ava_mask, ks_lim)
colorbar
set(gca, 'ytick', 1:length(eps_list), 'YTickLabel', eps_list)
title('\tau (median surrogate KS)')

end

function hFig = plotSumStat_deltaG(sum_stats_list, eta_vals)

hFig = makeMyFigure(8, 8);
hold on

for ii = 1:length(sum_stats_list)
    sum_stats = sum_stats_list{ii};
    deltag_pred = (sum_stats.alpha_values - 1)./(sum_stats.tau_values - 1) - sum_stats.gamma_values;
    % delta_g_pred5 = (sum_stats5.alpha_values - 1)./(sum_stats5.tau_values - 1) - sum_stats5.gamma_values;
    
    for jj = 1:length(eta_vals)
        eta_val = eta_vals(jj);
        try
            plot(sum_stats5.eps_list, deltag_pred(:, sum_stats.eta_list == eta_val),...
                'ko-', 'markersize', 10, 'linewidth', 1)
        catch
            y = squeeze(deltag_pred(:, sum_stats.eta_list == eta_val, :));
            n_rep = size(y, 2);
            errorbar(sum_stats.eps_list, nanmean(y, 2), nanstd(y, [], 2), 'o-', 'markersize', 10)
    %         plot(sum_stats.eps_list, squeeze(deltag_pred(:, sum_stats.eta_list == eta_val, :)),...
    %             'o-', 'linewidth', 1)
        end
    end
    ylabel('\gamma_{pred} - \gamma_{fit}')
    
end
title(['\eta = ' num2str(eta_vals)])

end

%% plot any field from sum_stats
function hFig = plotSumStat_field(sum_stats_list, eta_vals, field_name)

hFig = makeMyFigure(9, 9);
hold on

for ii = 1:length(sum_stats_list)
    sum_stats = sum_stats_list{ii};
    ss_field_vals = sum_stats.(field_name);
    
    %     deltag_pred = (sum_stats.alpha_values - 1)./(sum_stats.tau_values - 1) - sum_stats.gamma_values;
    % delta_g_pred5 = (sum_stats5.alpha_values - 1)./(sum_stats5.tau_values - 1) - sum_stats5.gamma_values;
    for jj = 1:length(eta_vals)
        eta_val = eta_vals(jj);
        try
            y = squeeze(ss_field_vals(:, sum_stats.eta_list == eta_val, :));
            n_rep = size(y, 2);
            errorbar(sum_stats.eps_list, nanmean(y, 2), nanstd(y, [], 2), 'o-', 'markersize', 10)
            %         plot(sum_stats.eps_list, squeeze(deltag_pred(:, sum_stats.eta_list == eta_val, :)),...
            %             'o-', 'linewidth', 1)
        catch
            
            
            plot(sum_stats.eps_list, ss_field_vals(:, sum_stats.eta_list == eta_val),...
                'ko-', 'markersize', 10, 'linewidth', 1)
        end
    end
    ylabel(field_name, 'interpreter', 'none')
    
end
title(['\eta = ' num2str(eta_vals)])
set(gca, 'color', 'none', 'fontsize', 12)

end

function h_eb = plotXYErrorbars(x, y, dx, dy)

    h_eb = errorbar(x, y, dy, dy, dx, dx, 'o');

end

%%
% %% plot summaries
% 
% plot_field = 'cell_FR'; %'field_samples';%
% num_rows = length(eps_list);
% num_cols = length(eta_list);
% fr_bins = linspace(0, 1, 501); %linspace(-3, 3, 101); %
% % n_samp = zeros(size(x_info));
% makeMyFigure(25, 25);
% for i_et = 1:length(eta_list)
%     for i_ep = 1:length(eps_list)
%         try
%             %%
%             fr_vals = cellfun(@(x) x.(plot_field), x_info(:, i_et, i_ep), 'UniformOutput', false);
%             subplot(num_rows, num_cols, i_et + (i_ep -1)*num_cols)
%             hold on
%             for ii = 1:length(fr_vals)
%                 
%                 histogram(fr_vals{ii}, fr_bins, 'normalization', 'cdf', 'DisplayStyle', 'stairs')
%             end
% %             legend(num2str(num_ava(:, i_et, i_ep), 'n_{ava} = %1.0f'))
%             axis([0 1 0 1])
%             set(gca, 'xtick', [], 'ytick', [])
%             if i_ep == length(eps_list)
%                 xlabel(plot_field, 'Interpreter', 'none')
%             end
%             if i_et == 1
%                 ylabel('cdf over cells')
%             end
%             title(['\epsilon = ' num2str(eps_list(i_ep)) ', \eta = ' num2str(eta_list(i_et))])
%         end
%     end
%     
% end
% not_empty = cellfun(@(x) ~isempty(x), x_info);
% n_samp = zeros(size(x_info));
% n_samp(not_empty) = cellfun(@(x) x.num_field_samples, x_info(not_empty));
% % n_samp = unique(cellfun(@(x) size(x.field_samples, 2), x_info));
% suph = suptitle({[plot_field ', number of field samples is ' num2str(unique(n_samp)')]; ... 
%     data_dir});
% suph.Interpreter = 'none';
% %% number of avalanches
% c_lim = [1 6];
% makeMyFigure(20, 10);
% 
% subplot(1, 2, 1)
% imagesc(eta_list, eps_list, log10(squeeze(num_ava)'), c_lim)
% ch = colorbar;
% ch.Ticks = c_lim(1):c_lim(2);
% ch.TickLabels = num2str(ch.Ticks', '10^%1.0f');
% ch.Label.String = 'avalanche count';
% set(gca, 'fontsize', 14)
% xlabel('\eta')
% ylabel('\epsilon')
% set(gca, 'ydir', 'normal')
% title(['Total simulation length: ' num2str(unique(nonzeros(n_samp(:)))', '%1.0g')])
% axis square
% num_ee = 20;
% % Xq = linspace(eta_list(1), eta_list(end), num_ee);
% % Yq = -linspace(eps_list(1), eps_list(end), 100);
% [Xq, Yq] = meshgrid(linspace(eta_list(1), eta_list(end), num_ee), ...
%     -linspace(eps_list(1), eps_list(end), num_ee));
% 
% X = eta_list;
% Y = -eps_list';
% V = squeeze(num_ava)';
% interp_count = (interp2(X, Y, (V), Xq, Yq));
% 
% subplot(1, 2, 2)
% imagesc(eta_list, eps_list, (log10(interp_count)), c_lim)
% hold on
% [x0, y0] = meshgrid(X, Y);
% plot(x0, -y0, 'w.')
% ch = colorbar;
% ch.Ticks = c_lim(1):c_lim(2);
% ch.TickLabels = num2str(ch.Ticks', '10^%1.0f');
% ch.Label.String = 'avalanche count';
% set(gca, 'fontsize', 14)
% xlabel('\eta')
% ylabel('\epsilon')
% set(gca, 'ydir', 'normal')
% title(['Interpolation'])
% axis square
% 
% print(gcf, '-dpdf', [plot_dir 'avalanche_counts'])