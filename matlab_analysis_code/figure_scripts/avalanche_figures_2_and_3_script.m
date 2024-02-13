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
fr_mats = cellfun(@(x) x.cell_FR/eff_bin_size, x_info_T5, 'UniformOutput',false);
figure()
tiledlayout(6, 8)
for ii = 1:length(fr_mats(:))
    nexttile
    histogram((fr_mats{ii}))
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
%% LOAD: 1-field fine sweep without replicates (F1): 
% rep_list = {'A', 'B', 'C', 'D', 'E'};
% eta_list_uf5 = [4.0, 6.0, 8.0];
% eps_list_uf5 = [ -6.0, -7.0, -8.0, -9.0, -10.0, -12.0] ;%-14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_f1, num_ava_f1] = loadSimulationRunResults('run_f1finesweep');
% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/fields1finesweep/';
rep_fn = dir([results_dir 'ava_decade_analysis*.mat']);
[~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');

x_reps_f1 = load([results_dir rep_fn(date_ord(1)).name]);
sum_stats_f1 = pullSummaryFitInfo({x_reps_f1}, x_reps_f1.x_var_list, x_reps_f1.y_var_list);
%%
%% LOAD: 5-field fine sweep without replicates (F): 
% rep_list = {'A', 'B', 'C', 'D', 'E'};
% eta_list_uf5 = [4.0, 6.0, 8.0];
% eps_list_uf5 = [ -6.0, -7.0, -8.0, -9.0, -10.0, -12.0] ;%-14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_f5, num_ava_f5] = loadSimulationRunResults('run_f5finesweep');
% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/fields5finesweep/';
rep_fn = dir([results_dir 'ava_decade_analysis*.mat']);
[~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');

x_reps_f5 = load([results_dir rep_fn(date_ord(1)).name]);
sum_stats_f5 = pullSummaryFitInfo({x_reps_f5}, x_reps_f5.x_var_list, x_reps_f5.y_var_list);

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
tiledlayout(h_Fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
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
dyn_tau = 0.5;
% set colors: 
f1f5_colors = lines(2);

%SOMETHING WAS BROKEN WHEN I "FIXED" THE GAMMA FITTING FUNCTION. MUST
%VERIFY WHY GAMMAFIT AND GAMMA PRED ARE SO DIFFERENT NOW. 
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
plot(log10(sum_stats_T1.rep_eta_vals), sum_stats_T1.tau_KS_min, 'o', 'color', f1f5_colors(1, :), 'LineWidth',1.5)
plot(log10(sum_stats_T5.rep_eta_vals), sum_stats_T5.tau_KS_min, 'o', 'color', f1f5_colors(2, :), 'LineWidth',1.5)

% plot(1.4, squeeze(sum_stats_T1_infTC.tau_KS_min), 'o', 'color', f1f5_colors(1, :), 'LineWidth',1.5)
% plot(1.4, squeeze(sum_stats_T5_infTC.tau_KS_min), 'o', 'color', f1f5_colors(2, :), 'LineWidth',1.5)


ylabel('KS statistic: \tau')
xlabel('\tau_F')
axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(1e4*10.^(-1:1:1)', '%1.0f'), 'color', 'none')


set(gcf, 'currentaxes', h_spQ(2));
hold on
plot(log10(sum_stats_T1.rep_eta_vals), sum_stats_T1.alpha_KS_min, 'o', 'color', f1f5_colors(1, :), 'LineWidth',1.5)
plot(log10(sum_stats_T5.rep_eta_vals), sum_stats_T5.alpha_KS_min, 'o', 'color', f1f5_colors(2, :), 'LineWidth',1.5)

% plot(1.4, squeeze(sum_stats_T1_infTC.tau_KS_min), 'o', 'color', f1f5_colors(1, :), 'LineWidth',1.5)
% plot(1.4, squeeze(sum_stats_T5_infTC.tau_KS_min), 'o', 'color', f1f5_colors(2, :), 'LineWidth',1.5)


ylabel('KS statistic: \alpha')
xlabel('\tau_F')
axis square
set(gca, 'xtick', -1:1:1, 'xticklabel', num2str(1e4*10.^(-1:1:1)', '%1.0f'), 'color', 'none')


% put infinite-time constant point on the gamma-gamma_pred plot
set(gcf, 'CurrentAxes', h_spG1(5))
hold on
plot(1.4, squeeze(sum_stats_T1_infTC.gamma_fit_values - sum_stats_T1_infTC.gamma_pred_values), 'o', 'color', f1f5_colors(1, :))
plot(1.4, squeeze(sum_stats_T5_infTC.gamma_fit_values - sum_stats_T5_infTC.gamma_pred_values), 'o', 'color', f1f5_colors(2, :))


try 
    % uses results of previous run, if they exist
    gamma_summary_T1T5 = plotTauAlphaGammaSummary({x_reps_T1, x_reps_T5}, {'1 field', '5 fields'}, ...
        h_Fig, h_spG1, gamma_summary_T1T5);
catch
% this plots in the first 5 subplot listed in h_spG1
    gamma_summary_T1T5 = plotTauAlphaGammaSummary({x_reps_T1, x_reps_T5}, {'1 field', '5 fields'}, ...
        h_Fig, h_spG1);
end
%

% labels
x_Lab = -.45;
text(h_sp1(1), x_Lab, 1, 'A', 'FontSize',16, 'Units','normalized')
text(h_sp1(2), x_Lab, 1, 'B', 'FontSize',16, 'Units','normalized')
text(h_sp1(3), x_Lab, 1, 'C', 'FontSize',16, 'Units','normalized')
text(h_sp1(4), x_Lab, 1, 'D', 'FontSize',16, 'Units','normalized')


text(h_sp5(1), x_Lab, 1, 'E', 'FontSize',16, 'Units','normalized')
text(h_sp5(2), x_Lab, 1, 'F', 'FontSize',16, 'Units','normalized')
text(h_sp5(3), x_Lab, 1, 'G', 'FontSize',16, 'Units','normalized')
text(h_sp5(4), x_Lab, 1, 'H', 'FontSize',16, 'Units','normalized')


text(h_spG1(1), x_Lab, 1, 'I', 'FontSize',16, 'Units','normalized')
text(h_spG1(2), x_Lab, 1, 'J', 'FontSize',16, 'Units','normalized')
text(h_spG1(3), x_Lab, 1, 'K', 'FontSize',16, 'Units','normalized')
text(h_spG1(4), x_Lab, 1, 'C', 'FontSize',16, 'Units','normalized')
text(h_spQ(1), x_Lab, 1, 'A', 'FontSize',16, 'Units','normalized')
text(h_spQ(2), x_Lab, 1, 'B', 'FontSize',16, 'Units','normalized')


text(h_spG1(5), x_Lab, 1, 'L', 'FontSize',16, 'Units','normalized')

% %% scratch
% 
% figure(), 
% errorbar(mean(gamma_summary_T1T5{1}.fit_gamma_range_start, 1), ...
%     std(gamma_summary_T1T5{1}.fit_gamma_range_start))
% hold on, 
% errorbar(mean(gamma_summary_T1T5{2}.fit_gamma_range_start, 1), ...
%     std(gamma_summary_T1T5{2}.fit_gamma_range_start))
% 
% % set unused plots to be invisible
% % set(h_sp1(5), 'Visible', 'off');
% % set(h_sp5(5), 'Visible', 'off');

print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/figure2_timescales_' datestr(now, 'yyyymmdd')])
print(gcf, '-dsvg', ['avalanches_matlab_code/plots/paper_figures/figure2_timescales_' datestr(now, 'yyyymmdd')])
%% Plot as a function of tau_F the gamma_fit vs decade

makeMyFigure(20, 10);
tiledlayout(1, 2);
d_val_g_T1 = cellfun(@(x) cellfun(@(y) y.x_cutoff(1), x'), x_reps_T1.all_gamma_pfit, 'UniformOutput', false);

fit_g_T1 = cellfun(@(x) cellfun(@(y) y.mhat, x'), x_reps_T1.all_gamma_pfit, 'UniformOutput', false);
tau_F_cmap = linspace(0, 1, size(d_val_g_T5, 2))'*[1 0.5 0];

nexttile
hold on
for ii = 1:length(fit_g_T1)
    try 
        y = cell2mat(fit_g_T1(:, ii));
        x = cell2mat(d_val_g_T1(:, ii));
    catch
    end

    errorbar(x(1, :), mean(y, 1), std(y, [], 1), 'linewidth', 1, 'color', tau_F_cmap(ii, :))
end
y_lim1 = ylim;
title({'\gamma_{fit} stability across decades'; 'one field'})
xlabel('lower bound of fit range (log_{10} D)')
ylabel('\gamma_{fit} (mean +/- SE across networks)')
set(gca, 'color', 'none')
axis square
tau_F_vals = 1e4*x_reps_T1.x_var_list;
legend(num2str(tau_F_vals', 'tau_F = %1.0g'), 'Location','eastoutside')

d_val_g_T5 = cellfun(@(x) cellfun(@(y) y.x_cutoff(1), x'), x_reps_T5.all_gamma_pfit, 'UniformOutput', false);

fit_g_T5 = cellfun(@(x) cellfun(@(y) y.mhat, x'), x_reps_T5.all_gamma_pfit, 'UniformOutput', false);
nexttile
hold on
for ii = 1:length(fit_g_T5)
    try 
        y = cell2mat(fit_g_T5(:, ii));
        x = cell2mat(d_val_g_T5(:, ii));
    catch
    end

    errorbar(x(1, :), mean(y, 1), std(y, [], 1), 'linewidth', 1, 'color', tau_F_cmap(ii, :))
end
ylim(y_lim1)
title({'\gamma_{fit} stability across decades'; 'five fields'})
xlabel('lower bound of fit range (log_{10} D)')
ylabel('\gamma_{fit} (mean +/- SE across networks)')
set(gca, 'color', 'none')
axis square
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/gamma_fit_fig2supp_' datestr(now, 'yyyymmdd')])

%% AGGREGATE METRICS FOR FIGURE 3

% compute a fit quality based on KS statistic
surr_test_fn = @(x, field_name) x.([field_name '_KS_min']) <= x.([field_name '_KS_p95']);
surr_z_fn = @(x, field_name) (x.([field_name '_KS_min']) - x.([field_name '_medSurrKS']))./x.([field_name '_stdSurrKS']);
%%
pass_surr_test_tau_UF1 = surr_test_fn(sum_stats_uf1, 'tau');
pass_surr_test_tau_UF5 = surr_test_fn(sum_stats_uf5, 'tau');
pass_surr_test_alpha_UF1 = surr_test_fn(sum_stats_uf1, 'alpha');
pass_surr_test_alpha_UF5 = surr_test_fn(sum_stats_uf5, 'alpha');

surr_z_tau_UF1 = surr_z_fn(sum_stats_uf1, 'tau'); 
surr_z_alpha_UF1 = surr_z_fn(sum_stats_uf1, 'alpha');
surr_z_tau_UF5 = surr_z_fn(sum_stats_uf5, 'tau'); 
surr_z_alpha_UF5 = surr_z_fn(sum_stats_uf5, 'alpha'); 

pass_surr_test_tau_F1 = surr_test_fn(sum_stats_f1, 'tau');
pass_surr_test_alpha_F1 = surr_test_fn(sum_stats_f1, 'alpha');
pass_surr_test_tau_F5 = surr_test_fn(sum_stats_f5, 'tau');
pass_surr_test_alpha_F5 = surr_test_fn(sum_stats_f5, 'alpha');


%%
surr_z_tau_F1 = surr_z_fn(sum_stats_f1, 'tau'); 
surr_z_alpha_F1 = surr_z_fn(sum_stats_f1, 'alpha');
surr_z_tau_F5 = surr_z_fn(sum_stats_f5, 'tau'); 
surr_z_alpha_F5 = surr_z_fn(sum_stats_f5, 'alpha'); 


%%
% lower cutoffs 
s_min_UF1 = sum_stats_uf1.tau_x_min;
s_min_UF5 = sum_stats_uf5.tau_x_min;
d_min_UF1 = sum_stats_uf1.alpha_x_min;
d_min_UF5 = sum_stats_uf5.alpha_x_min;


eta_UF1_vals = sum_stats_uf1.eta_list;
eta_UF1_colors = cool(length(eta_UF1_vals));
eta_UF5_vals = sum_stats_uf5.eta_list;
eta_UF5_colors = cool(length(eta_UF5_vals));

eps_UF1_vals = sum_stats_uf1.eps_list;
eps_UF5_vals = sum_stats_uf5.eps_list;

eta_F1_vals = sum_stats_f1.eta_list;
eta_F5_vals = sum_stats_f5.eta_list;

eps_F1_vals = sum_stats_f1.eps_list;
eps_F5_vals = sum_stats_f5.eps_list;


% this controls the maximum cutoffs allowed for avalanche sizes, durations
% and the maximum 'z-value' for the surrogate distribution analysis
max_ava_size = 30;
max_z_val = 4;

[ave_tau_vals_uf1, frac_tauSc_UF1] = extractExponentsGoodScaling(sum_stats_uf1, 'size', max_ava_size, max_z_val);
[ave_alpha_vals_uf1, frac_alphaSc_UF1] = extractExponentsGoodScaling(sum_stats_uf1, 'dur', max_ava_size, max_z_val);
[ave_tau_vals_uf5, frac_tauSc_UF5] = extractExponentsGoodScaling(sum_stats_uf5, 'size', max_ava_size, max_z_val);
[ave_alpha_vals_uf5, frac_alphaSc_UF5] = extractExponentsGoodScaling(sum_stats_uf5, 'dur', max_ava_size, max_z_val);

[ave_tau_vals_f1, frac_tauSc_F1] = extractExponentsGoodScaling(sum_stats_f1, 'size', max_ava_size, 2*max_z_val);
[ave_alpha_vals_f1, frac_alphaSc_F1] = extractExponentsGoodScaling(sum_stats_f1, 'dur', max_ava_size, 2*max_z_val);
[ave_tau_vals_f5, frac_tauSc_F5] = extractExponentsGoodScaling(sum_stats_f5, 'size', max_ava_size, 2*max_z_val);
[ave_alpha_vals_f5, frac_alphaSc_F5] = extractExponentsGoodScaling(sum_stats_f5, 'dur', max_ava_size, 2*max_z_val);


%% calculate fitted gammas using the new methods. 

gamma_summary_UF1 = plotTauAlphaGammaSummary(x_reps_uf1);
%%
gamma_summary_UF5 = plotTauAlphaGammaSummary(x_reps_uf5);
%%

gamma_summary_F1 = plotTauAlphaGammaSummary({x_reps_f1});

gamma_summary_F5 = plotTauAlphaGammaSummary({x_reps_f5});

%% 
gamma_UF1_fits = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma, -1), gamma_summary_UF1', 'UniformOutput',false));
gamma_UF1_range = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma_range, -1), gamma_summary_UF1', 'UniformOutput',false));

gamma_UF1_fits_r2 = gamma_UF1_fits;
gamma_UF1_fits_r2(gamma_UF1_range < 2) = nan;
gamma_UF1_fits_r2_ave = squeeze(mean(gamma_UF1_fits_r2, 1, 'omitnan'));
gamma_UF1_fits_r2_std = squeeze(std(gamma_UF1_fits_r2, [], 1, 'omitnan'));

gamma_UF5_fits = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma, -1), gamma_summary_UF5', 'UniformOutput',false));
gamma_UF5_range = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma_range, -1), gamma_summary_UF5', 'UniformOutput',false));

gamma_UF5_fits_r2 = gamma_UF5_fits;
gamma_UF5_fits_r2(gamma_UF5_range < 2) = nan;
gamma_UF5_fits_r2_ave = squeeze(mean(gamma_UF5_fits_r2, 1, 'omitnan'));
gamma_UF5_fits_r2_std = squeeze(std(gamma_UF5_fits_r2, [], 1, 'omitnan'));

% alternative - this fit_gamma is from the old decade processing
fit_gamma_uf1 = extractExponentsGoodScaling(sum_stats_uf1, 'gamma_fit', max_ava_size, max_z_val);
pred_gamma_uf1 = extractExponentsGoodScaling(sum_stats_uf1, 'gamma_pred', max_ava_size, max_z_val);
fit_gamma_uf5 = extractExponentsGoodScaling(sum_stats_uf5, 'gamma_fit', max_ava_size, max_z_val);
pred_gamma_uf5 = extractExponentsGoodScaling(sum_stats_uf5, 'gamma_pred', max_ava_size, max_z_val);


%% Figure 3: layout
h_Fig3 = makeMyFigure(2.54*6.5*1.6, 2.54*6.5*1.6);
tiledlayout(h_Fig3, 8, 5, 'TileSpacing', 'compact', 'Padding', 'tight');
h_ee = zeros(10, 1);
h_ee(1) = nexttile;%([2 1]);
h_ee(2) = nexttile([2 2]);%([2 1]);
h_ee(3) = nexttile([2 2]);%([2 1]);
h_ee(4) = nexttile;%([2 1]);
h_ee(5) = nexttile;%([2 1]);
h_ee(6) = nexttile;%([2 1]);
h_ee(7) = nexttile;%([2 1]);

%%
dg_lim = [-0.2 0.2];
n_dg_col = 4;
set(h_Fig3, 'currentaxes', h_ee(2))
imagesc(eta_UF1_vals, eps_UF1_vals, gamma_UF1_fits_r2_ave - pred_gamma_uf1, dg_lim);
colorbar;
colormap(h_ee(2), usa(n_dg_col));
axis square
set(gca, 'ydir', 'normal')
title('\gamma_{fit} - \gamma_{pred} (one field)')
xlabel('\eta')
ylabel('\epsilon') 

set(h_Fig3, 'currentaxes', h_ee(3))
imagesc(eta_UF5_vals, eps_UF5_vals, gamma_UF5_fits_r2_ave - pred_gamma_uf5, dg_lim);
colorbar;
colormap(h_ee(3), usa(n_dg_col));
axis square
set(gca, 'ydir', 'normal')
title('\gamma_{fit} - \gamma_{pred} (five fields)')
xlabel('\eta')
ylabel('\epsilon') 


%%

% these should all be in the same order: UF1, UF1, UF5, UF5. First tau,
% then alpha. In gamma, first the fit, then the prediction. The x- and
% y-values arrays are organized as such. 
ave_exp_arr = {ave_tau_vals_uf1, ave_alpha_vals_uf1, ave_tau_vals_uf5, ave_alpha_vals_uf5};
frac_sc_arr = {frac_tauSc_UF1, frac_alphaSc_UF1, frac_tauSc_UF5, frac_alphaSc_UF5};

ave_gamma_arr = {fit_gamma_uf1, pred_gamma_uf1, fit_gamma_uf5, pred_gamma_uf5};

x_vals_arr = {eta_UF1_vals, eta_UF1_vals, eta_UF5_vals, eta_UF5_vals};
y_vals_arr = {eps_UF1_vals, eps_UF1_vals, eps_UF5_vals, eps_UF5_vals};
gamma_title_str = {'\gamma_{fit}', '\gamma_{pred}', '\gamma_{fit}', '\gamma_{pred}'};
at_title_str = {'\tau', '\alpha', '\tau', '\alpha'};
title_str = {'1 field', '1 field', '5 fields', '5 fields'};
%%
% tau_scaling_UF5 = s_min_UF5 < max_ava_size & surr_z_tau_UF5 < max_z_val;
% tau_scaling_UF5(isnan(tau_scaling_UF5)) = false;
% frac_tauSc_UF5 = mean(tau_scaling_UF5, 3);
% 
% alpha_scaling_UF1 = d_min_UF1 < max_ava_size & surr_z_alpha_UF1 < max_z_val;
% alpha_scaling_UF1(isnan(alpha_scaling_UF1)) = false;
% frac_alphaSc_UF1 = mean(alpha_scaling_UF1, 3);
% 
% alpha_scaling_UF5 = d_min_UF5 < max_ava_size & surr_z_alpha_UF5 < max_z_val;
% alpha_scaling_UF5(isnan(alpha_scaling_UF5)) = false;
% frac_alphaSc_UF5 = mean(alpha_scaling_UF5, 3);

%% FIGURE 3 - SCALING DIAGRAM

hFig3 = makeMyFigure(2.54*6.5, 2.54*6);
nR = 3;
nC = 4;
tiledlayout(nR, nC, 'TileSpacing','compact', 'Padding','tight');

h_sp3 = zeros(nR*nC, 1);
for ii = 1:length(h_sp3)
    h_sp3(ii) = nexttile;
end
% h_sp3 = h_sp3([1 4 7 10 2 5 8 11 3 6 9 12]);
% 
for ii = 1:length(ave_exp_arr)
    set(hFig3, 'currentaxes', h_sp3(3*ii-2))
    imagesc(x_vals_arr{ii}, y_vals_arr{ii}, frac_sc_arr{ii}, [0 1])
    title([title_str{ii} ', ' at_title_str{ii} ', frac +'])
    colorbar

    set(hFig3, 'currentaxes', h_sp3(3*ii-1))
    imagesc(x_vals_arr{ii}, y_vals_arr{ii}, ave_exp_arr{ii}, [1.8 2.2])
    title([title_str{ii} ', ' at_title_str{ii}])
    colorbar
    
    set(hFig3, 'currentaxes', h_sp3(3*ii))
    imagesc(x_vals_arr{ii}, y_vals_arr{ii}, ave_gamma_arr{ii}, [1 1.5])
    title([title_str{ii} ', ' gamma_title_str{ii}])
    colorbar
end

%% FIGURE 3 - supplement (showing minimum cutoffs, etc. )
hFig3S = makeMyFigure(2.54*6.5, 2.54*6);

tiledlayout(3, 4, 'TileSpacing','tight', 'Padding','tight');

h_sp3S = zeros(12, 1);
for ii = 1:length(h_sp3S)
    h_sp3S(ii) = nexttile;
end

%% in progress: move to function. NOT SURE HOW TO PRESENT THESE ANALYSES
% WHAT IS THE MAIN POINT? 1) scaling over a range of parameters. 


% row 1: UF1 simulation, different eta curves
sp_ord = [1 2 5 6];

for i_at = 1:4
    set(gcf, 'CurrentAxes', h_sp3S(sp_ord(i_at)))

    if i_at == 1
        plot_val = s_min_UF1;
        pass_surr_test = pass_surr_test_tau_UF1;
        z_surr_test = surr_z_tau_UF1;
        eps_vals = eps_UF1_vals;
        eta_vals = eta_UF1_vals;
        eta_colors = eta_UF1_colors;
    elseif i_at == 2
        plot_val = d_min_UF1;
        pass_surr_test = pass_surr_test_alpha_UF1;
        z_surr_test = surr_z_alpha_UF1;
        eps_vals = eps_UF1_vals;
        eta_vals = eta_UF1_vals;
        eta_color = eta_UF1_colors;
    elseif i_at == 3
        plot_val = s_min_UF5;
        pass_surr_test = pass_surr_test_tau_UF5;
        z_surr_test = surr_z_tau_UF5;
        eps_vals = eps_UF5_vals;
        eta_vals = eta_UF5_vals;
        eta_colors = eta_UF5_colors;
    elseif i_at == 4
        plot_val = d_min_UF5;
        pass_surr_test = pass_surr_test_alpha_UF5;
        z_surr_test = surr_z_alpha_UF5;
        eps_vals = eps_UF5_vals;
        eta_vals = eta_UF5_vals;
        eta_colors = eta_UF5_colors;
    end

max_Z = 3;  % maximum z-score for "OK" fit
hold on
for i_eta = 1:length(eta_UF1_vals)
    
    y_vals = squeeze(plot_val(:, i_eta, :));
    

    plot(eps_vals, y_vals, '.', 'color', eta_colors(i_eta, :))

    % use an open circle for fits that are kind of good: in which the fit
    % KS score is within 3 SD of the surrogate average. 
    z_s_t = squeeze(z_surr_test(:, i_eta, :));
    y_z = y_vals;
    y_z(z_s_t > max_Z) = nan;
    plot(eps_vals, y_z, 'o', 'color', eta_colors(i_eta, :))

    % use star for fits that are good: passing the surrogate test at 5%
    p_s_t = squeeze(pass_surr_test(:, i_eta, :));
    y_p = y_vals;
    y_p(p_s_t) = nan;
    plot(eps_vals, y_p, '*', 'color', eta_colors(i_eta, :));
    ylabel('minimum cutoff')
end
end


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