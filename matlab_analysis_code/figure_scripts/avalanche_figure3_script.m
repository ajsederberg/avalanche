% Script for eta-epsilon plots (Figure 3 or 4)
dt_val = 0.008; % set effective time step to be 8 ms
                % 10000 steps = 80 s. Longest avalanches are ~8 s. 
                % 200 loops = 260 minutes -> about 4 hours recording
                % 1000 loops is about 20 hours recording 

                
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


%% Look at particular parameter combinations
% set parameters 

eta_f1 = 6;
eps_f1 = -18;

% generate multi-panel figure with all size/duration/size-duration scaling
% plots. 
h_Fig = makeMyFigure(2.54*6.5*1.6, 2.54*9*1.6);
tiledlayout(h_Fig, 6, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
h_sp1(1) = nexttile;%([2 1]);
h_sp1(2) = nexttile;%([2 1]);
h_sp1(3) = nexttile;%([2 1]);
h_sp1(4) = nexttile;%([2 1]);

for ii = 1:4
    h_sp2(ii) = nexttile;
end
for ii = 1:4
    h_sp3(ii) = nexttile;
end
for ii = 1:4
    h_sp4(ii) = nexttile;
end
for ii = 1:4
    h_sp5(ii) = nexttile;
end
for ii = 1:4
    h_sp6(ii) = nexttile;
end

f1f5_colors = lines(2);

plotSizeDurPDFScaling(x_reps_f1, eta_f1, eps_f1, [], [], h_Fig, h_sp1, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_uf1{1}, eta_f1, eps_f1, [], [], h_Fig, h_sp2, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_uf1{2}, eta_f1, eps_f1, [], [], h_Fig, h_sp3, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_uf1{3}, eta_f1, eps_f1, [], [], h_Fig, h_sp4, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_uf1{4}, eta_f1, eps_f1, [], [], h_Fig, h_sp5, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_uf1{5}, eta_f1, eps_f1, [], [], h_Fig, h_sp6, f1f5_colors(1, :));

plot_name = ['sim_length_compare_stim1_e' num2str(eps_f1) 'et' num2str(eta_f1)];

print(gcf, '-dpdf', ['avalanches_matlab_code/plots/sim_length_scaling_plots/' plot_name '.pdf'])
%% Insight: if lower cutoff is small (<10) in both size and duration, 
% scaling seems clear

% first check lower cutoff for size and for duration

d_min = sum_stats_f1.alpha_x_min;
s_min = sum_stats_f1.tau_x_min;
eps_vals = -sum_stats_f1.eps_list;
eta_vals = sum_stats_f1.eta_list;

figure()
tiledlayout(2,3)

nexttile
imagesc(eta_vals, eps_vals, s_min, [0 50])
title('minimum cutoff size')
xlabel('\eta')
ylabel('\epsilon')
nexttile
imagesc(eta_vals, eps_vals, d_min, [0 50])
title('minimum cutoff duration')

nexttile
imagesc(eta_vals, eps_vals, s_min <= 10 & d_min <= 10)
title('both are less than 10')
% calculate median, mean, and max firing rates across population
% go look at Hengen & Turrigiano 2016 for realistic firing
                % rate distributions -
has_entry = cellfun(@(x) ~isempty(x), x_info_f1);
fr_median = nan(size(has_entry));
fr_median(has_entry) = cellfun(@(x) median(x.cell_FR)/dt_val, x_info_f1(has_entry));
fr_median = squeeze(fr_median);

fr_max = nan(size(has_entry));
fr_max(has_entry) = cellfun(@(x) max(x.cell_FR)/dt_val, x_info_f1(has_entry));
fr_max = squeeze(fr_max);

nexttile
imagesc(eta_vals, eps_vals, log10(fr_median), [-2 1])
colorbar
title('median firing rate (log-scale)')

nexttile
imagesc(eta_vals, eps_vals, log10(fr_max), [-1 2])
colorbar
title('maximum firing rate (log-scale)')

%% Define eta - epsilon trajectory
eps_trace = [8 6 8];
eta_trace = [8 6 4];
eps_ind_trace = 0*eps_trace;
eta_ind_trace = 0*eta_trace;
for i_ee = 1:length(eps_trace)
    eps_ind_trace(i_ee) = find(eps_vals == eps_trace(i_ee));
    eta_ind_trace(i_ee) = find(eta_vals == eta_trace(i_ee));
end
eps_vals_uf = -x_reps_uf1{1}.y_var_list;
eta_vals_uf = x_reps_uf1{1}.x_var_list;
eps_ind_uf = 0*eps_trace;
eta_ind_uf = 0*eta_trace;
for i_ee = 1:length(eps_trace)
    eps_ind_uf(i_ee) = find(eps_vals_uf == eps_trace(i_ee));
    eta_ind_uf(i_ee) = find(eta_vals_uf == eta_trace(i_ee));
end
%%
% eps_ind_trace(2) = find(eps_vals == 6);
% eta_Ind_trace(2) = find(eta_vals == 4);
% 
% eps_ind_trace(3)= find(eps_vals == 8);
% eta_Ind_trace(3) = find(eta_vals == 4);

trace_fr_vals = zeros(1024, length(eps_ind_trace));
cmap = hsv(1024);
marker_list = {'o', 's', '+'};
figure()
hold on
for ii = 1:length(eps_ind_trace)
    
    trace_fr_vals(:, ii) = log10(x_info_f1{1, eps_ind_trace(ii), eta_ind_trace(ii)}.cell_FR/dt_val);
    if ii ==1 
        [~, ordC] = sort(trace_fr_vals(:, ii));

        [~, ord_ordC] = sort(ordC);
    end
    [~, ord] = sort(trace_fr_vals(:, ii));
    for jj = 1:1024
        plot(trace_fr_vals(ord(jj), ii), jj/1024,'marker', marker_list{ii}, 'Color',cmap(ord_ordC(ord(jj)), :))
    end
%     histogram(fr_vals, 'Normalization','cdf', 'numbins', 30, 'DisplayStyle','stairs', ...
%         'LineWidth', 1)
end

%%
h_Fig = makeMyFigure(2.54*6.5*1.6, 2.54*9*1.6);
tiledlayout(h_Fig, 6, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
h_sp1(1) = nexttile;%([2 1]);
h_sp1(2) = nexttile;%([2 1]);
h_sp1(3) = nexttile;%([2 1]);
h_sp1(4) = nexttile;%([2 1]);

for ii = 1:4
    h_sp2(ii) = nexttile;
end
for ii = 1:4
    h_sp3(ii) = nexttile;
end
for ii = 1:4
    h_sp4(ii) = nexttile;
end

for ii = 1:4
    h_sp5(ii) = nexttile;
end
for ii = 1:4
    h_sp6(ii) = nexttile;
end


plotSizeDurPDFScaling(x_reps_f1, eta_vals(eta_ind_trace(1)), ...
    -eps_vals(eps_ind_trace(1)), [], [], h_Fig, h_sp1, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_f1, eta_vals(eta_ind_trace(2)), ...
    -eps_vals(eps_ind_trace(2)), [], [], h_Fig, h_sp2, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_f1, eta_vals(eta_ind_trace(3)), ...
    -eps_vals(eps_ind_trace(3)), [], [], h_Fig, h_sp3, f1f5_colors(1, :));


plotSizeDurPDFScaling(x_reps_uf1{2}, eta_vals(eta_ind_trace(1)), ...
    -eps_vals(eps_ind_trace(1)), [], [], h_Fig, h_sp4, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_uf1{2}, eta_vals(eta_ind_trace(2)), ...
    -eps_vals(eps_ind_trace(2)), [], [], h_Fig, h_sp5, f1f5_colors(1, :));
plotSizeDurPDFScaling(x_reps_uf1{2}, eta_vals(eta_ind_trace(3)), ...
    -eps_vals(eps_ind_trace(3)), [], [], h_Fig, h_sp6, f1f5_colors(1, :));


%% gamma fits : very reliable across replicates
gamma_fit_arr = consolidateGammaFits(x_reps_uf1, eps_ind_trace, eta_ind_trace);

%% 
makeMyFigure(30, 10)

for jj = 1:length(eps_trace)
    nexttile
    hold on

    y_vals = cell2mat(cellfun(@(x) x.full_gamma_fit_mean, gamma_fit_arr(:, jj)', 'UniformOutput',false));
    dy_vals = cell2mat(cellfun(@(x) x.full_gamma_fit_se, gamma_fit_arr(:, jj)', 'UniformOutput',false));

    errorbar(y_vals(:, 1:end-1), dy_vals(:, 1:end-1), 'k')
    errorbar(y_vals(:, end), dy_vals(:, end), 'LineWidth',1)
    title(['\epsilon = ' num2str(eps_trace(jj)) ', \eta = ' num2str(eta_trace(jj))])
    xlabel('start point fit')
    ylabel('\gamma (SE fit)')
end
%% function to consolidate fits - UF
res_str = consolidateAlphaTauFits(x_reps_uf1, eps_trace(1), eta_trace(1), 'alpha');

all_alpha_str = cell(length(eps_vals_uf), length(eta_vals_uf));
all_tau_str = cell(length(eps_vals_uf), length(eta_vals_uf));
for ii = 1:length(eps_vals_uf)
    for jj = 1:length(eta_vals_uf)
        try
        all_alpha_str{ii, jj} = consolidateAlphaTauFits(x_reps_uf1, ...
            eps_vals_uf(ii), eta_vals_uf(jj), 'alpha');
        end
        try
        all_tau_str{ii, jj} = consolidateAlphaTauFits(x_reps_uf1, ...
            eps_vals_uf(ii), eta_vals_uf(jj), 'tau');
        end
    end
end

%% plots 
has_entry = cellfun(@(x) ~isempty(x), all_alpha_str);
% alpha information
alpha_vals = zeros(size(all_alpha_str));
alpha_vals(has_entry) = cellfun(@(x) x.exp_val_est, all_alpha_str(has_entry));

se_alpha_vals = zeros(size(all_alpha_str));
se_alpha_vals(has_entry) = cellfun(@(x) x.se_exp_val_est, all_alpha_str(has_entry));

ks_alpha_vals = zeros(size(all_alpha_str));
ks_alpha_vals(has_entry) = cellfun(@(x) x.ks_min, all_alpha_str(has_entry));

min_size_vals = zeros(size(all_alpha_str));
min_size_vals(has_entry) = cellfun(@(x) x.x_min_est, all_alpha_str(has_entry));

% tau information
tau_vals = zeros(size(all_tau_str));
tau_vals(has_entry) = cellfun(@(x) x.exp_val_est, all_tau_str(has_entry));

se_tau_vals = zeros(size(all_tau_str));
se_tau_vals(has_entry) = cellfun(@(x) x.se_exp_val_est, all_tau_str(has_entry));

ks_tau_vals = zeros(size(all_tau_str));
ks_tau_vals(has_entry) = cellfun(@(x) x.ks_min, all_tau_str(has_entry));

min_dur_vals = zeros(size(all_tau_str));
min_dur_vals(has_entry) = cellfun(@(x) x.x_min_est, all_tau_str(has_entry));


gamma_pred = (alpha_vals - 1)./(tau_vals - 1);
%% estimate gamma_pred_se
num_draws = 1000;
alpha_val_sur = alpha_vals(:) + repmat(se_alpha_vals(:), 1, num_draws).*randn(numel(alpha_vals), num_draws);
tau_val_sur = tau_vals(:) + repmat(se_tau_vals(:), 1, num_draws).*randn(numel(alpha_vals), num_draws);
gamma_pred_sur = (alpha_val_sur - 1)./(tau_val_sur - 1);

se_gamma_pred = reshape(std(gamma_pred_sur, [], 2), size(gamma_pred));
%% now get the values for gamma_fit

gamma_fit_summary = plotTauAlphaGammaSummary(x_reps_uf1);
%%
gamma_fit_vals = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma, -1), gamma_fit_summary', 'UniformOutput',false));
% set zeros to nan
gamma_fit_vals(gamma_fit_vals == 0) = nan;
gamma_fit_ave = squeeze(mean(gamma_fit_vals, 1, 'omitnan'));
gamma_fit_se = squeeze(std(gamma_fit_vals, [], 1, 'omitnan'));

figure()
nexttile
imagesc(eta_vals_uf, eps_vals_uf, tau_vals, [1 3])
xlabel('\eta')
ylabel('\epsilon')
title('\tau')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, log10(ks_tau_vals), [-3 -1])
xlabel('\eta')
ylabel('\epsilon')
title('log_{10} KS value, \tau fit')
ch = colorbar;
% ch.TickLabels = num2str((10.^ch.Ticks)', '%1.1d');

nexttile
imagesc(eta_vals_uf, eps_vals_uf, alpha_vals, [1 3])
xlabel('\eta')
ylabel('\epsilon')
title('\alpha')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, log10(ks_alpha_vals), [-3 -1])
xlabel('\eta')
ylabel('\epsilon')
title('log_{10} KS value, \alpha fit')
ch = colorbar;


nexttile
imagesc(eta_vals_uf, eps_vals_uf, gamma_pred, [1 1.6])
xlabel('\eta')
ylabel('\epsilon')
title('\gamma_{pred}')
colorbar
% nexttile
% imagesc(eta_vals_uf, eps_vals_uf, se_gamma_pred)
% xlabel('\eta')
% ylabel('\epsilon')
% title('\delta \gamma_{pred}')
% colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, gamma_fit_ave, [1 1.6])
xlabel('\eta')
ylabel('\epsilon')
title('\gamma_{fit}')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, gamma_fit_ave - gamma_pred, [-.2 .2])
xlabel('\eta')
ylabel('\epsilon')
title('\gamma_{fit} - \gamma_{pred}')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, (gamma_fit_ave - gamma_pred)./se_gamma_pred, [-3 3])
xlabel('\eta')
ylabel('\epsilon')
title('(\gamma_{fit} - \gamma_{pred})/(\delta \gamma_{pred})')
colorbar
%% compute rate distribution
n_neur = 1024;
J_i = randn(n_neur, 1);
prob_s_fun = @(h, J1, eta1, eps1) 1./(1 + exp(eta1*J1*h + eps1));


%% Compute, for fixed eta, epsilon, J, h_MLE and SE_h_MLE as h_true varies. 
% fixed parameters
% h_true_vals = linspace(-5, 5, 1001);
% dh = h_true_vals(2) - h_true_vals(1);

h_sig = 1;
prob_h = @(h) normpdf(h, 0, h_sig);
% prob_h = exp(-0.5*h_true_vals.^2/h_sig^2)/sqrt(2*pi*h_sig^2);
%%
r_i = zeros(n_neur, length(eps_vals), length(eta_vals));

for i_eps = 1:length(eps_vals)
    for i_eta = 1:length(eta_vals)
        for i_n = 1:n_neur
            int_fun = @(h) prob_s_fun(h, J_i(i_n), eta_vals(i_eta), eps_vals(i_eps)).*prob_h(h);
%         s_prob = prob_s_fun(h_true_vals, J_i, eta_vals(i_eta), eps_vals(i_eps));
%         r_i(:, i_eps, i_eta) = diag(dh*prob_h)*s_prob';
            r_i(i_n, i_eps, i_eta) = integral(@(x) int_fun(x), -5, 5);
        end
    end
end
%% compute for fixed eps, eta

r_i_test = computeFRatEtaEpsVal(1, 2, randn(128, 1));
histogram(r_i_test)
%% where is a realistic firing rate regime?
time_bin = 0.01;
eta1_vals = [6 5.5  5 4.5  4   4   4   4   4.5 5 5.5 6];
eps1_vals = [-10.8 -10.1 -9.4 -8.7 -8  -9 -10 -11 -11  -11 -11 -11];
t_h = [0 2 4 6 8 10 12 14 16 18 20 22];
% eta1_vals = [6.1 4.25];
% eps1_vals = [-8.1 -6.1];
r_i_ee = zeros(length(J_i), length(eta1_vals));
for i_ee = 1:length(eta1_vals)
    r_i_ee(:, i_ee) = computeFRatEtaEpsVal(eta1_vals(i_ee), -eps1_vals(i_ee), J_i);
end
% r_i2 = computeFRatEtaEpsVal(eta1_vals(2), -eps1_vals(2), J_i);
% r_i3 = computeFRatEtaEpsVal(eta1_vals(3), -eps1_vals(3), J_i);
[~, r_ord] = sort(abs(J_i), 'ascend');
% select a few values on the trajectory to highly in plots below
sel_ee = [1 5 8];
%% 

gamma_pred_err = gamma_pred - gamma_fit_ave;
[et_uf, ep_uf] = meshgrid(eta_vals_uf, eps_vals_uf);
[et1, ep2] = meshgrid(linspace(4, 8, 9), linspace(6, 14, 17));
gpe_int2_full = interp2(et_uf, ep_uf, gamma_pred_err, ...
    et1, ep2, 'nearest');

gpe_int2_traj = interp2(et_uf, ep_uf, gamma_pred_err, ...
    eta1_vals, -eps1_vals, 'linear');

%%
% this needs to match the length of sel_ee
etaeps_colors = [0 0.5 0; 1 1 0; 0.8 0.2 0];
etaeps_markers = {'ks', 'ko', 'k^'};

figure()

nexttile
contourf(eta_vals, -eps_vals, log10(squeeze(median(r_i, 1))/time_bin), ...
    -3:0.5:2)
colorbar
hold on
plot(eta1_vals, eps1_vals, 'wo-', 'linewidth', 2)
for i_ee = 1:length(sel_ee)
    plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
% plot(eta1_vals(2), eps1_vals(2), etaeps_markers{2}, 'markersize', 10,...
%     'linewidth', 2, 'MarkerFaceColor',etaeps_colors(2, :))
% plot(eta1_vals(3), eps1_vals(3), etaeps_markers{3}, 'markersize', 10,...
%     'linewidth', 2, 'MarkerFaceColor',etaeps_colors(3, :))

title('median firing rate')
ylim([-14 -2])

nexttile
imagesc(eta_vals_uf, -eps_vals_uf, gamma_pred_err)
colorbar
hold on
plot(eta1_vals, eps1_vals, 'wo-', 'linewidth', 2)
for i_ee = 1:length(sel_ee)
    plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
set(gca, 'ydir', 'normal')
ylim([-14 -2])


nexttile
hold on
for i_ee = 1:length(sel_ee)
    plot(log10(r_i_ee(r_ord, sel_ee(i_ee))/time_bin), linspace(0, 1, length(r_i1)), etaeps_markers{i_ee}, ...
        'MarkerFaceColor',etaeps_colors(i_ee, :), 'color', 'none')
end
% plot(log10(r_i2(r_ord)/time_bin), linspace(0, 1, length(r_i2)), etaeps_markers{2}, ...
%     'MarkerFaceColor',etaeps_colors(2, :), 'color', 'none')
% plot(log10(r_i3(r_ord)/time_bin), linspace(0, 1, length(r_i3)), etaeps_markers{3}, ...
%     'MarkerFaceColor',etaeps_colors(3, :), 'color', 'none')
axis tight
x_lims = ceil(xlim);
xticks = x_lims(1):0.5:x_lims(2);
set(gca, 'xtick', xticks, 'xticklabel', num2str(10.^xticks', '%1.0g'))
xlabel('firing rate')
ylabel('cdf over cells')

nexttile
plot(t_h, median(r_i_ee)/time_bin, 'linewidth', 1.5)
hold on
for i_ee = 1:length(sel_ee)
    plot(t_h(sel_ee(i_ee)), median(r_i_ee(:, sel_ee(i_ee)))/time_bin, etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
xlabel('time (hours?)')
ylabel('median firing rate, population')
yyaxis right
plot(t_h, gpe_int2_traj, 'linewidth', 1.5)
hold on
for i_ee = 1:length(sel_ee)
    plot(t_h(sel_ee(i_ee)), gpe_int2_traj(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
ylabel('\gamma prediction error')
%% variability across simulations : images
makeMyFigure(30, 8)

nexttile
imagesc(eta_vals_uf, eps_vals_uf, tau_vals, [1 3])
xlabel('\eta')
ylabel('\epsilon')
title('\tau average')
colorbar



nexttile
imagesc(eta_vals_uf, eps_vals_uf, alpha_vals, [1 3])
xlabel('\eta')
ylabel('\epsilon')
title('\alpha average')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf,(gamma_pred - gamma_fit_ave), [-0.2 0.2])
xlabel('\eta')
ylabel('\epsilon')
title('\gamma_{pred} - \gamma_{fit}')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, se_tau_vals, [0 0.1])
xlabel('\eta')
ylabel('\epsilon')
title('\tau SE across simulations')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, se_alpha_vals, [0 0.1])
xlabel('\eta')
ylabel('\epsilon')
title('\alpha SE across simulations')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, se_gamma_pred, [0 0.2])
xlabel('\eta')
ylabel('\epsilon')
title('\gamma_{pred} SE')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, se_tau_vals, [0 0.1])
xlabel('\eta')
ylabel('\epsilon')
title('\tau SE across simulations')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, min_size_vals)
xlabel('\eta')
ylabel('\epsilon')
title('min size')
colorbar

nexttile
imagesc(eta_vals_uf, eps_vals_uf, min_dur_vals)
xlabel('\eta')
ylabel('\epsilon')
title('min duration')
colorbar

% nexttile
% imagesc(eta_vals_uf, eps_vals_uf, se_tau_vals./ave_num_ava_uf1)
% xlabel('\eta')
% ylabel('\epsilon')
% title('\tau SE / # avalanches')
% colorbar
% 
% nexttile
% imagesc(eta_vals_uf, eps_vals_uf, se_alpha_vals.*sqrt(ave_num_ava_uf1))
% xlabel('\eta')
% ylabel('\epsilon')
% title('\alpha SE * sqrt ( # avalanches)')
% colorbar
%%

print(gcf, '-dpdf', 'avalanches_matlab_code/plots/paper_figures/longSim_exponent_variability.pdf')
%% variability: error bars
i_eta = 2;

makeMyFigure(30, 8)

nexttile
hold on
ph = plot(eps_vals_uf, tau_vals(:, i_eta) + se_tau_vals(:, i_eta)*[-2 0 2],'k');
ph(2).LineWidth = 1.5;
xlabel('\epsilon')
ylabel('\tau')
set(gca, 'fontsize', 12, 'color', 'none')


nexttile
ph = plot(eps_vals_uf, alpha_vals(:, i_eta) + se_alpha_vals(:, i_eta)*[-2 0 2],'k');
ph(2).LineWidth = 1.5;
xlabel('\epsilon')
ylabel('\alpha')
set(gca, 'fontsize', 12, 'color', 'none')
title(['estimated exponents and variation across simulations for \eta = ' ...
    num2str(eta_vals_uf(i_eta))])

nexttile
ph = plot(eps_vals_uf, gamma_pred(:, i_eta) + se_gamma_pred(:, i_eta)*[-2 0 2],'k');
ph(2).LineWidth = 1.5;
hold on
ph2 = plot(eps_vals_uf, gamma_fit_ave(:, i_eta)+ gamma_fit_se(:, i_eta)*[-2 0 2], '--m');
ph2(2).LineWidth = 1.5;

xlabel('\epsilon')
ylabel('\gamma_{pred}')
set(gca, 'fontsize', 12, 'color', 'none')

print(gcf, '-dpdf', ...
    ['avalanches_matlab_code/plots/paper_figures/longSim_exponent_variability_eta' ...
    num2str(eta_vals_uf(i_eta)) '.pdf'])
%% plot the error bars only
i_eta = 1;

makeMyFigure(30, 8)

nexttile
hold on
ph = plot(eps_vals_uf, se_tau_vals(:, i_eta),'k', 'LineWidth',1.5);
xlabel('\epsilon')
ylabel('\tau SD')
set(gca, 'fontsize', 12, 'color', 'none')


nexttile
ph = plot(eps_vals_uf, se_alpha_vals(:, i_eta),'k', 'linewidth', 1.5);
xlabel('\epsilon')
ylabel('\alpha SD')
set(gca, 'fontsize', 12, 'color', 'none')
title(['variation in exponents across simulations for \eta = ' ...
    num2str(eta_vals_uf(i_eta))])

nexttile
ph = plot(eps_vals_uf, se_gamma_pred(:, i_eta),'k', 'linewidth', 1.5);
hold on
ph2 = plot(eps_vals_uf, gamma_fit_se(:, i_eta), '--m', 'linewidth', 1.5);
xlabel('\epsilon')
ylabel('\gamma SD')
set(gca, 'fontsize', 12, 'color', 'none')
legend({'prediction', 'fit'})

print(gcf, '-dpdf', ...
    ['avalanches_matlab_code/plots/paper_figures/longSim_exponent_sd_eta' ...
    num2str(eta_vals_uf(i_eta)) '.pdf'])

%% plot as errorbars
figure()
nexttile
plot(eps_vals_uf,(gamma_fit_ave - gamma_pred), 'linewidth', 1.5 )
legend(num2str(eta_vals_uf', 'eta = %1.0f'))

nexttile
plot(eps_vals_uf,(gamma_fit_ave - gamma_pred)./se_gamma_pred, 'linewidth', 1.5 )
legend(num2str(eta_vals_uf', 'eta = %1.0f'))
%% check against the fine samples
% get new gamma_fit values for x_reps_f1
gamma_fit_summary_f1 = plotTauAlphaGammaSummary({x_reps_f1});
gamma_fit_summary_f1 = gamma_fit_summary_f1{1};
%%

makeMyFigure(30, 8);
nexttile
imagesc(eta_vals, eps_vals, sum_stats_f1.tau_values, [1 3])
xlabel('\eta')
ylabel('\epsilon')
title('\tau')
colorbar
% 
% nexttile
% imagesc(eta_vals, eps_vals, (sum_stats_f1.tau_KS_min), [-3 -1])
% xlabel('\eta')
% ylabel('\epsilon')
% title('log_{10} KS value, \tau fit')
% ch = colorbar;
% % ch.TickLabels = num2str((10.^ch.Ticks)', '%1.1d');

nexttile
imagesc(eta_vals, eps_vals, sum_stats_f1.alpha_values, [1 3])
xlabel('\eta')
ylabel('\epsilon')
title('\alpha')
colorbar
% 
% nexttile
% imagesc(eta_vals, eps_vals, (sum_stats_f1.alpha_KS_min), [-3 -1])
% xlabel('\eta')
% ylabel('\epsilon')
% title('log_{10} KS value, \alpha fit')
% ch = colorbar;


% nexttile
% imagesc(eta_vals, eps_vals, gamma_fit_summary_f1.gamma_prediction, [1 1.6])
% xlabel('\eta')
% ylabel('\epsilon')
% title('\gamma_{pred}')
% colorbar
% 
% nexttile
% imagesc(eta_vals, eps_vals, se_gamma_pred)
% xlabel('\eta')
% ylabel('\epsilon')
% title('\delta \gamma_{pred}')
% colorbar

nexttile
imagesc(eta_vals, eps_vals, gamma_fit_summary_f1.fit_gamma, [1 1.6])
xlabel('\eta')
ylabel('\epsilon')
title('\gamma_{fit}')
colorbar

nexttile
imagesc(eta_vals, eps_vals, gamma_fit_summary_f1.gamma_f_minus_p, [-.2 .2])
xlabel('\eta')
ylabel('\epsilon')
title('\gamma_{fit} - \gamma_{pred}')
colorbar

% nexttile
% imagesc(eta_vals, eps_vals, (gamma_fit_ave - gamma_pred)./se_gamma_pred, [-3 3])
% xlabel('\eta')
% ylabel('\epsilon')
% title('(\gamma_{fit} - \gamma_{pred})/(\delta \gamma_{pred})')
% colorbar

print(gcf, '-dpdf', 'avalanches_matlab_code/plots/paper_figures/shortSim_fine_onefield.pdf')

%% let's look at tau fits
tau_fit_arr = cell(length(x_reps_uf1)+1, length(eps_ind_trace));
for i_xra = 1:length(x_reps_uf1)+1
    for jj = 1:length(eta_ind_trace)
        if i_xra > length(x_reps_uf1)
            t_fit = x_reps_f1.all_tau_pfit{eps_ind_trace(jj), eta_ind_trace(jj)};
        else
            t_fit = x_reps_uf1{i_xra}.all_tau_pfit{eps_ind_uf(jj), eta_ind_uf(jj)};
        end
        
        tau_fit_arr{i_xra, jj}.x_min_vals = t_fit.x_min_vals;
        tau_fit_arr{i_xra, jj}.a_vals = t_fit.a_vals;
        tau_fit_arr{i_xra, jj}.se_a_vals = t_fit.se_a_vals;
        tau_fit_arr{i_xra, jj}.ks_stats = t_fit.ks_stats;

    end
end
%%
makeMyFigure(30, 20)

for jj = 1:length(eps_trace)
    nexttile
    hold on
    x_vals = cell2mat(cellfun(@(x) x.x_min_vals', tau_fit_arr(1:end-1, jj)', 'UniformOutput', false));
    y_vals = cell2mat(cellfun(@(x) x.a_vals', tau_fit_arr(1:end-1, jj)', 'UniformOutput',false));
    dy_vals = cell2mat(cellfun(@(x) x.se_a_vals', tau_fit_arr(1:end-1, jj)', 'UniformOutput',false));

    errorbar(x_vals, y_vals, dy_vals, 'k')
    errorbar(tau_fit_arr{end, jj}.x_min_vals, tau_fit_arr{end, jj}.a_vals, ...
        tau_fit_arr{end, jj}.se_a_vals, 'LineWidth',1)
    title(['\epsilon = ' num2str(eps_trace(jj)) ', \eta = ' num2str(eta_trace(jj))])
    xlabel('start point fit')
    ylabel('\tau (SE fit)')
end
for jj = 1:length(eps_trace)
    nexttile
    hold on
    x_vals = cell2mat(cellfun(@(x) x.x_min_vals', tau_fit_arr(1:end-1, jj)', 'UniformOutput', false));
    y_vals = cell2mat(cellfun(@(x) x.ks_stats', tau_fit_arr(1:end-1, jj)', 'UniformOutput',false));

    plot(x_vals, y_vals, 'k')
    plot(tau_fit_arr{end, jj}.x_min_vals, tau_fit_arr{end, jj}.ks_stats, ...
        'LineWidth',1)
    title(['\epsilon = ' num2str(eps_trace(jj)) ', \eta = ' num2str(eta_trace(jj))])
    xlabel('start point fit')
    ylabel('KS stat')
end
%% Alpha fits along the eta-eps trajectory
alpha_fit_arr = cell(length(x_reps_uf1)+1, length(eps_ind_trace));
for i_xra = 1:length(x_reps_uf1)+1
    for jj = 1:length(eta_ind_trace)
        if i_xra > length(x_reps_uf1)
            t_fit = x_reps_f1.all_alpha_pfit{eps_ind_trace(jj), eta_ind_trace(jj)};
        else
            t_fit = x_reps_uf1{i_xra}.all_alpha_pfit{eps_ind_uf(jj), eta_ind_uf(jj)};
        end

        alpha_fit_arr{i_xra, jj}.x_min_vals = t_fit.x_min_vals;
        alpha_fit_arr{i_xra, jj}.a_vals = t_fit.a_vals;
        alpha_fit_arr{i_xra, jj}.se_a_vals = t_fit.se_a_vals;
        alpha_fit_arr{i_xra, jj}.ks_stats = t_fit.ks_stats;

    end
end
%%
makeMyFigure(30, 20)

for jj = 1:length(eps_trace)
    nexttile
    hold on
    x_vals = cell2mat(cellfun(@(x) x.x_min_vals', alpha_fit_arr(1:end-1, jj)', 'UniformOutput', false));
    y_vals = cell2mat(cellfun(@(x) x.a_vals', alpha_fit_arr(1:end-1, jj)', 'UniformOutput',false));
    dy_vals = cell2mat(cellfun(@(x) x.se_a_vals', alpha_fit_arr(1:end-1, jj)', 'UniformOutput',false));

    errorbar(x_vals, y_vals, dy_vals, 'k')
    plot(x_vals(:, 1), mean(y_vals, 2), 'k', 'LineWidth', 1.5)
    errorbar(alpha_fit_arr{end, jj}.x_min_vals, alpha_fit_arr{end, jj}.a_vals, ...
        alpha_fit_arr{end, jj}.se_a_vals, 'LineWidth',1)
    title(['\epsilon = ' num2str(eps_trace(jj)) ', \eta = ' num2str(eta_trace(jj))])
    xlabel('start point fit')
    ylabel('\alpha (SE fit)')
end

for jj = 1:length(eps_trace)
    nexttile
    hold on
    x_vals = cell2mat(cellfun(@(x) x.x_min_vals', alpha_fit_arr(1:end-1, jj)', 'UniformOutput', false));
    y_vals = cell2mat(cellfun(@(x) x.ks_stats', alpha_fit_arr(1:end-1, jj)', 'UniformOutput',false));

    plot(x_vals, y_vals, 'k')
    plot(x_vals(:, 1), mean(y_vals, 2), 'k', 'LineWidth', 1.5)
    plot(alpha_fit_arr{end, jj}.x_min_vals, alpha_fit_arr{end, jj}.ks_stats, ...
        'LineWidth',1)
    title(['\epsilon = ' num2str(eps_trace(jj)) ', \eta = ' num2str(eta_trace(jj))])
    xlabel('start point fit')
    ylabel('KS stat')
end
% h_Fig = makeMyFigure(30, 30);
% tiledlayout(3,2)
% for iG = 1:5, h_spG(iG) = nexttile; end
% 
% %%\
% gamma_summary_UF1 = plotTauAlphaGammaSummary(x_reps_uf1, {'A', 'B', 'C', 'D', 'E'}, ...
%         h_Fig, h_spG);
