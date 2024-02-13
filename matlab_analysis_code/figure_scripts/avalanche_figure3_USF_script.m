%% Consolidate fits - SF
% res_str = consolidateAlphaTauFits(x_reps_sf1, eps_trace(1), eta_trace(1), 'alpha');
cmap_nan = [1 1 1 ; parula];
time_cmap_nan = [0.5*ones(1,3); parula(50)];

N_sf1 = 128;
%% LOAD: 1-field Ultra Small-Network Fine sweep with replicates (1USR): plotting firing rate, avalanche counts, and field values

rep_list = {'A', 'B', 'C', 'D', 'E', 'F'};
eta_list_sf1 = 1:2;
eps_list_sf1 = -2:-0.5:-5.5; % [ -6.0, -8.0, -10.0, -12.0, -14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[x_info_sf1, num_ava_sf1] = loadSimulationRunResults('run_f1smallultrafinesweep', rep_list);

% Summaries from avalanche analysis
results_dir = 'avalanches_matlab_code/analysis_results/fields1smallultrafinesweep/';
x_reps_sf1 = cell(length(rep_list), 1);
for ii = 1:length(x_reps_sf1)
    rep_dir = [results_dir 'rep' rep_list{ii} '/' ];
    rep_fn = dir([rep_dir 'ava_decade_analysis*.mat']);
    [~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');
    x_reps_sf1{ii} = load([rep_dir rep_fn(date_ord(1)).name]);
end
sum_stats_sf1 = pullSummaryFitInfo(x_reps_sf1, eta_list_sf1, eps_list_sf1);

%% Consolidate fits - SF
% res_str = consolidateAlphaTauFits(x_reps_sf1, eps_trace(1), eta_trace(1), 'alpha');

all_alpha_SF = cell(length(eps_list_sf1), length(eta_list_sf1));
all_tau_SF = cell(length(eps_list_sf1), length(eta_list_sf1));
for ii = 1:length(eps_list_sf1)
    for jj = 1:length(eta_list_sf1)
        try
        all_alpha_SF{ii, jj} = consolidateAlphaTauFits(x_reps_sf1, ...
            eps_list_sf1(ii), eta_list_sf1(jj), 'alpha');
        end
        try
        all_tau_SF{ii, jj} = consolidateAlphaTauFits(x_reps_sf1, ...
            eps_list_sf1(ii), eta_list_sf1(jj), 'tau');
        end
    end
end


%% plots 
has_entry = cellfun(@(x) ~isempty(x), all_alpha_SF);
% alpha information
alpha_vals = nan(size(all_alpha_SF));
alpha_vals(has_entry) = cellfun(@(x) x.exp_val_est, all_alpha_SF(has_entry));

se_alpha_vals = nan(size(all_alpha_SF));
se_alpha_vals(has_entry) = cellfun(@(x) x.se_exp_val_est, all_alpha_SF(has_entry));

ks_alpha_vals = nan(size(all_alpha_SF));
ks_alpha_vals(has_entry) = cellfun(@(x) median(x.ks_at_min), all_alpha_SF(has_entry));

min_size_vals = nan(size(all_alpha_SF));
min_size_vals(has_entry) = cellfun(@(x) x.x_min_est, all_alpha_SF(has_entry));

% tau information
tau_vals = nan(size(all_tau_SF));
tau_vals(has_entry) = cellfun(@(x) x.exp_val_est, all_tau_SF(has_entry));

se_tau_vals = nan(size(all_tau_SF));
se_tau_vals(has_entry) = cellfun(@(x) x.se_exp_val_est, all_tau_SF(has_entry));

ks_tau_vals = nan(size(all_tau_SF));
ks_tau_vals(has_entry) = cellfun(@(x) median(x.ks_at_min), all_tau_SF(has_entry));

min_dur_vals = nan(size(all_tau_SF));
min_dur_vals(has_entry) = cellfun(@(x) x.x_min_est, all_tau_SF(has_entry));

mult_ave_num_ava = squeeze(round(exp(mean(log(num_ava_sf1)))));
mult_ave_num_ava(mult_ave_num_ava == 0) = nan;
%% get surrogate KS stats

surr_ks_alpha = 0*alpha_vals;
surr_ks_alpha(has_entry) = cellfun(@(x) median(prctile(x.ks_surr, 95, 2)), ...
    all_alpha_SF(has_entry));

surr_ks_tau = 0*tau_vals;
surr_ks_tau(has_entry) = cellfun(@(x) median(prctile(x.ks_surr, 95, 2)), ...
    all_tau_SF(has_entry));

gamma_pred = (alpha_vals - 1)./(tau_vals - 1);
%% gamma pred from alpha, tau at ensemble min size
gamma_pred_at_min = cell(size(all_alpha_SF));
gamma_pred_at_min(has_entry) = cellfun(@(x, y) (x.exp_at_min-1)./(y.exp_at_min - 1), ...
    all_alpha_SF(has_entry), all_tau_SF(has_entry), 'UniformOutput',false);
gamma_pred_at_min(~has_entry) = cellfun(@(x) nan(10, 1), ...
    all_alpha_SF(~has_entry), 'UniformOutput',false);

mean_gamma_pred_atmin = cellfun(@mean, gamma_pred_at_min);
se_gamma_pred_atmin = cellfun(@std, gamma_pred_at_min);

% alpha at min
alpha_at_min = cell(size(all_alpha_SF));
alpha_at_min(has_entry) = cellfun(@(x) x.exp_at_min, ...
    all_alpha_SF(has_entry), 'UniformOutput',false);
alpha_at_min(~has_entry) = cellfun(@(x) nan(10, 1), ...
    all_alpha_SF(~has_entry), 'UniformOutput',false);

mean_alpha_atmin = cellfun(@mean, alpha_at_min);
se_alpha_atmin = cellfun(@std, alpha_at_min);

% tau at min
tau_at_min = cell(size(all_tau_SF));
tau_at_min(has_entry) = cellfun(@(x) x.exp_at_min, ...
    all_tau_SF(has_entry), 'UniformOutput',false);
tau_at_min(~has_entry) = cellfun(@(x) nan(10, 1), ...
    all_tau_SF(~has_entry), 'UniformOutput',false);

mean_tau_atmin = cellfun(@mean, tau_at_min);
se_tau_atmin = cellfun(@std, tau_at_min);


%% estimate gamma_pred_se
num_draws = 1000;
alpha_val_sur = alpha_vals(:) + repmat(se_alpha_vals(:), 1, num_draws).*randn(numel(alpha_vals), num_draws);
tau_val_sur = tau_vals(:) + repmat(se_tau_vals(:), 1, num_draws).*randn(numel(alpha_vals), num_draws);
gamma_pred_sur = (alpha_val_sur - 1)./(tau_val_sur - 1);

se_gamma_pred = reshape(std(gamma_pred_sur, [], 2), size(gamma_pred));
%% now get the values for gamma_fit
if ~exist('avalanches_matlab_code/gamma_fit_ultrasmallfine_20221212.mat', 'file')
    
gamma_fit_summary_SF = plotTauAlphaGammaSummary(x_reps_sf1);
% save gamma fit info

save(['avalanches_matlab_code/gamma_fit_ultrasmallfine_' ...
    datestr(now, 'yyyymmdd')], 'gamma_fit_summary_SF')
else
load('avalanches_matlab_code/gamma_fit_ultrasmallfine_20221212.mat', ...
    'gamma_fit_summary_SF')

end
%% Unnecessary plot block: unformatted eta-eps alpha and tau

gamma_fit_vals = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma, -1), gamma_fit_summary_SF', 'UniformOutput',false));
% set zeros to nan
gamma_fit_vals(gamma_fit_vals == 0) = nan;
gamma_fit_ave = squeeze(mean(gamma_fit_vals, 1, 'omitnan'));
gamma_fit_se = squeeze(std(gamma_fit_vals, [], 1, 'omitnan'));

figure()
nexttile
imagesc(eta_list_sf1, eps_list_sf1, tau_vals, [1 3])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\tau')
colorbar
% set(gca, 'ydir', 'normal')

nexttile
imagesc(eta_list_sf1, eps_list_sf1, log10(ks_tau_vals), [-3 -1])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('log_{10} KS value, \tau fit')
ch = colorbar;
% ch.TickLabels = num2str((10.^ch.Ticks)', '%1.1d');

nexttile
imagesc(eta_list_sf1, eps_list_sf1, alpha_vals, [1 3])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\alpha')
colorbar

nexttile
imagesc(eta_list_sf1, eps_list_sf1, log10(ks_alpha_vals), [-3 -1])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('log_{10} KS value, \alpha fit')
ch = colorbar;


nexttile
imagesc(eta_list_sf1, eps_list_sf1, gamma_pred, [1 1.6])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\gamma_{pred}')
colorbar
% nexttile
% imagesc(eta_vals_uf, eps_vals_uf, se_gamma_pred)
% xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
% title('\delta \gamma_{pred}')
% colorbar

nexttile
imagesc(eta_list_sf1, eps_list_sf1, gamma_fit_ave, [1 1.6])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\gamma_{fit}')
colorbar

nexttile
imagesc(eta_list_sf1, eps_list_sf1, gamma_fit_ave - gamma_pred, [-.2 .2])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\gamma_{fit} - \gamma_{pred}')
colorbar

nexttile
imagesc(eta_list_sf1, eps_list_sf1, (gamma_fit_ave - gamma_pred)./se_gamma_pred, [-2 2])
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('(\gamma_{fit} - \gamma_{pred})/(\delta \gamma_{pred})')
colorbar

%% 
eps_vals_sf1 = -eps_list_sf1;
eta_vals_sf1 = eta_list_sf1;
figure()
nexttile
ph = plot(eps_vals_sf1,(gamma_fit_ave - gamma_pred), 'linewidth', 1.5 );
assignColorsToLines(ph, parula(length(ph)))
legend(num2str(eta_vals_sf1', 'eta = %1.0f'))

nexttile
ph = plot(eps_vals_sf1,se_gamma_pred, 'linewidth', 1.5 );
assignColorsToLines(ph, parula(length(ph)))

legend(num2str(eta_vals_sf1', 'eta = %1.0f'))
ylabel('\gamma_{pred} uncertainty')
%% compute firing rates over eta-epsilon space for random J_i
% N.B. The simulation firing rates are in x_info_sf1
J_i = randn(128, 1);
n_neur = length(J_i);
eps_vals_full = 1.5:0.5:max(-eps_list_sf1);
eta_vals_full = 0:1:max(eta_list_sf1);
r_i = zeros(n_neur, length(eps_vals_full), length(eta_vals_full));

for i_eps = 1:length(eps_vals_full)
    for i_eta = 1:length(eta_vals_full)
        % value of epsilon should be positive here
        r_i(:, i_eps, i_eta) = computeFRatEtaEpsVal(eta_vals_full(i_eta), eps_vals_full(i_eps), J_i);

    end
end

%% Get curve for firing rate - avalanche count
mean_fr_ee = squeeze(mean(r_i(:, 2:end, 2:end), 1));

makeMyFigure(18, 16) 
nexttile
imagesc(eta_list_sf1, eps_list_sf1, mult_ave_num_ava)
set(gca, 'YDir', 'normal')
axis square
set(gca, 'fontsize', 14)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
ch = colorbar;
title('Number of avalanches')

nexttile
imagesc(eta_list_sf1, eps_list_sf1, mean_fr_ee)
set(gca, 'YDir', 'normal')
axis square
set(gca, 'fontsize', 14)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
ch = colorbar;
title('Mean population firing rate (FR)')

nexttile
ph = plot(mean_fr_ee, mult_ave_num_ava, 'o-', 'linewidth', 2);
xlabel('FR (per bin)')
ylabel('N_{ava}')
assignColorsToLines(ph, parula(length(ph)));
axis square
set(gca, 'fontsize', 14, 'color', 'none')
legend(ph, num2str(eta_list_sf1', 'eta = %1.0f'))

nexttile
scaled_sum_fr_ee = log10(mean_fr_ee*diag(1./eta_list_sf1));
scaled_num_ava = log10(mult_ave_num_ava*diag(eta_list_sf1));
ph = plot(scaled_sum_fr_ee, scaled_num_ava, 'o', 'linewidth', 2);
xlabel('log_{10}( FR / \eta )')
ylabel('log_{10} \eta N_{ava}')
assignColorsToLines(ph, parula(length(ph)));
% legend(ph, num2str(eta_list_sf1', 'eta = %1.0f'))
axis square
set(gca, 'fontsize', 14, 'color', 'none')



% print(gcf, '-dpdf', 'avalanches_matlab_code/plots/paper_figures/poster_figures/fr_ava_etaeps.pdf')
% print(gcf, '-dsvg', 'avalanches_matlab_code/plots/paper_figures/poster_figures/fr_ava_etaeps.svg')


%% Define trajectory in eta-epsilon space and calculate firing rates, etc. 
time_bin = 0.03;
% eta1_vals = [7    6   5  4  4   4   4   5  5.5   6  6.5  7];
% eps1_vals = [-10 -9  -8 -6 -8  -9 -10 -10 -10  -10 -10 -10];

% path 1: starting at low-information point
% eta1_vals = [ 8    7   5   4   4   4   5   6   7   8];
% eps1_vals = [-10  -9  -7  -6  -8  -10  -10  -10  -10  -10];

% path 2: starting at higher-information point
eta1_vals = [ 8   7    5     4   4     4   5   6   7   8];
eps1_vals = [-6  -5   -4    -3  -4.5  -6  -6  -6  -6  -6];

% eta1_vals = [ 8    7   6    5  4   4   4   5   6   7   8];
% eps1_vals = [-4 -3.5  -3 -2.5 -2  -3  -4  -4  -4  -4  -4];
t_h = 2*(1:length(eta1_vals)); %[0 2 4 6 8 10 12 14 16 18 20 22];
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
sel_ee = [1 4 6]; % path 1 and 2

%% 

gamma_pred_err = gamma_pred - gamma_fit_ave;
[et_sf, ep_sf] = meshgrid(eta_list_sf1, -eps_list_sf1);
[et1, ep2] = meshgrid(linspace(4, 8, 9), linspace(6, 14, 17));
gpe_int2_full = interp2(et_sf, ep_sf, gamma_pred_err, ...
    et1, ep2, 'nearest');

gpe_int2_traj = interp2(et_sf, ep_sf, gamma_pred_err, ...
    eta1_vals, -eps1_vals, 'linear');
%% binary: scaling or no [what happened with epsilon = 6?]
norm_gpe = double(abs(gamma_pred - gamma_fit_ave)./se_gamma_pred < 1);
% [et_sf, ep_sf] = meshgrid(eta_vals_sf1, eps_vals_sf1);
[et1, ep2] = meshgrid(linspace(1, 10, 10), linspace(2, 14, 7));
norm_gpe_int2_full = interp2(et_sf, ep_sf, norm_gpe, ...
    et1, ep2, 'nearest');

norm_gpe_int2_traj = interp2(et_sf, ep_sf, norm_gpe, ...
    eta1_vals, -eps1_vals, 'nearest');
%% load information calculation

i_calc = load('info_workspace_20220930.mat', 'eps_vals', 'eta_vals', 'mask_info_vals', 'r_i', 'j_vals', 'I_si_h10', 'J_i');
i_calc.r_i_1k = i_calc.r_i;
i_calc.r_i = zeros(length(i_calc.J_i), length(i_calc.eps_vals), length(i_calc.eta_vals));
for ii = 1:length(i_calc.eps_vals)
    for jj = 1:length(i_calc.eta_vals)
        i_calc.r_i(:, ii, jj) = computeFRatEtaEpsVal(i_calc.eta_vals(jj), i_calc.eps_vals(ii), i_calc.J_i);
    end
end
%%
i_calc.info_rate = i_calc.I_si_h10./squeeze(sum(i_calc.r_i, 1));
%% Compute log10(mean(r_i)/eta)
% [* for 10 cells, as in the info calculations, but think about this] 
scaled_FR_surf = (squeeze(mean(i_calc.r_i, 1))*diag(1./i_calc.eta_vals));

%% interpolate paths
[ic_eta, ic_eps] = meshgrid(i_calc.eta_vals, i_calc.eps_vals);
i_rate_path = interp2(ic_eta, ic_eps, i_calc.info_rate, eta1_vals, -eps1_vals);
info_path =  interp2(ic_eta, ic_eps, i_calc.I_si_h10, eta1_vals, -eps1_vals);
%% Figure 4: information images + avalanche boundaries
% Get curve for firing rate - avalanche count
mean_fr_ee = squeeze(mean(r_i(:, 2:end, 2:end), 1));

makeMyFigure(44, 16); 

%%%%%%%%% part 1: how the # of avalanches depends on eta, pop. fr
nexttile
imagesc(eta_list_sf1, eps_list_sf1, mult_ave_num_ava)
set(gca, 'YDir', 'normal')
axis square
set(gca, 'fontsize', 14, 'ytick', -14:4:-2)

xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('N_{ava}')
ch = colorbar;
ch.Label.String = 'Number of avalanches';

nexttile
imagesc(eta_list_sf1, eps_list_sf1, mean_fr_ee)
set(gca, 'YDir', 'normal')
axis square
set(gca, 'fontsize', 14, 'ytick', -14:4:-2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('FR')
ch = colorbar;
ch.Ticks = 0:0.1:0.3;
ch.Label.String = {'Pop. mean firing rate';'(per bin)'};

sh = nexttile;
ph = plot(mean_fr_ee, mult_ave_num_ava, 'o-', 'linewidth', 2);
xlabel('FR (mean, per bin)')
ylabel('N_{ava}')
assignColorsToLines(ph, parula(length(ph)));
axis square
set(gca, 'fontsize', 14, 'color', 'none')
% legend(ph, num2str(eta_list_sf1', 'eta = %1.0f'), 'Location','east')
colormap(sh, parula(length(eta_list_sf1)));
ch = colorbar;
ch.Ticks = [0 1];
ch.TickLabels = {'\eta = 1', '\eta = 10'};

nexttile
cv_levels = [0.5 1.5 2.5 3.5]; %N_sf1*10.^(-4:.4:-1.6); % contour levels
scaled_sum_fr_ee = N_sf1*(mean_fr_ee*diag(1./eta_list_sf1));
scaled_num_ava = (mult_ave_num_ava*diag(eta_list_sf1));
ph = plot(scaled_sum_fr_ee, scaled_num_ava, 'o', 'linewidth', 2);
% xlabel('log_{10}( FR / \eta )')
xlabel('\eta^{-1} N_{pop} FR')

ylabel('\eta N_{ava}')
assignColorsToLines(ph, parula(length(ph)));
y_lim = ylim;
hold on
plot([1; 1]*cv_levels, y_lim', 'k-')
set(gca, 'fontsize', 14, 'color', 'none')
% legend(ph, num2str(eta_list_sf1', 'eta = %1.0f'), 'Location','westoutside')
axis square
%%%%%%%%%%%%%%%%% part 2: scaling, information, cutoffs, etc. 

sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err), [-0.2 .2])

% imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err)./se_gamma_pred, [0 2])
hold on 
% [~, cr_h] = contour(i_calc.eta_vals, -i_calc.eps_vals, log10(med1Hz_bin), ...
%     dt_tics, 'LineColor',[1 0 0], 'LineWidth',1.5);
hold on
contour(i_calc.eta_vals, -i_calc.eps_vals, N_sf1*scaled_FR_surf, ...
    cv_levels, 'LineColor',[0 0 0], 'ShowText','on', ...
    'LineWidth', 1);

% contour(eta_vals_full, -eps_vals_full, log10(squeeze(mean(r_i, 1))/time_bin), ...
%     dt_tics, 'linecolor', [1 1 1]*0.8)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; usa(50).^.4];
colormap(sh, cmap_nan)
colorbar
% hold on
% plot(eta1_vals, eps1_vals, 'ko-', 'linewidth', 1)
% for i_ee = 1:length(sel_ee)
%     plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'fontsize', 14)
axis square
ylim([-15 -1])
title('\gamma_{pred} - \gamma_{fit}')
%%%%%%%%%%%%%%%% minimum size cutoff
sz_lim = [0 2];
sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(min_size_vals), sz_lim);
colormap(sh, time_cmap_nan)
hold on
contour(i_calc.eta_vals, -i_calc.eps_vals, N_sf1*scaled_FR_surf, ...
    cv_levels, 'LineColor',[0 0 0], 'ShowText','on', ...
    'LineWidth', 1);

ch_size = colorbar;
ch_size.Ticks = linspace(sz_lim(1), sz_lim(2), 3);
ch_size.TickLabels = num2str(10.^ch_size.Ticks', '%1.0f');
axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
set(gca, 'fontsize', 14)
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
ylim([-15 -1])

title('S_{min} (\tau fit)')

%%%%%%%%%%%%%%%% information
nexttile 
imagesc(i_calc.eta_vals, -i_calc.eps_vals, i_calc.I_si_h10/log(2), [0 2]);
ch = colorbar;
ch.Label.String = "I({s_i}, h) (bits)";
hold on
contour(i_calc.eta_vals, -i_calc.eps_vals, N_sf1*scaled_FR_surf, ...
    cv_levels, 'LineColor',[0 0 0], 'ShowText','on', ...
    'LineWidth', 1);
set(gca, 'YDir', 'normal')
axis square
set(gca, 'fontsize', 14, 'ytick', -14:4:-2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
ylim([-15 -1])
title('Information')
text(0, -18, '- log_{10} (FR/\eta)', 'FontWeight', 'normal', 'FontSize', 14)

%%%%%%%%%%%%%%%% information rate
nexttile 
infr_levels = 2.^(-1:1:2);
imagesc(i_calc.eta_vals, -i_calc.eps_vals, log2((i_calc.info_rate)/log(2)), [-1 2]);
ch = colorbar;
ch.Label.String = "I({s_i}, h)/spike";
ch.Ticks = log2(infr_levels);
ch.TickLabels = num2str(infr_levels');
hold on
contour(i_calc.eta_vals, -i_calc.eps_vals, N_sf1*scaled_FR_surf, ...
    cv_levels, 'LineColor',[0 0 0], 'ShowText','on', ...
    'LineWidth', 1);
contour(i_calc.eta_vals, -i_calc.eps_vals, ((i_calc.info_rate)/log(2)), ...
    infr_levels, 'LineColor',1 +[0 0 0], 'ShowText','on', ...
    'LineWidth', 1);
set(gca, 'YDir', 'normal')
axis square
set(gca, 'fontsize', 14, 'ytick', -14:4:-2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
ylim([-15 -1])
title('Information rate')
% text(0, -18, '- log_{10} (FR/\eta)', 'FontWeight', 'normal', 'FontSize', 14)



% % %%%%%%%%%%%%%%%%%%%% experiment length
% % exp_T = 1e7*time_bin/3600;
% % index_sim_ava = mult_ave_num_ava(6, 4);
% % sh = nexttile;
% % % time_tics = (-1:1:1) + log10(exp_T); % adaptive
% % time_tics = log10([3 30 300]); % good for time_bin = 0.01
% % imagesc(eta_list_sf1, eps_list_sf1, log10(exp_T*index_sim_ava./mult_ave_num_ava), time_tics([1 end]));
% % hold on
% % contour(i_calc.eta_vals, -i_calc.eps_vals, N_sf1*scaled_FR_surf, ...
% %     cv_levels, 'LineColor',[0 0 0], 'ShowText','on', ...
% %     'LineWidth', 1);
% % 
% % % contour(eta_vals_full, -eps_vals_full, log10(med1Hz_bin), ...
% % %     [0.5 1], 'LineColor',[1 0 0], 'LineWidth',1.5);
% % % [~, cr_h] = contour(i_calc.eta_vals, -i_calc.eps_vals, log10(med1Hz_bin), ...
% % %     dt_tics, 'LineColor',[1 0 0], 'LineWidth',1.5);
% % 
% % time_cmap_nan = [0.5*ones(1,3); parula(50)];
% % colormap(sh, time_cmap_nan)
% % ch_time = colorbar;
% % ch_time.Ticks = time_tics;
% % ch_time.TickLabels = num2str(10.^ch_time.Ticks', '%1.0f h');
% % axis square
% % set(gca, 'fontsize', 14)
% % set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
% % ylim([-15 -1])
% % 
% % title('Est. Req. Simulation Time')


% print(gcf, '-dpdf', 'avalanches_matlab_code/plots/paper_figures/poster_figures/info_fr_ava_etaeps.pdf')
% print(gcf, '-dsvg', 'avalanches_matlab_code/plots/paper_figures/poster_figures/info_fr_ava_etaeps.svg')
% print(gcf, '-dpng', 'avalanches_matlab_code/plots/paper_figures/poster_figures/info_fr_ava_etaeps.png')



%% Detailed trajectory plots
% this needs to match the length of sel_ee
etaeps_colors = cool(3); %[0 0.5 0; 1 1 0; 0.8 0.2 0];
etaeps_markers = {'ks', 'ko', 'k^'};

makeMyFigure(40, 18);
% makeMyFigure(32, 16);

nexttile
contourf(eta_vals_full, -eps_vals_full, log10(squeeze(mean(r_i, 1))/time_bin), ...
    -3:0.5:2)
colorbar
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
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

title('Mean FR')
axis square
set(gca, 'fontsize', 14)
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
ylim([-14 0])

%%%%%%%%%%%%%%%%%%%% number of avalanches

sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(mult_ave_num_ava), [3.5 6.5]);
ch = colorbar;
ch.Ticks = 4:6;
ch.TickLabels = num2str(ch.Ticks', '10^%1.0f');
axis square
set(gca, 'fontsize', 14)
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
ylim([-15 -1])

title('N_{ava}')

%%%%%%%%%%%%%%%% need some sort of goodness-of-alpha fit measure

%%%%%%%%%%%%%%%%%%%%% gamma_pred minus gamma_fit
sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err), [-0.2 .2])

% imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err)./se_gamma_pred, [0 2])
hold on 
contour(eta_vals_full, -eps_vals_full, log10(squeeze(median(r_i, 1))/time_bin), ...
    -3:0.5:2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_usa_nan = [0.5*ones(1,3) ; usa(50).^.4];
colormap(sh, cmap_usa_nan)
colorbar
hold on
plot(eta1_vals, eps1_vals, 'ko-', 'linewidth', 1)
for i_ee = 1:length(sel_ee)
    plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'fontsize', 14)
axis square
ylim([-15 -1])
title('\gamma_{pred} - \gamma_{fit}')

%%%%%%%%%%%%%%%%%%% std dev across runs for gamma_pred 
sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, se_gamma_pred, [0 .2])

% imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err)./se_gamma_pred, [0 2])
hold on 
contour(eta_vals_full, -eps_vals_full, log10(squeeze(median(r_i, 1))/time_bin), ...
    -3:0.5:2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
colormap(sh, cmap_nan)
colorbar
hold on
plot(eta1_vals, eps1_vals, 'wo-', 'linewidth', 2)
for i_ee = 1:length(sel_ee)
    plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'fontsize', 14)
axis square
ylim([-15 -1])
title('SE \gamma_{pred}')


nexttile
imagesc(i_calc.eta_vals, i_calc.eps_vals, i_calc.I_si_h10, [0 1.5])
ylim([0 14])
% hold on 
% contour(eta_vals_full, eps_vals_full, log10(squeeze(median(r_i, 1))/time_bin), ...
%     -3:0.5:2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('I(s_{1,.. 10}, h)')
colorbar
hold on
plot(eta1_vals, -eps1_vals, 'wo-', 'linewidth', 2)
for i_ee = 1:length(sel_ee)
    plot(eta1_vals(sel_ee(i_ee)), -eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
axis square
set(gca, 'fontsize', 14)

nexttile
[et1,ep1] = meshgrid(i_calc.eta_vals, i_calc.eps_vals);
[et1q,ep1q] = meshgrid(eta_vals_full, eps_vals_full);

i_per_spk = interp2(et1, ep1, i_calc.mask_info_vals, ...
    et1q, ep1q)./squeeze(mean(r_i, 1)/time_bin);
info_raw = interp2(et1, ep1, i_calc.mask_info_vals, ...
    et1q, ep1q);
info_10 = interp2(et1, ep1, i_calc.I_si_h10, ...
    et1q, ep1q);

imagesc(i_calc.eta_vals, i_calc.eps_vals, i_calc.info_rate, [0 4]), 
ylim([0 14])
colorbar
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('I_N / mean pop. FR')
hold on
plot(eta1_vals, -eps1_vals, 'wo-', 'linewidth', 2)
for i_ee = 1:length(sel_ee)
    plot(eta1_vals(sel_ee(i_ee)), -eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
axis square
set(gca, 'fontsize', 14)
nexttile
hold on
for i_ee = 1:length(sel_ee)
    plot(log10(r_i_ee(r_ord, sel_ee(i_ee))/time_bin), linspace(0, 1, n_neur), etaeps_markers{i_ee}, ...
        'MarkerFaceColor',etaeps_colors(i_ee, :), 'color', 'none')
end

% plot(log10(r_i2(r_ord)/time_bin), linspace(0, 1, length(r_i2)), etaeps_markers{2}, ...
%     'MarkerFaceColor',etaeps_colors(2, :), 'color', 'none')
% plot(log10(r_i3(r_ord)/time_bin), linspace(0, 1, length(r_i3)), etaeps_markers{3}, ...
%     'MarkerFaceColor',etaeps_colors(3, :), 'color', 'none')
axis tight
x_lims = ceil(xlim);
xticks = x_lims(1):0.5:x_lims(2);
set(gca, 'fontsize', 14, 'color', 'none', ...
    'xtick', xticks, 'xticklabel', num2str(10.^xticks', '%1.0g'))
axis square
xlabel('firing rate')
ylabel('cdf over cells')

sh = nexttile;
plot(t_h, median(r_i_ee)/time_bin, 'w', 'linewidth', 1.5)
y_lims = ylim;
y_lims(1) = 0;
hold on
imagesc(t_h, y_lims, [1 ; 1]*norm_gpe_int2_traj)
ylim(y_lims)
xlim(t_h([1 end]))
colormap(sh, parula(2));
ch = colorbar;
ch.Ticks = [0 1];
ch.TickLabels = {'no scaling', 'scaling'};
plot(t_h, median(r_i_ee)/time_bin,'r', 'linewidth', 1.5)

hold on
for i_ee = 1:length(sel_ee)
    plot(t_h(sel_ee(i_ee)), median(r_i_ee(:, sel_ee(i_ee)))/time_bin, etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
xlabel('time (hours?)')
ylabel('median firing rate, population')
axis square
yyaxis right
plot(t_h, info_path, 'linewidth', 1.5)
plot(t_h, i_rate_path, 'linewidth', 1.5)
ylabel({'Information rate';'(solid, per bin; dashed, per spike*)'})
% plot(t_h, gpe_int2_traj, 'linewidth', 1.5)
% hold on
% for i_ee = 1:length(sel_ee)
%     plot(t_h(sel_ee(i_ee)), gpe_int2_traj(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
% ylabel('\gamma prediction error')
axis square
set(gca, 'fontsize', 14)

print(gcf, '-dsvg', ['avalanches_matlab_code/plots/paper_figures/info_rate_scaling_plots_' num2str(yyyymmdd(datetime))])
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/info_rate_scaling_plots_' num2str(yyyymmdd(datetime)) '.pdf'])
%% Simple: information MLE (larger population)
i_calc.mask_info_vals
makeMyFigure(10, 8);
imagesc(i_calc.eta_vals, i_calc.eps_vals, i_calc.mask_info_vals)
ylim([0 14])
% hold on 
% contour(eta_vals_full, eps_vals_full, log10(squeeze(median(r_i, 1))/time_bin), ...
%     -3:0.5:2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('I_{est}(h_{MLE}, h)')
set(gca, 'fontsize', 14)
colorbar
print(gcf, '-dsvg', ['avalanches_matlab_code/plots/paper_figures/mle_info_rate_plots_' num2str(yyyymmdd(datetime))])


%% 
makeMyFigure(30, 10);
nexttile
ph = plot(i_per_spk(2:end, 2:end)', gamma_pred' - gamma_fit_ave', 'o', ...
    'linewidth', 2); 
xlabel('info/mean FR')
ylabel('\gamma_{pred} - \gamma_{fit}')
assignColorsToLines(ph, cool(length(ph)));
set(gca, 'color','none', 'fontsize', 12)
nexttile
ph = plot(info_raw(2:end, 2:end)', gamma_pred' - gamma_fit_ave', 'o', ...
    'linewidth', 2); 
xlabel('info (MLE)')
assignColorsToLines(ph, cool(length(ph)));
set(gca, 'color','none', 'fontsize', 12)
set(gca, 'color','none', 'fontsize', 12)
nexttile
ph = plot(info_10(2:end, 2:end)', gamma_pred' - gamma_fit_ave', 'o', ...
    'linewidth', 2); 
xlabel('I_{10} (exact)')
assignColorsToLines(ph, cool(length(ph)));
set(gca, 'color','none', 'fontsize', 12)
legend(ph, num2str(eps_list_sf1', 'eps = %1.0f'))
%%
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/info_crackling_' datestr(now, 'yyyymmdd')])
%% plots: Fig 3 / supplement? 
se_lim = [0 0.1];
at_lim = [1.9 2.8];
gam_lim = [1 1.4];

% eta_ind_ex = [2 2 5 7 8];
% eps_ind_ex = [3 7 4 2 7];

eta_ind_ex = [1 2 5 9];
eps_ind_ex = [2 3 5 7];

% eta_ind_ex2 = [1 3 5 7 9];
% eps_ind_ex2 = [4 4 4 4 4];
% 
ex_symbols = {'^', '*', 's', 'o', 'v'};
eps_colors = gray(length(eps_ind_ex));
eps_colors(:, 2) = 0;
makeMyFigure(30, 20);
tiledlayout(3, 4, 'Padding','tight', 'TileSpacing','tight')

% makeMyFigure(50, 14);
colormap(cmap_nan)


%%%%%%% alpha mean
nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, mean_alpha_atmin, at_lim)
hold on
for i_ex = 1:length(eta_ind_ex)
    for j_ex = 1:length(eps_ind_ex)
        plot(eta_vals_sf1(eta_ind_ex(i_ex)), eps_vals_sf1(eps_ind_ex(j_ex)), ex_symbols{i_ex}, ...
            'color', eps_colors(j_ex, :), 'linewidth', 1.5)
    end
end

% xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\alpha')
colorbar
% axis square
set(gca, 'fontsize', 14, 'ytick', 2:6:14, 'xtick', [])

nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, se_alpha_atmin, se_lim)
% xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
title('\alpha SE')
colorbar
% axis square
set(gca, 'fontsize', 14, 'ytick', [], 'xtick', [])




nexttile([2 2])
hold on
for i_ex = 1:length(eta_ind_ex)
    for j_ex = 1:length(eps_ind_ex)
        plot(tau_at_min{eps_ind_ex(j_ex), eta_ind_ex(i_ex)}, ...
            alpha_at_min{eps_ind_ex(j_ex), eta_ind_ex(i_ex)}, ...
            ex_symbols{i_ex}, ...
            'color', eps_colors(j_ex, :), 'linewidth', 1, 'MarkerSize', 8)
%     legent{i_ex} = ['\eta = ' num2str(eta_vals_sf1(eta_ind_ex(i_ex))) ...
%         ', \epsilon = ' num2str(eps_vals_sf1(eps_ind_ex(i_ex))) ]
    end
end
x_vals = xlim;
y_lims = ylim;
gf_min = min(gamma_fit_ave(:),[],  'omitnan');
gf_max = max(gamma_fit_ave(:), [], 'omitnan');
y1_vals = gf_min*x_vals + 1 - gf_min;
y2_vals = gf_max*x_vals + 1 - gf_max;
ph1 = plot(x_vals, y1_vals, 'k');
ph2 = plot(x_vals, y2_vals, 'k-.');
ylim(y_lims)
legend([ph1 ph2], num2str([gf_min; gf_max], 'gamma fit: %1.2f'), ...
    'Location','southeast', 'color', 'none')


xlabel('\tau')
ylabel('\alpha')
% axis square
title({'Exponents across network replicates';' for selected parameters'})
set(gca, 'FontSize', 16, 'color', 'none')
% legend(legent, 'location', 'northwest')

%%%%%% tau means
nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, mean_tau_atmin, at_lim)
% xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\tau')
colorbar
% axis square
set(gca, 'fontsize', 14, 'ytick', 2:6:14, 'xtick', [])

% SE plots


nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, se_tau_atmin, se_lim)
% xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
title('\tau SE')
colorbar
% axis square
set(gca, 'fontsize', 14, 'ytick', [], 'xtick', [])

%%%%%%%%%%%% gamma predicted, averages
nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, mean_gamma_pred_atmin, gam_lim)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('\gamma_{pred}')
colorbar
% axis square
set(gca, 'fontsize', 14, 'ytick', 2:6:14)

nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, se_gamma_pred_atmin, se_lim)
xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
title('\gamma_{pred} SE')
colorbar
% axis square
set(gca, 'fontsize', 14, 'ytick', [])

nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, gamma_fit_ave, gam_lim)
xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
title('\gamma_{fit}')

colorbar
% axis square
set(gca, 'fontsize', 14, 'ytick', [])




nexttile
imagesc(eta_vals_sf1, eps_vals_sf1, gamma_fit_se, se_lim)
xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
title('\gamma_{fit} SE')

colorbar
axis square
set(gca, 'fontsize', 14, 'ytick', [])

% suptitle('Average and SE of Exponents across N = 10 Simulation Replicates')




print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/figEE_suppl_exponents_' datestr(now, 'yyyymmdd')])
% 
% nexttile
% 
% x = cell2mat(alpha_at_min(:));
% y = cell2mat(tau_at_min(:));
% ph = plot(x, y, 'ko');
%% Figure 3 : scaling scaling everywhere but not a bit to read
% med1Hz_bin = squeeze(median(r_i, 1))/time_bin;
med1Hz_bin = i_calc.r_i_1k/time_bin;

% % this needs to match the length of sel_ee
% etaeps_colors = [0 0.5 0; 1 1 0; 0.8 0.2 0];
% etaeps_markers = {'ks', 'ko', 'k^'};

makeMyFigure(40, 18);
% makeMyFigure(32, 16);
%%%%%%%%%%%%%%%%%%%%% gamma_pred minus gamma_fit
sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err), [-0.2 .2])

% imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err)./se_gamma_pred, [0 2])
hold on 
% [~, cr_h] = contour(i_calc.eta_vals, -i_calc.eps_vals, log10(med1Hz_bin), ...
%     dt_tics, 'LineColor',[1 0 0], 'LineWidth',1.5);

% contour(eta_vals_full, -eps_vals_full, log10(squeeze(mean(r_i, 1))/time_bin), ...
%     dt_tics, 'linecolor', [1 1 1]*0.8)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; usa(50).^.4];
colormap(sh, cmap_nan)
colorbar
% hold on
% plot(eta1_vals, eps1_vals, 'ko-', 'linewidth', 1)
% for i_ee = 1:length(sel_ee)
%     plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'fontsize', 14)
axis square
ylim([-15 -1])
title('\gamma_{pred} - \gamma_{fit}')

%%%%%%%%%%%%%%%%%%%% minimum cutoff
sz_lim = [0.5 2];
sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(min_size_vals), sz_lim);
time_cmap_nan = [0.5*ones(1,3); parula(50)];
colormap(sh, time_cmap_nan)
ch_size = colorbar;
ch_size.Ticks = linspace(sz_lim(1), sz_lim(2), 3);
ch_size.TickLabels = num2str(10.^ch_size.Ticks', '%1.0f');
axis square
set(gca, 'fontsize', 14)
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
ylim([-15 -1])

title('Minimum Size Cutoff')

%%%%%%%%%%%%%%%%%%%% experiment length
exp_T = 1e7*time_bin/3600;
index_sim_ava = mult_ave_num_ava(6, 4);
sh = nexttile;
% time_tics = (-1:1:1) + log10(exp_T); % adaptive
time_tics = log10([3 30 300]); % good for time_bin = 0.01
imagesc(eta_list_sf1, eps_list_sf1, log10(exp_T*index_sim_ava./mult_ave_num_ava), time_tics([1 end]));
hold on
% contour(eta_vals_full, -eps_vals_full, log10(med1Hz_bin), ...
%     [0.5 1], 'LineColor',[1 0 0], 'LineWidth',1.5);
% [~, cr_h] = contour(i_calc.eta_vals, -i_calc.eps_vals, log10(med1Hz_bin), ...
%     dt_tics, 'LineColor',[1 0 0], 'LineWidth',1.5);

time_cmap_nan = [0.5*ones(1,3); parula(50)];
colormap(sh, time_cmap_nan)
ch_time = colorbar;
ch_time.Ticks = time_tics;
ch_time.TickLabels = num2str(10.^ch_time.Ticks', '%1.0f h');
axis square
set(gca, 'fontsize', 14)
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
ylim([-15 -1])

title('Est. Req. Simulation Time')

%%%%%%%%%%%%%%%%%%%% firing rates

nexttile
dt_tics = log10([0.5 1]); %-2:1:2;
% imagesc(eta_vals_full, -eps_vals_full, log10(med1Hz_bin), dt_tics([1 end]));
imagesc(i_calc.eta_vals, -i_calc.eps_vals, log10(med1Hz_bin), [-1.5 1.5])

hold on
% [~, cr_h] = contour(i_calc.eta_vals, -i_calc.eps_vals, log10(med1Hz_bin), ...
%     dt_tics, 'LineColor',[1 0 0], 'LineWidth',1.5);
ch_fr1 = colorbar;
ch_fr1.Ticks = 1.5*[-1 0 1];
ch_fr1.TickLabels = num2str(10.^ch_fr1.Ticks', '%1.0d');
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% hold on
% plot(eta1_vals, eps1_vals, 'wo-', 'linewidth', 2)
% for i_ee = 1:length(sel_ee)
%     plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
% plot(eta1_vals(2), eps1_vals(2), etaeps_markers{2}, 'markersize', 10,...
%     'linewidth', 2, 'MarkerFaceColor',etaeps_colors(2, :))
% plot(eta1_vals(3), eps1_vals(3), etaeps_markers{3}, 'markersize', 10,...
%     'linewidth', 2, 'MarkerFaceColor',etaeps_colors(3, :))

title('Median FR (spikes per s)')
axis square
set(gca, 'fontsize', 14)
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
ylim([-14 0])



%%%%%%%%%%%%%%%% need some sort of goodness-of-alpha fit measure




nexttile
imagesc(i_calc.eta_vals, -i_calc.eps_vals, i_calc.I_si_h10, [0 1.5])
% ylim([0 14])
ylim([-14 0])
% hold on 
% contour(eta_vals_full, eps_vals_full, log10(squeeze(median(r_i, 1))/time_bin), ...
%     -3:0.5:2)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('I(s_{1,.. 10}, h)')
colorbar
hold on
% contour(eta_vals_full, eps_vals_full, log10(med1Hz_bin), ...
%     [0.5 1], 'LineColor',[1 0 0], 'LineWidth',1.5);
[~, cr_h] = contour(i_calc.eta_vals, -i_calc.eps_vals, log10(med1Hz_bin), ...
    dt_tics, 'LineColor',[1 0 0], 'LineWidth',1.5);

% plot(eta1_vals, -eps1_vals, 'wo-', 'linewidth', 2)
% for i_ee = 1:length(sel_ee)
%     plot(eta1_vals(sel_ee(i_ee)), -eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
axis square
set(gca, 'fontsize', 14, 'YDir', 'normal', 'ytick', -14:4:-2)

nexttile
[et1,ep1] = meshgrid(i_calc.eta_vals, i_calc.eps_vals);
[et1q,ep1q] = meshgrid(eta_vals_full, eps_vals_full);

i_per_spk = interp2(et1, ep1, i_calc.mask_info_vals, ...
    et1q, ep1q)./squeeze(mean(r_i, 1)/time_bin);
info_raw = interp2(et1, ep1, i_calc.mask_info_vals, ...
    et1q, ep1q);
info_10 = interp2(et1, ep1, i_calc.I_si_h10, ...
    et1q, ep1q);

imagesc(i_calc.eta_vals, -i_calc.eps_vals, i_calc.info_rate, [0 4]), 
ylim([-14 0])
colorbar
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('I_N / mean pop. FR')
hold on
plot(eta1_vals, -eps1_vals, 'wo-', 'linewidth', 2)
for i_ee = 1:length(sel_ee)
    plot(eta1_vals(sel_ee(i_ee)), -eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end
axis square
set(gca, 'fontsize', 14, 'YDir', 'normal', 'ytick', -14:4:-2)


print(gcf, '-dsvg', ['avalanches_matlab_code/plots/paper_figures/fig3_time_scaling_plots_' num2str(yyyymmdd(datetime))])
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/fig3_time_scaling_plots_' num2str(yyyymmdd(datetime)) '.pdf'])

%% Figure 3: supplemental (fit quality metrics
makeMyFigure(30, 20);
tiledlayout(3, 4, 'TileSpacing','tight')
ks_lim = [-3 -1];
sz_lim = [log10(1) log10(100)];
%%%%%%%%%%%%%%%%%%%% tau fit quality
sh1 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(ks_tau_vals), ks_lim)

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh1, cmap_nan)
ch1 = colorbar;
ch1.Ticks = ks_lim(1):1:ks_lim(2);    
ch1.TickLabels = num2str(10.^ch1.Ticks', '%1.3f');
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])
title('KS statistic (P(S) ~ S^{\alpha})')

%%%%%%%%%%%%%%%%%%%% surrogate tau fit quality 
sh1 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(surr_ks_tau), ks_lim)

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh1, cmap_nan)
ch1 = colorbar;
ch1.Ticks = ks_lim(1):1:ks_lim(2);       
ch1.TickLabels = num2str(10.^ch1.Ticks', '%1.3f');

set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])
title('Surrogate KS, \tau (95%ile)')

%%%%%%%%%%%%%%%%%%%% minimum size values 
sh1 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(min_size_vals), sz_lim)

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh1, cmap_nan)
ch1 = colorbar;
ch1.Ticks = sz_lim; %sz_lim(1):1:sz_lim(2);       
ch1.TickLabels = num2str(10.^ch1.Ticks', '%1.0f');

set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])
title('S_{min} (\tau fit)')

%%%%%%%%%%%%%%%%%%%% alpha fit quality
sh2 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(ks_alpha_vals), ks_lim)

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh2, cmap_nan)
ch2 = colorbar;
ch2.Ticks = ks_lim(1):1:ks_lim(2);       
ch2.TickLabels = num2str(10.^ch2.Ticks', '%1.3f');

set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])
title('KS statistic (P(D) ~ D^{\alpha})')

%%%%%%%%%%%%%%%%%%%%% alpha fit surrogate KS stats
sh3 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(surr_ks_alpha), ks_lim);

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh3, cmap_nan)
ch1 = colorbar;
ch1.Ticks = ks_lim(1):1:ks_lim(2);    
ch1.TickLabels = num2str(10.^ch1.Ticks', '%1.3f');

set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])
title('Surrogate KS, \alpha (95%ile)')

%%%%%%%%%%%%%%%%%%%% minimum duration values 
sh1 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, log10(min_dur_vals), sz_lim)

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh1, cmap_nan)
ch1 = colorbar;
ch1.Ticks = sz_lim(1):1:sz_lim(2);       
ch1.TickLabels = num2str(10.^ch1.Ticks', '%1.0f');

set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])
title('D_{min} (\alpha fit)')
%%%%%%%%%%%%%%%%%%%%% gamma_fit SE

gfit_all_se = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma_se, -1), ... 
    gamma_fit_summary_SF', 'UniformOutput',false));

ave_gfit_se = squeeze(mean(gfit_all_se, 1));

sh3 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, ave_gfit_se)

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('SE \gamma_{fit}')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh3, cmap_nan)
colorbar

set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])

%%%%%%%%%%%%%%%%%%%%% gamma_fit range of scaling

gfit_range = cell2mat(cellfun(@(x) shiftdim(x.fit_gamma_range, -1), ... 
    gamma_fit_summary_SF', 'UniformOutput',false));

ave_range_gfit = squeeze(mean(gfit_range, 1));

sh3 = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, ave_range_gfit, [1 2.5])

axis square
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
title('Decades of scaling (S ~ D^{\gamma})')
% cmap_nan = [1 1 1 ; parula];
% colormap(sh, cmap_nan)
cmap_nan = [0.5*ones(1,3) ; parula(50)];
colormap(sh3, cmap_nan)
colorbar

set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
ylim([-15 -1])

%%%%%%%%%%%%%%%%%%% std dev across runs for gamma_pred 
sh = nexttile;
imagesc(eta_list_sf1, eps_list_sf1, se_gamma_pred, [0 .2])

% imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err)./se_gamma_pred, [0 2])
hold on 
contour(eta_vals_full, -eps_vals_full, log10(squeeze(median(r_i, 1))/time_bin), ...
    -3:0.5:2, 'linecolor', [1 1 1]*0.8)
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
cmap_nan = [1 1 1 ; parula];
colormap(sh, cmap_nan)
colorbar
% hold on
% plot(eta1_vals, eps1_vals, 'wo-', 'linewidth', 2)
% for i_ee = 1:length(sel_ee)
%     plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'fontsize', 14)
axis square
ylim([-15 -1])
title('SE \gamma_{pred}')


% %%%%%%%%%%%%%%%%%%%%% gamma_fit minus gamma_pred
% sh4 = nexttile;
% imagesc(eta_list_sf1, eps_list_sf1, abs(gamma_pred_err), [-0.2 .2])
% 
% axis square
% xlabel('\eta (gain)')
% ylabel('\epsilon (bias)')
% % cmap_nan = [1 1 1 ; parula];
% % colormap(sh, cmap_nan)
% cmap_nan = [0.5*ones(1,3) ; usa(50).^.4];
% colormap(sh4, cmap_nan)
% colorbar
% 
% set(gca, 'ydir', 'normal', 'ytick', -14:4:-2, 'Fontsize', 14)
% ylim([-15 -1])
% title('Crackling: \gamma_{pred} - \gamma_{fit}')
% 

%%%%%%%%%%%%%%%%%% contours of the firing rate plot (median fr = 1 Hz, dt =
%%%%%%%%%%%%%%%%%% xxx)

med_r_i_bin = squeeze(median(r_i, 1));
nexttile
dt_tics = -3.3:1:-0.3;
imagesc(eta_vals_full, -eps_vals_full, log10(med_r_i_bin), dt_tics([1 end]));
hold on
[~, cr_h] = contour(eta_vals_full, -eps_vals_full, log10(med_r_i_bin), ...
    dt_tics, 'LineColor',[1 0 0], 'LineWidth',1.5);
ch_fr1 = colorbar;
ch_fr1.Ticks = dt_tics;
ch_fr1.TickLabels = num2str(10.^dt_tics', '%1.0d');
xlabel('\eta (gain)')
ylabel('\epsilon (bias)')
% hold on
% plot(eta1_vals, eps1_vals, 'wo-', 'linewidth', 2)
% for i_ee = 1:length(sel_ee)
%     plot(eta1_vals(sel_ee(i_ee)), eps1_vals(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
% plot(eta1_vals(2), eps1_vals(2), etaeps_markers{2}, 'markersize', 10,...
%     'linewidth', 2, 'MarkerFaceColor',etaeps_colors(2, :))
% plot(eta1_vals(3), eps1_vals(3), etaeps_markers{3}, 'markersize', 10,...
%     'linewidth', 2, 'MarkerFaceColor',etaeps_colors(3, :))

title('Median FR (spikes per bin)')
axis square
set(gca, 'fontsize', 14)
set(gca, 'ydir', 'normal', 'ytick', -14:4:-2)
ylim([-14 0])
nexttile
ph = plot(mean_fr_ee, mult_ave_num_ava, 'o-', 'linewidth', 2);
xlabel('FR (per bin)')
ylabel('N_{ava}')
assignColorsToLines(ph, parula(length(ph)));
axis square
set(gca, 'fontsize', 14, 'color', 'none')
% legend(ph, num2str(eta_list_sf1', 'eta = %1.0f'), ...
%     'location', 'eastoutside')
ch = colorbar;
ch.Ticks = [0 1];
ch.TickLabels = {'\eta = 1', '\eta = 10'};

nexttile
scaled_sum_fr_ee = log(mean_fr_ee*diag(1./eta_list_sf1));
scaled_num_ava = (mult_ave_num_ava*diag(eta_list_sf1));
ph = plot(scaled_sum_fr_ee, scaled_num_ava, 'o', 'linewidth', 2);
xlabel('log( FR / \eta )')
ylabel('\eta N_{ava}')
assignColorsToLines(ph, parula(length(ph)));
% legend(ph, num2str(eta_list_sf1', 'eta = %1.0f'))
axis square
set(gca, 'fontsize', 14, 'color', 'none')




print(gcf, '-dsvg', ['avalanches_matlab_code/plots/paper_figures/fig3_suppl_plots_' num2str(yyyymmdd(datetime))])
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/fig3_suppl_plots_' num2str(yyyymmdd(datetime))])
%% Trajectory plots

figure

nexttile
hold on
for i_ee = 1:length(sel_ee)
    plot(log10(r_i_ee(r_ord, sel_ee(i_ee))/time_bin), linspace(0, 1, n_neur), etaeps_markers{i_ee}, ...
        'MarkerFaceColor',etaeps_colors(i_ee, :), 'color', 'none')
end

% plot(log10(r_i2(r_ord)/time_bin), linspace(0, 1, length(r_i2)), etaeps_markers{2}, ...
%     'MarkerFaceColor',etaeps_colors(2, :), 'color', 'none')
% plot(log10(r_i3(r_ord)/time_bin), linspace(0, 1, length(r_i3)), etaeps_markers{3}, ...
%     'MarkerFaceColor',etaeps_colors(3, :), 'color', 'none')
axis tight
x_lims = ceil(xlim);
xticks = x_lims(1):0.5:x_lims(2);
set(gca, 'fontsize', 14, 'color', 'none', ...
    'xtick', xticks, 'xticklabel', num2str(10.^xticks', '%1.0g'))
axis square
xlabel('firing rate')
ylabel('cdf over cells')


sh = nexttile;
plot(t_h, median(r_i_ee)/time_bin, 'w', 'linewidth', 1.5)
y_lims = ylim;
y_lims(1) = 0;
hold on
imagesc(t_h, y_lims, [1 ; 1]*norm_gpe_int2_traj)
ylim(y_lims)
xlim(t_h([1 end]))
colormap(sh, parula(2));
ch = colorbar;
ch.Ticks = [0 1];
ch.TickLabels = {'no scaling', 'scaling'};
plot(t_h, median(r_i_ee)/time_bin,'r', 'linewidth', 1.5)

hold on
for i_ee = 1:length(sel_ee)
    plot(t_h(sel_ee(i_ee)), median(r_i_ee(:, sel_ee(i_ee)))/time_bin, etaeps_markers{i_ee}, 'markersize', 10,...
        'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
end 
xlabel('time (hours?)')
ylabel('median firing rate, population')
axis square
yyaxis right
plot(t_h, info_path, 'linewidth', 1.5)
plot(t_h, i_rate_path, 'linewidth', 1.5)
ylabel({'Information rate';'(solid, per bin; dashed, per spike*)'})
% plot(t_h, gpe_int2_traj, 'linewidth', 1.5)
% hold on
% for i_ee = 1:length(sel_ee)
%     plot(t_h(sel_ee(i_ee)), gpe_int2_traj(sel_ee(i_ee)), etaeps_markers{i_ee}, 'markersize', 10,...
%         'linewidth', 2, 'MarkerFaceColor',etaeps_colors(i_ee, :))
% end
% ylabel('\gamma prediction error')
axis square
set(gca, 'fontsize', 14)
%%