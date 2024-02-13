function [h_Fig, h_sp, fit_info] = plotSizeDurPDFScaling(x_reps_str, eta_val, eps_val, tau_val, rep_string, h_Fig, h_sp, data_color, make_subplot)
% Function for plotting size and duration PDF subplots, along with size vs.
% duration scaling plot
% Input strut x_reps_str is loaded from the results directory; provide the
% desired parameter values through eta_val, eps_val, tau_val, and
% rep_string. x_reps_str should be a struct, not a cell array of structs. 
% Can selectively plot subplots by setting make_subplot to false for the plots you don't want. 
% determine the ii_y index of this simulation 

if nargin < 9
    make_subplot = true(length(h_sp));
end

if strcmp(x_reps_str.y_lab_sum, 'replicate')
    ii_y = find(cellfun(@(x) strcmp(x, rep_string), x_reps_str.y_var_list));
    y_val = rep_string;
else
    
    ii_y = find(x_reps_str.y_var_list == eps_val);
    y_val = eps_val;
end

% determine the ii_x index of this simulation 
if strcmp(x_reps_str.x_lab_sum, 'time constant')
    ii_x = find(x_reps_str.x_var_list == tau_val);
    x_val = tau_val;
else
    ii_x = find(x_reps_str.x_var_list == eta_val);
    x_val = eta_val;
end

x_lab_sum = x_reps_str.x_lab_sum;
y_lab_sum = x_reps_str.y_lab_sum;

ava_name_fn = x_reps_str.ava_name_fn;

% set the pdf y limits
pdf_y_lim = [-8 0];
% set size x_lim
size_x_lim = [0 5];
% set duration x_lim
dur_x_lim = [0 4];
d_msz = 10;
%%
if nargin < 6
    h_Fig = makeMyFigure(30, 25);
    h_sp(1) = subplot(3, 3, 1);
    h_sp(2) = subplot(3, 3, 2);
    h_sp(3) = subplot(3, 3, 3);
    data_color = lines(1);
end

try
    x_ava = load(ava_name_fn(ii_y, ii_x));
%     x_fr = load(fr_name_fn(ii_y, ii_x));
    
    [n_cts, size_bins, dur_bins] = calculateJointHistogram(x_ava.sizes, x_ava.durations);
    [ave_sizes, se_sizes, dur_bins, num_obs] = calculateAveSizeByDuration(x_ava.sizes, x_ava.durations);
           
    % require at least one avalanche of size 500
    if max(x_ava.sizes) < 5e2
        return
    end
    %%
    avalanche_count = length(x_ava.sizes);
    pfit_S = x_reps_str.all_tau_pfit{ii_y, ii_x};
    pfit_D = x_reps_str.all_alpha_pfit{ii_y, ii_x};
    % load the selected x_min for Duration and Size
    % distributions
    x_min_D = pfit_D.x_min;
    x_min_S = pfit_S.x_min;
    % load surrogate stats
    ks_stats_Surr_S = pfit_S.ks_surrogate;
    surr_a_S = pfit_S.a_est_surr;
    ks_stats_Surr_D = pfit_D.ks_surrogate;
    surr_a_D = pfit_D.a_est_surr;
    %% extract predicted gamma and fit gamma

    g_pred_all = x_reps_str.all_gamma_pred{ii_y, ii_x}.mhat;
    g_predSE = x_reps_str.all_gamma_pred{ii_y, ii_x}.mSE;
    pfit_G = x_reps_str.all_gamma_pfit{ii_y, ii_x};
    
    %% adaptive fit for gamma
    plot_diagnostics = true;
    [gamma_fit, gamma_range] = fitGammaRange(dur_bins, ave_sizes, plot_diagnostics);
    %% 
    x_min_vals = x_reps_str.x_min_vals;
    par_title_str = x_reps_str.par_title_str;
    
catch
    return
end

% %% reset current figure
% figure(h_Fig);

%% get slopes and error estimates
% title({['\tau = ' num2str(pfit_S.a_hat, '%1.3f') ' ' num2str(pfit_S.se_a_hat, ' (%1.3f SE)')]; ...
%     ['s_{min} = ' num2str(pfit_S.x_min) ', n_{tail} = ' num2str(n_tail_S)]})

fit_info.size_tau_best = pfit_S.a_hat;
% error in tau : estimate from changing the cutoff by one position
[~, ind] = find(pfit_S.x_min_vals == pfit_S.x_min);
nb_ind = ind + [-1 0 1];
nb_ind = nb_ind(nb_ind > 0 & nb_ind <= length(pfit_S.x_min_vals));
err_cutoff = range(pfit_S.a_vals(nb_ind))/(length(nb_ind)-1);
fit_info.size_tau_err = sqrt(err_cutoff^2 + pfit_S.se_a_hat^2);

% error in tau : how much does the slope estimate vary across cut-offs with
% ks-stats within X% of the minimum (where X% is relative error from the
% surrogate distribution
rel_err_frac = std(pfit_S.ks_surrogate)/mean(pfit_S.ks_surrogate);
[~, rerfr_inds] = find(pfit_S.ks_stats <= (1 + rel_err_frac)*min(pfit_S.ks_stats));
rerfr_err = std(pfit_S.a_vals(rerfr_inds));
fit_info.small_KSstat_tau_se = rerfr_err;

fit_info.dur_alpha_best = pfit_D.a_hat;
% error in tau : estimate from changing the cutoff by one position
[~, ind] = find(pfit_D.x_min_vals == pfit_D.x_min);
nb_ind = ind + [-1 0 1];
nb_ind = nb_ind(nb_ind > 0 & nb_ind <= length(pfit_D.x_min_vals));
err_cutoffD = range(pfit_D.a_vals(nb_ind))/(length(nb_ind)-1);
fit_info.dur_alpha_err = sqrt(err_cutoffD^2 + pfit_D.se_a_hat^2);

fit_info.gamma_pred = (fit_info.dur_alpha_best - 1)/(fit_info.size_tau_best - 1);

% propagate the error to get uncertainty in gamma predicted
[gamma_err, gamma_0] = propagateErrorDxi(@(x) (x(1)-1)/(x(2)-1), ...
    [fit_info.dur_alpha_best fit_info.size_tau_best], ...
    [fit_info.dur_alpha_err fit_info.size_tau_err]);
fit_info.gamma_pred_err = gamma_err;
fit_info.gamma_pred0 = gamma_0;

% find the uncertainty in gamma fit 
[gamma_fit_mean, gamma_fit_se, gamma_range1] = extractGammaFitRegion(pfit_G);

%% Subplot 1
% if h_sp is not empty, then make plot
if make_subplot(1)
set(gcf, 'CurrentAxes', h_sp(1));
freq_cts_size = sum(n_cts, 2);
log_size_bins = unique(round(logspace(0, 5, 51)));  % use unique(round()) to avoid having 1, 1.25, 1.58, etc.
size_pdf = histcounts(x_ava.sizes, bins2edges(log_size_bins), 'Normalization', 'pdf');
size_cts = histcounts(x_ava.sizes, bins2edges(log_size_bins), 'Normalization', 'count');
rel_err = 1./sqrt(size_cts-1);

size_pdf_SE = rel_err.*size_pdf; % get counting statistics error bars on PDF estimates

% ERRORBAR PLOT
errorbar(log10(log_size_bins), log10(size_pdf), ...
    log10(size_pdf) - log10(size_pdf-size_pdf_SE), ...  % lower errorbar
    log10(size_pdf + size_pdf_SE) - log10(size_pdf), '.', ...  % upper errorbar
    'markersize', d_msz, 'capsize', 0, 'color', data_color)

hold on
s_cutoff = [log10(pfit_S.x_min) max(log10(size_bins))];
use_bins = log10(log_size_bins) >= (s_cutoff(1)) & isfinite(log10(size_pdf));
b0 = mean((pfit_S.a_hat*log10(log_size_bins(use_bins)) + log10(size_pdf(use_bins))));

% FIT PLOT
phFit = plot(s_cutoff,b0 - pfit_S.a_hat*s_cutoff, 'k-', 'linewidth', 1);
plot([0 s_cutoff(1)],b0 - pfit_S.a_hat*[0 s_cutoff(1)], ':', ...
    'Color', phFit.Color, 'linewidth', 2);

% CUTOFF LINE
plot(log10(pfit_S.x_min)*[ 1 1], pdf_y_lim, 'color', 0.5*[1 1 1])
% xlim([0 dur_bins(end-1)])
n_tail_S = sum(freq_cts_size(size_bins >= pfit_S.x_min));
%                 ylim([-1 b0])
ylabel('size pdf (log_{10})')
xlabel('avalanche size')
title(['\tau = ' num2str(fit_info.size_tau_best, '%1.2f') ' ' ...
    num2str(fit_info.size_tau_err, ' (%1.2f SE)')])

set(gca, 'color', 'none', 'xtick', 0:1:max(log10(size_bins)), ... 
    'xticklabel', num2str((0:1:max(log10(size_bins)))', '10^{%1.0f}'))
axis square
xlim(size_x_lim)
ylim(pdf_y_lim)
end
%% subplot 2
if make_subplot(2)
set(gcf, 'CurrentAxes', h_sp(2));

log_dur_bins = unique(round(logspace(0, 5, 51)));  % use unique(round()) to avoid having 1, 1.25, 1.58, etc.
dur_pdf = histcounts(x_ava.durations, bins2edges(log_dur_bins), 'Normalization', 'pdf');
dur_cts = histcounts(x_ava.durations, bins2edges(log_dur_bins), 'Normalization', 'count');
rel_err_dur = 1./sqrt(dur_cts-1);

dur_pdf_SE = rel_err.*dur_pdf; % get counting statistics error bars on PDF estimates

% ERRORBAR PLOT
errorbar(log10(log_dur_bins), log10(dur_pdf), ...
    log10(dur_pdf) - log10(dur_pdf-dur_pdf_SE), ...  % lower errorbar
    log10(dur_pdf + dur_pdf_SE) - log10(dur_pdf), '.', ...  % upper errorbar
    'markersize', d_msz, 'capsize', 0, 'color', data_color)
hold on
d_cutoff = [log10(pfit_D.x_min) max(log10(dur_bins))];
use_bins = log10(log_dur_bins) >= (d_cutoff(1)) & isfinite(log10(dur_pdf));
b0 = mean((pfit_D.a_hat*log10(log_dur_bins(use_bins)) + log10(dur_pdf(use_bins))));
%

% FIT LINE PLOT
phFit2 = plot(d_cutoff, b0 - pfit_D.a_hat*d_cutoff, 'k-', 'linewidth', 1);
plot([0 d_cutoff(1)], b0 - pfit_D.a_hat*[0 d_cutoff(1)], ':', ...
    'Color', phFit2.Color, 'linewidth', 2);


% LOWER CUTOFF LINE
plot(log10(pfit_D.x_min)*[ 1 1], pdf_y_lim, 'color', 0.5*[1 1 1])

% TITLE AND LABEL
title(['\alpha = ' num2str(fit_info.dur_alpha_best, '%1.2f') ' ' ...
    num2str(fit_info.dur_alpha_err, ' (%1.2f SE)')])

ylabel('duration pdf (log_{10})')
xlabel('avalanche duration')
set(gca, 'color', 'none', 'xtick', 0:1:max(log10(dur_bins)), ... 
    'xticklabel', num2str((0:1:max(log10(dur_bins)))', '10^{%1.0f}'))
axis square
xlim(dur_x_lim)
ylim(pdf_y_lim)
end
%% subplot 3
if make_subplot(3)
% get the predicted gamma at the min-KS value for size and duration
g_pred = g_pred_all(x_min_vals == x_min_S,x_min_vals == x_min_D);
set(gcf, 'CurrentAxes', h_sp(3));
hold on
% ERRORBAR PLOT (DATA)
errorbar(log10(dur_bins), log10(ave_sizes),log10(ave_sizes) - log10(ave_sizes - se_sizes), ...
    log10(ave_sizes + se_sizes)- log10(ave_sizes), '.', ...
    'markersize', d_msz, 'capsize', 0, 'color', data_color)
% plot([0 d_cutoff(2)], pfit_G{use_G_fit}.bhat + pfit_G{use_G_fit}.mhat*[0 d_cutoff(2)], '-', 'linewidth', 1);
% FIT PLOT
plot([0 d_cutoff(2)], gamma_fit.p2 + gamma_fit.p1*[0 d_cutoff(2)], 'k-', 'linewidth', 1.5)
% plot([0 d_cutoff(2)], gamma_fit_mean*[0 d_cutoff(2)], 'k-', 'linewidth', 1.5);
% PREDICTION PLOT
pred_color = [0.9290    0.6940    0.1250];
plot([0 d_cutoff(2)], g_pred*[0 d_cutoff(2)], '-', ...
    'color', pred_color, 'linewidth', 1.5);

plot(gamma_range, 0*gamma_range, '-', 'color', 0.5*[1 1 1 ], 'linewidth', 3)
p1_CI = confint(gamma_fit);
title({['\gamma_{fit} = ' num2str(gamma_fit.p1, '%1.2f') ' ' ...
    num2str(p1_CI(:, 1)', '(%1.2f, %1.2f CI)') ', ' num2str(diff(gamma_range)) ' decades']; ...
    ['\gamma_{pred} = ' num2str(gamma_0, '%1.2f') ' ' num2str(gamma_err, '(%1.2f SE)')]});

xlabel('duration')
ylabel('average size')
set(gca, 'color', 'none', 'xtick', 0:1:max(log10(dur_bins)), ... 
    'xticklabel', num2str((0:1:max(log10(dur_bins)))', '10^{%1.0f}'))
lh = legend({'average size', 'fit', 'predicted', 'range'},'Location', 'northwest', 'Box','off', 'Color','none');

% lh = legend({'average size', 'fit', 'predicted', 'range'},'Location', 'EastOutside');
% lh.Position(1) = lh.Position(1)-0.06;
% lh.Position(2) = lh.Position(2)-0.03;
axis square
xlim(dur_x_lim)
ylim(size_x_lim)
end
%% subplot 4
if make_subplot(4)
set(gcf, 'CurrentAxes', h_sp(4))
% axes('Position',[.7 .2 .2 .2])
hold on
fit_error = log10(ave_sizes) - gamma_fit_mean*log10(dur_bins);
pred_error = log10(ave_sizes) - gamma_0*log10(dur_bins);
plot(log10(dur_bins), pred_error, '.', 'markersize', d_msz, 'color', data_color); %, ...
plot(log10(dur_bins), 0*log10(dur_bins), '-', 'color', pred_color, 'linewidth', 1.5)
% errorbar(log10(dur_bins), fit_error, log10(se_sizes), '.', ...
%     'capsize', 0, 'color', 0.2*[1 1 1])
ylim([-0.5 1])
axis square
xlabel('duration (log_{10})')
ylabel('prediction error')
set(gca, 'color', 'none', 'xtick', 0:1:max(log10(dur_bins)), ... 
    'xticklabel', num2str((0:1:max(log10(dur_bins)))', '10^{%1.0f}'))
end
% try
%     suptitle([par_title_str x_lab_sum ' = ' num2str(x_val) ', ' y_lab_sum y_val])
% catch
%     suptitle([par_title_str x_lab_sum ' = ' num2str(x_val) ', ' y_lab_sum num2str(y_val)])
%     
% end



% subplot(3, 3, 4)
% errorbar(log10(pfit_S.x_min_vals), pfit_S.a_vals, pfit_S.se_a_vals)
% ylabel('\tau_{MLE}')
% xlabel('s_{min}')
% yyaxis right
% plot(log10(pfit_S.x_min_vals), pfit_S.ks_stats)
% frac_tail_S = cellfun(@(x) mean(x_ava.sizes > x), num2cell(x_min_vals));
% %                 ks_statS_raw = pfit_S.ks_stats./frac_tail_S;
% % hold on
% % plot(log10(pfit_S.x_min_vals),ks_statS_raw)
% ylabel('KS statistic')
% set(gca, 'color', 'none')
% 
% subplot(3, 3, 5)
% errorbar(log10(pfit_D.x_min_vals), pfit_D.a_vals, pfit_D.se_a_vals)
% ylabel('\alpha_{MLE}')
% xlabel('d_{min}')
% yyaxis right
% plot(log10(pfit_D.x_min_vals), pfit_D.ks_stats)
% frac_tail_D = cellfun(@(x) mean(x_ava.durations > x), num2cell(x_min_vals));
% %                 ks_statD_raw = pfit_D.ks_stats./frac_tail_D;
% % hold on
% % plot(log10(pfit_D.x_min_vals),ks_statD_raw)
% ylabel('KS statistic')
% set(gca, 'color', 'none')
% 
% subplot(3, 3, 6)
% d_low_cut = cellfun(@(x) x.x_cutoff(1), pfit_G);
% m_hat_vs_x0 = cellfun(@(x) x.mhat, pfit_G);
% se_m_hat_vs_x0 = cellfun(@(x) x.mSE, pfit_G);
% x0 = cellfun(@(x) x.x_cutoff(1), pfit_G);
% hold on
% errorbar(x0, m_hat_vs_x0, se_m_hat_vs_x0, 'color', 0.2*[1 1 1])
% 
% % also plot line for gamma predicted
% %                 ['\gamma_{pred} = ' num2str(g_pred_all(x_min_vals == x_min_D,x_min_vals == x_min_S), '%1.3f')...
% %                     ' ' num2str(g_predSE(x_min_vals == x_min_D,x_min_vals == x_min_S), '(%1.3f SE)')]
% % prior to 5/5/22, these entries were transposed (S adn D
% % switched) which selected the wrong value.
% g_pred_ave = g_pred_all(x_min_vals == x_min_S,x_min_vals == x_min_D);
% g_pred_SE = g_predSE(x_min_vals == x_min_S,x_min_vals == x_min_D);
% plot(x0([1 end]), g_pred_ave*[ 1 1],'k', 'linewidth', 1)
% plot(x0([1 end]), g_pred_ave*[ 1 1] + g_pred_SE, 'k', 'linewidth', 0.5)
% plot(x0([1 end]), g_pred_ave*[ 1 1] - g_pred_SE, 'k', 'linewidth', 0.5)
% xlabel('minimum duration')
% ylabel('\gamma_{fit} for 1 decade')
% set(gca, 'color', 'none')
% 
% subplot(3,3,7)
% cla;
% ks_bins_S = logspace(log10(min(ks_stats_Surr_S)), log10(max(ks_stats_Surr_S)), 21);
% hold on
% histogram(log10(ks_stats_Surr_S), log10(ks_bins_S), 'EdgeColor', 'none')
% y_lims = ylim;
% ks_stat_hat = pfit_S.ks_stats(pfit_S.x_min_vals == pfit_S.x_min);
% plot(log10(ks_stat_hat)*[1 1], y_lims)
% xlabel('KS statistic')
% ylabel(['pdf (' num2str(length(ks_stats_Surr_S)) ' surrogates)'])
% title(['power law for s > ' num2str(pfit_S.x_min) ...
%     ', \tau = ' num2str(pfit_S.a_hat)])
% set(gca, 'color', 'none')
% 
% subplot(3,3,8)
% cla;
% ks_bins_D = logspace(log10(min(ks_stats_Surr_D)), log10(max(ks_stats_Surr_D)), 21);
% 
% hold on
% histogram(log10(ks_stats_Surr_D), log10(ks_bins_D), 'EdgeColor', 'none')
% 
% % histogram(ks_stats_Surr_D,  'EdgeColor', 'none')
% y_lims = ylim;
% ks_stat_hat = pfit_D.ks_stats(pfit_D.x_min_vals == pfit_D.x_min);
% plot(log10(ks_stat_hat)*[1 1], y_lims)
% xlabel('KS statistic')
% ylabel(['pdf (' num2str(length(ks_stats_Surr_D)) ' surrogates)'])
% title(['power law for s > ' num2str(pfit_D.x_min) ...
%     ', \alpha = ' num2str(pfit_D.a_hat)])
% set(gca, 'color', 'none')
% 
% subplot(3, 3, 9)
% g_xy = 1:length(pfit_S.x_min_vals);
% imagesc(g_xy, g_xy, g_pred_all)
% hold on
% sm_g_pred_all = imfilter(g_pred_all, fspecial('gaussian', [5 5], 2), 'replicate');
% 
% ave_fit_g = mean(m_hat_vs_x0(d_low_cut > log10(x_min_D) & d_low_cut < 2));
% se_fit_g = std(m_hat_vs_x0(d_low_cut > log10(x_min_D) & d_low_cut < 2));
% if isfinite(se_fit_g)
%     contour(g_xy, g_xy, sm_g_pred_all, ...
%         ave_fit_g + se_fit_g*[-2 2], 'color', [1 1 1])
% end
% plot(find(x_min_S == x_min_vals), find(x_min_D == x_min_vals), 'r*')
% colorbar
% set(gca, 'xtick', g_xy(1:12:end), 'XTickLabel', x_min_vals(1:12:end))
% set(gca, 'ytick', g_xy(1:8:end), 'YTickLabel', x_min_vals(1:8:end))
% set(gca, 'ydir', 'normal')
% axis square
% xlabel('d_{min}')
% ylabel('s_{min}')
% title('\gamma_{pred}')


end