function plotScalingVerification(x_reps_T1, eta_T1_infTC, eps_T1_infTC, T1_tau_vals, T1_taurep_vals, gamma_fit_summary, ii_rep, jj_t, t1_codes)


h_Fig = makeMyFigure(30, 10);
tiledlayout(1, 4, 'TileSpacing','compact');
h_sp = zeros(4, 1);
for ii =1 : 4
    h_sp(ii) = nexttile;
end
data_color = [0.4940    0.1840    0.5560];

if ~isempty(T1_tau_vals)
    % this plots results for selected tau, replicate
    plotSizeDurPDFScaling(x_reps_T1, eta_T1_infTC, eps_T1_infTC, T1_tau_vals(jj_t),...
        T1_taurep_vals{ii_rep}, h_Fig, h_sp, data_color);
else
    % this plots results for selected eta, epsilon
    plotSizeDurPDFScaling(x_reps_T1, eta_T1_infTC(jj_t), eps_T1_infTC(ii_rep), T1_tau_vals,...
        T1_taurep_vals, h_Fig, h_sp, data_color);
end
title_str = "Simulation " + t1_codes(ii_rep, jj_t);
suptitle(title_str)
% set(h_Fig, 'CurrentAxes', h_sp(3))

set(h_Fig, 'CurrentAxes', h_sp(4))
cla;
%% add to this plot the gamma fit values vs start point for fit
gamma_fit_arr = x_reps_T1.all_gamma_pfit{ii_rep, jj_t};
gamma_val = cellfun(@(x) x.mhat, gamma_fit_arr);
gamma_se = cellfun(@(x) x.mSE, gamma_fit_arr);
x_start = cellfun(@(x) x.x_cutoff(1), gamma_fit_arr);
errorbar(x_start, gamma_val, gamma_se, 'color', 0.3*[1 1 1],'LineWidth',1);
ylabel('gamma (fit) +/- SE')
title('gamma fit vals by decade')
hold on
plot(gamma_fit_summary.fit_gamma_range_start(ii_rep, jj_t) + [0 gamma_fit_summary.fit_gamma_range(ii_rep, jj_t)-1], ...
    gamma_fit_summary.fit_gamma(ii_rep, jj_t)*[1 1], 'k-', 'LineWidth', 2)
axis tight
%%
filedir = 'avalanches_matlab_code/plots/scaling_verification/';
filename = ['simulation_' num2str(t1_codes(ii_rep, jj_t), '%06.0f')];
print(gcf, '-dpdf', [filedir filename])
print(gcf, '-dpng', [filedir filename])

%% now clear axis and add back in text 
cla;
set(gca, 'Visible', 'off')
text(0, 1, title_str)

name_str = x_reps_T1.param_str(ii_rep, jj_t);
text(0, 0.75, name_str, 'Interpreter','none')

text(0, 0.5, ['Total Avalanche count: ' num2str(x_reps_T1.avalanche_count(ii_rep, jj_t))])

% tau fit z-scores
tau_KS = x_reps_T1.all_tau_pfit{ii_rep, jj_t}.min_KS_stat;
tau_KS_surr = x_reps_T1.all_tau_pfit{ii_rep, jj_t}.ks_surrogate;
z_tau = abs(tau_KS - mean(tau_KS_surr))/std(tau_KS_surr);
text(0, 0.25, ['KS val (z-score) tau: ' num2str(tau_KS, "%1.4f") num2str(z_tau, ' (%1.1f)')])

% alpha fit z-scores
alpha_KS = x_reps_T1.all_alpha_pfit{ii_rep, jj_t}.min_KS_stat;
alpha_KS_surr = x_reps_T1.all_alpha_pfit{ii_rep, jj_t}.ks_surrogate;
z_alpha = abs(alpha_KS - mean(alpha_KS_surr))/std(alpha_KS_surr);
text(0, 0, ['KS val (z-score) alpha: ' num2str(alpha_KS, "%1.4f") num2str(z_alpha, ' (%1.1f)')])

% gamma fit info
text(0, -0.25, ['gamma range = ' ...
    num2str(gamma_fit_summary.fit_gamma_range(ii_rep, jj_t), '%1.2f') ...
    ' start at ' num2str(gamma_fit_summary.fit_gamma_range_start(ii_rep, jj_t), '%1.2f')]);

print(gcf, '-dpdf', [filedir 'key_' filename '_' name_str '.pdf'])