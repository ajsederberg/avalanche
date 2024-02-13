function res_str = consolidateAlphaTauFits(p_fit_arr, eps_val, eta_val, exp_field)
% res_str has fields [alpha_est, se_alpha_est, x_min_est, ks_min] for all
% simulations with specified eps_val and eta_val
% consolidation: looks at all k_min vs x_min, then picks a single x_min.
% Pulls the fitted exponents from each run. Average and SE to get overall
% estimate. 


%%

% first pull out fitting info
if iscell(p_fit_arr)
    exp_fits = cellfun(@(x) x.(['all_' exp_field '_pfit'])...
        {abs(x.y_var_list) == abs(eps_val), abs(x.x_var_list) == abs(eta_val)}, ...
        p_fit_arr, 'UniformOutput', false);
else
    exp_fits = p_fit_arr.(['all_' exp_field '_pfit'])...
        {abs(p_fit_arr.y_var_list) == abs(eps_val), abs(p_fit_arr.x_var_list) == abs(eta_val)};
end

%% get rid of empty entries (missing simulation runs)
empty_row = cellfun(@(x) isempty(x), exp_fits);
exp_fits = exp_fits(~empty_row);

%%

x_vals = log10(cell2mat(cellfun(@(x) x.x_min_vals, exp_fits, 'UniformOutput', false)));
exp_vals = cell2mat(cellfun(@(x) x.a_vals, exp_fits, 'UniformOutput',false));
se_exp_vals = cell2mat(cellfun(@(x) x.se_a_vals, exp_fits, 'UniformOutput',false));
ks_vals = cell2mat(cellfun(@(x) x.ks_stats, exp_fits, 'UniformOutput',false));

% also grab the surrogate values
x_surr = cellfun(@(x) x.x_min, exp_fits);
ks_surr = cell2mat(cellfun(@(x) x.ks_surrogate', exp_fits, 'UniformOutput',false));

% check that all x_vals are the same
same_x = std(x_vals, [], 1) < 1e-4;
if ~all(same_x)
    disp('Check x_vals. They appear not to match across entries in p_fit_arr.')
end

%% choose x_min...

figure()

nexttile
plot(x_vals', ks_vals')
hold on
plot(x_vals(1, :)', mean(ks_vals, 1)', 'k', 'LineWidth',1.5)
xlabel('minimum cutoff')
ylabel('KS statistic')
set(gca, 'XTick', log10([3 10 30 100]), 'XTickLabel', [3 10 30 100])
[ks_min, ind_min] = min(exp(mean(log(ks_vals), 1)));
x_min_est = x_vals(1, ind_min);
plot(x_min_est, ks_min, 'r*', 'LineWidth',1.5)

if any(abs(x_surr - 10^x_min_est) < 2)
    use_surr = abs(x_surr - 10^x_min_est) < 2;

    nexttile
    histogram(log10(ks_surr(use_surr, :)))
    hold on
    plot(log10(ks_min)*[1 1], ylim)
    title(['surrogate set KS vals at x_{min} = ' num2str(10^x_min_est)])
    xlabel('log_{10} (KS)')
    res_str.ks_surr = ks_surr(use_surr, :);
else
    res_str.ks_surr = [];
end

% compute mean value of exponent at selected xmin
exp_val_est = mean(exp_vals(:, ind_min));
se_exp_val_est = std(exp_vals(:, ind_min));

nexttile
errorbar(x_vals', exp_vals', se_exp_vals')
hold on
errorbar(mean(x_vals, 1), mean(exp_vals, 1), std(exp_vals), 'k', 'linewidth', 1)
xlabel('minimum cutoff')
ylabel(exp_field)
title(['average at overall min cutoff: ' num2str(exp_val_est, '%1.2f')...
    '\pm ' num2str(se_exp_val_est, '%1.2f')])

nexttile
plot(x_vals(1, :), std(exp_vals, [], 1))
xlabel('minimum cutoff')
ylabel('SD across runs')

% compute KS stat for each run at this min value
ks_at_min = cellfun(@(x) x.ks_stats(ind_min), exp_fits);

% save output
res_str.ks_at_min = ks_at_min;
res_str.exp_at_min = cellfun(@(x) x.a_vals(ind_min), exp_fits);

res_str.x_min_est = 10^x_min_est;
res_str.exp_val_est = exp_val_est;
res_str.se_exp_val_est = se_exp_val_est;