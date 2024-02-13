function sum_stats = pullSummaryFitInfo(x_reps, eta_list, eps_list)


num_eta = length(eta_list);
num_eps = length(eps_list);
num_reps = length(x_reps);


tau_values = nan(num_eps, num_eta, num_reps);
tau_x_min = nan(num_eps, num_eta, num_reps);
tau_KS_min = nan(num_eps, num_eta, num_reps);
tau_KS_p95 = nan(num_eps, num_eta, num_reps);
tau_stdSurrKS = nan(num_eps, num_eta, num_reps);
tau_maxSurrKS = nan(num_eps, num_eta, num_reps);
tau_medSurrKS = nan(num_eps, num_eta, num_reps);
se_tau_values = nan(num_eps, num_eta, num_reps);

alpha_values = nan(num_eps, num_eta, num_reps);
alpha_x_min = nan(num_eps, num_eta, num_reps);
alpha_KS_min = nan(num_eps, num_eta, num_reps);
alpha_KS_p95 = nan(num_eps, num_eta, num_reps);
alpha_stdSurrKS = nan(num_eps, num_eta, num_reps);
alpha_maxSurrKS = nan(num_eps, num_eta, num_reps);
alpha_medSurrKS = nan(num_eps, num_eta, num_reps);
se_alpha_values = nan(num_eps, num_eta, num_reps);

gamma_fit_values = nan(num_eps, num_eta, num_reps);
se_gamma_fit_values = nan(num_eps, num_eta, num_reps);

gamma_pred_values = nan(num_eps, num_eta, num_reps);
se_gamma_pred_values = nan(num_eps, num_eta, num_reps);

get_g_val = @(g_arr) cellfun(@(x) x.mhat, g_arr);



for ii = 1:num_reps
    
    % get the indices of eta, epsilon in eta_list, eps_list: ensure that
    % across replications, the eta and epsilon values in the queried list
    % are in the expected locations based on eta_list and eps_list
    rep_eta_vals = x_reps{ii}.x_var_list;
    rep_eps_vals = x_reps{ii}.y_var_list;
    [~, eta_inds, rep_eta_inds] = intersect(eta_list, rep_eta_vals);
    [~, eps_inds, rep_eps_inds] = intersect(eps_list, rep_eps_vals);
    
    
    alpha_values(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat', rep_eps_inds, rep_eta_inds);
    se_alpha_values(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'se_a_hat', rep_eps_inds, rep_eta_inds);
%     alpha_KS_min(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'min_KS_stat', rep_eps_inds, rep_eta_inds);
    alpha_x_min(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'x_min', rep_eps_inds, rep_eta_inds);
    alpha_KS_min(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'min_KS_stat', @(x) log(x), rep_eps_inds, rep_eta_inds);
    alpha_medSurrKS(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'ks_surrogate', @(x) median(log(x)), rep_eps_inds, rep_eta_inds);
    alpha_KS_p95(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'ks_surrogate', @(x) prctile(log(x), 95), rep_eps_inds, rep_eta_inds);
    alpha_maxSurrKS(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'ks_surrogate', @(x) max(log(x)), rep_eps_inds, rep_eta_inds);
    alpha_stdSurrKS(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'ks_surrogate', @(x) std(log(x)), rep_eps_inds, rep_eta_inds);
    

    tau_values(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_tau_pfit', 'a_hat', rep_eps_inds, rep_eta_inds);
    se_tau_values(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_tau_pfit', 'se_a_hat', rep_eps_inds, rep_eta_inds);
%     tau_KS_min(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_tau_pfit', 'min_KS_stat', rep_eps_inds, rep_eta_inds);
    tau_x_min(eps_inds, eta_inds, ii) = pullFieldskipEmpty(x_reps{ii}, 'all_tau_pfit', 'x_min', rep_eps_inds, rep_eta_inds);
    tau_KS_min(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_tau_pfit', 'min_KS_stat', @(x) log(x), rep_eps_inds, rep_eta_inds);
    tau_medSurrKS(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_tau_pfit', 'ks_surrogate', @(x) median(log(x)), rep_eps_inds, rep_eta_inds);
    tau_KS_p95(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_tau_pfit', 'ks_surrogate', @(x) prctile(log(x), 95), rep_eps_inds, rep_eta_inds);
    tau_maxSurrKS(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_tau_pfit', 'ks_surrogate', @(x) max(log(x)), rep_eps_inds, rep_eta_inds);
    tau_stdSurrKS(eps_inds, eta_inds, ii) =  applyCellfunskipEmpty(x_reps{ii}, 'all_tau_pfit', 'ks_surrogate', @(x) std(log(x)), rep_eps_inds, rep_eta_inds);

    gamma_pred_values(eps_inds, eta_inds, ii) = (alpha_values(eps_inds, eta_inds, ii) - 1)./(tau_values(eps_inds, eta_inds, ii) - 1);
    rel_err = sqrt((se_alpha_values(eps_inds, eta_inds, ii)./alpha_values(eps_inds, eta_inds, ii)).^2 + ...
        (se_tau_values(eps_inds, eta_inds, ii)./tau_values(eps_inds, eta_inds, ii)).^2);
    se_gamma_pred_values(eps_inds, eta_inds, ii) = rel_err.*gamma_pred_values(eps_inds, eta_inds, ii);
    %     alpha_vals = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat');
    %     alpha_vals = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat');
    %     alpha_vals = pullFieldskipEmpty(x_reps{ii}, 'all_alpha_pfit', 'a_hat');
    use_gamma = x_reps{ii}.all_gamma_pfit(rep_eps_inds, rep_eta_inds);
    has_entry = cellfun(@(x) ~isempty(x), use_gamma);
    
    gamma_arr = cell(length(rep_eps_inds), length(rep_eta_inds));
    gamma_arr(has_entry) = cellfun(@(x) get_g_val(x), use_gamma(has_entry), ...
        'UniformOutput', false);
    gamma_vals = nan(length(rep_eps_inds), length(rep_eta_inds));
    se_gamma_vals = nan(length(rep_eps_inds), length(rep_eta_inds));
    gamma_vals(has_entry) = cellfun(@(x) mean(x(2:5)), gamma_arr(has_entry) );
    se_gamma_vals(has_entry) = cellfun(@(x) std(x(2:5)), gamma_arr(has_entry) );
    
    
    gamma_fit_values(eps_inds, eta_inds, ii) = gamma_vals;
    se_gamma_fit_values(eps_inds, eta_inds, ii) = se_gamma_vals;
end
clear gamma_vals
clear se_gamma_vals
clear ii
save('dummy_file.mat')

sum_stats = load('dummy_file.mat');

end