function gamma_fit_arr = consolidateGammaFits(x_reps_uf1, eps_ind_trace, eta_ind_trace)

gamma_fit_arr = cell(length(x_reps_uf1)+1, length(eps_ind_trace));
for i_xra = 1:length(x_reps_uf1)+1
    for jj = 1:length(eta_ind_trace)
        if i_xra > length(x_reps_uf1)
            g_fit = x_reps_f1.all_gamma_pfit{eps_ind_trace(jj), eta_ind_trace(jj)};
        else
            g_fit = x_reps_uf1{i_xra}.all_gamma_pfit{eps_ind_uf(jj), eta_ind_uf(jj)};
        end

        try 

        [gamma_fit_mean, gamma_fit_se, gamma_range, full_res] = extractGammaFitRegion(g_fit);
        gamma_fit_arr{i_xra, jj}.gamma_fit_mean = gamma_fit_mean;
        gamma_fit_arr{i_xra, jj}.gamma_fit_se = gamma_fit_se;
        gamma_fit_arr{i_xra, jj}.gamma_range = gamma_range;
        gamma_fit_arr{i_xra, jj}.full_gamma_fit_mean = full_res.mean_gamma_over_range;
        gamma_fit_arr{i_xra, jj}.full_gamma_fit_se = full_res.se_gamma_over_range;
        gamma_fit_arr{i_xra, jj}.full_gamma_end_ranges = full_res.gamma_end_ranges;
        catch
            gamma_fit_arr{i_xra, jj}.gamma_fit_mean = nan*gamma_fit_mean;
            gamma_fit_arr{i_xra, jj}.gamma_fit_se = nan*gamma_fit_se;
            gamma_fit_arr{i_xra, jj}.gamma_range = nan*gamma_range;
            gamma_fit_arr{i_xra, jj}.full_gamma_fit_mean = nan*full_res.mean_gamma_over_range;
            gamma_fit_arr{i_xra, jj}.full_gamma_fit_se = nan*full_res.se_gamma_over_range;
            gamma_fit_arr{i_xra, jj}.full_gamma_end_ranges = nan*full_res.gamma_end_ranges;
        end

    end
end
end