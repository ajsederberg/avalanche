function h_MLE_pars = compute_hMLEpars_forsim(h_vals, prob_s_fun, J_i, eta_val, eps_val, n_obs)
% simulates the system H = eta_val*J_i*h*s_i + eps_val*s_i, drawing n_obs
% samples for the "true" h values provided in h_vals. From samples,
% computes the MLE parameters (estimate and SE of estimate). 
% prob_s_fun is (1 + exp(eta*Jh + eps))^{-1}. 

    n_neur = length(J_i);   % length of J_i is number of neurons
    h_MLE_pars = zeros(length(h_vals), 3);

    for i_htr = 1:length(h_vals)
        h_true = h_vals(i_htr);
        p_s_i = prob_s_fun(h_true, J_i, eta_val, eps_val);
        if ~isinf(n_obs)
            % average over observations
            s_i_obs = mean(rand(n_neur, n_obs) < repmat(p_s_i, 1, n_obs), 2);
        else
            s_i_obs = p_s_i;
        end
        h_MLE = find_hMLE(s_i_obs, J_i, eps_val, eta_val);
        h_MLE_pars(i_htr, 1) = h_MLE;
        h_MLE_pars(i_htr, 2) = compute_h_unc(h_MLE, J_i, eta_val, eps_val, n_obs);
        % compute sigma_h at h_true as well
        h_MLE_pars(i_htr, 3) = compute_h_unc(h_true, J_i, eta_val, eps_val, n_obs);


    end

end