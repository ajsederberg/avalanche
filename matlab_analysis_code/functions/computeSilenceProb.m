function p_silence = computeSilenceProb(eta_val, eps_val, J_i, h_vals)

%%
    prob_s_fun = @(h, J1) 1./(1 + exp(eta_val*J1*h - eps_val));
    
    % vectorize
    J_i = J_i(:);
    h_vals = h_vals(:)';

    p_active = prob_s_fun(h_vals, J_i);

    p_silence = exp(sum(log(1 - p_active), 1));
    

    
end