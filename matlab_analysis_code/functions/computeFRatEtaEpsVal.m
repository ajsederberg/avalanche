function r_i = computeFRatEtaEpsVal(eta_val, eps_val, J_i)

    n_neur = length(J_i);
%     J_i = randn(n_neur, 1);
    prob_s_fun = @(h, J1) 1./(1 + exp(eta_val*J1*h + eps_val));
    
    
    %% Compute, for fixed eta, epsilon, J, h_MLE and SE_h_MLE as h_true varies. 
    % fixed parameters
    % h_true_vals = linspace(-5, 5, 1001);
    % dh = h_true_vals(2) - h_true_vals(1);
    
    h_sig = 1;
    prob_h = @(h) normpdf(h, 0, h_sig);
    % prob_h = exp(-0.5*h_true_vals.^2/h_sig^2)/sqrt(2*pi*h_sig^2);
    %%
    r_i = zeros(n_neur, 1);

    for i_n = 1:n_neur
        int_fun = @(h) prob_s_fun(h, J_i(i_n)).*prob_h(h);
        r_i(i_n) = integral(@(x) int_fun(x), -5, 5);
    end

end