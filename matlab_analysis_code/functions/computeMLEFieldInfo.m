function info_calc = computeMLEFieldInfo(eps_vals, eta_vals, J_i, n_obs)
% Computes the infromation between the MLE estiamte of the field value h
% and the true value h for conditionally independent spins coupled to h
% through the hamiltonian H = eta*J*s + eps*s. 
% % Inputs:
%       eps_vals: vector of epsilon values
%       eta_vals: vector of eta values
%       J_i: vector of coupling constants
%       n_obs: number of observations. Can be 1 or more. Integer. 
% % Output:
%       info_calc: struct with the following fields
%           info_asymp_nats: MLE infomration in nats, using numerical
%           integration over h values with resolution 0.01
%           info_nats_halfres: same, but dh = 0.02. 
% Note that negative information values occur, for instance, when the MLE
% algorithm fails. 


% fixed parameters
h_true_vals = linspace(-5, 5, 1001);
dh = h_true_vals(2) - h_true_vals(1);

h_sig = 1;
prob_h = exp(-0.5*h_true_vals.^2/h_sig^2)/sqrt(2*pi*h_sig^2);

info_asymp_nats = zeros(length(eps_vals), length(eta_vals));
info_asymp2_nats = zeros(length(eps_vals), length(eta_vals));



h_MLE_sigma_arr = cell(length(eps_vals), length(eta_vals));
for i_eps = 1:length(eps_vals)
    eps_val = eps_vals(i_eps); %6;
    for i_phi = 1:length(eta_vals)
        eta_val = eta_vals(i_phi); %7;

%         h_MLE_pars = compute_hMLEpars_forsim(h_true_vals, prob_s_fun, J_i, eps_val, eta_val, n_obs);

        %% compute MI, assuming h ~ N(0, 1): < 0.5 * log(1/(sigma_h^2))>_{h}
        sigma2_hMLE = compute_h_unc(h_true_vals, J_i, eta_val, eps_val, n_obs);

%%
        % Worry about behavior of sigma2_hMLE for values of h near zero
        % when the eps_val is large : if you have effectively no activity
        % because eps_val is so large, then MLE approx is no longer valid
        % and exact calculation must be used. 
        info_hhbar = dh*prob_h*log(1./sigma2_hMLE');
        info_asymp_nats(i_eps, i_phi) = info_hhbar;
        info_hhbar2 = 2*dh*prob_h(1:2:end)*log(1./sigma2_hMLE(1:2:end)');
        info_asymp2_nats(i_eps, i_phi) = info_hhbar2;
        
                
        
        h_MLE_sigma_arr{i_eps, i_phi} = sigma2_hMLE;
        
    end
end

info_calc.info_nats = info_asymp_nats;
info_calc.info_nats_halfres = info_asymp2_nats;
end