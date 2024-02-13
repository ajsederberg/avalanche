function [h_MLE_unc, f_info] = compute_h_unc(h_val, J_i, eta_val, eps_val, T)
% Estimates the SE of the MLE for h at this value of epsilon, phi, and J, 
% Also returns the Fisher information (f_info = 1/sigma^2)
% Computed from the second derivative w.r.t. h of the likelihood, which is
% -1/sigma^2, the uncertainty around the max-likelihood estimate of h
% this depends only on eps_val, phi_val, and J_i
% Inputs:
%       h_val: MLE value of h. if passed as a row vector, h_MLE_unc is
%       calculated at each value of h. 
%       J_i: num_neurons by 1 coupling coefficients between latent field
%           and neural observations
%       eps_val: scalar, bias toward silence/firing 
%       phi_val: scalar, multiplies Jhs term
%
% Output: 
%       h_MLE_unc: scalar, sigma. 2nd derivative of p({s_i}|h) is -1/sigma^2 
%       f_info: scalar, if h = hMLE, then this is the Fisher information

    f_info = T*sum((eta_val*J_i*ones(size(h_val))).^2./(4*(cosh(0.5*(eta_val*J_i*h_val + eps_val)).^2)), 1);
    
    h_MLE_unc = 1./sqrt(f_info);

end