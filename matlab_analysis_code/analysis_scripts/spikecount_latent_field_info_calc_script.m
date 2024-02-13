% Information calculation for population spike count and maximum likelihood
% estimator across phi, epsilon values. 

n_neur = 128; 

J_i = randn(n_neur/2, 1);

J_i = [-J_i ; J_i];
prob_s_fun = @(h, J1, eta1, eps1) 1./(1 + exp(eta1*J1*h + eps1));

%% Compute, for fixed eta, epsilon, J, h_MLE and SE_h_MLE as h_true varies. 
% fixed parameters
h_true_vals = linspace(-5, 5, 1001);
dh = h_true_vals(2) - h_true_vals(1);

h_sig = 1;
prob_h = exp(-0.5*h_true_vals.^2/h_sig^2)/sqrt(2*pi*h_sig^2);


n_obs = 1; 
eps_vals = 0:0.5:14 ; % eps values, range from SF1 simulations
eta_vals = 0:0.25:10; % eta values, range from SF1 simulations

info_asymp_nats = zeros(length(eps_vals), length(eta_vals));
info_asymp2_nats = zeros(length(eps_vals), length(eta_vals));

h_MLE_arr = cell(length(eps_vals), length(eta_vals));
err_MLE  = zeros(length(eps_vals), length(eta_vals));
for i_eps = 1:length(eps_vals)
    eps_val = eps_vals(i_eps); %6;
    parfor i_eta = 1:length(eta_vals)
        eta_val = eta_vals(i_eta); %7;
%%
        h_MLE_pars = compute_hMLEpars_forsim(h_true_vals, prob_s_fun, J_i, eps_val, eta_val, n_obs);

        %% compute MI, assuming h ~ N(0, 1): < 0.5 * log(1/(sigma_h^2))>_{h}
        sigma2_hMLE = compute_h_unc(h_true_vals, J_i, eta_val, eps_val, n_obs);

        %% This approximation sometimes fails and will return negative numbers that should be ignored
        info_hhbar = dh*prob_h*log(1./sigma2_hMLE');
        info_asymp_nats(i_eps, i_eta) = info_hhbar;
        info_hhbar2 = 2*dh*prob_h(1:2:end)*log(1./sigma2_hMLE(1:2:end)');
        info_asymp2_nats(i_eps, i_eta) = info_hhbar2;
        
        %% compute average difference between h_MLE and h_true
        err_MLE(i_eps, i_eta) = mean((h_MLE_pars(:, 1) - h_true_vals(:)).^2);
        
        h_MLE_arr{i_eps, i_eta} = {sigma2_hMLE, h_MLE_pars};
        
    end
end
%% 
wtd_err_MLE = cellfun(@(x) sum(dh*prob_h(:).*((x{2}(:, 1) - h_true_vals(:)).^2)), ...
    h_MLE_arr);

%% For fixed eta, epsilon, h: compute the spike count distribution and information
% This will be smaller than the MLE h* information when the MLE approach
% works, so there it is an underestimate of the total information. Where the
% MLE breaks down, the population spike count information (which can be
% calculated directly, without approximation) is a better measure of how
% well a latent variable could be inferred. 

plot_dir = 'avalanches_matlab_code/plots/information_spikecount/';
spk_ct_at_h = @(h, eta1, eps1, num_obs) sum(rand(n_neur, num_obs) < repmat(prob_s_fun(h, J_i, eta1, eps1), 1, num_obs), 1);
num_obs = 10000;
file_name = ['nNeur' num2str(n_neur) '_Nsamp' num2str(num_obs)];
save(['avalanches_matlab_code/analysis_results/spike_count_info_' file_name '.mat']);
make_plots = false;

info_arr_SC = cell(length(eps_vals), length(eta_vals));
%%
for i_eta = 1:length(eta_vals)
    for i_eps = 1:length(eps_vals)
        eps0 = eps_vals(i_eps);
        eta0 = eta_vals(i_eta);
        plot_tag = ['_eta' num2str(eta0*10, '%01.0f') '_eps' num2str(eps0*10, '%01.0f') ];

        sc_v_h = zeros(length(h_true_vals), num_obs);

        sc_bins = bins2edges(0:n_neur);
        hist_sc_v_h = zeros(length(h_true_vals), length(sc_bins) - 1);

        sc_v_h_approx = zeros(length(h_true_vals), length(sc_bins) - 1);


        parfor i_h = 1:length(h_true_vals)
            sc = spk_ct_at_h(h_true_vals(i_h), eta0, eps0, num_obs);
            sc_v_h(i_h, :) = sc;
            hist_sc_v_h(i_h, :) = histcounts(sc_v_h(i_h, :), sc_bins, 'Normalization','pdf');
            mu_est = mean(sc);
            sigma_est = std(sc);
            sc_v_h_approx(i_h, :) = normpdf(0:n_neur, mu_est, sigma_est);
        end

        % compute information: could probably improve drastically by
        % implmeneting a semi-parametric approximation (extract N, p for
        % each value of h, then evaluate the binomial distribution in the
        % information calculation. 
        p_SC_tot = dh*prob_h*hist_sc_v_h;

        H_SC = sum(-p_SC_tot.*log2(p_SC_tot + eps));

        H_SC_v_h = sum(-hist_sc_v_h.*log2(hist_sc_v_h + eps), 2);
        H_SC_v_h_approx = sum(-sc_v_h_approx.*log2(sc_v_h_approx + eps), 2);

        info_est_ee = H_SC - dh*prob_h*H_SC_v_h;
        info_est_approx_ee = H_SC - dh*prob_h*H_SC_v_h_approx;

        info_arr_SC{i_eps, i_eta}.sampled_hist_sc_v_h = hist_sc_v_h;
        info_arr_SC{i_eps, i_eta}.approx_hist_sc_v_h = sc_v_h_approx;
        info_arr_SC{i_eps, i_eta}.sampled_p_SC_tot = p_SC_tot;
        info_arr_SC{i_eps, i_eta}.sampled_H_SC = H_SC;
        info_arr_SC{i_eps, i_eta}.sampled_H_SC_v_h = H_SC_v_h;
        info_arr_SC{i_eps, i_eta}.approx_H_SC_v_h = H_SC_v_h_approx;
        info_arr_SC{i_eps, i_eta}.sampled_info_est_ee = info_est_ee;
        info_arr_SC{i_eps, i_eta}.sampled_info_approx_ee = info_est_approx_ee;

        if make_plots
            h_Fig = makeMyFigure(20, 20);
            nexttile([2 2])
            imagesc([], h_true_vals, hist_sc_v_h)
            xlabel('spike count')
            ylabel('h')
            title(['P(SC | h), \epsilon = ' num2str(eps0) ', \eta = ' num2str(eta0) ', I = ' num2str(info_est_ee, '%1.2f')])

            nexttile([2 1])
            hold on
            plot(H_SC_v_h, h_true_vals)
            plot(H_SC_v_h_approx, h_true_vals)
            xlabel('est. entropy (SC | h)')
            ylabel('h')

            nexttile([1 2])
            bar(0:n_neur, p_SC_tot, 'BarWidth',1)
            xlabel('spike count')
            ylabel('p(SC)')
            title(['Estimated entropy: ' num2str(H_SC, '%1.2f')])

            print(h_Fig,'-dpdf', [plot_dir file_name plot_tag]);
        end
        
    end
%     close all
    % periodic save
    save(['avalanches_matlab_code/analysis_results/spike_count_info_' file_name '.mat'], ...
        'info_arr_SC', '-append');
end



%% Compute p_ava, p_silence: (n_neur set at top of script, should be 128)

eta_list_ultra = linspace(0, 10, 101);
eps_list_ultra = linspace(0, -14, 141);
[p_sil_vs_ee, p_ava_vs_ee, ave_p_silence, ave_p_ava, h_vals] = ...
    computeProbSilenceAvalanche(eta_list_ultra, eps_list_ultra, n_neur);
% %% Compute p_silence and P_ava for 8-cell groups
% num_draws = 2000;
% [p_sil_vs_ee8, p_ava_vs_ee8, ave_p_silence8, ave_p_ava8] = computeProbSilenceAvalanche(eta_list_ultra, eps_list_ultra, 8, num_draws);
% %% Compute p_silence and P_ava for 1024-cell groups
% num_draws1k = 8;
% [p_sil_vs_ee1024, p_ava_vs_ee1024, ave_p_silence1024, ave_p_ava1024] = computeProbSilenceAvalanche(eta_list_ultra, eps_list_ultra, 1024, num_draws1k);

%% interpolate
p_silence_ee = interp2(eta_list_ultra', eps_list_ultra, ave_p_silence, eta_vals_sf1', eps_vals_sf1);
p_ava_ee = interp2(eta_list_ultra', eps_list_ultra, ave_p_ava, eta_vals_sf1', eps_vals_sf1);

%% Compute variance in p_ava
prob_h = normpdf(h_vals);
dh = h_vals(2) - h_vals(1);
p_ava_squared = sum((p_ava_vs_ee.^2).*repmat(reshape(prob_h*dh, [1 1 length(prob_h)]), [size(p_ava_vs_ee, [1 2]) 1]), 3);
p_ava_average = sum((p_ava_vs_ee).*repmat(reshape(prob_h*dh, [1 1 length(prob_h)]), [size(p_ava_vs_ee, [1 2]) 1]), 3);

sigma_P = sqrt(p_ava_squared - p_ava_average.^2);
%% Save everything, but smaller info array

info_arr_SC = cellfun(@(x) rmfield(x, {'approx_hist_sc_v_h', 'sampled_hist_sc_v_h', 'sampled_p_SC_tot'}), info_arr_SC);
save(['avalanches_matlab_code/analysis_results/spike_count_probAva_info_' file_name '.mat'], '-v7.3')


%% helper functions

function [h_MLE, lh_expr, rh_val] = find_hMLE(s_i_vals, J_i, eta_val, eps_val)
% semi-analytic calculation of MLE value of h:
% given model parameters and observation(s) s_i_vals, 
% returns maximum likelihood value of h. 
% 
% found by solving 
% \sum_i(\frac{\phi J_i}{1 + \exp(\phi J_i h + \epsilon)} = \sum_i \phi s_i J_i
% for h. 
% 
% Inputs:
%       s_i_vals: num_neurons by num_observations matrix (in interval
%           [0,1]; single observation will be binary but function also
%           works for average over observations. 
%       J_i: num_neurons by 1 coupling coefficients between latent field
%           and neural observations
%       eps_val: scalar, bias toward silence/firing 
%       phi_val: scalar, multiplies Jhs term
%
% Output: 
%       h_MLE: scalar, max-likelihood estimate of h under the single
%           latent field model 
%%
    lh_expr = @(h) sum(eta_val*J_i./(1 + exp(eta_val*J_i*h + eps_val)));
    rh_val = sum(eta_val*mean(s_i_vals, 2).*J_i);
    
    mle_fun = @(h) lh_expr(h) - rh_val;
    
    % finds where mle_fun(h) = 0, start from initial value of h = 0. 
    h_MLE = fsolve(mle_fun, 0);
        

end

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