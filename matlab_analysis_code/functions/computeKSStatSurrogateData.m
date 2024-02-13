function [ks_stats_D, surr_a_est] = computeKSStatSurrogateData(a_hat, x_min, x_max_val, n_total, n_tail, emp_bins, num_surr)
% generate surrogate datasets and measure the KS statistic between model
% and surrogate dataset
% ks_stats_D: over power law only (eqn 3.9)
% ks_stats_Surr: over entire region

% tail_bins are those greater than the lower cutoff x_min
    tail_bins = emp_bins(emp_bins >= x_min);
    powerLaw_ccdf = computeCDFpowerlawdistribution(tail_bins, a_hat, x_min);

    %%%%%%%%%%%%%%% THis is hack to avoid out of bound problems
%     tail_bins_eval = tail_bins;
% %     x_max_val = 1e4;
%     tail_bins_eval(tail_bins_eval > x_max_val) = x_max_val;    
%     %         cdf_model = 1 - computeCDFpowerlawdistribution(x_bins_eval, a_val, x_min);
% %%
% 
%     powerLaw_ccdf = computeCDFpowerlawdistribution(tail_bins_eval, a_hat, x_min);
%     pe4_val = computeCDFpowerlawdistribution(x_max_val, a_hat, x_min);
%     % super dumb bad hack to interpolate from 1e4 to max(tail_bins) - should at
%     % least do this in log-space to preserve power law 
%     powerLaw_ccdf(tail_bins > x_max_val) = interp1([x_max_val max(tail_bins)], [pe4_val 0], tail_bins(tail_bins > x_max_val));
%%
    % tail_bins_approx = unique(round(logspace(log10(tail_bins(1)), log10(tail_bins(end)), 201)));
    %%
%     scale_fac = 1 - empirical_cdf(emp_bins == x_min);
    %%
%     cdf_model(emp_bins >= x_min) = (1 - scale_fac) + scale_fac*(1 - powerLaw_ccdf);
    %%

%     cdf_preTail = empirical_cdf(emp_bins < x_min);
%     cdf_preTail = cdf_preTail/cdf_preTail(end);

%%
%     ks_stats_Surr = zeros(num_surr, 1);
    ks_stats_D = zeros(num_surr, 1);
    surr_a_est = zeros(num_surr, 1);
    tic
    for ii = 1:num_surr
%         rng('shuffle')
        % with probability n_tail/n_total, each surrogate observation is assigned
        % to the tail
        num_tail_surr = sum(rand(n_total, 1) < n_tail/n_total);

        % those not in the tail are in the "pre-tail"
    %     num_preTail_surr = n_total - num_tail_surr;
%%

        tail_draws = rand(num_tail_surr, 1);
        tail_x_i_draws = round(interp1(1 - powerLaw_ccdf, tail_bins, tail_draws, 'previous'));
        
     %%   
        % drop tail draws that are larger than the x_max_val
        tail_x_i_draws = tail_x_i_draws(tail_x_i_draws <= x_max_val);
        
        a_s_hat = powerLawMLEClauset(tail_x_i_draws, x_min);
        surr_a_est(ii) = a_s_hat;
        ks_stats_D(ii) = computeKSStatPowerLaw(tail_x_i_draws, a_s_hat, x_min, x_max_val); 

    %     pretail_draws = rand(num_preTail_surr, 1);
    %     preTail_x_i_draws = interp1([0; cdf_preTail], [0 pretail_bins], pretail_draws, 'next');
    % 
    %     surr_cdf = histcounts([preTail_x_i_draws; tail_x_i_draws], [emp_bins Inf], ... 
    %         'Normalization', 'cdf');
    %     ks_stats_Surr(ii) = max(abs(surr_cdf - cdf_model'));
    %     
    %     % over powerlaw only; interpolating from 100 log-spaced bins in the tail
    %     % is a numerical error ov about 10^-5 (for values checked, xmin = 3, a
    %     % = 1.9)
    %     surr_a_est(ii) = a_s_hat;
    %     
    %     est_powerLaw_ccdf = 1 - computeCDFpowerlawdistribution(tail_bins_eval, a_s_hat, x_min);
    % 
    % %     pL_ccdf_approx = 1 - computeCDFpowerlawdistribution(tail_bins_approx, a_s_hat, x_min);
    % %     est_powerLaw_ccdf = interp1(tail_bins_approx, pL_ccdf_approx, tail_bins, 'spline');
    %     
    % %     powerLaw_ccdf = 1 - computeCDFpowerlawdistribution(tail_bins, a_s_hat, x_min);
    %     surr_PL_cdf = histcounts(tail_x_i_draws, [0 tail_bins-.5], 'Normalization', 'cdf');
    %     ks_stats_D(ii) = max(abs(surr_PL_cdf - est_powerLaw_ccdf));

    if mod(ii, 100) == 0
        disp(['shuffle ' num2str(ii) ' done'])
    end
    end
end
