function [ks_D, n_obs_tail, ks_stat] = computeKSStatPowerLaw(x_i, a_val, x_min, x_max_vals)
% Computes the KS statistic (maximum absolute difference between
% cumulative density function) for model and observations
% ks_D is eqn 3.9 from Clauset: this is the KS statistic for x > x_min ONLY

%   

    max_n_possible = sum(x_i >= x_min & x_i <= max(x_max_vals));
    if max_n_possible < 500
        ks_stat = nan(size(x_max_vals));
        ks_D = ks_stat;
        n_obs_tail = zeros(size(x_max_vals));
        for ii = 1:length(x_max_vals)
            n_obs_tail(ii) = sum(x_i >= x_min & x_i <=x_max_vals(ii));
        end
        
    else
        % compute the model CDF at 50 log-spaced point, then interpolate
        % up to x_bins - otherwise very inefficient. 
        x_bins = x_min:max(x_max_vals);
        x_bins_eval = x_bins;
%         x_bins_eval(x_bins > 1e4) = 1e4;    
        cdf_model = 1 - computeCDFpowerlawdistribution(x_bins_eval, a_val, x_min);
        % This ^^^ is a hack for counts > 1e4; should adjust to have
        % power-law behavior if it matters 
    %     %% with 50 points, error relative to exact is < 1e-6. OK. 
    %     tic
    %     x_bin_log = unique(round(logspace(log10(min(x_bins)), log10(max(x_bins)), 50)));
    %     cdf_model_eff = computeCDFpowerlawdistribution(x_bin_log, a_val, x_min);
    %     toc

    %     cdf_model = 1 - interp1(x_bin_log, cdf_model_eff, x_bins, 'spline');
        %% construct the semi-parametric model for all values of x_max in
        % x_max_vals
        ks_D = nan(length(x_max_vals), 1);
        ks_stat = nan(length(x_max_vals), 1);
        n_obs_tail = zeros(length(x_max_vals), 1);

        for ii = 1:length(x_max_vals)
            x_max = x_max_vals(ii);

            % compute the empirical CDF of data over the entire range
            x_full = 1:x_max;
            x_efull = bins2edges(x_full);
    %     x_ebins = bins2edges(x_bins);
            cdf_empirical = histcounts(x_i, x_efull, 'normalization', 'cdf');
            cdf_emp_powerlaw = histcounts(x_i(x_i >= x_min & x_i <= x_max), x_efull, 'Normalization', 'cdf');

            n_obs_tail(ii) = sum(x_i >= x_min & x_i <= x_max);
            frac_tail = n_obs_tail(ii)/sum(x_i <= x_max);
            if n_obs_tail(ii)  > 500
                % semi-parametric model is the empirical CDF up to x_min, then
                % the power-law cdf above x_min
                cdf_SP_model = 0*cdf_empirical;
                cdf_SP_model(x_full < x_min) = cdf_empirical(x_full < x_min);
                % assign semi-parametric model cdf values to the power-law prediction;
                % (1-frac_tail) and frac_tail ensure continuity with cdf up to x_min
                % (i.e. CDF ends at 1)
                cdf_SP_model(x_full >= x_min & x_full < x_max) = ...
                    (1-frac_tail) + frac_tail*cdf_model(x_bins >= x_min & x_bins < x_max);
                % set last two entries fo cdf_SP_model equal
                cdf_SP_model(end) = cdf_SP_model(end-1);

                ks_stat(ii) = max(abs(cdf_empirical - cdf_SP_model));

                % compute the ks_D, the ks-statistic over the power-law portion
                % only. note that first entry of cdf_model is 0, so start with
                % second (probability of having x_min)
                ks_D(ii) = max(abs(cdf_model(x_bins > x_min & x_bins <= x_max) - ...
                    cdf_emp_powerlaw(x_full >= x_min & x_full < x_max)));


            end
        end
    end
end