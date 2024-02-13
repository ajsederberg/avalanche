function [gamma_fit_mean, gamma_fit_se, gamma_range, full_res] = extractGammaFitRegion(pfit_G)
% starting at each point, how many decades of scaling can you have before
% the standard deviation of the mean gamma fits (scaled by sqrt(4), because
% we have overlapping quarters) is > 3 SE of an typical individual gamma
% SE? 

% skip if pfit_G is empty / return NaN
if isempty(pfit_G)
    gamma_fit_mean = nan;
    gamma_fit_se = nan;
    gamma_range = [nan nan];
    return;
end



g_means = cellfun(@(x) x.mhat, pfit_G);
g_ses = cellfun(@(x) x.mSE, pfit_G);

mean_gamma_over_range = 0*g_means; % if starting from this index, what is the average gamma
se_gamma_over_range = 0*g_ses; % if starting from this index, what is the uncertainty in gamma? 
gamma_end_ranges = 0*g_ses;    % if starting from each index, how many decades?
% start at 2 because pfit_G{1} is starting from a duration of 1 
for ii = 2:length(pfit_G)
    
    init_g = g_means(ii);
    init_se = g_ses(ii);
    
    end_range = ii;
    while end_range < length(pfit_G) && 2*std(init_g) < 3*init_se
        end_range = end_range + 1;
        init_g = g_means(ii:end_range);
        init_se = sqrt(mean(g_ses(ii:end_range).^2));
    end
    
    gamma_end_ranges(ii) = end_range;
    mean_gamma_over_range(ii) = mean(g_means(ii:end_range));
    se_gamma_over_range(ii) = std(g_means(ii:end_range));
end

% find the starting point from which you have the most decades of scaling
[~, ind_start] = max(gamma_end_ranges - (1:length(gamma_end_ranges))');

gamma_range = [pfit_G{ind_start}.x_cutoff(1) pfit_G{gamma_end_ranges(ind_start)}.x_cutoff(2)];
gamma_fit_mean = mean_gamma_over_range(ind_start(1));
gamma_fit_se = se_gamma_over_range(ind_start(1)); 

full_res.mean_gamma_over_range = mean_gamma_over_range;
full_res.se_gamma_over_range = se_gamma_over_range;
full_res.gamma_end_ranges = gamma_end_ranges;

%%%%%%%%%%%%%%%%%%%%% complicated version - but kind of works?
% % %%
% % g_fit_mat = zeros(length(pfit_G));
% % g_fit_var_mat = ones(length(pfit_G));
% % g_fit_frac_SE_mat = zeros(length(pfit_G));
% % 
% % g_means = cellfun(@(x) x.mhat, pfit_G);
% % g_ses = cellfun(@(x) x.mSE, pfit_G);
% % 
% % g_range = zeros(length(pfit_G));
% % 
% % for ii = 2:length(pfit_G)
% %     for jj = 2:ii
% %         
% %         % question for morning brain: what does it make sense to compare
% %         % here? the idea is that the variation in the means is comparable
% %         % to the error in any given mean. 
% %         
% %         % compute the mean over this interval and the average standard
% %         % error at each point over this interval
% %         ave_g = mean(g_means(jj:ii));
% %         mean_g_var = mean(g_ses(jj:ii).^2);
% %         
% %         % compute the variance across the mean fits
% %         var_g = var(g_means(jj:ii));
% %         
% %         g_fit_mat(ii,jj) = ave_g;
% %         g_fit_mat(jj,ii) = ave_g;
% %         
% %         g_fit_var_mat(ii, jj) = var_g + mean_g_var;
% %         g_fit_var_mat(jj, ii) = g_fit_var_mat(ii, jj);
% %         % compute fraction of total variation across segments 
% %         % that comes from single-point variability rather than cross-point
% %         % variability 
% %         g_fit_frac_SE_mat(ii, jj) = mean_g_var/g_fit_var_mat(ii, jj);
% %         g_fit_frac_SE_mat(jj, ii) = mean_g_var/g_fit_var_mat(ii, jj);
% %         
% %         % range in decades
% %         g_range(ii, jj) = 1 + abs(ii-jj)*0.25;
% %     end
% % end
% % 
% % % % first strategy: set criteria below
% % % %% criteria:
% % % % require overall SE of < 0.1 (g_fit_se_mat) and fraction of variance
% % % % attributed to fit error of greater then 50%
% % % % then maximize the range
% % % 
% % % allowed_points = g_fit_frac_SE_mat > 0.5 & g_fit_frac_SE_mat < 1 & g_fit_se_mat < 0.1;
% % % 
% % % max_range = max(g_range(allowed_points));
% % % 
% % % if ~isempty(max_range)
% % %     [i_max, j_max] = find(g_range.*allowed_points == max_range);
% % % else
% % %     % if no points satisfy the criterion, take the 80th percentile of
% % %     % g_fit_frac_SE_mat
% % %     g_min = prctile(g_fit_frac_SE_mat(g_fit_frac_SE_mat< 1), 80);
% % %     allowed_points = g_fit_frac_SE_mat > g_min & g_fit_se_mat < 0.1;
% % % 
% % %     max_range = max(g_range(allowed_points));
% % %     [i_max, j_max] = find(g_range.*allowed_points == max_range);
% % % 
% % % end
% % 
% % %% alternative strategy: for each possible range, find the highest
% % % 'fit_frac' value. Select the longest with where that is > 0.5. 
% %  
% % g_rng_vals = unique(nonzeros(g_range));
% % se_frac_rng = 0*g_rng_vals;
% % se_val_rng = 0*g_rng_vals;
% % 
% % for ii = 1:length(g_rng_vals)
% %     [se_frac_rng(ii), i_frac_max] = max(g_fit_frac_SE_mat(g_range == g_rng_vals(ii)));
% %     se_vals = g_fit_var_mat(g_range == g_rng_vals(ii));
% %     se_val_rng(ii) = se_vals(i_frac_max);
% % 
% % end
% % i_rng = find(se_frac_rng > 0.25, 1, 'last');
% % [i_best, j_best] = find(g_fit_frac_SE_mat == se_frac_rng(i_rng), 1);
% % % keyboard
% % % half_frac_range = 
% % %%
% % if ~isempty(i_best)
% % %     if length(i_max) > 1
% % %         keyboard
% % %     end
% %     % take first element of i_max, j_max; 
% % %     i_max = i_max(1);
% % %     j_max = j_max(1);
% %     
% %     gamma_fit_mean = g_fit_mat(i_best, j_best);
% %     gamma_fit_se = g_fit_var_mat(i_best, j_best);
% %     gamma_range = minmax( [pfit_G{i_best}.x_cutoff pfit_G{j_best}.x_cutoff]);
% % else
% %     gamma_fit_mean = nan;
% %     gamma_fit_se = nan;
% %     gamma_range = [0 0];
% % end