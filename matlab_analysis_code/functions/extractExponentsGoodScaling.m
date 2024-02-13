function [ave_exp_vals, frac_scale, has_scaling] = extractExponentsGoodScaling(sum_stats, size_or_dur, max_ava_size, max_z_val)
% Takes the sum_stats struct and filters the desired expoennts (tau or
% alpha) by the quality of scaling metrics, using the criteria provided 
% for "size" returns tau values; for "dur", returns alpha values. 

    if strcmp(size_or_dur, 'size')
        s_min = sum_stats.tau_x_min;
        surr_z = (sum_stats.tau_KS_min - sum_stats.tau_medSurrKS)./sum_stats.tau_stdSurrKS;

        exp_vals = sum_stats.tau_values;

    elseif strcmp(size_or_dur, 'dur')
        s_min = sum_stats.alpha_x_min;
        surr_z = (sum_stats.alpha_KS_min - sum_stats.alpha_medSurrKS)./sum_stats.alpha_stdSurrKS;
        exp_vals = sum_stats.alpha_values;

    elseif strcmp(size_or_dur(1:5), 'gamma')
        
        s_min = (sum_stats.alpha_x_min + sum_stats.tau_x_min)/2;
        surr_z_t = (sum_stats.tau_KS_min - sum_stats.tau_medSurrKS)./sum_stats.tau_stdSurrKS;
        surr_z_a = (sum_stats.alpha_KS_min - sum_stats.alpha_medSurrKS)./sum_stats.alpha_stdSurrKS;
        surr_z = (surr_z_a + surr_z_t)/2;
        if strcmp(size_or_dur, 'gamma_pred')
            exp_vals = sum_stats.gamma_pred_values;
        elseif strcmp(size_or_dur, 'gamma_fit')
            exp_vals = sum_stats.gamma_fit_values;
        end
    else
        ave_exp_vals = [];
        frac_scale = [];
        disp('size_or_dur must be "size", "dur", or "gamma"')
        return;
    end

    has_scaling = s_min < max_ava_size & surr_z < max_z_val;
    has_scaling(isnan(has_scaling)) = false;
    frac_scale = mean(has_scaling, 3);
    
    exp_vals(~has_scaling) = nan;
    ave_exp_vals = mean(exp_vals, 3, 'omitnan');

end