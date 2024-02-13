function [ave_sizes, se_sizes, dur_bins, num_obs] = calculateAveSizeByDuration(sizes, durations)
% Calculates the average avalanche size by duration of avalanche. Allows
% creating duration bins larger than 1 (default dur_bin_size is 1). The
% min_obs_cutoff drops all duration bins with fewer than the specified
% number of observations. 
%
    % calculated joint counts, returns bin centers. bin width is 1. 
    [n_cts, size_bins, dur_bins] = calculateJointHistogram(sizes, durations);

    % normalize to get p( size | duration)
    n_cts_norm = n_cts*diag(1./sum(n_cts, 1));
        
    % calculate < size p(size | duration) > 
    ave_sizes = size_bins*n_cts_norm;
    
    % number of observations of each duration
    num_obs = sum(n_cts, 1);

    % calculate expectations of (s_i - mu)^2 by evaluating sum(
    % (size - ave)^2 * p(size|duration)) (summation over duration values. 
    var_sizes_bias = sum(bsxfun(@(x,y) (x-y).^2, size_bins', ave_sizes).*n_cts_norm, 1);
    
    % denominator N-1 not N
    var_sizes = var_sizes_bias.*num_obs./(num_obs-1);
    se_sizes = sqrt(var_sizes);
end