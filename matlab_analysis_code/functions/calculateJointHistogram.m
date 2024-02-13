function [n_cts, size_bins, dur_bins] = calculateJointHistogram(sizes, durations)
% This calls histcounts2 using bins for discrete data: bins are all unique
% sizes and durations from the observed data. 

    % set bins
    dur_bins = unique(double(durations));
    dur_edges = [dur_bins - 0.5 dur_bins(end)+0.5];
    size_bins = unique(sizes);
    size_edges = [size_bins - 0.5 size_bins(end)+0.5];
    
    % tally observations
    [n_cts, ~, ~] = histcounts2(sizes, double(durations), ...
        size_edges, dur_edges, 'normalization', 'count');

end