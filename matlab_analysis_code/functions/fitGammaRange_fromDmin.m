function [gamma_fit, gamma_range] = fitGammaRange_fromDmin(dur_bins, ave_sizes, d_min)
% Function for adaptively fitting the size ~ duration scaling curve
% Input strut x_reps_str is loaded from the results directory; provide the
% desired parameter values through eta_val, eps_val, tau_val, and
% rep_string. x_reps_str should be a struct, not a cell array of structs.

f_opt = fitoptions('Method', 'linearLeastSquares');
% goal is to find how far a fit that initiates at d_min (provided in linear,
% not log, units) predicts the data
%%
x = log10(dur_bins);
y = log10(ave_sizes);

x_bins = 0.5:0.5:ceil(max(x));


%%

% Calculate the slope and range if we start at d_min
xMin = log10(d_min);
minStart_end_x = findRangeGoodFitFromStart(x, y, xMin);
minS_fit_range = x > xMin & x <= minStart_end_x;
%%
fit_minS = fit(x(minS_fit_range)', y(minS_fit_range)', 'poly1', f_opt);
gamma_fit = fit_minS;
gamma_range = [xMin minStart_end_x];



end

function [new_end_x, ind_end, fit0] = findRangeGoodFitFromStart(x, y, x_start)
% Force a start at x_start; keep a fit based on the first decade. What is
% the upper limit on the fit range? 
%%
f_opt = fitoptions('Method', 'linearLeastSquares');
x_end_vals = (x_start+0.5):0.05:max(y);
first_fit_range = x > x_start & x <= x_start+1;

[fit0, gof] = fit(x(first_fit_range)', y(first_fit_range)', 'poly1', f_opt);
%     e_m_val(i_e) = fit0.p1;
%     e_b_val(i_e) = fit0.p2;
%     ci_fit0 = confint(fit0);
%     e_m_CI(i_e, :) = ci_fit0(:, 1)';
%     e_b_CI(i_e, :) = ci_fit0(:, 2)';


% keep track of the original fraction of points in the prediction
% interval
new_end_pred_frac = zeros(length(x_end_vals), 1);
min_pred_frac = zeros(length(x_end_vals), 1);
% e_m_val = zeros(length(x_end_vals), 1);
% e_m_CI = zeros(length(x_end_vals), 2);
% e_b_val = zeros(length(x_end_vals), 1);
% e_b_CI = zeros(length(x_end_vals), 2);
%%
% fit error shift
e_resid_dec1 = zeros(length(x_end_vals), 2);
for i_e = 1:length(x_end_vals)
    extended_range = x > (x_end_vals(i_e)-0.5) & x <= x_end_vals(i_e);

    % calculate prediction on last half decade of extended range
    y_pred = predint(fit0, x(extended_range)');

    % calculate if value is in the predicted interval
    in_predint = y(extended_range)' > y_pred(:, 1) & y(extended_range)' < y_pred(:, 2);

    % compute fraction of observations in range 
    new_end_pred_frac(i_e) = mean(in_predint);
    min_pred_frac(i_e) = floor(0.8*sum(extended_range))./sum(extended_range);

end
%%

% now take maximum range: most (95%?) of data within the 95% CI -- too stringent?
% try 80%
[new_end_x, ind_end] = max(x_end_vals(new_end_pred_frac > min_pred_frac));

end