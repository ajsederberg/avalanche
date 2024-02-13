function [gamma_fit, gamma_range, minS_gamma_fit, minS_gamma_range] = fitGammaRange(dur_bins, ave_sizes, plot_diagnostics, plot_name)
% Function for adaptively fitting the size ~ duration scaling curve
% Input strut x_reps_str is loaded from the results directory; provide the
% desired parameter values through eta_val, eps_val, tau_val, and
% rep_string. x_reps_str should be a struct, not a cell array of structs. 


f_opt = fitoptions('Method', 'linearLeastSquares');
% goal is to find the broadest range over which 95% of the data falls
% within the 95% CI of the fit; also want to know the overall error in
% slope estimate
%%
if nargin < 3
    plot_diagnostics = false;
elseif plot_diagnostics
    makeMyFigure(30, 20);
    tiledlayout(2,3);

    ex_cmap = lines(2);
    if nargin < 4
        plot_name = datestr(now, 'yyyymmdd_HHMMSS');
    end
end

%%
    x = log10(dur_bins);
    y = log10(ave_sizes);

    x_starts = 0.2:0.1:(max(x)-1);
    x_bins = 0.5:0.5:ceil(max(x));


    % first, get the fit over each decade starting at x_ct(ii)
    fit_rmse = zeros(length(x_starts), 1);
    adj_rsq = zeros(length(x_starts), 1);
    m_val = zeros(length(x_starts), 1);
    b_val = zeros(length(x_starts), 1);
    m_CI = zeros(length(x_starts), 2);
    b_CI = zeros(length(x_starts), 2);
    invm_CI = zeros(length(x_starts), 2);
    pred_frac = nan(length(x_starts), length(x_bins));

    % consistent with fitting inverse (x = f(y))
    cons_inv = zeros(length(x_starts), 1);

    for ii = 1 : length(x_starts)
            
            fit_range = x <= (1 + x_starts(ii)) & x > x_starts(ii);
if sum(fit_range) > 10
            [fit0, gof] = fit(x(fit_range)', y(fit_range)', 'poly1', f_opt);

            fit_rmse(ii) = gof.rmse;
            adj_rsq(ii) = gof.adjrsquare;
            m_val(ii) = fit0.p1;
            b_val(ii) = fit0.p2;
            ci_fit0 = confint(fit0);
            m_CI(ii, :) = ci_fit0(:,1)';
            b_CI(ii, :) = ci_fit0(:, 2)';

            % check inverse fit
            [fitI, gofI] = fit(y(fit_range)', x(fit_range)', 'poly1', f_opt);
            % compute 1/m from inverse fit
            ci_fitI = 1./confint(fitI);
            invm_CI(ii, :) = ci_fitI(:,1)';

            cons_inv(ii) = ~ ( min(ci_fitI(:, 1)) > max(ci_fit0(:, 1)) || ... 
                max(ci_fitI(:, 1)) < min(ci_fit0(:, 1)));

            % calculate prediction on entire range
            y_pred = predint(fit0, x');

            % calculate if value is in the predicted interval
            in_predint = y' > y_pred(:, 1) & y' < y_pred(:, 2);

            % calculate fraction inside predicted interval 
            z = binYbyX(x', double(in_predint), x_bins, false, true);

            % enter in pred_frac
            pred_frac(ii, :) = z;
            
end
    end



%%
%     % start with the decade with the longest 95%-CI prediction interval
%     % agreement for which the fit and inverse fit have consistent slopes


    % "length" is how many SUBSEQUENT values of x_start have fit coefficients within CI bounds for this
    % value
    l_95 = 0*cons_inv;
    for ii = 1:(length(l_95)-1)
        this_CI = m_CI(ii, :);
        kk = 1;
        
        % fit criterion is that m, b from the linear fits is consistent
        % over the range
        fit_crit = (~(m_CI(ii+kk, 2) < m_CI(ii, 1) || m_CI(ii+kk, 1) > m_CI(ii, 2))) && ...
            (~(b_CI(ii+kk, 2) < b_CI(ii, 1) || b_CI(ii+kk, 1) > b_CI(ii, 2))); 
        while ii + kk < length(m_CI) && fit_crit
            kk = kk + 1;

            % update fit_crit
            fit_crit = (~(m_CI(ii+kk, 2) < m_CI(ii, 1) || m_CI(ii+kk, 1) > m_CI(ii, 2))) && ...
                (~(b_CI(ii+kk, 2) < b_CI(ii, 1) || b_CI(ii+kk, 1) > b_CI(ii, 2))); 
    
        end
        l_95(ii) = kk;

    end
%%
    % set to 0 if the inverse fit was not consistent
    % This was an idea for a consistency check, but in practice, it
    % unfairly penalizes ranges where there's more spread in y at fixed x
    % (creating a bias in the inverse fit). Could be fixed by choosing
    % ranges differently, but for now let's drop it. (9/20/22)
%     l_95(~cons_inv) = 0;
    %%
    if any(l_95 > 0)
        [max_L, ind_start] = max(l_95);
        x_end = x_starts(ind_start(1)) + 1;
        x_start = x_starts(ind_start(1));

    else
        [~, ind_start] = max(sum(pred_frac, 2));
        x_end = x_starts(ind_start) + 1;
        x_start = x_starts(ind_start);
    end
%% methods plots

    if plot_diagnostics
        
        % plot a couple of fits
        nexttile
        %%
        fit_range1 =  x <= x_end & x > x_start;
        if x_start > 1.8
            x_ex2 = 0.5;
        else
            x_ex2 = x_start + 1.5;
        end
        [~, ind2] = min(abs(x_starts - x_ex2));
        fit_range2 = x <= (1 + x_starts(ind2)) & x > x_starts(ind2);
        lh = zeros(5,1);

        hold on
        lh(1) = plot(x, y, '.', 'Color', 0.5*[1 1 1]);
        lh(2) = plot(x(fit_range1), y(fit_range1), '.', 'Color',ex_cmap(1, :));
        lh(3) = plot(x(fit_range2), y(fit_range2), '.', 'Color',ex_cmap(2, :));

        lh(4) = plot(x([1 end]), m_val(ind_start)*x([1 end]) + b_val(ind_start), ...
            'color', ex_cmap(1, :));
        lh(5) = plot(x([1 end]), m_val(ind2)*x([1 end]) + b_val(ind2), ...
            'color', ex_cmap(2, :));
        xlabel('duration (log scale)')
        ylabel('average size (log scale)')
        set(gca, 'Color', 'none')
        axis square
        legend(lh, {'all durations', 'ex. decade 1', 'ex. decade 2', ... 
            'fit to dec. 1', 'fit to dec. 2'}, 'Location','northoutside', 'NumColumns', 2)
        title('Part 1: select lower cutoff (d_{min}) for S ~ D^{\gamma}')
%%
        nexttile
        %% plot m and conf-int
        cla
%         yyaxis left
        hold on
        plot(x_starts, m_val, 'k')
        plot(x_starts, m_CI, 'color', 0.5*[1 1 1])
        xlabel('d_{min} (lower bound on fitted decade)')
        ylabel('slope (\gamma, 95% CI)')
        set(gca, 'color', 'none')
        axis square

        yyaxis right
        hold on
        plot(x_starts, b_val, 'k-.')
        plot(x_starts, b_CI, ':', 'color', 0.5*[1 1 1])
        ylabel('intercept (b, 95% CI)')
        title('fit parameters for single-decade fits')

        %% plot "length" 
        nexttile
        plot(x_starts, l_95)
        xlabel('d_{min}')
        ylabel('number of consistent start points')
        title('fit range starting from each d_{min}')
        set(gca, 'color', 'none')
        axis square
        suptitle('Estimate of \gamma value and power-law range')
    end 
%%
% or start with the fit criterion? 
% [~, ind_start] = min(fit_rmse);
% x_start = x_starts(ind_start);
% x_end = x_start + 1;
% or from where the confidence intervals on the slope parameter are
% consistent? (overlapping x values, though)
%% now extend the range as much as possible
x_end_vals = x_end:0.05:max(x);
[new_end_x, ind_end, e_m_val, e_m_CI, e_b_val, e_b_CI, e_resid_dec1, fits_start] = extendFitRangeFromStart(x, y, x_start, x_end_vals, x_bins);


    
%     % keep track of the original fraction of points in the prediction
%     % interval
%     x_end_vals = x_end:0.05:max(x);
%     new_end_pred_frac = zeros(length(x_end_vals), length(x_bins));
%     e_m_val = zeros(length(x_end_vals), 1);
%     e_m_CI = zeros(length(x_end_vals), 2);
%     e_b_val = zeros(length(x_end_vals), 1);
%     e_b_CI = zeros(length(x_end_vals), 2);
%         
%     % end extend
% %     e_invm_CI = zeros(length(x_end_vals), 2);
% %     e_invm_CI = zeros(length(x_end_vals), 2);
% %     e_cons_inv = zeros(length(x_end_vals), 1);
% % 
% 
%     % fit error shift
%     e_resid_dec1 = zeros(length(x_end_vals), 2);
%     for i_e = 1:length(x_end_vals) 
%         new_fit_range = x > x_start & x <= x_end_vals(i_e);
%     
%         [fit0, gof] = fit(x(new_fit_range)', y(new_fit_range)', 'poly1', f_opt);
% %         [fitI, gofI] = fit(y(new_fit_range)', x(new_fit_range)', 'poly1', f_opt);
%     
%         e_m_val(i_e) = fit0.p1;
%         e_b_val(i_e) = fit0.p2;
%         ci_fit0 = confint(fit0);
%         e_m_CI(i_e, :) = ci_fit0(:, 1)';
%         e_b_CI(i_e, :) = ci_fit0(:, 2)';
% 
%          % compute 1/m from inverse fit
% %         ci_fitI = 1./confint(fitI);
% %         invm_CI(i_e, :) = ci_fitI(:,1)';
% 
% %         e_cons_inv(i_e) = ~ ( min(ci_fitI(:, 1)) > max(ci_fit0(:, 1)) || ... 
% %             max(ci_fitI(:, 1)) < min(ci_fit0(:, 1)));
% % 
% 
%         % calculate prediction on entire range
%         y_pred = predint(fit0, x');
%     
%         % calculate if value is in the predicted interval
%         in_predint = y' > y_pred(:, 1) & y' < y_pred(:, 2);
%     
%         % calculate fraction inside predicted interval 
%         z = binYbyX(x', double(in_predint), x_bins, false, true);
% 
%         % enter in pred_frac
%         new_end_pred_frac(i_e, :) = z;
% 
%         %% check for shift in fit error over the first decade of the fit
%         fit_resid = y(new_fit_range)' - fit0(x(new_fit_range)');
%         first_decade = x(new_fit_range) <= x_start + 1;
%         e_resid_dec1(i_e, 1) = mean(fit_resid(first_decade));
%         e_resid_dec1(i_e, 2) = std(fit_resid(first_decade))/sqrt(sum(first_decade) -1 );
%     end
% 
%     %% now choose longest range 
% 
%     % count up observations in each bin
%     num_x_obs = histcounts(x, [0 x_bins]);
% 
%     min_pred_frac = floor(0.95*num_x_obs)./num_x_obs;
%     range_x = sum(new_end_pred_frac > repmat(min_pred_frac, size(new_end_pred_frac, 1), 1), 2);
% 
%     % require fit at low end of range 
%     fits_start = abs(e_resid_dec1(:, 1))./e_resid_dec1(:, 2) < 3;
%     range_x(~fits_start) = 0;
% 
%     % if inverse fit is not consistent, set to 0
% %     range_x(e_cons_inv == 0) = 0;
% 
%     % now take maximum range
%     [new_end_x, ind_end] = max(x_end_vals(range_x == max(range_x)));

    %% fit over this range
    final_fit_range = x > x_start & x <= new_end_x;
    %%
    [fit_final, gof_final] = fit(x(final_fit_range)', y(final_fit_range)', 'poly1', f_opt);
    gamma_fit = fit_final;
    gamma_range = [x_start new_end_x];
    %%

    % also calculate the slope and range if we start at fixed lower bound (0.5)
    xMin = 0.5;
    minStart_end_x = extendFitRangeFromStart(x, y, xMin, (1+xMin):.05:max(x), x_bins);
    minS_fit_range = x > xMin & x <= minStart_end_x;
    %%
    fit_minS = fit(x(minS_fit_range)', y(minS_fit_range)', 'poly1', f_opt);
    minS_gamma_fit = fit_minS;
    minS_gamma_range = [xMin minStart_end_x];

    if plot_diagnostics

        %% 
        nexttile
        %%
        fit_range1 =  x <= new_end_x & x > x_start;
        if new_end_x > x_start + 1.5
            x_ex2 = x_start + 1;
        else
            x_ex2 = min(max(x), x_start + 1.2);
        end
        [~, ind2] = min(abs(x_end_vals - x_ex2));
        fit_range2 = x <= (x_end_vals(ind2)) & x > x_start;
        lh = zeros(6,1);

        hold on
        lh(1) = plot(x, y, '.', 'Color', 0.5*[1 1 1]);
        lh(2) = plot(x(fit_range1), y(fit_range1), '.', 'Color',ex_cmap(1, :));
        lh(3) = plot(x(fit_range2), y(fit_range2), '.', 'Color',ex_cmap(2, :));

        lh(4) = plot(x([1 end]), e_m_val(ind_end)*x([1 end]) + e_b_val(ind_end), ...
            'color', ex_cmap(1, :));
        lh(5) = plot(x([1 end]), e_m_val(ind2)*x([1 end]) + e_b_val(ind2), ...
            'color', ex_cmap(2, :));
        axis square
        lh(6) = plot(x_start*[1 1], ylim, 'k-');
        xlabel('duration (log scale)')
        ylabel('average size (log scale)')
        set(gca, 'Color', 'none')
        legend(lh, {'all durations', ['(i) d_{max} = ' num2str(new_end_x)],...
            ['(ii) d_{max} = ' num2str(x_ex2)],...
            'fit to (i)', 'fit to (ii)', 'd_{min}'}, ...
            'Location', 'NorthOutside','NumColumns', 2, 'color', 'none')

        title('Part 2: select upper cutoff (d_{max}) for S ~ D^{\gamma}')
        nexttile
        %%
        cla;    
        hold on
        errorbar(x_end_vals, e_resid_dec1(:, 1), e_resid_dec1(:, 2), ...
            'color', 0.5*[1 1 1])
        errorbar(x_end_vals(fits_start), e_resid_dec1(fits_start, 1), e_resid_dec1(fits_start, 2), ...
            'color', 0*[1 1 1], 'linewidth', 1)
        plot(new_end_x*[1 1], ylim, 'k')
        legend({'all end points', 'within 95% CI', 'selected d_{max}'}, 'Location','northoutside', 'NumColumns',2)
        xlabel('d_{max}')
        ylabel(['ave error on [' num2str([x_start x_start+1]) '] ( +/- 1 SE)'])
        title('residuals for fits ending at d_{max}')
        set(gca, 'color', 'none')
        axis square

        nexttile
        %%
        lh = zeros(3,1);
        cla
        hold on
        lh(1) = plot(x, y, '.', 'Color', 0.5*[1 1 1]);
        lh(2) = plot(x(final_fit_range), y(final_fit_range), '.', 'Color',ex_cmap(1, :));
        lh(3) = plot(x([1 end]), fit_final.p1*x([1 end]) + fit_final.p2, ...
            'color', ex_cmap(1, :));
        set(gca, 'color', 'none')
        xlabel('duration (log scale)')
        ylabel('average size (log scale)')
        axis square
        legend(lh, {'all', 'final fit range', 'final fit'}, 'Color', 'none', 'Location','northoutside', 'NumColumns', 2  )
        title(['\gamma = ' num2str(fit_final.p1, '%1.2f') ', range: ' num2str(diff(gamma_range))])

        print(gcf, '-dpdf', ['avalanches_matlab_code/plots/gamma_fit_diagnostics/' ... 
            'gamma_fit_' plot_name])
        close(gcf);
    end

end

function [new_end_x, ind_end, e_m_val, e_m_CI, e_b_val, e_b_CI, e_resid_dec1, fits_start] = extendFitRangeFromStart(x, y, x_start, x_end_vals, x_bins)
    
    f_opt = fitoptions('Method', 'linearLeastSquares');

    % keep track of the original fraction of points in the prediction
    % interval
    new_end_pred_frac = zeros(length(x_end_vals), length(x_bins));
    e_m_val = zeros(length(x_end_vals), 1);
    e_m_CI = zeros(length(x_end_vals), 2);
    e_b_val = zeros(length(x_end_vals), 1);
    e_b_CI = zeros(length(x_end_vals), 2);
        
    % end extend
%     e_invm_CI = zeros(length(x_end_vals), 2);
%     e_invm_CI = zeros(length(x_end_vals), 2);
%     e_cons_inv = zeros(length(x_end_vals), 1);
% 

    % fit error shift
    e_resid_dec1 = zeros(length(x_end_vals), 2);
    for i_e = 1:length(x_end_vals) 
        new_fit_range = x > x_start & x <= x_end_vals(i_e);
    
        [fit0, gof] = fit(x(new_fit_range)', y(new_fit_range)', 'poly1', f_opt);
%         [fitI, gofI] = fit(y(new_fit_range)', x(new_fit_range)', 'poly1', f_opt);
    
        e_m_val(i_e) = fit0.p1;
        e_b_val(i_e) = fit0.p2;
        ci_fit0 = confint(fit0);
        e_m_CI(i_e, :) = ci_fit0(:, 1)';
        e_b_CI(i_e, :) = ci_fit0(:, 2)';

         % compute 1/m from inverse fit
%         ci_fitI = 1./confint(fitI);
%         invm_CI(i_e, :) = ci_fitI(:,1)';

%         e_cons_inv(i_e) = ~ ( min(ci_fitI(:, 1)) > max(ci_fit0(:, 1)) || ... 
%             max(ci_fitI(:, 1)) < min(ci_fit0(:, 1)));
% 

        % calculate prediction on entire range
        y_pred = predint(fit0, x');
    
        % calculate if value is in the predicted interval
        in_predint = y' > y_pred(:, 1) & y' < y_pred(:, 2);
    
        % calculate fraction inside predicted interval 
        z = binYbyX(x', double(in_predint), x_bins, false, true);

        % enter in pred_frac
        new_end_pred_frac(i_e, :) = z;

        %% check for shift in fit error over the first decade of the fit
        fit_resid = y(new_fit_range)' - fit0(x(new_fit_range)');
        first_decade = x(new_fit_range) <= x_start + 1;
        e_resid_dec1(i_e, 1) = mean(fit_resid(first_decade));
        e_resid_dec1(i_e, 2) = std(fit_resid(first_decade))/sqrt(sum(first_decade) -1 );
    end

    %% now choose longest range 

    % count up observations in each bin
    num_x_obs = histcounts(x, [0 x_bins]);

    min_pred_frac = floor(0.95*num_x_obs)./num_x_obs;
    range_x = sum(new_end_pred_frac > repmat(min_pred_frac, size(new_end_pred_frac, 1), 1), 2);

    % require fit at low end of range 
    fits_start = abs(e_resid_dec1(:, 1))./e_resid_dec1(:, 2) < 3;
    range_x(~fits_start) = 0;

    % if inverse fit is not consistent, set to 0
%     range_x(e_cons_inv == 0) = 0;

    % now take maximum range
    [new_end_x, ind_end] = max(x_end_vals(range_x == max(range_x)));

end