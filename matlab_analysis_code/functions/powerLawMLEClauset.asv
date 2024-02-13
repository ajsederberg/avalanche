function a_MLE = powerLawMLEClauset(x_i, x_min)
% returns the MLE of a power law exponent for the observations x_i,
% restricted to observations larger than x_min. If x_min is empty,
% estimates x_min using the Clauset methods
% x_i are assumed to be DISCRETE

% UPDATE 1/18/22: uses a look-up table (as a global variable), 
% glu_dhurwitzZeta(1, a, x_min) and
% glu_hurwitzZeta(a, x_min) to evaluate these functions .
    if x_min == 1 || isempty(x_min)
        
        meanLogX = -mean(log(x_i (x_i >= 1)));
        
        alpha_fun = @(a) zeta(1, a) - meanLogX*zeta(a);
        [a_MLE, fvalF, exitflag] = fsolve(@(x) alpha_fun(x), 1.5);
    elseif x_min <= 6
        %%
        [glu_hurwitzZeta, glu_dhurwitzZeta] = loadhurwitzZetaLookup();


        meanLogX = -mean(log(x_i (x_i >= x_min)));

%         alpha_fun = @(a) hurwitzZeta(1, a, x_min) - meanLogX*hurwitzZeta(a, x_min);
        alpha_fun = @(a) glu_dhurwitzZeta(a, x_min) - meanLogX*glu_hurwitzZeta(a, x_min);
        [a_MLE, fvalF, exitflag] = fsolve(@(x) alpha_fun(x), 1.5);
    else
        % use the approximate function for x_min > 6 (faster)
        x_i_samples = x_i(x_i >= x_min);
        a_MLE = 1 + 1/mean(log(x_i_samples./(x_min - 0.5)));
        
    end

end