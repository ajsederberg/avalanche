function b_pdf = binopdf_approx(n, T, p)
    % shape variables 
    rep_n = repmat(reshape(n, numel(n), 1), 1, numel(p));
    rep_p = repmat(reshape(p, 1, numel(p)), numel(n), 1);
    
    b_pdf = 0*rep_p;
    % condition for normal approximation to be good
    n_a_inds = rep_n > 9*(1 - rep_p)./rep_p & rep_n > 9*(rep_p./(1 - rep_p));
    
    b_pdf(n_a_inds) = normpdf(rep_n(n_a_inds) - 0.5, rep_p(n_a_inds)*T, sqrt(T*(eps + rep_p(n_a_inds)).*(eps + 1 - rep_p(n_a_inds))));
    % if normal approximation is not accurate, use exact binomial
    b_pdf(~n_a_inds) = binopdf(rep_n(~n_a_inds), T, rep_p(~n_a_inds));

    % cut off (set to 0) anything less than eps
    b_pdf(b_pdf < eps) = 0;


end