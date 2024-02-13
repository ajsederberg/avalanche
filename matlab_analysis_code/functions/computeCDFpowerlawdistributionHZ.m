function cdf_Paxm = computeCDFpowerlawdistributionHZ(x_vals, a_hat, x_min)
% computes the complementary CDF of a discrete power law distribution using
% the generalized zeta function (Hurwitz zeta)
%     cdf_Paxm = hurwitzZeta(a_hat, x_vals)/hurwitzZeta(a_hat, x_min);

%     [glu_hurwitzZeta, ~] = loadhurwitzZetaLookup();

%     global glu_hurwitzZeta
%     cdf_Paxm = glu_hurwitzZeta(a_hat, x_vals)/glu_hurwitzZeta(a_hat, x_min);
cdf_Paxm = 0*x_vals;

cdf_Paxm(1) = hurwitzZeta(a_hat, x_vals(1))/hurwitzZeta(a_hat, x_min);
if length(x_vals) > 1
    %     if any(x_vals > 1e4)
    % evaluate at 1e3, well within interpolation bounds
    % this finds constants for log(ccdf) = a0 - (a_hat-1)*log(x)
    %
    a0 = (a_hat - 1)*log(x_vals(1)) + log(cdf_Paxm(1));

%     logxH = log(x_vals(2:end));

%     cdf_Paxm(2:end) = exp(a0 - (a_hat - 1)*log(x_vals(2:end)));
    cdf_Paxm = exp(a0 - (a_hat - 1)*log(x_vals));
    %         warning('Note: Hurwitz Zeta evaluated outside of bounds (x > 1e4).')

end
end