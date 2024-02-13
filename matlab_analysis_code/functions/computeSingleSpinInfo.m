function info_r_onespin = computeSingleSpinInfo(eta_vals, eps_vals, J_i)
% J_i is the coupling strength and is fixed for each neuron, this in
% information from single sample of single spin
%
%% Compute information for single neuron, using rates.
r_h = @(h, eta_h, eps_h) (1 + exp(-eta_h*h*J_i + eps_h)).^(-1);

info_r_onespin = zeros(length(eps_vals), length(eta_vals));
r_bar_onespin = zeros(length(eps_vals), length(eta_vals));
for ii = 1:length(eps_vals)
    for jj = 1:length(eta_vals)

        eta0 = eta_vals(jj);
        eps0 = eps_vals(ii);
        r_bar = integral(@(x) r_h(x, eta0, eps0).*exp(-x.^2/2)/sqrt(2*pi), -5, 5);
        r_bar_onespin(ii, jj) = r_bar;
        int_fun = @(h) (1 - r_h(h, eta0, eps0)).*log(eps + ((1 - r_h(h, eta0, eps0))/(1-r_bar))) + ...
            (r_h(h, eta0, eps0)).*log(eps + ((r_h(h, eta0, eps0))/(r_bar)));

        info_r = integral(@(x) normpdf(x).*int_fun(x), -5, 5);
        info_r_onespin(ii, jj) = info_r;
    end
end
end