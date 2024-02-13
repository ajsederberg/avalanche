function [p_sil_vs_ee, p_ava_vs_ee, ave_p_silence, ave_p_ava, h_vals] = computeProbSilenceAvalanche(eta_list_ultra, eps_list_ultra, N_Neur, num_draws)
% Calculates the probability of silence for values of eta and epsilon in
% provided list, for N_Neur neurons and num_draws draws of J_i coupling
% values. Uses that to compute the probability of avalanche onset occuring.
% 

h_vals = linspace(-5, 5, 201);


% lets average over a bunch of J_i's
if nargin < 4
    num_draws = 100;
end

p_sil_vs_ee = zeros(length(eps_list_ultra),length(eta_list_ultra), length(h_vals));
for jj =1 :length(eta_list_ultra)
    this_eta_val = eta_list_ultra(jj);

    for ii = 1:length(eps_list_ultra)
        this_eps_val = eps_list_ultra(ii);
        p_sil_vs_ee_draw = zeros(length(h_vals), num_draws);

        parfor i_draw = 1:num_draws
            p_sil_vs_ee_draw(:, i_draw) = computeSilenceProb(this_eta_val, this_eps_val, randn(N_Neur, 1), h_vals);
        end
        p_sil_vs_ee(ii, jj, :) = mean(p_sil_vs_ee_draw, 2);
    end
    
%     nexttile
%     imagesc(h_vals, eps_list_ultra, squeeze(p_sil_vs_ee(:, jj, :)), [0 1])
%     title(['\eta = ' num2str(eta_list_sf1(jj))])
%     yyaxis right
%     plot(h_vals, normpdf(h_vals), 'linewidth', 1.5)
end

%% Compute probability of avalanches
p_ava_vs_ee = p_sil_vs_ee.*(1 - p_sil_vs_ee);
dh = h_vals(2) - h_vals(1);
ave_p_silence = dh*sum(p_sil_vs_ee.*repmat(reshape(normpdf(h_vals), [1 1 length(h_vals)]), [size(p_sil_vs_ee, [1 2]) 1]), 3);
ave_p_ava = dh*sum(p_sil_vs_ee.*(1 - p_sil_vs_ee).*repmat(reshape(normpdf(h_vals), [1 1 length(h_vals)]), [size(p_sil_vs_ee, [1 2]) 1]), 3);
