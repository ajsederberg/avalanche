% Information calculation for model across phi, epsilon values. 
% QUESTION for size scaling: can I calculate the probability of silence for
% populations of size N? and does this help interpret the finite size
% effects of the eta-epsilon plane? 

n_neur = 1000; 

J_i = randn(n_neur, 1);
prob_s_fun = @(h, J1, eta1, eps1) 1./(1 + exp(eta1*J1*h + eps1));


%% Compute, for fixed eta, epsilon, J, h_MLE and SE_h_MLE as h_true varies. 
% fixed parameters
h_true_vals = linspace(-5, 5, 1001);
dh = h_true_vals(2) - h_true_vals(1);

h_sig = 1;
prob_h = exp(-0.5*h_true_vals.^2/h_sig^2)/sqrt(2*pi*h_sig^2);

% n_neur = 1000;

J_i = j_vals(1:10); %randn(n_neur, 1);
n_obs = 1; % think about firing rates over all and how many spikes to collect: what is reasonable value of n_obs?
%  vary these parameters
eps_vals = 0:0.5:20; %0:0.5:15;
eta_vals = 0.25:0.25:10; %1:1:100; %logspace(-1, 2, 20);

info_asymp_nats = zeros(length(eps_vals), length(eta_vals));
info_asymp2_nats = zeros(length(eps_vals), length(eta_vals));



h_MLE_arr = cell(length(eps_vals), length(eta_vals));
for i_eps = 1:length(eps_vals)
    eps_val = eps_vals(i_eps); %6;
    for i_phi = 1:length(eta_vals)
        eta_val = eta_vals(i_phi); %7;

%         h_MLE_pars = compute_hMLEpars_forsim(h_true_vals, prob_s_fun, J_i, eps_val, eta_val, n_obs);

        %% compute MI, assuming h ~ N(0, 1): < 0.5 * log(1/(sigma_h^2))>_{h}
        sigma2_hMLE = compute_h_unc(h_true_vals, J_i, eta_val, eps_val, n_obs);

%%
        % worried about behavior of sigma2_hMLE for values of h near zero
        % when the eps_val is large : if you have effectively no activity
        % because eps_val is so large, then information about the latent
        % field should go to zero. how do we set up this integral then? 
        info_hhbar = dh*prob_h*log(1./sigma2_hMLE');
        info_asymp_nats(i_eps, i_phi) = info_hhbar;
        info_hhbar2 = 2*dh*prob_h(1:2:end)*log(1./sigma2_hMLE(1:2:end)');
        info_asymp2_nats(i_eps, i_phi) = info_hhbar2;
        
                
        
        h_MLE_arr{i_eps, i_phi} = sigma2_hMLE;
        
    end
end

%%

figure()
plot(h_true_vals, r_h(h_true_vals, 4, 8))


%% compute the average rate

r_i = zeros(length(eps_vals), length(eta_vals));

for i_eps = 1:length(eps_vals)
    for i_phi = 1:length(eta_vals)
        s_prob = prob_s_fun(h_true_vals, J_i, eta_vals(i_phi), eps_vals(i_eps));
        r_i(i_eps, i_phi) = dh*prob_h*mean(s_prob, 1)';
    end
end
%% dependence on N
% 4/28/22 to-do: calculate these curves across eta-epsilon - the parts that
% are independent of N and T. 
% h_true_vals = linspace(-5, 5, 301);
h_vals_Ndep = linspace(-5, 5, 101);
% n_neur_vals = round(logspace(3, 6, 7));
n_neur_vals = round(logspace(3, 5, 5));
sigma_v_N =  zeros(length(n_neur_vals), length(h_vals_Ndep));
eta_val = 100;
eps_val = 0;
for ii = 1:length(n_neur_vals)
    sigma_v_N(ii, :) = compute_h_unc(h_vals_Ndep, randn(n_neur_vals(ii), 1), eta_val, eps_val, 1);
end
%%
makeMyFigure(15, 12);
subplot(1, 4, [1 2 3])
ph = plot(h_vals_Ndep, diag(sqrt(n_neur_vals))*(sigma_v_N), 'linewidth', 1);
assignColorsToLines(ph, parula(length(ph)))
set(gca, 'color', 'none')
xlabel('h')
ylabel('\sigma_0 (h)')
title(['\eta = ' num2str(eta_val) ', \epsilon = ' num2str(eps_val)])
legend(num2str(n_neur_vals', 'N_{neur} = %1.0g'), 'Location', 'eastoutside')
%% For small FRs, where h_{MLE} is not well determined, compute information 
% for a single neuron then scale up. Idea is that, in this range, only a
% handful of the units are contributing, and those will tend to be highly
% correlated. 
j_vals_eq = linspace(-5, 5, 21);
info_r_cell = zeros(length(eps_vals), length(eta_vals), length(j_vals_eq));
for ii = 1:length(j_vals_eq)
    info_r_cell(:, :, ii) = computeSingleSpinInfo(eta_vals, eps_vals, j_vals_eq(ii)); 
end
%% average over values of the spin
dj = j_vals_eq(2)-j_vals_eq(1);

ave_info_one_cell = dj*sum(info_r_cell.*...
    repmat(reshape(normpdf(j_vals_eq'), [1 1 21]), [size(info_r_cell, [1 2]), 1]), 3);
%%

%%
%  The information between the latent field and a single spin is bounded by
%  the entropy of the spin, and this is not dependent on T.  
%  The information between the latent field and the max-likelihood 
%  estimator increases with T and N. T and N are not fully interchangeable
%  because coupling values can vary from one cell to the other. 

%  The real question that I have is how to transition from the multi-neuron
%  estimator to the maximum likelihood estimator. 
%  

% this takes, at each eta/epsilon, the larger of the one-spin info and MLE info 
i_j_sampmax = 17;   % match to index of typical extreme value of a sample of N J's 
merge_info = max(info_r_cell(:, :, i_j_sampmax), info_asymp_nats);

figure()
tiledlayout(2, 4)
nexttile
imagesc(eta_vals, eps_vals, merge_info)


eps_inds = [5 10 15];
for ii_ei = 1:length(eps_inds)
nexttile
i_eps1 = eps_inds(ii_ei);
hold on
plot(eta_vals, merge_info(i_eps1, :), ':', 'linewidth', 2)
plot(eta_vals, info_r_cell(i_eps1, :, i_j_sampmax))
plot(eta_vals, info_asymp_nats(i_eps1, :))
title(['\epsilon = ' num2str(eps_vals(i_eps1))])
xlabel('\eta')
ylabel("information")
ylim([0 2])
legend({'merged', 'one-cell', 'MLE'})
end

J_inds = [1 8 15 21];
% for ii_ei = 1:length(eps_inds)
nexttile
% i_j_val = J_inds(ii_ei);
hold on
% plot(eta_vals, merge_info(i_eps1, :), ':', 'linewidth', 2)
plot(eta_vals, squeeze(info_r_cell(i_eps1, :, J_inds)))
% plot(eta_vals, info_asymp_nats(i_eps1, :), 'k')
% title(['J_i = ' num2str(j_vals_eq(i_j_val))])
xlabel('\eta')
ylabel("information")
ylim([0 2])
legend(num2str(j_vals_eq(J_inds)', 'J = %1.1f'))
% end
%%

j_1 = max(J_i);
j_2 = min(J_i);
info_pair = computePairSpinsInfo(eta_vals, eps_vals, j_1, j_2);
%%
info_pair2 = computePairSpinsInfo(eta_vals, eps_vals, j_vals(1), j_vals(end));
%%
tic;
[I_si_h2, I_ind_i2, I_multi2] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals(1:2));
disp(['2 neurons done in ' num2str(toc)])
tic
[I_si_h3, I_ind_i3, I_multi3] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals(1:3));
disp(['3 neurons done in ' num2str(toc)])
tic
[I_si_h4, I_ind_i4, I_multi4] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals(1:4));
disp(['4 neurons done in ' num2str(toc)])
tic
[I_si_h5, I_ind_i5, I_multi5] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals(1:5));
disp(['5 neurons done in ' num2str(toc)])

tic
[I_si_h6, I_ind_i6, I_multi6] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals(1:6));
disp(['6 neurons done in ' num2str(toc)])

tic
[I_si_h8, I_ind_i8, I_multi8] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals(1:8));
disp(['8 neurons done in ' num2str(toc)])
%%
tic
[I_si_h10, I_ind_i10, I_multi10] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals(1:10));
disp(['10 neurons done in ' num2str(toc)])


%% Are we to a regime where pair-info approaches MLE-info? 
figure()
tiledlayout(2,2)
nexttile
imagesc(info_asymp_nats, [0 2])
colorbar

nexttile
imagesc(I_si_h10, [0 2])
colorbar

nexttile
imagesc(rect(info_asymp2_nats) - I_si_h10)
colorbar

nexttile
hold on
eps_val0 = 5;
[~, eps_val_ind] = min(abs(eps_vals - eps_val0));
plot(eta_vals, info_asymp_nats(eps_val_ind, :))
plot(eta_vals, I_si_h2(eps_val_ind, :))
plot(eta_vals, I_si_h4(eps_val_ind, :))
plot(eta_vals, I_si_h6(eps_val_ind, :))
plot(eta_vals, I_si_h8(eps_val_ind, :))
plot(eta_vals, I_si_h10(eps_val_ind, :))
xlabel('\eta')
ylabel('information (nats)')
legend({'h_{MLE}, N =10; T = 1', 'two cells','four cells', 'six cells', 'eight cells', 'ten cells'}, 'Location','best')
title(['\epsilon = ' num2str(eps_val)])
%%
sigma2_hMLE = compute_h_unc(h_true_vals, J_i, 1, 10, 1);

%%
% mask_info_vals = info_vals_nats/log(2);
mask_info_vals = info_asymp_nats/log(2);
mask_info_vals(r_i < 0.01) = 0;

info_interp = @(x_eta, y_eps) interp2(eta_vals, eps_vals, mask_info_vals, x_eta, y_eps);
rate_interp = @(x_eta, y_eps) interp2(eta_vals, eps_vals, r_i, x_eta, y_eps);

%% load : eta_list, eps_list: list of eta, epsilon values for which there 
% are simulation results

rep_list = {'A', 'B', 'C', 'D', 'E'};
eta_list = [4.0, 6.0, 8.0];
eps_list = [ -6.0,  -8.0, -10.0, -12.0, -14.0, -16.0, -18.0, -20.0];
% eps_list = [ -10.0, -12.0];
[xu1_info, num_ava_u1] = loadSimulationRunResults('run_f1ultrafinesweep', rep_list);
results_dir = 'avalanches_matlab_code/analysis_results/fields1ultrafinesweep/';
x_reps = cell(length(rep_list), 1);
for ii = 1:length(x_reps)
    rep_dir = [results_dir 'rep' rep_list{ii} '/' ];
    rep_fn = dir([rep_dir 'ava_decade_analysis*.mat']);
    [~, date_ord] = sort(arrayfun(@(x) x.datenum, rep_fn), 'descend');
    x_reps{ii} = load([rep_dir rep_fn(date_ord(1)).name]);
end
%%
sum_stats_uf1 = pullSummaryFitInfo(x_reps, eta_list, eps_list);




%% this function is vry close: I am waiting on more simulation results to ill in 
% though we could use the abbreviated run data as well, for better
% coverage. 
plotInfoSpaceTrajectory(eta_vals, eps_vals, mask_info_vals, r_i, sum_stats_uf1)

%% Number of avalanches [over simulations of the same length] versus information
% Thinking is: yes, there is scaling over a large region of parameter
% space, but much of this falls where there are few avalanches. 
mean_num_ava = squeeze(mean(num_ava_u1, 1));
[eta_list_mesh, eps_list_mesh] = meshgrid(eta_list, -eps_list);
info_for_sims = info_interp(eta_list_mesh, eps_list_mesh);
rates_for_sims = rate_interp(eta_list_mesh, eps_list_mesh);
star_eps = 2;
% TIME to "collect" avalanches for quantifying power laws is 1/num_ava (all
% simulations had the same length)
makeMyFigure(27, 18);
tiledlayout(2, 3)
nexttile
plot(1./(mean_num_ava), info_for_sims, 'o-'),
hold on
plot(1./(mean_num_ava(star_eps, :)), info_for_sims(star_eps, :), 'k*'),

set(gca, 'color', 'none', 'fontsize', 12, 'box', 'off')
axis square
ylabel('information (times log(NT))'), 
xlabel('avalanche rate')
% xlabel('log_{10} N avalanches per sim')
legend(num2str(eta_list', 'eta = %1.0f'), 'location', 'southeast')
% xlim([5.3 6])

nexttile
gamma_pred = (sum_stats_uf1.alpha_values - 1)./(sum_stats_uf1.tau_values - 1);
delta_g = gamma_pred - sum_stats_uf1.gamma_fit_values;
ave_delta_g = squeeze(nanmedian(abs(delta_g), 3));
se_delta_g = squeeze(nanstd(abs(delta_g), [], 3));
hold on
plot(ave_delta_g, info_for_sims, 'o-'),
plot(ave_delta_g(star_eps, :), info_for_sims(star_eps, :), 'k*'),
set(gca, 'color', 'none', 'fontsize', 12, 'box', 'off')
axis square
ylabel('information (times log(NT))'), 
xlabel('abs (gamma - gamma_{pred})')
legend(num2str(eta_list', 'eta = %1.0f'), 'location', 'southeast')


nexttile
plot(rates_for_sims, ave_delta_g, 'o-'),
set(gca, 'color', 'none', 'fontsize', 12, 'box', 'off')
axis square
xlabel('rate (pop spikes per bin)'), 
ylabel('abs (gamma - gamma_{pred})')
% legend(num2str(eta_list', 'eta = %1.0f'), 'location', 'southeast')

nexttile
[~, inds] = intersect(eta_vals, eta_list);
plot(r_i(:, inds), mask_info_vals(:, inds), 'o-'),
hold on
[~, ind_star_eps] = min(abs(eps_vals + eps_list(star_eps)));
plot(r_i(ind_star_eps, inds), mask_info_vals(ind_star_eps, inds), 'k*'),

set(gca, 'color', 'none', 'fontsize', 12, 'box', 'off')
axis square
xlabel('rate (pop spikes per bin)'), 
ylabel('information (times log(NT))'), 


nexttile
d_eps = eps_vals(2)-eps_vals(1);
hold on
plot(edges2bins(eps_vals), diff(mask_info_vals(:, inds), [], 1)/d_eps)
plot(eps_vals([ 1 end]), [0 0], 'k')
set(gca, 'color', 'none', 'fontsize', 12, 'box', 'off')
ylabel('dI / d\epsilon')
xlabel('\epsilon')
ylim([-0.5 0.5])
xlim([0 10])
print(gcf, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/info_deltaGamma_avalanche_rate_' datestr(now, 'yyyymmdd')])



%% When does the information calculation break down? 
% prob_s_fun = @(h, J1, eta1, eps1) 1./(1 + exp(eta1*J1*h + eps1));
figure()

eta_big = 4;
test_J_i = randn(10, 1);
h_eta_big = [-3 -2 -1 0 1 2 3];
cmap = lines(length(h_eta_big));
for i_h = 1:length(h_eta_big)
s_i_vals = prob_s_fun(h_eta_big(i_h), test_J_i, eta_big, 0);

[h_MLE, lh_expr, rh_val] = find_hMLE(s_i_vals, test_J_i, eta_big, 0);
sigma_hMLE = compute_h_unc(h_MLE, test_J_i, eta_big, 0, 1);
sigma_hTrue = compute_h_unc(h_eta_big(i_h), test_J_i, eta_big, 0, 1);

%%
subplot(2,2,1)
hold on
plot(h_true_vals, lh_expr(h_true_vals), 'k')
plot(h_true_vals([1 end]), rh_val*[1 1], 'color', cmap(i_h, :), 'LineWidth', 1)
plot(h_eta_big(i_h)*[1 1], ylim, 'color', cmap(i_h, :), 'LineWidth', 1)

subplot(2, 2, 3)
hold on
plot(h_eta_big(i_h), h_MLE, 'o', 'MarkerFaceColor', cmap(i_h, :))
plot(h_eta_big, h_eta_big, 'k')
xlabel('field value (h)')
ylabel('estimated h')

subplot(2, 2, 4)
hold on
plot(h_eta_big(i_h), sigma_hMLE, 'o', 'MarkerFaceColor', cmap(i_h, :))
plot(h_eta_big(i_h), sigma_hTrue, '^', 'MarkerFaceColor', cmap(i_h, :))
% plot(h_eta_big, h_eta_big, 'k')
xlabel('field value (h)')
ylabel({'uncertainty in estimator';'(o eval at MLE, ^ at true value)'})

subplot(2, 2, 2)
hold on
histogram(s_i_vals, linspace(0, 1, 21), 'normalization', 'cdf', ...
    'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', cmap(i_h, :))
end
suptitle(['\eta = ' num2str(eta_big) ', N = ' num2str(length(test_J_i))])
%%
makeMyFigure(30,20);
sh = subplot(2,2,1);

% imagesc((1:length(phi_vals)), eps_vals, info_vals, [0 max(info_vals(:))])
imagesc(eta_vals, eps_vals, mask_info_vals, [0 max(mask_info_vals(:))])
colormap(sh, gray)
hold on
% contour((1:length(phi_vals)), eps_vals, r_i)
fr_contours = [0.05 0.1 0.2 0.3];
contour(eta_vals, eps_vals, r_i, fr_contours, 'linewidth', 1)
% contour(eta_vals, eps_vals, mask_info_vals, 0.5:.2:2.3, 'color', 'k', 'linewidth', 1.5)
colorbar
% set(gca, 'xtick', 1:2:20, 'XTickLabel', round(10*phi_vals(1:2:20))/10)
xlabel('field multiplier \eta')
ylabel('bias \epsilon')
title('I(h, h_{est})')
set(gca, 'fontsize' , 12)
hold on
% ylim([0 20])


init_eta_val = 6; % pre-MD eta value
init_eps_val = 8;
r_init = rate_interp(init_eta_val, init_eps_val);
% "recovery epsilon": find the values of epsilon that keeps constant FR
% when eta drops by 50% (MD experiment)
init_rate_vs_eps = rate_interp(init_eta_val/2, eps_vals)/r_init;
% recovery epsilon value
rec_eps_val = interp1(init_rate_vs_eps, eps_vals, 1, 'spline');

% recovery eta value
rec_eta_val = init_eta_val*2;


md_steps = 2;
rec_steps = 10;
rescale_steps = 10;
num_steps = md_steps + rec_steps + rescale_steps;
% 
% % program eta/eps trajectories step by step
% eps_trajectory = zeros(1, num_steps);
% eta_trajectory = zeros(1, num_steps);
% 
% eps_trajectory(1) = [init_eps_val];
% eta_trajectory(1) = [init_eta_val];
% eps_trajectory(2:md_steps) = rec_eps_val;
% eta_trajectory(2:md_steps) = init_eta_val/2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % program eta/eps trajectory by hand
eps_trajectory = [init_eps_val ones(1, md_steps-1)*rec_eps_val linspace(rec_eps_val, init_eps_val, rec_steps) init_eps_val*ones(1, rescale_steps)];
eta_trajectory = [init_eta_val*ones(1, md_steps) init_eta_val*ones(1, rec_steps) linspace(init_eta_val, rec_eta_val, rescale_steps) ];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


inp_scale_trajectory = [1 0.5*ones(1, md_steps - 1 + rec_steps + rescale_steps)];

effective_eta_traj = eta_trajectory.*inp_scale_trajectory;

plot(effective_eta_traj, eps_trajectory, 'wo-', 'LineWidth', 1)
subplot(4, 2, 2)
hold on
time_vec = 1:length(eps_trajectory);
plot(time_vec, info_interp(effective_eta_traj, eps_trajectory),'k', 'linewidth', 1.5)

ylabel('I(h, h_{est})')
yyaxis right
plot(time_vec, rate_interp(effective_eta_traj, eps_trajectory), 'linewidth', 1.5)
ylabel('population average firing rate')

subplot(4, 2, 4)
hold on
time_vec = 1:length(eps_trajectory);
plot(time_vec, eps_trajectory, 'k', 'linewidth', 1.5)
ylabel('\epsilon ("I-cell rate")')

yyaxis right

plot(time_vec, eta_trajectory,'linewidth', 1.5)

ylabel('"synaptic scaling" (\eta)')

% subplot(223)
% imagesc(eta_vals, eps_vals, mask_info_vals - beta*(r_i - r_init).^2)
% colorbar
% % set(gca, 'xtick', 1:2:20, 'XTickLabel', round(10*phi_vals(1:2:20))/10)
% xlabel('field multiplier \eta')
% ylabel('bias \epsilon')
% title(['I(h, h_{est}) - \beta (<r_i> - r_0)^2  (\beta = ' num2str(beta) ')'])
% set(gca, 'fontsize' , 12)
% % ylim([0 20])

subplot(223)
imagesc(eta_vals, eps_vals, r_i)
colorbar
hold on
contour(eta_vals, eps_vals, r_i, fr_contours, 'color', 0.5*[1 1 1], 'linewidth', 1)

% set(gca, 'xtick', 1:2:20, 'XTickLabel', round(10*phi_vals(1:2:20))/10)
xlabel('field multiplier \eta')
ylabel('bias \epsilon')
title('firing rate')
set(gca, 'fontsize' , 12)
% ylim([0 20])

%%
print(gcf, '-dpdf', ['write_up/MD_trajectory_ic' num2str(eps_trajectory(1)) '_' datestr(now, 'yyyymmdd')])

%% save info calcs
save(['avalanches_matlab_code/analysis_results/info_calculations_onefield_' datestr(now, 'yyyymmdd')])
%% load avalanche analysis
ava_results = load('avalanches_matlab_code/analysis_results/fields1finesweep/ava_decade_analysis_20220304.mat');
% load review_notes
load('avalanches_matlab_code/analysis_results/fields1finesweep/ava_decade_review_20220307.mat')
%%
x_var_list = ava_results.x_var_list;
y_var_list = ava_results.y_var_list;
%%

has_scaling = arrayfun(@(x) x.size_scaling | x.duration_scaling, review_notes, 'UniformOutput', false);
has_crackling = arrayfun(@(x) x.crackling, review_notes, 'UniformOutput', false);

rev_note_entries = arrayfun(@(x) ~isempty(x.crackling), review_notes);
has_crackling_flat = false(size(has_scaling));
has_crackling_flat(rev_note_entries) = arrayfun(@(x) x.crackling, review_notes(rev_note_entries));

dcc_est = 0*has_crackling_flat;
has_gamma_fit = cellfun(@(x) ~isempty(x), ava_results.all_gamma_pfit);
i_dec_start = 3;    % start fit at x_cutoff = 0.5
fit_gamma = 0*has_gamma_fit;
fit_gamma(has_gamma_fit) = cellfun(@(x) x{2}.mhat, ava_results.all_gamma_pfit(has_gamma_fit));
pred_gamma = 0*has_gamma_fit;
pred_gamma(has_gamma_fit) = cellfun(@(x) x.mhat(3, 3), ava_results.all_gamma_pred(has_gamma_fit));


sim_eta_vals = repmat(x_var_list, length(y_var_list), 1);
sim_eps_vals = repmat(-y_var_list, 1, length(x_var_list));

sim_crackles_eta_list = sim_eta_vals(has_crackling_flat);
sim_crackles_eps_list = sim_eps_vals(has_crackling_flat);
subplot(2, 2,1)
hold on
plot(sim_crackles_eta_list, sim_crackles_eps_list, 'm*')


subplot(4, 2, 6);
yyaxis left
cla;

predG_interp = interp2(x_var_list, -y_var_list, pred_gamma, eta_trajectory, eps_trajectory);
fitG_interp = interp2(x_var_list, -y_var_list, fit_gamma, eta_trajectory, eps_trajectory);
plot(predG_interp, 'LineWidth', 1)
hold on
plot(fitG_interp, 'LineWidth', 1)
ylabel('\gamma')

yyaxis right
cla;
dcc_traj = abs(fitG_interp - predG_interp);
plot(dcc_traj, 'k', 'linewidth', 1)
ylabel('DCC (|\gamma_{fit} - \gamma_{pred}|', 'Color', [0 0 0])

%%  MULTIPLE FIELDS SIMULATIONS AND MLE



%%
% % calcualte joint probability over h, h_MLE
% % each row corresponds to an h_MLE (retained in 'sort_h_MLE')
% fine_h_vals = linspace(-10, 10, 501);
% joint_p_hh = zeros(size(h_MLE_pars, 1), length(fine_h_vals));
% 
% % step through the MLE h values in order:
% [sort_h_MLE, mle_ord] = sort(h_MLE_pars(:, 1));
% 
% for i_hh = 1:length(h_true_vals)
%     h_mle = sort_h_MLE(i_hh);
%     % take the sig_h and h_true corresponding to this order:
%     sig_h = h_MLE_pars(mle_ord(i_hh), 2);
%     h_true = h_true_vals(mle_ord(i_hh));
%     
%     % put into the joint prob matrix in position i_hh (maps to h_MLE)
%     joint_p_hh(i_hh, :) = exp(-(fine_h_vals - h_mle).^2/(2*sig_h^2)); %.*exp(-fine_h_vals.^2/2)/(2*pi*sig_h);
% end
% %%
% % computing marginals (AND FAILING)
% dhMLE = diff(sort_h_MLE);
% dh = diff(fine_h_vals);
% prob_hMLE = joint_p_hh(2:end, :)'*dhMLE;
% prob_htrue = joint_p_hh(:, 2:end)*dh';
% 
% figure()
% subplot(2, 2, 1)
% imagesc(fine_h_vals, sort_h_MLE, log10(joint_p_hh), [-10 0])
% colorbar
% hold on
% plot(sort_h_MLE, sort_h_MLE, 'r.')
% xlabel('h')
% ylabel('h_{MLE}')
% 
% subplot(2, 2, 2)
% plot(sort_h_MLE, prob_htrue)
% hold on
% plot(fine_h_vals, normpdf(fine_h_vals, 0, 1))
% xlabel('h')
% ylabel('p(h)')
% 
% subplot(2, 2, 3)
% ph = plot(fine_h_vals, log(joint_p_hh(1:20:end, :)));
% ylim([-50 1])
% assignColorsToLines(ph, parula(length(ph)))
% xlabel('h')
% ylabel('p(h_{MLE}|h_{true})')
% 
% subplot(2, 2, 4)
% plot(fine_h_vals, prob_hMLE)
% hold on
% plot(fine_h_vals, normpdf(fine_h_vals, 0, 1))
% 
% xlabel('h_{MLE}')
% ylabel('p(h_{MLE})')
%%
% 
% figure()
% subplot(2, 2, 1)
% hold on
% plot(h_vals, h_vals, 'k')
% plot(h_vals', h_MLE_pars(:, 1))
% xlabel('true h')
% ylabel('h_{MLE}')
% title(['n_{neur} = ' num2str(n_neur) ', n_{obs} = ' num2str(n_obs)])
% 
% subplot(2, 2, 2)
% plot(h_vals, compute_h_unc(h_vals, J_i, eps_val, phi_val))
% xlabel('h_{MLE}')
% ylabel('SE_{hMLE}')


%% some checks on the model

%% Takes a somewhat long time to run: how does data-dependence of h_{MLE} determine h_{SE}

h_vals = linspace(-5, 5, 101);
n_neur = 10;
eps_val = 10;
eta_val = 7;
n_obs = 100;
J_i = randn(n_neur, 1);

h_MLE_pars = compute_hMLEpars_forsim(h_vals, prob_s_fun, J_i, eps_val, eta_val, n_obs);

% run the simulation 20 times with only a single obervation, compute
% h_{MLE} and SE_{h_{MLE}} and then plot. 
h_MLE_pars_short = zeros([size(h_MLE_pars) 20]);
for ii = 1:20
    h_MLE_pars_short(:, :, ii) = compute_hMLEpars_forsim(h_vals, prob_s_fun, J_i, eps_val, eta_val, 1);
end

mean_data_hMLE = mean(h_MLE_pars_short, 3);
se_data_hMLE = std(h_MLE_pars_short, [], 3);

figure()
subplot(2, 2, 1)
hold on
plot(h_vals, h_vals, 'k')
errorbar(h_vals, mean_data_hMLE(:, 1), se_data_hMLE(:, 1), 'linewidth', 1)
errorbar(h_vals, h_MLE_pars_short(:, 1, 1), h_MLE_pars_short(:, 2, 1))
ylim([-5 5])
xlabel('h_{true}')
ylabel('h_{MLE}')
legend({'true', '20-sample ave', '1-sample est'})


subplot(2, 2,3)
hold on
plot(h_vals, squeeze(h_MLE_pars_short(:, 2, :)), 'o')
plot(h_vals, h_MLE_pars(:, 2), 'k', 'linewidth', 2)
xlabel('h_{true}')
ylabel('SE_{h_{MLE}(s_i|h_{true})}')
ylim([0 10])
subplot(2, 2, 4)
plot(h_vals, compute_h_unc(h_vals, J_i, eta_val, eps_val, n_obs))
xlabel('h_{MLE}')
ylabel('SE_{h_{MLE}}')
% plot([0 10], [0 10], 'k')


%% helper functions
function logL = computeLikelihood(h_vals, s_i_bar, J_i, eta_val, eps_val)
% Computes (1/T) log L \; T is the number of samples used to compute s_i_bar 
    logL = 0*h_vals;
    
    l_fun = @(h) -eta_val*sum(J_i.*s_i_bar)*h - sum(eps_val*s_i_bar) - ...
        sum(log(1 + exp(-eta_val*J_i*h - eps_val)));
    for i_h = 1:length(h_vals)
        logL(i_h) = l_fun(h_vals(i_h));

    
    end
end

%%
function [h_MLE, lh_expr, rh_val] = find_hMLE(s_i_vals, J_i, eta_val, eps_val)
% semi-analytic calculation of MLE value of h:
% given model parameters and observation(s) s_i_vals, 
% returns maximum likelihood value of h. 
% 
% found by solving 
% \sum_i(\frac{\phi J_i}{1 + \exp(\phi J_i h + \epsilon)} = \sum_i \phi s_i J_i
% for h. 
% 
% Inputs:
%       s_i_vals: num_neurons by num_observations matrix (in interval
%           [0,1]; single observation will be binary but function also
%           works for average over observations. 
%       J_i: num_neurons by 1 coupling coefficients between latent field
%           and neural observations
%       eps_val: scalar, bias toward silence/firing 
%       phi_val: scalar, multiplies Jhs term
%
% Output: 
%       h_MLE: scalar, max-likelihood estimate of h under the single
%           latent field model 
%%
    lh_expr = @(h) sum(eta_val*J_i./(1 + exp(eta_val*J_i*h + eps_val)));
    rh_val = sum(eta_val*mean(s_i_vals, 2).*J_i);
    
    mle_fun = @(h) lh_expr(h) - rh_val;
    
    % finds where mle_fun(h) = 0, start from initial value of h = 0. 
    h_MLE = fsolve(mle_fun, 0);
        

end

function [h_MLE_unc, f_info] = compute_h_unc(h_val, J_i, eta_val, eps_val, T)
% Estimates the SE of the MLE for h at this value of epsilon, phi, and J, 
% Also returns the Fisher information (f_info = 1/sigma^2)
% Computed from the second derivative w.r.t. h of the likelihood, which is
% -1/sigma^2, the uncertainty around the max-likelihood estimate of h
% this depends only on eps_val, phi_val, and J_i
% Inputs:
%       h_val: MLE value of h. if passed as a row vector, h_MLE_unc is
%       calculated at each value of h. 
%       J_i: num_neurons by 1 coupling coefficients between latent field
%           and neural observations
%       eps_val: scalar, bias toward silence/firing 
%       phi_val: scalar, multiplies Jhs term
%
% Output: 
%       h_MLE_unc: scalar, sigma. 2nd derivative of p({s_i}|h) is -1/sigma^2 
%       f_info: scalar, if h = hMLE, then this is the Fisher information

    f_info = T*sum((eta_val*J_i*ones(size(h_val))).^2./(4*(cosh(0.5*(eta_val*J_i*h_val + eps_val)).^2)), 1);
    
    h_MLE_unc = 1./sqrt(f_info);

end

function h_MLE_pars = compute_hMLEpars_forsim(h_vals, prob_s_fun, J_i, eta_val, eps_val, n_obs)
% simulates the system H = eta_val*J_i*h*s_i + eps_val*s_i, drawing n_obs
% samples for the "true" h values provided in h_vals. From samples,
% computes the MLE parameters (estimate and SE of estimate). 
% prob_s_fun is (1 + exp(eta*Jh + eps))^{-1}. 

    n_neur = length(J_i);   % length of J_i is number of neurons
    h_MLE_pars = zeros(length(h_vals), 3);

    for i_htr = 1:length(h_vals)
        h_true = h_vals(i_htr);
        p_s_i = prob_s_fun(h_true, J_i, eta_val, eps_val);
        if ~isinf(n_obs)
            % average over observations
            s_i_obs = mean(rand(n_neur, n_obs) < repmat(p_s_i, 1, n_obs), 2);
        else
            s_i_obs = p_s_i;
        end
        h_MLE = find_hMLE(s_i_obs, J_i, eps_val, eta_val);
        h_MLE_pars(i_htr, 1) = h_MLE;
        h_MLE_pars(i_htr, 2) = compute_h_unc(h_MLE, J_i, eta_val, eps_val, n_obs);
        % compute sigma_h at h_true as well
        h_MLE_pars(i_htr, 3) = compute_h_unc(h_true, J_i, eta_val, eps_val, n_obs);


    end

end

function h_MLE_pars = bootstrap_hMLEpars_forMultiHsim(h1_vals, J_i, K_i, eta_val, eps_val, n_obs)
% estimates the MLE parameters (estimate and SE of estimate) for simulations
% with two fields (h1 and h2) and fixed J, K, eps, phi and evaluated 
% for s_i drawn n_obs times. Uses the single-field MLE estimator. Calling
% repeated for the same h1 values allows bootstrapping SE in the
% multi-field estimator. 
% NOTE: written for eta, epsilon, to be current with simulations from
% 11/15/21 onward. 
% Hamiltonian of the system is H = eta*(J_i h_1 + K_i h_2) s_i + eps*s_i
%   h1_vals: "true" values of h1, used for each "simulation" 
%   J_i: coupling to h1
%   K_i: coupling to secondary field, h2 ~ N(0, 1)
%   eps_val: bias
%   eta_val: scale for latent field coupling (same for J and K, i.e.
%   global)
%   n_obs: number of samples to draw

    % reshape coupling to be column vector
    J_i = J_i(:);
    K_i = K_i(:);

    prob_s_i_fun_h1h2 = @(h1, h2) 1./(1 + exp(eta_val*(repmat(J_i*h1, [1 length(h2)]) + K_i*h2)+eps_val));

    n_neur = length(J_i);   % length of J_i is number of neurons
    h_MLE_pars = zeros(length(h1_vals), 3);
%     h_MLE_noK = zeros(length(h1_vals), 2);
    for i_htr = 1:length(h1_vals)
        h_true = h1_vals(i_htr);
        h2_obs = randn(1, n_obs);
        p_s_i = prob_s_i_fun_h1h2(h_true, h2_obs);
        
            % average over observations
        s_i_obs = mean(rand(n_neur, n_obs) < p_s_i, 2);
%         s_i_obs'
%         h_MLE = find_h1MLE_twofields(s_i_obs, J_i, K_i, eta_val, eps_val);
        
        h_MLE = find_hMLE(s_i_obs, J_i, eta_val, eps_val);
        
        h_MLE_pars(i_htr, 1) = h_MLE;

%         compute_h_unc_latenth2(h1_val, J_i, K_i, eta_val, eps_val, s_i_ave, T)
%         h_MLE_pars(i_htr, 2) = compute_h_unc_latenth2(h_MLE, J_i, K_i, eta_val, eps_val, s_i_obs, n_obs);
        h_MLE_pars(i_htr, 2) = compute_h_unc(h_MLE, J_i, eta_val, eps_val, n_obs);
%         % compute sigma_h at h_true as well
        h_MLE_pars(i_htr, 3) = compute_h_unc(h_true, J_i, eta_val, eps_val, n_obs);
% 
%         h_MLE_noK(i_htr, 2) = compute_h_unc(h_MLE_noK(i_htr, 1), J_i, eta_val, eps_val, n_obs);

    end

end

function [h_MLE_pars, h_MLE_noK] = compute_hMLEpars_forMultiHsim(h1_vals, J_i, K_i, eta_val, eps_val, n_obs)
% computes the MLE parameters (estimate and SE of estimate) for simulations
% with two fields (h1 and h2) and fixed J, K, eps, phi and evaluated 
% for s_i drawn n_obs times. 
% NOTE: written for eta, epsilon, to be current with simulations from
% 11/15/21 onward. 
% Hamiltonian of the system is H = eta*(J_i h_1 + K_i h_2) s_i + eps*s_i
%   h1_vals: "true" values of h1, used for each "simulation" 
%   J_i: coupling to h1
%   K_i: coupling to secondary field, h2 ~ N(0, 1)
%   eps_val: bias
%   eta_val: scale for latent field coupling (same for J and K, i.e.
%   global)
%   n_obs: number of samples to draw

    % reshape coupling to be column vector
    J_i = J_i(:);
    K_i = K_i(:);

    prob_s_i_fun_h1h2 = @(h1, h2) 1./(1 + exp(eta_val*(repmat(J_i*h1, [1 length(h2)]) + K_i*h2)+eps_val));

    n_neur = length(J_i);   % length of J_i is number of neurons
    h_MLE_pars = zeros(length(h1_vals), 3);
    h_MLE_noK = zeros(length(h1_vals), 2);
    for i_htr = 1:length(h1_vals)
        h_true = h1_vals(i_htr);
        h2_obs = randn(1, n_obs);
        p_s_i = prob_s_i_fun_h1h2(h_true, h2_obs);
        
            % average over observations
        s_i_obs = mean(rand(n_neur, n_obs) < p_s_i, 2);
        s_i_obs'
        h_MLE = find_h1MLE_twofields(s_i_obs, J_i, K_i, eta_val, eps_val);
        h_MLE_pars(i_htr, 1) = h_MLE;
        
        h_MLE_noK(i_htr, 1) = find_hMLE(s_i_obs, J_i, eta_val, eps_val);
        
%         compute_h_unc_latenth2(h1_val, J_i, K_i, eta_val, eps_val, s_i_ave, T)
        h_MLE_pars(i_htr, 2) = compute_h_unc_latenth2(h_MLE, J_i, K_i, eta_val, eps_val, s_i_obs, n_obs);
%         % compute sigma_h at h_true as well
%         h_MLE_pars(i_htr, 3) = compute_h_unc(h_true, J_i, eps_val, eta_val, n_obs);
% 
        h_MLE_noK(i_htr, 2) = compute_h_unc(h_MLE_noK(i_htr, 1), J_i, eta_val, eps_val, n_obs);

    end

end

function [h_MLE, mle_lhs, mle_rhs] = find_h1MLE_twofields(s_i_vals, J_i, K_i, eta_val, eps_val)
% semi-analytic calculation of MLE value of h:
% given model parameters and observation(s) s_i_vals, 
% returns maximum likelihood value of h. 
% 
% found by solving 
% \int_{h2} P(h2) P(s_i | {h1,h2}) \sum_i(\frac{\phi J_i}{1 + \exp(\phi J_i h + \epsilon)} = \sum_i \phi s_i J_i
% for h. 
% 
% Inputs:
%       s_i_vals: num_neurons by num_observations matrix (in interval
%           [0,1]; single observation will be binary but function also
%           works for average over observations. 
%       J_i: num_neurons by 1 coupling coefficients between latent field
%           and neural observations
%       eps_val: scalar, bias toward silence/firing 
%       eta_val: scalar, multiplies (Jh1 + Kh2)*s term
%       prob_h2: P(h2), the second field. (normpdf(x)) *** note this is
%       hard-coded for now; maintain consistency with compute_hMLEpars_forMultiHsim
%           if different pdfs are used. 
%
% Output: 
%       h_MLE: scalar, max-likelihood estimate of h under the single
%           latent field model 
%%
    % computes for each cell how often s_i = 1 in samples
    s_i_ave = mean(s_i_vals, 2);
    
    % function: prob of s_i given J,K:
    prob_s_i_vh1h2 = @(h1, h2, s_i) 1./(1 + exp(eta_val*(h1*J_i.*s_i + h2*K_i.*s_i) + eps_val*s_i));

    num_mesh = 1000;
    h1_mesh = linspace(-3.5, 3.5, num_mesh);
    h2_mesh = h1_mesh;
    dh2 = h2_mesh(2)-h2_mesh(1);
    prob_h2 = @(x) normpdf(x);
    
    mle_lhs = 0*h1_mesh;
    
    mle_rhs = 0*h1_mesh;
    
    lhs_vals = zeros(1, num_mesh);
    rhs_vals = zeros(1, num_mesh);
    %%
    figure(), 
    legents = cell(10, 1);
    cmap = parula(10);
    num_mod = ceil(num_mesh/10);

    ctr = 0;
    for i_h1 = 1:num_mesh
        intgrndLHS = @(h2) sum(((1-s_i_ave).*J_i.*prob_s_i_vh1h2(h1_mesh(i_h1), h2, 0).*prob_s_i_vh1h2(h1_mesh(i_h1), h2, 1) + ...
            s_i_ave.*J_i.*prob_s_i_vh1h2(h1_mesh(i_h1), h2, 1).*prob_s_i_vh1h2(h1_mesh(i_h1), h2, 1))*prob_h2(h2));
        intgrndRHS = @(h2) sum(((1-s_i_ave).*J_i.*prob_s_i_vh1h2(h1_mesh(i_h1), h2, 0).*s_i_ave + ...
            s_i_ave.*J_i.*prob_s_i_vh1h2(h1_mesh(i_h1), h2, 1).*s_i_ave)*prob_h2(h2));
            
            for i_h2 = 1:num_mesh
                lhs_vals(i_h2) = dh2*intgrndLHS(h2_mesh(i_h2));
                rhs_vals(i_h2) = dh2*intgrndRHS(h2_mesh(i_h2));
            end
            
        if mod(i_h1, num_mod) == 1
            ctr = ctr + 1;
            subplot(2, 2, 1)
            hold on
            plot(h2_mesh, lhs_vals, 'color', cmap(ctr, :), 'linewidth', 1)
            legents{ctr} = (['h_1 = ' num2str(h1_mesh(i_h1))]);

            subplot(2, 2, 2)
            hold on
            plot( h2_mesh, rhs_vals, 'color', cmap(ctr, :), 'linewidth', 1)
            
        end
        mle_lhs(i_h1) = sum(lhs_vals); %integral(intgrndLHS, -3.5, 3.5);


        mle_rhs(i_h1) = sum(rhs_vals); %integral(intgrndRHS, -3.5, 3.5);
    end
       
    subplot(221)
    legend(legents)
    ylabel('sum_i (J_i p(s_i = 1) p(s_i|h1,h2) p(h2)')
    xlabel('h_2')
    subplot(222)
    ylabel('sum_i J_i <s>_i p(s_i | h1, h2)*p(h2)')
    xlabel('h_2')
    %%
    subplot(2, 1, 2)
    hold on
    plot(h1_mesh, mle_lhs, 'linewidth', 2)
    plot(h1_mesh, mle_rhs, 'linewidth', 2)
    xlabel('h1')
%     % evaluate the probability of 
%     lh_expr = @(h) sum(phi_val*J_i./(1 + exp(phi_val*J_i*h + eps_val)));
%     
    diff_fun = @(x) interp1(mle_lhs - mle_rhs, h1_mesh, x, 'linear');
    
%     mle_fun_vals = mle_lhs - mle_rhs;
    h_MLE = diff_fun(0);
%     rh_val = sum(s_i_ave.*J_i);
    
%     mle_fun = @(h) lh_expr(h) - rh_val;
    
    % finds where mle_fun(h) = 0, start from initial value of h = 0. 
%     h_MLE = fsolve(mle_fun, 0);
        

end

function [h_MLE_unc, f_info] = compute_h_unc_latenth2(h1_val, J_i, K_i, eta_val, eps_val, s_i_ave, T)
% Estimates the SE of the MLE for h at this value of epsilon, phi, and J, 
% Also returns the Fisher information (f_info = 1/sigma^2)
% Computed from the second derivative w.r.t. h of the likelihood, which is
% -1/sigma^2, the uncertainty around the max-likelihood estimate of h
% this depends only on eps_val, phi_val, and J_i
% Inputs:
%       h_val: MLE value of h. if passed as a row vector, h_MLE_unc is
%       calculated at each value of h. 
%       J_i: num_neurons by 1 coupling coefficients between latent field
%           and neural observations
%       eps_val: scalar, bias toward silence/firing 
%       phi_val: scalar, multiplies Jhs term
%
%
% Output: 
%       h_MLE_unc: scalar, sigma. 2nd derivative of p({s_i}|h) is -1/sigma^2 
%       f_info: scalar, if h = hMLE, then this is the Fisher information
%
% Methods:
%   Evaluates : Sum [ X''(s_i, h1) / X(s_i, h1) - (X'(s_i))/X(s_i))^2] *
%   Prod(X(s_i, h1)) 
%   where X(s_i, h1) = integral ( p(s_i | h1, h2) p(h2) dh2) is the
%   likelihood of a single data point after integrating out the second
%   field


%   need to reshape inputs (control input shape: row vectors for J, K, s
    J_i = reshape(J_i, 1, numel(J_i));
    K_i = reshape(K_i, 1, numel(K_i));
    s_i_ave = reshape(s_i_ave, 1, numel(s_i_ave));
    

    % control shape: column vector for h2_mesh
    num_mesh = 300;
    h2_mesh = linspace(-3.5, 3.5, num_mesh)';
    dh2 = h2_mesh(2)-h2_mesh(1);
    prob_h2 = @(x) normpdf(x);
    %% useful function that shows up : eta*(J_i*h1 + K_i*h_2), 
    %     has dimensions of h2 by K_i
    etaKH_sq = @(h2) eta_val*(bsxfun(@plus, h1_val*J_i, h2*K_i));
    %% function for (1 + exp(eta(Jh1 + Kh2) + eps)^(-1)
    prob_s_i1_vh1h2 = @(h2) 1./(1 + exp(etaKH_sq(h2) + eps_val));

    %% p(h2) p(s | h1,h2): has dimensions of h2 by s_i
    %  one term for each in time average, <s> P(s = 1 |h1,h2) + (1 - <s>) P(s = 0|h1,h2)
%     int_prefac_s0 = @(h2) diag(prob_h2(h2))*(prob_s_i1_vh1h2(h2)*diag(s_i_ave));
%     int_prefac_s1 = @(h2) diag(prob_h2(h2))*((1 - prob_s_i1_vh1h2(h2))*diag(1-s_i_ave));
%     
% TOTAL CLUSTER BELOW
    %% First term X (expanded in h2) for s_i = 0 and s_i = 1
    X_h2_s_i_1 = @(h2) diag(prob_h2(h2))*(prob_s_i1_vh1h2(h2));
    X_h2_s_i_0 = @(h2) diag(prob_h2(h2))*((1 - prob_s_i1_vh1h2(h2)));
    
    %% Second term X' (expanded in h2) for s_i = 0 and s_i = 1
    Xp_h2_s_i_1 = @(h2) (X_h2_s_i_1(h2)*diag(eta_val*J_i)).*(prob_s_i1_vh1h2(h2) - 1);
    Xp_h2_s_i_0 = @(h2) (X_h2_s_i_0(h2)*diag(eta_val*J_i)).*(prob_s_i1_vh1h2(h2));
    
    %% second derivative terms: etaJ^2*(1/1+exp - s_i)^2 - (etaJ)^2/4cosh^2(JK/2)
    %% term 1: (eta J)^2 / (4 cosh^2((etaJK + eps )/ 2))
    term1_fn = @(h2) (1./(4*cosh((etaKH_sq(h2) + eps_val)/2).^2))*diag(eta_val^2*J_i.^2);
    
%     %% term 2: s_i_ave*(eta J)^2 / (1 + exp(etaJK + eps))
    
%     term2_fn = @(h2) (1./(1 + exp(etaKH_sq(h2) + eps_val)))*diag(eta_val^2*J_i.^2.*s_i_ave);
    term2_fn0 = @(h2) (prob_s_i1_vh1h2(h2).^2)*diag(eta_val^2*J_i.^2);
    term2_fn1 = @(h2) ((prob_s_i1_vh1h2(h2) - 1).^2)*diag(eta_val^2*J_i.^2);

    Xpp_h2_s_i_1 = @(h2) Xp_h2_s_i_1(h2).*(-term1_fn(h2) + term2_fn1(h2));
    Xpp_h2_s_i_0 = @(h2) Xp_h2_s_i_0(h2).*(-term1_fn(h2) + term2_fn0(h2));
    
    %% 
    Xpp_over_X = sum(T*s_i_ave.*(sum(Xpp_h2_s_i_1(h2_mesh), 1)./sum(X_h2_s_i_1(h2_mesh), 1)) + ... 
        T*(1 - s_i_ave).*(sum(Xpp_h2_s_i_0(h2_mesh), 1)./sum(X_h2_s_i_0(h2_mesh), 1)));
    
    Xp_over_X = sum(T*s_i_ave.*(sum(Xp_h2_s_i_1(h2_mesh), 1).^2./sum(X_h2_s_i_1(h2_mesh), 1).^2) + ...
        T*(1 - s_i_ave).*(sum(Xp_h2_s_i_0(h2_mesh), 1).^2./sum(X_h2_s_i_0(h2_mesh), 1).^2));
    
    %% total cluster
    logprodX = sum(T*(1 -s_i_ave).*log(sum(Xp_h2_s_i_0(h2_mesh),1 )) + ... 
        T*s_i_ave.*log(sum(Xp_h2_s_i_1(h2_mesh), 1)));
    %%
%     %% term 3 : (eta J)^2 / (1 + exp(etaJK + eps))^2
%     term3_fn = @(h2) ((1./(1 + exp(etaKH_sq(h2) + eps_val))).^2)*diag(eta_val^2*J_i.^2);

%     noH2_sig_v2 = int_prefac(h2_mesh).*term1_fn(h2_mesh);
    % term1 + term2 - term3
    sig_v_h2 = int_prefac(h2_mesh).*(term1_fn(h2_mesh) + term2_fn(h2_mesh) - term3_fn(h2_mesh));
% %     f_info = T*sum((eta_val*J_i*ones(size(h_val))).^2./(4*(cosh(0.5*(eta_val*J_i*h_val + eps_val)).^2)), 1);
    f_info = T*sum(sum(sig_v_h2));
    h_MLE_unc = 1./sqrt(f_info);

end

function hFig = plotInfoSpaceTrajectory(eta_vals, eps_vals, mask_info_vals, r_i, sim_delta_gamma)

hFig = makeMyFigure(30,20);
sh = subplot(2,2,1);

% imagesc((1:length(phi_vals)), eps_vals, info_vals, [0 max(info_vals(:))])
imagesc(eta_vals, eps_vals, mask_info_vals, [0 0.9*max(mask_info_vals(:))])
colormap(sh, gray)
hold on
% contour((1:length(phi_vals)), eps_vals, r_i)
fr_contours = [0.05 0.1 0.2 0.3];
contour(eta_vals, eps_vals, r_i, fr_contours, 'linewidth', 1)
% contour(eta_vals, eps_vals, mask_info_vals, 0.5:.2:2.3, 'color', 'k', 'linewidth', 1.5)
colorbar
% set(gca, 'xtick', 1:2:20, 'XTickLabel', round(10*phi_vals(1:2:20))/10)
xlabel('field multiplier \eta')
ylabel('bias \epsilon')
title('I(h, h_{est})')
set(gca, 'fontsize' , 12)
hold on
% ylim([0 20])


info_interp = @(x_eta, y_eps) interp2(eta_vals, eps_vals, mask_info_vals, x_eta, y_eps);
rate_interp = @(x_eta, y_eps) interp2(eta_vals, eps_vals, r_i, x_eta, y_eps);

init_eta_val = 6; % pre-MD eta value
init_eps_val = 8;
r_init = rate_interp(init_eta_val, init_eps_val);
% "recovery epsilon": find the values of epsilon that keeps constant FR
% when eta drops by 50% (MD experiment)
init_rate_vs_eps = rate_interp(init_eta_val/2, eps_vals)/r_init;
% recovery epsilon value
rec_eps_val = interp1(init_rate_vs_eps, eps_vals, 1, 'spline');

% recovery eta value
rec_eta_val = init_eta_val*2;

%%
seg_steps = 20;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % program eta/eps trajectory by hand
eps_trajectory = [init_eps_val ones(1, seg_steps-1)*rec_eps_val;  linspace(rec_eps_val, init_eps_val, seg_steps);  init_eps_val*ones(1, seg_steps)];
eta_trajectory = [init_eta_val*ones(1, seg_steps) ; init_eta_val*ones(1, seg_steps); linspace(init_eta_val, rec_eta_val, seg_steps) ];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_colors = hsv(size(eps_trajectory, 1));

inp_scale_trajectory = 0.5*ones(size(eps_trajectory));
inp_scale_trajectory(1,1) = 1;
time_vec = 1:3*seg_steps;

for ii = 1:size(eps_trajectory, 1)
    effective_eta_traj = eta_trajectory(ii, :).*inp_scale_trajectory(ii, :);
    eff_eps_traj = eps_trajectory(ii, :);
    t_vec = (ii-1)*seg_steps + (1:seg_steps);
    % plot on the parameter diagram
    subplot(2, 2, 1)
    hold on
    plot(effective_eta_traj, eff_eps_traj, '-', 'color', step_colors(ii, :), 'LineWidth', 1)
    
    % plot rates etc vs time
    subplot(4, 2, 2)
    hold on
    plot(t_vec, info_interp(effective_eta_traj, eff_eps_traj),'-', 'color', step_colors(ii, :), 'linewidth', 1.5)

    ylabel('I(h, h_{est})')
    
    subplot(4, 2, 4)
    hold on
    plot(t_vec, rate_interp(effective_eta_traj, eff_eps_traj),'-', 'color', step_colors(ii, :), 'linewidth', 1.5)
    ylabel('population average firing rate')

    subplot(4, 2, 6)
    hold on
    plot(t_vec, eff_eps_traj, '-', 'color', step_colors(ii, :), 'linewidth', 1.5)
    ylabel('\epsilon ("threshold")')

    subplot(4, 2, 8)
    hold on
    plot(t_vec, effective_eta_traj, '-', 'color', step_colors(ii, :), 'linewidth', 1.5)
    ylabel('"synaptic scaling" (\eta)')

end
%%
gamma_pred = (sim_delta_gamma.alpha_values  - 1)./(sim_delta_gamma.tau_values -1);
gamma_fit = sim_delta_gamma.gamma_fit_values;
delta_gamma = squeeze(nanmean(abs(gamma_pred - gamma_fit), 3));
eta_sim = sim_delta_gamma.eta_list;
eps_sim = sim_delta_gamma.eps_list;
subplot(2, 2, 3)
imagesc(eta_sim, -eps_sim, (delta_gamma), [0 0.2])
xlabel('field multiplier \eta')
ylabel('bias \epsilon')
colorbar

xlim(eta_vals([1 end]))
ylim(eps_vals([1 end]))
end


function info_r_onespin = computeSingleSpinInfo(eta_vals, eps_vals, J_i)
% J_i is the coupling strength and is fixed for each neuron
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

function [info_pair, I_1, I_2, I_12] = computePairSpinsInfo(eta_vals, eps_vals, j_1, j_2)
% computes the infromation that the pair of spins {s1, s2} has about the
% latent field

I_1 = computeSingleSpinInfo(eta_vals, eps_vals, j_1);
I_2 = computeSingleSpinInfo(eta_vals, eps_vals, j_2);

% compute information s1 and s2 have about each other
r_h = @(h, eta_h, eps_h, J_i) (1 + exp(-eta_h*h*J_i + eps_h)).^(-1);

I_12 = 0*I_1;

for ii = 1:length(eps_vals)
    for jj = 1:length(eta_vals)

        eta0 = eta_vals(jj);
        eps0 = eps_vals(ii);

        r_bar1 = integral(@(x) r_h(x, eta0, eps0, j_1).*exp(-x.^2/2)/sqrt(2*pi), -5, 5);
        r_bar2 = integral(@(x) r_h(x, eta0, eps0, j_2).*exp(-x.^2/2)/sqrt(2*pi), -5, 5);

        % compute P(s1, s2) : integral (over h) of p(h) p(s1 | h) p(s2 | h)
        P_s1s2 = zeros(2);
        % 0,0
        r00_fun = @(h) (1 - r_h(h, eta0, eps0, j_1)).*(1 - r_h(h, eta0, eps0, j_2));
        P_s1s2(1,1) = integral(@(x) normpdf(x).*r00_fun(x), -5, 5);

        % 0,1
        r01_fun = @(h) (1 - r_h(h, eta0, eps0, j_1)).*(r_h(h, eta0, eps0, j_2));
        P_s1s2(1,2) = integral(@(x) normpdf(x).*r01_fun(x), -5, 5);

        % 1,0
        r10_fun = @(h) (r_h(h, eta0, eps0, j_1)).*(1 - r_h(h, eta0, eps0, j_2));
        P_s1s2(2,1) = integral(@(x) normpdf(x).*r10_fun(x), -5, 5);

        % 1,1
        r11_fun = @(h) (r_h(h, eta0, eps0, j_1)).*(r_h(h, eta0, eps0, j_2));
        P_s1s2(2,2) = integral(@(x) normpdf(x).*r11_fun(x), -5, 5);

        %
        % compute information 
        P_ind_s1s2 = [(1 - r_bar1)*(1 - r_bar2) (1-r_bar1)*(r_bar2); ...
            (r_bar1)*(1-r_bar2) (r_bar1)*(r_bar2)];
       
        % 
%         r_bar1_onespin(ii, jj) = r_bar1;
%         r_bar2_onespin(ii, jj) = r_bar2;
        
%         info_r = integral(@(x) normpdf(x).*info_fun(x), -5, 5);
        I_12(ii, jj) = sum(P_s1s2.*log(P_s1s2./P_ind_s1s2), 'all');
    end
end

info_pair = I_1 + I_2 - I_12;




end

function [I_si_h, I_ind_i, I_multi] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals)
% computes the information that the set of spins {couplings j_vals} has about the
% latent field
% warning: number of terms computed goes at 2^length(j_vals)

if length(j_vals) > 15
    I_si_h = nan;
    I_ind_i = nan;
    I_multi = nan;
    disp('will not attempt for length(j_vals) > 8')
    return
end
I_ind_i = zeros(length(eps_vals), length(eta_vals), length(j_vals));
for ii = 1:length(j_vals)
    I_ind_i(:, :, ii) = computeSingleSpinInfo(eta_vals, eps_vals, j_vals(ii));
end
I_ind = sum(I_ind_i, 3);
% compute multi information among spins with coupling values j_vals

% this gives the rate function for a spin with coupling J_i as a function
% of field h
r_h = @(h, eta_h, eps_h, J_i) (1 + exp(-eta_h*h*J_i + eps_h)).^(-1);

I_multi = 0*I_ind;

%% construct binary state vecs
num_states = 2^length(j_vals);
s_i_include = zeros(num_states, length(j_vals));
for ii = 1:length(j_vals)
    s_i_include(:, ii) = mod(floor(((1:num_states)-1)/2^(ii-1)), 2);
end
%%
for ii = 1:length(eps_vals)
    for jj = 1:length(eta_vals)

        eta0 = eta_vals(jj);
        eps0 = eps_vals(ii);

        % compute <s_i> for each j_i
        r_bar_i = 0*j_vals;
        for kk = 1:length(j_vals)
            r_bar_i(kk) = integral(@(x) r_h(x, eta0, eps0, j_vals(kk)).*exp(-x.^2/2)/sqrt(2*pi), -5, 5);
        end

        % compute joint distribution : integral (over h) of p(h) p(s1 | h) p(s2 | h)
        P_si_joint = zeros(num_states, 1);
        % 0,0
        P_ind_s_i = ones(num_states, 1);
        for i_state = 1:num_states
            %% compose the function P(s_i = 0,1 | h) for all 0,1 combinations
            r_fun = @(h) 1;
            for kk = 1:length(j_vals)
                % construct joint probabiltiy term
                new_r_fun = @(h) (1-s_i_include(i_state, kk)) + ...
                    (2*s_i_include(i_state, kk) - 1)*r_h(h, eta0, eps0, j_vals(kk));
                r_fun = @(h) new_r_fun(h).*r_fun(h);

                % calculate independent probability term
                P_ind_s_i(i_state) = P_ind_s_i(i_state)*((1-s_i_include(i_state, kk)) + ...
                    (2*s_i_include(i_state, kk) - 1)*r_bar_i(kk));
            end
%%
            %% evaluate integral over h
            P_si_joint(i_state) = integral(@(x) normpdf(x).*r_fun(x), -5, 5);

        end
        
        %
        % compute information 
%         P_ind_s1s2 = [(1 - r_bar1)*(1 - r_bar2) (1-r_bar1)*(r_bar2); ...
%             (r_bar1)*(1-r_bar2) (r_bar1)*(r_bar2)];
       
        % 
%         r_bar1_onespin(ii, jj) = r_bar1;
%         r_bar2_onespin(ii, jj) = r_bar2;
        
%         info_r = integral(@(x) normpdf(x).*info_fun(x), -5, 5);
        I_multi(ii, jj) = sum(P_si_joint.*log(eps + P_si_joint./P_ind_s_i));
    end
end

I_si_h = I_ind - I_multi;




end