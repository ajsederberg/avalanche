% Script to compute the branching ratio 
% for running locally, can call 'setSederbergAvaProjPaths' instead of
% setting data_dir, print_dir here. 
% data_dir is the path to the sample data (can be relative path from cwd)
data_dir = 'sample_simulation_data/';

% where you want to print outputs
print_dir = '';
data_sub_dir = 'fields1smallfinesweep/';

rep_string = 'repA';
file_name = 'sweep_smallfine_repA_stim1e-8.0et3.0ph1.0p1.0tstat';

load([data_dir data_sub_dir file_name 'fr_dynamics.mat']);
% This countains 'pop_FR,' the number of events across the population in
% each timestep of the simulation and 'fields,' the value of the latent
% variable over time. For the infinite time constant simulation, this will
% be constant for 10000 steps, and then change. 
%% For infinite time-constant simulation, this block reshapes the pop_FR 
% vector into a matrix with each column corresponding to a different field
% value. 

% Find the number of times the field changed
num_field_changes = length(find(diff([0 fields])~= 0));
% Interval between changes was constant so reshape will work
seg_fr_resps = reshape(pop_FR, length(pop_FR)/num_field_changes, num_field_changes);


%% now find a_{t+k} , a_t 
max_k = 25;
[m_k_vals, k_vals] = fitKOffsetAtAt2(seg_fr_resps, max_k);

%% 
hFig = plotMKAnalysis(m_k_vals, k_vals, seg_fr_resps);
%%
exportgraphics(hFig, [print_dir data_sub_dir rep_string '/branching_ratio_' file_name '.pdf'])

%% Apply to finite time constant time series. 

% load parameters
load([data_dir 'fields5timesweep/timeB_stim5e-12.0et4.0ph1.0p1.0t0.5fr_stats.mat'])
%%
N_F = size(J, 1);
N_neur = size(J, 2);
J_if = J'/sqrt(N_F);
eps_val = eta*epsilon;
eta_val = eta;

etaeps_tag = ['eta' num2str(eta_val, '%1.0f') 'eps' num2str(-eps_val, '%1.0f')];
%% simulate observations
num_loops = 1000;
x_max = 10000;
tau_val = 0.5*x_max;
pop_counts = zeros(x_max, num_loops);

for i_loop = 1:num_loops
T = x_max;
dt = 1;
N_T = ceil(T/dt);
tau_F = tau_val + 0*logspace(-1, 1, N_F);   % this controls the timescales of latent variables
h_F = zeros(2*N_T, N_F);
for i_f = 1:N_F
    h_F(:, i_f) = generateOUProcess(tau_F(i_f), dt, 2*T);
end
% drop initial burn-in period
h_F = h_F(T+1:end, :);
b_fac_t = exp(-eta_val*J_if*h_F' + eps_val);
prob_fire_t = b_fac_t./(1 + b_fac_t);


s_i = double(rand(N_neur, N_T) < prob_fire_t);
pop_counts(:, i_loop) = sum(s_i, 1)';
end
%%
max_k = 25;
drop_zeros = true;
[m_k_vals_tauF, k_vals2] = fitKOffsetAtAt2(pop_counts, max_k, drop_zeros);


hFig2 = plotMKAnalysis(m_k_vals_tauF, k_vals2, pop_counts);
%%
exportgraphics(hFig2, [print_dir data_sub_dir rep_string '/branching_ratio_tau' num2str(tau_val) '_5field_' etaeps_tag '.pdf'])

%%
function [m_k_vals, k_vals] = fitKOffsetAtAt2(seg_fr_resps, max_k, drop_zero)

if nargin < 3
    drop_zero = false;
end
    k_vals = 1:max_k;
    m_k_vals = 0*k_vals;
    
    for i_k = 1:length(k_vals)
    
        x_vals = seg_fr_resps(1:end-i_k, :);
        y_vals = seg_fr_resps(1+i_k:end, :);

        if drop_zero
            use_xy = x_vals~=0 & y_vals ~= 0;
            x_vals = x_vals(use_xy);
            y_vals = y_vals(use_xy);
        end
        
        p = polyfit(x_vals(:), y_vals(:), 1);
        m_k_vals(i_k) = p(1);
    end
end

function hFig = plotMKAnalysis(m_k_vals, k_vals, seg_fr_resps)

hFig = makeMyFigure(20, 10);
tiledlayout(1, 2)
nexttile
x1_vals = seg_fr_resps(1:end-1, :);
y1_vals = seg_fr_resps(2:end, :);
hc2 = histcounts2(y1_vals(:), x1_vals(:), 'BinMethod','integers', 'Normalization','probability');
% plot(x1_vals(:), y1_vals(:), '.')
imagesc(log10(hc2+eps), [-10 0]);
xlabel('a_t')
ylabel('a_{t+1}')
hc = colorbar;
% hc.Limits = [-7 0];
hold on
p1 = polyfit(x1_vals(:), y1_vals(:), 1);
u1_vals = unique(y1_vals);
plot(u1_vals, polyval(p1, u1_vals), '-', 'LineWidth',1.5);
set(gca, 'ydir', 'normal', 'fontsize', 12)
title({'color: density of points (log10 scale)'; ...
    'line: <a_{t+1}|a_{t}> = m_1 a_t + b'; ...
    ['m_1 = ' num2str(p1(1), '%1.2f')]})


nexttile
plot(k_vals, log10(m_k_vals), 'o', 'linewidth', 1)
ylim([-1 0])

p2 = polyfit(k_vals, log10(m_k_vals), 1);
hold on
plot(k_vals, polyval(p2, k_vals), '-', 'linewidth', 1.5)
title({'k-step estimate of branching ratio'; ...
    ['slope = ' num2str(p2(1), '%1.2g')];...
    ['branching ratio: ' num2str(exp(p2(1)))]})
xlabel('time offset k')
ylabel('log_{10} (m_k) (from <a_{t+k} | a_{t}>)')
set(gca, 'color', 'none', 'fontsize', 12)
legend({'m_k', 'y = k log(m_k)'}, 'location', 'southeast')

end