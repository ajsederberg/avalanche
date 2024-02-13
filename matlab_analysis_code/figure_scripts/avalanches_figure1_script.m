% Visualization of spiking activity in our model 

%% generate latent process, random coupled spins, covariance matrix

N_F = 5;
N_neur = 128;
J_if = randn(N_neur, N_F)/sqrt(N_F);
% 
T = 150;
dt = 0.003; % in seconds
N_T = ceil(T/dt); % 
tau_F = 15 + 0*logspace(-1, 1, N_F); % in seconds
h_F = zeros(N_T, N_F);
for i_f = 1:N_F
    h_F(:, i_f) = generateOUProcess(tau_F(i_f), dt, T);
end
%%

%%
eps_val = -10 ; %2.67;
eta_val = 8;

b_fac_t = exp(-eta_val*J_if*h_F' + eps_val);
prob_fire_t = b_fac_t./(1 + b_fac_t);


s_i = rand(N_neur, N_T) < prob_fire_t;

%% embed in 2 dimensions
y_pos = tsne(prob_fire_t);
% blur a little
% y_pos = y_pos + 0.5*(rand(size(y_pos))-.5);
%% and cluster
idx = kmeans(prob_fire_t, 50);
[~, k_ord] = sort(idx);

%% find an avalanches
num_active = sum(s_i, 1);

size_thr = prctile(num_active, 80);
% find an avlanche with at least the 80th %ile active; rounding and taking mode finds
% a long one
% bigS = mode(round(find(num_active > size_thr), -2));
% or pick one
bigS = 16333;
% find the last 0 before it
has0 = find(num_active == 0);
last0 = sum(has0 <= bigS);
start_point = has0(last0);
end_point = has0(last0+1);

%%
makeMyFigure(33, 10);
nexttile
hold on
imagesc((1:N_T)*dt, [], s_i)
colormap(flipud(gray))
plot(dt*[start_point end_point end_point start_point start_point], N_neur*[0 0 1 1 0], 'g-', 'linewidth', 3)
axis tight
xlabel('time (s)')
ylabel('neurons')
title(['N_F = ' num2str(N_F) ', \tau_F = ' num2str(tau_F(1))])
set(gca, 'fontsize', 16)
print(gcf, '-dsvg', 'avalanches_matlab_code/plots/paper_figures/poster_figures/avalanche_raster_example.svg')
%%
makeMyFigure(20, 10)

for ii = 1:N_F
    subplot(N_F, 1, ii)
    plot(h_F(:, ii), 'k')
    set(gca,'color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', ...
        'YColor', 'none', 'fontsize', 14)
    xlim([N_T - 15000 N_T] )
end

print(gcf, '-dsvg', 'avalanches_matlab_code/plots/paper_figures/poster_figures/latent_fields_example.svg')

%%
makeMyFigure(33, 5);
% time_steps = [7153 2150 2300 2450 2520]
time_steps = round(linspace(start_point, end_point, 7));
for i_t = 1:length(time_steps)
    nexttile
    hold on
    on_cells = s_i(:, time_steps(i_t));
    plot(y_pos(~on_cells, 1), y_pos(~on_cells, 2), 'o', 'linewidth', 1, ...
        'Color',0.5*[1 1 1]);
    
    plot(y_pos(on_cells, 1), y_pos(on_cells, 2), 'ko', 'linewidth', 1, ...
        'MarkerFaceColor',[0 1 0 ]);
    set(gca,'color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', ...
        'YColor', 'none', 'fontsize', 14)
    axis square
    title(['t = ' num2str(time_steps(i_t)*dt, '%1.1f')], 'FontWeight','normal')
end
print(gcf, '-dsvg', 'avalanches_matlab_code/plots/paper_figures/poster_figures/avalanche_progression.svg')

%%
function y = generateOUProcess(tau, dt, T)

    N = ceil(T/dt);
    y = zeros(N, 1);
    
    for i_t = 2:N
        y(i_t) = (1 - dt/tau)*y(i_t-1) + sqrt(2*dt/tau)*randn;
    end

end