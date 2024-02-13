% Sub-sampling script: analysis of avalanche size distributions with fewer
% samples for a simulation with scaling at high sample count

% example parameters: N_F = 1, eps = -18, eta = 4
param_str = 'e-6.0et6.0';
uf_files = dir(['avalanches_code/data/fields1ultrafinesweep/sweep*' param_str '*tstat.mat']);
uffr_files = dir(['avalanches_code/data/fields1ultrafinesweep/sweep*' param_str '*fr_stats.mat']);

f_file = ['avalanches_code/data/fields1finesweep/envfixedJh_sweep_fine_av_stim1' param_str 'ph1.0p1.0tstat.mat'];
%% check out the time simulations
param_str = 't0.5';
uf_files = dir(['avalanches_code/data/fields1timesweep/time*' param_str '.mat']);
uffr_files = dir(['avalanches_code/data/fields1timesweep/time*' param_str 'fr_stat.mat']);




%%
uf_ava_array = cell(length(uf_files), 1);

for i_uf = 1:length(uf_files)
    uf_ava_array{i_uf} = load([uf_files(i_uf).folder filesep uf_files(i_uf).name]);
end

%% plot size pdfs
makeMyFigure(50, 10);
n_decimate = 5;
dec_map = hot(n_decimate*2);
dec_map = flipud(dec_map(1:n_decimate, :));

tiledlayout(1, 5)

for i_rep = 1:length(uf_ava_array)
    all_sizes = uf_ava_array{i_rep}.sizes;
    L_tot = length(all_sizes);
    L_dec = round(logspace(-1.5, 0, n_decimate)*L_tot);
    sh =nexttile;
    eh = zeros(n_decimate, 1);
    for i_dec = 1:n_decimate
        
        eh(i_dec) = plotLogscalePDF(sh, all_sizes(randperm(L_tot, L_dec(i_dec))));

%         eh(i_dec) = plotLogscalePDF(sh, all_sizes(1:L_dec(i_dec)));
    end
    assignColorsToLines(eh, dec_map);
    set(gca, 'color', 'none')
    xlabel('size (binned)')
    ylabel('pdf')
    legend(num2str(L_dec', 'N_{ava} = %1.0f'))
end
suptitle(param_str)

%%

for ii = 1:length(rep_list)
    x_info{ii} = load([uffr_files(ii).folder filesep uffr_files(ii).name]);
end

%%
J_mats = cellfun(@(x) x.J, x_info, 'UniformOutput', false);
[~, ord1] = sort(J_mats{1}(1, :));
figure()
for ii = 1:length(J_mats)
    subplot(length(J_mats), 2, 2*ii - 1)
    imagesc(J_mats{ii}(:, ord1))
    
    subplot(length(J_mats), 2, 2*ii )
    histogram(J_mats{ii})
end

%% Firing rate distribution
eff_bin_size = 0.003;    % bin size in seconds
fr_mats = cellfun(@(x) x.cell_FR/eff_bin_size, x_info, 'UniformOutput',false);
figure()
tiledlayout(3, 2)
for ii = 1:length(fr_mats)
    nexttile
    histogram(fr_mats{ii})
end
%% 
function eh = plotLogscalePDF(h_sub, x)
% plot pdf(x) on log scale
d_msz = 8;
data_color = [0 0 0];
log_size_bins = unique(round(logspace(0, 5, 51)));  % use unique(round()) to avoid having 1, 1.25, 1.58, etc.
size_pdf = histcounts(x, bins2edges(log_size_bins), 'Normalization', 'pdf');
size_cts = histcounts(x, bins2edges(log_size_bins), 'Normalization', 'count');
rel_err = 1./sqrt(size_cts-1);

size_pdf_SE = rel_err.*size_pdf; % get counting statistics error bars on PDF estimates
set(gcf, 'currentaxes', h_sub)
hold on
% ERRORBAR PLOT
eh = errorbar(log10(log_size_bins), log10(size_pdf), ...
    log10(size_pdf) - log10(size_pdf-size_pdf_SE), ...  % lower errorbar
    log10(size_pdf + size_pdf_SE) - log10(size_pdf), '.', ...  % upper errorbar
    'markersize', d_msz, 'capsize', 0, 'color', data_color)
end