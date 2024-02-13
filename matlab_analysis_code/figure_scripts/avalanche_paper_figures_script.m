%% Requires these files, which can be found in 'avalanches_matlab_code/':
%       % analysis_results/fields1finesweep/ava_decade_analysis*
%       % analysis_results/fields5finesweep/ava_decade_analysis*
%       % analysis_results/fields20finesweep/ava_decade_analysis*
%       % analysis_results/timesweep_stim1_e12e04/ava_decade_analysis*
%       % analysis_results/timesweep_stim5_e12e04/ava_decade_analysis*
%       % analysis_results/timesweep_stim20_e12e04/ava_decade_analysis*

% And in avalanches_code/data:
%       fr_stats files, tstat files

%% Script to plot FR statistics across simulations

ss_str.run_f20_e3010_timesweep = false;       % pull data from 20-field, time sweep
ss_str.run_f5_e1204_timesweep = false;       % pull data from 5-field, eps = -12, eta = 4 time sweep
ss_str.run_f1_e1204_timesweep = false;       % pull data from 1-field, eps = -12, eta = 4 time sweep

ss_str.run_f1finesweep = false;       % pull data from 1-field, fine sweep of eta, eps
ss_str.run_f5finesweep = false;       % pull data from 5-field, fine sweep
ss_str.run_f20finesweep = true;       % pull data from 20-field, fine sweep





%% Get the firing rate distribution, coupling parameters, etc for each sim. 
sim_set = {'run_f1finesweep', 'run_f5finesweep', 'run_f20finesweep'};
sim_set_colors = lines(length(sim_set));
fr_out = cell(length(sim_set), 1);
for kk = 1:length(sim_set)
    
    ss_str = switchSimulationSet(ss_str, sim_set{kk});
    
    % set only one of the above to true/false, then run this script
    setLoadingFunctionNames
    
    clear fr_stats_array
    fr_stats_array(length(y_var_list), length(x_var_list)) = struct('J', [], ...
        'epsilon', [], 'eta', [], 'cell_FR', [], 'num_field_samples', []);

    for ii = 1:length(y_var_list)
        for jj = 1:length(x_var_list)
            try
                fr_stats_array(ii, jj) = load(fr_name_fn(ii, jj), ...
                    'J', 'epsilon', 'eta', 'cell_FR', 'num_field_samples');
            catch
                disp(['failed for i,j ' num2str([ii jj], '%1.0f , %1.0f')])
            end
        end
    end

    has_fr_stat = arrayfun(@(x) ~isempty(x.J), fr_stats_array);
    fr_out{kk} = fr_stats_array;
end
%%
nR = length(y_var_list);
nC = length(x_var_list);
h_fig1 = makeMyFigure(3*nC, 4*nR);


%%
figure(h_fig1);
for kk = 1:length(sim_set)
    fr_stats_array = fr_out{kk};
    for ii = 1:nR
        for jj = 1:nC
    %         if has_fr_stat(ii,jj)
%             try
                subplot(nR, nC, (ii-1)*nC + jj)
                hold on
    %             histogram(fr_stats_array(ii,jj).cell_FR)
                logRank = log10(1:length(fr_stats_array(ii,jj).cell_FR));

                plot(logRank, sort(fr_stats_array(ii,jj).cell_FR, 'descend'))
    %             title(['\eta = ' num2str(fr_stats_array(ii,jj).eta) ...
    %                 ', \epsilon = ' num2str(fr_stats_array(ii,jj).epsilon)])
%             end
        end
    end
end
%%
figure(h_fig1);
for ii = 1:nR
    for jj = 1:nC
        subplot(nR, nC, (ii-1)*nC + jj)
        set(gca, 'color', 'none')
        ylim([0. 0.4])
        xlim([0 3])
        if ii == nR
            xlabel([x_lab_sum ' = ' num2str(x_var_list(jj))])
        else
            set(gca, 'xtick', [])
        end
        if jj == 1
            ylabel([y_lab_sum ' = ' num2str(y_var_list(ii))])
        else
            set(gca, 'ytick', [])
        end
    end
end
%%

%% Plot avalanche analyses : where is there scaling? 

avalanche_analysis_files = {'fields1finesweep/ava_decade_analysis_20220217.mat', ...
    'fields5finesweep/ava_decade_analysis_20220215.mat', ...
    'fields20finesweep/ava_decade_analysis_20220216.mat'};

ava_out = cell(length(avalanche_analysis_files), 1);
for ii =1 :length(avalanche_analysis_files)
    
    ava_out{ii} = load(['avalanches_matlab_code/analysis_results/' avalanche_analysis_files{ii}], ...
        'all_alpha_pfit', 'all_gamma_pfit', 'all_tau_pfit');
end

%% find where size durations (tau) scale

tau_scaling = nan([size(fr_stats_array) length(ava_out)]);
alpha_scaling = nan([size(fr_stats_array) length(ava_out)]);
gamma_scaling = nan([size(fr_stats_array) length(ava_out)]);

tau_values = nan([size(fr_stats_array) length(ava_out)]);
alpha_values = nan([size(fr_stats_array) length(ava_out)]);
gamma_values = nan([size(fr_stats_array) length(ava_out)]);
tau_se_values = nan([size(fr_stats_array) length(ava_out)]);
alpha_se_values = nan([size(fr_stats_array) length(ava_out)]);
gamma_se_values = nan([size(fr_stats_array) length(ava_out)]);


% number of SD from mean of surrogate data for KS stat
tau_KS_nSD = nan([size(fr_stats_array) length(ava_out)]);
alpha_KS_nSD = nan([size(fr_stats_array) length(ava_out)]);


for ii = 1:length(ava_out)
    % find the percentile of the actual KS statistic relative to the
    % surrogate data set 
    has_afits = cellfun(@(x) ~isempty(x), ava_out{ii}.all_alpha_pfit);
    has_tfits = cellfun(@(x) ~isempty(x), ava_out{ii}.all_tau_pfit);
    
    apv = nan(size(has_afits));
    apv(has_afits) = cellfun(@(x) mean(x.min_KS_stat > x.ks_surrogate), ... 
        ava_out{ii}.all_alpha_pfit(has_afits));
    alpha_scaling(:, :, ii) = apv;
    avs = nan(size(has_tfits));
    avs(has_afits) = cellfun(@(x) x.a_hat, ava_out{ii}.all_alpha_pfit(has_afits));
    alpha_values(:, :, ii) = avs;
    avs(has_afits) = cellfun(@(x) x.se_a_hat, ava_out{ii}.all_alpha_pfit(has_afits));
    alpha_se_values(:, :, ii) = avs;
    
    % find the number of SD from the mean of the surrogate statistic for KS
    % stat
    apv(has_afits) = cellfun(@(x) (x.min_KS_stat - mean(x.ks_surrogate))/std(x.ks_surrogate), ... 
        ava_out{ii}.all_alpha_pfit(has_afits));
    alpha_KS_nSD(:, :, ii) = apv;
    
    
    tpv = nan(size(has_tfits));
    tpv(has_tfits) = cellfun(@(x) mean(x.min_KS_stat > x.ks_surrogate), ... 
        ava_out{ii}.all_tau_pfit(has_tfits));
    tau_scaling(:, :, ii) = tpv;
    tvs = nan(size(has_tfits));
    tvs(has_tfits) = cellfun(@(x) x.a_hat, ava_out{ii}.all_tau_pfit(has_tfits));
    tau_values(:, :, ii) = tvs;
    tvs(has_tfits) = cellfun(@(x) x.se_a_hat, ava_out{ii}.all_tau_pfit(has_tfits));
    tau_se_values(:, :, ii) = tvs;

    % find the number of SD from the mean of the surrogate statistic for KS
    % stat
    tpv(has_tfits) = cellfun(@(x) (x.min_KS_stat - mean(x.ks_surrogate))/std(x.ks_surrogate), ... 
        ava_out{ii}.all_tau_pfit(has_tfits));
    tau_KS_nSD(:, :, ii) = tpv;
    
    %% get the gamma vs cutoff values
    has_gfits = cellfun(@(x) numel(x) > 1, ava_out{ii}.all_gamma_pfit);
    gpv = cell(size(has_gfits));
    gpv_se = gpv;
    gpv(has_gfits) = cellfun(@(x) cellfun(@(y) y.mhat, x), ...
        ava_out{ii}.all_gamma_pfit(has_gfits), 'UniformOutput', false);
    gpv_se(has_gfits) = cellfun(@(x) cellfun(@(y) y.mSE, x), ...
        ava_out{ii}.all_gamma_pfit(has_gfits), 'UniformOutput', false);
    
    % scaling quality is SD of mean fit parameters 
    %%
    gpv_scq = nan(size(has_tfits));
    gpv_scq(has_gfits) = cellfun(@(x) std(x(2:min(5, end))) , ... 
        gpv(has_gfits));
    gamma_scaling(:, :, ii) = gpv_scq;
    gvs = nan(size(has_gfits));
    gvs(has_gfits) = cellfun(@(x) mean(x(2:min(5, end))) , ... 
        gpv(has_gfits));
    gamma_values(:, :, ii) = gvs;
end
%%
gamma_prediction = (alpha_values - 1)./(tau_values - 1);

crackle_err = abs(gamma_prediction - gamma_values)./gamma_scaling;


%%


figure()
for ii = 1:3
    subplot(3, 2, (ii-1)*2 + 1)
    has_scaling = tau_KS_nSD(:, :, ii) ;
    imagescWNans(has_scaling, [0 5])
    colorbar
    
    subplot(3, 2, (ii-1)*2 + 2)
    has_scaling = alpha_KS_nSD(:, :, ii);
    imagescWNans(has_scaling, [0 5])
    colorbar
end

%% OVERLAY on fig 1 where we see scaling, by run

figure(h_fig1);
for ii = 1:nR
    for jj = 1:nC
        subplot(nR, nC, (ii-1)*nC + jj)
        hold on
        for kk = 1:length(sim_set)
            if alpha_KS_nSD(ii, jj, kk) < 3
                text(kk, 0.2, '\alpha', 'Color', sim_set_colors(kk, :))
            end
            
            if tau_KS_nSD(ii, jj, kk) < 3
                text(kk, 0.3, '\tau', 'Color', sim_set_colors(kk, :))
            end
            
            if crackle_err(ii, jj, kk) < 1
                text(kk, 0.4, 'c', 'Color', sim_set_colors(kk, :), ...
                    'FontWeight', 'bold')
            elseif crackle_err(ii, jj, kk) < 2
                text(kk, 0.4, 'c', 'color', sim_set_colors(kk, :))
            end
            
        end
    end
end
%%
print(h_fig1, '-dpdf', ['avalanches_matlab_code/plots/paper_figures/fr_overview_' ...
    datestr(now, 'yyyymmdd')], '-painters') 
