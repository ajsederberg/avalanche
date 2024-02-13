
if ss_str.run_f20_e3010_timesweep
    data_dir = 'avalanches_code/data/fields20timesweep/';
    
    % file prefix string: 20 fields, eps = -30, eta = 10
    name_str = '_stim20e-30.0et10.0ph1.0p1.0';
    dir_name_str ='timesweep_stim20_e30et10';
    par_title_str = '20 fields, \epsilon = -30, \eta = 10, ';
    % time constant list and replicate list
    x_var_list = [0.1, 0.3, 1.0, 3.0, 10.0];
    y_var_list = {'A', 'B'};   % two replicates of the simulation
    % generate parameter string from parameter values
    param_str = @(i_rep, i_time) ['time' y_var_list{i_rep} name_str ...
        't' num2str(x_var_list(i_time), '%1.1f')];
    
    % x and y labels of summary plots
    x_lab_sum = 'time constant';
    y_lab_sum = 'replicate';
elseif ss_str.run_f5_e1204_timesweep
    data_dir = 'avalanches_code/data/fields5timesweep/';
    
    % file prefix string: 20 fields, eps = -30, eta = 10
    name_str = '_stim5e-12.0et4.0ph1.0p1.0';
    dir_name_str ='timesweep_stim5_e12et04';
    par_title_str = '5 fields, \epsilon = -12, \eta = 4, ';
    % time constant list and replicate list
    x_var_list = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0];
    y_var_list = {'A', 'B', 'C', 'D', 'E'};   % five replicates of the simulation
    % generate parameter string from parameter values
    param_str = @(i_rep, i_time) ['time' y_var_list{i_rep} name_str ...
        't' num2str(x_var_list(i_time), '%1.1f')];
    
    % x and y labels of summary plots
    x_lab_sum = 'time constant';
    y_lab_sum = 'replicate';
elseif ss_str.run_f1_e1204_timesweep
    data_dir = 'avalanches_code/data/fields1timesweep/';
    
    % file prefix string: 20 fields, eps = -30, eta = 10
    name_str = '_stim1e-8.0et4.0ph1.0p1.0';
    dir_name_str ='timesweep_stim1_e12et04';
    par_title_str = '1 fields, \epsilon = -8, \eta = 4, ';
    % time constant list and replicate list
    x_var_list = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0];
    y_var_list = {'A', 'B', 'C', 'D', 'E'};   % five replicates of the simulation
    
    % generate parameter string from parameter values
    param_str = @(i_rep, i_time) ['time' y_var_list{i_rep} name_str ...
        't' num2str(x_var_list(i_time), '%1.1f')];

    
    % x and y labels of summary plots
    x_lab_sum = 'time constant';
    y_lab_sum = 'replicate';
elseif ss_str.run_f1finesweep
    data_dir = 'avalanches_code/data/fields1finesweep/';

    % file prefix string
    name_str = 'envfixedJh_sweep_fine_av_stim1' ; 

    dir_name_str ='fields1finesweep';
    par_title_str = '1 field, infinite time constant, ';
    
    % parameter values
    y_var_list = -2:-2:-20;
    x_var_list = 2:2:10;    
    
    % generate parameter string from parameter values
    param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
        num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];
    
    % x and y labels of summary plots
    x_lab_sum = '\eta';
    y_lab_sum = '\epsilon';
    
    ava_results_file = 'ava_decade_analysis_20220304.mat'; %'ava_decade_analysis_20220217.mat';
elseif ss_str.run_f5finesweep
    data_dir = 'avalanches_code/data/fields5finesweep/';

    % file prefix string
    name_str = 'envfixedJh_sweep_fine_av_stim5' ; 
    dir_name_str ='fields5finesweep';
    par_title_str = '5 fields, infinite time constant, ';
   
    % parameter values
    y_var_list = -2:-2:-20;
    x_var_list = 2:2:10; 
    
    % generate parameter string from parameter values
    param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
        num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];
    
    % x and y labels of summary plots
    x_lab_sum = '\eta';
    y_lab_sum = '\epsilon';
    
    ava_results_file = 'ava_decade_analysis_20220226.mat'; %'ava_decade_analysis_20220217.mat';

elseif ss_str.run_f20finesweep
    data_dir = 'avalanches_code/data/fields20finesweep/';

    % file prefix string
    name_str = 'envfixedJh_sweep_fine_av_stim20' ; 
    dir_name_str ='fields20finesweep';
    par_title_str = '20 fields, infinite time constant, ';
      
    % parameter values
    y_var_list = -2:-2:-20;
    x_var_list = 2:2:10;
    
    % generate parameter string from parameter values
    param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
        num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];

    % x and y labels of summary plots
    x_lab_sum = '\eta';
    y_lab_sum = '\epsilon';
elseif ss_str.run_f5ultrafinesweep
    % ss_str will have a replicate field
    try 
        rep_str = ['rep' ss_str.replicate];
        
        data_dir = 'avalanches_code/data/fields5ultrafinesweep/';

    % file prefix string
        name_str = ['sweep_ultrafine_' rep_str '_stim5']; 
        dir_name_str =['fields5ultrafinesweep/' rep_str '/'];
        par_title_str = ['5 fields, infinite time constant, ' rep_str ', '];
   
    % parameter values
    %     y_var_list = -2:-2:-20;
    %     x_var_list = 2:2:10; 
        y_var_list = [-6 -7 -8 -9 -10 -12 -14 -16 -18 -20];
        x_var_list = [4 6 8];    
        % generate parameter string from parameter values
        param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
            num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];

        % x and y labels of summary plots
        x_lab_sum = '\eta';
        y_lab_sum = '\epsilon';
    catch
        disp('MUST PROVIDE REPLICATE: ss_str.replicate = "A"; etc.')
    end
elseif ss_str.run_f1ultrafinesweep
    % ss_str will have a replicate field
    try 
        rep_str = ['rep' ss_str.replicate];
        
        data_dir = 'avalanches_code/data/fields1ultrafinesweep/';

    % file prefix string
        name_str = ['sweep_ultrafine_' rep_str '_stim1']; 
        dir_name_str =['fields1ultrafinesweep/' rep_str '/'];
        par_title_str = ['1 field, infinite time constant, ' rep_str ', '];
   
    % parameter values
    %     y_var_list = -2:-2:-20;
    %     x_var_list = 2:2:10; 
        y_var_list = [-6 -8 -10 -12 -14 -16 -18 -20];
        x_var_list = [4 6 8];    
        % generate parameter string from parameter values
        param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
            num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];

        % x and y labels of summary plots
        x_lab_sum = '\eta';
        y_lab_sum = '\epsilon';
        
%         ava_results_file = 'ava_decade_analysis_20220426.mat'; %'ava_decade_analysis_20220217.mat';
    catch
        disp('MUST PROVIDE REPLICATE: ss_str.replicate = "A"; etc.')
    end
elseif ss_str.run_f1smallfinesweep
    % ss_str will have a replicate field
    try 
        rep_str = ['rep' ss_str.replicate];
        
        data_dir = 'avalanches_code/data/fields1smallfinesweep/';

    % file prefix string
        name_str = ['sweep_smallfine_' rep_str '_stim1']; 
        dir_name_str =['fields1smallfinesweep/' rep_str '/'];
        par_title_str = ['1 field, small N, infinite time constant, ' rep_str ', '];
   
    % parameter values
        y_var_list = -2:-2:-14;
        x_var_list = 1:10; 
%         y_var_list = [-6 -8 -10 -12 -14 -16 -18 -20];
%         x_var_list = [4 6 8];    
        % generate parameter string from parameter values
        param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
            num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];

        % x and y labels of summary plots
        x_lab_sum = '\eta';
        y_lab_sum = '\epsilon';
        
%         ava_results_file = 'ava_decade_analysis_20220426.mat'; %'ava_decade_analysis_20220217.mat';
    catch
        disp('MUST PROVIDE REPLICATE: ss_str.replicate = "A"; etc.')
    end

elseif ss_str.run_f1smallultrafinesweep
    % ss_str will have a replicate field
    try 
        rep_str = ['rep' ss_str.replicate];
        
        data_dir = 'avalanches_code/data/fields1smallultrafinesweep/';

    % file prefix string
        name_str = ['sweep_smallfine_' rep_str '_stim1']; 
        dir_name_str =['fields1smallultrafinesweep/' rep_str '/'];
        par_title_str = ['1 field, small N, infinite time constant, ' rep_str ', '];
   
    % parameter values
        y_var_list = -2:-.5:-5.5;
        x_var_list = 1:2; 
%         y_var_list = [-6 -8 -10 -12 -14 -16 -18 -20];
%         x_var_list = [4 6 8];    
        % generate parameter string from parameter values
        param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
            num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];

        % x and y labels of summary plots
        x_lab_sum = '\eta';
        y_lab_sum = '\epsilon';
        
%         ava_results_file = 'ava_decade_analysis_20220426.mat'; %'ava_decade_analysis_20220217.mat';
    catch
        disp('MUST PROVIDE REPLICATE: ss_str.replicate = "A"; etc.')
    end

elseif ss_str.run_morrell1field
    data_dir = 'avalanches_code/data/avalanche_data_Morrell/';

    % file prefix string
    name_str = 'envfixedfieldJ_md2dsweep_etaepssweep_coarse_av_stim1' ; 
% e-5.0et5.0ph1.0p1.0tstat
    dir_name_str ='fields1morrell';
    par_title_str = '1 field, infinite time constant, ';
    
    % parameter values
    y_var_list = -5:-2.5:-25;
    x_var_list = 5:5:50;    
    
    % generate parameter string from parameter values
    param_str = @(i_eps, i_eta) [name_str 'e' num2str(y_var_list(i_eps), '%1.1f') 'et' ...
        num2str(x_var_list(i_eta), '%1.1f') 'ph1.0p1.0tstat'];
    
    % x and y labels of summary plots
    x_lab_sum = '\eta';
    y_lab_sum = '\epsilon';



end

plot_dir = ['avalanches_matlab_code/plots/' dir_name_str ];
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

save_dir = ['avalanches_matlab_code/analysis_results/' dir_name_str ];
if ~strcmp(save_dir(end), '/')
    save_dir = [save_dir '/'];
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

% look for the avalanche analysis files
file_list = dir([save_dir 'ava_decade_analysis_*.mat']);
if ~isempty(file_list)
    % if there are multiple results files here, take the most recent
    date_nums = arrayfun(@(x) x.datenum, file_list);
    [~, ord] = sort(date_nums, 'descend');
    ava_results_file = file_list(ord(1)).name;
end


% end of file name  string
ava_str_end = '.mat';
fr_str_end = 'fr_stats.mat';

% this function generates the file name for the specified parameter INDICES
ava_name_fn = @(i_y, i_x) [data_dir param_str(i_y, i_x) ava_str_end];
fr_name_fn = @(i_y, i_x) [data_dir param_str(i_y, i_x) fr_str_end];
