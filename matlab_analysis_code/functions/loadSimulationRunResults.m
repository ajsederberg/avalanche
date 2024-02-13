function [x_info, num_ava] = loadSimulationRunResults(sim_string, rep_list)

    if nargin == 1
        rep_list = {''};
        ss_str = selectSimulation(sim_string);
    else 
        ss_str = selectSimulation(sim_string, rep_list{1});

    end
%     rep_list = {''};
    setLoadingFunctionNames

    x_info = cell(length(rep_list), length(y_var_list), length(x_var_list));
    num_ava = zeros(length(rep_list), length(y_var_list), length(x_var_list));
    for iR = 1:length(rep_list)
        ss_str = selectSimulation(sim_string, rep_list{iR});
        setLoadingFunctionNames

    %     data_dir = 'avalanches_code/data/fields5finesweep/';
    % 
    %     % file prefix string
    %     name_str = 'envfixedJh_sweep_fine_av_stim5' ; 
    %     dir_name_str ='fields5finesweep';
    %     par_title_str = '5 fields, infinite time constant, '
    %    
    %     % parameter values
    %     y_var_list = -2:-2:-20;
    %     x_var_list = 2:2:10; 

    % eta_list = [4.0, 6.0, 8.0];
    % eps_list = [-4.0,-5.0, -6.0, -7.0, -8.0, -9.0, -10.0, -12.0];
    % eps_list = [ -10.0, -12.0];

        for i_y = 1:length(y_var_list)
            for i_x = 1:length(x_var_list)

                try
                    x_ava = load(ava_name_fn(i_y, i_x));
                    num_ava(iR, i_y, i_x) = length(x_ava.sizes);
                end
                try
                    x_info{iR, i_y, i_x} = load(fr_name_fn(i_y, i_x));
                end
            end
        end
    end

end
