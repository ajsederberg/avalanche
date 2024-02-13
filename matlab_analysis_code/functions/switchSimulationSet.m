function ss_str = switchSimulationSet(ss_str, field_name)
% sets the specified field to true, all others to false. 
% for selecting which of the simulation sets to load
% options may include:

% run_f20_e3010_timesweep = false;       % pull data from 20-field, time sweep
% run_f5_e1204_timesweep = false;       % pull data from 5-field, eps = -12, eta = 4 time sweep
% run_f1_e1204_timesweep = false;       % pull data from 1-field, eps = -12, eta = 4 time sweep
% 
% run_f1finesweep = false;       % pull data from 1-field, fine sweep of eta, eps
% run_f5finesweep = false;       % pull data from 5-field, fine sweep
% run_f20finesweep = true;       % pull data from 20-field, fine sweep


all_fields = fieldnames(ss_str);
for ii = 1:length(all_fields)
    ss_str.(all_fields{ii}) = false;
end

ss_str.(field_name) = true;