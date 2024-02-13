function ss_str = selectSimulation(true_field, rep_field)
% sets the specified field to true in s_str. All other fields are false.
% These match field names in setLoadingFunctionNames and correspond to
% specific simulation runs. 


% initialize ss_str
ss_str.run_morrell1field = false;       % analyze the envfield_fixedJ sims from Mia
ss_str.run_f20_e3010_timesweep = false;       % pull data from 20-field, time sweep
ss_str.run_f5_e1204_timesweep = false;       % pull data from 5-field, eps = -12, eta = 4 time sweep
ss_str.run_f1_e1204_timesweep = false;       % pull data from 1-field, eps = -12, eta = 4 time sweep
ss_str.run_f1finesweep = false;       % pull data from 1-field, fine sweep of eta, eps
ss_str.run_f5finesweep = false;       % pull data from 5-field, fine sweep
ss_str.run_f20finesweep = false;       % pull data from 20-field, fine sweep
ss_str.run_f5ultrafinesweep = false;       % pull data from 5-field, ultrafine sweep
ss_str.run_f1ultrafinesweep = false;       % pull data from 1-field, ultrafine sweep
ss_str.run_f1smallfinesweep = false;       % pull data from 1-field, small network, fine sweep
ss_str.run_f1smallultrafinesweep = false;       % pull data from 1-field, small network, ultrafine limited sweep

% set the "true_field" to true
ss_str.(true_field) = true;

if nargin > 1
    ss_str.replicate = rep_field;
end
