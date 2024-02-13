#Save results of sweeps over time constants for 20 fields. This accompanies 'avalanche_timesweep_run20fields_20211216.py'
#in future versions, the saving functions should be run at the same time as simulation. 
#Date written and executed: 2021-12-17
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *
from placerg.avafunc import *


N0 = 4096 # number of cells
N = 1024 # number of cells after removing silent ones

percell= 1.0 # probability that each field is accepted

phi=1.0 # multiply by this constant to adjust overall activity of 
# nonplace cells

# when using 'runsim_avalanche_env_fixedJh', the hamiltonian is H = eta*Jh+epsilon
# when using 'runsim_avalanche_env', the hamiltonian is H = eta*(Jh+epsilon)
epsilon = -3.0
# timed out , but also ran for [-10.0, -20.0]. Rerunning SLURM with more time. 
times = [0.1, 0.3, 1.0, 3.0, 10.0]
time_list = ['0.1', '0.3', '1.0', '3.0', '10.0']
eta = 10.0 # multiply by this constant to increase overall activity of 
# network
#30


nstim = 20 #np.arange(2,21) # number of nonplace stimuli

datadir_name = '/home/sede0018/sede0018/avalanche_project/data/'
save_dir = '/home/sede0018/sede0018/avalanche_project/data/fields20timesweep/'

# file_prefix = 'env_noreset_av_stim20'
fileA_prefix = 'timeA_stim20e-30.0et10.0ph1.0p1.0t'
fileB_prefix = 'timeB_stim20e-30.0et10.0ph1.0p1.0t'


# set a keyword for this run: 
file_keyword = ''
for i in range(len(time_list)):
#     time = times[i]
#     runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel = 'timeA')
#     runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel = 'timeB')
#     fileA_name = fileA_prefix + time_list[i]
#     try:
#         save_sim_stats(fileA_name, datadir_name, save_dir)
#         save_avalanches(fileA_name, datadir_name, save_dir)
#         print(fileA_name)
#     except:
#         print(fileA_name + ' simulation files missing (likely too few events)')
        
    fileB_name = fileB_prefix + time_list[i]    
    try:
        save_sim_stats(fileB_name, datadir_name, save_dir)
        save_avalanches(fileB_name, datadir_name, save_dir)
        print(fileB_name)
    except:
        print(fileB_name + ' simulation files missing (likely too few events)')