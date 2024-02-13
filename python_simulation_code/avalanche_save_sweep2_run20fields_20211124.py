#Save results of sweeps over eps, eta for 20 fields. This accompanies 'avalanche_sweep2_run20fields_20211124.py'
#in future versions, the saving functions should be run at the same time as simulation. 
#Date written and executed: 2021-12-17
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *
from placerg.avafunc import *

N0 = 2048 # number of cells
N = 512 # number of cells after removing silent ones

percell= 1.0 # probability that each field is accepted

phi=1.0 # multiply by this constant to adjust overall activity of 
# nonplace cells


# timed out , but also ran for [-10.0, -20.0]. Rerunning SLURM with more time. 
time = 'stationary'

#these should match (if running simulations and saving script together) - be aware of hamiltonian def. 
# when using 'runsim_avalanche_env_fixedJh', the hamiltonian is H = eta*Jh+epsilon
epsilons= [-10.0, -20.0, -30.0, -40.0, -50.0]
eps_list = ['-10.0', '-20.0', '-30.0', '-40.0', '-50.0']#['-15.0', '-5.0', '-22.5', '-7.5', '-30.0', '-10.0']

etas= [5.0, 10.0, 15.0, 20.0, 25.0] # multiply by this constant to increase overall activity of 
eta_list  = ['5.0', '10.0', '15.0', '20.0', '25.0'] #[  '5.0',  '5.0',   '7.5',  '7.5',  '10.0', '10.0']

# network
#30

nstim = 20 #np.arange(2,21) # number of nonplace stimuli

datadir_name = '/home/sede0018/sede0018/avalanche_project/data/'
save_dir = '/home/sede0018/sede0018/avalanche_project/data/fields20sweep/'

# file_prefix = 'env_noreset_av_stim20'
file_prefix = 'envfixedJh_sweep__av_stim20'
# runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=100, keyword = 'noreset'):
# runsim_avalanche_env_fixedfieldJ(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset'):
# set a keyword for this run: 
file_keyword = ''
for ii in range(len(epsilons)):
    for jj in range(len(etas)):
        file_name = file_prefix + 'e' + eps_list[ii] + 'et' + eta_list[jj] + 'ph1.0p1.0tstat'
        try:
            save_sim_stats(file_name, datadir_name, save_dir)
            save_avalanches(file_name, datadir_name, save_dir)
            print(file_name)
        except:
            print(file_name + ' simulation files missing (likely too few events)')