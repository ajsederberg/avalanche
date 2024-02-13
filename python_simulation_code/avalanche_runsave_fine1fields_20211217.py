#Run and save results of sweeps over eps, eta for 5 fields. 
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


# timed out , but also ran for [-10.0, -20.0]. Rerunning SLURM with more time. 
time = 'stationary'

#these should match (if running simulations and saving script together) - be aware of hamiltonian def. 
# when using 'runsim_avalanche_env_fixedJh', the hamiltonian is H = eta*Jh+epsilon
#epsilons = [-20.0, -18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0]
#eps_list = ['-20.0', '-18.0', '-16.0', '-14.0', '-12.0', '-10.0', '-8.0', '-6.0', '-4.0', '-2.0']
epsilons = [-6.0, -4.0, -2.0]
eps_list = ['-6.0', '-4.0', '-2.0']
etas= [2.0, 4.0, 6.0, 8.0, 10.0] # multiply by this constant to increase overall activity 
eta_list  = ['2.0', '4.0', '6.0', '8.0', '10.0'] 


nstim = 1 #np.arange(2,21) # number of nonplace stimuli

datadir_name = '/home/sede0018/sede0018/avalanche_project/data/'
save_dir = '/home/sede0018/sede0018/avalanche_project/data/fields1finesweep/'

file_prefix = 'envfixedJh_sweep_fine_av_stim1'

# set a keyword for this run: 
file_keyword = 'fine'
for ii in range(len(epsilons)):
    for jj in range(len(etas)):
        
        runsim_avalanche_env_fixedJh(N0, nstim, percell, time, phi, etas[jj], epsilons[ii], keyword = file_keyword)

        file_name = file_prefix + 'e' + eps_list[ii] + 'et' + eta_list[jj] + 'ph1.0p1.0tstat'
        try:
            save_sim_stats(file_name, datadir_name, save_dir)
            save_avalanches(file_name, datadir_name, save_dir)
            print(file_name)
        except:
            print(file_name + ' simulation files missing (likely too few events)')
