#Run and save results of sweeps over eps, eta for 5 fields. 
#Date written and executed: 2022-04-26
#Runs multiple replicates at each set of parameter values
# ran test run with two replicates at one set of parameter values
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *
from placerg.avafunc import *


N0 = 128 # number of cells
N = 128 # number of cells after removing silent ones

percell= 1.0 # probability that each field is accepted

phi=1.0 # multiply by this constant to adjust overall activity of 
# nonplace cells


# timed out , but also ran for [-10.0, -20.0]. Rerunning SLURM with more time. 
time = 'stationary'

#these should match (if running simulations and saving script together) - be aware of hamiltonian def. 
# when using 'runsim_avalanche_env_fixedJh', the hamiltonian is H = eta*Jh+epsilon
# epsilons = [ -12.0,   -10.0,   -9.0,   -8.0,  -7.0,   -6.0]
# eps_list = ['-12.0', '-10.0', '-9.0', '-8.0','-7.0', '-6.0']
epsilons = [ -6.0]
eps_list = ['-6.0']

# etas = [3.0, 2.0, 5.0]
# eta_list = ['3.0', '2.0', '5.0']


etas= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] # multiply by this constant to increase overall activity 
eta_list  = ['1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0'] 

rep_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
# rep_list = ['A', 'B', 'C', 'D', 'E']


nstim = 1 #np.arange(2,21) # number of nonplace stimuli

datadir_name = '/Users/sederberglab/Dropbox/Projects_Python/avalanches/data/'
save_dir = '/Users/sederberglab/Dropbox/Projects_Python/avalanches/data/fields1smallfinesweep/'

# provide to runsim_avalanche_env_fixedJh
file_name_start = 'sweep_smallfine'

for ii in range(len(epsilons)):
    for jj in range(len(etas)):
        for kk in range(len(rep_list)):
            # set a keyword for this run: 
            file_keyword = 'rep' + rep_list[kk]
           
            runsim_avalanche_env(N0, nstim, percell, time, phi, etas[jj], epsilons[ii], inputlabel = file_name_start +'_' + file_keyword, loop = 200)
                
            file_name = file_name_start +'_' + file_keyword + '_stim1' + 'e' + eps_list[ii] + 'et' + eta_list[jj] + 'ph1.0p1.0tstat'
            #try:
            print(file_name)

            save_sim_stats(file_name, datadir_name, save_dir)

            try: 
                save_avalanches(file_name, datadir_name, save_dir)
            except: 
                print(file_name + ' avalanche analysis failed') 
            #except:
                #print(file_name + ' simulation files missing (likely too few events)')
