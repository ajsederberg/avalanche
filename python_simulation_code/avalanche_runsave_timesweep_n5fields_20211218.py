#Runs and saves results of sweeps over time constants for 5 fields. 
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

# SET EPS and ETA here
# when using 'runsim_avalanche_env_fixedJh', the hamiltonian is H = eta*Jh+epsilon
# when using 'runsim_avalanche_env', the hamiltonian is H = eta*(Jh+epsilon)
epsilon = -3.0
# timed out , but also ran for [-10.0, -20.0]. Rerunning SLURM with more time. 
times = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0,  3.0, 5.0, 10.0, 20.0]
time_list = ['0.1', '0.2', '0.3', '0.5', '1.0', '2.0', '3.0', '5.0', '10.0', '20.0']
eta = 4.0 # multiply by this constant to increase overall activity of network
#MAKE SURE EPS and ETA match the file prefix, and nstim 
file_prefix = '_stim5e-12.0et4.0ph1.0p1.0t'


# MAKE SURE THIS MATCHES THE SAVE DIRECTORY
nstim = 5 #np.arange(2,21) # number of nonplace stimuli

datadir_name = '/home/sede0018/sede0018/avalanche_project/data/'
save_dir = '/home/sede0018/sede0018/avalanche_project/data/fields5timesweep/'

# replicate list
rep_list = ['A','B','C','D','E']

# set a keyword for this run: 
file_keyword = ''
for i in range(len(time_list)):
    for j in range(len(rep_list)):
        time = times[i]
        inp_string = 'time' + rep_list[j]
        runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel = inp_string)
        file_name = inp_string + file_prefix + time_list[i]
        try:
            save_sim_stats(file_name, datadir_name, save_dir)
            save_avalanches(file_name, datadir_name, save_dir)
            print(file_name)
        except:
            print(file_name + ' simulation files missing (likely too few events)')