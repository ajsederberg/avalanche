#Run sweeps over number of fields
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *


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
    
eta = 10.0 # multiply by this constant to increase overall activity of 
# network
#30


nstim = 20 #np.arange(2,21) # number of nonplace stimuli

# runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=100, keyword = 'noreset'):
# runsim_avalanche_env_fixedfieldJ(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset'):
# set a keyword for this run: 
file_keyword = ''
for i in range(len(times)):
    time = times[i]
#         runsim_avalanche_env_fixedJh(N0, nstim, percell, time, phi, etas[j], epsilons[i], keyword = file_keyword)
    runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel = 'timeA')
    runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel = 'timeB')