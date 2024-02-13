#Run sweeps over number of fields
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *


N0 = 2048 # number of cells
N = 512 # number of cells after removing silent ones

percell= 1.0 # probability that each field is accepted

phi=1.0 # multiply by this constant to adjust overall activity of 
# nonplace cells

# when using 'runsim_avalanche_env_fixedJh', the hamiltonian is H = eta*Jh+epsilon
epsilons= [ -30.0, -40.0, -50.0]
# timed out , but also ran for [-10.0, -20.0]. Rerunning SLURM with more time. 
time = 'stationary'
    
etas= [5.0, 10.0, 15.0, 20.0, 25.0] # multiply by this constant to increase overall activity of 
# network
#30


nstim = 20 #np.arange(2,21) # number of nonplace stimuli

# runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=100, keyword = 'noreset'):
# runsim_avalanche_env_fixedfieldJ(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset'):
# set a keyword for this run: 
file_keyword = ''
for i in range(len(epsilons)):
    for j in range(len(etas)):
#         runsim_avalanche_env(N0, nstim, percell, time, phi, etas[j], epsilons[i], [], N, file_keyword)
#         runsim_avalanche_env_fixedfieldJ(N0, nstim, percell, time, phi, eta, epsilon, keyword = keyword)
# runsim_avalanche_env_fixedJh(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset'):
        runsim_avalanche_env_fixedJh(N0, nstim, percell, time, phi, etas[j], epsilons[i], keyword = file_keyword)
