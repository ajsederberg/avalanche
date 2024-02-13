#Run sweeps over number of fields
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *


N0 = 1024 # number of cells
N = 256 # number of cells after removing silent ones

percell= 1.0 # probability that each field is accepted

phi=1.0 # multiply by this constant to adjust overall activity of 
# nonplace cells

epsilons= [-3.0, -1.0]

time = 'stationary'
    
etas= [5.0, 15./2., 10.0] # multiply by this constant to increase overall activity of 
# network
#30


nstim = 20 #np.arange(2,21) # number of nonplace stimuli

# runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=100, keyword = 'noreset'):
# set a keyword for this run: 
file_keyword = 'fields20limswp'
for i in range(len(epsilons)):
    for j in range(len(etas)):
        runsim_avalanche_env(N0, nstim, percell, time, phi, etas[j], epsilons[i], [], N, file_keyword)