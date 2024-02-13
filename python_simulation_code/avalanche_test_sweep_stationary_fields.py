#Run sweeps over number of fields
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *


N0 = 500 # number of cells

percell= 1.0 # probability that each field is accepted

phi=1.0 # multiply by this constant to adjust overall activity of 
# nonplace cells

epsilon= -3.0

time = 'stationary'
    
eta= 15./2. # multiply by this constant to increase overall activity of 
# network
#30


nstims = [20] #np.arange(2,21) # number of nonplace stimuli


for i in range(len(nstims)):
    nstim = nstims[i]
    runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon)


