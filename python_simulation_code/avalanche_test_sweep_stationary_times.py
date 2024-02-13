#Run sweeps over time constants
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.runfunc import *


N0 = 1024*4 # number of cells

percell= 1.0 # probability that each field is accepted

phi=1.0 # multiply by this constant to adjust overall activity of 
# nonplace cells

epsilon= -3.
        
#times = [0.1, 0.25, 0.5, 0.75, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, \
#        110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 200.0, 0.05, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.35, 0.4,\
#        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, \
#        5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5 , 0.055, 0.065, 0.075, 0.085, 0.095, 0.125, 0.15, 0.225, 0.325, 0.425, 0.525, 0.625, 0.725, 0.825, 0.925, 1.25, 2.25, 3.25, 4.25, \
#        5.25, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, \
#        47.5, 52.5, 57.5]

#times = [ 0.05, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.35, 0.4,\
#        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, \
#        5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5 , 0.055, 0.065, 0.075, 0.085, 0.095, 0.125, 0.15, 0.225, 0.325, 0.425, 0.525, 0.625, 0.725, 0.825, 0.925, 1.25, 2.25, 3.25, 4.25, \
#        5.25, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, \
#        47.5, 52.5, 57.5]

times = ['stationary']

#times = [42.5, 47.5, 52.5]

#times = [50.0, 55.0]

#times = [60.0, 65.0, 70.0, 75.0]
#times = [75.0, 80.0]


    
eta= 15./2. # multiply by this constant to increase overall activity of 
# network
#30


nstim = 1 # number of nonplace stimuli


for i in range(len(times)):
    time = times[i]
    runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon)


