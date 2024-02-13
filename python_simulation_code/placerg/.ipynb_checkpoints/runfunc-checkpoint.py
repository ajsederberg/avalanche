#runfunc.py
import numpy as np
from placerg.funcs import *
from placerg.funcsrg import *
from placerg.objects import *
from placerg.funcsboot import *
import glob

def return_target_directory():
    # target_dir = '/home/sede0018/sede0018/avalanche_project/data/'
    # target_dir = '/Users/sederberglab/Dropbox/Projects_Python/avalanches/data/'
    target_dir = '/Users/sederberglab/Dropbox/Projects_Python/LVM/data/'
    return target_dir

    
def runsim_avalanche_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=128, keyword = 'noreset', loop = 40):
    """
    Run a simulation. fields and couplings are selected randomly. same as runsim_avalanche_env_sweep but with different sim filenaming
    --------------------------------------------------------------
    Inputs:
    --------------------------------------------------------------
    N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
    N: number of cells to save
    nstim: number of latent fields
    percell: probability that each field couples to a given cell
    time: timeconstant of latent field
    phi: field multiplier
    eta: variance multiplier
    epsilon: bias
    inputlabel: label for sims filenames
    keyword: another label for sims filenames
    """
    N0 = N0 # number of cells

    N = N # number of cells after removing silent ones
 
    # loop = 40 #1600 # number of track runs; changed to an input onn 4/26/22

    # for the avalanche function to work, this MUST be equal to 10000. (4/26/22)
    xmax = 10000 #9600 # 50 track length

    dt=1. # increment for measurement locations

    # set up network structure

    nstim = nstim # number of nonplace stimuli

    percell= percell # probability that each field is accepted

    
    if nstim==0:
        npprob=np.array([1,0])
    else:
        npprob=np.array([1-percell, percell]) 
    # probability that cell is coupled to percell fields
    

    # now initialize distribution parameters of fields and couplings
    # here the place field couplings will be gamma distributed
    # the nonplace field couplings will be normally distributed

    vj = 0. # mean of couplings for nonplace fields
  
   

    sj = 1. # standard deviation of couplings for nonplace fields


    time=time
    if time == 'stationary':
        timelabel = 'stat'
        timeconst=time
    elif type(time)!=float:
        print('running with non-uniform time constants')
        timeconst=time
        timelabel=inputlabel
        print(timelabel)
        
        
    else: 
        timeconst = np.full((nstim,), time) # mean length of stochastic process in track lengths
        timelabel=time

        # the call 'avalanche_runsim' uses a Hamiltonian H = eta*(Jh + epsilon). However, we have 
        # changed our convention to H = eta*Jh + epsilon, and this is reflected in the name_env variable below. 
    phi=phi # multiply by this constant to adjust overall activity of 
    # nonplace cells

    eta = eta # multiply by this constant to increase overall activity of 
    # network

    epsilon= epsilon


        #name_env='/home/mia/OneDrive/simsrg/env_av_stim{}e{}et{}ph{}p{}t{}'.format(nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel) #reset all every track run
    
#     name_env='/home/mia/OneDrive/simsrg/env_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
    
#     name_env=return_target_directory()+'env_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
    name_env=return_target_directory()+'{}_stim{}e{}et{}ph{}p{}t{}'.format(inputlabel, nstim, np.round(epsilon, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
    
#     name_env=return_target_directory()+'{}_stim{}e{}et{}ph{}p{}t{}'.format(inputlabel, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
    
    
    print(name_env)
    
    avalanche_runsim(N0=N0, N=N, loop=loop, xmax=xmax, dt=dt, nstim=nstim,\
    percell=percell, \
    npprob=npprob, vj=vj,  sj=sj,\
    timeconst=timeconst,\
    phi=phi, eta=eta, epsilon=epsilon, name = name_env, num_pages = 5)

def runsim_DLV_env(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=128, keyword = 'noreset', loop = 40):
    """
    Run a dynamic latent variable simulation. fields and couplings are selected randomly. 
    same as runsim_avalanche_env but with different sim filenaming
    --------------------------------------------------------------
    Inputs:
    --------------------------------------------------------------
    N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
    N: number of cells to save
    nstim: number of latent fields
    percell: probability that each field couples to a given cell
    time: timeconstant of latent field
    phi: field multiplier
    eta: variance multiplier
    epsilon: bias
    inputlabel: label for sims filenames
    keyword: another label for sims filenames
    """
    N0 = N0 # number of cells

    N = N # number of cells after removing silent ones
 
    # loop = 40 #1600 # number of track runs; changed to an input onn 4/26/22

    # for the avalanche function to work, this MUST be equal to 10000. (4/26/22)
    xmax = 1000 # track length

    dt=1. # increment for measurement locations

    # set up network structure

    nstim = nstim # number of nonplace stimuli

    percell= percell # probability that each field is accepted

    
    if nstim==0:
        npprob=np.array([1,0])
    else:
        npprob=np.array([1-percell, percell]) 
    # probability that cell is coupled to percell fields
    

    # now initialize distribution parameters of fields and couplings
    # here the place field couplings will be gamma distributed
    # the nonplace field couplings will be normally distributed

    vj = 0. # mean of couplings for nonplace fields
  
   

    sj = 1. # standard deviation of couplings for nonplace fields


    time=time
    if time == 'stationary':
        timelabel = 'stat'
        timeconst=time
    elif type(time)!=float:
        print('running with non-uniform time constants')
        timeconst=time
        timelabel=inputlabel
        print(timelabel)
        
        
    else: 
        timeconst = np.full((nstim,), time) # mean length of stochastic process in track lengths
        timelabel=time

        # the call 'avalanche_runsim' uses a Hamiltonian H = eta*(Jh + epsilon). However, we have 
        # changed our convention to H = eta*Jh + epsilon, and this is reflected in the name_env variable below. 
    phi=phi # multiply by this constant to adjust overall activity of 
    # nonplace cells

    eta = eta # multiply by this constant to increase overall activity of 
    # network

    epsilon= epsilon

    name_env=return_target_directory()+'{}_N{}stim{}e{}et{}ph{}p{}t{}'.format(inputlabel, N, nstim, np.round(epsilon, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
    
    
    print(name_env)
    
    avalanche_runsim(N0=N0, N=N, loop=loop, xmax=xmax, dt=dt, nstim=nstim,\
    percell=percell, \
    npprob=npprob, vj=vj,  sj=sj,\
    timeconst=timeconst,\
    phi=phi, eta=eta, epsilon=epsilon, name = name_env, num_pages = 5)

    return name_env

    
def runsim_avalanche_env_sweep(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='Nsweep',  N=1024):
    """
    Run a simulation. fields and couplings are selected randomly. same as runsim_avalanche_env but with different sim filenaming
    --------------------------------------------------------------
    Inputs:
    --------------------------------------------------------------
    N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
    N: number of cells to save
    nstim: number of latent fields
    percell: probability that each field couples to a given cell
    time: timeconstant of latent field
    phi: field multiplier
    eta: variance multiplier
    epsilon: bias
    inputlabel: label for sims filenames
    """
    N0 = N0 # number of cells

    N = N # number of cells after removing silent ones
 
    loop = 40 #1600 # number of track runs

    xmax = 10000 #9600 # 50 track length

    dt=1. # increment for measurement locations

    # set up network structure

    nstim = nstim # number of nonplace stimuli

    percell= percell # probability that each field is accepted

    
    if nstim==0:
        npprob=np.array([1,0])
    else:
        npprob=np.array([1-percell, percell]) 
    # probability that cell is coupled to percell fields
    

    # now initialize distribution parameters of fields and couplings
    # here the place field couplings will be gamma distributed
    # the nonplace field couplings will be normally distributed

    vj = 0. # mean of couplings for nonplace fields
  
   

    sj = 1. # standard deviation of couplings for nonplace fields


    time=time
    if time == 'stationary':
        timelabel = 'stat'
        timeconst=time
    elif type(time)!=float:
        print('running with non-uniform time constants')
        timeconst=time
        timelabel=inputlabel
        print(timelabel)
        
        
    else: 
        timeconst = np.full((nstim,), time) # mean length of stochastic process in track lengths
        timelabel=time

    phi=phi # multiply by this constant to adjust overall activity of 
    # nonplace cells

    eta = eta # multiply by this constant to increase overall activity of 
    # network

    epsilon= epsilon


    #name_env='/home/mia/OneDrive/simsrg/env_av_stim{}e{}et{}ph{}p{}t{}'.format(nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel) #reset all every track run
    
#     name_env='/home/mia/OneDrive/simsrg/env_{}_av_N{}stim{}e{}et{}ph{}p{}t{}'.format(inputlabel, N, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
    name_env=return_target_directory()+'env_{}_av_N{}stim{}e{}et{}ph{}p{}t{}'.format(inputlabel, N, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
    
    avalanche_runsim(N0=N0, N=N, loop=loop, xmax=xmax, dt=dt, nstim=nstim,\
    percell=percell, \
    npprob=npprob, vj=vj,  sj=sj,\
    timeconst=timeconst,\
    phi=phi, eta=eta, epsilon=epsilon, name = name_env, num_pages = 5)

    

def runsim_avalanche_env_fixedJ(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset'):
    """
    Run a simulation. fields are selected randomly. couplings are fixed.
    --------------------------------------------------------------
    Inputs:
    --------------------------------------------------------------
    N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
    N: number of cells to save
    nstim: number of latent fields
    percell: probability that each field couples to a given cell
    time: timeconstant of latent field
    phi: field multiplier
    eta: variance multiplier
    epsilon: bias
    inputlabel: label for sims filenames
    keyword: another label for sims filenames
    """
    N0 = N0 # number of cells

    N = N # number of cells after removing silent ones
 
    loop = 40 #1600 # number of track runs

    xmax = 10000 #9600 # 50 track length

    dt=1. # increment for measurement locations

    # set up network structure

    nstim = nstim # number of nonplace stimuli

    percell= percell # probability that each field is accepted

    
    if nstim==0:
        npprob=np.array([1,0])
    else:
        npprob=np.array([1-percell, percell]) 
    # probability that cell is coupled to percell fields
    

    # now initialize distribution parameters of fields and couplings
    # here the place field couplings will be gamma distributed
    # the nonplace field couplings will be normally distributed

    vj = 0. # mean of couplings for nonplace fields
  
   

    sj = 1. # standard deviation of couplings for nonplace fields


    time=time
    if time == 'stationary':
        timelabel = 'stat'
        timeconst=time
    elif type(time)!=float:
        print('running with non-uniform time constants')
        timeconst=time
        timelabel=inputlabel
        print(timelabel)
        
        
    else: 
        timeconst = np.full((nstim,), time) # mean length of stochastic process in track lengths
        timelabel=time

    phi=phi # multiply by this constant to adjust overall activity of 
    # nonplace cells

    eta = eta # multiply by this constant to increase overall activity of 
    # network

    epsilon= epsilon


    #name_env='/home/mia/OneDrive/simsrg/env_av_stim{}e{}et{}ph{}p{}t{}'.format(nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel) #reset all every track run
    
    #name_env='/home/mia/OneDrive/simsrg/envfixedJ_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
#     name_env='/home/mia/OneDrive/simsrg/envfixedJ_md2dsweep_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  #FOR UNCOUPLED MODEL

    name_env=return_target_directory()+'envfixedJ_md2dsweep_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  #FOR UNCOUPLED MODEL
    
    avalanche_runsim_fixedJ(N0=N0, N=N, loop=loop, xmax=xmax, dt=dt, nstim=nstim,\
    percell=percell, \
    npprob=npprob, vj=vj,  sj=sj,\
    timeconst=timeconst,\
    phi=phi, eta=eta, epsilon=epsilon, name = name_env, num_pages = 5) 
    
def runsim_avalanche_env_fixedfield(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset'):
    """
    Run a simulation. couplings are selected randomly. fields are fixed.
    --------------------------------------------------------------
    Inputs:
    --------------------------------------------------------------
    N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
    N: number of cells to save
    nstim: number of latent fields
    percell: probability that each field couples to a given cell
    time: timeconstant of latent field
    phi: field multiplier
    eta: variance multiplier
    epsilon: bias
    inputlabel: label for sims filenames
    keyword: another label for sims filenames
    """
    N0 = N0 # number of cells

    N = N # number of cells after removing silent ones
 
    loop = 40 #1600 # number of track runs

    xmax = 10000 #9600 # 50 track length

    dt=1. # increment for measurement locations

    # set up network structure

    nstim = nstim # number of nonplace stimuli

    percell= percell # probability that each field is accepted

    
    if nstim==0:
        npprob=np.array([1,0])
    else:
        npprob=np.array([1-percell, percell]) 
    # probability that cell is coupled to percell fields
    

    # now initialize distribution parameters of fields and couplings
    # here the place field couplings will be gamma distributed
    # the nonplace field couplings will be normally distributed

    vj = 0. # mean of couplings for nonplace fields
  
   

    sj = 1. # standard deviation of couplings for nonplace fields


    time=time
    if time == 'stationary':
        timelabel = 'stat'
        timeconst=time
    elif type(time)!=float:
        print('running with non-uniform time constants')
        timeconst=time
        timelabel=inputlabel
        print(timelabel)
        
        
    else: 
        timeconst = np.full((nstim,), time) # mean length of stochastic process in track lengths
        timelabel=time

    phi=phi # multiply by this constant to adjust overall activity of 
    # nonplace cells

    eta = eta # multiply by this constant to increase overall activity of 
    # network

    epsilon= epsilon


    #name_env='/home/mia/OneDrive/simsrg/env_av_stim{}e{}et{}ph{}p{}t{}'.format(nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel) #reset all every track run
    
    #name_env='/home/mia/OneDrive/simsrg/envfixedfield_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # FOR COUPLED HAMILTONIAN MODEL
#     name_env='/home/mia/OneDrive/simsrg/envfixedfield_md2dsweep_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  # FOR UNCOUPLED HAMILTONIAN MODEL
    
    name_env=return_target_directory()+'envfixedfield_md2dsweep_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  # FOR UNCOUPLED HAMILTONIAN MODEL
    
    avalanche_runsim_fixedfield(N0=N0, N=N, loop=loop, xmax=xmax, dt=dt, nstim=nstim,\
    percell=percell, \
    npprob=npprob, vj=vj,  sj=sj,\
    timeconst=timeconst,\
    phi=phi, eta=eta, epsilon=epsilon, name = name_env, num_pages = 5) 


def runsim_avalanche_env_fixedfieldJ(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset'):
    """
    Run a simulation. fields are fixed. couplings are fixed.
    --------------------------------------------------------------
    Inputs:
    --------------------------------------------------------------
    N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
    N: number of cells to save
    nstim: number of latent fields
    percell: probability that each field couples to a given cell
    time: timeconstant of latent field
    phi: field multiplier
    eta: variance multiplier
    epsilon: bias
    inputlabel: label for sims filenames
    keyword: another label for sims filenames
    """
    N0 = N0 # number of cells

    N = N # number of cells after removing silent ones
 
    loop = 40 #1600 # number of track runs

    xmax = 10000 #9600 # 50 track length

    dt=1. # increment for measurement locations

    # set up network structure

    nstim = nstim # number of nonplace stimuli

    percell= percell # probability that each field is accepted

    
    if nstim==0:
        npprob=np.array([1,0])
    else:
        npprob=np.array([1-percell, percell]) 
    # probability that cell is coupled to percell fields
    

    # now initialize distribution parameters of fields and couplings
    # here the place field couplings will be gamma distributed
    # the nonplace field couplings will be normally distributed

    vj = 0. # mean of couplings for nonplace fields
  
   

    sj = 1. # standard deviation of couplings for nonplace fields


    time=time
    if time == 'stationary':
        timelabel = 'stat'
        timeconst=time
    elif type(time)!=float:
        print('running with non-uniform time constants')
        timeconst=time
        timelabel=inputlabel
        print(timelabel)
        
        
    else: 
        timeconst = np.full((nstim,), time) # mean length of stochastic process in track lengths
        timelabel=time

    phi=phi # multiply by this constant to adjust overall activity of 
    # nonplace cells

    eta = eta # multiply by this constant to increase overall activity of 
    # network

    epsilon= epsilon


    #name_env='/home/mia/OneDrive/simsrg/env_av_stim{}e{}et{}ph{}p{}t{}'.format(nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel) #reset all every track run
    
    #name_env='/home/mia/OneDrive/simsrg/envfixedfieldJ_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
#     name_env='/home/mia/OneDrive/simsrg/envfixedfieldJ_md2dsweep_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  #FOR UNCOUPLED MODEL
    
    name_env=return_target_directory()+'envfixedfieldJ_md2dsweep_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  #FOR UNCOUPLED MODEL
    
    avalanche_runsim_fixedfieldJ(N0=N0, N=N, loop=loop, xmax=xmax, dt=dt, nstim=nstim,\
    percell=percell, \
    npprob=npprob, vj=vj,  sj=sj,\
    timeconst=timeconst,\
    phi=phi, eta=eta, epsilon=epsilon, name = name_env, num_pages = 5) 

    
def runsim_avalanche_env_fixedJh(N0, nstim, percell, time, phi, eta, epsilon, inputlabel='mixed',  N=1024, keyword = 'noreset', file_prefix='envfixedJh_sweep'):
    """
    Run a simulation. fields are fixed. couplings are fixed. Modified 23Nov2021 by AJS. 
    --------------------------------------------------------------
    Inputs:
    --------------------------------------------------------------
    N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
    N: number of cells to save
    nstim: number of latent fields
    percell: probability that each field couples to a given cell
    time: timeconstant of latent field
    phi: field multiplier
    eta: variance multiplier
    epsilon: bias
    inputlabel: label for sims filenames (time constant label)
    keyword: another label for sims filenames
    """
    N0 = N0 # number of cells

    N = N # number of cells after removing silent ones
 
    loop = 100 #1600 # number of track runs # was 40 in the first runs 

    xmax = 10000 # this MUST be set to 10000 or avalanche counting will not work
                 # note that this sets hard avalanche upper cutoff in durations

    dt=1. # increment for measurement locations

    # set up network structure

    nstim = nstim # number of nonplace stimuli

    percell= percell # probability that each field is accepted

    
    if nstim==0:
        npprob=np.array([1,0])
    else:
        npprob=np.array([1-percell, percell]) 
    # probability that cell is coupled to percell fields
    

    # now initialize distribution parameters of fields and couplings
    # here the place field couplings will be gamma distributed
    # the nonplace field couplings will be normally distributed

    vj = 0. # mean of couplings for nonplace fields
  
   

    sj = 1. # standard deviation of couplings for nonplace fields


    time=time
    if time == 'stationary':
        timelabel = 'stat'
        timeconst=time
    elif type(time)!=float:
        print('running with non-uniform time constants')
        timeconst=time
        timelabel=inputlabel
        print(timelabel)
        
        
    else: 
        timeconst = np.full((nstim,), time) # mean length of stochastic process in track lengths
        timelabel=time

    phi=phi # multiply by this constant to adjust overall activity of 
    # nonplace cells

    eta = eta # multiply by this constant to increase overall activity of 
    # network

    epsilon= epsilon


    #name_env='/home/mia/OneDrive/simsrg/env_av_stim{}e{}et{}ph{}p{}t{}'.format(nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel) #reset all every track run
    
    #name_env='/home/mia/OneDrive/simsrg/envfixedfieldJ_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon*eta, 1), np.round(eta,1), phi, percell, timelabel)  # DO NOT reset all every track run
#     name_env='/home/mia/OneDrive/simsrg/envfixedfieldJ_md2dsweep_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  #FOR UNCOUPLED MODEL
    
    name_env=return_target_directory()+ file_prefix + '_{}_av_stim{}e{}et{}ph{}p{}t{}'.format(keyword, nstim, np.round(epsilon, 2), np.round(eta,2), phi, percell, timelabel)  #FOR UNCOUPLED MODEL
    
    avalanche_runsim_fixedfieldJ(N0=N0, N=N, loop=loop, xmax=xmax, dt=dt, nstim=nstim,\
    percell=percell, \
    npprob=npprob, vj=vj,  sj=sj,\
    timeconst=timeconst,\
    phi=phi, eta=eta, epsilon=epsilon, name = name_env, num_pages = 5) 
