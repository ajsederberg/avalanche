#avafunc.py

# import stuff
from placerg.funcs import *
from placerg.funcsrg import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from placerg.funcsall import *
from scipy.optimize import least_squares
from scipy.sparse import csr_matrix
from scipy.io import savemat

# 
def calc_avalanche_inf(act):
    """
    calculate avalanche sizes and durations for infinite time constant simulations.
    -------------------------------------------------------------------------------------------
    INPUTS:
    -----------------------------------------------------------------------------------------
    act: sparse activity array of shape (number of cells, track length*number of track runs)
    -----------------------------------------------------------------------------------------------
    Results:
    -----------------------------------------------------------------------------------------------
    sizes and durations of avalanches
    """
    theta = 0 #0.5*np.median(np.sum(env.pmat, axis=0))
    #pmat_tot = np.sum(act, axis=0)
    pmat_tot = np.array(act.sum(axis=0))[0]
    pmat_tot_bool = pmat_tot != theta
    avalanches = [] # location of avalanche start, duration of avalanche
    min_avalanche = 1
    state = False
    a = 0 # length of current avalanche
    counter = 0 #track length counter
    track_len = 10000    # this is why loop must be set to 10000; to fix this magic number issue,
                         # this function needs to take loop as an input. 
    for i in range(len(pmat_tot_bool)):
        if counter <= track_len-1:
            if state == True: # avalanche state
                if pmat_tot_bool[i] == True: # if avalanche continues
                    a += 1
                    
                else: # if avalanche ends
                    if a >= min_avalanche:
                        avalanches.append((i - a, a))
                    a = 0
                    state = False
            else:
                if pmat_tot_bool[i] == True: # avalanche starts
                    state = True
                    a += 1
                else:
                    pass
            counter += 1
        else:
            counter = 0
            state = False
            a = 0
            
            if pmat_tot_bool[i] == True: # avalanche starts
                state = True
                a += 1
            else:
                pass
            counter += 1
        
    # do integral
    area = []
    for avalanche in avalanches:
        sum_area = np.sum(pmat_tot[avalanche[0]:avalanche[0]+avalanche[1]])
        area.append(sum_area)

    durations = np.stack(avalanches)[:,1]
    sizes = area
    
    return sizes, durations

def calc_av_means(sizes, durations):
    """
    calculate mean size per duration of avalanches
    -------------------------------------------------------------------------------------------
    INPUTS:
    -----------------------------------------------------------------------------------------
    sizes: sizes of avalches
    durations: corresponding durations of avalanches
    -----------------------------------------------------------------------------------------------
    Results:
    -----------------------------------------------------------------------------------------------
    mean size and corresponding duration of avalanches
    """
    crits = np.stack((durations, sizes)).T

    duration_mean = []
    size_mean = []
    uniques = np.unique(crits[:,0])
    for unique in uniques:
        wh = np.where(crits[:,0]==unique)
        mean = np.mean(crits[:,1][wh])
        duration_mean.append(unique)
        size_mean.append(mean)
    duration_mean = np.array(duration_mean)
    size_mean = np.array(size_mean)

    return size_mean, duration_mean

def calc_avalanche(act):
    """
    calculate avalanche sizes and durations for finite time constant simulations.
    -------------------------------------------------------------------------------------------
    INPUTS:
    -----------------------------------------------------------------------------------------
    act: sparse activity array of shape (number of cells, track length*number of track runs)
    -----------------------------------------------------------------------------------------------
    Results:
    -----------------------------------------------------------------------------------------------
    sizes and durations of avalanches
    """
    theta = 0 #0.5*np.median(np.sum(env.pmat, axis=0))
    pmat_tot = np.array(act.sum(axis=0))[0]
    pmat_tot_bool = pmat_tot != theta
    avalanches = [] # location of avalanche start, duration of avalanche
    min_avalanche = 1
    state = False
    a = 0 # length of current avalanche
    #counter = 0 #track length counter
    #track_len = 1600
    for i in range(len(pmat_tot_bool)):
        if state == True: # avalanche state
            if pmat_tot_bool[i] == True: # if avalanche continues
                a += 1
            else: # if avalanche ends
                if a >= min_avalanche:
                    avalanches.append((i - a, a))
                a = 0
                state = False
        else:
            if pmat_tot_bool[i] == True: # avalanche starts
                state = True
                a += 1
            else:
                pass
        #counter += 1
            
    # do integral
    area = []
    for avalanche in avalanches:
        sum_area = np.sum(pmat_tot[avalanche[0]:avalanche[0]+avalanche[1]])
        area.append(sum_area)

    durations = np.stack(avalanches)[:,1]
    sizes = area

    return sizes, durations

def save_avalanches(file_name, data_dir, save_dir):
    """
    Extracts avalanches from data in envname and saves to external file
    --------------------------------------------------------------------------------------------------------
    Arguments:
    --------------------------------------------------------------------------------------------------------
    envname: name of simulation files
     --------------------------------------------------------------------------------------------------------
    Returns:
    --------------------------------------------------------------------------------------------------------
    
    """
    sizes = []
    durations = []
    size_means = []
    duration_means = []
    act = []
    
    # data is in directory data_dir
    envname = data_dir + file_name
    for i in range(5):
        pmat = sparse.load_npz(envname+'_'+str(i)+'.npz')
        act.append(pmat)
        del pmat

    act = sparse.hstack([act[i] for i in range(len(act))]).tocsr()
    print('converted to csr')
    sizes, durations = calc_avalanche_inf(act)
    print(sizes[0:5])
    print(durations[0:5])
    
    mdic = {"sizes": sizes, "durations": durations}
    save_filename = save_dir+file_name+'.mat'
    print(save_filename)
    savemat(save_filename, mdic)
    
def save_sim_stats(file_name, data_dir, save_dir):
    """
    Extracts avalanches from data in envname and saves to external file
    --------------------------------------------------------------------------------------------------------
    Arguments:
    --------------------------------------------------------------------------------------------------------
    envname: name of simulation files
     --------------------------------------------------------------------------------------------------------
    Returns:
    --------------------------------------------------------------------------------------------------------
    
    """
#     sizes = []
#     durations = []
#     size_means = []
#     duration_means = []
    act = []
    
    # data is in directory data_dir
    envname = data_dir + file_name
    for i in range(5):
        pmat = sparse.load_npz(envname+'_'+str(i)+'.npz')
        act.append(pmat)
        del pmat

    act = sparse.hstack([act[i] for i in range(len(act))]).tocsr()
    print('converted to csr')
#     sizes, durations = calc_avalanche_inf(act)
#     print(sizes[0:5])
#     print(durations[0:5])

#     FR_mean = np.mean(act)
    #print('FR_mean = '+str(FR_mean))
#     act_var = np.var(act.toarray())
    #print('act. var.  = '+str(act_var))
    cell_FR = np.mean(act, axis=1)
#     rate_var = np.var(rates)

    savefile = envname    
    sf = h5py.File(savefile+'.h5', 'r+')
    sg = sf.get('dataset_0')
    J = np.array(sg['J'])
    epsilon  = np.array(sg['epsilon'])
    eta  = np.array(sg['eta'])
    xmax = np.array(sg['xmax'])
    loop = np.array(sg['loop'])
    dt = np.array(sg['dt'])
    
    print(sg.keys())

    
    numpages = 5  # magic number, hard-coded. this is why processes go 0 to 4. 
    process_0 = np.array(sg['process_0'])
    process_1 = np.array(sg['process_1'])
    process_2 = np.array(sg['process_2'])
    process_3 = np.array(sg['process_3'])
    process_4 = np.array(sg['process_4'])
    print(process_0.shape)
    sf.close()
    
    # fields is HUGE, downsample to one per loop
    fields = np.vstack([process_0, process_1, process_2, process_3, process_4]).T
    num_field_samples = fields.size
    print(xmax)
    print(loop)
    print(dt)
    field_samples = fields[:, np.arange(0, int(loop*xmax*dt*numpages),int(xmax))]
    
    mdic = {"J": J, "epsilon": epsilon, "eta": eta, "cell_FR": cell_FR, "xmax": xmax, "loop": loop,\
            "num_field_samples": num_field_samples, "field_samples": field_samples}
    
    save_filename = save_dir+file_name+'fr_stats.mat'
    print(save_filename)
    savemat(save_filename, mdic)


#     fig, ax = plt.subplots(1,3,figsize=(2*27/4,27/(4*1.5)))
    
#     plt.subplot(1, 3, 1)
#     cts2, bins2 = np.histogram(fields)
#     plt.hist(bins2[:-1], bins2, weights=cts2)
    
#     plt.subplot(1, 3, 3)
#     counts, bins = np.histogram(cell_FR)
#     plt.hist(bins[:-1], bins, weights=counts)