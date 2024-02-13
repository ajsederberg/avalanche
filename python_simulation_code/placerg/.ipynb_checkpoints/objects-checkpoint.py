#objects.py

import numpy as np
import copy
from placerg.funcs import *
#from placerg.funcsrg import *
import placerg._dice6 as _dice6
from scipy import sparse
import h5py
import os


def avalanche_runsim(N0, N, loop, xmax, dt, nstim, percell, \
                  npprob, vj,  sj, \
                  timeconst, phi, eta, epsilon, name, num_pages):
        """
        Run a simulation. Couplings and fields are selected randomly. Use functions from runfunc.py to call this function and run sims
        --------------------------------------------------------------
        Inputs:
        --------------------------------------------------------------
        N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
        N: number of cells in saved simulation
        loop: number of track loops
        xmax: track length
        dt: time step (1)
        nstim: number of latent fields
        percell: probability that each field couples to a given cell
        npprob: [1-percell, percell] 
        vj: latent field mean
        sj: latent field standard deviation
        timeconst: timeconstant of latent field
        phi: field multiplier
        eta: variance multiplier
        epsilon: bias
        name: name for sims
        num_pages: number of memory pages
        """
                  
        savefile = name
        sf = h5py.File(savefile+'.h5', 'w')
        sg = sf.create_group('dataset_0')
        
        sg.create_dataset('N', data = N)
        sg.create_dataset('loop', data=loop)
        sg.create_dataset('xmax', data = xmax)
        sg.create_dataset('dt', data = dt)
        sg.create_dataset('nstim', data=nstim)
        sg.create_dataset('percell', data = percell)
        sg.create_dataset('vj', data = vj)
        sg.create_dataset('sj', data = sj)
        sg.create_dataset('phi', data = phi)
        sg.create_dataset('eta', data = eta)
        sg.create_dataset('epsilon', data = epsilon)
        
        
        
        choice=np.array([0,1])
        
        sigmas= np.full((nstim,), 1.) 
                            # standard deviation of stochastic process. 
                            # shape: (nstim,)
                            
        J=fillJ(np.zeros((nstim, N0)), N0,\
                vj,\
                sj, nstim, npprob,\
                choice, phi) 
                
                
        mean = 0        
        for i in range(num_pages):             
            if timeconst == 'stationary':
                process= stim_stationary(nstim, loop, int(loop*xmax*dt)).T 
            else:                  
                taus= timeconst*xmax                   
                if i==0:
                    sg.create_dataset('timeconst', data = timeconst)
                    process= stim(taus, sigmas, dt, \
                        int(loop*xmax*dt), xmax).T 
                else:
                    # edit (12/16/21), remove commented-out if you want each track run to start from previous state
                    # process= stim(taus, sigmas, dt, \
                        # int(loop*xmax*dt), xmax, state_input = prev_state).T 
                    process= stim(taus, sigmas, dt, \
                        int(loop*xmax*dt), xmax).T 
                    
                    
                            # this makes an array of shape 
                            # (loop*dt*xmax, nstim)
                            
            sg.create_dataset('process_'+str(i), data = process)
                    
#                             # computeh computes the Hamiltonian H = eta*(Jh + epsilon)
#             pmat=np.vstack(_dice6.dice6(computeP(computeh(fillfields(N0, \
#                 np.zeros((loop, xmax, nstim)),\
#                 process, loop) , J, eta, epsilon) ))).T
                              # computeh_uncoupled computes the Hamiltonian H = eta*Jh + epsilon
                              # this is more commonly used starting late 2021
            pmat=np.vstack(_dice6.dice6(computeP(computeh_uncoupled(fillfields(N0, \
                np.zeros((loop, xmax, nstim)),\
                process, loop) , J, eta, epsilon) ))).T     
    
            prev_state = process[-1,:]
               
            del process
       
            mean_i = np.mean(pmat, axis=1)
            print(savefile)
                      
            sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat))
        
            del pmat
            
            mean += mean_i
        
        inds=nonzerocell(mean, N)
        
        for i in range(num_pages):
            pmat = sparse.load_npz(savefile+'_'+str(i)+'.npz')
            pmat = pmat.toarray()[inds, :]
            sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat))       
        
        del pmat
        
        J=J[:, inds]
        sg.create_dataset('J', data = J)
        sf.close()
                  
        
def avalanche_runsim_fixedJ(N0, N, loop, xmax, dt, nstim, percell, \
                  npprob, vj,  sj, \
                  timeconst, phi, eta, epsilon, name, num_pages):
        """
        Run a simulation. fields are selected randomly. Couplings are fixed. Use functions from runfunc.py to call this function and run sims
        --------------------------------------------------------------
        Inputs:
        --------------------------------------------------------------
        N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
        N: number of cells in saved simulation
        loop: number of track loops
        xmax: track length
        dt: time step (1)
        nstim: number of latent fields
        percell: probability that each field couples to a given cell
        npprob: [1-percell, percell] 
        vj: latent field mean
        sj: latent field standard deviation
        timeconst: timeconstant of latent field
        phi: field multiplier
        eta: variance multiplier
        epsilon: bias
        name: name for sims
        num_pages: number of memory pages
        """
                  
        savefile = name
        sf = h5py.File(savefile+'.h5', 'w')
        sg = sf.create_group('dataset_0')
        
        sg.create_dataset('N', data = N)
        sg.create_dataset('loop', data=loop)
        sg.create_dataset('xmax', data = xmax)
        sg.create_dataset('dt', data = dt)
        sg.create_dataset('nstim', data=nstim)
        sg.create_dataset('percell', data = percell)
        sg.create_dataset('vj', data = vj)
        sg.create_dataset('sj', data = sj)
        sg.create_dataset('phi', data = phi)
        sg.create_dataset('eta', data = eta)
        sg.create_dataset('epsilon', data = epsilon)
        
        
        
        choice=np.array([0,1])
        
        sigmas= np.full((nstim,), 1.) 
                            # standard deviation of stochastic process. 
                            # shape: (nstim,)
                            
        J=fillJ_fixed(np.zeros((nstim, N0)), N0,\
                vj,\
                sj, nstim, npprob,\
                choice, phi) 
                
                
        mean = 0        
        for i in range(num_pages):             
            if timeconst == 'stationary':
                process= stim_stationary(nstim, loop, int(loop*xmax*dt)).T 
            else:                  
                taus= timeconst*xmax                   
                if i==0:
                    sg.create_dataset('timeconst', data = timeconst)
                    process= stim(taus, sigmas, dt, \
                        int(loop*xmax*dt), xmax).T 
                else:
                    process= stim(taus, sigmas, dt, \
                        int(loop*xmax*dt), xmax, state_input = prev_state).T 
                    
                    
                            # this makes an array of shape 
                            # (loop*dt*xmax, nstim)
                            
            sg.create_dataset('process_'+str(i), data = process)
                    
                            # computeh_uncoupled computes the Hamiltonian H = eta*Jh + epsilon
            pmat=np.vstack(_dice6.dice6(computeP(computeh_uncoupled(fillfields(N0, \
                np.zeros((loop, xmax, nstim)),\
                process, loop) , J, eta, epsilon) ))).T
            
            prev_state = process[-1,:]
               
            del process
       
            mean_i = np.mean(pmat, axis=1)
                      
            sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat))
        
            del pmat
            
            mean += mean_i
        
        inds=nonzerocell_fixedJ(mean, N)
        
        for i in range(num_pages):
            pmat = sparse.load_npz(savefile+'_'+str(i)+'.npz')
            pmat = pmat.toarray()[inds, :]
            sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat))       
        
        del pmat
        
        J=J[:, inds]
        sg.create_dataset('J', data = J)
        sf.close()
        
        
def avalanche_runsim_fixedfield(N0, N, loop, xmax, dt, nstim, percell, \
                  npprob, vj,  sj, \
                  timeconst, phi, eta, epsilon, name, num_pages):
        """
        Run a simulation. couplings are selected randomly. fields are fixed. Use functions from runfunc.py to call this function and run sims
        --------------------------------------------------------------
        Inputs:
        --------------------------------------------------------------
        N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
        N: number of cells in saved simulation
        loop: number of track loops
        xmax: track length
        dt: time step (1)
        nstim: number of latent fields
        percell: probability that each field couples to a given cell
        npprob: [1-percell, percell] 
        vj: latent field mean
        sj: latent field standard deviation
        timeconst: timeconstant of latent field
        phi: field multiplier
        eta: variance multiplier
        epsilon: bias
        name: name for sims
        num_pages: number of memory pages
        """
                  
        savefile = name
        sf = h5py.File(savefile+'.h5', 'w')
        sg = sf.create_group('dataset_0')
        
        sg.create_dataset('N', data = N)
        sg.create_dataset('loop', data=loop)
        sg.create_dataset('xmax', data = xmax)
        sg.create_dataset('dt', data = dt)
        sg.create_dataset('nstim', data=nstim)
        sg.create_dataset('percell', data = percell)
        sg.create_dataset('vj', data = vj)
        sg.create_dataset('sj', data = sj)
        sg.create_dataset('phi', data = phi)
        sg.create_dataset('eta', data = eta)
        sg.create_dataset('epsilon', data = epsilon)
        
        
        
        choice=np.array([0,1])
        
        sigmas= np.full((nstim,), 1.) 
                            # standard deviation of stochastic process. 
                            # shape: (nstim,)
                            
        J=fillJ(np.zeros((nstim, N0)), N0,\
                vj,\
                sj, nstim, npprob,\
                choice, phi) 
                
                
        mean = 0        
        for i in range(num_pages):             
            if timeconst == 'stationary':
                process= stim_stationary_fixedfield(nstim, loop, int(loop*xmax*dt), seed = i).T 
            else:     
                raise NotImplementedError('dynamic fields not implemented')              
                #taus= timeconst*xmax                   
                #if i==0:
                #    sg.create_dataset('timeconst', data = timeconst)
                #    process= stim(taus, sigmas, dt, \
                #        int(loop*xmax*dt), xmax).T 
                #else:
                #    process= stim(taus, sigmas, dt, \
                #        int(loop*xmax*dt), xmax, state_input = prev_state).T 
                    
                    
                            # this makes an array of shape 
                            # (loop*dt*xmax, nstim)
                            
            sg.create_dataset('process_'+str(i), data = process)
                    
        
            pmat=np.vstack(_dice6.dice6(computeP(computeh_uncoupled(fillfields(N0, \
                np.zeros((loop, xmax, nstim)),\
                process, loop) , J, eta, epsilon) ))).T
            
            prev_state = process[-1,:]
               
            del process
       
            mean_i = np.mean(pmat, axis=1)
                      
            sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat))
        
            del pmat
            
            mean += mean_i
        try:
            inds=nonzerocell(mean, N)
            for i in range(num_pages):
                pmat = sparse.load_npz(savefile+'_'+str(i)+'.npz')
                pmat = pmat.toarray()[inds, :]
                sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat)) 
            del pmat
        
            J=J[:, inds]
            sg.create_dataset('J', data = J)
            sf.close()
        
        except:
            print('Not enough active cells: sim deleted')
            sf.close()
            os.remove(savefile+'.h5')
            for i in range(num_pages): 
                os.remove(savefile+'_'+str(i)+'.npz')
            
            
            
        
    
        
        
        
def avalanche_runsim_fixedfieldJ(N0, N, loop, xmax, dt, nstim, percell, \
                  npprob, vj,  sj, \
                  timeconst, phi, eta, epsilon, name, num_pages):
        """
        Run a simulation. fields are fixed. Couplings are fixed. Use functions from runfunc.py to call this function and run sims
        --------------------------------------------------------------
        Inputs:
        --------------------------------------------------------------
        N0: number of initial cells to simulate. we simulate N0 cells, then remove the silent cells and select N cells to save. recommended 4 times the numer opf desired cells, N
        N: number of cells in saved simulation
        loop: number of track loops
        xmax: track length
        dt: time step (1)
        nstim: number of latent fields
        percell: probability that each field couples to a given cell
        npprob: [1-percell, percell] 
        vj: latent field mean
        sj: latent field standard deviation
        timeconst: timeconstant of latent field
        phi: field multiplier
        eta: variance multiplier
        epsilon: bias
        name: name for sims
        num_pages: number of memory pages
        """
                  
        savefile = name
        sf = h5py.File(savefile+'.h5', 'w')
        sg = sf.create_group('dataset_0')
        
        sg.create_dataset('N', data = N)
        sg.create_dataset('loop', data=loop)
        sg.create_dataset('xmax', data = xmax)
        sg.create_dataset('dt', data = dt)
        sg.create_dataset('nstim', data=nstim)
        sg.create_dataset('percell', data = percell)
        sg.create_dataset('vj', data = vj)
        sg.create_dataset('sj', data = sj)
        sg.create_dataset('phi', data = phi)
        sg.create_dataset('eta', data = eta)
        sg.create_dataset('epsilon', data = epsilon)
        
        
        
        choice=np.array([0,1])
        
        sigmas= np.full((nstim,), 1.) 
                            # standard deviation of stochastic process. 
                            # shape: (nstim,)
                            
        J=fillJ_fixed(np.zeros((nstim, N0)), N0,\
                vj,\
                sj, nstim, npprob,\
                choice, phi) 
                
                
        mean = 0        
        for i in range(num_pages):             
            if timeconst == 'stationary':
                process= stim_stationary_fixedfield(nstim, loop, int(loop*xmax*dt), seed = i).T 
            else:       
                raise NotImplementedError('dynamic fields not implemented')            
                #taus= timeconst*xmax                   
                #if i==0:
                #    sg.create_dataset('timeconst', data = timeconst)
                #    process= stim(taus, sigmas, dt, \
                #        int(loop*xmax*dt), xmax).T 
                #else:
                #    process= stim(taus, sigmas, dt, \
                #        int(loop*xmax*dt), xmax, state_input = prev_state).T 
                    
                    
                            # this makes an array of shape 
                            # (loop*dt*xmax, nstim)
                            
            sg.create_dataset('process_'+str(i), data = process)
                    
        
            pmat=np.vstack(_dice6.dice6(computeP(computeh_uncoupled(fillfields(N0, \
                np.zeros((loop, xmax, nstim)),\
                process, loop) , J, eta, epsilon) ))).T
            
            prev_state = process[-1,:]
               
            del process
       
            mean_i = np.mean(pmat, axis=1)
                      
            sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat))
        
            del pmat
            
            mean += mean_i
        
        try:
            inds=nonzerocell_fixedJ(mean, N)
            for i in range(num_pages):
                pmat = sparse.load_npz(savefile+'_'+str(i)+'.npz')
                pmat = pmat.toarray()[inds, :]
                sparse.save_npz(savefile+'_'+str(i)+'.npz', sparse.csr_matrix(pmat)) 
            del pmat
        
            J=J[:, inds]
            sg.create_dataset('J', data = J)
            sf.close()
        
        except:
            print('Not enough active cells: sim deleted')
            sf.close()
            os.remove(savefile+'.h5')
            for i in range(num_pages): 
                os.remove(savefile+'_'+str(i)+'.npz')
