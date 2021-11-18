import numpy as np
import sys
import copy


try:
    TT = float(sys.argv[1])
except:
    TT = 10000.
    print("Time set to defaultvalue 20000")

try:
    tstep = float(sys.argv[2])
except:
    tstep = 0.5
    print("Timestep set to defaultvalue 1.")    

try:
    seed = float(sys.argv[3])
except:
    seed = 0.
    print("Seed set to defaultvalue 0")

try:
    deme = int(sys.argv[4])
except:
    deme = 200
    print("demes set to defaultvalue 200")

try:
    alph = float(sys.argv[5])
except:
    alph = 0.1
    print("alph set to defaultvalue 0.1")

try:
    rr = float(sys.argv[6])
except:
    rr = 0.3
    print("r set to defaultvalue 0.3")

try:
    S = int(sys.argv[7])
except:
    S = 100
    print("Species set to defaultvalue 100")

try:
    logdiff = float(sys.argv[8])
except:
    logdiff = -1.
    print("Log[lambda] set to defaultvalue -1")

try:
    kk = float(sys.argv[9])
except:
    kk = 10.
    print("K set to defaultvalue 10.")



def laplacian(Z):
    Zleft = Z[:,0:-2]
    Zright = Z[:,2:]
    Zcenter = Z[:,1:-1]
    return ( Zleft  + Zright - 2. * Zcenter)/2. 



def identical_competition_mat(species, alpha):
    mat  = alpha*np.ones((species,species),dtype=float)
    mat2 = mat + (1.-alpha)*np.identity(species) #- beta * np.diag(np.diag(mat))
    return mat2


def random_competition_mat_Gauss(species, alpha,alphaspreadio):
    matrand=np.random.normal(alpha, alphaspreadio/np.sqrt(species), size=(species,species))
    matzeroONdiag  = np.ones((species,species),dtype=float)-np.identity(species)
    mat2 = matzeroONdiag*matrand #+np.identity(species)
    return mat2



def func(state,r,rmat,alphamat,k):

    state = state * (rmat *(1.-np.matmul(np.identity(len(state)),state)/k)-r*np.matmul(alphamat,state)/k)
    
    return state



def update(species,demes,f,state,diff,dt,args=()):
        
        #For Euler forward
        k1=f(state[:,1:-1],*args) 
        
        #For 4th order RungeKutta (RK4)
        #k2=f(state[:,1:-1]+k1*dt/2.,*args)
        #k3=f(state[:,1:-1]+k2*dt/2.,*args)
        #k4=f(state[:,1:-1]+k3*dt,*args)

        #Demographic fluctuations
        noise = np.sqrt(dt)*(np.random.poisson(state[:,1:-1])-state[:,1:-1])
        
        #Calculate contributions from interactions
        #For Euler:
        LV = dt * k1 
        #For RK4:
        #LV=(dt/6.)* (k1+2.*k2+2*k3+k4)         
        

        #Add interactions:
        state[:,1:-1] = state[:,1:-1] + LV
        
        if(np.min(state[:,1:-1])<0.):
            print("negative density from interaction contributions!")
            state[:,1:-1]=-1.

        else: 

            #periodic boundary conditions:
            for i in range(len(state)): 
                state[i,0] = state[i,-2]
                state[i,-1] = state[i,1]

            
            #Short-range migration
            migration = diff * laplacian(state) #Diffusive migration
        
            #... or global migration
            #migration = diff * np.transpose(state[:,1:-1].mean(axis=1)-np.transpose(state[:,1:-1])) #Global migration
            
            #Add migration:
            state[:,1:-1]= state[:,1:-1] + dt * migration 

            if(np.min(state)<0.):
                print("negative density from migration!")
                state[:,1:-1]=-1.


            #Add noise:
            state[:,1:-1]= state[:,1:-1] + noise


        return state



def CoupledCompetitiveLV(T,species,demes,r,rstd,alpha,alphastd,LogDiff,Diffstd,seed,k,timestep):

    Nall=open('N_%f_%f_%f.csv' % (species,LogDiff,alpha),'wb')

    Diff=10.**LogDiff
    

    print("Log[diff] is %f, nr of species is %f, and the interaction strength is %f" % (LogDiff,species,alpha))


    np.random.seed(int(seed))

    
    #Initialize interaction matrix
    mat = random_competition_mat_Gauss(species,alpha,alphastd)
    
    #Initialize fitnesses
    rvec= np.random.normal(r, rstd, size=(species)) #Gauss
    rmat= np.transpose(np.tile(rvec,(demes,1)))
    
    #Initialize migration rates
    Diffvec=np.random.normal(Diff,Diffstd, size=(species))
    Diffvec[Diffvec < 0.] = 0.
    Diffmat=np.transpose(np.tile(Diffvec,(demes,1)))
    

    #computation parameters
    dt=timestep
    dtoriginal=copy.deepcopy(dt)
    savespace=1                                  #How many demes are getting saved
    savespecies=1                                #How many species are getting saved
    savesteps=50.                                #How many timesteps are getting saved
    Minimaltime=1.                            #First saved time point
    samplingtime=(T-Minimaltime)/savesteps
    totaltime=0.
    it=0.
    
    

    #initial state
    perturbation=0.05
    state = k*(1.+perturbation*(2.*np.random.random(size=(species,demes+2))-1.))
    state = np.round(state)
    

    #Periodic boundary conditions
    for i in range(len(state)):
            state[i,0] = state[i,-2]
            state[i,-1] = state[i,1]

    #Show mean and save current abundances
    print(np.mean(state[:,1:-1:savespace]/k))
    np.savetxt(Nall,state[::savespecies,1:-1:savespace]/k, delimiter=",")
    

    while (totaltime < T):
            
            if (np.max(state)>10.**(-15.)):
                
                statec=copy.deepcopy(state)
                state = update(species,demes,func,state,Diffmat,dt,args=(r,rmat,mat,k))

                if(np.min(state)<0.):
                    dt=0.5*dt           #Half the time step
                    print('Violation: densities negative! Adapted timestep to ' + repr(dt))
                    state = statec
                    if dt<10.**-8.:
                        print("not converged")
                        break    


                else:

                    totaltime=totaltime+dt 
                        
                    if (totaltime >= Minimaltime and totaltime >= (Minimaltime+samplingtime*it)):
                        print('time=' + repr(totaltime))
                        
                        #Show mean and save current abundances
                        np.savetxt(Nall,state[::savespecies,1:-1:savespace]/k, delimiter=",")
                        print(np.mean(state[:,1:-1:savespace]/k))
                        dt=dtoriginal
                        it=it+1

            else:
                
                print(totaltime)
                print("extinct")
                totaltime=T+dt 


                    
    print('time=' + repr(totaltime))
    

    #Show mean and save current abundances
    np.savetxt(Nall,state[::savespecies,1:-1:savespace]/k, delimiter=",")
    print(np.mean(state[:,1:-1:savespace]/k))
    #print(state/k)
    Nall.close()
    



rstd=0.
alphstd=0.
diffstd=0.

CoupledCompetitiveLV(TT,S,deme,rr,rstd,alph,alphstd,logdiff,diffstd,seed,kk,tstep)








