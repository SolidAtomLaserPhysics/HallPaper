import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


'''
changes the submit script to 5 Nodes for the calculation
'''
def createDirectories(path):
    if not os.path.exists(path):                                  #make directory if not exists already
        os.mkdir(path)
    if not os.path.exists(path + "/hdf5BackupFolder"):
        os.mkdir(path + "/hdf5BackupFolder")


def createParameters(u, mu, beta, q, t, tPri, tPriPri, kSteps, QMCDraftDirectory, QMCCalculationDirectory):
    if q == 3:
        with open(QMCDraftDirectory + '/Parameters_Q3.in', 'r', encoding='utf-8') as file:
            data = file.readlines()
    
            data[2] = "NAt             = 3 #q des Magnetfeldes p/q      \n"
            data[5] = "beta            = {} #-> In Einheiten 1/t.     \n".format(beta)  
            data[7] = "mu              = {}         \n".format(mu) 
            data[10] = "#fileold         = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_kSteps{}DMFT.hdf5 #-> Hier wird jede Iteration gespeichert. Von dort kann ich weiterrechnen.      \n".format(u, mu, beta, q, t, tPri, tPriPri, kSteps)
            data[11] = "DMFTsteps       = 50 #-> Anzahl der DMFT loops     \n"
            data[14] = "FileNamePrefix  = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT   \n".format(u, mu, beta, q, t, tPri, tPriPri, kSteps)
            data[23] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[27] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[31] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[36] = "Nwarmups        = 1e7   \n" 
            data[37] = "Nmeas           = 1e6   \n"
            data[38] = "NCorr           = 30 #-> Wie viele Schritte zwischen zwei Messungen. \n"
        with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Parameters.in'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps), 'w', encoding='utf-8') as file:
            file.writelines(data)
    
    if q == 4:
        with open(QMCDraftDirectory + '/Parameters_Q4.in', 'r', encoding='utf-8') as file:
            data = file.readlines()
    
            data[2] = "NAt             = 4 #q des Magnetfeldes p/q      \n"
            data[5] = "beta            = {} #-> In Einheiten 1/t.     \n".format(beta)  
            data[7] = "mu              = {}         \n".format(mu) 
            data[10] = "#fileold         = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_DMFT.hdf5 #-> Hier wird jede Iteration gespeichert. Von dort kann ich weiterrechnen.      \n".format(u, mu, beta, q, t, tPri, tPriPri)
            data[11] = "DMFTsteps       = 100 #-> Anzahl der DMFT loops     \n"
            data[14] = "FileNamePrefix  = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_DMFT   \n".format(u, mu, beta, q, t, tPri, tPriPri)
            data[23] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[27] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[31] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[35] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[41] = "Nmeas           = 5e5   \n"
            data[42] = "NCorr           = 30 #-> Wie viele Schritte zwischen zwei Messungen. \n"
        with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Parameters.in'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps), 'w', encoding='utf-8') as file:
            file.writelines(data)

    if q == 9:
        with open(QMCDraftDirectory + '/Parameters_Q9.in', 'r', encoding='utf-8') as file:
            data = file.readlines()
    
            data[2] = "NAt             = 9 #q des Magnetfeldes p/q      \n"
            data[5] = "beta            = {} #-> In Einheiten 1/t.     \n".format(beta)  
            data[7] = "mu              = {}         \n".format(mu) 
            data[10] = "#fileold         = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_DMFT.hdf5 #-> Hier wird jede Iteration gespeichert. Von dort kann ich weiterrechnen.      \n".format(u, mu, beta, q, t, tPri, tPriPri)
            data[11] = "DMFTsteps       = 5 #-> Anzahl der DMFT loops     \n"
            data[14] = "FileNamePrefix  = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_DMFT   \n".format(u, mu, beta, q, t, tPri, tPriPri)
           # data[23] = "Udd             = {} #-> Hubbard U     \n".format(u)
            #data[27] = "Udd             = {} #-> Hubbard U     \n".format(u)
            #data[31] = "Udd             = {} #-> Hubbard U     \n".format(u)
            #data[35] = "Udd             = {} #-> Hubbard U     \n".format(u)
            #data[41] = "Nmeas           = 1e6   \n"
            #data[42] = "NCorr           = 30 #-> Wie viele Schritte zwischen zwei Messungen. \n"
        with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Parameters.in'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps), 'w', encoding='utf-8') as file:
            file.writelines(data)





    #now copy the FortranCode and change what to change
def createHamiltonian(q, QMCDraftDirectory, QMCCalculationDirectory,
                      beta, u, mu,t,tPri,tPriPri,kSteps):
    with open(QMCDraftDirectory + '/hofstadter_hamiltonian.f90', 'r', encoding='utf-8') as file:           #calculate ED results on /scratch
        dataFortran = file.readlines()

        dataFortran[5] = "  integer, parameter :: q={}   \n".format(q)

        dataFortran[6] = "  integer, parameter :: ksteps={}  \n".format(kSteps)

        dataFortran[8] = "  real(kind=8), parameter :: t={}d0   \n".format(t)
        dataFortran[9] = "  real(kind=8), parameter :: t1={}d0  \n".format(tPri)
        dataFortran[10] = "  real(kind=8), parameter :: t2={}d0  \n".format(tPriPri)
        
    with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/hofstadter_hamiltonian.f90'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps), 'w', encoding='utf-8') as file:
        file.writelines(dataFortran)




# can change running time and number of cores here
def createSubmit(QMCDraftDirectory, QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, numberOfCores):
    with open(QMCDraftDirectory + '/submit_dmft_script', 'r', encoding='utf-8') as file:
        dataSubmit = file.readlines()
    
        dataSubmit[1] = "#SBATCH -t 12:00:00           \n"
        dataSubmit[2] = "#SBATCH -n {}   \n".format(numberOfCores)
        dataSubmit[4] = "#SBATCH -p standard96     \n"                                                      
        dataSubmit[12] = "mpirun -np {} /home/hhpnhytt/w2dynamics/DMFT.py Parameters.in      \n".format(numberOfCores)
    
    with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/submit_dmft_script'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps), 'w', encoding='utf-8') as file:
        file.writelines(dataSubmit)




#just copy the run.sh
def createRun(QMCDraftDirectory, QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps):
    with open(QMCDraftDirectory + '/run.sh', 'r', encoding='utf-8') as file:
        dataRun = file.readlines()
    
    with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/run.sh'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps), 'w', encoding='utf-8') as file:
        file.writelines(dataRun)

    







