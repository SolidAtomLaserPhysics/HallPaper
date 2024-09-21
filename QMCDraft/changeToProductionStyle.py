import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


'''
changes the submit script to 5 Nodes (=96*5 = 480 cores) or so for the calculation
'''
def createSubmit(QMCDraftDirectory, QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,numberOfCores):
    with open(QMCDraftDirectory + '/submit_dmft_script', 'r', encoding='utf-8') as file:
        dataSubmit = file.readlines()
    
        dataSubmit[1] = "#SBATCH -t 9:00:00           \n"
        dataSubmit[2] = "#SBATCH -n {}   \n".format(numberOfCores)
        dataSubmit[4] = "#SBATCH -p standard96     \n"   
        dataSubmit[5] = "#SBATCH -A hhpnhytt \n"                                                   
        dataSubmit[12] = "mpirun -np {} /home/hhpnhytt/w2dynamics/DMFT.py Parameters.in      \n".format(numberOfCores)
    
    with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/submit_dmft_script'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps), 'w', encoding='utf-8') as file:
        file.writelines(dataSubmit)



def createParameters(u, mu, beta, q, t, tPri, tPriPri, kSteps, QMCDraftDirectory, QMCCalculationDirectory, nCorr):
    if q == 3:
        with open(QMCDraftDirectory + '/Parameters_Q3.in', 'r', encoding='utf-8') as file:
            data = file.readlines()
    
            data[2] = "NAt             = 3 #q des Magnetfeldes p/q      \n"
            data[5] = "beta            = {} #-> In Einheiten 1/t.     \n".format(beta)  
            data[7] = "mu              = {}         \n".format(mu) 
            data[9] = "readold         = -1 #-> Im ersten Schritt auskommentieren, sagt mir wo ich bei Abbruch weiterrechnen soll. -1 ist z.B. die letzte Iteration wenn die Rechnung durchgelaufen ist.   \n  "
            data[10] = "fileold         = U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT.hdf5 #-> Hier wird jede Iteration gespeichert. Von dort kann ich weiterrechnen.      \n".format(u, mu, beta, q, t, tPri, tPriPri, kSteps)
            data[11] = "DMFTsteps       = 10 #-> Anzahl der DMFT loops     \n"          #usually set it to 12
            data[14] = "FileNamePrefix  = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT   \n".format(u, mu, beta, q, t, tPri, tPriPri, kSteps)
            data[23] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[27] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[31] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[36] = "Nwarmups        = 1e7   \n" 
            data[37] = "Nmeas           = 1e7   \n"
            data[38] = "NCorr           = {} #-> Wie viele Schritte zwischen zwei Messungen. \n".format(nCorr)
        with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Parameters.in'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps), 'w', encoding='utf-8') as file:
            file.writelines(data)
    
    if q == 4:
        with open(QMCDraftDirectory + '/Parameters_Q4.in', 'r', encoding='utf-8') as file:
            data = file.readlines()
    
            data[2] = "NAt             = 4 #q des Magnetfeldes p/q      \n"
            data[5] = "beta            = {} #-> In Einheiten 1/t.     \n".format(beta)  
            data[7] = "mu              = {}         \n".format(mu) 
            data[9] = "readold         = -1 #-> Im ersten Schritt auskommentieren, sagt mir wo ich bei Abbruch weiterrechnen soll. -1 ist z.B. die letzte Iteration wenn die Rechnung durchgelaufen ist.    \n "
            data[10] = "fileold         = U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT.hdf5 #-> Hier wird jede Iteration gespeichert. Von dort kann ich weiterrechnen.      \n".format(u, mu, beta, q, t, tPri, tPriPri, kSteps)
            data[11] = "DMFTsteps       = 30 #-> Anzahl der DMFT loops     \n"
            data[14] = "FileNamePrefix  = U{}_mu{}_B{}_L{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT   \n".format(u, mu, beta, q, t, tPri, tPriPri, kSteps)
            data[23] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[27] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[31] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[35] = "Udd             = {} #-> Hubbard U     \n".format(u)
            data[36] = "Nwarmups        = 1e7   \n" 
            data[37] = "Nmeas           = 1e7   \n"
            data[42] = "NCorr           = {} #-> Wie viele Schritte zwischen zwei Messungen. \n".format(nCorr)
    
        with open(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Parameters.in'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps), 'w', encoding='utf-8') as file:
            file.writelines(data)









