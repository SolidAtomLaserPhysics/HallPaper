import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


'''
calls a Julia script in the Draft folder,
that reads out the hdf5 file of each calculation and then prints out the self energies
MatFreq     Re(Sigma)       Imag(Sigma)         sigma Re(Sigma)         sigma Imag(Sigma)
'''
def extractSelfEnergy(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, version):
    os.chdir(QMCDraftDirectory)                  
    os.system('julia readQMCSelfEnergyResults.jl {} {} {} {} {} {} {} {} {} {}'.format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version, QMCCalculationDirectory))


'''
basically the same, but here we only look out for the different Iterations(Array[n]) to compare convergence
'''
def extractSelfEnergyIterations(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, version, iteration):
    os.chdir(QMCDraftDirectory)  
    os.system('julia readQMCSelfEnergyResultsIterations.jl {} {} {} {} {} {} {} {} {} {} {}'.format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version, QMCCalculationDirectory, iteration))
        


'''
calls a Julia script in the Draft folder,
that reads out the hdf5 file of each calculation and then prints out density
'''
def extractDensity(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, version):
    os.chdir(QMCDraftDirectory)                      
    os.system('julia readQMCDensityResults.jl {} {} {} {} {} {} {} {} {} {}'.format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version, QMCCalculationDirectory))


'''
looks for number of done iterations in slurm
'''
def extractNumberOfIter(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, version):
    pass
    #TODO: do it

