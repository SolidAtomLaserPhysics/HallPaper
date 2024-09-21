import numpy as np
import os
import math
import cmath
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

#for Antrag
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)



'''
Here you find a lot of plots I did with QMC data

This is how an inset works
    fig, ax1 = plt.subplots()

    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ax2 = plt.axes([0,0,1,1])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.3,0.4,0.4,0.35])          #left bottom width height
    ax2.set_axes_locator(ip)
'''



#TODO: nochmal uieberarbeiten mit error
def plotSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, version,points):
    dataSelfEnergy = np.loadtxt(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/self-en_wim_QMC_{}.csv'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, version))
    MatsubaraFreq = dataSelfEnergy[:, 0]
    RealPart = dataSelfEnergy[:, 1]
    ImagPart = dataSelfEnergy[:, 2]

    plt.plot(MatsubaraFreq[0:points], ImagPart[0:points], label = "U=1")                            #only plot until Matsubara frequency is the 20s value to see more structure
    plt.xlabel(r'$\nu$')
    plt.ylabel(r'Im($\Sigma$)')
    plt.legend()
    plt.grid()
    plt.savefig(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/selfEnergyPlot_QMC-U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}.png'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, u,beta,q,mu,t,tPri,tPriPri,kSteps,version))
    plt.clf()



def plotAllSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, Versions, points):
    for version in Versions:
        dataSelfEnergy = np.loadtxt(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/self-en_wim_QMC_{}.csv'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, version))
        MatsubaraFreq = dataSelfEnergy[:, 0]
        RealPart = dataSelfEnergy[:, 1]
        ImagPart = dataSelfEnergy[:, 2]
        RealError = dataSelfEnergy[:, 3]
        ImagError = dataSelfEnergy[:, 4]


        #plt.plot(MatsubaraFreq[0:points], ImagPart[0:points], label = 'Version {}'.format(version), marker = 'x', markersize = 0.2, linestyle = 'None')                            #only plot until Matsubara frequency is the 20s value to see more structure
        plt.errorbar(MatsubaraFreq[0:points], ImagPart[0:points], yerr = RealError[0:points], label = 'Version {}'.format(version), marker = 'x', linestyle = 'None')
        #plt.errorbar(MatsubaraFreq[0:points], ImagPart[0:points], yerr = 0.5, label = 'Version {}'.format(version), marker = 'x', linestyle = 'None')
        plt.xlabel(r'$\nu$')
        plt.ylabel(r'Im($\Sigma$)')
        plt.grid()
        plt.legend()
    plt.savefig(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/selfEnergyPlotAllVersions_QMC-points{}.png'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, points))
    plt.clf()
 





'''
plots the density n with mu on the x-axis
'''
#TODO: maybe error bars as well
def plotDensityQMC(QMCCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps,version):
    allDensities = np.zeros(len(Mu))
    for i_mu in range(len(Mu)):
        dataDensity = np.loadtxt(QMCCalculationDirectory + '''/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/densityValues_QMC_{}.csv'''.format(u,beta,q,Mu[i_mu],t,tPri,tPriPri,kSteps, version), dtype=float)
        allDensities[i_mu] = (dataDensity[0] + dataDensity[3])/2

    
    plt.plot(Mu, allDensities, marker = 'x', label = "Density n")                         
    plt.xlabel(r'$\mu$')
    plt.ylabel('n')
    plt.grid()
    plt.legend()
    #since not mu dependent anymore have this in the general directory, not in the calculation directory itself
    plt.savefig(QMCCalculationDirectory + '/density_QMC-U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}.png'.format(u,beta,q,t,tPri,tPriPri,kSteps,version))
    plt.clf()






'''
plots the density n of All Versions of QMC with mu on the x-axis
'''
#TODO: maybe error bars as well
def plotAllDensityQMC(QMCCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps,Versions):
    for version in Versions:
        allDensities = np.zeros(len(Mu))
        for i_mu in range(len(Mu)):
            dataDensity = np.loadtxt(QMCCalculationDirectory + '''/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/densityValues_QMC_{}.csv'''.format(u,beta,q,Mu[i_mu],t,tPri,tPriPri,kSteps, version), dtype=float)
            allDensities[i_mu] = (dataDensity[0] + dataDensity[3])/2

        
        plt.plot(Mu, allDensities, marker = 'x', label = 'Version {}'.format(version))                         
        plt.xlabel(r'$\mu$')
        plt.ylabel('n')
        plt.grid()
        plt.legend()
    #since not mu dependent anymore have this in the general directory, not in the calculation directory itself
    plt.savefig(QMCCalculationDirectory + '/density_QMC-U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}.png'.format(u,beta,q,t,tPri,tPriPri,kSteps))
    plt.clf()


'''
QMC and ED comparison for thesis
-------------------------------------------------
'''


def plotDensityQMCED(QMCCalculationDirectory, NevalinnaCalculationDirectory, EDFilling,QMCFilling, u,beta,q,Mu,t,tPri,tPriPri, kSteps, version,points, bathSites):
    plt.plot(Mu, EDFilling, marker = 'x', label = "Density n from ED")        
    plt.plot(Mu, QMCFilling, marker = 'o', label = "Density n from QMC")                       
    plt.xlabel(r'$\mu$')
    plt.ylabel('n')
    plt.grid()
    plt.legend()
    #since not mu dependent anymore have this in the general directory, not in the calculation directory itself
    plt.savefig(NevalinnaCalculationDirectory + '/density_QMC-ED-U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.png'.format(u,beta,q,t,tPri,tPriPri,kSteps,version, bathSites))
    plt.clf()




#plot self energy of correct ED, QMC and Julia code
def plotSigmaQMCED(QMCCalculationDirectory, NevalinnaCalculationDirectory, dataSelfEnergyED,dataSelfEnergyQMC, u,beta,q,mu,t,tPri,tPriPri, kSteps, version,points, bathSites):
    
    MatsubaraFreq = dataSelfEnergyQMC[:, 0]
    RealPart = dataSelfEnergyQMC[:, 1]
    ImagPart = dataSelfEnergyQMC[:, 2]
    plt.plot(MatsubaraFreq[0:points], ImagPart[0:points], label = "QMC", marker = "x")                            #only plot until Matsubara frequency is the 20s value to see more structure
    
    MatsubaraFreqED = dataSelfEnergyED[:, 0]
    RealPartED = dataSelfEnergyED[:, 1]
    ImagPartED = dataSelfEnergyED[:, 2]
    plt.plot(MatsubaraFreqED[0:points], ImagPartED[0:points], label = "ED", marker = "+")                            #only plot until Matsubara frequency is the 20s value to see more structure

    plt.ylabel(r'Im($\Sigma(i\nu_n)$)', fontsize=25, rotation=90)
    plt.xlabel(r'$\nu_{n}$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.legend(fontsize=15)
    plt.grid()
    plt.tight_layout()
    



    plt.savefig(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/selfEnergyComaprisonQMCED.pdf'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.clf()





'''
only testing
-------------------------------------------------
'''

def plotSigmaAllIterationsQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, Versions, points, Iterations):
    for version in Versions:
        for iteration in Iterations:
            dataSelfEnergy = np.loadtxt(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/self-en_wim_QMC_{}_iteration{}.csv'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, version, iteration))
            MatsubaraFreq = dataSelfEnergy[:, 0]
            RealPart = dataSelfEnergy[:, 1]
            ImagPart = dataSelfEnergy[:, 2]
            RealError = dataSelfEnergy[:, 3]
            ImagError = dataSelfEnergy[:, 4]


            #plt.plot(MatsubaraFreq[0:points], ImagPart[0:points], label = 'Iteration{}'.format(iteration), marker = 'x', markersize = 0.2, linestyle='None)                            #only plot until Matsubara frequency is the 20s value to see more structure
            plt.plot(MatsubaraFreq[0:points], ImagPart[0:points], label = 'Iteration{}'.format(iteration), marker = 'x')
            #plt.errorbar(MatsubaraFreq[0:points], ImagPart[0:points], yerr = RealError[0:points], label = 'Version {}'.format(version), marker = 'x', linestyle = 'None')
            #plt.errorbar(MatsubaraFreq[0:points], ImagPart[0:points], yerr = 0.5, label = 'Version {}'.format(version), marker = 'x', linestyle = 'None')
            plt.xlabel(r'$\nu$')
            plt.ylabel(r'Im($\Sigma$)')
            plt.legend()
        plt.savefig(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/selfEnergyPlotAllIterations_QMC-points{}_version{}.png'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, points,version))
        plt.clf()
 


'''
specific plot for the Antrag
ALSO: here you can find how to input an inset
'''
def plotAntragSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, version,points):
    fig, ax1 = plt.subplots()

    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ax2 = plt.axes([0,0,1,1])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.3,0.4,0.4,0.35])          #left bottom width height
    ax2.set_axes_locator(ip)

    u = 1.0
    mu = 0.5
    dataSelfEnergy = np.loadtxt(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/self-en_wim_QMC_{}.csv'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, version))
    MatsubaraFreq = dataSelfEnergy[:, 0]
    RealPart = dataSelfEnergy[:, 1]
    ImagPart = dataSelfEnergy[:, 2]
    ax1.plot(MatsubaraFreq[0:points], ImagPart[0:points], label = "U=1.0",  marker = 'x', markersize = 5)                            #only plot until Matsubara frequency is the 20s value to see more structure
    u = 2.0
    mu = 1.0
    dataSelfEnergy = np.loadtxt(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/self-en_wim_QMC_{}.csv'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, version))
    MatsubaraFreq = dataSelfEnergy[:, 0]
    RealPart = dataSelfEnergy[:, 1]
    ImagPart = dataSelfEnergy[:, 2]
    halfArray = int((len(MatsubaraFreq)/2))
    ax1.plot(MatsubaraFreq[0:points], ImagPart[0:points], label = "U=2.0", marker = 'x', markersize = 5) 
    ax1.set_xlabel(r'$\nu$', fontsize = 16)
    ax1.set_ylabel(r'Im($\Sigma$)', fontsize = 16)
    ax1.legend(fontsize="16")
    ax1.grid()




    u = 3.0
    mu = 1.5
    dataSelfEnergy = np.loadtxt(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/self-en_wim_QMC_{}.csv'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, version))
    MatsubaraFreq = dataSelfEnergy[:, 0]
    RealPart = dataSelfEnergy[:, 1]
    ImagPart = dataSelfEnergy[:, 2]
    ax2.plot(MatsubaraFreq[0:points, ImagPart[0:points]], label = "U=3.0", marker = 'x', markersize = 5, markerfacecolor='green', color='green') 
    ax2.set_xlabel(r'$\nu$', fontsize = 16)
    ax2.set_ylabel(r'Im($\Sigma$)', fontsize = 16)
    ax2.legend(fontsize="16")
    ax2.grid() 
    
    #TODO: lesbar
    fig.savefig(QMCCalculationDirectory + '/AntragSelfEnergiesCompared.pdf'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, u,beta,q,mu,t,tPri,tPriPri,kSteps,version))
    fig.clf()



