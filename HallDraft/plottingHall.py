import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess
import math
import cmath

from scipy import linalg


def getHallValueDirect(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallDirect/runHall.out'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'r', encoding='utf-8') as file:
        data = file.readlines()
        Hall = float(data[9])       #value always in third row of run.out

        return Hall

          

def plotHallDirect(NevalinnaCalculationDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, HallArray):
    x = Mu[:]       #all mu on x-axis
    y = HallArray[:]       
    plt.plot(x, y, label = "Hall from gm_wre", marker = 'x')     

    plt.ylabel(r'$\sigma_{xy}$')
    plt.xlabel(r'$\mu$')
    plt.legend()
    plt.savefig(NevalinnaCalculationDirectory + '/HallResultDirect_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.png'.format(u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.clf()




def getHallValueEta(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, eta):
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}/runHall.out'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), 'r', encoding='utf-8') as file:
        data = file.readlines()
        Hall = float(data[len(data) - 1])       #value always last row of run.out

    return Hall


def plotHallEta(NevalinnaCalculationDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, HallArray, eta):

    x = Mu[:]       #all mu on x-axis
    y = HallArray[:]       
    plt.plot(x, y, label = "$\eta =$ {}".format(eta), marker = 'x')     

    plt.ylabel(r'$\sigma_{xy}$', fontsize=25, rotation=90)
    plt.xlabel(r'$\mu$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.legend(fontsize=15)
    plt.grid()
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/HallResultEta{}_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.pdf'.format(eta, u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.savefig(NevalinnaCalculationDirectory + '/HallResultEta{}_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.png'.format(eta, u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.clf()




def plotHallAllEta(NevalinnaCalculationDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, HallArray, Eta):
    #x = Mu[:]       #all mu on x-axis
    #y = HallArray[:]       
    #plt.plot(x, y, label = "Hall from gm_wre", marker = 'x')   
    
    HallArrayEta = np.zeros(len(Mu))
    for eta in Eta:
        #load data
        for m, mu in enumerate(Mu):
            HallArrayEta[m] = getHallValueEta(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, eta)          
        x = Mu[:]       #all mu on x-axis
        y = HallArrayEta[:]   
        if (u == 0):    
            plt.plot(x, y, label = "$\delta =$ {}".format(eta), marker = 'x') 
        else:
            plt.plot(x, y, label = "$\eta =$ {}".format(eta), marker = 'x')   
            
    #plt.axvline(x = 0.65, label = "half filling", linestyle= "dashed")

    plt.ylabel(r'$\sigma_{xy}\left[\frac{e^2}{2\pi\hbar}\right]$', fontsize=25, rotation=90)
    plt.xlabel("$\mu$", fontsize=25, rotation=0)
    plt.ylim = (0,4.5)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.legend(fontsize=14, loc='upper left', bbox_to_anchor=(-0.03, 1.15), ncol=2, fancybox=False, shadow=False)  #eig 1.35 
    #plt.legend(fontsize = 15)
    plt.grid()
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/HallResultAllEta_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.pdf'.format(u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.savefig(NevalinnaCalculationDirectory + '/HallResultAllEta_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.png'.format(u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.clf()





#-----------------------------------------------------------------------------------------------------------
#last plots for colloquium
def plotHallEtaColloq(NevalinnaCalculationDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, HallArray, eta):

    x = Mu[:]       #all mu on x-axis
    y = HallArray[:]       
    plt.plot(x, y, marker = 'x', label = "$\sigma_{xy}$")     

    plt.ylabel(r'$\sigma_{xy}\left[\frac{e^2}{2\pi\hbar}\right]$', fontsize=25, rotation=90)
    plt.xlabel("$\mu$", fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    #plt.axvline(x = 1.05, label = "half filling", linestyle= "dashed", color="red")
    plt.grid()
    plt.tight_layout()
    plt.legend(fontsize=14, loc='upper right', ncol=1, fancybox=False, shadow=False)
    plt.savefig(NevalinnaCalculationDirectory + '/HallResultEta{}_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}Colloq.pdf'.format(eta, u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.savefig(NevalinnaCalculationDirectory + '/HallResultEta{}_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}Colloq.png'.format(eta, u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.clf()