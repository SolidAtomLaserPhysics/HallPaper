import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import subprocess


'''
Hall by Nevanlinna with different eta
'''

def copyHallCodeEta(HallDraftDirectory,NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta):
    with open(HallDraftDirectory + '/hall_conductivity.cpp', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[21] = "const double t = {}; \n".format(t)
        data[22] = "const double tPrime = {}; \n".format(tPri)
        data[23] = "const double tPrimePrime = {}; \n".format(tPriPri)
        path = NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta)
        data[451] = "    fstream Sigma_file(\"{}/selfEnergyHall.dat\",fstream::in); \n".format(path)
        data[579] = "  param.open(\"{}/parameters.dat\", fstream::in); \n".format(path)

    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}/hall_conductivity.cpp'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), 'w', encoding='utf-8') as file:
        file.writelines(data)

def copyHallHeaderCodeEta(HallDraftDirectory,NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta):
    with open(HallDraftDirectory + '/hall_conductivity.h', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[14] = "const int L ={}; \n".format(q)

    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}/hall_conductivity.h'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), 'w', encoding='utf-8') as file:
        file.writelines(data)


def setParametersHallEta(NevalinnaCalculationDirectory,HallDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, filling, wPoints,wMin,wMax, eta):
    p = 1           #set p = 1 fixed since will not do anything else
    with open(HallDraftDirectory + '/parameters.dat', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[1] = "> {}    {}     {}     {} \n".format(u, mu, beta, filling)
        data[3] = "> {}    {}  \n".format(p, q)
        data[5] = "> {}     {}     {} \n".format(wPoints, wMin, wMax)
        kPoints = kSteps
        data[7] = "> {}       {} \n".format(kPoints, kPoints)
        #if U=0 then eta is not the Nevanlinna eta, but Eta in the Hall Code
        if (u==0):
            data[10] = "> 0.001 {} 5 \n".format(eta)
        else:
            data[10] = "> 0.001 0.01 5 \n"

    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}/parameters.dat'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), 'w', encoding='utf-8') as file:
        file.writelines(data)


def hallSubmitEta(NevalinnaCalculationDirectory,HallDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, numberOfCores, eta):
    with open(HallDraftDirectory + '/submit_hall', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[1] = "#SBATCH -t 6:00:00 \n"
        data[2] = "#SBATCH -n {} \n".format(numberOfCores)
        data[3] = "#SBATCH --ntasks-per-node 96 \n"
        data[4] = "#SBATCH -p standard96   \n"         #need large memory for calculation
        data[5] = "#SBATCH -A hhpnhytt  \n"
        data[7] = "export OMP_NUM_THREADS={} \n".format(numberOfCores)

        #one > means overwriting, >> means anhaengen
        data[10] = "./hall_conductivity > runHall.out"  # 2> runHall.err"

    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}/submit_hall'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), 'w', encoding='utf-8') as file:
        file.writelines(data)



'''
since only self energy looking like Greens function continued
recalculate it to normal looking self energy (with hartree term when HartreeTerm = True)
HartreeTerm = False only useful when added later (e.g. in Hall Code)
'''
def normalizeNevalinnaOutputSelfEnergy(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, n, HartreeTerm, eta):
    sigmaValuesOutput = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/output.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
    for i in range(len(sigmaValuesOutput)):
        sigmaValuesOutput[i,2] = sigmaValuesOutput[i,2] * (u*u *n/2 * (1-n/2))
        if HartreeTerm:
            pass
        else:
            sigmaValuesOutput[i,1] = sigmaValuesOutput[i,1] * (u*u *n/2 * (1-n/2))
            np.savetxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}/selfEnergyHall.dat'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), sigmaValuesOutput)

