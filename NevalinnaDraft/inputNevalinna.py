import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import subprocess





'''
get filling of ED, search in run.out for densimp
'''
def getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/run.out'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'r', encoding='utf-8') as file:
        data = file.readlines()
        for line in range(len(data)):
            #go through each line from bottom to top
            content = data[len(data) - 1 - line]
            #if densimp in that line
            if (content.find("densimp") != -1):
                #find all numbers in that line with densimp
                match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
                final_list = [float(x) for x in re.findall(match_number, content)]

                #if found densimp the first time, stop for
                return float(final_list[0])
                continue


'''
get self energy function self-en_wim from ED and make it look like Greens function
'''
def NevalinnaSelfEnergyData(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, n, eta):
    #load self energy
    sigmaValues = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/self-en_wim'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites))
    greensFunctionValues = np.zeros((len(sigmaValues), 3))
    #frequency values stay the same
    for i in range(len(sigmaValues)):
        greensFunctionValues[i,0] = sigmaValues[i,0]
        #real part
        greensFunctionValues[i,1] = (sigmaValues[i,1] - u*n/2)/(u*u *n/2 * (1-n/2))
        #imag part
        greensFunctionValues[i,2] = sigmaValues[i,2]/(u*u *n/2 * (1-n/2))
    np.savetxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/data.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), greensFunctionValues)




'''
build input.txt of C Nevalinna Code
ifile imag_num ofile
TODO: find a good imag_num number
'''
def NevalinnaInputTxt(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, eta, imag_num):
    with open(NevalinnaDraftDirectory + '/input.txt', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[0] = NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/data.txt '.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta) + str(imag_num) + ' ' + NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/output.txt '.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta)
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/input.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), 'w', encoding='utf-8') as file:
        file.writelines(data)



def NevalinnaSubmitAll(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, Eta):
    with open(NevalinnaDraftDirectory + '/submit_nevalinna', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[1] = "#SBATCH -t 1:00:00 \n"
        data[2] = "#SBATCH -N 4 \n"
        data[5] = "#SBATCH -A hhpnhytt  \n"

        for m, mu in enumerate(Mu):
            for e, eta in enumerate(Eta):
                path = NevalinnaCalculationDirectory + "/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}".format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites,eta)
                #data[9 + 2*e*len(Mu) + 2*m] = "g++ -o {}/nevanlinna {}/nevanlinna.cpp -I /home/hhpnhytt/eigen-3.4.0 -lgmp -lgmpxx \n".format(path, path)
                data[10 + 2*e*len(Mu) + 2*m] = "{}/nevanlinna <  ".format(path) + '{}/input.txt    &   \n'.format(path)
        data[16 + 4*len(Mu)*len(Eta)] = "wait \n"

    with open(NevalinnaCalculationDirectory + '/submit_nevalinna_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}'.format(u,beta,q,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)





'''
copy the Code into ED directory
Also changes values in nevanlinna.h file
'''
def NevalinnaCode(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, wPoints,wMin,wMax, eta):
    os.system('cp ' + NevalinnaDraftDirectory + '/nevanlinna.cpp ' + NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/nevanlinna.cpp'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
    os.system('cp ' + NevalinnaDraftDirectory + '/schur.h ' + NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/schur.h'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
    with open(NevalinnaDraftDirectory + '/nevanlinna.h', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[73] = "    real_domain_data (std::string ofile) : ofs(ofile), N_real_({}), omega_min({}), omega_max({}), eta({}) {} \n".format(wPoints, wMin, wMax, eta,"{")     #python cant have { in string since thinks } forgotten
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/nevanlinna.h'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), 'w', encoding='utf-8') as file:
        file.writelines(data)




def DOSScript(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, Eta):
    
    

    for mu in Mu:
        for eta in Eta:
            with open(NevalinnaDraftDirectory + '/calcDOS.py', 'r', encoding='utf-8') as file:
                data = file.readlines()
                data[6] = "kStepsFile = {} \n".format(kSteps)
                data[9] = "p = 1 \n"
                data[10] = "q = {} \n".format(q)
                data[11] = "kSteps = {} \n".format(120)                 #kpoints for calculation of DOS
                data[12] = "u = {} \n".format(u)
                data[13] = "beta = {} \n".format(beta)
                data[14] = "version = {} \n".format(version)
                data[15] = "bathSites = {} \n".format(bathSites)
                data[16] = "t = {} \n".format(t)
                data[17] = "tPri = {} \n".format(tPri)
                data[18] = "tPriPri = {} \n".format(tPriPri)
                data[19] = "mu = {} \n".format(mu)
                data[20] = "eta = {} \n".format(eta)
                data[21] = "NevalinnaCalculationDirectory =  \"{}\" \n".format(NevalinnaCalculationDirectory)


            with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/calcDOS.py'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites,eta), 'w', encoding='utf-8') as file:
                file.writelines(data)



    with open(NevalinnaDraftDirectory + '/submit_DOS', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[1] = "#SBATCH -t 12:00:00 \n"
        data[2] = "#SBATCH -N 1 \n"


        for m, mu in enumerate(Mu):
            for e, eta in enumerate(Eta):
                path = NevalinnaCalculationDirectory + "/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}".format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites,eta)
                data[10 + e*len(Mu) + m] = ' python3 {}/calcDOS.py  &   \n'.format(path)
        data[16 + len(Mu)*len(Eta)] = "wait \n"

    with open(NevalinnaCalculationDirectory + '/submit_DOS_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}'.format(u,beta,q,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)


    