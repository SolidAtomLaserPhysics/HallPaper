import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess




'''
Use Julians Script to get bath parameters for ED Code
take the hdf5 from QMC
and output a hubb.andpar to the Nevalinna Calculation directory immediately
'''
def getBathParameter(QMCCalculationDirectory,NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    hdf5PathName = '''/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/hdf5BackupFolder/U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT_Version{}.hdf5'''.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, u,mu,beta,q,t,tPri,tPriPri,kSteps,version)
    NevalinnaPathName = '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites)
    process = subprocess.run(['julia /scratch/projects/hhp00048/codes/scripts/LadderDGA_utils/fitAndersonParams_2.jl ' +
                                QMCCalculationDirectory + '{} {} {} {} {} '.format(hdf5PathName,beta,u,mu,bathSites) +
                                NevalinnaCalculationDirectory + '{}/EDCodeCorrect/hubb.andpar_bak'.format(NevalinnaPathName)], capture_output = True, shell = True)
    output = process.stdout.decode("utf-8")



'''
creating directories to prepare for calculations 
'''
def createDirectories(path, Eta):
    if not os.path.exists(path):                                  #make directory if not exists already
        os.mkdir(path)
    if not os.path.exists(path + "/EDCodeCorrect"):
        os.mkdir(path + "/EDCodeCorrect")
    for eta in Eta:
        if not os.path.exists(path + "/NevalinnaEta{}".format(eta)):
            os.mkdir(path + "/NevalinnaEta{}".format(eta))
        if not os.path.exists(path + "/HallEta{}".format(eta)):
            os.mkdir(path + "/HallEta{}".format(eta))







'''
der korrekte ED Code ab Ende November
'''
def EDHubbDatCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, fullyNewED):
    #test iteration number
    if fullyNewED:
        iterations = 1000
    else:
        iterations = 0
    with open(NevalinnaDraftDirectory + '/EDCodeCorrect/hubb.dat', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[1] = " {}d0,  0.d0, 0.0d0 \n".format(u)
        data[3] = " {}, -12.0, 12.0, 0.01 \n".format(beta)
        data[5] = " {}, 0, 0.d0, {},  1.d-9 \n".format(bathSites + 1, iterations)
    #write into that file
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/hubb.dat'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
                file.writelines(data)

def EDHubbAndParCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, fullyNewED):
    if fullyNewED:
        if (bathSites == 6):
            with open(NevalinnaDraftDirectory + '/EDCodeCorrect/hubb.andpar', 'r', encoding='utf-8') as file:
                data = file.readlines()
                data[22] = "  {}          #chemical potential \n".format(mu)
            #write into that file
            with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/hubb.andpar_bak'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
                file.writelines(data)

def EDTPriCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    with open(NevalinnaDraftDirectory + '/EDCodeCorrect/tpri.dat', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[0] = "        t={}d0 \n".format(t)
        data[1] = "        t1={}d0 \n".format(tPri)
        data[2] = "        t2={}d0 \n".format(tPriPri)
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/tpri.dat'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)

#just copy the init.h
def InitHCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    with open(NevalinnaDraftDirectory + '/EDCodeCorrect/init.h', 'r', encoding='utf-8') as file:
        data = file.readlines()
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/init.h'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)

def EDSubmitCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    with open(NevalinnaDraftDirectory + '/EDCodeCorrect/submit_EDdmft_script', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[1] = "#SBATCH -t 10:00:00 \n"
        data[2] = "#SBATCH -n 960 \n"               #10 nodes if produktion, 1 for testing
        data[5] = "#SBATCH -A hhpnhytt  \n"
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/submit_EDdmft_script'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)

#run is not used at the moment
def EDRunCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    with open(NevalinnaDraftDirectory + '/EDCodeCorrect/run.sh', 'r', encoding='utf-8') as file:
        data = file.readlines()
        data[6] = "mpirun -np {} ./run.x > run.out \n".format((bathSites+2)*(bathSites+2))
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/run.sh'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)

def EDCodeCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, fullyNewED):
    with open(NevalinnaDraftDirectory + '/EDCodeCorrect/ed_dmft_parallel_frequencies_corrected_and_commented_betterSampling_november.f', 'r', encoding='utf-8') as file:
        data = file.readlines()
        if fullyNewED:
            data[165] = "      symm=.false.  \n"                #use it for ED Calculations, when I am in Georgs system and this line si allowed
        data[17] = "      parameter (nss={})  \n".format(bathSites+1)
        #data[18] = "      parameter (prozessoren={})       \n".format((bathSites+2) * (bathSites+2))
        data[18] = "      parameter (prozessoren={})  \n".format(960)  #only testing 64, mostly 960 (10 nodes) when annealing
        data[23] = "      parameter (ksteps={}) \n".format(kSteps)
        data[26] = "      parameter (L={})   \n".format(q)

    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/ed_dmft_parallel_frequencies_corrected_and_commented_betterSampling_november.f'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)


def EDAnnealingScript(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites):
    with open(NevalinnaDraftDirectory + '/submit_dmft_scan_script', 'r', encoding='utf-8') as file:
        data = file.readlines()

        data[18] = "EDdir=\"{}\" \n".format(NevalinnaCalculationDirectory)
        data[19] = "muArray=("
        for imu in range(len(Mu) - 1):
            data[19] = data[19] + "{} ".format(Mu[imu])
        data[19] = data[19] + "{})\n".format(Mu[len(Mu)-1])
        data[21] = "U={} \n".format(u)
        data[22] = "B={} \n".format(beta)
        data[23] = "q={} \n".format(q)
        data[24] = "t={} \n".format(t)
        data[25] = "tPri={} \n".format(tPri)
        data[26] = "tPriPri={} \n".format(tPriPri)
        data[27] = "kSteps={} \n".format(kSteps)
        data[28] = "version={} \n".format(version)
        data[29] = "bathSites={} \n".format(bathSites)


    with open(NevalinnaCalculationDirectory + '/submit_dmft_scan_script_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}'.format(u,beta,q,t,tPri,tPriPri,kSteps,version,bathSites), 'w', encoding='utf-8') as file:
        file.writelines(data)

