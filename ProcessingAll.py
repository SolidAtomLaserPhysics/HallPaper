import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess



import QMCDraft.createDirectories as CreateDir
import QMCDraft.changeToProductionStyle as changeToProd
import QMCDraft.extractQMCResults as extract
import QMCDraft.PlottingQMC as plot
import QMCDraft.createMatrix as createMat

import NevalinnaDraft.inputED as inputED
import NevalinnaDraft.inputNevalinna as inputNevalinna
import NevalinnaDraft.plottingNevalinna as plotNevalinna

import HallDraft.inputHall as inputHall
import HallDraft.plottingHall as plotHall





if __name__ == "__main__":
#test of MaxEnt
    U = [2.0]
    #full for U0 but less
    #Mu = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
    #      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    #full U0
    #Mu = [-1.2,-1.15,-1.1,-1.05,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,
    #      0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
    #full U1
    #Mu = [-1.2,-1.15,-1.1,-1.05,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,
    #      0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,
    #      1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0] 
    #full U2
    Mu = [-1.2,-1.15,-1.1,-1.05,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,
          0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
          #1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0,
          #2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,3.0]
    #full U3
    #Mu = [-1.2,-1.15,-1.1,-1.05,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,
    #      0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,
    #      1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0,
    #      2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,3.0,
    #      3.05,3.1,3.15,3.2,3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9]    #3.95,4
    #full for U1 but less
    #Mu = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
    #      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
    #      1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
    #full for U2 but less
    #Mu = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
    #      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
    #      1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
    #      2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
    #full for U3 but less
    #Mu = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
    #      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
    #      1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
    #      2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
    #      3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]
    #single mu
    #Mu = [-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]          #2.2, 1.0, 1.8
    #Mu = [0.5, 1.0]

    Beta = [50.0]
    Q = [4]
    T = [0.25]
    TPrime = [-0.125]
    TPrimePrime = [0.0]
    KSteps = [240]
    wPoints = 2000
    wMin = -5.0             #usually -5 to 5, for U =0 also -1,1, for U1 and U2 at least -4,4, U3 -3,3          #Achtung U=3 gerade -6,6
    wMax = 5.0
    #imag_num is the number of Matsubara frequencies I use in Nevanlinna
    imag_num = 40
    #Eta U0
    #Eta = [0.01, 0.05, 0.1]      #no 0.001
    #Eta U1
    #Eta = [0.05, 0.1, 0.15, 0.2]                   #0.005 und 0.01 eh zu schlecht und SC hat die nicht
    #Eta U2
    Eta = [0.005, 0.01, 0.05, 0.1, 0.15, 0.2]      #
    #Eta U3
    #Eta = [0.005, 0.01, 0.05, 0.1, 0.15]           #U3 SC hat 0.2 nicht aber correlation verschmiert eh zu sehr??
    #Eta = [0.05]

    numberOfCores = 96       #for Hall code set this to 1 node = 96   #change to 480 for better statistik after precalculation


#QMC
    createDirectories = False
    preCalculation = False
    renameHDF5 = False
    backUp = False
    copyHDF5 = False
    version = 1                                         #Version number which I put behind the backup hdf5 (precalc = 0)
    prepareForAndDoNextCalculation = False
#Plotting/Extracting
    extractResults = False
    plotResults = False
    VersionToInvestigate = [1]
    points = 60
    #Iterations = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35, 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]
    Iterations = [1,2,3,4,5,6,7,8,9,10]
    plotQMCEDComparison = False


#Nevalina
    prepareDirectories = False
    prepareEdInput = False
    calculateEdData = False
    EDAnnealing = False
    prepareNevalinnaInput = False
    submitNevalinna = False
    submitAllNevalinna = False
    bathSites = 6


#Plotting und so
    plotDensityED = False
    calculateDOS = False           #Only Maxent
    plotDOS = False
    prepareDOSNevalinna = False
    calculateDOSScript = False
    plotDOSNevalinna = False
    plotGreensFunction = False

#Hall
    prepareHall = False
    submitHall = True
    showHall = False


#paths for different cases
    QMCDraftDirectory = '/scratch/usr/hhpnhytt/QMCDraft'
    MaxEntDraftDirectory = '/scratch/usr/hhpnhytt/MaxEntDraft'
    NevalinnaDraftDirectory = '/scratch/usr/hhpnhytt/NevalinnaDraft'
    HallDraftDirectory = '/scratch/usr/hhpnhytt/HallDraft'

    #QMCDraftDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/QMCDraft'
    #MaxEntDraftDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/MaxEntDraft'
    #NevalinnaDraftDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/NevalinnaDraft'


    #QMCCalculationDirectory = '/scratch/usr/hhpnhytt/finalQMC/finalQMCResults/QMCSCValues/lastCodes'     #originalDMFTpy
    QMCCalculationDirectory = '/scratch/projects/hhp00048/NicoQMCData/finalQMC/finalQMCResults/QMCSCValues/lastCodes'

    MaxEntCalculationDirectory = '/scratch/usr/hhpnhytt/MaxEntTests/QMCMaxEntHighAccuracy/8-11-23'
    #MaxEntCalculationDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/MaxEntValues/20-10-23'


    NevalinnaCalculationDirectory = '/scratch/usr/hhpnhytt/Results/q4/allFlat/Calculation/U2'
    #NevalinnaCalculationDirectory = '/scratch/usr/hhpnhytt/Results/newMySystem/U3Annealing'



    #On Physnet
    #QMCCalculationDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/23-11-23'
    #QMCDraftDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/QMCDraft'
    #MaxEntDraftDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/MaxEntDraft'
    #MaxEntCalculationDirectory = '/afs/physnet.uni-hamburg.de/users/th1_ro/nhyttrek/Desktop/MA/MACode/QMC/MaxEntValues/7-11-23'





    '''
    Order of steps:
        create Directories
        make the pre Calculation

        rename the HDF5 to back it up
        make a backup of the precalculation
        (can copy the hdf5 from the iteration before, only copy since need it for next iteration) with version = 1
        do the next calculation

        rename hdf5 for the next calculation (and overwritte of version before in main folder)
        backup with version = 1
        (can copy the hdf5 from the iteration before, only copy since need it for next iteration) with version = 2
        do the next calculation

        rename hdf5 for the next calculation (and overwritte of version before in main folder)
        backup with version = 2
        (can copy the hdf5 from the iteration before, only copy since need it for next iteration) with version = 3

        NEXT is calculation of version = 4
    '''
    #QMC
    for u in U:
        for beta in Beta:
            for q in Q:
                for t in T:
                    for tPri in TPrime:
                        for tPriPri in TPrimePrime:
                            for kSteps in KSteps:
                                for mu in Mu:
                                    '''
                                    all directories and files will be generated and put into its folder
                                    after this, everything is ready to start the Calculation
                                    '''
                                    if createDirectories:
                                        CreateDir.createDirectories(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps))
                                        print("dir created")
                                        CreateDir.createParameters(u, mu, beta, q, t, tPri, tPriPri, kSteps, QMCDraftDirectory, QMCCalculationDirectory)
                                        print("parameters.in created")
                                        CreateDir.createSubmit(QMCDraftDirectory, QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, numberOfCores)
                                        print("submit created")
                                        createMat.writeHKFile(1, q, t, tPri, tPriPri, kSteps, QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Hk_Hofstadter.dat'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps))
                                        print("epsilon Matrix calculated")

                                    '''
                                    start the first Calculation to get starting point of nCorr
                                    '''
                                    if preCalculation:
                                        os.chdir(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps))
                                        os.system('sbatch submit_dmft_script' )

                                    '''
                                    for that you rename the latest hdf5 file to an easy name
                                    '''
                                    if renameHDF5:
                                        os.chdir(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps))
                                        #rename hdf5 from last calc to an easier name
                                        hdf5Name = 'U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT.hdf5'.format(u,mu,beta,q,t,tPri,tPriPri,kSteps)
                                        #renaming itself
                                        try:
                                            os.system('mv U{}_mu{}_B{}_L* {}'.format(u,mu,beta, hdf5Name))
                                        except:
                                            print("no new hdf5 file")

                                    '''
                                    Here I manually backup
                                    I backup after the calculation of Version n and after renaming the hdf5
                                    I save hdf5, parameters.in and submit with its Version into the backup Folder.
                                    With this I can resolve all calculations
                                    '''
                                    if backUp:
                                        os.chdir(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps))
                                        hdf5Name = 'U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT.hdf5'.format(u,mu,beta,q,t,tPri,tPriPri,kSteps)
                                        #copy old/renamed hdf5 to a safe place and remove the old version then
                                        try:
                                            if (version == 0):
                                                os.system('cp {} hdf5BackupFolder/U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT_Version{}.hdf5'.format(hdf5Name, u,mu,beta,q,t,tPri,tPriPri,kSteps, version))
                                            else:
                                                #move full .hdf5 since this is the .hdf5 of the old Version and was needed for the current calculation
                                                os.system('cp {} hdf5BackupFolder/U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT_Version{}.hdf5'.format(hdf5Name, u,mu,beta,q,t,tPri,tPriPri,kSteps, version))
                                        except:
                                            print("no hdf5 file to backup")
                                        #move the Parameters.in which belongs to that Version into the Backup Folder
                                        try:
                                            os.system('mv Parameters.in hdf5BackupFolder/Parameters_Version{}.in'.format(version))
                                        except:
                                            print("no Parameters.in file to backup")
                                        #move the submit_dmft_script which belongs to that Version into the Backup Folder
                                        try:
                                            os.system('mv submit_dmft_script hdf5BackupFolder/submit_dmft_script_Version{}'.format(version))
                                        except:
                                            print("no submit file to backup")
                                        #move the slurm which belongs to that Version into the Backup Folder
                                        try:
                                            os.system('mv slurm* hdf5BackupFolder/slurm_Version{}'.format(version))
                                        except:
                                            print("no slurm file to backup")

                                    '''
                                    can manually copy hdf5 ready for next calculation into the hdf5BackupFolder
                                    With that I can investigate it allready, although it was not backuped
                                    This happens AFTER renaming the hdf5 file AFTER backup
                                    '''
                                    if copyHDF5:
                                        os.chdir(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps))
                                        hdf5Name = 'U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT.hdf5'.format(u,mu,beta,q,t,tPri,tPriPri,kSteps)
                                        try:
                                            os.system('cp {} hdf5BackupFolder/U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT_Version{}.hdf5'.format(hdf5Name, u,mu,beta,q,t,tPri,tPriPri,kSteps, version))
                                        except:
                                            print("no hdf5 file to copy")

                                    '''
                                    prepare and do the next calculation
                                    '''
                                    if prepareForAndDoNextCalculation:
                                        hdf5Name = 'U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT.hdf5'.format(u,mu,beta,q,t,tPri,tPriPri,kSteps)
                                        #write nCorr into Parameters.in so that it is ready for real calculation
                                        process = subprocess.run(['julia  /scratch/projects/hhp00048/codes/scripts/LadderDGA_utils/ncorr.jl ' + QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/{}'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, hdf5Name)]
                                                               ,capture_output = True, shell = True)
                                        print(process.stdout.decode("utf-8"))
                                        nCorr = int(process.stdout.decode("utf-8"))
                                        print(nCorr)
                                        changeToProd.createSubmit(QMCDraftDirectory, QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,numberOfCores)
                                        changeToProd.createParameters(u, mu, beta, q, t, tPri, tPriPri, kSteps, QMCDraftDirectory, QMCCalculationDirectory, nCorr)
                                        os.chdir(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps))
                                        os.system('sbatch submit_dmft_script' )


                                    if extractResults:
                                        '''
                                        will extract results(density and self energy) out of the hdf5 file
                                        '''
                                        for version in VersionToInvestigate:
                                            #extract.extractSelfEnergy(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version)
                                            #extract.extractDensity(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version)
                                            try:
                                                for iteration in Iterations:
                                                    extract.extractSelfEnergyIterations(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, version, iteration)
                                            except:
                                                pass





    '''
    Next there are all steps of Nevalina
    '''
    for u in U:
        for beta in Beta:
            for q in Q:
                for t in T:
                    for tPri in TPrime:
                        for tPriPri in TPrimePrime:
                            for kSteps in KSteps:         
                                '''
                                just preparing directories 
                                '''
                                for mu in Mu:
                                    if prepareDirectories:
                                        inputED.createDirectories(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites), Eta)

                                '''
                                ED stuff 
                                '''
                                for mu in Mu: 
                                    if prepareEdInput:
                                        inputED.EDHubbDatCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, True)
                                        inputED.EDHubbAndParCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, True)
                                        inputED.EDTPriCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                        inputED.InitHCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                        inputED.EDSubmitCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                        inputED.EDRunCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                        inputED.EDCodeCorrect(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, True)
                                        #inputED.getBathParameter(QMCCalculationDirectory,NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)

                                if EDAnnealing:
                                    inputED.EDAnnealingScript(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                    os.chdir(NevalinnaCalculationDirectory)
                                    os.system('sbatch submit_dmft_scan_script_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}'.format(u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
                                #if some mu but annealing all mu is better
                                for mu in Mu:
                                    if calculateEdData:
                                        os.chdir(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect'.format(u,
                                                                                    beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites))
                                        os.system('sbatch submit_EDdmft_script' )



                                '''
                                Navanlinna stuff
                                '''
                                for mu in Mu:
                                    if prepareNevalinnaInput:
                                        filling = inputNevalinna.getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                        print("filling {}".format(filling))
                                        for eta in Eta:
                                            inputNevalinna.NevalinnaSelfEnergyData(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, filling, eta)
                                            inputNevalinna.NevalinnaInputTxt(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, eta,imag_num)
                                            inputNevalinna.NevalinnaCode(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, wPoints,wMin,wMax, eta)
                                            #compile code I justed copied to that directory and build submit
                                            os.chdir(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}'.format(u,
                                                                                        beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
                                            #has to load gcc/13.2.0 to have g++ work        module load gcc/13.2.0
                                            os.system('g++ -o nevanlinna nevanlinna.cpp -I /home/hhpnhytt/eigen-3.4.0 -lgmp -lgmpxx' )
                                if prepareNevalinnaInput:
                                    #build Nevalinna submit for many calcs
                                    inputNevalinna.NevalinnaSubmitAll(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites,Eta)
                                #submit nevanlinna
                                if submitAllNevalinna:
                                    os.chdir(NevalinnaCalculationDirectory)
                                    os.system('sbatch submit_nevalinna_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}'.format(u,beta,q,t,tPri,tPriPri,kSteps,version,bathSites))


                                '''
                                Hall stuff 
                                '''
                                for mu in Mu:
                                    if prepareHall:
                                        #get filling from ED results
                                        if (u != 0):
                                            filling = inputNevalinna.getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                        else:
                                            filling = 0
                                        #all steps for eta values
                                        for eta in Eta:
                                            #if U=0 have self energy 0
                                            if (u != 0):
                                                inputHall.normalizeNevalinnaOutputSelfEnergy(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, filling, False, eta)
                                            #copy codes
                                            os.system("cp {}/Makefile {}/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}/Makefile".format(HallDraftDirectory,
                                                                                    NevalinnaCalculationDirectory,u, beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
                                            inputHall.copyHallCodeEta(HallDraftDirectory,NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta)
                                            inputHall.copyHallHeaderCodeEta(HallDraftDirectory,NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta)
                                            inputHall.setParametersHallEta(NevalinnaCalculationDirectory,HallDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, filling, wPoints,wMin,wMax, eta)
                                            os.chdir(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}'.format(u,
                                                                                    beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
                                            os.system("make")
                                            inputHall.hallSubmitEta(NevalinnaCalculationDirectory,HallDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, numberOfCores, eta)

                                    if submitHall:
                                        for eta in Eta:
                                            os.chdir(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/HallEta{}'.format(u,
                                                                                    beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
                                            print(eta)
                                            os.system('sbatch submit_hall')





                                '''
                                Now some Plotting
                                '''
                                if plotResults:
                                    #plot self energy
                                    #for mu in Mu:
                                    #    plot.plotSigmaQMCED(QMCCalculationDirectory, NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, version,points, bathSites)
                                    #TODO: test plot density
                                    for mu in Mu:
                                        #    plot.plotAllSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, VersionToInvestigate, points)
                                        #    for version in VersionToInvestigate:
                                        #        plot.plotSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version,30)
                                        plot.plotSigmaAllIterationsQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, VersionToInvestigate, points, Iterations)

                                    #plot densities
                                    #plot.plotAllDensityQMC(QMCCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps,VersionToInvestigate)
                                    #for version in VersionToInvestigate:
                                    #    plot.plotDensityQMC(QMCCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps, version)
                                    #        plot.plotAntragSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, version,20)



                                if plotQMCEDComparison:
                                    EDFilling = np.zeros(len(Mu))
                                    QMCFilling = np.zeros(len(Mu))
                                    versionED = 1
                                    versionQMC = 1
                                    #for m, mu in enumerate(Mu):
                                        #there is only 1 ED version, but you have to see which
                                        #EDFilling[m] = inputNevalinna.getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, versionED, bathSites)
                                        #version of QMC you want to compare with
                                        #QMCFilling[m] = inputMaxEnt.getFillingQMC(QMCCalculationDirectory,QMCDraftDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, versionQMC, bathSites)
                                    #plot.plotDensityQMCED(QMCCalculationDirectory, NevalinnaCalculationDirectory, EDFilling,QMCFilling, u,beta,q,Mu,t,tPri,tPriPri, kSteps, versionED,points, bathSites)
                                    for mu in Mu:
                                        dataSelfEnergyED = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/self-en_wim'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,versionED,bathSites))
                                        extract.extractSelfEnergy(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,versionQMC)
                                        dataSelfEnergyQMC = np.loadtxt(QMCCalculationDirectory + '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/self-en_wim_QMC_{}.csv'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps, versionQMC))
                                        plot.plotSigmaQMCED(QMCCalculationDirectory, NevalinnaCalculationDirectory, dataSelfEnergyED,dataSelfEnergyQMC, u,beta,q,mu,t,tPri,tPriPri, kSteps, versionED,points, bathSites)



                                if plotDensityED:
                                    fillingArray = np.zeros(len(Mu))
                                    for m, mu in enumerate(Mu):
                                        fillingArray[m] = inputNevalinna.getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                    plotNevalinna.plotDensityED(NevalinnaCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps,version, bathSites, fillingArray)



                                if calculateDOSScript:
                                    inputNevalinna.DOSScript(NevalinnaCalculationDirectory,NevalinnaDraftDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, Eta)
                                    os.chdir(NevalinnaCalculationDirectory)
                                    #os.system('sbatch submit_DOS_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}'.format(u,
                                    #                                                beta,q,t,tPri,tPriPri,kSteps,version,bathSites))

                                #if prepareDOSNevalinna:
                                    #Matrix = plotNevalinna.createMatrix(p, q, t, tPri, tPriPri, kSteps)
                                for mu in Mu:
                                    '''
                                    calculates DOS like in Georgs paper using the selfEnergy and energy matrix
                                    '''
                                    if prepareDOSNevalinna:
                                        print(mu)
                                        filling = inputNevalinna.getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites)
                                        for eta in Eta:
                                            plotNevalinna.normalizeNevalinnaOutputSelfEnergy(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, filling, True, eta)
                                            #SelfEnergy = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/selfEnergyHartreeForDOS.txt'.format(u,
                                            #                                        beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
                                            #NNevalinna = len(SelfEnergy)
                                            #A = plotNevalinna.calculateDOS(q,kSteps,mu, Matrix,SelfEnergy,NNevalinna)
                                            #np.savetxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/DOS.dat'.format(u,
                                            #                                        beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), A)


                                    if plotDOSNevalinna:
                                        noninteracting = np.loadtxt(NevalinnaCalculationDirectory + "/DOSNonInteracting/DensityOfStates.dat")
                                        plotNevalinna.plotDOSNoninteracting(NevalinnaCalculationDirectory,noninteracting)
                                        #plots Sigma(nu) directly after Nevanlinna
                                        #plotNevalinna.plotSigmaRenormalizedAll(NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites, Eta)
                                        #plotNevalinna.plotDOSAll(NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites, Eta)

                                        #for eta in Eta:
                                        #    A = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/DOS.dat'.format(u,
                                        #                                                beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
                                        #    SelfEnergyForFrequencies = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/output.txt'.format(u,
                                        #                                                beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
                                        #    plotNevalinna.plotDOS(NevalinnaCalculationDirectory, A,SelfEnergyForFrequencies, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites, eta)



                                if showHall:
                                    #plot of Hall with All Eta as comparison
                                    plotHall.plotHallAllEta(NevalinnaCalculationDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, HallArrayDirect, Eta)

                                    #plot Hall of each eta individually
                                    HallArrayEta = np.zeros(len(Mu))
                                    for eta in Eta:
                                        for m, mu in enumerate(Mu):
                                            print(mu)
                                            HallArrayEta[m] = plotHall.getHallValueEta(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, eta)
                                        plotHall.plotHallEtaColloq(NevalinnaCalculationDirectory, u,beta,Mu,q,t,tPri,tPriPri,kSteps, version, bathSites, HallArrayEta, eta)
