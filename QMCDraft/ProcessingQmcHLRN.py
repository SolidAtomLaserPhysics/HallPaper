import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess

import createDirectories as CreateDir
import changeToProductionStyle as changeToProd
import extractQMCResults as extract
import PlottingQMC as plot
import createMatrix as createMat





if __name__ == "__main__":
#test of MaxEnt
    U = [1.0]
    Mu = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] 
    Beta = [50.0]
    Q = [3]
    T = [0.25]
    TPrime = [0.0]
    TPrimePrime = [0.0]                                               
    KSteps = [240]

    numberOfCores = 96          #change to 480 for better statistik after precalculation


#True or False dependent on calculations and what I need for that
    createDirectories = False
    preCalculation = True
    renameHDF5 = False
    backUp = False
    copyHDF5 = False
    Version = 0                                         #Version number which I put behind the backup hdf5 (precalc = 0)
    prepareForAndDoNextCalculation = False


#True or False on what we want to Extract/Plot
#also which versions to do
    extractResults = False
    plotResults = False
    VersionToInvestigate = [0]
    #points to plot
    points = 10
    Iterations = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]


#paths for different cases
    QMCDraftDirectory = '/scratch/usr/hhpnhytt/finalQMC/QMCDraft'

    QMCCalculationDirectory = '/scratch/usr/hhpnhytt/finalQMC/QMCMaxentTestHighAccuracy/24-10-23'
    #QMCCalculationDirectory = '/scratch/usr/hhpnhytt/finalQMC/finalQMCResults/QMCSCValues/10-10-23'


    #in future, I will create the directories locally on my computer (with this path)
    #QMCCalculationDirectory = '/afs/physnet.uni-hamburg.de/users/th1_km/nhyttrek/Desktop/MA/MACode/QMC/QMCDraft/24-10-23'
    #QMCDraftDirectory = '/afs/physnet.uni-hamburg.de/users/th1_km/nhyttrek/Desktop/MA/MACode/QMC/QMCDraft'


    '''
    Order of steps:
        create Directories
        make the pre Calculation

        rename the HDF5 to back it up
        make a backup of the precalculation
        (can copy the hdf5 from the iteration before, only copy since need it for next iteration) with Version = 1
        do the next calculation

        backup with Version = 1
        rename hdf5 for the next calculation
        (can copy the hdf5 from the iteration before, only copy since need it for next iteration) with Version = 2
        do the next calculation

        backup with Version = 2
        rename hdf5 for the next calculation
        (can copy the hdf5 from the iteration before, only copy since need it for next iteration) with Version = 3

        NEXT is calculation of version = 4
    '''
                                    


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
                                            if (Version == 0):
                                                os.system('cp {} hdf5BackupFolder/U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT_Version{}.hdf5'.format(hdf5Name, u,mu,beta,q,t,tPri,tPriPri,kSteps, Version))   
                                            else:
                                                #move full .hdf5 since this is the .hdf5 of the old Version and was needed for the current calculation
                                                os.system('mv {} hdf5BackupFolder/U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT_Version{}.hdf5'.format(hdf5Name, u,mu,beta,q,t,tPri,tPriPri,kSteps, Version - 1))      
                                        except:
                                            print("no hdf5 file to backup")
                                        #move the Parameters.in which belongs to that Version into the Backup Folder
                                        try:
                                            os.system('mv Parameters.in hdf5BackupFolder/Parameters_Version{}.in'.format(Version))         
                                        except:
                                            print("no Parameters.in file to backup")
                                        #move the submit_dmft_script which belongs to that Version into the Backup Folder
                                        try:
                                            os.system('mv submit_dmft_script hdf5BackupFolder/submit_dmft_script_Version{}'.format(Version))         
                                        except:
                                            print("no submit file to backup")
                                        #move the slurm which belongs to that Version into the Backup Folder
                                        try:
                                            os.system('mv slurm* hdf5BackupFolder/slurm_Version{}'.format(Version))         
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
                                            os.system('cp {} hdf5BackupFolder/U{}_mu{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_DMFT_Version{}.hdf5'.format(hdf5Name, u,mu,beta,q,t,tPri,tPriPri,kSteps, Version))   
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
                                    







                                    for version in VersionToInvestigate:
                                        '''
                                        will extract results(density and self energy) out of the hdf5 file
                                        '''
                                        if extractResults:
                                            #extract.extractSelfEnergy(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version)
                                            #extract.extractDensity(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version)
                                            for iteration in Iterations:
                                                extract.extractSelfEnergyIterations(QMCCalculationDirectory, QMCDraftDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, version, iteration)





                                if plotResults:
                                    #plot self energy
                                    #TODO: test plot density
                                    for mu in Mu:
                                        #pass
                                        #    plot.plotAllSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps, VersionToInvestigate, points)
                                        #    for version in VersionToInvestigate:
                                        #        plot.plotSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version,30)
                                        plot.plotSigmaAllIterationsQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, VersionToInvestigate, points, Iterations)
                                    
                                    #plot densities
                                    #plot.plotAllDensityQMC(QMCCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps,VersionToInvestigate)
                                    #for version in VersionToInvestigate:
                                    #    plot.plotDensityQMC(QMCCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps, version)
                                    #        plot.plotAntragSigmaQMC(QMCCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri, kSteps, version,20)
                                    







