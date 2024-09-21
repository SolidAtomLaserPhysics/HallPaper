import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess
import math
import cmath

from scipy import linalg







'''
OUT: epsilepsMatrix Array[k,k, q,q]
'''
def createMatrix(p, q, t, tPri, tPriPri, kSteps):


    '''
    set variables
    '''
    j = complex(0,1)                                                #gives the imaginary number i
    B = p/q
    epsMatrix = np.zeros((kSteps, kSteps, q, q), dtype=complex)
    TMatrix = np.zeros((kSteps, kSteps, q, q), dtype=complex)
    TPriMatrix = np.zeros((kSteps, kSteps, q, q), dtype=complex)
    TPriPriMatrix = np.zeros((kSteps, kSteps, q, q), dtype=complex)

    '''
    build kPoints
    '''
    dkx = (2*math.pi)/kSteps
    dky = (2*math.pi)/(kSteps*q)

    '''
    build values 
    '''
    for ikx in range(kSteps):
        for iky in range(kSteps):
            kx = -math.pi + ikx * dkx
            ky = -math.pi/q + iky * dky

            if (q == 3):
                '''
                T Matrix
                '''
                TMatrix[ikx, iky, 0,0] = 2*t*math.cos(kx)
                TMatrix[ikx, iky, 0,1] = t
                TMatrix[ikx, iky, 0,2] = t*cmath.exp(j*ky*q)

                TMatrix[ikx, iky, 1,0] = t
                TMatrix[ikx, iky, 1,1] = 2*t*math.cos(kx + (2*math.pi*B))
                TMatrix[ikx, iky, 1,2] = t

                TMatrix[ikx, iky, 2,0] = t*cmath.exp(-j*ky*q)
                TMatrix[ikx, iky, 2,1] = t
                TMatrix[ikx, iky, 2,2] = 2*t*math.cos(kx + (4*math.pi*B))

                '''
                TPrime Matrix
                '''
                TPriMatrix[ikx, iky, 0,0] = 0
                TPriMatrix[ikx, iky, 0,1] = 2*tPri * math.cos(kx + (math.pi*B))
                TPriMatrix[ikx, iky, 0,2] = 2*tPri * math.cos(kx - (math.pi*B)) * cmath.exp(j*ky*q)

                TPriMatrix[ikx, iky, 1,0] = 2*tPri * math.cos(kx + (math.pi*B))
                TPriMatrix[ikx, iky, 1,1] = 0
                TPriMatrix[ikx, iky, 1,2] = 2*tPri * math.cos(kx + ((3)*math.pi*B))

                TPriMatrix[ikx, iky, 2,0] = 2*tPri * math.cos(kx - (math.pi*B)) * cmath.exp(-j*ky*q)
                TPriMatrix[ikx, iky, 2,1] = 2*tPri * math.cos(kx + ((3)*math.pi*B))
                TPriMatrix[ikx, iky, 2,2] = 0


                '''
                TPrimePrime Matrix
                '''
                TPriPriMatrix[ikx, iky, 0,0] = 2*tPriPri * math.cos(2*kx)
                TPriPriMatrix[ikx, iky, 0,1] = tPriPri * cmath.exp(j*ky*q)
                TPriPriMatrix[ikx, iky, 0,2] = tPriPri

                TPriPriMatrix[ikx, iky, 1,0] = tPriPri * cmath.exp(-j*ky*q)
                TPriPriMatrix[ikx, iky, 1,1] = 2*tPriPri * math.cos(2*kx + (4*math.pi*B))
                TPriPriMatrix[ikx, iky, 1,2] = tPriPri * cmath.exp(j*ky*q)

                TPriPriMatrix[ikx, iky, 2,0] = tPriPri
                TPriPriMatrix[ikx, iky, 2,1] = tPriPri * cmath.exp(-j*ky*q)
                TPriPriMatrix[ikx, iky, 2,2] = 2*tPriPri * math.cos(2*kx + (8*math.pi*B))

                epsMatrix = TMatrix + TPriMatrix + TPriPriMatrix

    return epsMatrix




'''
calculate and plot DOS which is A(w) = 
IN:
    Matrix Array[k,k, q,q]
    SelfEnergy Array[NNevalinna, 5]
    NNevalinna number of real frequencies omega
'''
def calculateDOS(q,kSteps,mu, Matrix,SelfEnergy,NNevalinna):
    '''
    initialize
    '''
    denominator = np.zeros((q,q), dtype=complex)
    A = np.zeros(NNevalinna, dtype=float)
    imag = complex(0,1)


    for w in range(NNevalinna):  
        if (SelfEnergy[w, 0] > -5) and (SelfEnergy[w, 0] < 5):  
            for ikx in range(kSteps):
                for iky in range(kSteps):  
                    for i in range(q):
                        for j in range(q):
                            #TODO: sit das wirklich richtig so?
                            denominator[i,j] = (SelfEnergy[w, 0] + mu - (SelfEnergy[w, 1] + imag*SelfEnergy[w,2]))*np.identity(q)[i,j] - Matrix[ikx,iky, i,j]          
                    insideMatrix = np.linalg.inv(denominator) * 1/(kSteps*kSteps * 2*math.pi)                 #Matrix inversion and prefactors to 

                    A[w] = A[w] + np.imag(np.trace(insideMatrix))
            #print("one w done and w= " + w)
            #print(A[w])

    return (A/q)

                         

def plotDOSNoninteracting(NevalinnaCalculationDirectory,noninteracting):
    x = noninteracting[:,0]        #frequency taken from that file
    y = noninteracting[:,1]
    plt.plot(x[800:1200], y[800:1200])  #no marker since 2000 points too many for marker   
    #plt.plot(x[400:1600], y2[400:1600], label = "noninteracting", marker = '+')   

    plt.ylabel(r'$A(\omega)$', fontsize=25, rotation=90)
    plt.xlabel(r'$\omega$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(x[0::60], fontsize = 15, rotation = 0)
    plt.legend(fontsize=15)
    plt.grid()
    plt.xlim(-0.8,0.8)
    plt.ylim(0, 15)
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/DOSNonInteracting/DOS.pdf')
    plt.clf()


def plotDOS(NevalinnaCalculationDirectory, A,SelfEnergy, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites, eta):
    x = SelfEnergy[:, 0]        #frequency taken from that file
    y = -A[:]/math.pi
    plt.plot(x[400:1600], y[400:1600])   

    plt.ylabel(r'$A(\omega)$', fontsize=25, rotation=90)
    plt.xlabel(r'$\omega$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.legend(fontsize=15)
    plt.grid()
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/DOS.pdf'.format(u,
                                                                                    beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
    plt.clf()


def plotDOSAll(NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites, Eta):
    for eta in Eta:
        SelfEnergy = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/selfEnergyHartreeForDOS.txt'.format(u,
                                                                beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
        A = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/DOS.dat'.format(u,
                                                                beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
        x = SelfEnergy[:, 0]        #frequency taken from that file
        y = -A[:]/math.pi
        plt.plot(x[400:1600], y[400:1600], label = "$\eta$ = {}".format(eta))   

    plt.ylabel(r'$A(\omega)$', fontsize=25, rotation=90)
    plt.xlabel(r'$\omega$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.ylim(0, 0.3)
    plt.legend(fontsize=15, loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=3, fancybox=False, shadow=False)
    plt.grid()
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/DOSAll.pdf'.format(u,
                                                                                    beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
    plt.clf()



def plotSigmaRenormalized(NevalinnaCalculationDirectory,SelfEnergy, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites, eta):
    x = SelfEnergy[:, 0]        #frequency taken from that file
    y = SelfEnergy[:, 2]
    plt.plot(x[400:1600], y[400:1600])   

    plt.ylabel(r'$Im(\Sigma(\omega))$', fontsize=25, rotation=90)
    plt.xlabel(r'$\omega$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.legend(fontsize=15)
    plt.grid()
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/SigmaRenormalized.pdf'.format(u,
                                                                                    beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
    plt.clf()



def plotSigmaRenormalizedAll(NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites, Eta):
    for eta in Eta:
        SelfEnergy = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/selfEnergyHartreeForDOS.txt'.format(u,
                                                                beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
        x = SelfEnergy[:, 0]        #frequency taken from that file
        y = SelfEnergy[:, 2]
        plt.plot(x[400:1600], y[400:1600], label = "$\eta$ = {}".format(eta))   

    plt.ylabel(r'$Im(\Sigma(\omega))$', fontsize=25, rotation=90)
    plt.xlabel(r'$\omega$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.ylim(-10, 0)
    plt.legend(fontsize=15, loc='upper center', bbox_to_anchor=(0.5, 1.35), ncol=3, fancybox=False, shadow=False)
    plt.grid()
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/SigmaRenormalizedAll.pdf'.format(u,
                                                                                    beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta))
    plt.clf()


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
            sigmaValuesOutput[i,1] = sigmaValuesOutput[i,1] * (u*u *n/2 * (1-n/2))  +  u*n/2       #last term the Hartree term
            np.savetxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/selfEnergyHartreeForDOS.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites, eta), sigmaValuesOutput)
        else:
            pass




'''
plots density n of ED calculations with mu as x
'''
def plotDensityED(NevalinnaCalculationDirectory, u,beta,q,Mu,t,tPri,tPriPri,kSteps,version, bathSites, fillingArray):
    x = Mu[:]       #all mu on x-axis
    y = fillingArray[:]       
    plt.plot(x, y, label = "filling n", marker = 'x')     

    plt.ylabel('n', fontsize=25, rotation=90)
    plt.xlabel(r'$\mu$', fontsize=25, rotation=0)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 0)
    plt.legend(fontsize=15)
    plt.grid()
    plt.tight_layout()
    plt.savefig(NevalinnaCalculationDirectory + '/DensityED_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.pdf'.format(u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.savefig(NevalinnaCalculationDirectory + '/DensityED_U{}_B{}_q{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}.png'.format(u,
                                                                                    beta,q,t,tPri,tPriPri,kSteps,version,bathSites))
    plt.clf()












'''
have different DOS in one graph
insert parameters inside function
'''
def plotAllDOS(MaxEntCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version):
    p=1
    Matrix = createMatrix(p, q, t, tPri, tPriPri, kSteps)

    u = 1.0
    mu = 0.5
    SelfEnergy = np.loadtxt(MaxEntCalculationDirectory + '/MaxEnt_U{}_B{}_q{}_mu{}_t{}_tPri{}_tpriPri{}_kSteps{}_version{}/maxEnt.out.avspec_self.dat'
                                                                            .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    NMaxEnt = len(SelfEnergy)
    A = calculateDOS(q,kSteps,mu, Matrix,SelfEnergy,NMaxEnt)
    np.savetxt(MaxEntCalculationDirectory + "/MaxEnt_U{}_B{}_q{}_mu{}_t{}_tPri{}_tpriPri{}_kSteps{}_version{}/DOS.dat".format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version), A)
    plt.plot(SelfEnergy[:, 0], A[:], label = "U = 1.0")

    u = 2.0
    mu = 1.0
    SelfEnergy = np.loadtxt(MaxEntCalculationDirectory + '/MaxEnt_U{}_B{}_q{}_mu{}_t{}_tPri{}_tpriPri{}_kSteps{}_version{}/maxEnt.out.avspec_self.dat'
                                                                            .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    NMaxEnt = len(SelfEnergy)
    A = calculateDOS(q,kSteps,mu, Matrix,SelfEnergy,NMaxEnt)
    np.savetxt(MaxEntCalculationDirectory + "/MaxEnt_U{}_B{}_q{}_mu{}_t{}_tPri{}_tpriPri{}_kSteps{}_version{}/DOS.dat".format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version), A)
    plt.plot(SelfEnergy[:, 0], A[:], label = "U = 2.0")

    u = 3.0
    mu = 1.5
    SelfEnergy = np.loadtxt(MaxEntCalculationDirectory + '/MaxEnt_U{}_B{}_q{}_mu{}_t{}_tPri{}_tpriPri{}_kSteps{}_version{}/maxEnt.out.avspec_self.dat'
                                                                            .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    NMaxEnt = len(SelfEnergy)
    A = calculateDOS(q,kSteps,mu, Matrix,SelfEnergy,NMaxEnt)
    np.savetxt(MaxEntCalculationDirectory + "/MaxEnt_U{}_B{}_q{}_mu{}_t{}_tPri{}_tpriPri{}_kSteps{}_version{}/DOS.dat".format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version), A)
    plt.plot(SelfEnergy[:, 0], A[:], label = "U = 3.0")


    plt.ylabel(r'$A(\omega)$')
    plt.xlabel(r'$\omega$')
    plt.legend()
    plt.savefig(MaxEntCalculationDirectory + "/AllDOS.png")
    plt.clf()



'''
plots Greens function of nevalinna and different maxent priors for Georgs system
Nevalinna without QMC is the "correct" versions, since only Georgs ED Code and Nevalinna
'''
def plotSpectralAll(QMCCalculationDirectory, NevalinnaCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites):
    MaxEntGreenFlat = np.loadtxt(QMCCalculationDirectory + 
                                        '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/aw_dmft-last_0_0_0.dat'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps))
    plt.plot(MaxEntGreenFlat[:, 0], MaxEntGreenFlat[:, 1], label = "MaxEnt new", linewidth = 0.5)
    NevalinnaGreenWithQMC = np.loadtxt(NevalinnaCalculationDirectory + '/withQMC/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/output.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites)
                                                                            .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    plt.plot(NevalinnaGreenWithQMC[:, 0], NevalinnaGreenWithQMC[:, 1], label = "Nevalinna with QMC", linewidth = 0.5)
    NevalinnaGreenWithoutQMC = np.loadtxt(NevalinnaCalculationDirectory + '/withoutQMC/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/output.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites)
                                                                            .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    plt.plot(NevalinnaGreenWithoutQMC[:, 0], NevalinnaGreenWithoutQMC[:, 1], label = "Nevalinna without QMC", linewidth = 0.5, linestyle = "dashed")
    plt.ylabel(r'$A(\omega)$')
    plt.xlabel(r'$\omega$')
    plt.ylim((0,2))
    plt.legend()
    plt.savefig("/scratch/usr/hhpnhytt/AllSpectralGeorgsSystemU{}_beta{}_5.png".format(u,beta))
    plt.clf()


def plotGreensAll(QMCCalculationDirectory, NevalinnaCalculationDirectory,MaxEntDraftDirectory,MaxEntCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites):
    os.chdir(MaxEntDraftDirectory)                  
    os.system('julia readQMCResultsToGreensMaxEnt.jl {} {} {} {} {} {} {} {} {} {} {}'
                .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version, QMCCalculationDirectory, MaxEntCalculationDirectory))
    #MaxEntGreenFlat = np.loadtxt(QMCCalculationDirectory + 
    #                                    '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Greens.txt'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps))
    #plt.plot(MaxEntGreenFlat[1000:1020, 0], MaxEntGreenFlat[1000:1020, 3], label = "QMC", linewidth = 0.5, linestyle = "dashed", marker = "x")
    #NevalinnaGreenWithQMC = np.loadtxt(NevalinnaCalculationDirectory + '/withQMC/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCode/gm_wim'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites)
    #                                                                        .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    #plt.plot(NevalinnaGreenWithQMC[0:20, 0], NevalinnaGreenWithQMC[0:20, 2], label = "ED with QMC", linewidth = 0.5, linestyle = "None", marker = "x")
    #NevalinnaGreenWithoutQMC = np.loadtxt(NevalinnaCalculationDirectory + '/withoutQMC/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCode/gm_wim'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites)
    #                                                                        .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    #plt.plot(NevalinnaGreenWithoutQMC[0:20, 0], NevalinnaGreenWithoutQMC[0:20, 2], label = "ED without QMC", linewidth = 0.5, linestyle = "dashed", marker = "x")
    #plt.ylabel(r'$G(\nu)$')
    #plt.xlabel(r'$\nu$')
    #plt.legend()
    #plt.savefig("/scratch/usr/hhpnhytt/AllGreensGeorgsSystemU{}_beta{}_5.png".format(u,beta))
    #plt.clf()


def plotGreensDiff(QMCCalculationDirectory, NevalinnaCalculationDirectory,MaxEntDraftDirectory,MaxEntCalculationDirectory, u,beta,q,mu,t,tPri,tPriPri,kSteps,version, bathSites):
    #os.chdir(MaxEntDraftDirectory)                  
    #os.system('julia readQMCResultsToMaxEnt.jl {} {} {} {} {} {} {} {} {} {} {}'
    #            .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version, QMCCalculationDirectory, MaxEntCalculationDirectory))
    MaxEntGreenFlat = np.loadtxt(QMCCalculationDirectory + 
                                        '/finalQMC_U{}_B_{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}/Greens.txt'.format(u,beta,q,mu,t,tPri,tPriPri, kSteps))
    NevalinnaGreenWithQMC = np.loadtxt(NevalinnaCalculationDirectory + '/withQMC/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCode/gm_wim'.format(u,beta,q,mu,t,tPri,tPriPri,kSteps,version,bathSites)
                                                                            .format(u, beta, q, mu, t, tPri, tPriPri, kSteps, version))
    plt.plot(MaxEntGreenFlat[1000:1050, 0], MaxEntGreenFlat[1000:1050, 1] -  NevalinnaGreenWithQMC[0:50, 1], label = "QMC", linewidth = 0.5)

    plt.ylabel(r'$G(\nu)$')
    plt.xlabel(r'$\nu$')
    plt.legend()
    plt.savefig("/scratch/usr/hhpnhytt/DiffGreensGeorgsSystemU{}__beta{}_5.png".format(u,beta))
    plt.clf()
