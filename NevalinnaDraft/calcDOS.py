import numpy as np
import math
import cmath
import re


kStepsFile = 240


p = 1
q = 3
kSteps = 240
u = 
beta = 
version = 
bathSites = 
t = 
tPri = 
tPriPri = 
mu = 
eta = 
NevalinnaCalculationDirectory = 


'''
get filling of ED, search in run.out for densimp
'''
def getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, kStepsFile):
    with open(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/EDCodeCorrect/run.out'.format(u,beta,q,mu,t,tPri,tPriPri,kStepsFile,version,bathSites), 'r', encoding='utf-8') as file:
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

    print("drin")
    for w in range(NNevalinna):  
        if (SelfEnergy[w, 0] > -3) and (SelfEnergy[w, 0] < 3):  
            for ikx in range(kSteps):
                for iky in range(kSteps):  
                    for i in range(q):
                        for j in range(q):
                            #TODO: sit das wirklich richtig so?
                            denominator[i,j] = (SelfEnergy[w, 0] + mu - (SelfEnergy[w, 1] + imag*SelfEnergy[w,2]))*np.identity(q)[i,j] - Matrix[ikx,iky, i,j]          
                    insideMatrix = np.linalg.inv(denominator) * 1/(kSteps*kSteps * 2*math.pi)                 #Matrix inversion and prefactors to 

                    A[w] = A[w] + np.imag(np.trace(insideMatrix))
            print("mu = {} with w= {} \n".format(mu,w))
            #print(A[w])
        else: 
            pass

    return (A/q)



'''
since only self energy looking like Greens function continued
recalculate it to normal looking self energy (with hartree term when HartreeTerm = True)
HartreeTerm = False only useful when added later (e.g. in Hall Code)
'''
def normalizeNevalinnaOutputSelfEnergy(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, n, HartreeTerm, eta, kStepsFile):
    sigmaValuesOutput = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/output.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kStepsFile,version,bathSites, eta))    
    for i in range(len(sigmaValuesOutput)):
        sigmaValuesOutput[i,2] = sigmaValuesOutput[i,2] * (u*u *n/2 * (1-n/2))
        if HartreeTerm:
            sigmaValuesOutput[i,1] = sigmaValuesOutput[i,1] * (u*u *n/2 * (1-n/2))  +  u*n/2       #last term the Hartree term
            np.savetxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/selfEnergyHartreeForDOS.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kStepsFile,version,bathSites, eta), sigmaValuesOutput)
        else:
            pass


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




















Matrix = createMatrix(p, q, t, tPri, tPriPri, kSteps)

filling = getFillingED(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, kStepsFile)

normalizeNevalinnaOutputSelfEnergy(NevalinnaCalculationDirectory, u,beta,mu,q,t,tPri,tPriPri,kSteps, version, bathSites, filling, True, eta, kStepsFile)
SelfEnergy = np.loadtxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/selfEnergyHartreeForDOS.txt'.format(u,beta,q,mu,t,tPri,tPriPri,kStepsFile,version,bathSites, eta))
NNevalinna = len(SelfEnergy)
A = calculateDOS(q,kSteps,mu, Matrix,SelfEnergy,NNevalinna)
np.savetxt(NevalinnaCalculationDirectory + '/Nevalinna_U{}_B{}_q{}_mu{}_t{}_tPri{}_tPriPri{}_kSteps{}_version{}_bathSites{}/NevalinnaEta{}/DOS.dat'.format(u,beta,q,mu,t,tPri,tPriPri,kStepsFile,version,bathSites, eta), A)









