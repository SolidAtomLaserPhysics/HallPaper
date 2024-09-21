from telnetlib import TM
import numpy as np
import math
import cmath




def writeHKFile(p, q, t, tPri, tPriPri, kSteps, filename):
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
    build HK.dat file
    '''
    with open(filename, "w") as file:
        file.write("{}   {}   {} \n".format( kSteps*kSteps, q, q))






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


                



    with open(filename, "a+") as file:
        if (q == 3):
            for ikx in range(kSteps):
                for iky in range(kSteps):
                    kx = -math.pi + ikx * dkx
                    ky = -math.pi/q + iky * dky

                    file.write("   {} {} {} \n".format(kx, ky, 0))
                    file.write("{} {} {} {} {} {} \n".format(np.real(epsMatrix[ikx, iky, 0,0]), np.imag(epsMatrix[ikx, iky, 0,0]), np.real(epsMatrix[ikx, iky, 0,1]), np.imag(epsMatrix[ikx, iky, 0,1]), np.real(epsMatrix[ikx, iky, 0,2]), np.imag(epsMatrix[ikx, iky, 0,2])))
                    file.write("{} {} {} {} {} {} \n".format(np.real(epsMatrix[ikx, iky, 1,0]), np.imag(epsMatrix[ikx, iky, 1,0]), np.real(epsMatrix[ikx, iky, 1,1]), np.imag(epsMatrix[ikx, iky, 1,1]), np.real(epsMatrix[ikx, iky, 1,2]), np.imag(epsMatrix[ikx, iky, 1,2])))
                    file.write("{} {} {} {} {} {} \n".format(np.real(epsMatrix[ikx, iky, 2,0]), np.imag(epsMatrix[ikx, iky, 2,0]), np.real(epsMatrix[ikx, iky, 2,1]), np.imag(epsMatrix[ikx, iky, 2,1]), np.real(epsMatrix[ikx, iky, 2,2]), np.imag(epsMatrix[ikx, iky, 2,2])))




    '''
    build values 
    '''
    if (q == 4):
        for ikx in range(kSteps):
            for iky in range(kSteps):
                kx = -math.pi + ikx * dkx
                ky = -math.pi/q + iky * dky

                ''''
                same for a 4x4 Matrix == q = 4
                '''
            
                '''
                T Matrix
                '''
                TMatrix[ikx, iky, 0,0] = 2*t*math.cos(kx)
                TMatrix[ikx, iky, 0,1] = t
                TMatrix[ikx, iky, 0,2] = 0
                TMatrix[ikx, iky, 0,3] = t*cmath.exp(j*ky*q)

                TMatrix[ikx, iky, 1,0] = t
                TMatrix[ikx, iky, 1,1] = 2*t*math.cos(kx + (2*math.pi*B))
                TMatrix[ikx, iky, 1,2] = t
                TMatrix[ikx, iky, 1,3] = 0

                TMatrix[ikx, iky, 2,0] = 0
                TMatrix[ikx, iky, 2,1] = t
                TMatrix[ikx, iky, 2,2] = 2*t*math.cos(kx + (4*math.pi*B))
                TMatrix[ikx, iky, 2,3] = t

                TMatrix[ikx, iky, 3,0] = t*cmath.exp(-j*ky*q)
                TMatrix[ikx, iky, 3,1] = 0
                TMatrix[ikx, iky, 3,2] = t
                TMatrix[ikx, iky, 3,3] = 2*t*math.cos(kx + (6*math.pi*B))

                '''
                TPrime Matrix
                '''
                TPriMatrix[ikx, iky, 0,0] = 0
                TPriMatrix[ikx, iky, 0,1] = 2*tPri * math.cos(kx + (math.pi*B))
                TPriMatrix[ikx, iky, 0,2] = 0
                TPriMatrix[ikx, iky, 0,3] = 2*tPri * math.cos(kx - (math.pi*B)) * cmath.exp(j*ky*q)

                TPriMatrix[ikx, iky, 1,0] = 2*tPri * math.cos(kx + (math.pi*B))
                TPriMatrix[ikx, iky, 1,1] = 0
                TPriMatrix[ikx, iky, 1,2] = 2*tPri * math.cos(kx + ((3)*math.pi*B))
                TPriMatrix[ikx, iky, 1,3] = 0

                TPriMatrix[ikx, iky, 2,0] = 0
                TPriMatrix[ikx, iky, 2,1] = 2*tPri * math.cos(kx + ((3)*math.pi*B))
                TPriMatrix[ikx, iky, 2,2] = 0
                TPriMatrix[ikx, iky, 2,3] = 2*tPri * math.cos(kx + ((5)*math.pi*B))

                TPriMatrix[ikx, iky, 3,0] = 2*tPri * math.cos(kx - (math.pi*B)) * cmath.exp(-j*ky*q)
                TPriMatrix[ikx, iky, 3,1] = 0
                TPriMatrix[ikx, iky, 3,2] = 2*tPri * math.cos(kx + ((5)*math.pi*B))
                TPriMatrix[ikx, iky, 3,3] = 0


                '''
                TPrimePrime Matrix
                '''
                TPriPriMatrix[ikx, iky, 0,0] = 2*tPriPri*math.cos(2*kx)
                TPriPriMatrix[ikx, iky, 0,1] = 0
                TPriPriMatrix[ikx, iky, 0,2] = tPriPri + tPriPri*cmath.exp(j*ky*q)
                TPriPriMatrix[ikx, iky, 0,3] = 0

                TPriPriMatrix[ikx, iky, 1,0] = 0
                TPriPriMatrix[ikx, iky, 1,1] = 2*tPriPri*math.cos(2*kx + (4*math.pi*B))
                TPriPriMatrix[ikx, iky, 1,2] = 0
                TPriPriMatrix[ikx, iky, 1,3] = tPriPri + tPriPri*cmath.exp(j*ky*q)

                TPriPriMatrix[ikx, iky, 2,0] = tPriPri + tPriPri*cmath.exp(-j*ky*q)
                TPriPriMatrix[ikx, iky, 2,1] = 0
                TPriPriMatrix[ikx, iky, 2,2] = 2*tPriPri*math.cos(2*kx + (8*math.pi*B))
                TPriPriMatrix[ikx, iky, 2,3] = 0

                TPriPriMatrix[ikx, iky, 3,0] = 0
                TPriPriMatrix[ikx, iky, 3,1] = tPriPri + tPriPri*cmath.exp(-j*ky*q)
                TPriPriMatrix[ikx, iky, 3,2] = 0
                TPriPriMatrix[ikx, iky, 3,3] = 2*tPriPri*math.cos(2*kx + (12*math.pi*B))

                epsMatrix = TMatrix + TPriMatrix + TPriPriMatrix


                



    with open(filename, "a+") as file:
        if (q == 4):
            for ikx in range(kSteps):
                for iky in range(kSteps):
                    kx = -math.pi + ikx * dkx
                    ky = -math.pi/q + iky * dky

                    file.write("   {} {} {} \n".format(kx, ky, 0))
                    file.write("{} {} {} {} {} {} {} {} \n".format(np.real(epsMatrix[ikx, iky, 0,0]), np.imag(epsMatrix[ikx, iky, 0,0]), np.real(epsMatrix[ikx, iky, 0,1]), np.imag(epsMatrix[ikx, iky, 0,1]), np.real(epsMatrix[ikx, iky, 0,2]), np.imag(epsMatrix[ikx, iky, 0,2]), np.real(epsMatrix[ikx, iky, 0,3]), np.imag(epsMatrix[ikx, iky, 0,3])))
                    file.write("{} {} {} {} {} {} {} {} \n".format(np.real(epsMatrix[ikx, iky, 1,0]), np.imag(epsMatrix[ikx, iky, 1,0]), np.real(epsMatrix[ikx, iky, 1,1]), np.imag(epsMatrix[ikx, iky, 1,1]), np.real(epsMatrix[ikx, iky, 1,2]), np.imag(epsMatrix[ikx, iky, 1,2]), np.real(epsMatrix[ikx, iky, 1,3]), np.imag(epsMatrix[ikx, iky, 1,3])))
                    file.write("{} {} {} {} {} {} {} {} \n".format(np.real(epsMatrix[ikx, iky, 2,0]), np.imag(epsMatrix[ikx, iky, 2,0]), np.real(epsMatrix[ikx, iky, 2,1]), np.imag(epsMatrix[ikx, iky, 2,1]), np.real(epsMatrix[ikx, iky, 2,2]), np.imag(epsMatrix[ikx, iky, 2,2]), np.real(epsMatrix[ikx, iky, 2,3]), np.imag(epsMatrix[ikx, iky, 2,3])))
                    file.write("{} {} {} {} {} {} {} {} \n".format(np.real(epsMatrix[ikx, iky, 3,0]), np.imag(epsMatrix[ikx, iky, 3,0]), np.real(epsMatrix[ikx, iky, 3,1]), np.imag(epsMatrix[ikx, iky, 3,1]), np.real(epsMatrix[ikx, iky, 3,2]), np.imag(epsMatrix[ikx, iky, 3,2]), np.real(epsMatrix[ikx, iky, 3,3]), np.imag(epsMatrix[ikx, iky, 3,3])))




