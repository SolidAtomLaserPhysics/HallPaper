using HDF5
using DelimitedFiles



U = ARGS[1]
Beta = ARGS[2]
q = ARGS[3]
mu = ARGS[4]
t = ARGS[5]
tPrime = ARGS[6]
tPrimePrime = ARGS[7]
kSteps = ARGS[8]
version = ARGS[9]
QMCCalculationDirectory = ARGS[10]
it = ARGS[11]

fileName = string(QMCCalculationDirectory, "/finalQMC_U$(U)_B_$(Beta)_q$(q)_mu$(mu)_t$(t)_tPri$(tPrime)_tPriPri$(tPrimePrime)_kSteps$(kSteps)/hdf5BackupFolder/U$(U)_mu$(mu)_B$(Beta)_q$(q)_t$(t)_tPri$(tPrime)_tPriPri$(tPrimePrime)_kSteps$(kSteps)_DMFT_Version$(version).hdf5")



f = h5open(fileName, "r")
if parse(Int64, it) < 10
    Values = read(f["dmft-00$(it)/ineq-001/siw-full/value"])
    Errors = read(f["dmft-00$(it)/ineq-001/siw-full/error"])
else
    Values = read(f["dmft-0$(it)/ineq-001/siw-full/value"])
    Errors = read(f["dmft-0$(it)/ineq-001/siw-full/error"])
end

Frequencies = read(f[".axes/iw"])

siw = Values[:,1,1,1,1]
error = Errors[:,1,1,1,1]

#only reads out 1001 till 1201 which are the first 200 points with posiitve Matsubara frequency
RealPart = [real(siw[1001:1201])]                                #apparently Real Part is 2D array but only 1 element in that one direction. In reality vector
ImagPart = [imag(siw[1001:1201])]
RealError = [real(error[1001:1201])]
ImagError = [imag(error[1001:1201])]
MatsubaraFreq = Frequencies[1001:1201]




vv = [MatsubaraFreq RealPart[1] ImagPart[1] RealError[1] ImagError[1]]

outputName = string(QMCCalculationDirectory, "/finalQMC_U$(U)_B_$(Beta)_q$(q)_mu$(mu)_t$(t)_tPri$(tPrime)_tPriPri$(tPrimePrime)_kSteps$(kSteps)/self-en_wim_QMC_$(version)_iteration$(it).csv")

writedlm(outputName, vv)

