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



fileName = string(QMCCalculationDirectory, "/finalQMC_U$(U)_B_$(Beta)_q$(q)_mu$(mu)_t$(t)_tPri$(tPrime)_tPriPri$(tPrimePrime)_kSteps$(kSteps)/hdf5BackupFolder/U$(U)_mu$(mu)_B$(Beta)_q$(q)_t$(t)_tPri$(tPrime)_tPriPri$(tPrimePrime)_kSteps$(kSteps)_DMFT_Version$(version).hdf5")
f = h5open(fileName, "r")
Values = read(f["dmft-last/ineq-001/occ/value"])    
Error = read(f["dmft-last/ineq-001/occ/error"])            







occupation1 = Values[:,:,1,1]
error1 = Error[:,:,1,1]


outputNameValues = string(QMCCalculationDirectory, "/finalQMC_U$(U)_B_$(Beta)_q$(q)_mu$(mu)_t$(t)_tPri$(tPrime)_tPriPri$(tPrimePrime)_kSteps$(kSteps)/densityValues_QMC_$(version).csv")
outputNameErrors = string(QMCCalculationDirectory, "/finalQMC_U$(U)_B_$(Beta)_q$(q)_mu$(mu)_t$(t)_tPri$(tPrime)_tPriPri$(tPrimePrime)_kSteps$(kSteps)/densityErrors_QMC_$(version).csv")
writedlm(outputNameValues, Values)
writedlm(outputNameErrors, Error)

