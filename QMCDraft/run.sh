#!/bin/bash

module load gcc/8.3.0
module load openmpi/gcc.9/4.1.4

mpif90 hofstadter_hamiltonian.f90 -o hofstadter_hamiltonian.x -g -ffixed-line-length-0
mpirun --oversubscribe ./hofstadter_hamiltonian.x
