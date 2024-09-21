#!/bin/bash

for U in 1.0 2.0 3.0
do
   for beta in 30.0 50.0 70.0 100.0
   do 
      for q in 3 4
      do
         for tPrime 0.25
            cd /home/hhpnhytt/QMC/QMC_U($U)_B($beta)_P1_L($q)_tPrime($tPrime)
            #find file of first calc and find Ncorr with hdf5
            filename=$(find -name '*.hdf5')
            echo $filename
            Ncorr=$(julia /scratch/projects/hhp00048/codes/scripts/LadderDGA_utils/ncorr.jl $filename)

            #change Parameters.in to production calculation style (more DMFT steps and Nmeas = 10^6)
            #change the submit file to production calculation style (5 nodes and 480 cores)
            #also change Ncorr found before
            python3 change Ncorr.py $Ncorr
            python3 changeToProductionStyle.py

            #last step is to submit it
            sbatch submit_dmft_script 
         done
      done
   done
done



