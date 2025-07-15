#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=ADcreep
#SBATCH -o %j.o  # stdout; %j expands to jobid
#SBATCH --time=13:00:00
#SBATCH --partition=hbm
#SBATCH --ntasks-per-node=100
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --wckey=nrc
# Load the necessary modules
module purge
module load use.moose moose-dev-openmpi/2025.04.22

# Run your application and use tee to capture the output
mpiexec -n 100 moose-dev-exec /home/huyngd/projects/moose_AD/modules/xfem/xfem-opt -i ADCCG_enriched_creep_crack_h0p25.i ####--distributed-mesh  ####  --split-mesh 112 --split-file ct_h0p4 #
