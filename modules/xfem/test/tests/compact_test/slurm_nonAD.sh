#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=creep
#SBATCH -o %j.o  # stdout; %j expands to jobid
#SBATCH --time=04:30:00
#SBATCH --partition=general
#SBATCH --ntasks-per-node=80
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --wckey=nrc
# Load the necessary modules
module purge
module load use.moose moose-dev-openmpi/2025.04.22

# Run your application and use tee to capture the output
mpiexec -n 80 moose-dev-exec /home/huyngd/projects/moose_AD/modules/xfem/xfem-opt -i CCG_enriched_creep_crack.i ####--distributed-mesh  ####  --split-mesh 112 --split-file ct_h0p4 #
