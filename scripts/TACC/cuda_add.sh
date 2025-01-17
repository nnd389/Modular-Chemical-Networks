#!/bin/bash

# This is a SLURM submission script for cuda_add.jl

#SBATCH -J cuda_job       # Job name
#SBATCH -o job_output.txt # Output file
#SBATCH -n 1              # Tasks
#SBATCH -N 1              # Nodes
#SBATCH -p normal         # Parition: Normal queue
#SBATCH -t 00:05:00       # Time

# Load the required modules
module load cuda

# Run your Julia script
julia cuda_add.jl

