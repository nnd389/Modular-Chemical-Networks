#!/bin/bash

#SBATCH -J umist            # Job name
#SBATCH -e umist_out_err_test1.txt
#SBATCH -o umist_out_test1.txt # Name of stdout output file
#SBATCH -p normal           # Queue (partition) name
#SBATCH -N 1                # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 00:20:00         # Run time (hh:mm:ss)
#SBATCH --mail-type=all     # Send email at begin and end of job
#SBATCH -A AST24033         # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=haninadlt@gmail.com

# Any other commands must follow all #SBATCH directives...


# Launch serial code...
apptainer exec julia_latest.sif julia umist_test1.jl      # Do not use ibrun or any other MPI launcher
