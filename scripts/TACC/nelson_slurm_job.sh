#!/bin/bash
#SBATCH -J Nina_Nelson_job
#SBATCH -o output.txt
#SBATCH -t 00:05:00 	# Time limit (hh:mm:ss)
#SBATCH -n 1       	# Number of tasks
#SBATCH -c 4 		# Number of CPUs per task
#SBATCH -p normal    	# Queue name (partition)
#SBATCH -A CDA23007	# Allocation
#SBATCH -N 2 		# Number of nodes


#Load the julia module?
#module load julia

# Run your Julia script
~/.juliaup/bin/julia nelson.jl
