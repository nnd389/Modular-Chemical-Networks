#!/bin/bash
#SBATCH -J Nina_Nelson_ODEs_job
#SBATCH -o nelson_catalyst.txt
#SBATCH -t 00:02:00 	# Time limit (hh:mm:ss)
#SBATCH -n 1       	# Number of tasks
#SBATCH -c 4 		# Number of CPUs per task
#SBATCH -p normal    	# Queue name (partition)
#SBATCH -A AST24033	# Allocation
#SBATCH -N 1 		# Number of nodes


#Load any modules

# Run your Julia script
apptainer exec julia_latest.sif julia nelson_catalyst.jl
