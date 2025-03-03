#!/bin/bash
#SBATCH -J Nina_Nelson_ODEs_job
#SBATCH -o nelson_odes.txt
#SBATCH -t 00:09:59 	# Time limit (hh:mm:ss)
#SBATCH -n 1       	# Number of tasks
#SBATCH -c 4 		# Number of CPUs per task
#SBATCH -p gpu-a100    	# Queue name (partition, either use normal or gpu-a100)
#SBATCH -A AST24033	# Allocation
#SBATCH -N 1 		# Number of nodes


#Load any modules
#module load tacc-apptainer

# Run your Julia script
apptainer exec julia_latest.sif julia nelson_odes.jl
