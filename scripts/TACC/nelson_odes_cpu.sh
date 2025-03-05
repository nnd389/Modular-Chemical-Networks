#!/bin/bash
#SBATCH -J Nelson_odes_cpu
#SBATCH -o nelson_odes_cpu.txt
#SBATCH -t 00:04:49 	# Time limit (hh:mm:ss)
#SBATCH -n 1       	# Number of tasks
#SBATCH -c 4 		# Number of CPUs per task
#SBATCH -p normal    	# Queue name (partition, either use normal or gpu-a100)
#SBATCH -A AST24033	# Allocation
#SBATCH -N 1 		# Number of nodes


#Load any modules
module load tacc-apptainer

#Set Julia threading
export JULIA_NUM_THREADS=4

# Run your Julia script
apptainer exec julia_latest.sif julia nelson_odes_cpu.jl
