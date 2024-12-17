#!/bin/bash -e
#SBATCH --account massey03345 
#SBATCH --job-name=GPUJob   # job name (shows up in the queue)
#SBATCH --time=00-00:10:00  # Walltime (DD-HH:MM:SS)
#SBATCH --gpus-per-node=1   # GPU resources required per node
#SBATCH --cpus-per-task=2   # number of CPUs per task (1 by default)
#SBATCH --mem=512MB         # amount of memory per node (1 by default)

# load CUDA module
module purge
module load CUDA/11.0.2

# display information about the available GPUs
nvidia-smi

# check the value of the CUDA_VISIBLE_DEVICES variable
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
