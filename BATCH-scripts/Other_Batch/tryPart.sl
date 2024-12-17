#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryPart_v1
#SBATCH --time 4:00:00
#SBATCH --mem 2GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 2
#SBATCH -e tryPart_v1.err
#SBATCH -o tryPart_v1.out
#SBATCH --export NONE

threads=32

module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler partition --reads ../*fq --cluster_dirs ./cluster_*
