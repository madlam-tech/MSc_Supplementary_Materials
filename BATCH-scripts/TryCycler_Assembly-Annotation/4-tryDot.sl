#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryDots
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryDots.err
#SBATCH -o tryDots.out
#SBATCH --export NONE

threads=32
module purge

module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler dotplot  --cluster_dir ./cluster_001
trycycler dotplot  --cluster_dir ./cluster_002
