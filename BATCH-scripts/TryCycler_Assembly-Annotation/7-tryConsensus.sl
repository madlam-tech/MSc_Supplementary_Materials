#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryConsensus
#SBATCH --time 1:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e 7-tryCon%j.err
#SBATCH --export NONE

threads=32


module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler consensus --cluster_dir ./clusters/cluster_001
trycycler consensus --cluster_dir ./clusters/cluster_002
