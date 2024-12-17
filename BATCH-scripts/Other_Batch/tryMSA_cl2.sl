#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryMSA_v1
#SBATCH --time 8:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryMSA_v1.err
#SBATCH -o tryMSA_v1.out
#SBATCH --export NONE

threads=32


module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

#trycycler msa --cluster_dir ./cluster_001
trycycler msa --cluster_dir ./cluster_002
#trycycler msa --cluster_dir ./cluster_003
#trycycler msa --cluster_dir ./cluster_004
#trycycler msa --cluster_dir ./cluster_005
