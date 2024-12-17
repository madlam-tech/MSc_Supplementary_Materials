#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryClu_v2
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryClu_v2.err
#SBATCH -o tryClu_v2.out
#SBATCH --export NONE

module purge

module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler cluster --distance 0.04 --threads 32 -—assemblies assemblies/*.fasta --reads ma1263_1.fq --out_dir tryclust_v2
