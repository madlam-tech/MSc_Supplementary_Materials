#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryClust 
#SBATCH --time 04:00:00
#SBATCH --mem 24GB
#SBATCH --cpus-per-task 32
#SBATCH -e 3-tryClust_%j.err
#SBATCH -o 3-tryClust_%j.out


mkdir -p tryCluster
module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler cluster --threads $SLURM_CPUS_PER_TASK --distance 0.02 --assemblies assemblies/*.fasta --reads ../filtered/*.fq --out_dir tryCluster
