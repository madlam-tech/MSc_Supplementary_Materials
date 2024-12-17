#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryClust 
#SBATCH --time 4:00:00
#SBATCH --mem 12GB
#SBATCH --cpus-per-task 16
#SBATCH -e tryClust_%j.err
#SBATCH -o tryClust_%j.out


module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler cluster --threads $SLURM_CPUS_PER_TASK --assemblies assemblies/*.fasta --reads *.fq --out_dir tryClust
