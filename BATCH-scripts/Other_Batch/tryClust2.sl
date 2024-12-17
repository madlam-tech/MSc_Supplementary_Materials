#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryClust 
#SBATCH --time 4:00:00
#SBATCH --mem 12GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryClust.err
#SBATCH -o tryClust.out
#SBATCH --mail-user=madlam@me.com
#SBATCH --mail-type=FAIL
#SBATCH --export NONE

threads=32

module purge
module load Mash/2.3-GCC-9.2.0
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler cluster --assemblies assemblies/*.fasta --reads *.fq --out_dir tryClust2

