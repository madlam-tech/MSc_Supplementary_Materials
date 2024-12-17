#!/bin/bash -e
#SBATCH --account        massey03345
#SBATCH --job-name       flye08_220916
#SBATCH --mem            1G
#SBATCH --cpus-per-task  2
#SBATCH --time           00:10:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e ss07.err
#SBATCH -o ss07.out


module purge
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0

trycycler subsample --thread=16 --genome_size 5.5m --reads test.fq --out_dir read_subsets


