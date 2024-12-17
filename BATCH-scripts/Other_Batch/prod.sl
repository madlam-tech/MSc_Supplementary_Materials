#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal
#SBATCH --time          00:20:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1
#SBATCH --error         prodigal.err
#SBATCH --output        prodigal.out

module purge
module load prodigal/2.6.3-GCC-11.3.0


    prodigal -i BC04_2022-09-23.fasta -p single \
             -d predictions/BC04.genes.fna \
             -a predictions/BC04.genes.faa \
             -o predictions/BC04.genes.gbk
