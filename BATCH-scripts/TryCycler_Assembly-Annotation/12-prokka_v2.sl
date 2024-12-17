#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 00:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka_v2%j.err
#SBATCH -o 12-prokka_v2%j.log
mkdir -p prokka 

module purge
module load prokka/1.14.5-GCC-9.2.0

prokka final.fna \
--force \
--auto \
--prefix prokka \
--addgenes \
--addmrna \
--coverage 70 \

