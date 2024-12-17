#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 00:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e canu-prokka%j.err
#SBATCH -o canu-prokka%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0

prokka BC11_canu.contigs.fasta
--outdir ./prokka \
--force \
--prefix prokka \
--addgenes \
--addmrna \
--coverage 70 \
