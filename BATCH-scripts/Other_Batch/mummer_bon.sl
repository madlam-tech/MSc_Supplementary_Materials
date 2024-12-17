#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer_bin
#SBATCH --time 00:45:00
#SBATCH --mem 6gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer.err
#SBATCH -o nucmer.out


module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0


for f in BC*; do

nucmer -p P_bon.fasta ${f}


done
