#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer_bo
#SBATCH --time 00:10:00
#SBATCH --mem 1gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer_bon.err
#SBATCH -o nucmer_bon.out
#SBATCH --export NONE

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0

for f in BC*.*



nucmer -p ${f} ./Pb_bo_ref.fa  ./${f}

show-coords -rcl ${f}.bon.delta > ${f}.bon.coords

delta-filter -q ${f}.bon.delta > ${f}.bon.filter


done
