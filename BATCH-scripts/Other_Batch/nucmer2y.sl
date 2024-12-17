#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer_ag
#SBATCH --time 00:10:00
#SBATCH --mem 1gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer.err
#SBATCH -o nucmer.out
#SBATCH --export NONE

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0

for f in BC*.*;do


nucmer -p ${f} ./P_ag.fa  ./${f}

show-coords -rcl ${f}.delta > ${f}.coords

delta-filter -q ${f}.delta > ${f}.filter

show-tiling -q ${f} > ${f}.tiling

done
