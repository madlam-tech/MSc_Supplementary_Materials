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


nucmer -p ${f} ./Pb_ag_ref.fa  ./${f}

show-coords -rcl ./${f}.delta > ./${f}.coords

delta-filer -q ./${f}.delta > ./${f}.filter


done
