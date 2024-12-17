#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer3
#SBATCH --time 00:10:00
#SBATCH --mem 1gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer3.err
#SBATCH -o nucmer3.out
#SBATCH --export NONE

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0

for f in BC*.*;do


nucmer -p ./nucmer_out/${f} ./Pb_ag_ref.fa  ./${f}

show-coords -rcl ./nucmer_out/${f}.delta > ./nucmer_out/${f}.coords

delta-filter -q ./nucmer_out/${f}.delta > ./nucmer_out/${f}.filter


done
