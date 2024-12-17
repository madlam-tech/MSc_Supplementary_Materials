#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer_hay
#SBATCH --time 00:10:00
#SBATCH --mem 1gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer_hay.err
#SBATCH -o nucmer_hay.out
#SBATCH --export NONE

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0

for f in BC*.*;do


nucmer -p ./nucmer_out/${f}.hay ./Pb_ha_ref.fa  ./${f}.hay

show-coords -rcl ./nucmer_out/${f}.hay.delta > ./nucmer_out/${f}.hay.coords

delta-filter -q ./nucmer_out/${f}.hay.delta > ./nucmer_out/${f}.hay.filter


done
