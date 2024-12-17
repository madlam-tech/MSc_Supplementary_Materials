#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer_agvbo
#SBATCH --time 00:10:00
#SBATCH --mem 1gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer2_ag_bo.err
#SBATCH -o nucmer_ag_bo.out
#SBATCH --export NONE

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0


nucmer -p ./nucmer_out/ag_v_bo ./Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fa  ./Pb_bo_GCF_009455625.1_ASM945562v1_genomic.fa

show-coords -rcl ./nucmer_out/ag_v_bo.delta > ./nucmer_out/ag_v_bo.coords

delta-filter -q ./nucmer_out/ag_v_bo.delta > ./nucmer_out/ag_v_bo.filter

