#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer_self
#SBATCH --time 00:10:00
#SBATCH --mem 1gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer.err
#SBATCH -o nucmer.out
#SBATCH --export NONE

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0


nucmer -p self ./Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fna  ./Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fna

show-coords -rcl ./nucmer_out/Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fna.delta > ./nucmer_out/Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fna.coords

delta-filer -q ./Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fna.delta > ./nucmer_out/./Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fna.filter

