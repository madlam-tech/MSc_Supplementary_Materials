#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mummer_ag
#SBATCH --time 00:45:00
#SBATCH --mem 6gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e nucmer.err
#SBATCH -o nucmer.out
#SBATCH --export NONE

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0


mummer -maxmatch -l 20 -b -n -k 3 -threads 3 P_ag.fa  BC.fa
