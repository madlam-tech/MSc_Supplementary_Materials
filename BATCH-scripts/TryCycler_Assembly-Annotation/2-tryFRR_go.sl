#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      2-tryFRR_go
#SBATCH --time          00:02:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 2
#SBATCH --error         2-tryFRR_go%j.err


mkdir -p assemblies


sbatch 2-tryF_set1A.sl 
sbatch 2-tryF_set1B.sl
sbatch 2-tryRR_set2.sl 
sbatch 2-tryRR_set3.sl 
