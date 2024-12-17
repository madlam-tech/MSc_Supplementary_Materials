#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      annotate_DRAM
#SBATCH --time          5:00:00
#SBATCH --mem           30Gb
#SBATCH --cpus-per-task 24
#SBATCH --error         14-DRAM_%x_%A.log
#SBATCH --output        14-DRAM_%x_%A.out

# Load modules
module purge
module load DRAM/1.3.5-Miniconda3

# Run DRAM
DRAM.py annotate -i final.fna \
                 --checkm_quality ./checkM/checkm.txt \
                 -o dram_annotations --threads $SLURM_CPUS_PER_TASK
