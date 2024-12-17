#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      distill_DRAM
#SBATCH --time          00:10:00
#SBATCH --mem           2Gb
#SBATCH --cpus-per-task 24
#SBATCH --error         14-DRAM_dist_%x_%A.log
#SBATCH --output        14-DRAM_dist_%x_%A.out

# Load modules
module purge
module load DRAM/1.3.5-Miniconda3



DRAM.py distill -i dram_annotations/annotations.tsv -o dram_distillation --trna_path dram_annotations/trnas.tsv --rrna_path dram_annotations/rrnas.tsv
