#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      progressive_mauve
#SBATCH --time          18:00:00
#SBATCH --mem           8GB
#SBATCH --cpus-per-task 2
#SBATCH --output 	prog_mauve_Pb_%j.out
#SBATCH --error		prog_mauve_Pb_%j.err


module purge


./progressiveMauve --output ma_Pb_mauve.xmfa --output-guide-tree=ma_Pb_guide_tree.txt --seed-family /nesi/nobackup/massey03345/fnas_4_mauve/ma*.fna
