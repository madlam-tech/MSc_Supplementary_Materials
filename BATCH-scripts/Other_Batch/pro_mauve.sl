#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      progressive_mauve
#SBATCH --time          48:00:00
#SBATCH --mem           35GB
#SBATCH --cpus-per-task 4
#SBATCH --output 	prog_mauve_%j.out
#SBATCH --error		prog_mauve_%j.err


module purge


./progressiveMauve --output Pb_mauve.xmfa --output-guide-tree=Pb_guide_tree.txt --seed-family /nesi/nobackup/massey03345/fnas_4_mauve/*.fna
