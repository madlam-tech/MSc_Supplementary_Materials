#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      progressive_mauve
#SBATCH --time          36:00:00
#SBATCH --mem           35GB
#SBATCH --cpus-per-task 4
#SBATCH --output        /nesi/nobackup/massey03345/fnas_4_mauve/fnas_4_mauve/prog_mauve_%j.out
#SBATCH --error         /nesi/nobackup/massey03345/fnas_4_mauve/fnas_4_mauve/prog_mauve_%j.err

module purge

# Define the directory containing the input files
input_dir="/nesi/nobackup/massey03345/fnas_4_mauve/fnas_4_mauve"

# Run progressiveMauve with specified input files
./progressiveMauve --output "$input_dir/Pb_mauve.xmfa" --output-guide-tree="$input_dir/Pb_guide_tree.txt" --seed-family "$input_dir"/*.fna

