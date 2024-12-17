#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J pseudo
#SBATCH --time 02:30:00
#SBATCH --mem 4GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 22-pseudo_%j.log


mkdir -p psuedogenes 

module purge
module load Python/3.9.9-gimkl-2020a
module load DIAMOND/2.1.1-GCC-11.3.0


# Define input and output directories
input_dir="./prokka-2023-05-06_v9"
output_dir="./pseudogenes"

# Loop through input files and run PseudoFinder on each one
for file in "${input_dir}"/*.gbk; do
    # Get filename without extension
    filename=$(basename "${file%.*}")

pseudofinder.py annotate --diamond 
