#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J pseudofinder
#SBATCH --time 8:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e pseudofinder_%j.err
#SBATCH -o pseudofinder_%j.log

# Load required modules
module load DIAMOND/2.1.6-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-04
module load Python/3.10.5-gimkl-2022a

# Set the input directory
input_dir="/nesi/nobackup/massey03345/fnas_4_annotation/MLST_genes/prokka"

# Create a datestamped output directory
timestamp=$(date +%Y%m%d_%H%M%S)
output_dir="/nesi/nobackup/massey03345/fnas_4_annotation/MLST_genes/pseudogene_output_$timestamp"

# Ensure the directory structure is clean before starting
mkdir -p "$output_dir"

# Find .gbk files in the specified directory and its subdirectories
find "$input_dir" -type f -name "*.gbk" | while read -r filename; do
    # Extract basename without extension
    basename=$(basename "$filename" .gbk)
    
    # Create a directory for each file processed
    mkdir -p "$output_dir/$basename"

    # Run pseudofinder.py for the current file
    ./pseudofinder.py annotate --diamond --skip_makedb -g "$filename" -db /nesi/nobackup/massey03345/uniprot/uniprot_sprot.dmnd -op "$output_dir/$basename"
done

