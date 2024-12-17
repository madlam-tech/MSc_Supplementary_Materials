#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J pseudofinder
#SBATCH --time 16:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e pseudofinder_%j.err
#SBATCH -o pseudofinder_%j.log

module load DIAMOND/2.1.6-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-04
module load Python/3.10.5-gimkl-2022a

# Ensure the directory structure is clean before starting
rm -rf pseudogene_output

# Create a directory to store output for all pseudogene files
mkdir pseudogene_output

# Find .gbk files in the specified directory and its subdirectories
find /nesi/nobackup/massey03345/fnas_4_annotation/prokka_4_pseudo -type f -name "*.gbk" | while read -r filename; do
    # Extract basename without extension
    basename=$(basename "$filename" .gbk)
    
    # Create a directory for each file processed
    mkdir -p "pseudogene_output/$basename"

    # Run pseudofinder.py for the current file
    ./pseudofinder.py annotate --diamond --skip_makedb -g "$filename" -db /nesi/nobackup/massey03345/uniprot/uniprot_sprot.dmnd -op "pseudogene_output/$basename"
done

