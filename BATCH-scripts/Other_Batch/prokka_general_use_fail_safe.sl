#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 05:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e prokka_fs_%j.err
#SBATCH -o prokka_fs_%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0

# Output directory
output_dir="./prokka"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Find all .fna files in the parent directories and process them
find ../* -type f -name "*.fna" | while read -r fna_file; do
    # Check if the file is readable and non-empty
    if [[ ! -r "$fna_file" || ! -s "$fna_file" ]]; then
        echo "Skipping unreadable or empty file: $fna_file"
        continue
    fi

    # Extract basename of the file (without extension)
    file_basename=$(basename "$fna_file" .fna)
    
    # Create a directory for this file's Prokka output in the designated output directory
    mkdir -p "$output_dir/$file_basename"
    
    # Run Prokka on the fna file
    prokka "$fna_file" \
    --force \
    --outdir "$output_dir/$file_basename" \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70 || {
        echo "Prokka failed for $fna_file, continuing with next file."
        continue
    }
done

