#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 05:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka_v14%j.err
#SBATCH -o 12-prokka_v14%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0

# Directory containing the fna files (current directory in this case)
input_dir="."

# Create directory for Prokka output on NeSI
output_dir="./prokka_output"
mkdir -p "$output_dir"

# Find all .fna files in the current directory
find "$input_dir" -maxdepth 1 -type f -name '*.fna' | while read -r fna_file; do
    # Extract basename of the file (without path and extension)
    file_basename=$(basename -- "$fna_file" .fna)

    # Create directory for this file's Prokka output in the designated output directory
    mkdir -p "$output_dir/$file_basename"

    # Run Prokka on the fna file
    prokka "$fna_file" \
    --force \
    --outdir "$output_dir/$file_basename" \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70
done

