#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 03:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka_v2%j.err
#SBATCH -o 12-prokka_v2%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0

# Directory containing the fna files
input_dir="/nesi/nobackup/massey03345/fnas_4_annotation/staph"

# Create directory for Prokka output on NeSI
output_dir="/nesi/nobackup/massey03345/fnas_4_annotation/prokka_staph"
mkdir -p "$output_dir"

# Navigate to the directory containing fna files
cd "$input_dir"

# Find all genomic.fna files across subdirectories
find . -type f -name '*_genomic.fna' -exec basename {} \; | while read fna_file; do
    # Extract basename of the file (without extension)
    file_basename=$(basename -- "$fna_file" .fna)
    
    # Create directory for this file's Prokka output in the designated output directory
    mkdir -p "$output_dir/$file_basename"
    
    # Run Prokka on the fna file
    prokka "$input_dir/$file_basename/$fna_file" \
    --force \
    --outdir "$output_dir/$file_basename" \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70
done

