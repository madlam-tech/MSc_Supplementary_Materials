#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 00:30:00
#SBATCH --mem 40GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka%j.err
#SBATCH -o 12-prokka%j.log

# Define the file to process
fna_file="GCF_000006765.1_ASM676v1_genomic.fna"

# Define the output directory
output_dir="prokka"

# Ensure the prokka directory exists
mkdir -p "$output_dir"

# Load the necessary module
module purge
module load prokka/1.14.5-GCC-9.2.0

# Check if the file exists
if [ -f "$fna_file" ]; then
    # Extract basename without extension
    file_basename=$(basename "$fna_file" .fna)

    # Create directory for this file's output
    mkdir -p "$output_dir/$file_basename"

    # Run Prokka on the .fna file
    prokka "$fna_file" \
    --outdir "$output_dir/$file_basename" \
    --force \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70
else
    echo "File $fna_file does not exist."
    exit 1
fi

# Feedback
echo "Processing complete. Output saved to $output_dir/$file_basename."
