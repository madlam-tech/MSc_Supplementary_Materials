#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 01:00:00
#SBATCH --mem 40GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka%j.err
#SBATCH -o 12-prokka%j.log

mkdir -p prokka

module purge
module load prokka/1.14.5-GCC-9.2.0

# Array of specific filenames to process
files_to_process=("ma102_1b.fna")

# Loop through each specified file
for fna_file in "${files_to_process[@]}"; do
    # Check if the file exists in the current directory
    if [ -f "$fna_file" ]; then
        # Extract basename without extension
        file_basename=$(basename "$fna_file" .fna)

        # Create directory for this file's output
        mkdir -p prokka/$file_basename

        # Run Prokka on the .fna file
        prokka "$fna_file" \
        --outdir prokka/ma102_1b \
        --force \
        --prefix ma102_1b \
        --addgenes \
        --addmrna \
        --coverage 70
    else
        echo "File not found: $fna_file"
    fi
done
