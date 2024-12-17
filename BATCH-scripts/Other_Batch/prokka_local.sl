#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 00:30:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka%j.err
#SBATCH -o 12-prokka%j.log

mkdir -p prokka_2

module purge
module load prokka/1.14.5-GCC-9.2.0

# Function to process .fna files in the specified directory
process_directory() {
    local dir="$1"

    # Loop through .fna files in the current directory
    for fna_file in "$dir"/*.fna; do
        # Check if the file exists
        if [ -f "$fna_file" ]; then
            # Extract basename without extension
            file_basename=$(basename "$fna_file" .fna)

            # Create directory for this file's output
            mkdir -p prokka_2/$file_basename

            # Run Prokka on the .fna file
            prokka --compliant "$fna_file" \
            --outdir prokka_2/$file_basename \
            --force \
            --prefix $file_basename \
            --addgenes \
            --addmrna \
            --coverage 70
        fi
    done
}

# Start processing from the main directory (current directory)
process_directory "$(pwd)"
