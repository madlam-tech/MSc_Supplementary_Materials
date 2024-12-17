#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 02:30:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka%j.err
#SBATCH -o 12-prokka%j.log

mkdir -p prokka

module purge
module load prokka/1.14.5-GCC-9.2.0

# Function to process directories recursively
process_directory() {
    local dir="$1"

    # Loop through .fna files in the current directory
    for fna_file in "$dir"/*.fna; do
        # Check if the file exists
        if [ -f "$fna_file" ]; then
            # Extract basename without extension
            file_basename=$(basename "$fna_file" .fna)

            # Check if the output directory for this file already exists and has the expected output
            if [ -d "prokka/$file_basename" ] && [ -f "prokka/$file_basename/$file_basename.gff" ]; then
                echo "Output for $fna_file already exists, skipping..."
                continue
            fi

            # Create directory for this file's output
            mkdir -p prokka/$file_basename

            # Run Prokka on the .fna file
            prokka "$fna_file" \
            --outdir prokka/$file_basename \
            --force \
            --prefix $file_basename \
            --addgenes \
            --addmrna \
            --coverage 70
        fi
    done

    # Recursively process subdirectories
    for subdir in "$dir"/*/; do
        if [ -d "$subdir" ]; then
            process_directory "$subdir"
        fi
    done
}

# Start processing from the main directory
process_directory "$(pwd)"
