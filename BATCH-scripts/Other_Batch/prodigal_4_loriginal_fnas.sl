#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          01:00:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1
#SBATCH --error         prod%j.err
#SBATCH --output        prod%j.out

# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

# Function to process directories recursively
process_directory() {
    local dir="$1"

    # Loop through .fna files in the current directory
    for filename in "$dir"/*.fna; do
        # Check if the file exists
        if [ -f "$filename" ]; then
            # Extract basename without extension
            basename=$(basename "$filename" .fna)

            # Create a directory for each file processed
            mkdir -p "$dir/prodigal_original"

            # Run Prodigal for the current file
            prodigal -i "$filename" \
                     -o "$dir/prodigal_original/${basename}_genes.fa" \
                     -a "$dir/prodigal_original/${basename}_proteins.faa" \
                     -p single
        fi
    done

    # Recursively process subdirectories
    for subdir in "$dir"/*/; do
        process_directory "$subdir"
    done
}

# Start processing from the specified directory
process_directory "/nesi/nobackup/massey03345/fnas_4_annotation/fnas"

