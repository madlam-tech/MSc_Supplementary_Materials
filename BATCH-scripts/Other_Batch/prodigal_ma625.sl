#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          00:20:00
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

    # Loop through .fasta files in the current directory
    for filename in "$dir"/*.fasta; do
        # Check if the file exists
        if [ -f "$filename" ]; then
            # Extract basename without extension
            basename=$(basename "$filename" .fasta)

            # Create a directory for each file processed
            mkdir -p "$dir/prodigal"

            # Run Prodigal for the current file
            prodigal -i "$filename" \
                     -o "$dir/prodigal/${basename}_genes.fa" \
                     -a "$dir/prodigal/${basename}_proteins.faa" \
                     -p single
        fi
    done

    # Recursively process subdirectories
    for subdir in "$dir"/*/; do
        process_directory "$subdir"
    done
}

# Start processing from the current directory
process_directory "$(pwd)"


