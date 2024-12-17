#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          00:20:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1

# Get the current date in YYYYMMDD format
current_date=$(date +%Y%m%d)

#SBATCH --error         prod_${current_date}_%j.err
#SBATCH --output        prod_${current_date}_%j.out

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

            # Create a predictions directory for each file processed
            mkdir -p "$dir/predictions"

            # Run Prodigal for the current file
            prodigal -i "$filename" \
                     -o "$dir/predictions/${basename}_genes.fa" \
                     -a "$dir/predictions/${basename}_proteins.faa" \
                     -p single
        fi
    done

    # Recursively process subdirectories
    for subdir in "$dir"/*/; do
        if [ -d "$subdir" ]; then
            process_directory "$subdir"
        fi
    done
}

# Start processing from the current directory
process_directory "$(pwd)"

