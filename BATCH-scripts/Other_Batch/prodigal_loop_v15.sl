#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          04:00:00
#SBATCH --mem           10GB
#SBATCH --cpus-per-task 4
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

            # Create base directory for this file's output
            base_dir="$dir/$basename"
            mkdir -p "$base_dir/predictions"
            mkdir -p "$base_dir/checkM"

            # Run Prodigal on the .fna file
            prodigal -i "$filename" \
                     -o "$base_dir/predictions/${basename}_genes.gbk" \
                     -a "$base_dir/predictions/${basename}_proteins.faa" \
                     -d "$base_dir/predictions/${basename}_nucl.fnn" \
                     -p meta

            # Additional processing for checkM can be added here if needed
            # Placeholder for checkM processing
            # checkm command example: checkm analyze -o 2 -t 4 -x fna "$filename" "$base_dir/checkM"
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
