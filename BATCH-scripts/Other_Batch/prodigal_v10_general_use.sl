#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v10
#SBATCH --time          01:00:00
#SBATCH --mem          	1GBB
#SBATCH --cpus-per-task 1
#SBATCH --error         prod%j.err
#SBATCH --output        prod%j.out

# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

# Function to process directories
process_directory() {
    local dir="$1"

    # Loop through .fna files in the current directory
    for filename in "$dir"/*.fna; do
        # Check if the file exists
        if [ -f "$filename" ]; then
            # Extract basename without extension
            basename=$(basename "$filename" .fna)

            # Create the output directory in the same directory where the .fna file is found
            output_dir="$dir/predictions"
            mkdir -p "$output_dir"

            # Run Prodigal for the current file
            prodigal -i "$filename" \
                     -o "$output_dir/${basename}_genes.fa" \
                     -a "$output_dir/${basename}_proteins.faa" \
                     -p single
        fi
    done
}

# Start processing from the current directory
current_dir="$(pwd)"
process_directory "$current_dir"

# Process one level below the current directory
for subdir in "$current_dir"/*/; do
    [ -d "$subdir" ] && process_directory "$subdir"
done

