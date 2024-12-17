#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 05:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka_v14%j.err
#SBATCH -o 12-prokka_v14%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0

# Function to process directories
process_directory() {
    local dir="$1"
    local current_dir="$2"

    # Loop through .fna files in the current directory
    for filename in "$dir"/*.fna; do
        # Check if the file exists
        if [ -f "$filename" ]; then
            # Extract basename of the file (without extension)
            file_basename=$(basename "$filename" .fna)
            
            # Create directory for this file's Prokka output in the designated output directory
            output_dir="$current_dir/prokka/$file_basename"
            mkdir -p "$output_dir"
            
            # Run Prokka on the fna file
            prokka "$filename" \
            --force \
            --outdir "$output_dir" \
            --prefix "$file_basename" \
            --addgenes \
            --addmrna \
            --coverage 70
        fi
    done
}

# Start processing from the current directory
current_dir="$(pwd)"
process_directory "$current_dir" "$current_dir"

# Process one level below the current directory
for subdir in "$current_dir"/*/; do
    [ -d "$subdir" ] && process_directory "$subdir" "$current_dir"
done

