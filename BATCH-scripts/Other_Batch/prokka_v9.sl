#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 04:30:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka%j.err
#SBATCH -o 12-prokka%j.log

# Load modules
module purge
module load prokka/1.14.5-GCC-9.2.0

# Function to process each .fna file found
process_file() {
    local file="$1"
    local dir=$(dirname "$file")
    local basename=$(basename "$file" .fna)

    # Create a directory for Prokka output
    mkdir -p "$dir/prokka"

    # Run Prokka on the file
    prokka "$file" \
        --outdir "$dir/prokka" \
        --force \
        --prefix "$basename" \
        --addgenes \
        --addmrna \
        --coverage 70
}

# Export the function to be used by find's exec
export -f process_file

# Find and process each .fna file
find . -type f -name "*.fna" -exec bash -c 'process_file "$0"' {} \;

echo "Prokka annotation completed."

