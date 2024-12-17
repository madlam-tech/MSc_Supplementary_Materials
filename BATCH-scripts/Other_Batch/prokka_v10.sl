#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka_v5
#SBATCH --time 03:30:00
#SBATCH --mem 40GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka_v10%j.err
#SBATCH -o 12-prokka_v10%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0
module load RNAmmer/1.2-GCC-9.2.0-Perl-5.30.1

# Date stamp for output directory
date_stamp=$(date +%Y%m%d)
output_dir="prokka_output_$date_stamp"
mkdir -p "$output_dir"

# Function to process directories recursively
process_directory() {
    local dir="$1"

    # Loop through .fna files in the current directory
    for fna_file in "$dir"/*.fna; do
        # Check if the file exists
        if [ -f "$fna_file" ]; then
            # Extract basename without extension
            file_basename=$(basename "$fna_file" .fna)

            # Create directory for this file's Prokka output
            mkdir -p "$output_dir/$file_basename"

            # Run Prokka with RNAmmer on the .fna file
            prokka "$fna_file" \
            --force \
            --outdir "$output_dir/$file_basename" \
            --prefix "$file_basename" \
            --addgenes \
            --addmrna \
            --coverage 70 \
            --rnammer
        fi
    done

    # Recursively process subdirectories
    for subdir in "$dir"/*/; do
        process_directory "$subdir"
    done
}

# Start processing from the main directory
main_dir="/nesi/nobackup/massey03345/fnas_4_annotation"
process_directory "$main_dir"

