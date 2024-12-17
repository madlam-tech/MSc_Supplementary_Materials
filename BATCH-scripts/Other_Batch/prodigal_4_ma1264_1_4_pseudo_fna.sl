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

    # Only process the specific directory
    if [ "$dir" = "/nesi/nobackup/massey03345/fnas_4_annotation" ]; then
        # Specify the input file
        input_file="$dir/ma1263_1_4_pseudo.fna"

        # Check if the file exists
        if [ -f "$input_file" ]; then
            # Extract basename without extension
            basename=$(basename "$input_file" .fna)

            # Create a directory for the output
            mkdir -p "$dir/prodigal"

            # Run Prodigal for the current file
            prodigal -i "$input_file" \
                     -o "$dir/prodigal/${basename}_genes.fa" \
                     -a "$dir/prodigal/${basename}_proteins.faa" \
                     -p single
        fi
    fi

    # No need to recursively process subdirectories
}

# Start processing from the specified directory
process_directory "/nesi/nobackup/massey03345/fnas_4_annotation"

