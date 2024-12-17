#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          06:00:00
#SBATCH --mem           10GB
#SBATCH --cpus-per-task 4
#SBATCH --error         prod%j.err
#SBATCH --output        prod%j.out

# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

# Define the file to process
fna_file="GCF_000006765.1_ASM676v1_genomic.fna"

# Define the output directory
output_dir="prokka"

# Check if the file exists
if [ -f "$fna_file" ]; then
    # Extract basename without extension
    file_basename=$(basename "$fna_file" .fna)

    # Create base directory for this file's output
    base_dir="$output_dir/$file_basename"
    mkdir -p "$base_dir/predictions"
    mkdir -p "$base_dir/checkM"

    # Run Prodigal on the .fna file
    prodigal -i "$fna_file" \
             -o "$base_dir/predictions/${file_basename}_genes.gbk" \
             -a "$base_dir/predictions/${file_basename}_proteins.faa" \
             -d "$base_dir/predictions/${file_basename}_nucl.fnn" \
             -p meta

    # Additional processing for checkM can be added here if needed
    # Placeholder for checkM processing
    # checkm command example: checkm analyze -o 2 -t 4 -x fna "$fna_file" "$base_dir/checkM"
else
    echo "File $fna_file does not exist."
    exit 1
fi

# Feedback
echo "Processing complete. Output saved to $base_dir."
