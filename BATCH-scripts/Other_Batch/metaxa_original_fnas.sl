#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      metaxa2
#SBATCH --time          05:00:00
#SBATCH --mem           10GB
#SBATCH --cpus-per-task 4
#SBATCH --error         metaxa_%j.err
#SBATCH --output        metaxa_%j.out

# Load modules
module purge
module load Metaxa2/2.2.3-gimkl-2022a

# Define input directory
input_dir="/nesi/nobackup/massey03345/fnas_4_annotation"

# Define output directory
output_dir="/nesi/nobackup/massey03345/metaxa_from_original_fnas"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .fna file in the input directory
for fna_file in "$input_dir"/*.fna; do
    # Extract the filename without extension
    filename=$(basename "$fna_file" .fna)
    # Create a new directory for Metaxa2 output
    mkdir -p "$output_dir/${filename}_metaxa"
    
    # Run Metaxa2 for each ribosome type
    for ribosome_type in ssu lsu; do
        metaxa2 -g "$ribosome_type" --mode genome \
                -i "$fna_file" \
                -o "$output_dir/${filename}_metaxa/${filename}_genes.fa.${ribosome_type}"
    done
done

