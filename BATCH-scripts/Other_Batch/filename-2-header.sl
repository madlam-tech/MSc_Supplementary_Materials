#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH --job-name      filename-2-header
#SBATCH --time          00:05:00
#SBATCH --mem           1G
#SBATCH -e filename-2-header_%j

mkdir -p filename-2-header

# Input directory containing FASTA files
input_directory="./sequences"

# Output directory for modified FASTA files
output_directory="./filename-2-header"

# Process each FASTA file in the input directory
for input_file in "$input_directory"/*.fasta; do
    # Get the file name without extension
    filename=$(basename "$input_file")
    header_name="${filename%.*}"

    # Construct output file path
    output_file="$output_directory/$filename"

    # Add the file name to the header in the output file
    awk -v header="$header_name" '/^>/{print ">" header; next} {print}' "$input_file" > "$output_file"
done
