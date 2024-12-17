#!/bin/bash

# Define the input and output directories
INPUT_DIR="./ordered_fasta_files"
OUTPUT_DIR="./concatenated_fasta_files"
mkdir -p "$OUTPUT_DIR"

# Loop through all .fasta files in the ordered directory
for file in "$INPUT_DIR"/*.fasta; do
    # Ensure the file exists
    [ -e "$file" ] || continue

    # Extract the base name up to the second underscore
    base_name=$(basename "$file" .fasta | awk -F'_' '{print $1"_"$2}')

    # Define the new concatenated output file
    new_file="${OUTPUT_DIR}/${base_name}_concatenated.fasta"

    # Prepend the header and concatenate sequences
    echo ">${base_name}" > "$new_file"

    # Extract all sequences and concatenate them into a single line
    awk '
    /^>/ {next}                # Skip headers
    {gsub(/\n/, ""); print}    # Remove newlines and print sequences
    ' "$file" | tr -d '\n' >> "$new_file"  # Concatenate without any breaks

    echo "" >> "$new_file"  # Ensure a final newline for proper formatting

    echo "Processed: $file -> $new_file"
done

echo "All ordered FASTA files have been concatenated. Results are in: $OUTPUT_DIR"
