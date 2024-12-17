#!/bin/bash

# Loop through all .fasta files in the current directory
for file in *.fasta; do
    # Extract the basename up to the second underscore
    # This assumes that there are at least two underscores in the filename
    base_name=$(echo $(basename "$file" .fasta) | cut -d'_' -f1-2)

    # New file name
    new_file="${base_name}_concatenated.fasta"

    # Process each file
    # Remove headers, concatenate sequences, and prepend new header
    echo ">${base_name}" > "$new_file"
    awk 'BEGIN {RS=">"; ORS=""} 
        NR > 1 {
            header = $1; 
            sub(header, ""); 
            gsub("\n", ""); 
            print
        }' "$file" >> "$new_file"

    echo "" >> "$new_file" # Ensure the file ends with a newline
done

echo "Processing complete."


