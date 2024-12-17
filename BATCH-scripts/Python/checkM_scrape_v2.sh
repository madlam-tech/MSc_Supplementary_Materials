#!/bin/bash

# Define the output file
output_file="checkm_summary_v2.tsv"

# Initialize the output file and write the header
echo -e "Filename\tDirectory\tBin Id\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity" > "$output_file"

# Loop through all subdirectories containing 'checkm.txt'
find . -type f -name "checkm.txt" -print0 | while IFS= read -r -d $'\0' file; do
    # Extract the directory name of the checkm.txt file
    dir_name=$(dirname "$file")
    
    # Extract the GCF number from the grandparent directory
    grandparent_dir=$(basename "$(dirname "$dir_name")")

    # Read the second line from checkm.txt (data line)
    data_line=$(tail -n +2 "$file" | head -n 1)

    # Append the directory name and data line to the output file
    echo -e "${grandparent_dir}\t${grandparent_dir}\t${data_line}" >> "$output_file"
done

echo "Processing complete. Output saved to $output_file."

