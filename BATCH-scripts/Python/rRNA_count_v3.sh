#!/bin/bash

# Output file
output_file="barrnap_rRNA_counts.tsv"
echo -e "Directory\tFile\t16S\t23S\t5S" > "$output_file"

# Loop through directories starting with 'barrnap_output_'
for dir in barrnap_output_*; do
    # Find all *_rrna.fa files up to 2 directories deep
    find "$dir" -maxdepth 2 -type f -name "*_rrna.fa" | while read -r file; do
        # Initialize counters
        count_16S=0
        count_23S=0
        count_5S=0

        # Count occurrences of each rRNA type
        count_16S=$(grep -c "16S_rRNA" "$file")
        count_23S=$(grep -c "23S_rRNA" "$file")
        count_5S=$(grep -c "5S_rRNA" "$file")

        # Extract the base name of the current directory and file
        base_name=$(basename "$dir")
        file_name=$(basename "$file")

        # Write to the output file
        echo -e "$base_name\t$file_name\t$count_16S\t$count_23S\t$count_5S" >> "$output_file"
    done
done

# Display the result (optional, remove this line if not needed)
cat "$output_file"
