#!/bin/bash

# Function to extract the first 16S gene from the rrna file, clean the header, and append it to the genes file
append_cleaned_first_16S_gene() {
    local rrna_file="$1"
    local genes_file="$2"
    local output_file="$3"

    echo "Checking files: $rrna_file and $genes_file"

    # Check if both files exist
    if [[ -f "$rrna_file" && -f "$genes_file" ]]; then
        # Extract the first 16S gene sequence from the rrna file and clean the header
        awk '/^>.*16S/ {if (found) exit; found=1} found {print}' "$rrna_file" | sed '1s/.*/>16S/' > first_16S_gene.fa

        # Append the cleaned first 16S gene to the genes file and save to the output file
        cat "$genes_file" first_16S_gene.fa > "$output_file"
        echo "Appended and cleaned first 16S gene from $rrna_file to $genes_file and saved to $output_file"

        # Clean up temporary file
        rm first_16S_gene.fa
    else
        echo "Files $rrna_file and/or $genes_file do not exist."
    fi
}

# Directory to start searching (current directory)
start_directory="."

# Find all *_16S_rRNA.fa files in all subdirectories and append the cleaned first 16S gene to their corresponding *_11_genes.gbk files
find "$start_directory" -type f -name "*_16S_rRNA.fa" | while read -r rrna_file; do
    # Extract the base name without the _16S_rRNA.fa suffix
    base_name=$(basename "$rrna_file" _16S_rRNA.fa)

    echo "Processing $base_name"

    # Find the corresponding genes file in all subdirectories
    genes_file=$(find "$start_directory" -type f -name "${base_name}_11_genes.gbk" | head -n 1)

    if [[ -n "$genes_file" ]]; then
        # Define the output file
        output_file="$(dirname "$rrna_file")/${base_name}_combined.fasta"

        # Append the cleaned first 16S gene to the genes file
        append_cleaned_first_16S_gene "$rrna_file" "$genes_file" "$output_file"
    else
        echo "Corresponding genes file ${base_name}_11_genes.gbk not found for $rrna_file"
    fi
done

# Specific handling for ma1263 files
for ma_file in $(find "$start_directory" -type f -name "ma1263*16S.rrna.fa"); do
    base_name=$(basename "$ma_file" _16S.rrna.fa)
    echo "Processing specific ma1263 file: $ma_file"

    genes_file=$(find "$start_directory" -type f -name "${base_name}_11_genes.gbk" | head -n 1)
    if [[ -n "$genes_file" ]]; then
        output_file="$(dirname "$ma_file")/${base_name}_combined.fasta"
        append_cleaned_first_16S_gene "$ma_file" "$genes_file" "$output_file"
    else
        echo "Corresponding genes file ${base_name}_11_genes.gbk not found for $ma_file"
    fi
done

echo "File processing completed."

