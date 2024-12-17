#!/bin/bash

# Directories
barrnap_dir="./barrnap_output_20241217"    # Directory with *_rrna.fa files
genes_dir="./MLST_intermediate_genes"      # Directory with *_11_genes.gbk files
output_dir="./MLST_intermediate_genes/combined_files" # Combined files output directory

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Function to append cleaned first 16S gene
append_cleaned_first_16S_gene() {
    local rrna_file="$1"
    local genes_file="$2"
    local output_file="$3"

    echo "Processing: 16S file = $rrna_file, genes file = $genes_file"

    if [[ -f "$rrna_file" && -f "$genes_file" ]]; then
        # Extract the first 16S sequence
        awk '/^>.*16S/ {if (found) exit; found=1} found {print}' "$rrna_file" | \
            sed '1s/.*/>16S/' > temp_16S.fa
        
        if [[ -s temp_16S.fa ]]; then
            # Combine files
            cat "$genes_file" temp_16S.fa > "$output_file"
            echo "Saved combined file: $output_file"
            rm temp_16S.fa
        else
            echo "ERROR: No 16S rRNA sequence found in $rrna_file"
            rm temp_16S.fa
        fi
    else
        echo "ERROR: One or both files do not exist: $rrna_file, $genes_file"
    fi
}

# Step 1: Find and combine files
echo "Step 1: Combining 16S rRNA and gene files..."
find "$barrnap_dir" -type f -name "*_rrna.fa" | while read -r rrna_file; do
    # Extract base name
    base_name=$(basename "$rrna_file" _rrna.fa)
    
    # Corresponding genes file
    genes_file="${genes_dir}/${base_name}_11_genes.gbk"
    output_file="${output_dir}/${base_name}_combined.fasta"

    # Combine files
    append_cleaned_first_16S_gene "$rrna_file" "$genes_file" "$output_file"
done

# Step 2: Specific handling for ma1263 files
echo "Step 2: Handling specific ma1263 files..."
find "$barrnap_dir" -type f -name "ma1263*rrna.fa" | while read -r ma_file; do
    base_name=$(basename "$ma_file" _rrna.fa)
    genes_file="${genes_dir}/${base_name}_11_genes.gbk"
    output_file="${output_dir}/${base_name}_combined.fasta"

    append_cleaned_first_16S_gene "$ma_file" "$genes_file" "$output_file"
done

echo "File processing completed. Combined files are in $output_dir"
