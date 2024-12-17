#!/bin/bash

# Define the directory containing combined FASTA files
INPUT_DIR="./MLST_intermediate_files"
OUTPUT_DIR="./ordered_fasta_files"
mkdir -p "$OUTPUT_DIR"

# Define the order of genes
GENES=("16S" "atpD" "gltB" "gyrB" "infB" "lepA" "recA" "rpoB" "trpB")

# Placeholder sequence for missing genes
PLACEHOLDER_SEQUENCE="NNNNNNNNNN"

# Function to extract a gene sequence from a file
extract_gene_sequence() {
    local fasta_file="$1"
    local gene="$2"
    local gene_sequence=""

    gene_sequence=$(awk -v gene="$gene" '/^>/{f=0} /^>/{if($0 ~ gene){f=1}} f' "$fasta_file")
    echo "$gene_sequence"
}

# Process each combined FASTA file
for fasta_file in "$INPUT_DIR"/*_combined.fasta; do
    base_name=$(basename "$fasta_file" _combined.fasta)
    output_file="$OUTPUT_DIR/${base_name}_ordered.fasta"

    # Initialize the output file
    > "$output_file"

    # Extract and append each gene in the defined order
    for gene in "${GENES[@]}"; do
        gene_sequence=$(extract_gene_sequence "$fasta_file" "$gene")
        if [[ -n "$gene_sequence" ]]; then
            echo "$gene_sequence" >> "$output_file"
        else
            echo ">${gene}" >> "$output_file"
            echo "$PLACEHOLDER_SEQUENCE" >> "$output_file"
            echo "Gene $gene not found in $fasta_file, adding placeholder."
        fi
    done

    echo "Generated ordered FASTA file: $output_file"
done

echo "All files processed. Ordered FASTA files are in: $OUTPUT_DIR"

