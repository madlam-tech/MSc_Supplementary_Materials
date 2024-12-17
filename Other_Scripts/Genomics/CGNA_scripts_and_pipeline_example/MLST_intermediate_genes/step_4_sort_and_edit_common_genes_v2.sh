#!/bin/bash

# Define the directory containing combined FASTA files and the output directory
INPUT_DIR="."
OUTPUT_DIR="./ordered_fasta_files"
PRESENCE_ABSENCE_FILE="genes_presence_absence_transposed.tsv"
mkdir -p "$OUTPUT_DIR"

# Extract common genes from the presence-absence matrix
common_genes=$(awk '
NR == 1 { for (i = 2; i <= NF; i++) headers[i] = $i }
NR > 1 {
  for (i = 2; i <= NF; i++) {
    if ($i == 0) delete headers[i]
  }
}
END {
  for (i in headers) print headers[i]
}' "$PRESENCE_ABSENCE_FILE")

# Convert the common genes into an array
IFS=$'\n' read -r -d '' -a COMMON_GENES <<< "$common_genes"

# Print common genes for debugging
echo "Common genes: ${COMMON_GENES[*]}"

# Function to extract a gene sequence from a file
extract_gene_sequence() {
    local fasta_file="$1"
    local gene="$2"
    local gene_sequence=""

    gene_sequence=$(awk -v gene=">$gene" '$0 ~ gene {f=1; print; next} /^>/ {f=0} f {print}' "$fasta_file")
    echo "$gene_sequence"
}

# Process each combined FASTA file
for fasta_file in "$INPUT_DIR"/*_combined.fasta; do
    base_name=$(basename "$fasta_file" _combined.fasta)
    output_file="$OUTPUT_DIR/${base_name}_ordered.fasta"

    # Initialize the output file
    > "$output_file"

    # Extract and append each common gene in the defined order
    for gene in "${COMMON_GENES[@]}"; do
        gene_sequence=$(extract_gene_sequence "$fasta_file" "$gene")
        
        # Debugging output
        echo "Processing gene $gene in file $fasta_file"
        if [[ -n "$gene_sequence" ]]; then
            echo "Found sequence for gene $gene in file $fasta_file"
            echo "$gene_sequence" >> "$output_file"
        else
            echo "Gene $gene not found in file $fasta_file"
        fi
    done

    echo "Generated ordered FASTA file: $output_file"
done

echo "All files processed. Ordered FASTA files are in: $OUTPUT_DIR"
