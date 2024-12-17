#!/bin/bash

# Directories
INPUT_DIR="./MLST_intermediate_genes/combined_files"
OUTPUT_DIR="./ordered_fasta_files"
PRESENCE_ABSENCE_FILE="genes_presence_absence_transposed.tsv"
BEFORE_TABLE="genes_summary_before.tsv"
AFTER_TABLE="genes_summary_after.tsv"

# Create output directories
mkdir -p "$OUTPUT_DIR"

# Step 1: Extract common genes
echo "Identifying common genes..."
common_genes=$(awk '
NR == 1 { for (i = 2; i <= NF; i++) genes[i] = $i }  # Capture all gene names
NR > 1 { 
  for (i = 2; i <= NF; i++) {
    if ($i == 0) delete genes[i]                   # Remove genes with any absence (0)
  }
}
END {
  for (i in genes) print genes[i]                  # Print common genes
}' "$PRESENCE_ABSENCE_FILE")

# Check if any common genes were found
if [[ -z "$common_genes" ]]; then
    echo "Error: No common genes found across all files. Check your input matrix."
    exit 1
fi

IFS=$'\n' read -r -d '' -a COMMON_GENES <<< "$common_genes"
echo "Common genes: ${COMMON_GENES[*]}"

# Step 2: Generate "BEFORE" table (gene presence/absence)
echo "Generating 'before' table..."
echo -e "File\tGenes_Present" > "$BEFORE_TABLE"
for fasta_file in "$INPUT_DIR"/*_combined.fasta; do
    base_name=$(basename "$fasta_file" _combined.fasta)
    genes=$(grep "^>" "$fasta_file" | sed 's/^>//' | paste -sd "," -)
    echo -e "${base_name}\t${genes}" >> "$BEFORE_TABLE"
done

# Step 3: Reorder FASTA files based on common genes
echo "Reordering FASTA files..."
for fasta_file in "$INPUT_DIR"/*_combined.fasta; do
    base_name=$(basename "$fasta_file" _combined.fasta)
    output_file="$OUTPUT_DIR/${base_name}_ordered.fasta"

    # Initialize output file
    > "$output_file"

    # Extract each common gene in order
    for gene in "${COMMON_GENES[@]}"; do
        awk -v gene=">$gene" '
        $0 ~ gene {f=1; print; next} 
        /^>/ {f=0} 
        f {print}
        ' "$fasta_file" >> "$output_file"
    done
    echo "Generated ordered FASTA file: $output_file"
done

# Step 4: Generate "AFTER" table (gene order in new FASTA files)
echo "Generating 'after' table..."
echo -e "File\tGene_Order" > "$AFTER_TABLE"
for fasta_file in "$OUTPUT_DIR"/*_ordered.fasta; do
    base_name=$(basename "$fasta_file" _ordered.fasta)
    gene_order=$(grep "^>" "$fasta_file" | sed 's/^>//' | paste -sd "," -)
    echo -e "${base_name}\t${gene_order}" >> "$AFTER_TABLE"
done

echo "Processing complete!"
echo "Tables generated:"
echo " - Before editing table: $BEFORE_TABLE"
echo " - After editing table: $AFTER_TABLE"
echo "Ordered FASTA files are in: $OUTPUT_DIR"
