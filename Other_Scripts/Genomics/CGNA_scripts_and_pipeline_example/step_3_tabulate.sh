#!/bin/bash

# Set directories
ROOT_DIR="."                       # Root directory (current directory)
SEARCH_DIR="${ROOT_DIR}"           # Directory to search for *_combined.fasta files
TEMPDIR="${ROOT_DIR}/temp_genes"   # Temporary directory
OUTFILE="${ROOT_DIR}/genes_presence_absence_transposed.tsv" # Output file

# Create necessary directories
mkdir -p "$TEMPDIR"

echo "Step 1: Extracting gene names from *_combined.fasta files..."

# Extract gene names from all *_combined.fasta files in subdirectories
find "$SEARCH_DIR" -type f -name "*_combined.fasta" | while read -r fasta_file; do
    base_name=$(basename "$fasta_file" _combined.fasta)
    grep "^>" "$fasta_file" | sed 's/^>//' | awk '{print $1}' > "$TEMPDIR/${base_name}.txt"
    echo "Processed: $fasta_file -> $TEMPDIR/${base_name}.txt"
done

# Step 2: Generate a list of all unique genes
echo "Step 2: Generating unique list of genes..."
awk '{print $1}' $TEMPDIR/*.txt | sort | uniq > "$TEMPDIR/all_genes.txt"

# Step 3: Write the header row for the output file
echo -ne "File" > "$OUTFILE"
awk '{printf "\t%s", $1}' "$TEMPDIR/all_genes.txt" >> "$OUTFILE"
echo "" >> "$OUTFILE"

# Step 4: Generate presence/absence matrix
echo "Step 3: Generating presence/absence matrix..."
find "$SEARCH_DIR" -type f -name "*_combined.fasta" | while read -r fasta_file; do
    base_name=$(basename "$fasta_file" _combined.fasta)
    echo -ne "$base_name" >> "$OUTFILE"
    while read -r gene; do
        if grep -q "^${gene}$" "$TEMPDIR/${base_name}.txt"; then
            echo -ne "\t1" >> "$OUTFILE"
        else
            echo -ne "\t0" >> "$OUTFILE"
        fi
    done < "$TEMPDIR/all_genes.txt"
    echo "" >> "$OUTFILE"
done

# Cleanup temporary directory
rm -rf "$TEMPDIR"

echo "Gene presence/absence matrix generated in: $OUTFILE"
