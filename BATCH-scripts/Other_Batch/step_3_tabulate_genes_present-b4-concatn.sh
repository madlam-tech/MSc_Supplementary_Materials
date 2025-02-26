#!/bin/bash

# Directory containing FASTA files
DIR="."

# Temporary and output directory setup
TEMPDIR="./temp_genes"
OUTFILE="./genes_presence_absence_transposed.tsv"
mkdir -p "$TEMPDIR"

# Extract gene names from each file and store in temporary files
for fasta_file in "$DIR"/*.fasta; do
    base_name=$(basename "$fasta_file" .fasta)
    grep "^>" "$fasta_file" | sed 's/>//' > "$TEMPDIR/${base_name}.txt"
done

# Create a list of all unique genes
awk '{print $1}' $TEMPDIR/*.txt | sort | uniq > "$TEMPDIR/all_genes.txt"

# Print header for the output file
echo -ne "File\t" > "$OUTFILE"
while read -r gene; do
    echo -ne "$gene\t" >> "$OUTFILE"
done < "$TEMPDIR/all_genes.txt"
echo "" >> "$OUTFILE"

# Initialize an array to hold the sums for each gene
declare -A gene_sums
while read -r gene; do
    gene_sums["$gene"]=0
done < "$TEMPDIR/all_genes.txt"

# Fill the table with 1s and 0s indicating presence/absence of genes
for fasta_file in "$DIR"/*.fasta; do
    base_name=$(basename "$fasta_file" .fasta)
    echo -ne "$base_name\t" >> "$OUTFILE"
    while read -r gene; do
        if grep -q "^$gene$" "$TEMPDIR/${base_name}.txt"; then
            echo -ne "1\t" >> "$OUTFILE"
            ((gene_sums["$gene"]++))
        else
            echo -ne "0\t" >> "$OUTFILE"
        fi
    done < "$TEMPDIR/all_genes.txt"
    echo "" >> "$OUTFILE"
done

# Add the bottom row with the sum of times each gene is present
echo -ne "Sum\t" >> "$OUTFILE"
while read -r gene; do
    echo -ne "${gene_sums["$gene"]}\t" >> "$OUTFILE"
done < "$TEMPDIR/all_genes.txt"
echo "" >> "$OUTFILE"

# Cleanup temporary directory
rm -rf "$TEMPDIR"

echo "Transposed gene presence/absence matrix with sums generated in: $OUTFILE"

