#!/bin/bash

# Directory containing the fasta files
DIR="./"

# Loop through each fasta file in the directory
for file in "$DIR"/*.fasta; do
    # Use sed to remove the '_1' suffix and save to a temporary file
    sed 's/_1$//' "$file" > "${file%.fasta}_no_suffix.fasta"
done

