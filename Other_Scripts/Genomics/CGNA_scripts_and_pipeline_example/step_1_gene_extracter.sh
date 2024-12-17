#!/bin/bash

# Define the input directory (prokka_7 and its subdirectories)
input_dir="./prokka_7"

# Define the output directory
output_dir="./MLST_intermediate_genes"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Define gene names to be extracted, ensuring correct names and separators
genes="atpD|lepA|gltB|16S|gyrB|phaC|rpoB|recA|infB|trpB"

# Find all .gbk files in the specified directory and subdirectories
find "$input_dir" -type f -name "*.gbk" | while read -r gbk_file; do
    # Extract the base name for output naming
    base_name=$(basename "${gbk_file}" .gbk)

    # Define the output file path
    output_file="${output_dir}/${base_name}_11_genes.gbk"

    # Call Python script to extract the genes
    # Assuming 'gene_extractor.py' is a script that accepts these parameters and handles gene extraction
    python3 gene_extractor.py "$gbk_file" "$output_file" "$genes"

    echo "Processed $gbk_file and output to $output_file"
done

echo "All files processed. Output stored in $output_dir"
