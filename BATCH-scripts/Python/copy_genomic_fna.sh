#!/bin/bash

# Define the destination directory
destination_dir="./fnas_complete_set"

# Create the destination directory if it does not exist
mkdir -p "$destination_dir"

# Find and copy all *genomic.fna files in all subdirectories
find ./ -type f -name '*genomic.fna' -exec cp {} "$destination_dir" \;

echo "All *genomic.fna files have been copied to $destination_dir."
