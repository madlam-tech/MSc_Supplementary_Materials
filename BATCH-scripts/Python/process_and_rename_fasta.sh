#!/bin/bash

# Directory containing your FASTA files and the TXT file
fasta_dir="."

# TXT file in the current directory
txt_file="./species_names.txt"

# Name of your conda environment
conda_env="myenv"

# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$conda_env"

# Python script to read TXT and create a mapping
python_script="
import csv
import sys
import os

# Read TXT file to create a mapping from accession number to species name
def read_species_mapping(txt_file):
    mapping = {}
    with open(txt_file, 'r') as txtfile:
        for line in txtfile:
            if line.strip() and not line.startswith('Accession Number') and not line.startswith('==='):
                acc, species = line.split(None, 1)
                mapping[acc] = species.strip()
    return mapping

# Replace header in FASTA file with species name
def replace_header(fasta_file, header_name):
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
    if lines:
        lines[0] = f'>{header_name}\n'
    with open(fasta_file, 'w') as file:
        file.writelines(lines)

if __name__ == '__main__':
    species_mapping = read_species_mapping('$txt_file')
    fasta_files = sys.argv[1:]
    for fasta_file in fasta_files:
        acc = os.path.basename(fasta_file).split('_concatenated')[0]
        if acc.startswith('ma'):
            header_name = acc.upper()
        else:
            species_name = species_mapping.get(acc, 'Unknown_Species')
            genus, species = species_name.split()[:2]
            header_name = f'{genus[0]}.{species}'
        print(f'Processing file: {fasta_file}')
        print(f'Extracted accession number: {acc}')
        print(f'Mapped species name: {header_name}')
        replace_header(fasta_file, header_name)
"

# Find all fasta files in the current directory and apply the Python script
find "$fasta_dir" -maxdepth 1 -name "*.fasta" -exec python -c "$python_script" '{}' +

# Deactivate the conda environment
conda deactivate

echo "FASTA headers updated successfully."

