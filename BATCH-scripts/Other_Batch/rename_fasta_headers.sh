#!/bin/bash

# Directory containing your FASTA files
fasta_dir="/Users/matthewadlam/Desktop/NeSI-2024-04-30/Staph/9_genes"
csv_file="/Users/matthewadlam/Desktop/NeSI-2024-04-30/Staph/summaries/summaries/species_names.csv"

# Python script to read CSV and create a mapping
python_script="
import csv, sys
mapping = {}
with open('$csv_file', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        mapping[row['Accession Number']] = row['Species Name']
for fasta in sys.argv[1:]:
    acc = fasta.split('/')[-1].split('_')[0]
    species = mapping.get(acc, 'Unknown_Species')
    with open(fasta, 'r') as file:
        data = file.read().replace('>Assembly', f'>{species}')
    with open(fasta, 'w') as file:
        file.write(data)
"

# Find all fasta files and apply the Python script
find "$fasta_dir" -name "*.fasta" -exec python -c "$python_script" '{}' +

