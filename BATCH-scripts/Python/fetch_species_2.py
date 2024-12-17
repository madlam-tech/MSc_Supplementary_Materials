#!/usr/bin/env python3

from Bio import Entrez
import os

# Set the email address to be used with Entrez
Entrez.email = "madlam@me.com"

# Get the list of accession numbers from the files in the ./stripped_concatenated directory
input_dir = "./stripped_concatenated"
accession_numbers = [filename.split('_concatenated')[0] for filename in os.listdir(input_dir) if filename.endswith('.fasta')]

def fetch_species_name(accession):
    handle = Entrez.esearch(db="assembly", term=accession)
    record = Entrez.read(handle)
    handle.close()
    if record["IdList"]:
        assembly_id = record["IdList"][0]
        summary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary_record = Entrez.read(summary_handle)
        summary_handle.close()
        try:
            return summary_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        except (IndexError, KeyError):
            return "Unknown species"
    return "Unknown species"

# Fetch species names for each accession number
species_names = {acc: fetch_species_name(acc) for acc in accession_numbers}

# Write the results to a TSV file
with open("species_names.tsv", "w") as tsv_file:
    tsv_file.write(f"{'Accession Number'}\t{'Species'}\n")
    for acc, species in species_names.items():
        tsv_file.write(f"{acc}\t{species}\n")
