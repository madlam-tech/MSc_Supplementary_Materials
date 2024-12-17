import os
import sys
from Bio import SeqIO
from Bio import Entrez
import pandas as pd

# Set the email address to be used with Entrez
Entrez.email = "madlam@me.com"  # Replace with your email

def fetch_species_name_and_data(accession):
    try:
        handle = Entrez.esearch(db="assembly", term=accession)
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            assembly_id = record["IdList"][0]
            summary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            species_name = summary_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
            # Split the species name to get genus and species separately
            genus, species = species_name.split()[:2]
            return accession, genus, species
    except Exception as e:
        print(f"Error fetching species name for {accession}: {e}")
    return accession, "Unknown", "Unknown"

def update_headers_and_create_table(input_dir, output_file):
    data = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            accession = "_".join(filename.split('_')[:2])
            accession, genus, species = fetch_species_name_and_data(accession)
            data.append([accession, genus, species])
    
    df = pd.DataFrame(data, columns=["Accession", "Genus", "Species"])
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    input_directory = sys.argv[1]
    output_file = sys.argv[2]
    update_headers_and_create_table(input_directory, output_file)
