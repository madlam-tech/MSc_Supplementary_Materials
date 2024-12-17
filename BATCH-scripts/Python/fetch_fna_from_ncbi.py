import os
import sys
from Bio import Entrez
from Bio import SeqIO

# Set the email address to be used with Entrez
Entrez.email = "your_email@example.com"  # Replace with your email

# Define the input TSV file and output directory
species_file = "species_names.txt"
output_dir = "./fnas_complete_set"

# Create the output directory if it does not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def fetch_fna_file(accession):
    try:
        # Search for the nucleotide entry
        handle = Entrez.esearch(db="nucleotide", term=accession)
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            nucleotide_id = record["IdList"][0]
            
            # Fetch the nucleotide sequence
            fetch_handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
            sequence_data = fetch_handle.read()
            fetch_handle.close()
            
            # Save the sequence data to a file
            output_file = os.path.join(output_dir, f"{accession}.fna")
            with open(output_file, "w") as f:
                f.write(sequence_data)
            print(f"Successfully fetched and saved {accession}.fna")
        else:
            print(f"No results found for {accession}")
    except Exception as e:
        print(f"Error fetching data for {accession}: {e}")

# Read the species file and extract accession numbers
with open(species_file, "r") as f:
    lines = f.readlines()[3:]  # Skip the first three lines (header)
    for line in lines:
        accession = line.split()[0]
        fetch_fna_file(accession)

print("Processing complete. All .fna files have been fetched and saved to", output_dir)
