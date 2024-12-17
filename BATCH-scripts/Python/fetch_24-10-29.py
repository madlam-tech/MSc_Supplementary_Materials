import os
import sys
import pandas as pd
from Bio import Entrez

# Set your email for NCBI Entrez
Entrez.email = "your_email@example.com"  # Replace with your email

def fetch_species(accession):
    try:
        handle = Entrez.esummary(db="assembly", term=accession)
        record = Entrez.read(handle)
        handle.close()
        
        if record["DocumentSummarySet"]["DocumentSummary"]:
            docsum = record["DocumentSummarySet"]["DocumentSummary"][0]
            return docsum.get("Organism", "Unknown Species")
        else:
            return "Species Not Found"
    except Exception as e:
        print(f"Error fetching species for {accession}: {e}")
        return "Error"

def main(input_directory, output_file):
    results = []
    for filename in os.listdir(input_directory):
        if filename.endswith(".fna"):
            # Extract accession number from filename
            accession = filename.split("_")[1]
            species = fetch_species(accession)
            results.append({"Accession": accession, "Species": species})
            print(f"Fetched {species} for {accession}")

    # Save results to a TSV file
    df = pd.DataFrame(results)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fetch_species_name.py <input_directory> <output_file>")
    else:
        input_directory = sys.argv[1]
        output_file = sys.argv[2]
        main(input_directory, output_file)
