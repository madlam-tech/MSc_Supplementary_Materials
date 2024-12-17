import os
import sys
import pandas as pd
from Bio import Entrez

# Set your email here
Entrez.email = "your_email@example.com"

def fetch_species_data(accession):
    try:
        handle = Entrez.esearch(db="assembly", term=accession)
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            assembly_id = record["IdList"][0]
            summary_handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            doc_summary = summary_record['DocumentSummarySet']['DocumentSummary'][0]
            data = {
                "Accession": accession,
                "Title": doc_summary.get("Title", "N/A"),
                "Organism": doc_summary.get("Organism", "N/A"),
                "TaxId": doc_summary.get("Taxid", "N/A"),
                "Bioproject": doc_summary.get("BioprojectId", "N/A"),
                "Biosample": doc_summary.get("BiosampleId", "N/A"),
                "AssemblyName": doc_summary.get("AssemblyName", "N/A"),
                "SubmissionDate": doc_summary.get("SubmissionDate", "N/A")
            }
            return data
        else:
            print(f"No records found for accession: {accession}")
    except Exception as e:
        print(f"Error fetching data for {accession}: {e}")
    return None

def main(input_file, output_tsv):
    data = []
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            filename = os.path.basename(line)
            if filename.endswith(".fna"):
                parts = filename.split('_')
                if len(parts) >= 2 and parts[0].startswith("GCF"):
                    accession = parts[0] + '_' + parts[1]
                    print(f"Processing accession: {accession}")  # Debugging statement
                    species_data = fetch_species_data(accession)
                    if species_data:
                        data.append(species_data)
                else:
                    print(f"Skipping invalid filename: {filename}")
    
    output_df = pd.DataFrame(data)
    output_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Saved data to {output_tsv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fetch_ncbi_data.py <input_file> <output_tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_tsv = sys.argv[2]
    main(input_file, output_tsv)
