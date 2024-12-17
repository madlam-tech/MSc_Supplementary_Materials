import os
import sys
import time
import pandas as pd
from Bio import Entrez

# Set the email address to be used with Entrez
Entrez.email = "madlam@me.com"  # Replace with your email

def fetch_species_data(accession, retries=3, delay=5):
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="assembly", term=accession, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                assembly_id = record["IdList"][0]
                summary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
                summary_record = Entrez.read(summary_handle)
                summary_handle.close()
                docsum = summary_record['DocumentSummarySet']['DocumentSummary'][0]
                data = {
                    'AssemblyAccession': docsum.get('AssemblyAccession', 'N/A'),
                    'Bioproject': docsum.get('Bioproject', 'N/A'),
                    'Biosample': docsum.get('Biosample', 'N/A'),
                    'WGS': docsum.get('WGS', 'N/A'),
                    'RefSeq_category': docsum.get('RefSeq_category', 'N/A'),
                    'Taxid': docsum.get('Taxid', 'N/A'),
                    'SpeciesTaxid': docsum.get('SpeciesTaxid', 'N/A'),
                    'Organism': docsum.get('Organism', 'N/A'),
                    'Infraspecific_name': docsum.get('Infraspecific_name', 'N/A'),
                    'Isolate': docsum.get('Isolate', 'N/A'),
                    'VersionStatus': docsum.get('VersionStatus', 'N/A'),
                    'AssemblyLevel': docsum.get('AssemblyLevel', 'N/A'),
                    'ReleaseType': docsum.get('ReleaseType', 'N/A'),
                    'GenbankAssemblyAccession': docsum.get('GenbankAssemblyAccession', 'N/A'),
                    'SubGroup': docsum.get('SubGroup', 'N/A')
                }
                return data
        except Exception as e:
            print(f"Error fetching data for {accession}, attempt {attempt+1}/{retries}: {e}")
            time.sleep(delay)
    return None

def main(input_directory, output_tsv):
    data = []
    for filename in os.listdir(input_directory):
        if filename.endswith(".fasta"):
            basename = "_".join(filename.split('_')[:2])
            fetched_data = fetch_species_data(basename)
            if fetched_data:
                data.append(fetched_data)
    
    columns = [
        'AssemblyAccession', 'Bioproject', 'Biosample', 'WGS', 'RefSeq_category',
        'Taxid', 'SpeciesTaxid', 'Organism', 'Infraspecific_name', 'Isolate',
        'VersionStatus', 'AssemblyLevel', 'ReleaseType', 'GenbankAssemblyAccession',
        'SubGroup'
    ]
    
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python fetch_ref_data.py <input_directory> <output_tsv>")
    else:
        input_directory = sys.argv[1]
        output_tsv = sys.argv[2]
        main(input_directory, output_tsv)

