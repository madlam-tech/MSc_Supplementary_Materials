import os
from Bio import Entrez
from Bio import SeqIO

# Set your email for NCBI Entrez
Entrez.email = "your_email@example.com"

# Input and output directories
INPUT_DIR = "./concatenated_fasta_files"
OUTPUT_DIR = "./reheaded_fasta_files"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to fetch genus and species from NCBI Assembly database
def fetch_genus_species(accession):
    try:
        print(f"Fetching genus and species for accession: {accession}")
        
        # Query NCBI Assembly database
        handle = Entrez.esearch(db="assembly", term=accession, retmode="xml")
        search_result = Entrez.read(handle)
        handle.close()

        if not search_result["IdList"]:
            print(f"Warning: No result found for {accession}")
            return None

        # Fetch assembly summary
        assembly_id = search_result["IdList"][0]
        summary_handle = Entrez.esummary(db="assembly", id=assembly_id, retmode="xml")
        summary_result = Entrez.read(summary_handle)
        summary_handle.close()

        # Extract organism name
        organism = summary_result['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        genus, species = organism.split()[:2]
        formatted_name = f"{genus[0]}.{species}"
        print(f"Fetched: {accession} -> {formatted_name}")
        return formatted_name
    except Exception as e:
        print(f"Error fetching data for {accession}: {e}")
        return None

# Reheading the FASTA files
def rehead_fasta_files():
    for fasta_file in os.listdir(INPUT_DIR):
        if fasta_file.endswith("_concatenated.fasta"):
            accession = fasta_file.split("_")[0] + "_" + fasta_file.split("_")[1]
            input_path = os.path.join(INPUT_DIR, fasta_file)
            output_path = os.path.join(OUTPUT_DIR, fasta_file)

            # Fetch genus and species
            new_header = fetch_genus_species(accession)
            if not new_header:
                new_header = accession  # Fallback to accession name

            # Rehead the FASTA file
            with open(input_path, "r") as infile, open(output_path, "w") as outfile:
                for record in SeqIO.parse(infile, "fasta"):
                    record.id = new_header
                    record.description = ""
                    SeqIO.write(record, outfile, "fasta")
            print(f"Processed: {fasta_file} -> {output_path}")

if __name__ == "__main__":
    print("Fetching genus and species from NCBI and reheading FASTA files...")
    rehead_fasta_files()
    print("Processing complete. Reheaded files are in:", OUTPUT_DIR)
