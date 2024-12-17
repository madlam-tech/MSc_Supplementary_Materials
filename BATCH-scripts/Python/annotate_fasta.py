from Bio import Entrez
from Bio import SeqIO
import os

# Set the email address to be used with Entrez
Entrez.email = "your_email@example.com"

def fetch_species_name(gene_id):
    handle = Entrez.efetch(db="nucleotide", id=gene_id, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    try:
        return records[0]["GBSeq_organism"]
    except IndexError:
        return "Unknown species"

def annotate_fasta(fasta_file):
    output_file = os.path.splitext(fasta_file)[0] + "_annotated.fasta"
    with open(fasta_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            gene_id = record.id.split("|")[1]
            species_name = fetch_species_name(gene_id)
            new_description = f"{record.description} {species_name}"
            record.description = new_description
            SeqIO.write(record, output_handle, "fasta")

# Directory containing your fasta files
input_directory = "./fasta_files"
output_directory = "./annotated_fasta_files"

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Process each fasta file in the input directory
for fasta_file in os.listdir(input_directory):
    if fasta_file.endswith(".fasta"):
        input_path = os.path.join(input_directory, fasta_file)
        output_path = os.path.join(output_directory, os.path.splitext(fasta_file)[0] + "_annotated.fasta")
        annotate_fasta(input_path)

print("Annotation complete.")

