 cat fetch_species.py
from Bio import Entrez

# Set the email address to be used with Entrez
Entrez.email = "madlam@me.com"

# List of accession numbers
accession_numbers = [
    "GCA_000013425.1", "GCA_001611955.1", "GCA_002208805.2", "GCA_003571725.1",
    "GCA_003812505.1", "GCA_006094375.1", "GCA_016599795.1", "GCF_000006765.1",
    "GCF_000018225.1", "GCF_000063585.1", "GCF_000195955.2", "GCF_000225955.1",
    "GCF_000246755.1", "GCF_000426725.1", "GCF_000709415.1", "GCF_003387165.1",
    "GCF_004359515.1", "GCF_007814115.1", "GCF_008330805.1", "GCF_009731575.1",
    "GCF_011604685.1", "GCF_013030075.1", "GCF_013374215.1", "GCF_017821535.1",
    "GCF_025272815.1", "GCF_900143275.1", "GCF_900215245.1", "GCF_900458255.1",
    "ma1263_1"
]

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


