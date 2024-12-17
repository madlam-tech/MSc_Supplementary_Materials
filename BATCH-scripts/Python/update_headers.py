import csv
import sys

# Read TSV file to create a mapping from accession number to species name
def read_species_mapping(tsv_file):
    mapping = {}
    with open(tsv_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            acc_number = row["File"].replace("_ordered", "")
            species_name = row["Species"]  # Ensure the correct column is used
            mapping[acc_number] = species_name
    return mapping

# Replace header in FASTA file with species name
def replace_header(fasta_file, species_name, output_file):
    with open(fasta_file, "r") as file:
        lines = file.readlines()
    if lines:
        if species_name == "Unknown_Species":
            new_header = f">Unknown_Species\n"
        else:
            genus, species = species_name.split(maxsplit=1)
            genus_initial = genus[0]
            new_header = f">{genus_initial}.{species.replace(' ', '_')}\n"
        lines[0] = new_header
    with open(output_file, "w") as file:
        file.writelines(lines)

if __name__ == "__main__":
    species_mapping = read_species_mapping(sys.argv[1])
    fasta_files = sys.argv[2:]
    for fasta_file in fasta_files:
        acc = fasta_file.split("/")[-1].split("_concatenated")[0]
        species_name = species_mapping.get(acc, "Unknown_Species")
        if acc.startswith("ma"):
            species_name = acc.upper()
        output_file = sys.argv[3] + "/" + fasta_file.split("/")[-1]
        print(f"Processing file: {fasta_file}")
        print(f"Extracted accession number: {acc}")
        print(f"Mapped species name: {species_name}")
        try:
            replace_header(fasta_file, species_name, output_file)
        except IndexError:
            print(f"Error: Unable to replace header for {fasta_file}")
