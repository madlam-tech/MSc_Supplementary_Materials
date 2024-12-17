import sys
from Bio import SeqIO

def extract_genes(gbk_path, output_path, genes_to_find):
    with open(gbk_path, 'r') as input_handle, open(output_path, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "gene" and "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                    if gene_name in genes_to_find:
                        # Extract the sequence
                        seq = feature.location.extract(record).seq
                        # Write to output
                        output_handle.write(f">{gene_name}\n{str(seq)}\n")

if __name__ == "__main__":
    gbk_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    genes_to_look_for = set(sys.argv[3].split('|'))  # genes passed as "gene1|gene2|gene3"
    extract_genes(gbk_file_path, output_file_path, genes_to_look_for)

