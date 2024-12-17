import os
import glob
from Bio import SeqIO
import pandas as pd

# Define the path pattern to search for GBK files
path_pattern = "./*/**/*.gbk"

# Initialize a list to store the results
results = []

# Function to count rRNA genes in a GBK file
def count_rrna_genes(gbk_file):
    count_16S = 0
    count_23S = 0
    count_5S = 0
    print(f"Processing file: {gbk_file}")
    for record in SeqIO.parse(gbk_file, "genbank"):
        print(f"  Processing record: {record.id}")
        for feature in record.features:
            if feature.type == "rRNA":
                product = feature.qualifiers.get("product", [""])[0]
                print(f"    Found rRNA feature: {product}")
                if "16S" in product:
                    count_16S += 1
                    print("      Incrementing 16S count")
                elif "23S" in product:
                    count_23S += 1
                    print("      Incrementing 23S count")
                elif "5S" in product:
                    count_5S += 1
                    print("      Incrementing 5S count")
    return count_16S, count_23S, count_5S

# Find all GBK files and process each
print("Searching for GBK files...")
for gbk_file in glob.glob(path_pattern, recursive=True):
    count_16S, count_23S, count_5S = count_rrna_genes(gbk_file)
    results.append({
        "File": os.path.basename(gbk_file),
        "16S rRNA": count_16S,
        "23S rRNA": count_23S,
        "5S rRNA": count_5S
    })

# Create a DataFrame from the results
df = pd.DataFrame(results)

# Save the DataFrame to a CSV file
output_csv = "rRNA_gene_counts.csv"
df.to_csv(output_csv, index=False)

print("Results saved to", output_csv)
print(df)
