import os
import subprocess
from glob import glob

# Define directories
input_dir = '/nesi/nobackup/massey03345/fnas_complete_set/pseudogene_output_2024-06-13/'
output_dir = '/nesi/nobackup/massey03345/fnas_complete_set/merged_output/'
sorted_dir = os.path.join(output_dir, 'sorted')
merged_dir = os.path.join(output_dir, 'merged')

# Create directories if they don't exist
os.makedirs(sorted_dir, exist_ok=True)
os.makedirs(merged_dir, exist_ok=True)

# Function to extract basename
def extract_basename(filepath):
    return os.path.basename(filepath).split('_genomic')[0]

# Function to preprocess GFF files to ensure they are tab-delimited and have the correct number of columns
def preprocess_gff(gff_file, cleaned_file):
    with open(gff_file, 'r') as infile, open(cleaned_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                columns = line.strip().split('\t')
                if len(columns) >= 9:  # Ensure the line has at least 9 columns (standard GFF format)
                    outfile.write('\t'.join(columns) + '\n')

# Get all GFF files
gff_files = glob(os.path.join(input_dir, '*.gff'))

# Sort and merge GFF files by basename
basename_dict = {}

for gff_file in gff_files:
    basename = extract_basename(gff_file)
    if basename not in basename_dict:
        basename_dict[basename] = []
    basename_dict[basename].append(gff_file)

# Sort and merge files
for basename, files in basename_dict.items():
    cleaned_files = []
    sorted_files = []
    for file in files:
        cleaned_file = os.path.join(sorted_dir, os.path.basename(file) + '.cleaned')
        cleaned_files.append(cleaned_file)
        preprocess_gff(file, cleaned_file)
        
        sorted_file = cleaned_file + '.sorted'
        sorted_files.append(sorted_file)
        # Sort the file using bedtools
        subprocess.run(['bedtools', 'sort', '-i', cleaned_file, '-o', sorted_file])
    
    merged_file = os.path.join(merged_dir, f'{basename}_merged.gff')
    # Merge sorted files using bedtools
    subprocess.run(['bedtools', 'merge', '-i'] + sorted_files + ['-o', merged_file])

print("Sorting and merging completed.")
