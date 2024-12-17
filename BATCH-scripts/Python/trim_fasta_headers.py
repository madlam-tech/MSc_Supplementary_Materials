#!/usr/bin/env python3

import sys
import os

def trim_header(fasta_file, output_file):
    with open(fasta_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                new_header = line.split("_")[0] + "\n"
                outfile.write(new_header)
            else:
                outfile.write(line)

if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for fasta_file in os.listdir(input_dir):
        if fasta_file.endswith(".fasta"):
            input_path = os.path.join(input_dir, fasta_file)
            output_path = os.path.join(output_dir, fasta_file)
            print(f"Processing file: {input_path}")
            trim_header(input_path, output_path)
