#!/usr/bin/env python3
import os
import sys
import pandas as pd

def count_sequences(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        return len(lines)

def main(input_directory, output_tsv):
    results = []
    for filename in os.listdir(input_directory):
        if filename.endswith("_5s_rRNA_results.txt"):
            file_path = os.path.join(input_directory, filename)
            sequence_count = count_sequences(file_path)
            results.append({
                "Filename": filename,
                "5S_rRNA_Count": sequence_count
            })
            print(f"Counted {sequence_count} 5S rRNA sequences in {filename}")

    # Save results to a TSV file
    df = pd.DataFrame(results)
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Results saved to {output_tsv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python count_5s_rRNA_sequences.py <input_directory> <output_tsv>")
    else:
        input_directory = sys.argv[1]
        output_tsv = sys.argv[2]
        main(input_directory, output_tsv)
