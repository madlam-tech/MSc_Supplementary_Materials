import os
import re
import csv

def parse_xmfa(xmfa_file):
    block_id = None
    sequence_id = None
    sequence = ""
    start = None
    end = None
    data = []

    with open(xmfa_file, 'r') as xmfa:
        for line in xmfa:
            line = line.strip()
            if line.startswith(">"):
                if block_id and sequence:
                    data.append([block_id, sequence_id, start, end, sequence])
                block_info = line.split()
                block_id = re.search(r'\d+', block_info[0]).group()  # Extract block ID
                coords = re.search(r'\d+:\d+-\d+', block_info[0]).group().split(":")[1]
                start, end = map(int, coords.split("-"))
                sequence_id = block_info[-1].split('/')[-1]
                sequence = ""
            elif not line.startswith("#") and line != "=":
                sequence += line
        if block_id and sequence:
            data.append([block_id, sequence_id, start, end, sequence])

    return data

def save_to_csv(data, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Block ID", "Sequence ID", "Start", "End", "Sequence"])
        writer.writerows(data)

# Directory containing XMFA files
directory = './'  # Replace with the path to your directory if different

# Loop through all files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".xmfa"):
        xmfa_file = os.path.join(directory, filename)
        parsed_data = parse_xmfa(xmfa_file)
        csv_filename = os.path.splitext(filename)[0] + "_parsed.csv"
        csv_file = os.path.join(directory, csv_filename)
        save_to_csv(parsed_data, csv_file)
        print(f"Parsed and saved: {csv_file}")
