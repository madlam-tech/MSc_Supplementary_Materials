import os
import re
import csv

def parse_xmfa(xmfa_file):
    blocks = []
    block_id = None
    sequence_id = None
    sequence_data = []

    print(f"Opening file: {xmfa_file}")
    with open(xmfa_file, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            line = line.strip()

            print(f"Processing line {line_num}: {line}")

            if line.startswith('>'):
                if sequence_data:
                    print(f"Appending block: Block ID = {block_id}, Sequence ID = {sequence_id}, Sequence Length = {len(''.join(sequence_data))}")
                    blocks.append((block_id, sequence_id, ''.join(sequence_data)))
                    sequence_data = []
                
                block_info = re.findall(r'>\s*(\d+):', line)
                if block_info:
                    block_id = block_info[0]
                    print(f"Found Block ID: {block_id}")
                else:
                    block_id = "unknown"  # Handle cases where block ID might not be present
                    print("Warning: Block ID not found, assigning 'unknown'")

                sequence_id = line.split()[-1]  # Extract Sequence ID from the last part of the line
                print(f"Found Sequence ID: {sequence_id}")
            elif line and not line.startswith('#') and not line.startswith('='):
                sequence_data.append(line)
                print(f"Appending sequence data, current length: {len(''.join(sequence_data))}")

        if sequence_data:
            print(f"Appending final block: Block ID = {block_id}, Sequence ID = {sequence_id}, Sequence Length = {len(''.join(sequence_data))}")
            blocks.append((block_id, sequence_id, ''.join(sequence_data)))

    print(f"Total blocks parsed: {len(blocks)}")
    return blocks

def save_to_csv(blocks, output_file):
    print(f"Saving parsed data to CSV file: {output_file}")
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Block ID", "Sequence ID", "Sequence"])
        for block in blocks:
            print(f"Writing Block ID = {block[0]}, Sequence ID = {block[1]}, Sequence Length = {len(block[2])}")
            writer.writerow(block)
    print(f"CSV file {output_file} has been created successfully.")

# Loop through all the XMFA files in the directory
directory = './'
print(f"Scanning directory: {directory} for XMFA files")
for filename in os.listdir(directory):
    if filename.endswith(".xmfa"):
        xmfa_file = os.path.join(directory, filename)
        output_file = xmfa_file.replace('.xmfa', '_parsed.csv')
        
        print(f"Processing XMFA file: {xmfa_file}")
        parsed_data = parse_xmfa(xmfa_file)
        
        if parsed_data:
            save_to_csv(parsed_data, output_file)
            print(f"CSV file generated: {output_file}")
        else:
            print(f"No data parsed from {xmfa_file}")

print("Processing complete.")
