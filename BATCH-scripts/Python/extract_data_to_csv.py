import csv
import sys

def parse_xmfa(file_path):
    blocks = []
    current_block = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('='):
                if current_block:
                    blocks.append(current_block)
                    current_block = []
            else:
                current_block.append(line.strip())
        if current_block:
            blocks.append(current_block)
    return blocks

def write_csv(blocks, output_path):
    with open(output_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Block ID', 'Sequence ID', 'Sequence'])
        block_id = 1
        for block in blocks:
            seq_id = None
            for line in block:
                if line.startswith('>'):
                    seq_id = line.split()[0][1:]
                elif seq_id is not None:
                    seq = line
                    csv_writer.writerow([block_id, seq_id, seq])
            block_id += 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_data_to_csv.py <input_xmfa_file> <output_csv_file>")
        sys.exit(1)
    
    input_xmfa_file = sys.argv[1]
    output_csv_file = sys.argv[2]
    
    blocks = parse_xmfa(input_xmfa_file)
    write_csv(blocks, output_csv_file)

    print(f"CSV file generated: {output_csv_file}")
