import csv

def parse_xmfa(xmfa_file):
    """
    Parse the XMFA file and return a list of dictionaries with the parsed data.
    """
    data = []
    block_id = 1  # Initialize block_id counter
    with open(xmfa_file, 'r') as file:
        for line in file:
            if line.startswith('='):
                block_id += 1  # Increment block_id for each new block
            elif line.startswith('>'):
                parts = line.strip().split()
                seq_id = parts[3]  # Sequence ID is always the 4th part
                # Handle different possible positions for start and end
                if ':' in parts[1]:
                    start, end = map(int, parts[1].split(':')[1].split('-'))
                else:
                    start, end = map(int, parts[1].split('-'))
                
                strand = parts[2] if parts[2] in ['+', '-'] else parts[3]

                current_block = {
                    "Sequence_ID": seq_id,
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "LCB_ID": block_id  # Use 1-based indexing for LCB_ID
                }
                data.append(current_block)
    return data

def summarize_lcb_data(lcb_data):
    """
    Summarize the LCB data.
    """
    summarized_data = []
    for entry in lcb_data:
        summarized_entry = {
            "Sequence_ID": entry['Sequence_ID'],
            "Start": entry['Start'],
            "End": entry['End'],
            "Strand": entry['Strand'],
            "LCB_ID": entry['LCB_ID']  # LCB_ID is already 1-based from parsing
        }
        summarized_data.append(summarized_entry)
    return summarized_data

def write_summary_to_csv(summarized_data, output_file):
    """
    Write the summarized LCB data to a CSV file.
    """
    with open(output_file, 'w', newline='') as file:
        fieldnames = ["Sequence_ID", "Start", "End", "Strand", "LCB_ID"]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        
        writer.writeheader()
        for entry in summarized_data:
            writer.writerow(entry)

if __name__ == "__main__":
    input_xmfa_file = "./combined.xmfa"
    output_summary_file = "./Paraburkholderia_lcbs_summary.csv"

    # Parse the XMFA file
    lcb_data = parse_xmfa(input_xmfa_file)
    
    # Summarize the LCB data
    summarized_data = summarize_lcb_data(lcb_data)
    
    # Write the summarized data to a CSV file
    write_summary_to_csv(summarized_data, output_summary_file)
    
    print(f"Summarized data written to {output_summary_file}")
