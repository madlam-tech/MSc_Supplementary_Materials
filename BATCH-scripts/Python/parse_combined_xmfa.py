def parse_xmfa(file_path, output_csv):
    with open(file_path, 'r') as infile:
        block_id = None
        sequence_id = None
        sequence = ""
        results = []
        
        for line in infile:
            line = line.strip()
            
            # Identify the start of a block
            if line.startswith(">"):
                # If we're in the middle of a block, save it
                if block_id is not None and sequence:
                    results.append([block_id, sequence_id, sequence])
                
                # Start a new block
                block_info = line.split()[0].replace(">", "")
                block_id = block_info.split(":")[0]
                sequence_id = line.split()[-1].split('/')[-1]
                sequence = ""  # Reset sequence for the new block
                
            elif not line.startswith("#"):
                # Accumulate sequence data
                sequence += line
                
        # Don't forget to add the last block
        if block_id is not None and sequence:
            results.append([block_id, sequence_id, sequence])
        
    if not results:
        print("No blocks were found. Please check the format of the file.")
    else:
        # Write to CSV
        import csv
        with open(output_csv, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(["Block ID", "Sequence ID", "Sequence"])
            writer.writerows(results)

        print(f"CSV file generated: {output_csv}")

# Run the parser
parse_xmfa('./combined.xmfa', 'parsed_xmfa.csv')
