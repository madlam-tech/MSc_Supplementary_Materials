def parse_xmfa(file_path, output_csv):
    with open(file_path, 'r') as infile:
        block_id = None
        sequence_id = None
        sequence = ""
        results = []
        
        for line in infile:
            line = line.strip()
            
            # Debug information: print each line
            print(f"Processing line: {line}")

            if line.startswith(">"):
                # If we're in the middle of a block, save it
                if block_id is not None and sequence:
                    print(f"Adding block: {block_id}, sequence_id: {sequence_id}")
                    results.append([block_id, sequence_id, sequence])
                
                # Start a new block
                block_id_part = line.split()[0]
                block_id = block_id_part.replace(">", "").split(":")[0]
                sequence_id = line.split()[-1].split('/')[-1]
                sequence = ""  # Reset sequence for the new block
                
            elif not line.startswith("#"):
                # Accumulate sequence data
                sequence += line
                
        # Don't forget to add the last block
        if block_id is not None and sequence:
            print(f"Adding final block: {block_id}, sequence_id: {sequence_id}")
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
parse_xmfa('./combined.xmfa', 'parsed_xmfa_debug.csv')
