def read_xmfa(file):
    with open(file, 'r') as f:
        return f.readlines()

def merge_xmfa(files, output_file):
    merged_data = []
    for file in files:
        xmfa_data = read_xmfa(file)
        for line in xmfa_data:
            if line.startswith('>'):
                # Modify the header if necessary
                merged_data.append(line)
            elif line.startswith('='):
                # End of block, add separator
                merged_data.append(line)
            else:
                # Sequence data
                merged_data.append(line)
    
    with open(output_file, 'w') as f:
        f.writelines(merged_data)

if __name__ == "__main__":
    import glob
    xmfa_files = glob.glob('*.xmfa')
    output_file = 'combined_alignment.xmfa'
    merge_xmfa(xmfa_files, output_file)
    print(f"Merged {len(xmfa_files)} XMFA files into {output_file}")
