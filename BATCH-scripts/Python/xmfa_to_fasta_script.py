import sys

def xmfa_to_fasta(xmfa_file):
    with open(xmfa_file, 'r') as f:
        lines = f.readlines()
    
    fasta_output = []
    current_seq = ""
    current_header = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                fasta_output.append(current_seq)
                current_seq = ""
            current_header = line
            fasta_output.append(current_header)
        elif not line.startswith('#') and line != '=':
            current_seq += line

    if current_seq:
        fasta_output.append(current_seq)

    return "\n".join(fasta_output)

if __name__ == "__main__":
    xmfa_file = sys.argv[1]
    fasta_output = xmfa_to_fasta(xmfa_file)
    print(fasta_output)
