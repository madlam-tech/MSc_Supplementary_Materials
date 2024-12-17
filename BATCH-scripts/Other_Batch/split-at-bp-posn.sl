#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J split
#SBATCH --time 00:05:00
#SBATCH --mem 2GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e split%j.err
#SBATCH --export NONE


# Input FASTA file
input_file="ma626_1_1_rev.fasta"

# Base position to split the sequence
split_position=2353418

# Output file names
output_file1="ma626_1_1_rev-5p.fasta"
output_file2="ma626_1_1_rev-3p.fasta"

# Variables for sequence data and current position
sequence=""
position=0

# Read the input FASTA file
while IFS= read -r line; do
    if [[ $line =~ ^[>] ]]; then
        # Process sequence identifier line
        if [[ -n $sequence ]]; then
            # Split the sequence at the desired position
            substring1=${sequence:0:split_position}
            substring2=${sequence:split_position}
            
            # Write to output files
            echo "$header" >> "$output_file1"
            echo "$substring1" >> "$output_file1"
            
            echo "$header" >> "$output_file2"
            echo "$substring2" >> "$output_file2"
        fi
        
        # Save the new sequence header
        header="$line"
        sequence=""
    else
        # Concatenate sequence lines
        sequence+="$line"
    fi
done < "$input_file"

# Process the last sequence
if [[ -n $sequence ]]; then
    substring1=${sequence:0:split_position}
    substring2=${sequence:split_position}
    
    echo "$header" >> "$output_file1"
    echo "$substring1" >> "$output_file1"
    
    echo "$header" >> "$output_file2"
    echo "$substring2" >> "$output_file2"
fi
