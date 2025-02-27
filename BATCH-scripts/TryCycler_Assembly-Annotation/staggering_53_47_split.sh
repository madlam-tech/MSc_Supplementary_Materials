#!/bin/bash


# This shell script will remove 53% of the seqeunce at a predetermined start point (ion % of total script size
# this is to split the doublet contigs at differening start points


# Function to display usage information
usage() {
    echo "Usage: $0 <fasta_file> <start_percentage>"
    exit 1
}

# Check if the correct number of arguments was provided
if [ $# -ne 2 ]; then
    usage
fi

fasta_file=$1
start_percentage=$2

# Check if the file exists
if [ ! -f "$fasta_file" ]; then
    echo "File not found!"
    exit 1
fi

# Check if start_percentage is a valid number between 0 and 100
if ! [[ "$start_percentage" =~ ^[0-9]+$ ]] || [ "$start_percentage" -lt 0 ] || [ "$start_percentage" -gt 100 ]; then
    echo "Error: start_percentage must be a number between 0 and 100."
    exit 1
fi

# Extract the sequence from the FASTA file
sequence=$(grep -v '^>' "$fasta_file" | tr -d '\n')

# Get the length of the sequence
sequence_length=${#sequence}

# Calculate the starting point and new length
start_point=$(echo "$sequence_length * $start_percentage / 100" | bc)
new_length=$(echo "$sequence_length * 0.53 / 1" | bc)

# Extract the new sequence starting from the calculated start point
new_sequence=${sequence:$start_point:$new_length}

# Extract the header
header=$(grep '^>' "$fasta_file")

# Write the new sequence to a new file
output_file="split_${fasta_file}"
echo "$header" > "$output_file"
echo "$new_sequence" >> "$output_file"

echo "The new sequence has been saved to $output_file"

