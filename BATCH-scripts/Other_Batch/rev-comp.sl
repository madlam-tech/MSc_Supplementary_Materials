#!/bin/bash

input_file="ma626_1_1_.fasta"
output_file="ma626_1_1_rev.fasta"

# Function to compute the reverse complement of a DNA sequence
revcomp() {
  echo "$1" | rev | tr "ACTGactg" "TGACtgac"
}

# Initialize variables
header=""
sequence=""

# Read the input FASTA file line by line
while IFS= read -r line; do
  # Check if the line is a header line
  if [[ $line =~ ^\> ]]; then
    # Process the previous sequence if it exists
    if [[ ! -z $sequence ]]; then
      echo -e "$header"
      echo -e "$(revcomp "$sequence")"
    fi

    # Store the new header and reset the sequence
    header="$line"
    sequence=""
  else
    # Concatenate the sequence lines
    sequence+="$line"
  fi
done < "$input_file"

# Process the last sequence
if [[ ! -z $sequence ]]; then
  echo -e "$header"
  echo -e "$(revcomp "$sequence")"
fi > "$output_file"

