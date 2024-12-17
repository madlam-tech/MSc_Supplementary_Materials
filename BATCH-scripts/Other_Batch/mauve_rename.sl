#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      mauve_rename
#SBATCH --time          1:00:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 2
#SBATCH --output 	mauve_rename_%j.out
#SBATCH --error		mauve_rename_%j.err

# Define source directory
source_dir="/nesi/nobackup/massey03345/fnas_4_mauve_test"

# Create a text file to store the mappings
mapping_file="file_mapping.txt"

# Initialize counter
counter=1

# Loop through each file in the source directory
for file in "$source_dir"/*.fna; do
    # Extract the base name of the file
    base_name=$(basename "$file")

    # Generate the new file name
    new_name="sequence_${counter}.fna"

    # Rename the file
    mv "$file" "$source_dir/$new_name"

    # Store the mapping in the text file
    echo "$base_name -> $new_name" >> "$source_dir/$mapping_file"

    # Increment counter
    ((counter++))
done

