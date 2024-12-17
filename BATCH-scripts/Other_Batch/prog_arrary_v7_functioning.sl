#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      john_testing_Mauve
#SBATCH --time          0:25:00
#SBATCH --array         1-253  # 23 sequences = 253 pairs
#SBATCH --mem           5GB
#SBATCH --cpus-per-task 2
#SBATCH --output        prog_mauve_array_%j.out
#SBATCH --error         prog_mauve_array_%j.err

module purge 2>/dev/null

# Define the directory containing the sequence files
directory="/nesi/nobackup/massey03345/fnas_4_mauve/renamed_originals"

# Define the output directory
output_dir="/nesi/nobackup/massey03345/fnas_4_mauve/renamed_originals/output"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Get the list of sequence files
files=("$directory"/*.fna)

# Initialize the total comparison count
total_comparisons=0

# Get the sequence file indices from the SLURM array task ID
index1=$((SLURM_ARRAY_TASK_ID % ${#files[@]}))
index2=$((SLURM_ARRAY_TASK_ID / ${#files[@]}))

# Extract the corresponding sequence file paths
sequence_file1="${files[$index1]}"
sequence_file2="${files[$index2]}"

# Extract the sequence names from the file paths
sequence_name1=$(basename "$sequence_file1" .fna)
sequence_name2=$(basename "$sequence_file2" .fna)

# Define the unique temporary directory for this comparison
temp_dir=$(mktemp -d)

# Copy the sequence files to the temporary directory
cp "$sequence_file1" "$sequence_file2" "$temp_dir"

# Copy the progressiveMauve directory to the temporary directory
cp -r ./progressiveMauve "$temp_dir"

# Change directory to the temporary directory
cd "$temp_dir" || exit

# Define the output file names using sequence names
output_file="${sequence_name1}_${sequence_name2}.xmfa"
guide_tree_file="${sequence_name1}_${sequence_name2}_tree.txt"

# Run progressiveMauve
./progressiveMauve  --output="$output_file" --skip-gapped-alignment --output-guide-tree="$guide_tree_file" --seed-family "$sequence_file1" "$sequence_file2"

# Move the relevant output files to the output directory
mv "$output_file" "$guide_tree_file" "$output_dir"

# Change back to the original directory
cd - || exit

# Remove the temporary directory
rm -r "$temp_dir"

# Display the total number of comparisons
echo "Total number of comparisons: $SLURM_ARRAY_TASK_COUNT"

