#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      rerun_missing_Mauve
#SBATCH --time          0:35:00
#SBATCH --mem           8GB
#SBATCH --cpus-per-task 2
#SBATCH --output        prog_mauve_rerun_%j.out
#SBATCH --error         prog_mauve_rerun_%j.err
#SBATCH --array         1-<MAX_ARRAY_TASK_ID> # This will be set dynamically based on the number of missing combinations

module purge 2>/dev/null

# Define the directory containing the sequence files
directory="/nesi/nobackup/massey03345/fnas_4_mauve/renamed_originals"

# Define the output directory
output_dir="/nesi/nobackup/massey03345/fnas_4_mauve/renamed_originals/output-24-05-27"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# File containing the missing combinations
missing_combinations="/nesi/nobackup/massey03345/fnas_4_mauve/renamed_originals/output-24-05-2/missing_combinations.txt"

# Get the total number of lines in the missing combinations file
total_lines=$(wc -l < "$missing_combinations")

# Adjust SLURM_ARRAY_TASK_ID to handle combinations in the missing_combinations.txt file
SLURM_ARRAY_TASK_ID=$((SLURM_ARRAY_TASK_ID % total_lines + 1))

# Read the specific combination for this task ID
combination=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$missing_combinations")
sequence_name1=$(echo $combination | cut -d' ' -f1)
sequence_name2=$(echo $combination | cut -d' ' -f2)

# Define the sequence file paths
sequence_file1="$directory/${sequence_name1}.fna"
sequence_file2="$directory/${sequence_name2}.fna"

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
