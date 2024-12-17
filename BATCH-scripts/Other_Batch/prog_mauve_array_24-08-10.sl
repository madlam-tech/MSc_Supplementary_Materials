#!/bin/bash -e

#SBATCH --account       uoa04147
#SBATCH --job-name      progressiveMauve
#SBATCH --time          01:00:00
#SBATCH --mem           8G
#SBATCH --cpus-per-task 2
#SBATCH --output        ./job_%A_%a.out
#SBATCH --error         ./job_%A_%a.err
#SBATCH --array         1-36  # Set the correct range for your sequences

module purge 2>/dev/null
module load Java

# Define the input directory
input_dir="./renamed"

# Get the list of sequence files in the input directory
files=("$input_dir"/sequence_*.fna)
num_seq_files=${#files[@]}

# Dynamically calculate the number of combinations for the array range
if (( num_seq_files < 2 )); then
    echo "ERROR: Not enough sequence files. At least 2 are required."
    exit 1
fi

# Calculate the total number of pairwise comparisons (choose 2 from n)
total_comparisons=$((num_seq_files * (num_seq_files - 1) / 2))

# Ensure the SLURM array is properly sized
if (( SLURM_ARRAY_TASK_ID > total_comparisons )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID exceeds the number of comparisons."
    exit 1
fi

# Convert the array task ID to a pair of indices
task_id=0
for (( i = 0; i < num_seq_files - 1; i++ )); do
    for (( j = i + 1; j < num_seq_files; j++ )); do
        task_id=$((task_id + 1))
        if (( task_id == SLURM_ARRAY_TASK_ID )); then
            sequence_file1="${files[$i]}"
            sequence_file2="${files[$j]}"
            break 2
        fi
    done
done

# Extract sequence names from the file paths
sequence_name1=$(basename "$sequence_file1" .fna)
sequence_name2=$(basename "$sequence_file2" .fna)

# Define the output file names using sequence names
output_file="alignment_${sequence_name1}_${sequence_name2}.xmfa"
output_guide_file="guide_tree_${sequence_name1}_${sequence_name2}.txt"

# Run progressiveMauve
./progressiveMauve --output="$output_file" --output-guide-tree="$output_guide_file" --seed-family "$sequence_file1" "$sequence_file2"

echo "Finished comparison of $sequence_name1 and $sequence_name2."
