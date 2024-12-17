#!/bin/bash -e
#SBATCH --account=massey03345
#SBATCH --job-name=pseudofinder
#SBATCH --time=0:20:00
#SBATCH --mem=4GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --error=pseudofinder_%j_%a.err
#SBATCH --output=pseudofinder_%j_%a.log
#SBATCH --array=0-91

# Load required modules
module load DIAMOND/2.1.6-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-04
module load Python/3.10.5-gimkl-2022a

# Define input and output directories
input_dir="./prokka_2"
output_dir="pseudogene_output"

# Ensure the directory structure is clean before starting
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    rm -rf "$output_dir"
    mkdir -p "$output_dir"
fi

# Wait for the first task to create the output directory
sleep 60

# Get the list of .gbk files
all_files=($(find "$input_dir" -type f -name "*.gbk"))
total_files=${#all_files[@]}

# Number of array jobs
num_jobs=92

# Calculate the portion of files for this task
files_per_task=$(( (total_files + num_jobs - 1) / num_jobs ))

# Calculate the start and end index for this task
start_index=$(( SLURM_ARRAY_TASK_ID * files_per_task ))
end_index=$(( start_index + files_per_task - 1 ))

# Ensure the end index does not exceed the total number of files
if [ $end_index -ge $total_files ]; then
    end_index=$(( total_files - 1 ))
fi

# Process the files assigned to this task
for (( I = start_index; I <= end_index; I++ )); do
    if [ $I -ge $total_files ]; then
        break
    fi
    file=${all_files[$I]}
    base_name=$(basename "$file" .gbk)
    output_file="$output_dir/${base_name}_pseudogenes.txt"
    
    # Run pseudofinder or any other required processing on the file
    pseudofinder "$file" > "$output_file"
done

echo "Task $SLURM_ARRAY_TASK_ID completed."
