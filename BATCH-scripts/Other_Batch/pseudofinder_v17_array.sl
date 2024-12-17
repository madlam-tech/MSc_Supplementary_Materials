#!/bin/bash -e
#SBATCH --account=massey03345
#SBATCH --job-name=pseudofinder
#SBATCH --time=1:20:00
#SBATCH --mem=1GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --error=pseudofinder_%j_%a.err
#SBATCH --output=pseudofinder_%j_%a.log
#SBATCH --array=0-5

# Load required modules
module load DIAMOND/2.1.6-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-04
module load Python/3.10.5-gimkl-2022a

# Input directory directly set to the specified path
input_dir="./prokka_4"

# Get the current date in YYYY-MM-DD format
current_date=$(date +%F)

# Define the output directory with date-stamped
output_dir="./pseudogene_output_${current_date}"

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
num_jobs=21

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
    filename="${all_files[$I]}"
    # Extract basename without extension
    basename=$(basename "$filename" .gbk)

    # Create a directory for each file processed
    mkdir -p "$output_dir/$basename"

    # Run pseudofinder.py for the current file
    ./pseudofinder.py annotate --diamond --skip_makedb --compliant -g "$filename" -db /nesi/nobackup/massey03345/uniprot/uniprot_sprot.dmnd -op "$output_dir/$basename"
done

echo "Task $SLURM_ARRAY_TASK_ID completed."
