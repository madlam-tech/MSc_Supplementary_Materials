#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J pseudofinder
#SBATCH --time 1:20:00
#SBATCH --mem 2GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e pseudofinder_%j_%a.err
#SBATCH -o pseudofinder_%j_%a.log
#SBATCH --array 0-10

# Load required modules
module load DIAMOND/2.1.6-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-04
module load Python/3.10.5-gimkl-2022a

# Input directory directly set to the specified path
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

# Calculate the portion of files for this task
files_per_task=$(( (total_files + 20 - 1) / 20 ))

# Calculate the start and end index for this task
start_index=$(( SLURM_ARRAY_TASK_ID * files_per_task ))
end_index=$(( start_index + files_per_task - 1 ))

# Ensure the end index does not exceed the total number of files
if [ $end_index -ge $total_files ]; then
    end_index=$(( total_files - 1 ))
fi

# Process the files assigned to this task
for (( i = start_index; i <= end_index; i++ )); do
    filename="${all_files[$i]}"
    # Extract basename without extension
    basename=$(basename "$filename" .gbk)

    # Create a directory for each file processed
    mkdir -p "$output_dir/$basename"

    # Run pseudofinder.py for the current file
    ./pseudofinder.py annotate --diamond --skip_makedb --compliant -g "$filename" -db /nesi/nobackup/massey03345/uniprot/uniprot_sprot.dmnd -op "$output_dir/$basename"
done
