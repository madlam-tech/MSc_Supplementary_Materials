#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka_array
#SBATCH --time 06:30:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --array=1-89%15  # Adjust the array size and limit to 15 concurrent jobs
#SBATCH -e 12-prokka_%A_%a.err
#SBATCH -o 12-prokka_%A_%a.log

# Load required modules
module purge
module load prokka/1.14.5-GCC-9.2.0

# Get the file corresponding to this array task
fna_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fna_files.txt)

# Extract basename without extension
file_basename=$(basename "$fna_file" .fna)

# Define the output directory for this file
output_dir="prokka_3/$file_basename"

# Check if the output directory already exists
if [ -d "$output_dir" ]; then
    echo "Output for $fna_file already exists. Skipping..."
else
    # Create directory for this file's output
    mkdir -p "$output_dir"

    # Run Prokka on the .fna file
    prokka --compliant "$fna_file" \
    --outdir "$output_dir" \
    --force \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70
fi

echo "Completed Prokka annotation for $fna_file"
