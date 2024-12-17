#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka_array
#SBATCH --time 01:30:00
#SBATCH --mem 3GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --array=1-6%6  # Adjust the array size and limit to 30 concurrent jobs
#SBATCH -e 12-prokka_%A_%a.err
#SBATCH -o 12-prokka_%A_%a.log

# Load required modules
module purge
module load prokka/1.14.5-GCC-9.2.0
#module load barrnap/0.9-GCC-9.2.0

# Get the file corresponding to this array task
fna_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fna_files.txt)

# Extract basename without extension
file_basename=$(basename "$fna_file" .fna)

# Truncate the locus tag and prefix after the second underscore
truncated_basename=$(echo "$file_basename" | awk -F'_' '{print $1 "_" $2}')

# Ensure the prefix does not exceed 37 characters
if [ ${#truncated_basename} -gt 20 ]; then
    truncated_basename=${truncated_basename:0:20}
fi

# Define the output directory for this file
output_dir="prokka_4/$file_basename"

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
           --prefix "$truncated_basename" \
           --addgenes \
           --addmrna \
           --coverage 70 \
           --locustag "$truncated_basename" \
           --rfam

    # Run Barrnap on the .fna file
    # barrnap "$fna_file" --outseq "$output_dir/${file_basename}_rrna.fa" --compliant  > "$output_dir/${file_basename}_barrnap.gff"
fi

echo "Completed Prokka annotation for $fna_file"
