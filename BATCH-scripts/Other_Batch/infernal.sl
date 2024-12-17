#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J infernal_run
#SBATCH --time 04:00:00
#SBATCH --mem 2GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e infernal_run%j.err
#SBATCH -o infernal_run%j.log

module purge
module load Infernal/1.1.4-GCC-9.2.0

# Path to the rRNA covariance model, e.g., bacteria
rrna_model="/path/to/bacteria.cm"

# Date stamp for output directory
date_stamp=$(date +%Y%m%d)
output_dir="infernal_output_$date_stamp"
mkdir -p "$output_dir"

# Find all .fna files in the specified directory
for fna_file in /nesi/nobackup/massey03345/fnas_4_annotation/all_fnas/*.fna; do
    file_basename=$(basename -- "$fna_file" .fna)
    
    # Create directory for this file's Infernal output
    mkdir -p "$output_dir/$file_basename"
    
    # Run Infernal's cmsearch on the .fna file
    cmsearch --cpu 4 --tblout "$output_dir/$file_basename/$file_basename"_infernal.tbl "$rrna_model" "$fna_file" > "$output_dir/$file_basename/$file_basename"_infernal.out
done

