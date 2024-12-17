#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J barrnap_run
#SBATCH --time 01:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e barrnap_run%j.err
#SBATCH -o barrnap_run%j.log

module purge
module load barrnap/0.9-GCC-9.2.0

# Date stamp for output directory
date_stamp=$(date +%Y%m%d)
output_dir="barrnap_output_$date_stamp"
mkdir -p "$output_dir"

# Find all .fna files in the specified directory
find ./prokka_4 -type f -name "*.fna" | while read fna_file; do
    file_basename=$(basename -- "$fna_file" .fna)

    # Create directory for this file's Barrnap output
    mkdir -p "$output_dir/$file_basename"

    # Run Barrnap on the .fna file
    barrnap "$fna_file" --outseq "$output_dir/$file_basename/${file_basename}_rrna.fa" > "$output_dir/$file_basename/${file_basename}_barrnap.gff"
done
