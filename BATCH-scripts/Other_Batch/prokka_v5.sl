#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka_v5
#SBATCH --time 03:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka_v2%j.err
#SBATCH -o 12-prokka_v2%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0
module load RNAmmer/1.2-GCC-9.2.0-Perl-5.30.1

# Date stamp for output directory
date_stamp=$(date +%Y%m%d)
output_dir="prokka_output_$date_stamp"
mkdir -p "$output_dir"

# Find all .fna files in the specified directory
for fna_file in /nesi/nobackup/massey03345/fnas_4_annotation/all_fnas/*.fna; do
    file_basename=$(basename -- "$fna_file" .fna)
    
    # Create directory for this file's Prokka output
    mkdir -p "$output_dir/$file_basename"
    
    # Run Prokka with RNAmmer on the .fna file
    prokka "$fna_file" \
    --force \
    --outdir "$output_dir/$file_basename" \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70 \
    --rnammer
done

