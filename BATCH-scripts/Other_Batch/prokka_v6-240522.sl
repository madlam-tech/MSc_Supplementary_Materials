#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka_v5
#SBATCH --time 03:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8

# Date stamp for error and output files
date_stamp=$(date +%Y%m%d)
#SBATCH -e 12-prokka_v2${date_stamp}%j.err
#SBATCH -o 12-prokka_v2${date_stamp}%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0
module load RNAmmer/1.2-GCC-9.2.0-Perl-5.30.1

# Find all .fna files in the specified directory and below
find /nesi/nobackup/massey03345/fnas_4_annotation/alt_genomes -name '*.fna' | while read -r fna_file; do
    parent_dir=$(dirname "$fna_file")
    file_basename=$(basename -- "$fna_file" .fna)

    # Run Prokka with RNAmmer on the .fna file
    prokka "$fna_file" \
    --force \
    --outdir "$parent_dir" \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70 \
    --rnammer
done

