#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 03:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 12-prokka_v2%j.err
#SBATCH -o 12-prokka_v2%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0

# Create directory for Prokka output
#mkdir -p prokka_output

# Find all _4_pseudo.fna files in the current directory
for fna_file in *.fna; do
    # Extract basename of the file (without extension)
    file_basename=$(basename -- "$fna_file" .fna)
    
    # Create directory for this file's Prokka output
    mkdir -p "prokka_output/$file_basename"
    
    # Run Prokka on the fna file
    prokka "$fna_file" \
    --force \
    --outdir "prokka_output/$file_basename" \
    --prefix "$file_basename" \
    --addgenes \
    --addmrna \
    --coverage 70
done

