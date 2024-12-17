#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 1:00:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e prodigal-prokka%j.err
#SBATCH -o prodigal-prokka%j.log

# Load Prokka module
module purge
module load prokka/1.14.5-GCC-9.2.0

# Define paths to files
faa_file="/nesi/nobackup/massey03345/fnas_4_annotation/prodigal/ma1263_1_4_pseudo_genes.faa"
fna_file="/nesi/nobackup/massey03345/fnas_4_annotation/old_header_fixed_header_fnas/ma1263_1_4_pseudo.fna"

# Run Prokka for the specified files
prokka "$faa_file" \
    --outdir "/nesi/nobackup/massey03345/fnas_4_annotation/prodigal/ma1263_1_4_pseudo_prokka" \
    --force \
    --prefix prokka \
    --addgenes \
    --addmrna \
    --proteins "$faa_file"

prokka "$fna_file" \
    --outdir "/nesi/nobackup/massey03345/fnas_4_annotation/prodigal/ma1263_1_4_pseudo_prokka" \
    --force \
    --prefix prokka \
    --addgenes \
    --addmrna

