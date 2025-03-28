#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 10:00:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e prodigal-prokka%j.err
#SBATCH -o prodigal-prokka%j.log

# Load Prokka module
module purge
module load prokka/1.14.5-GCC-9.2.0

# Directory containing .faa files generated by Prodigal
faa_dir="/nesi/nobackup/massey03345/fnas_4_mauve/prodigal/"

# Loop through each .faa file
for faa_file in "$faa_dir"/*.faa; do
    # Extract file name without extension
    file_name=$(basename "$faa_file" .faa)
    
    # Create a directory using the file name
    mkdir -p "$file_name"_prokka
    
    # Run Prokka for the current .faa file
    prokka "$faa_file" \
        --outdir ./prodigal/"$file_name"_prokka \
        --force \
        --prefix prokka \
        --addgenes \
        --addmrna \
        --proteins "$faa_file"  # Specify the input protein sequences
done

