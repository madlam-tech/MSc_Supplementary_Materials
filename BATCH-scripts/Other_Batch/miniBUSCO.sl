#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      BUSCO
#SBATCH --time          12:00:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 4
#SBATCH --output        BUSCO_%j.out
#SBATCH --error         BUSCO_%j.err

# Load modules
module purge
module load miniBUSCO/0.2.1-gimkl-2022a
module load miniprot/0.11-GCC-11.3.0

# Directory to store BUSCO output directories
output_dir="/nesi/nobackup/massey03345/fnas_4_mauve"

# Loop through each directory matching the patterns you provided
for dir in /nesi/nobackup/massey03345/fnas_4_mauve/*/prokka; do
    # Check if the directory contains .faa files
    if [ -n "$(find "$dir" -maxdepth 1 -name '*.faa')" ]; then
        # Extract the strain name from the directory path
        strain_name=$(basename "$(dirname "$dir")")
        
        # Extract the parent directory (e.g., GCA_000300095.1) from the directory path
        parent_dir=$(basename "$(dirname "$(dirname "$dir")")")

        # Create the output directory for BUSCO results for the strain within the parent directory
        mkdir -p "$output_dir/$parent_dir/busco/$strain_name"

        # Run BUSCO for the current strain
        busco -i "$dir"/*.faa -o "$strain_name" -l ./betaproteobacteria_odb10 -c "$SLURM_CPUS_PER_TASK" -m proteins --out_path "$output_dir/$parent_dir/busco/$strain_name/"
    fi
done

