#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH --job-name checkm_run
#SBATCH --time 04:00:00
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 8
#SBATCH --error checkm_run%j.err
#SBATCH --output checkm_run%j.out

module purge
module load CheckM/1.2.1-gimkl-2022a-Python-3.10.5

# Define the main directory
main_dir="/nesi/nobackup/massey03345/fnas_4_annotation/staph/all_fnas"

# Navigate to the main directory
cd "$main_dir"

# Loop through .fna files in the current directory
for fna_file in *.fna; do
    echo "Processing file: $fna_file"
    
    # Extract basename from the .fna file
    file_basename=$(basename "$fna_file" .fna)

    # Check if the .faa file exists in the same directory
    faa_file="${file_basename}_proteins.faa"
    if [[ -f "$faa_file" ]]; then
        echo "Found .faa file: $faa_file"
        
        # Create the output directory for CheckM
        checkm_dir="$main_dir/checkM/${file_basename}/output"
        mkdir -p "$checkm_dir"
        
        # Run CheckM lineage_wf
        checkm lineage_wf --genes -t 8 --pplacer_threads 8 -x faa --tab_table -f "$checkm_dir/checkm.txt" "$main_dir" "$checkm_dir"
        echo "CheckM analysis completed for $fna_file"
    else
        echo "No .faa file found for $file_basename"
    fi
done

