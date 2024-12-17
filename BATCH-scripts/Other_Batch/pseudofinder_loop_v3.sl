#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J pseudofinder
#SBATCH --time 16:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e pseudofinder_%j.err
#SBATCH -o pseudofinder_%j.log

module load DIAMOND/2.1.6-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-04
module load Python/3.10.5-gimkl-2022a

# Ensure the directory structure is clean before starting
rm -rf pseudogene_output

# Create a directory to store output for all pseudogene files
mkdir pseudogene_output

# Loop through each GenBank file (*.gbk)
for filename in ../*/*.gbk; do
    # Check if the file exists
    if [ -f "$filename" ]; then
        # Extract basename without extension
        basename=$(basename "$filename" .gbk)
        
        # Create a directory for each file processed
        mkdir -p "pseudogene_output/$basename"

        # Run pseudofinder.py for the current file
        pseudofinder.py annotate --diamond --skip_makedb -g "$filename" -db /nesi/nobackup/massey03345/uniprotreference.dmnd -op "pseudogene_output/$basename"
    fi
done

# Move FNA files to the appropriate directory
mv ../*/*.fna pseudogene_output/

