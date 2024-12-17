#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tRNAscan_run
#SBATCH --time 02:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e tRNAscan_run%j.err
#SBATCH -o tRNAscan_run%j.log

module purge
module load NeSI
module load tRNAscan-SE/2.0.7-GCC-9.2.0

# Directory containing .fna files
cd /nesi/nobackup/massey03345/fnas_4_annotation/all_fnas

# Loop through all .fna files in the directory
for fna_file in *.fna; do
    # Generate output file name by replacing .fna with _trnascan.gff
    output_gff="${fna_file%.fna}_trnascan.gff"

    # Run tRNAscan-SE on the file
    tRNAscan-SE -G -o "$output_gff" "$fna_file"

    echo "Finished processing $fna_file"
done

echo "All files processed."
