#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      metaxa2
#SBATCH --time          03:00:00
#SBATCH --mem           10GB
#SBATCH --cpus-per-task 4
#SBATCH --error         metaxa_%j.err
#SBATCH --output        metaxa_%j.out

# Load modules
module purge
module load Metaxa2/2.2.3-gimkl-2022a

# Create output directory if it doesn't exist
mkdir -p ribosomes/

# Loop through each .fna file in the current directory
for fna_file in *.fna; do
    # Extract the filename without extension
    filename=$(basename "$fna_file" .fna)
    # Extract the directory name based on the filename pattern (assuming filename is in the format "GCA_000300095.1_ASM30009v1_genomic")
    dirname=$(echo "$filename" | cut -d '_' -f 1,2)
    # Create a new directory based on the extracted directory name
    mkdir -p "$dirname/metaxa"
    
    # Run Metaxa2 for each ribosome type
    for ribosome_type in ssu lsu; do
        metaxa2 -g "$ribosome_type" --mode genome \
                -i "$fna_file" \
                -o "$dirname/metaxa/${filename}_genes.fa.${ribosome_type}"
    done
done

