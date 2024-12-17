#!/bin/bash -e
#SBATCH --job-name=merge_annotations
#SBATCH --account=massey03345
#SBATCH --time=01:00:00
#SBATCH --mem=4GB
#SBATCH --cpus-per-task=4
#SBATCH --output=merge_annotations_%j.log
#SBATCH --error=merge_annotations_%j.err

# Load necessary modules
module purge
module load BEDTools/2.30.0-GCC-11.3.0
module load Python/3.8.2-gimkl-2020a

# Define directories
PROKKA_DIR="/nesi/nobackup/massey03345/fnas_complete_set/prokka/"
BARRNAP_DIR="/nesi/nobackup/massey03345/fnas_complete_set/barrnap_output_20240605/"
PSEUDOFINDER_DIR="/nesi/nobackup/massey03345/fnas_complete_set/pseudogene_output_2024-06-13/"
MERGED_DIR="/nesi/nobackup/massey03345/fnas_complete_set/merged_output/"

# Create merged output directory
mkdir -p $MERGED_DIR

# Find and extract annotations from Prokka, Barrnap, and Pseudofinder GFF files
find $PROKKA_DIR -name "*.gff" -exec grep -v "^#" {} \; > ${MERGED_DIR}/prokka_annotations.gff
find $BARRNAP_DIR -name "*.gff" -exec grep -v "^#" {} \; > ${MERGED_DIR}/barrnap_annotations.gff
find $PSEUDOFINDER_DIR -name "*.gff" -exec grep -v "^#" {} \; > ${MERGED_DIR}/pseudofinder_annotations.gff

# Merge annotations using bedtools
bedtools intersect -a ${MERGED_DIR}/prokka_annotations.gff -b ${MERGED_DIR}/barrnap_annotations.gff -wa -wb > ${MERGED_DIR}/prokka_barrnap_merged.gff
bedtools intersect -a ${MERGED_DIR}/prokka_barrnap_merged.gff -b ${MERGED_DIR}/pseudofinder_annotations.gff -wa -wb > ${MERGED_DIR}/merged_annotations.gff

# Run the Python script for further processing if needed
# python3 merge_annotations.py

echo "Final annotations merged and saved to ${MERGED_DIR}/merged_annotations.gff"
