#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryMM2
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e  tryMM%j.err
#SBATCH --export NONE

module purge 
module load minimap2/2.24-GCC-9.2.0 
# Input files
consensus_fasta="final_chr3.fasta"
raw_reads_fastq="BC09.fq"

# Output files
mapping_sam="mapping_chr3.sam"

# Mapping parameters
minimap2_params="-ax map-ont -t 8"

# Map consensus to raw reads
minimap2 $minimap2_params $consensus_fasta ../filtered/$raw_reads_fastq > $mapping_sam
