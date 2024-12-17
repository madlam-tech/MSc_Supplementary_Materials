#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryMM2
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryMM%j.err
#SBATCH --export NONE

module purge 
module load minimap2/2.24-GCC-9.2.0 
# Input files
consensus_fasta="CP_+90.fasta"
raw_reads_fastq="BC10.fq"

# Output files
mapping_sam="mapping.sam"

# Mapping parameters
minimap2_params="-ax map-ont -t 8"

# Map consensus to raw reads
minimap2 $minimap2_params $consensus_fasta $raw_reads_fastq > $mapping_v_CP_+90.sam
