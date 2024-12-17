#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryMM2
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryMM%j.err

module purge 
module load minimap2/2.24-GCC-9.2.0 
# Input files
consensus_fasta="2022-06-09-BC08-ma102_a_cl3.fasta"
raw_reads_fastq="BC08.fq"

# Output files
mapping_paf="mapping.paf"

# Mapping parameters
minimap2_params="-ax map-ont -t 8"

# Map consensus to raw reads
minimap2 $minimap2_params -k 27 --sam-hit-only $consensus_fasta $raw_reads_fastq > $mapping_paf
