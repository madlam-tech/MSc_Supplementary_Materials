#!/bin/bash -e
#SBATCH --account uoa04147
#SBATCH --job-name BLAST_5S_rRNA
#SBATCH --time 00:05:00  # Adjust time based on expected workload
#SBATCH --mem 10G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e BLASTN%j.err
#SBATCH --array=1-98

module purge
module load BLAST/2.10.0-GCC-9.2.0
module load BLASTDB/2024-07

# Directory containing GCF*.fna files
DIR="./"

# Define the 5S rRNA sequence query file
QUERY_FILE="./Pb_ag_5S.fa"

# Find all GCF*.fna files and store them in an array
FASTAS=($(find $DIR -name "GCF*.fna"))

# Get the current FASTA file for this array job
TARGET=${FASTAS[$SLURM_ARRAY_TASK_ID-1]}

# Output format and BLAST options
FORMAT="6 qseqid qstart qend qseq sseqid sgi sacc sstart send staxids sscinames stitle length evalue bitscore"
BLASTOPTS="-evalue 0.05 -max_target_seqs 10"
BLASTAPP=blastn

# Create a BLAST database for the current target file
makeblastdb -in $TARGET -dbtype nucl -out ${TARGET%.*}_db

# Run BLAST using the 5S rRNA sequence as the query against the target genome
$BLASTAPP $BLASTOPTS -db ${TARGET%.*}_db -query $QUERY_FILE -outfmt "$FORMAT" \
    -out ${TARGET%.*}_5s_rRNA_results.txt -num_threads $SLURM_CPUS_PER_TASK

# Optional: Clean up BLAST database files to save space
rm ${TARGET%.*}_db.*
