#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH --job-name      BLAST
#SBATCH --time          00:40:00  # ~10 CPU minutes / MB blastn query vs nt
#SBATCH --mem           30G
#SBATCH -e 16-blastn_%j


module purge
module load BLAST/2.9.0-gimkl-2018b
module load BLASTDB/2023-04


# This script takes one argument, the FASTA file of query sequences.
FORMAT="6 qseqid qstart qend qseq sseqid sgi sacc sstart send staxids sscinames stitle length evalue bitscore"
BLASTOPTS="-evalue 0.05 -max_target_seqs 10"
BLASTAPP=blastn
DB=nt
#BLASTAPP=blastx
#DB=nr

$BLASTAPP $BLASTOPTS -db $DB -query cluster_001_consensus.fasta -outfmt "$FORMAT" \
    -out ma1263_1_1.$DB.$BLASTAPP -num_threads $SLURM_CPUS_PER_TASK
