#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH --job-name      BLAST
#SBATCH --time          00:30:00  # ~10 CPU minutes / MB blastn query vs nt
#SBATCH --mem           3G
#SBATCH -e 34-blastn_%j

module purge
module load BLAST/2.9.0-gimkl-2018b

BLASTOPTS="-evalue 0.5 -max_target_seqs 30"
for filename in *.fasta
do 
        blastn $BLASTOPTS -query ${filename} \
	 -subject ma200_1_flye.fna \
	 -out ${filename}-v-flye-ma200.txt \
	 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
done

