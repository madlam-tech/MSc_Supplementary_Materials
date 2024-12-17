#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH --job-name      BLAST
#SBATCH --time          01:00:00  # ~10 CPU minutes / MB blastn query vs nt
#SBATCH --mem           3G
#SBATCH -e 31-blastn_%j

mkdir -p ./blastn-clone-v-clone

module purge
module load BLAST/2.9.0-gimkl-2018b


# Directory containing FASTA files
fasta_directory="./blastn-clone-v-clone"

# Output directory for BLAST results
output_directory="./blastn-clone-v-clone"



# Perform BLAST comparisons for each pair of FASTA files
fasta_files=("$fasta_directory"/*.fasta)

for (( i=0; i<${#fasta_files[@]}; i++ )); do
    for (( j=i+1; j<${#fasta_files[@]}; j++ )); do
        query="${fasta_files[i]}"
        subject="${fasta_files[j]}"
        output="${output_directory}/blast_result_${i}_${j}.txt"

        blastn -query "$query" -subject "$subject" -out "$output" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    done
done

