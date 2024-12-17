#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryMDK
#SBATCH --time 00:30:00
#SBATCH --mem 10GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e 10-tryMDK%j.err
#SBATCH --export NONE

threads=32


module purge
module load medaka/1.6.0-Miniconda3-4.12.0


for c in ./cluster_*; do
    medaka_consensus -i "$c"/4_reads.fastq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r941_min_sup_g507 -t 12
    mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
    rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi  # clean up
done

