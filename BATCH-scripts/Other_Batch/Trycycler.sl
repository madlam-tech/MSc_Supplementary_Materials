#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryMDK_v1
#SBATCH --time 4:00:00
#SBATCH --mem 20GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryMDK_v1.err
#SBATCH -o tryMDK_v1.out
#SBATCH --export NONE

threads=32


module purge
module load medaka/1.6.0-Miniconda3-4.12.0

for c in ./cluster_*; do
    medaka_consensus -i ../*fq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r941_min_hac_g507 -t 2
    mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
    rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi  # clean up
done


