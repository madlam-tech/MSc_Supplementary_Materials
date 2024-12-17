#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      quast
#SBATCH --time          01:00:00
#SBATCH --mem           20GB
#SBATCH --cpus-per-task 20
#SBATCH --error         quast.err
#SBATCH --output        quast.out

module purge
module load QUAST/5.2.0-gimkl-2022a

    quast.py --nanopore barcode07.fastq
               -r Pb_ag_GCF_009455635.1_ASM945563v1_genomic.fna, Pb_bo_GCF_009455625.1_ASM945562v1_genomic.fna, Pb_hay_GCF_009455685.1_ASM945568v1_genomic.fna \
               -g quast_test_data/BC07genes.gff

