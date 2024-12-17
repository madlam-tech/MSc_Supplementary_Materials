#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J diamondDB
#SBATCH --time 10:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e makedb_%j.log

module load DIAMOND/2.1.6-GCC-11.3.0
module load BLASTDB/2023-04
module load BLAST/2.13.0-GCC-11.3.0
module load Python/3.10.5-gimkl-2022a


diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot

