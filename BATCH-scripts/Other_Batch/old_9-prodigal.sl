#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          00:20:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1
#SBATCH --error         9-prod%j.log



# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

for dir in clusters/cluster_00*;

do

mkdir -p ./${dir}/predictions

prodigal -i ./${dir}/8_medaka.fasta -o ./${dir}/predictions/8_medaka_genes.fa -a ./${dir}/predictions/8_medaka.proteins.faa -p single

done
