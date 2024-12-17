#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          00:20:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1
#SBATCH --error         prod%j.err



# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

for dir in ./clusters/cluster_00*;

do

mkdir -p ./clusters//${dir}/predictions

prodigal -i ./clusters/${dir}/8_medaka.fasta -o ./clusters/${dir}/predictions/8_medaka_genes.fa -a ./clusters/${dir}/predictions/8_medaka.proteins.faa -p single

done

