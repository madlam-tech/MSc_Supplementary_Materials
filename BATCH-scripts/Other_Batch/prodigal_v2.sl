#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          00:10:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1
#SBATCH --error         prod.err
#SBATCH --output       prod.out

# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

mkdir -p predictions/

for filename in *.fna

do 

prodigal -i ${filename} -o ./predictions/${filename}_genes.fa -a ./predictions/${filename}_proteins.faa -p single

done
