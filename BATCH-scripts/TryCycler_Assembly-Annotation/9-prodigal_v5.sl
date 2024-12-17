#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_v2
#SBATCH --time          00:30:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1
#SBATCH --error         9-prod%j.log

# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

for filepath in ./*fasta
do
    filename=$(basename "$filepath")
    filename="${filename%.*}"
    mkdir -p  "./$filename"
    mkdir -p  "./$filename"/predictions
    cp $filepath ./filename

prodigal -i ./$filepath -o ./$filename/predictions/$filename.fa -a ./$filename/predictions/$filename.faa -p single

done
