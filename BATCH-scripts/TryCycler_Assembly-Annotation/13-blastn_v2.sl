#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J blastn
#SBATCH --time 01:00:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 13-blastn_%j.log
#SBATCH -o 13-blastn_%j.log2

module purge
module load BLAST/2.9.0-gimkl-2018b

for filename in *.fna

do

blastn -taxids 2 -query ${filename} -db /opt/nesi/db/blast/2023-04 -out ${filename}.blast -outfmt 0

done


