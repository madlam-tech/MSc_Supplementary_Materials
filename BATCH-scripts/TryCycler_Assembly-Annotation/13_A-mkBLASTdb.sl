#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mkBLASTdb
#SBATCH --time 01:00:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 13-mkBlastDB_%j.log
#SBATCH -o mkBlastDB_%j.log2

module purge
module load BLAST/2.9.0-gimkl-2018b

for filename in *.fna

do

makeblastdb -in ${filename} -dbtype nucl -parse_seqids -out ${filename}.db

done









