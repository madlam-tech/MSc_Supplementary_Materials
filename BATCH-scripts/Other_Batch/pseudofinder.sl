#!/bin/bash -e

#SBATCH --account massey03345
#SBATCH -J PsFinder-diamondDB
#SBATCH --time 08:00:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 23-pseudogenes.log

module load DIAMOND/2.1.6-GCC-11.3.0
module load BLASTDB/2024-04
module load BLAST/2.13.0-GCC-11.3.0
module load Python/3.10.5-gimkl-2022a

for filename in *gbk

do

basename=$(basename ${filename} .gbk)
mkdir ${basename}

./pseudofinder/pseudofinder.py annotate --diamond --skip_makedb -g ${filename} -db  /nesi/nobackup/massey03345/uniprot -op ${basename}



done
