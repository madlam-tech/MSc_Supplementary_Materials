#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J spades
#SBATCH --time 00:15:00
#SBATCH --mem 15GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e sp-ass.err
#SBATCH -o sp-ass.out

module purge
module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5

for filename in *.fastq

do

spades.py --nanopore -k auto --threads 12 ./${filename} -o ./assemblies/${filename}spades.fasta

done

