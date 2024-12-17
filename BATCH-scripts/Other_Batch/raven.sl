#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J raven
#SBATCH --time 00:25:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e tc-ass.err
#SBATCH -o tc-ass.out

module purge
module load Raven/1.5.0-GCC-9.2.0


for filename in *.fastq

do

raven --threads 16 ./${filename} > ./assemblies/${filename}raven.fasta && rm raven.cereal

done

