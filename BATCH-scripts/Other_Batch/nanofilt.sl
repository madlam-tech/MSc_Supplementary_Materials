#!/bin/bash -e
#SBATCH --account        massey03345
#SBATCH --job-name       nanofilt
#SBATCH --mem            1G
#SBATCH --cpus-per-task  4
#SBATCH --time           00:10:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e nanofilt07_2022-08-19.err
#SBATCH -o nanofilt07_2022-08-19.out


module purge
module load nanofilt/2.6.0-gimkl-2020a-Python-3.8.2

for filename in *.fastq
do

filtlong --min_length 1000 --keep_percent 95 ${filename}.fastq > ${filename}.filt.fastq

done


