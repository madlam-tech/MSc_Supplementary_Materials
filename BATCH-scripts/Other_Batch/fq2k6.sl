#!/bin/sh -e
#SBATCH --account massey03345
#SBATCH -J k6
#SBATCH --time 03:00:00
#SBATCH --mem 65GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e k2_1.err
#SBATCH -o k1_1.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0
for filename in $(ls *paired_1.fastq | sed 's/paired_1.fastq//')

do 
    kraken2 --db /opt/nesi/db/Kraken2/standard-2022-07 --report ${filename}_k2_huplas_min3.tax --report-minimizer-data --minimum-hit-groups 3 --paired ${filename}paired_1.fastq ${filename}paired_2.fastq --output ${filename}_k2_hu_plas_min3.txt
done
