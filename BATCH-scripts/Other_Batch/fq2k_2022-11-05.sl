#!/bin/sh -e
#SBATCH --account massey03345
#SBATCH -J k2_2022-11-05
#SBATCH --time 05:00:00
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
    kraken2 --db /opt/nesi/db/Kraken2/standard-2022-07 --threads 8 --report ${filename}.k2report \
	--report-minimizer-data --minimum-hits-groups 3 \
	${filename}paired_1.fastq ${filename}paired_2.fastq \
	--output ${filename}_k2.txt
done
