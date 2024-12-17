#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J filtlong07
#SBATCH --time 00:15:00
#SBATCH --mem 1GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 1
#SBATCH -e FL2022-09-22.err
#SBATCH -o FL2022_09_22.out

module purge
module load Filtlong/0.2.0

for filename in *.fastq.gz

do

filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ${filename} | gzip > ../filtlong/${filename}output.gzip

done
