#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 1-flye2racon.sl
#SBATCH --time 8:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e flye2racon%j.err
#SBATCH -o flye2racon%j.out
#SBATCH --export NONE

mkdir -p filtlong

module load Filtlong/0.2.0
module load minimap2/2.24-GCC-9.2.0

for filename in ../BC06_07_09/*.fq

filtlong --min_length 1000 --keep_percent 95 ${filename} > ./filtlong/${filename}

minimap2 –x ava-ont \
 ./filtered/*fq \
| gzip -1 > ./minimap.paf.gz


racon ./filtered/*fq
./minimap.paf.gz \
 ./assembly.fasta \ 
> ./miniasm.racon.consensus.fasta
