#!/bin/bash -e
#SBATCH --account        massey03345
#SBATCH --job-name       flye12_i2
#SBATCH --mem            16G
#SBATCH --ntasks	2
#SBATCH --cpus-per-task  8
#SBATCH --time           6:00:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e flye%j.err


module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

mkdir -p flye_i5

flye --nano-hq   *.fastq \
	-o ./flye_i5 \
	-i 5 \
	-m 1000 \
