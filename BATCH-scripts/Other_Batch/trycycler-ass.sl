#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J trycycler_ass
#SBATCH --time 02:00:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e tc-ass.err
#SBATCH -o tc-ass.out

module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5
module load Raven/1.5.0-GCC-9.2.0

mkdir barcode01_flye
mkdir barcode01_raven

for filename in *.fastq

flye --nano-hq   ${filename}.fastq —threads 16 -o ./barcode01_flye\${filename} -i 2 -m 1000

raven --threads 16 ${filename} > ./barcode01_raven\${filename} 
