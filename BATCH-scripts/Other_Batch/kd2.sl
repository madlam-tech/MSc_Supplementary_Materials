#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kneaddata
#SBATCH --time 0:30:00
#SBATCH --mem 12GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e kd2.err
#SBATCH -o kd2.out
#SBATCH --export NONE

module purge
module load Python/3.7.3-gimkl-2018b
module load Bowtie2/2.4.1-GCC-9.2.0
module load FastQC/0.11.9

kneaddata --input1 EO1219unmapped.bamR1.fq --input2 EO1219unmapped.bamR2.fq \
	--reference-db /nesi/nobackup/massey03345/bowtieDB \
	--bypass-trf \ 
	--run-fastqc-start \
	--run-fastqc-end \ 
	-o kd2b_out
