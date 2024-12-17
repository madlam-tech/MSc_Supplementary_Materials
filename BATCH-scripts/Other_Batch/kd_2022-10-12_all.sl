#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kneaddata_EO5xx
#SBATCH --time 02:00:00
#SBATCH --mem 12GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e kd_EO5xx.err
#SBATCH -o kd_EO5xx.out
#SBATCH --export NONE

module purge
module load Python/3.7.3-gimkl-2018b
module load Bowtie2/2.4.1-GCC-9.2.0
module load FastQC/0.11.9

for filename *.fq

do

kneaddata --input1 /nesi/nobackup/massey03345/EO5xxPE/${filename}.fq \
	--input2 /nesi/nobackup/massey03345/EO5xxPE/${filename}.fq \
	--reference-db /nesi/nobackup/massey03345/bowtieDB \
	--bypass-trf \ 
	--run-fastqc-start \
	--run-fastqc-end \ 
	-o kd_2022-10-12_hu_plasmid_out

done
