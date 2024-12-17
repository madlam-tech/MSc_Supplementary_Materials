#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kneaddata
#SBATCH --time 0:30:00
#SBATCH --mem 12GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e bt2.err
#SBATCH -o bt2.out
#SBATCH --export NONE

module purge
module load Python/3.7.3-gimkl-2018b
module load Bowtie2/2.4.1-GCC-9.2.0
module load FastQC/0.11.9

kneaddata --unpaired BC10.fq \
	--reference-db ./bowtieDB \
	--bypass-trf \ 
	--run-fastqc-start \
	--run-fastqc-end \ 
	--output BC10_decon_out_bt
	--output-prefix BC10_decon_bt2
