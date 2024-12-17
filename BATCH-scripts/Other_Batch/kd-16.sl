#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kneaddata
#SBATCH --time 06:00:00
#SBATCH --mem 30GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e kd2-16.err
#SBATCH -o kd2-16.out
#SBATCH --export NONE

module purge
module load Python/3.7.3-gimkl-2018b
module load Bowtie2/2.4.1-GCC-9.2.0
module load FastQC/0.11.9

kneaddata --unpaired BC10.fq --reference-db ./bowtieDB --output BC10_decon_out_kd --bypass-trim --bypass-trf --threads 2
