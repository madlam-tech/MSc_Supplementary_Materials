#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kneaddata_EO5xx_hu_plas
#SBATCH --time 02:00:00
#SBATCH --mem 12GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e kd_EO5xx_hu_plas.err
#SBATCH -o kd_EO5xx_hu_plas.out


module purge
module load Python/3.7.3-gimkl-2018b
module load Bowtie2/2.4.1-GCC-9.2.0
module load FastQC/0.11.9

for SAMPLE in $(ls *.bamR1.fq | sed 's/.bamR1.fq//')
do

kneaddata --input1 ${SAMPLE}.bamR1.fq --input2 ${SAMPLE}.bamR2.fq \
	--bypass-trf \
	-o kd2b_hu_plas_out \
	-db /nesi/nobackup/massey03345/bowtieDB \
	--run-fastqc-start \
	--run-fastqc-end  

done



