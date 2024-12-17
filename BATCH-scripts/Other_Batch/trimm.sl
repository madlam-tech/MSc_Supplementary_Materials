#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J trimmomatic
#SBATCH --time 0:15:00
#SBATCH --mem 2GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e trimmomatic.err
#SBATCH -o trimmomatic.out

module purge
module load Trimmomatic/0.39-Java-1.8.0_144

for filename in *_r1.fastq 
do 
    # Remove the _r1.fastq suffix and store the base name
    base_name="${filename%_r1.fastq}"
    
    trimmomatic PE -threads 2 -phred33 ${filename} ${base_name}_r2.fastq \
        ${base_name}_r1.qc.fastq ${base_name}_s1.qc.fastq ${base_name}_r2.qc.fastq.gz ${base_name}_s2.qc.fastq.gz \
        HEADCROP:10 SLIDINGWINDOW:4:30 MINLEN:100
done



