#!/bin/bash -e

#SBATCH --account       massey03345      
#SBATCH --job-name      fq2fa_interleave      
#SBATCH --time          00:30:00          
#SBATCH --mem           2G           
#SBATCH --cpus-per-task 2 
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out 

module load SAMtools/1.9-GCC-7.4.0

for filename in _r1.fastq 

do 

    # Remove the _r1.fastq suffix and store the base name
    base_name="${filename%_r1.fastq}"

samtools merge ${base_name}_r1.fastq ${base_name}_r2.fastq ${filename}.fna

done


