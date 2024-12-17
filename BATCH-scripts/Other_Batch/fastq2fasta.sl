#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J fastq2fasta
#SBATCH --time 00:015:00
#SBATCH --mem 15GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e fastq2fasta.err
#SBATCH -o fastq2fasta.out
#SBATCH --export NONE

module purge
module load seqmagick/0.7.0-gimkl-2018b-Python-3.7.3

for filename in *.fq

do 

sed -n '1~4s/^@/>/p;2~4p' ${filename} > ${filename}.fasta
    
done
