#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J fastqc check
#SBATCH --time 01:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e fastqc.err
#SBATCH -o fastqc.out
#SBATCH --export NONE



module purge
module load FastQC/0.11.7

do

fastqc mock_R1.good.fastq.gz mock_R2.good.fastq.gz

done