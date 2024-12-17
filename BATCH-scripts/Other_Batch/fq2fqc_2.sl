#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J fastqc
#SBATCH --time 01:00:00
#SBATCH --mem 12GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e 2_fastqc.err
#SBATCH -o 2_fastqc.out

module load FastQC/0.12.1

# Use find to locate all paired-end fastq files in the specified directory
for filename in /nesi/nobackup/massey03345/CD_clones/*_R1.fastq
do
  # Get the file name without the path and extension
  base_filename=$(basename -- "$filename")
  filename_noext="${base_filename%_R1.fastq}"
  
  # Run FastQC on the R1 and R2 files
  fastqc -o /nesi/nobackup/massey03345/CD_clones/"${filename_noext}.report" "$filename"
  fastqc -o /nesi/nobackup/massey03345/CD_clones/"${filename_noext}.report" "/nesi/nobackup/massey03345/CD_clones/${filename_noext}_R2.fastq"
done

