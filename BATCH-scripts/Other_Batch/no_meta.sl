#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      meta_strip
#SBATCH --time          00:20:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 1
#SBATCH --error         No_meta.err
#SBATCH --output        No_meta.out


for pred_file in predictions/*.fna;
do
    file_base=$(basename ${pred_file} .fna)

    cut -f1 -d ' ' predictions/${file_base}.fna > predictions/${file_base}.no_metadata.fna
    cut -f1 -d ' ' predictions/${file_base}.faa > predictions/${file_base}.no_metadata.faa
done
