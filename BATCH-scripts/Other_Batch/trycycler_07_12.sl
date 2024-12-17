#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J trycycler-ss
#SBATCH --time 00:10:00
#SBATCH --mem 512mb
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e tc-ss07_12_10_5.2.err
#SBATCH -o tc-ss.out
#SBATCH --export NONE

module purge
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0

trycycler subsample --genome_size 5.2m --count 12 --min_read_depth 10 --reads barcode07.fq --out_dir read_subsets
