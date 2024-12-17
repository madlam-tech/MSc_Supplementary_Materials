#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J trycycler-ss-07
#SBATCH --time 00:20:00
#SBATCH --mem 40gb
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e tc_07_06_4_gen-ss.err
#SBATCH -o tc-ss.out
#SBATCH --export NONE

module purge
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0

trycycler subsample  --count 12 --genome_size 4.2 --min_read_depth 15 --reads barcode07.fastq --out_dir read_subsets
