#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J trycycler-ss-03
#SBATCH --time 00:20:00
#SBATCH --mem 6gb
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e tc_07_06_no_gen-ss.err
#SBATCH -o tc-ss.out
#SBATCH --export NONE

module purge
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0

trycycler subsample  --count 12 --min_read_depth 10 --reads barcode09.fastq --out_dir read_subsets
