#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J try12_NG
#SBATCH --time 01:00:00
#SBATCH --mem 16GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e tc-ss_NG.err
#SBATCH -o tc-ss_NG.out
#SBATCH --export NONE

module purge
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0

trycycler subsample  --count 12 --min_read_depth 10 --reads ./*fq --out_dir read_subsets_NG
