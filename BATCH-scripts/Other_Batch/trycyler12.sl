#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J trycy12
#SBATCH --time 00:15:00
#SBATCH --mem 2GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e tc-ss12.err
#SBATCH -o tc-ss12.out
#SBATCH --export NONE

module purge
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0

trycycler subsample --genome_size 5m --count 12 --min_read_depth 10 --reads ./*fq --out_dir read_subsets
