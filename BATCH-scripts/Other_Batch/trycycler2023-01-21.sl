#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J mtest
#SBATCH --time 04:30:00
#SBATCH --mem 24GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e mtest.err
#SBATCH -o mtest.out
#SBATCH --export NONE

module purge
module load Filtlong/0.2.0
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0


mkdir filtered

for filename in *.fq

do


filtlong --min_length 1000 --keep_percent 95 ${filename} > ./filtered/${filename}

trycycler subsample  --count 10 --min_read_depth 10 --reads ./filtered/${filename} --out_dir ${filename}_read_subsets


done
