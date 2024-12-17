#!/bin/sh
#SBATCH --account massey03345
#SBATCH -J trycycler_
#SBATCH --time 00:30:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 1-trySRS_%j.log

module purge
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0
module load Filtlong/0.2.0


cat ./../*.fastq >> ./../cat.fq
mkdir -p ./../filtered
mkdir -p ./../read_subsets


for filename in ../filtered/*.fq

do

filtlong --min_length 1000 --keep_percent 95 ./../${filename} > ./../filtered/${filename}

trycycler subsample  --count 18 --min_read_depth 120 --reads ./../filtered/${filename} --out_dir ./../read_subsets

done
