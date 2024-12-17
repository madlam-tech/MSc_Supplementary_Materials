#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J trycycler-ss-03
#SBATCH --time 00:20:00
#SBATCH --mem 2gb
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e tc_09.err
#SBATCH -o tc_09.out
#SBATCH --export NONE

module purge
module load Filtlong/0.2.0
module load Trycycler/0.4.2-gimkl-2020a-Python-3.8.2
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0


filtlong --min_length 1000 --keep_percent 95 BC09.fq > BC09_fl.fq


trycycler subsample  -r BC09_fl.fq -o ../BC09_read_subsets \
		--genome_size 7m \
		--count 12 \
		--min_read_depth 20 \
		-t 12
