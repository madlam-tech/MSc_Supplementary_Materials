#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J snakemake-ONT
#SBATCH --time 10:00:00
#SBATCH --mem 30GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 4
#SBATCH -e sm..err
#SBATCH -o sm.out

module purge
module load Miniconda3/4.9.2
module load snakemake/7.6.2-gimkl-2020a-Python-3.9.9

snakemake -s /nesi/nobackup/massey03345/nanopore/2022-09-23/2022-09-23/fastq/files/barcode1_snakemake/ont-assembly-snake/Snakefile \
	--use-mamba \
	--cores 20
	--config medaka_model=r941_min_high__g621
	filtlong_min_read_length=500
