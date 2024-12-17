#!/bin/sh
#SBATCH --account        massey03345
#SBATCH --job-name       ont-guppy-gpu_job
#SBATCH --gpus-per-node  P100:1
#SBATCH --mem            12G
#SBATCH --cpus-per-task  4
#SBATCH --time           10:00:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e guppy.err
#SBATCH -o guppy.out


module purge
module load medaka/1.6.0-Miniconda3-4.12.0


medaka_consensus -t 14 -m r941_min_hac_g507 \
	-i /nesi/nobackup/massey03345/nanopore/2022-09-23/2022-09-23/fastq/files/barcode01/barcode01.fq \
	-d /nesi/nobackup/massey03345/nanopore/2022-09-23/2022-09-23/fastq/files/barcode01_flye_i5_gen_5.5/assembly.fasta \
	-o /nesi/nobackup/massey03345/nanopore/2022-09-23/2022-09-23/fastq/files/barcode01/medaka
