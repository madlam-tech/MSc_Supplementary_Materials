#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J bowtie-check
#SBATCH --time 00:10:00
#SBATCH --mem 2GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e b2-check.err
#SBATCH -o b2-check.out
#SBATCH --export NONE

module purge

module load Bowtie2/2.4.5-GCC-11.3.0

	bowtie2-inspect -a -n -s -o bt2_DB_inspect /nesi/nobackup/massey03345/nanopore/bowtieDB
