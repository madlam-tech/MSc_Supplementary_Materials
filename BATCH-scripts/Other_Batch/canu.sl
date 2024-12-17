#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J canu_v2
#SBATCH --mem 20GB
#SBATCH --nodes 1
#SBTACH --cpus-per-task 2
#SBATCH --ntasks-per-node 24
#SBATCH --time 12:00:00
#SBATCH --output canu%j.log

module purge
module load Canu/2.2-GCC-9.2.0

canu \
  -p BC11_canu \
  -d BC11_canu_out \
  -genomeSize=8.8m \
  -nanopore-raw ./*.fq \
  -minReadLength=5000 \
  -minOverlapLength=1000 \
  -corOutCoverage=200 \
  -maxThreads=24 \
  -useGrid=false \
  -maxMemory=20G
