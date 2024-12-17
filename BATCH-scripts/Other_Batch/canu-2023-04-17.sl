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

for filepath in ./*fastq; do
    filename=$(basename "$filepath")
    filename="${filename%.*}"
    mkdir -p  "./$filename"

canu \
  -p $filename \
  -d ./$filename \
  -genomeSize=5.3m \
  -nanopore-raw ./${filename}.fastq \
  -minReadLength=1000 \
  -minOverlapLength=1000 \
  -corOutCoverage=200 \
  -maxThreads=24 \
  -useGrid=false \
  -maxMemory=20G
done

