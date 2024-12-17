#!/bin/bash -e
#SBATCH --account        massey03345
#SBATCH --job-name       Bowtie_build
#SBATCH --mem            10G
#SBATCH --cpus-per-task  2
#SBATCH --time           05:00:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e BT_build.err
#SBATCH -o BT_build.out


module purge
module load Bowtie2/2.4.5-GCC-11.3.0



bowtie2-build --large-index -f ./ref_gen/*.fna ./bowtieDB


