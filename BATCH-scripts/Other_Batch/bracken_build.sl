#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J bracken-build
#SBATCH --time 00:02:00
#SBATCH --mem 512mB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e BrBuild.err
#SBATCH -o BrBuild.out
#SBATCH --export NONE

module purge
module load Bracken/2.7-GCC-11.3.0
module load Kraken2/2.1.2-GCC-9.2.0

bracken-build -d /opt/nesi/db/Kraken2/standard-2022-07 -t 1 -k 35 -l 150
