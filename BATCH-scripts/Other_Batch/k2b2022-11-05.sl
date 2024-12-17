#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J bracken
#SBATCH --time 00:20:00
#SBATCH --mem 5GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 1
#SBATCH -e bracken.err
#SBATCH -o bracken.out
#SBATCH --export NONE

module purge
module load Bracken/2.6.2-GCCcore-9.2.0

for filename in *.k2

do 

    bracken -d /opt/nesi/db/Kraken2/standard-2018-09 -i ${filename}.k2 -o ${filename}.bracken -r 150 -l S -t 1
    
done 

