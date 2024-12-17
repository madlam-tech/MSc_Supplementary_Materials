#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 28-roary
#SBATCH --time 02:00:00
#SBATCH --mem 1GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 8
#SBATCH -e 28-roary-alone_%j.log


module purge
module load Roary/3.13.0-gimkl-2020a

roary -r -s -p 8 -i 70 *.gff
