#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 29-roary-2
#SBATCH --time 20:00:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e 29-roary-2_%j.log


module purge
module load Roary/3.13.0-gimkl-2020a

roary -f ./ -e -n -v ./*.gff
