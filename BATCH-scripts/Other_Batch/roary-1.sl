#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 28-roary-1
#SBATCH --time 2:00:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e roary-1_%j.log


module purge
module load Roary/3.13.0-gimkl-2020a

roary -r -s -p 8 -i 70 *.gff
