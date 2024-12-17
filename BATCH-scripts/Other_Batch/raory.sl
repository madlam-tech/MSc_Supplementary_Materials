#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J pan-roary
#SBATCH --time 10:00:00
#SBATCH --mem 12GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 12
#SBATCH -e 29-pan-roary_%j.log


module purge
module load Roary/3.13.0-gimkl-2020a

roary -f ./prokka_4_pseudo -e -n -v ./gffs_roary/*gff
