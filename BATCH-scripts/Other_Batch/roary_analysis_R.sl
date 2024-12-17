#!/bin/bash -e
#SBATCH --account=massey03345
#SBATCH --job-name=R-roary
#SBATCH --time=2:00:00
#SBATCH --mem=48GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --error=R_roary_%j.err
#SBATCH --output=R_roary_%j.log

module purge
module load R/4.3.2-foss-2023a
module load R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1

Rscript R_pseudo_roary.R
