#!/bin/bash
#SBATCH --account=uoa04147
#SBATCH --partition=milan
#SBATCH --job-name=PhyloseqAnalysis
#SBATCH --time=2:00:00            
#SBATCH --cpus-per-task=10     
#SBATCH --mem=40G                  # Increased from 1G to 4G
#SBATCH --output=stegen_%j.out
#SBATCH --error=stegen_%j.err

module purge
module load R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1

Rscript stegen_NTI.Rscript
