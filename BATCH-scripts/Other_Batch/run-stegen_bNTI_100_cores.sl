#!/bin/bash
#SBATCH --account=uoa04147
#SBATCH --partition=milan
#SBATCH --job-name=PhyloseqAnalysis
#SBATCH --time=2:00:00            
#SBATCH --cpus-per-task=100     
#SBATCH --mem=10G                  # Increased from 1G to 10G
#SBATCH --output=phyloseq_analysis_%j.out
#SBATCH --error=phyloseq_analysis_%j.err

module purge
module load R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1

Rscript stegen_NTI_100_cores_v2.Rscript

