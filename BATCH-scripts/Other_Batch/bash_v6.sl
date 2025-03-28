#!/bin/sh
#SBATCH --account=massey03345
#SBATCH --partition=milan
#SBATCH --job-name=PhyloseqAnalysis
#SBATCH --time=2:00:00            
#SBATCH --cpus-per-task=10     
#SBATCH --mem=25G                  # Increased from 4G to 25G
#SBATCH --output=phyloseq_analysis_%j.out
#SBATCH --error=phyloseq_analysis_%j.err

module purge
module load R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1

Rscript create_null_bNTI_v6.Rscript

