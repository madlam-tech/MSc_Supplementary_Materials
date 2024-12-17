#!/bin/bash

#SBATCH --account=uoa04147
#SBATCH --partition=milan
#SBATCH --job-name=PhyloseqAnalysis
#SBATCH --time=2:00:00            
#SBATCH --mail-type=BEGIN,END,FAIL  # Notifications for job start, end, and failure
#SBATCH --mail-user=madlam@me.com
#SBATCH --cpus-per-task=50     
#SBATCH --mem=80G                  # Increased from 1G to 4G
#SBATCH --output=phyloseq_analysis_50_cores_%j.out
#SBATCH --error=phyloseq_analysis_50_cores_%j.err

module purge
module load R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1

Rscript create_null_bNTI_v9.Rscript
