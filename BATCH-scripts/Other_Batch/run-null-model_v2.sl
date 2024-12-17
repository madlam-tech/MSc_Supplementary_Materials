#!/bin/bash

# SLURM job submission script for running Phyloseq Analysis

#SBATCH --account=massey03345       # Specify the account code
#SBATCH --partition=milan           # Specify the partition or queue
#SBATCH --job-name=PhyloseqAnalysis # Set a name for the job, visible in `squeue`
#SBATCH --time=2:00:00              # Set maximum time for the job hh:mm:ss
#SBATCH --cpus-per-task=100         # Number of cores to allocate
#SBATCH --mem=3G                    # Memory per node (all cores), adjust according to your job's requirements
#SBATCH --output=phyloseq_analysis_%j.out  # Standard output file
#SBATCH --error=phyloseq_analysis_%j.err   # Standard error file

# Clear all loaded modules to avoid conflicts
module purge

# Load the required module for the job
module load R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1

# Execute the R script
Rscript create_null_bNTI.Rscript

