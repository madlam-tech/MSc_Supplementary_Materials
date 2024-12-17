#!/bin/bash
#SBATCH --account=massey03345
#SBATCH --job-name=clustal-omega-job
#SBATCH --output=clustal-omega-output_%j.txt
#SBATCH --error=clustal-omega-error_%j.txt
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=6G



# Load the Clustal Omega module
module purge
module load Clustal-Omega/1.2.4-gimkl-2020a



# Run Clustal Omega
clustalo --auto  -i cat.fasta -o Pseudomonas.aln --outfmt=clu --guidetree-out=Pseudomonas_guidetree.dnd --clustering-out=Pseudomonas_cluster_file --force

echo "Clustal Omega alignment completed."

