#!/bin/bash
#SBATCH --account=massey03345
#SBATCH --job-name=clustal-omega-job
#SBATCH --output=clustal-omega-output_%j.txt
#SBATCH --error=clustal-omega-error_%j.txt
#SBATCH --time=1:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=5G



# Load the Clustal Omega module
module purge
module load Clustal-Omega/1.2.4-gimkl-2020a



# Run Clustal Omega
clustalo --auto  -i cat.fasta -o Staph.aln --outfmt=clu --guidetree-out=Staph_guidetree.dnd --clustering-out=Staph_cluster_file --force

echo "Clustal Omega alignment completed."

