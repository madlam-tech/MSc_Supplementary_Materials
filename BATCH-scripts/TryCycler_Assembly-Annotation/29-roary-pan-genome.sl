#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 29-pan-roary
#SBATCH --time 01:00:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e 29-pan-roary_%j.log


module purge
module load Roary/3.13.0-gimkl-2020a

query_pan_genome -a intersection *gff

mkdir paraburk_core

sed 's/:\s/\t/g' pan_genome_results > paraburk_core/core_genes.table
cut -f 1 paraburk_core/core_genes.table > paraburk_core/core_genes.list
