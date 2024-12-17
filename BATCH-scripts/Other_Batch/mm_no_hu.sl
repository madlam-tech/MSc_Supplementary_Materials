#!/bin/bash -e
#SBATCH --account        massey03345
#SBATCH --job-name       Minimap_human removal
#SBATCH --mem            2G
#SBATCH —ntasks 2
#SBATCH --cpus-per-task  8
#SBATCH --time           10:00:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e MM08.err
#SBATCH -o MM08.out


module purge
module load Nextflow/22.04.3


nextflow run phac-nml/ncov-dehoster -profile conda 
	--nanopore \ 
	--minimap2 \
	--fastq_directory /nesi/nobackup/massey03345/nanopore/2022-09-06/20220906_MA/20220906MA/20220906_1416_MN23427_FAV19227_07777e28/fastq/files/barcode08 \ 
	--human_ref /nesi/nobackup/massey03345/nanopore/2022-09-06/20220906_MA/20220906MA/20220906_1416_MN23427_FAV19227_07777e28/fastq/files/barcode08 \
	--run_name Barcode08_no_hu 
