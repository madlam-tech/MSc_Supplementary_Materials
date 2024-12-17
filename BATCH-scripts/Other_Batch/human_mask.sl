#!/bin/bash -e
#SBATCH -A massey03345
#SBATCH -J host_filt_bbmap_index
#SBATCH -J 2.qc_bbmap_ref
#SBATCH --time 00:20:00
#SBATCH --mem 23GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e host_filt_bbmap_index.err
#SBATCH -o host_filt_bbmap_index.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

cd /nesi/nobackup/massey03345/nanopore/2022-09-06/20220906_MA/20220906MA/20220906_1416_MN23427_FAV19227_07777e28/fastq/BBMask_human_reference/

# Load BBMap module
module purge
module load BBMap/38.90-gimkl-2020a

# Build indexed reference file via BBMap


srun bbmap.sh -usejni=t ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz -Xmx23
