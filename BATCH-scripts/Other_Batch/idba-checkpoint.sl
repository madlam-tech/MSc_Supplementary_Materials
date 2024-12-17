#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J idbaud_assembly
#SBATCH --time 00:35:00
#SBATCH --mem 4GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e idbaud_assembly.err
#SBATCH -o idbaud_assembly.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load IDBA/1.1.3-gimkl-2017a

for i in *.fna

do

	idba_ud --num_threads 8 --mink 33 --maxk 99 --step 22 -r {i}.fna -o {i}ass/

done
