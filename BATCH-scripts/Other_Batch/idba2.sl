#!/bin/bash -e
#SBATCH -A massey03345
#SBATCH -J idbaud_assembly
#SBATCH --time 00:35:00
#SBATCH --mem 4GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e idbaud_assembly2.err
#SBATCH -o idbaud_assembly2.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load IDBA/1.1.3-gimkl-2017a

for filename in *.fna

do

	idba_ud --num_threads 8 --mink 33 --maxk 99 --step 22l -l ${filename} -o ${filename}_ass2

done





seqmagick convert --min-length 1000 EO1220unmapped_spades_assembly/EO1221_seq_magic.fna EO1220umnapped/EO1220_1000_seq_magic.fna
