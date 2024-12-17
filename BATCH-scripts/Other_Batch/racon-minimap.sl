#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J Rac-MMPsl
#SBATCH --time 00:20:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e Rac-MMP%j.err
#SBATCH --export NONE

threads=32

module purge
module load Python/3.10.5-gimkl-2022a
module load any2fasta/0.4.2-GCC-9.2.0
module load minimap2/2.24-GCC-9.2.0
module load Raven/1.5.0-GCC-9.2.0
module load Racon/1.5.0-GCC-11.3.0
module load miniasm/0.3-20191007-GCC-11.3.0

./miniasm_and_minipolish.sh ./BC10_filt.fq "$threads" > assemblies/BC10_filt.gfa && any2fasta assemblies/BC10_filt.gfa > assemblies/BC10_filt.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/BC10_filt.gfa ./BC10_filt.fq > assemblies/BC10_filt.fasta




