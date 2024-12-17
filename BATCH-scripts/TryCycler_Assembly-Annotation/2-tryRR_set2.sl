#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 2-tryRR_08-15.sl
#SBATCH --time 2:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e 2-tryRR_08-15%j.err
#SBATCH --export NONE

threads=32

module purge
module load Python/3.10.5-gimkl-2022a
module load any2fasta/0.4.2-GCC-9.2.0
module load minimap2/2.24-GCC-9.2.0
module load Raven/1.5.0-GCC-9.2.0
module load Racon/1.5.0-GCC-11.3.0
module load miniasm/0.3-20191007-GCC-11.3.0

./miniasm_and_minipolish.sh ./sample_08.fastq "$threads" > assemblies/assembly_08.gfa && any2fasta assemblies/assembly_08.gfa > assemblies/assembly_08.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_09.gfa ./sample_09.fastq > assemblies/assembly_09.fasta

./miniasm_and_minipolish.sh ./sample_10.fastq "$threads" > assemblies/assembly_10.gfa && any2fasta assemblies/assembly_10.gfa > assemblies/assembly_10.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_11.gfa ./sample_11.fastq > assemblies/assembly_11.fasta

./miniasm_and_minipolish.sh ./sample_12.fastq "$threads" > assemblies/assembly_12.gfa && any2fasta assemblies/assembly_12.gfa > assemblies/assembly_12.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_13.gfa ./sample_13.fastq > assemblies/assembly_13.fasta

./miniasm_and_minipolish.sh ./sample_14.fastq "$threads" > assemblies/assembly_14.gfa && any2fasta assemblies/assembly_14.gfa > assemblies/assembly_14.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_15.gfa ./sample_15.fastq > assemblies/assembly_15.fasta
