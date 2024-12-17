#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 2-tryRR_16-24.sl
#SBATCH --time 2:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e 2-tryRR_16-24%j.err
#SBATCH --export NONE

threads=32

module purge
module load Python/3.10.5-gimkl-2022a
module load any2fasta/0.4.2-GCC-9.2.0
module load minimap2/2.24-GCC-9.2.0
module load Raven/1.5.0-GCC-9.2.0
module load Racon/1.5.0-GCC-11.3.0
module load miniasm/0.3-20191007-GCC-11.3.0

./miniasm_and_minipolish.sh ./sample_16.fastq "$threads" > assemblies/assembly_16.gfa && any2fasta assemblies/assembly_16.gfa > assemblies/assembly_16.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_17.gfa ./sample_17.fastq > assemblies/assembly_17.fasta

./miniasm_and_minipolish.sh ./sample_18.fastq "$threads" > assemblies/assembly_18.gfa && any2fasta assemblies/assembly_18.gfa > assemblies/assembly_18.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_19.gfa ./sample_19.fastq > assemblies/assembly_19.fasta

./miniasm_and_minipolish.sh ./sample_20.fastq "$threads" > assemblies/assembly_20.gfa && any2fasta assemblies/assembly_20.gfa > assemblies/assembly_20.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_21.gfa ./sample_21.fastq > assemblies/assembly_21.fasta

./miniasm_and_minipolish.sh ./sample_22.fastq "$threads" > assemblies/assembly_22.gfa && any2fasta assemblies/assembly_22.gfa > assemblies/assembly_22.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_23.gfa ./sample_23.fastq > assemblies/assembly_23.fasta

./miniasm_and_minipolish.sh ./sample_24.fastq "$threads" > assemblies/assembly_24.gfa && any2fasta assemblies/assembly_24.gfa > assemblies/assembly_24.fasta


